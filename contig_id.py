"""Contig species identification and confirmed abundance for the eToL workflow.

After :mod:`assembler` builds contigs from each taxon's matched reads, the paper
(Hu, Haas & Lathe, *BMC Microbiology* 2022;22:317) does two things with them:

1. **Species identification** -- BLAST each contig against a comprehensive rRNA
   reference (the paper used "BLAST at NCBI"; COBLAST stays local and uses a
   registered reference DB, e.g. a SILVA SSU/LSU rRNA database) and take the
   closest homolog as the species call.
2. **Confirmed abundance** -- BLAST each contig against *that taxon's own*
   matched reads, keeping only reads that align at near-100% identity, and count
   them. This re-confirms abundance against the contig rather than the 64-mer
   probe.

The paper's third step, re-probing the library with key contigs, is not
implemented here (deferred by design).

Both BLAST calls reuse the local BLAST+ install. Species naming searches an
indexed reference DB (``-db``); confirmation is a small bl2seq-style search of a
taxon's reads against its handful of contigs (``-subject``, no index needed).
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
from typing import Any

from blast_runner import reads_to_fasta
from config import blast_exe
from database_registry import blast_safe_path
from human_filter import (
    HUMAN_BITSCORE_THRESHOLD,
    extract_reads,
    find_human_read_ids,
)


# Reads aligning to a contig at or above this percent identity are "confirmed"
# (the paper keeps "100% or near-100% identity"). 99% tolerates a single
# sequencing error over a short rRNA read.
DEFAULT_CONFIRM_IDENTITY_PCT = 99.0
# Species naming is a real homology search, so keep a modest e-value gate; the
# confirmation search is identity-gated in Python, so its e-value stays loose.
DEFAULT_NAME_EVALUE = "1e-5"
DEFAULT_TIMEOUT_SECONDS = 1800

# Re-probing (Hu, Haas & Lathe 2022, Box 3): use a taxon's most-abundant contig
# as a fresh probe to pull more reads from the same library, then re-assemble.
# The paper's example used the "two most abundant" contigs; COBLAST uses the
# single most-abundant one because each extra contig is a full-length, highly-
# conserved rRNA query against the whole patient library -- the dominant cost of
# this step. ponytail: 1 contig; raise toward 2 only if recovery sensitivity
# measurably suffers (after a permissive net first pass it recovers ~0 either way).
REPROBE_TOP_CONTIGS = 1
# Re-probe matches are gated on the same E-value as the net (drop E >= this), so
# re-probing is no more permissive than the first search that found the taxon.
REPROBE_EVALUE = 0.01
# A contig can match many reads; lift the default cap so read recovery isn't
# silently truncated. ponytail: 100k covers ordinary SRAs; raise only if a single
# contig legitimately matches more reads than this in one library.
REPROBE_MAX_TARGET_SEQS = "100000"


# Reference hit titles look like ``ACC.start.end Domain;Phylum;...;Species``; the
# species is the last ``;``-delimited rank (the leading accession lives in the
# first field, so it never reaches the last). Human reads from rRNA can slip the
# genome-level human filter and assemble into contigs the reference correctly
# calls Homo sapiens -- the homolog string is the only signal that catches them.
HUMAN_HOMOLOG_MARKER = "Homo sapiens"


def species_from_homolog(stitle: str) -> str:
    """Reduce a reference hit title to just its species name.

    Takes the last non-empty ``;``-delimited rank, so ``ACC.s.e Bacteria;...;
    Pseudomonas fluorescens`` becomes ``Pseudomonas fluorescens``. A title with
    no lineage (no ``;``) is returned as-is, so this never raises.
    """
    ranks = [part.strip() for part in (stitle or "").split(";") if part.strip()]
    return ranks[-1] if ranks else (stitle or "").strip()


def _is_human_homolog(stitle: str) -> bool:
    """True when a contig's closest homolog is human (checked over the lineage)."""
    return HUMAN_HOMOLOG_MARKER in (stitle or "")


def _best_homolog_per_query(tabular_text: str) -> dict[str, dict[str, Any]]:
    """Pick the highest-bitscore reference hit per contig from BLAST outfmt 6.

    Expects ``6 qseqid pident length bitscore stitle`` rows. ``stitle`` is the
    last field (it may contain spaces, never tabs), so the split keeps it whole.
    """
    best: dict[str, dict[str, Any]] = {}
    for line in tabular_text.splitlines():
        fields = line.split("\t", 4)
        if len(fields) < 5:
            continue
        qseqid = fields[0]
        try:
            pident = float(fields[1])
            length = int(fields[2])
            bitscore = float(fields[3])
        except ValueError:
            continue
        stitle = fields[4].strip()
        current = best.get(qseqid)
        if current is None or bitscore > current["bitscore"]:
            best[qseqid] = {
                "homolog": stitle,
                "pident": pident,
                "length": length,
                "bitscore": bitscore,
            }
    return best


def _confirmed_reads_from_tabular(
    tabular_text: str, identity_pct: float
) -> tuple[dict[str, int], set[str]]:
    """Count reads (qseqid) confirming each contig (sseqid) at >= identity_pct.

    Expects ``6 qseqid sseqid pident`` rows. Returns per-contig confirmed read
    counts and the set of all distinct confirmed reads (so a taxon-level total
    counts a read once even if it matches several of that taxon's contigs).
    """
    per_contig: dict[str, set[str]] = {}
    confirmed: set[str] = set()
    for line in tabular_text.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        read_id, contig_id, pident_str = fields[0], fields[1], fields[2]
        try:
            pident = float(pident_str)
        except ValueError:
            continue
        # Tiny epsilon so 99.0 passes a >=99 gate despite float formatting.
        if pident + 1e-9 >= identity_pct:
            per_contig.setdefault(contig_id, set()).add(read_id)
            confirmed.add(read_id)
    return {cid: len(reads) for cid, reads in per_contig.items()}, confirmed


def name_contigs(
    named_contigs: dict[str, str],
    reference_db_prefix: str,
    *,
    evalue: str = DEFAULT_NAME_EVALUE,
    num_threads: int | str | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
) -> dict[str, dict[str, Any]]:
    """BLAST contigs against the reference DB; return the best homolog per id.

    ``named_contigs`` maps a caller-chosen unique id to a contig sequence (the
    caller batches every taxon's contigs into one search and maps ids back). The
    returned dict is keyed by the same ids.
    """
    if not named_contigs:
        return {}
    with tempfile.TemporaryDirectory(prefix="contig_name_") as tmpdir:
        query_path = Path(tmpdir) / "contigs.fasta"
        query_path.write_text(reads_to_fasta(named_contigs), encoding="utf-8")
        command = [
            str(blast_exe("blastn")),
            "-task",
            "megablast",
            "-query",
            str(query_path),
            "-db",
            blast_safe_path(reference_db_prefix),
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            "5",
            "-outfmt",
            "6 qseqid pident length bitscore stitle",
        ]
        if num_threads:
            command += ["-num_threads", str(num_threads)]
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "Contig species-ID BLAST failed: "
            + ((completed.stderr or completed.stdout).strip() or "unknown error")
        )
    return _best_homolog_per_query(completed.stdout)


def confirm_contig_reads(
    contigs: list[dict[str, Any]],
    reads: dict[str, str],
    *,
    identity_pct: float = DEFAULT_CONFIRM_IDENTITY_PCT,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
) -> tuple[dict[str, int], set[str]]:
    """Count a taxon's reads that align to its contigs at >= identity_pct.

    Contigs are the BLAST subjects (a handful of short sequences) and the reads
    are the queries, so no makeblastdb is needed.
    """
    subjects = {
        str(contig.get("id", "")): str(contig.get("sequence", ""))
        for contig in contigs
        if contig.get("id") and contig.get("sequence")
    }
    if not subjects or not reads:
        return {}, set()
    with tempfile.TemporaryDirectory(prefix="contig_confirm_") as tmpdir:
        reads_path = Path(tmpdir) / "reads.fasta"
        contigs_path = Path(tmpdir) / "contigs.fasta"
        reads_path.write_text(reads_to_fasta(reads), encoding="utf-8")
        contigs_path.write_text(reads_to_fasta(subjects), encoding="utf-8")
        # ponytail: -subject runs single-threaded (BLAST ignores -num_threads
        # with a subject file). Fine for a few contigs x a few thousand reads;
        # build a real DB per taxon only if a taxon's read count makes this slow.
        command = [
            str(blast_exe("blastn")),
            "-task",
            "megablast",
            "-query",
            str(reads_path),
            "-subject",
            str(contigs_path),
            "-outfmt",
            "6 qseqid sseqid pident",
        ]
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "Confirmed-abundance BLAST failed: "
            + ((completed.stderr or completed.stdout).strip() or "unknown error")
        )
    return _confirmed_reads_from_tabular(completed.stdout, identity_pct)


def _representative_contig(contigs: list[dict[str, Any]]) -> dict[str, Any] | None:
    """The taxon's most-supported contig (most reads, then longest)."""
    named = [contig for contig in contigs if contig.get("closest_homolog")]
    pool = named or contigs
    if not pool:
        return None
    return max(pool, key=lambda c: (int(c.get("num_reads", 0)), int(c.get("length", 0))))


def identify_contigs(
    contigs_by_species: dict[str, list[dict[str, Any]]],
    reads_by_taxon: dict[str, list[str]],
    all_reads: dict[str, str],
    *,
    reference_db_prefix: str,
    identity_pct: float = DEFAULT_CONFIRM_IDENTITY_PCT,
    num_threads: int | str | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
) -> dict[str, dict[str, Any]]:
    """Annotate contigs with their closest homolog + confirmed read counts.

    Mutates each contig dict in ``contigs_by_species`` in place (adds
    ``closest_homolog``, ``homolog_pident``, ``confirmed_reads``) and returns a
    per-taxon summary ``{taxon: {closest_homolog, homolog_pident,
    confirmed_reads}}`` for the species table and exports.
    """
    # 1. Species naming: one batched BLAST of every contig against the reference
    #    DB (one DB open beats one search per taxon). Synthetic ids map back to
    #    (taxon, contig) because CAP3 restarts its "Contig1.." numbering per taxon.
    named: dict[str, str] = {}
    origin: dict[str, tuple[str, int]] = {}
    for taxon, contigs in contigs_by_species.items():
        for index, contig in enumerate(contigs):
            sequence = str(contig.get("sequence", ""))
            if not sequence:
                continue
            synthetic_id = f"c{len(named)}"
            named[synthetic_id] = sequence
            origin[synthetic_id] = (taxon, index)
    homologs = name_contigs(
        named,
        reference_db_prefix,
        num_threads=num_threads,
        timeout_seconds=timeout_seconds,
    )
    for synthetic_id, info in homologs.items():
        taxon, index = origin[synthetic_id]
        contig = contigs_by_species[taxon][index]
        contig["closest_homolog"] = info["homolog"]
        contig["closest_species"] = species_from_homolog(info["homolog"])
        contig["homolog_pident"] = info["pident"]

    # 1b. Drop human contigs: human rRNA reads can survive the genome-level human
    #     filter and assemble here; the reference homolog is the only signal that
    #     identifies them. Removing them from contigs_by_species also shrinks the
    #     contig_count, exports and FASTA download. ponytail: substring match on
    #     the homolog lineage; tighten to a taxonomy field only if it misfires.
    for taxon, contigs in contigs_by_species.items():
        contigs_by_species[taxon] = [
            contig
            for contig in contigs
            if not _is_human_homolog(str(contig.get("closest_homolog", "")))
        ]

    # 2. Confirmed abundance: per taxon, BLAST that taxon's reads against its own
    #    contigs and count reads at >= identity_pct.
    identification: dict[str, dict[str, Any]] = {}
    for taxon, contigs in contigs_by_species.items():
        reads = {
            read_id: all_reads[read_id]
            for read_id in reads_by_taxon.get(taxon, [])
            if read_id in all_reads
        }
        per_contig, confirmed = confirm_contig_reads(
            contigs, reads, identity_pct=identity_pct, timeout_seconds=timeout_seconds
        )
        for contig in contigs:
            contig["confirmed_reads"] = per_contig.get(str(contig.get("id", "")), 0)
        representative = _representative_contig(contigs)
        identification[taxon] = {
            "closest_homolog": (representative or {}).get("closest_homolog", ""),
            "closest_species": (representative or {}).get("closest_species", ""),
            "homolog_pident": (representative or {}).get("homolog_pident"),
            "confirmed_reads": len(confirmed),
        }
    return identification


# --- Re-probing (Box 3): extend contigs with more library reads ---------------

def _reprobe_reads_by_query(tabular_text: str) -> dict[str, set[str]]:
    """Group library read ids (sseqid) by the contig (qseqid) that recovered them.

    Expects ``6 qseqid sseqid evalue`` rows; keeps matches below ``REPROBE_EVALUE``
    (the net's gate), so re-probing is no more permissive than the first search.
    """
    by_query: dict[str, set[str]] = {}
    for line in tabular_text.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        query_id, read_id, evalue_str = fields[0], fields[1], fields[2]
        try:
            evalue = float(evalue_str)
        except ValueError:
            continue
        if evalue < REPROBE_EVALUE:
            by_query.setdefault(query_id, set()).add(read_id)
    return by_query


def _reprobe_hits(
    named_contigs: dict[str, str],
    patient_db_prefix: str,
    *,
    num_threads: int | str | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
) -> dict[str, set[str]]:
    """BLAST key contigs against the patient library; return reads per contig id."""
    if not named_contigs:
        return {}
    with tempfile.TemporaryDirectory(prefix="contig_reprobe_") as tmpdir:
        query_path = Path(tmpdir) / "key_contigs.fasta"
        query_path.write_text(reads_to_fasta(named_contigs), encoding="utf-8")
        command = [
            str(blast_exe("blastn")),
            "-task",
            "megablast",
            "-query",
            str(query_path),
            "-db",
            blast_safe_path(patient_db_prefix),
            "-evalue",
            str(REPROBE_EVALUE),
            "-max_target_seqs",
            REPROBE_MAX_TARGET_SEQS,
            "-outfmt",
            "6 qseqid sseqid evalue",
        ]
        if num_threads:
            command += ["-num_threads", str(num_threads)]
            # mt_mode 0 (let BLAST split by db), matching the eToL net search: the
            # reprobe is few long contigs vs a huge patient library, so the work is
            # db-bound. Forcing mt_mode 1 (split by query) flatlines it across
            # threads -- the regression blast_runner.EXACT_MATCH_MT_MODE documents.
            try:
                if int(num_threads) > 1:
                    command += ["-mt_mode", "0"]
            except (TypeError, ValueError):
                pass
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "Re-probing BLAST failed: "
            + ((completed.stderr or completed.stdout).strip() or "unknown error")
        )
    return _reprobe_reads_by_query(completed.stdout)


def _key_contigs(contigs: list[dict[str, Any]], top_contigs: int) -> list[dict[str, Any]]:
    """The taxon's most-supported contigs (most reads first), capped at top_contigs."""
    ranked = sorted(contigs, key=lambda c: int(c.get("num_reads", 0)), reverse=True)
    return [c for c in ranked[:top_contigs] if c.get("sequence")]


def reprobe_and_reassemble(
    contigs_by_species: dict[str, list[dict[str, Any]]],
    reads_by_taxon: dict[str, list[str]],
    all_reads: dict[str, str],
    *,
    patient_db_prefix: str,
    source_fasta_path: str,
    assembler: Any,
    top_contigs: int = REPROBE_TOP_CONTIGS,
    num_threads: int | str | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    human_db_prefix: str | None = None,
    human_bitscore_threshold: float = HUMAN_BITSCORE_THRESHOLD,
    assembly_pool: Any = None,
) -> dict[str, int]:
    """One round of contig re-probing (Hu, Haas & Lathe 2022, Box 3).

    Uses each taxon's top contigs as probes against the patient library, pulls the
    new reads they recover, optionally human-filters them (the paper filters
    re-probe matches too), then re-assembles each taxon from its original + new
    reads. Mutates ``contigs_by_species`` (replacing a taxon's contigs when
    re-assembly yields any), ``reads_by_taxon`` and ``all_reads`` (so the later
    confirmed-abundance step counts against the expanded read set) in place.

    ``assembly_pool`` is an optional ``Executor``: when given, the per-taxon CAP3
    re-assemblies fan out across it (each runs in its own temp dir, so they are
    independent); left ``None`` they run serially. The shared-dict mutations are
    always applied serially afterward.

    Returns ``{new_reads, reprobed_taxa, human_removed}`` for reporting.
    """
    # 1. One batched BLAST of every taxon's key contigs against the library.
    named: dict[str, str] = {}
    origin: dict[str, str] = {}
    for taxon, contigs in contigs_by_species.items():
        for contig in _key_contigs(contigs, top_contigs):
            synthetic_id = f"k{len(named)}"
            named[synthetic_id] = str(contig.get("sequence", ""))
            origin[synthetic_id] = taxon
    reads_by_query = _reprobe_hits(
        named, patient_db_prefix, num_threads=num_threads, timeout_seconds=timeout_seconds
    )

    # 2. Per taxon, the reads not already assigned to it (the genuinely new ones).
    existing_by_taxon = {taxon: set(ids) for taxon, ids in reads_by_taxon.items()}
    new_ids_by_taxon: dict[str, set[str]] = {}
    all_new_ids: set[str] = set()
    for synthetic_id, read_ids in reads_by_query.items():
        taxon = origin[synthetic_id]
        novel = read_ids - existing_by_taxon.get(taxon, set())
        if novel:
            new_ids_by_taxon.setdefault(taxon, set()).update(novel)
            all_new_ids.update(novel)
    if not all_new_ids:
        return {"new_reads": 0, "reprobed_taxa": 0, "human_removed": 0}

    # 3. Recover the new reads' sequences in one pass.
    new_reads, _method = extract_reads(patient_db_prefix, source_fasta_path, sorted(all_new_ids))

    # 4. Optional human filter on the new reads (paper filters re-probe matches).
    human_removed = 0
    if human_db_prefix and new_reads:
        human_ids = find_human_read_ids(
            new_reads,
            human_db_prefix,
            bitscore_threshold=human_bitscore_threshold,
            num_threads=num_threads,
            timeout_seconds=timeout_seconds,
        )
        if human_ids:
            human_removed = len(human_ids)
            new_reads = {rid: seq for rid, seq in new_reads.items() if rid not in human_ids}
    all_reads.update(new_reads)

    # 5. Re-assemble each re-probed taxon from its original + surviving new reads.
    # The CAP3 calls are independent (own temp dir each), so fan them out over the
    # shared pool when one is supplied; the mutations below stay on this thread.
    def _reassemble(taxon: str):
        novel_ids = new_ids_by_taxon[taxon]
        surviving = {rid for rid in novel_ids if rid in new_reads}
        if not surviving:
            return None
        combined_ids = list(existing_by_taxon.get(taxon, set())) + sorted(surviving)
        reads = {rid: all_reads[rid] for rid in combined_ids if rid in all_reads}
        return taxon, combined_ids, len(surviving), assembler.assemble(reads)

    taxa = list(new_ids_by_taxon)
    if assembly_pool is not None:
        results = list(assembly_pool.map(_reassemble, taxa))
    else:
        results = [_reassemble(taxon) for taxon in taxa]

    added_total = 0
    reprobed_taxa = 0
    for result in results:
        if not result:
            continue
        taxon, combined_ids, surviving_count, new_contigs = result
        if new_contigs:
            contigs_by_species[taxon] = [contig.to_dict() for contig in new_contigs]
            reads_by_taxon[taxon] = combined_ids
            added_total += surviving_count
            reprobed_taxa += 1
    return {
        "new_reads": added_total,
        "reprobed_taxa": reprobed_taxa,
        "human_removed": human_removed,
    }
