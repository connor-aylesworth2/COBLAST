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

from config import blast_exe
from database_registry import blast_safe_path


# Reads aligning to a contig at or above this percent identity are "confirmed"
# (the paper keeps "100% or near-100% identity"). 99% tolerates a single
# sequencing error over a short rRNA read.
DEFAULT_CONFIRM_IDENTITY_PCT = 99.0
# Species naming is a real homology search, so keep a modest e-value gate; the
# confirmation search is identity-gated in Python, so its e-value stays loose.
DEFAULT_NAME_EVALUE = "1e-5"
DEFAULT_TIMEOUT_SECONDS = 1800


def _to_fasta(records: dict[str, str]) -> str:
    """Render ``{id: sequence}`` as FASTA text."""
    return "".join(f">{rec_id}\n{sequence}\n" for rec_id, sequence in records.items())


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
        query_path.write_text(_to_fasta(named_contigs), encoding="utf-8")
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
        reads_path.write_text(_to_fasta(reads), encoding="utf-8")
        contigs_path.write_text(_to_fasta(subjects), encoding="utf-8")
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
        contig["homolog_pident"] = info["pident"]

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
            "homolog_pident": (representative or {}).get("homolog_pident"),
            "confirmed_reads": len(confirmed),
        }
    return identification
