"""Tests for contig species identification + confirmed abundance.

The BLAST subprocess calls are never invoked: the pure tabular parsers are
tested directly, and ``identify_contigs`` is tested with its two BLAST helpers
monkeypatched, so the orchestration (synthetic-id mapping, in-place annotation,
per-taxon summary) is exercised without a real blastn.
"""

from assembler import Contig

import contig_id
from contig_id import (
    _best_homolog_per_query,
    _confirmed_reads_from_tabular,
    _representative_contig,
    _reprobe_reads_by_query,
    identify_contigs,
    reprobe_and_reassemble,
    species_from_homolog,
)


# --- homolog -> species parsing ----------------------------------------------

def test_species_from_homolog_takes_last_lineage_rank():
    # Real reference-title format: accession in the first field, species last.
    full = ("CP031450.6146683.6148220 Bacteria;Pseudomonadota;Gammaproteobacteria;"
            "Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas fluorescens")
    assert species_from_homolog(full) == "Pseudomonas fluorescens"
    # No lineage (bare description) is returned as-is; trailing ';' tolerated.
    assert species_from_homolog("Thermotoga maritima") == "Thermotoga maritima"
    assert species_from_homolog("A;B;Hoeflea olei;") == "Hoeflea olei"
    assert species_from_homolog("") == ""


# --- tabular parsers ----------------------------------------------------------

def test_best_homolog_keeps_highest_bitscore_and_full_title():
    # qseqid pident length bitscore stitle  (stitle has spaces/semicolons)
    text = (
        "c0\t99.5\t450\t800.0\tBacteria;Thermotogae;Thermotoga maritima\n"
        "c0\t88.0\t300\t400.0\tBacteria;something else\n"
        "c1\t97.0\t480\t905.5\tArchaea;Methanocaldococcus jannaschii\n"
    )
    best = _best_homolog_per_query(text)
    assert best["c0"]["homolog"] == "Bacteria;Thermotogae;Thermotoga maritima"
    assert best["c0"]["pident"] == 99.5  # the 800-bit row wins over the 400-bit one
    assert best["c1"]["homolog"] == "Archaea;Methanocaldococcus jannaschii"


def test_best_homolog_skips_malformed_rows():
    text = "c0\tnot_a_number\t1\t2\ttitle\nc1\t99.0\t100\t150.0\tGood hit\n"
    best = _best_homolog_per_query(text)
    assert "c0" not in best
    assert best["c1"]["homolog"] == "Good hit"


def test_confirmed_reads_counts_distinct_reads_per_contig_at_cutoff():
    # qseqid(read) sseqid(contig) pident
    text = (
        "r1\tContig1\t100.0\n"
        "r2\tContig1\t99.0\n"
        "r2\tContig1\t99.4\n"   # same read, still one
        "r3\tContig1\t80.0\n"   # below 99 -> not confirmed
        "r4\tContig2\t99.0\n"
    )
    per_contig, confirmed = _confirmed_reads_from_tabular(text, identity_pct=99.0)
    assert per_contig == {"Contig1": 2, "Contig2": 1}  # r1,r2 ; r4
    assert confirmed == {"r1", "r2", "r4"}


def test_representative_contig_prefers_named_then_most_reads():
    contigs = [
        {"id": "Contig1", "num_reads": 2, "length": 300, "closest_homolog": "X"},
        {"id": "Contig2", "num_reads": 9, "length": 300},  # more reads, but unnamed
    ]
    assert _representative_contig(contigs)["id"] == "Contig1"


# --- identify_contigs orchestration -------------------------------------------

def test_identify_contigs_annotates_and_summarizes(monkeypatch):
    contigs_by_species = {
        "B0_Tmaritima_16S": [
            {"id": "Contig1", "sequence": "AAA", "num_reads": 5, "length": 3},
            {"id": "Contig2", "sequence": "CCC", "num_reads": 2, "length": 3},
        ],
        "A_Mjannaschii_16S": [
            {"id": "Contig1", "sequence": "GGG", "num_reads": 9, "length": 3},
        ],
    }
    reads_by_taxon = {"B0_Tmaritima_16S": ["r1", "r2"], "A_Mjannaschii_16S": ["r9"]}
    all_reads = {"r1": "AAA", "r2": "AAA", "r9": "GGG"}

    seq_homolog = {
        "AAA": "Thermotoga maritima",
        "CCC": "Some divergent clone",
        "GGG": "Methanocaldococcus jannaschii",
    }

    def fake_name(named, reference_db_prefix, **_kwargs):
        return {
            sid: {"homolog": seq_homolog[seq], "pident": 99.0, "length": 3, "bitscore": 99.0}
            for sid, seq in named.items()
        }

    def fake_confirm(contigs, reads, *, identity_pct=99.0, **_kwargs):
        if not contigs or not reads:
            return {}, set()
        return {contigs[0]["id"]: len(reads)}, set(reads)

    monkeypatch.setattr(contig_id, "name_contigs", fake_name)
    monkeypatch.setattr(contig_id, "confirm_contig_reads", fake_confirm)

    identification = identify_contigs(
        contigs_by_species,
        reads_by_taxon,
        all_reads,
        reference_db_prefix="ref",
    )

    # Per-contig annotation written in place, mapped back to the right contig.
    t = contigs_by_species["B0_Tmaritima_16S"]
    assert t[0]["closest_homolog"] == "Thermotoga maritima"
    assert t[0]["confirmed_reads"] == 2          # fake_confirm gave all reads to Contig1
    assert t[1]["closest_homolog"] == "Some divergent clone"
    assert t[1]["confirmed_reads"] == 0          # not in per_contig -> 0

    # Per-taxon summary: representative is the most-supported contig (Contig1, 5 reads),
    # confirmed total is the distinct read count.
    assert identification["B0_Tmaritima_16S"]["closest_homolog"] == "Thermotoga maritima"
    assert identification["B0_Tmaritima_16S"]["confirmed_reads"] == 2
    assert identification["A_Mjannaschii_16S"]["closest_homolog"] == "Methanocaldococcus jannaschii"
    assert identification["A_Mjannaschii_16S"]["confirmed_reads"] == 1


def test_identify_contigs_drops_human_contigs(monkeypatch):
    contigs_by_species = {
        "H3_Human_18S": [
            {"id": "Contig1", "sequence": "AAA", "num_reads": 9, "length": 3},
        ],
        "B0_T": [
            {"id": "Contig1", "sequence": "GGG", "num_reads": 5, "length": 3},
        ],
    }
    reads_by_taxon = {"H3_Human_18S": ["r1"], "B0_T": ["r9"]}
    all_reads = {"r1": "AAA", "r9": "GGG"}

    seq_homolog = {
        "AAA": "ACC.1.2 Eukaryota;Metazoa;Chordata;Mammalia;Homo sapiens (human)",
        "GGG": "ACC.3.4 Bacteria;Thermotogae;Thermotoga maritima",
    }

    def fake_name(named, reference_db_prefix, **_kwargs):
        return {
            sid: {"homolog": seq_homolog[seq], "pident": 99.0, "length": 3, "bitscore": 99.0}
            for sid, seq in named.items()
        }

    monkeypatch.setattr(contig_id, "name_contigs", fake_name)
    monkeypatch.setattr(contig_id, "confirm_contig_reads", lambda *a, **k: ({}, set()))

    identification = identify_contigs(
        contigs_by_species, reads_by_taxon, all_reads, reference_db_prefix="ref"
    )

    # Human contig removed from the taxon's list and its identification blanked.
    assert contigs_by_species["H3_Human_18S"] == []
    assert identification["H3_Human_18S"]["closest_species"] == ""
    # Non-human taxon keeps its contig and gets a parsed species name.
    assert [c["id"] for c in contigs_by_species["B0_T"]] == ["Contig1"]
    assert identification["B0_T"]["closest_species"] == "Thermotoga maritima"


def test_identify_contigs_skips_sequenceless_contigs(monkeypatch):
    contigs_by_species = {"T": [{"id": "Contig1", "sequence": "", "num_reads": 0, "length": 0}]}
    captured = {}

    def fake_name(named, reference_db_prefix, **_kwargs):
        captured["named"] = named
        return {}

    monkeypatch.setattr(contig_id, "name_contigs", fake_name)
    monkeypatch.setattr(contig_id, "confirm_contig_reads", lambda *a, **k: ({}, set()))

    identify_contigs(contigs_by_species, {"T": []}, {}, reference_db_prefix="ref")
    assert captured["named"] == {}  # empty-sequence contig is not sent to BLAST


# --- re-probing ---------------------------------------------------------------

def test_reprobe_reads_by_query_groups_and_evalue_gates():
    # qseqid(contig) sseqid(read) evalue
    text = (
        "k0\tr1\t1e-30\n"
        "k0\tr2\t0.0\n"
        "k0\tr3\t0.5\n"     # E >= 0.01 -> dropped
        "k1\tr9\t1e-5\n"
        "k1\tbad\tnope\n"   # unparseable evalue -> skipped
    )
    grouped = _reprobe_reads_by_query(text)
    assert grouped == {"k0": {"r1", "r2"}, "k1": {"r9"}}


def _fake_assembler(num_reads_label=None):
    """Assembler stand-in: collapses any reads into one labelled contig."""

    class _Fake:
        name = "fake"

        def assemble(self, reads):
            if not reads:
                return []
            return [Contig(id="Reassembled1", sequence="ACGTACGT", num_reads=len(reads))]

    return _Fake()


def test_reprobe_and_reassemble_pulls_reads_and_rebuilds(monkeypatch):
    contigs_by_species = {
        "B0_T": [
            {"id": "Contig1", "sequence": "AAA", "num_reads": 5, "length": 3},
            {"id": "Contig2", "sequence": "CCC", "num_reads": 2, "length": 3},
        ],
        "A_M": [{"id": "Contig1", "sequence": "GGG", "num_reads": 9, "length": 3}],
    }
    reads_by_taxon = {"B0_T": ["r1", "r2"], "A_M": ["r9"]}
    all_reads = {"r1": "AAA", "r2": "AAA", "r9": "GGG"}

    # Re-probe finds new reads keyed by the contig sequence (so it maps back to
    # whatever synthetic id the orchestrator assigned).
    def fake_reprobe_hits(named, patient_db_prefix, **_kwargs):
        seq_to_reads = {"AAA": {"rNEW1", "rNEW2"}, "GGG": {"r9b"}}
        return {sid: set(seq_to_reads[seq]) for sid, seq in named.items() if seq in seq_to_reads}

    def fake_extract(db_prefix, source_fasta, read_ids):
        seqs = {"rNEW1": "AAAACGT", "rNEW2": "AAAACGT", "r9b": "GGGACGT"}
        return {rid: seqs[rid] for rid in read_ids if rid in seqs}, "blastdbcmd"

    monkeypatch.setattr(contig_id, "_reprobe_hits", fake_reprobe_hits)
    monkeypatch.setattr(contig_id, "extract_reads", fake_extract)

    stats = reprobe_and_reassemble(
        contigs_by_species,
        reads_by_taxon,
        all_reads,
        patient_db_prefix="patient",
        source_fasta_path="",
        assembler=_fake_assembler(),
    )

    assert stats == {"new_reads": 3, "reprobed_taxa": 2, "human_removed": 0}
    # Contigs replaced by the re-assembly output.
    assert [c["id"] for c in contigs_by_species["B0_T"]] == ["Reassembled1"]
    # The taxon's read set now includes the recovered reads (order-independent).
    assert set(reads_by_taxon["B0_T"]) == {"r1", "r2", "rNEW1", "rNEW2"}
    assert "rNEW1" in all_reads  # new sequences merged in for confirmed-abundance


def test_reprobe_reassembles_in_parallel_when_pool_given(monkeypatch):
    from concurrent.futures import ThreadPoolExecutor

    contigs_by_species = {
        "B0_T": [{"id": "Contig1", "sequence": "AAA", "num_reads": 5, "length": 3}],
        "A_M": [{"id": "Contig1", "sequence": "GGG", "num_reads": 9, "length": 3}],
    }
    reads_by_taxon = {"B0_T": ["r1"], "A_M": ["r9"]}
    all_reads = {"r1": "AAA", "r9": "GGG"}

    def fake_reprobe_hits(named, patient_db_prefix, **_kwargs):
        seq_to_reads = {"AAA": {"rNEW1"}, "GGG": {"r9b"}}
        return {sid: set(seq_to_reads[seq]) for sid, seq in named.items() if seq in seq_to_reads}

    def fake_extract(db_prefix, source_fasta, read_ids):
        seqs = {"rNEW1": "AAAACGT", "r9b": "GGGACGT"}
        return {rid: seqs[rid] for rid in read_ids if rid in seqs}, "blastdbcmd"

    monkeypatch.setattr(contig_id, "_reprobe_hits", fake_reprobe_hits)
    monkeypatch.setattr(contig_id, "extract_reads", fake_extract)

    with ThreadPoolExecutor(max_workers=2) as pool:
        stats = reprobe_and_reassemble(
            contigs_by_species,
            reads_by_taxon,
            all_reads,
            patient_db_prefix="patient",
            source_fasta_path="",
            assembler=_fake_assembler(),
            assembly_pool=pool,
        )

    assert stats == {"new_reads": 2, "reprobed_taxa": 2, "human_removed": 0}
    assert [c["id"] for c in contigs_by_species["B0_T"]] == ["Reassembled1"]
    assert set(reads_by_taxon["A_M"]) == {"r9", "r9b"}


def test_reprobe_human_filters_new_reads(monkeypatch):
    contigs_by_species = {"B0_T": [{"id": "Contig1", "sequence": "AAA", "num_reads": 5, "length": 3}]}
    reads_by_taxon = {"B0_T": ["r1"]}
    all_reads = {"r1": "AAA"}

    monkeypatch.setattr(
        contig_id, "_reprobe_hits", lambda named, db, **k: {sid: {"rNEW1", "rNEW2"} for sid in named}
    )
    monkeypatch.setattr(
        contig_id,
        "extract_reads",
        lambda db, src, ids: ({rid: "AAAACGT" for rid in ids}, "blastdbcmd"),
    )
    # rNEW2 looks human and must be dropped before re-assembly.
    monkeypatch.setattr(contig_id, "find_human_read_ids", lambda reads, db, **k: {"rNEW2"})

    stats = reprobe_and_reassemble(
        contigs_by_species,
        reads_by_taxon,
        all_reads,
        patient_db_prefix="patient",
        source_fasta_path="",
        assembler=_fake_assembler(),
        human_db_prefix="human",
    )

    assert stats["human_removed"] == 1
    assert stats["new_reads"] == 1  # only rNEW1 survived
    assert set(reads_by_taxon["B0_T"]) == {"r1", "rNEW1"}
    assert "rNEW2" not in all_reads
