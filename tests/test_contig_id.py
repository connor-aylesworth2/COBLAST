"""Tests for contig species identification + confirmed abundance.

The BLAST subprocess calls are never invoked: the pure tabular parsers are
tested directly, and ``identify_contigs`` is tested with its two BLAST helpers
monkeypatched, so the orchestration (synthetic-id mapping, in-place annotation,
per-taxon summary) is exercised without a real blastn.
"""

import contig_id
from contig_id import (
    _best_homolog_per_query,
    _confirmed_reads_from_tabular,
    _representative_contig,
    identify_contigs,
)


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
