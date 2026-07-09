"""Unit tests for the single-walk SRA project helpers (no filesystem needed)."""

from pathlib import Path

import sra_workflow

from sra_workflow import (
    SraFileSummary,
    build_fetch_script_lines,
    find_blast_prefixes,
    find_fasta_files,
    find_sra_files,
    source_fasta_for_blast_prefix,
)


def test_fetch_script_indexes_each_run_into_a_blastdb():
    # The whole point of the fetch flow: a fetched run must end blast-ready, or
    # there is nothing to register. Guard the prefetch -> FASTA -> makeblastdb chain.
    steps = build_fetch_script_lines(
        ["SRR1"],
        Path("/data/sra"),
        Path("/t/prefetch"),
        Path("/t/fastq-dump"),
        Path("/t/makeblastdb"),
    )
    headers = [header for header, _ in steps]
    cmds = [cmd for _, cmd in steps]
    assert len(steps) == 3
    assert "--progress" in cmds[0]  # live download bar so the user sees progress
    assert "--max-size u" in cmds[0]  # no 20G default cap silently truncating a run
    assert "--split-spot" in cmds[1]  # mates stay separate, not chimeric
    assert "makeblastdb" in cmds[2] and "-parse_seqids" in cmds[2]  # id index for eToL
    assert headers == [  # each step announces its position before it runs
        "[run 1/1] SRR1 - step 1/3: downloading .sra",
        "[run 1/1] SRR1 - step 2/3: converting to FASTA",
        "[run 1/1] SRR1 - step 3/3: building BLAST database",
    ]


def test_find_fasta_files_filters_and_sorts():
    files = [Path("d/b.fasta"), Path("d/a.fa"), Path("d/x.txt"), Path("d/c.fna")]
    assert find_fasta_files(files) == sorted(
        [Path("d/b.fasta"), Path("d/a.fa"), Path("d/c.fna")]
    )


def test_find_sra_files_filters_case_insensitively():
    files = [Path("d/r1.sra"), Path("d/r2.fasta"), Path("d/r3.SRA")]
    assert find_sra_files(files) == sorted([Path("d/r1.sra"), Path("d/r3.SRA")])


def test_find_blast_prefixes_dedupes_and_skips_volume_files():
    nin = Path("db/patient.nin")
    nal = Path("db/alias.nal")
    volume = Path("db/big.00.nin")  # a numbered DB volume -> not a prefix
    prefixes = find_blast_prefixes([nin, Path("db/patient.nhr"), nal, volume])

    assert str(nin.with_suffix("")) in prefixes
    assert str(nal.with_suffix("")) in prefixes
    assert str(volume.with_suffix("")) not in prefixes
    # patient appears once even though several of its files were present.
    assert prefixes.count(str(nin.with_suffix(""))) == 1


def test_source_fasta_for_blast_prefix_prefers_same_stem():
    fasta_files = [
        SraFileSummary("project/fasta/SRR1_1.fasta", 10, "10 B"),
        SraFileSummary("project/fasta/SRR1.fasta", 20, "20 B"),
        SraFileSummary("project/fasta/SRR1_2.fasta", 10, "10 B"),
    ]

    source = source_fasta_for_blast_prefix("project/blastdb/SRR1", fasta_files)

    assert source == "project/fasta/SRR1.fasta"


def test_source_fasta_for_blast_prefix_does_not_guess_between_multiple_files():
    fasta_files = [
        SraFileSummary("project/fasta/reads_1.fasta", 10, "10 B"),
        SraFileSummary("project/fasta/reads_2.fasta", 10, "10 B"),
    ]

    assert source_fasta_for_blast_prefix("project/blastdb/SRR1", fasta_files) == ""


def test_register_sra_database_preserves_an_existing_source_path(monkeypatch):
    existing = type("Existing", (), {"source_fasta_path": "reads/SRR1.fasta"})()
    monkeypatch.setattr(sra_workflow, "get_database_by_prefix", lambda _prefix: existing)
    # Keep the test BLAST-free: the display name now comes from the DB's own title.
    monkeypatch.setattr(
        sra_workflow, "verify_database_prefix", lambda _prefix: {"database_title": ""}
    )
    captured = {}

    def register(**fields):
        captured.update(fields)
        return fields

    monkeypatch.setattr(sra_workflow, "register_existing_database", register)

    sra_workflow.register_sra_blast_database(
        accession="SRR1", db_prefix_path="blastdb/SRR1"
    )

    assert captured["source_fasta_path"] == "reads/SRR1.fasta"
