"""Tests for matched-read recovery used by the secondary human filter."""

from types import SimpleNamespace

import human_filter


def test_extract_reads_falls_back_for_ids_missing_from_blastdbcmd(monkeypatch):
    monkeypatch.setattr(
        human_filter,
        "extract_reads_via_blastdbcmd",
        lambda _prefix, _ids: {"read-1": "AAAA"},
    )
    monkeypatch.setattr(
        human_filter,
        "extract_reads_from_fasta",
        lambda _path, ids: {read_id: "CCCC" for read_id in ids},
    )

    reads, method = human_filter.extract_reads(
        "patient-db", "patient.fasta", ["read-1", "read-2"]
    )

    assert reads == {"read-1": "AAAA", "read-2": "CCCC"}
    assert method == "blastdbcmd+source_fasta"


def test_extract_reads_reports_none_when_no_recovery_method_succeeds(monkeypatch):
    monkeypatch.setattr(
        human_filter, "extract_reads_via_blastdbcmd", lambda _prefix, _ids: None
    )
    monkeypatch.setattr(
        human_filter, "extract_reads_from_fasta", lambda _path, _ids: {}
    )

    reads, method = human_filter.extract_reads(
        "patient-db", "patient.fasta", ["read-1"]
    )

    assert reads == {}
    assert method == "none"


def test_find_human_read_ids_requires_full_query_coverage(monkeypatch):
    captured_command = []

    def run(command, **_kwargs):
        captured_command.extend(command)
        return SimpleNamespace(
            returncode=0,
            stdout="full-read\t100\npartial-read\t99.5\ninvalid\tnot-a-number\n",
            stderr="",
        )

    monkeypatch.setattr(human_filter, "blast_exe", lambda name: name)
    monkeypatch.setattr(human_filter.subprocess, "run", run)

    human_ids = human_filter.find_human_read_ids(
        {"full-read": "ACGT", "partial-read": "TGCA"}, "human-db"
    )

    assert human_ids == {"full-read"}
    qcov_index = captured_command.index("-qcov_hsp_perc")
    assert captured_command[qcov_index + 1] == "100"
    outfmt_index = captured_command.index("-outfmt")
    assert captured_command[outfmt_index + 1] == "6 qseqid qcovhsp"
