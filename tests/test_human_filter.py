"""Tests for matched-read recovery used by the secondary human filter."""

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
