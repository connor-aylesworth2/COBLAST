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


def test_find_human_read_ids_uses_bitscore_threshold(monkeypatch):
    captured_command = []

    def run(command, **_kwargs):
        captured_command.extend(command)
        return SimpleNamespace(
            returncode=0,
            # 200.5 > 150 -> human; 150 is not strictly > 150 -> kept; 80 -> kept.
            stdout="human-read\t200.5\nborderline\t150\nweak-read\t80.0\ninvalid\tnope\n",
            stderr="",
        )

    monkeypatch.setattr(human_filter, "blast_exe", lambda name: name)
    monkeypatch.setattr(human_filter.subprocess, "run", run)

    human_ids = human_filter.find_human_read_ids(
        {"human-read": "ACGT", "borderline": "TTTT", "weak-read": "TGCA"}, "human-db"
    )

    assert human_ids == {"human-read"}
    outfmt_index = captured_command.index("-outfmt")
    assert captured_command[outfmt_index + 1] == "6 qseqid bitscore"
    # The bitscore is the sole criterion: no coverage filter is applied.
    assert "-qcov_hsp_perc" not in captured_command


def test_find_human_read_ids_threshold_is_configurable(monkeypatch):
    def run(_command, **_kwargs):
        return SimpleNamespace(
            returncode=0,
            stdout="read-a\t130\nread-b\t90\n",
            stderr="",
        )

    monkeypatch.setattr(human_filter, "blast_exe", lambda name: name)
    monkeypatch.setattr(human_filter.subprocess, "run", run)

    human_ids = human_filter.find_human_read_ids(
        {"read-a": "ACGT", "read-b": "TGCA"}, "human-db", bitscore_threshold=100.0
    )

    assert human_ids == {"read-a"}  # 130 > 100, 90 < 100
