"""Tests for the pluggable CAP3 contig assembler and its binary resolver.

These never invoke a real CAP3 binary: the executable resolver and
``subprocess.run`` are monkeypatched, and the fake process writes CAP3's
filename-based output files into the working directory the assembler chose.
"""

from pathlib import Path
from types import SimpleNamespace

import pytest

import assembler
import config
from assembler import Cap3Assembler, Contig, parse_cap3_contigs, _ace_read_counts
from etol_summary import group_read_ids_by_taxon


SAMPLE_CONTIGS = """>Contig1
ACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGT
>Contig2
TTTTGGGGCCCCAAAA
"""

SAMPLE_ACE = """AS 2 16
CO Contig1 40 12 3 U
ACGTACGT
CO Contig2 16 3 1 U
TTTTGGGG
"""


# --- Contig / parsing ---------------------------------------------------------

def test_contig_to_dict_and_length():
    contig = Contig(id="Contig1", sequence="ACGTACGT", num_reads=5)
    assert contig.length == 8
    assert contig.to_dict() == {
        "id": "Contig1",
        "sequence": "ACGTACGT",
        "num_reads": 5,
        "length": 8,
    }


def test_ace_read_counts_parses_contig_lines():
    assert _ace_read_counts(SAMPLE_ACE) == {"Contig1": 12, "Contig2": 3}


def test_parse_cap3_contigs_joins_wrapped_sequence_and_applies_counts():
    contigs = parse_cap3_contigs(SAMPLE_CONTIGS, {"Contig1": 12})
    assert [contig.id for contig in contigs] == ["Contig1", "Contig2"]
    assert contigs[0].sequence == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    assert contigs[0].length == 40
    assert contigs[0].num_reads == 12
    assert contigs[1].num_reads == 0  # missing from counts -> defaults to 0


# --- Cap3Assembler ------------------------------------------------------------

def test_is_available_reflects_cap3_exe(monkeypatch):
    monkeypatch.setattr(assembler, "cap3_exe", lambda: Path("cap3"))
    assert Cap3Assembler().is_available() is True

    def _missing():
        raise FileNotFoundError("nope")

    monkeypatch.setattr(assembler, "cap3_exe", _missing)
    assert Cap3Assembler().is_available() is False


def test_assemble_skips_when_fewer_than_two_reads(monkeypatch):
    calls = []
    monkeypatch.setattr(assembler, "cap3_exe", lambda: Path("cap3"))
    monkeypatch.setattr(
        assembler.subprocess,
        "run",
        lambda *args, **kwargs: calls.append(args) or SimpleNamespace(returncode=0),
    )
    assert Cap3Assembler().assemble({"r1": "ACGTACGT"}) == []
    assert calls == []  # CAP3 is never invoked for a lone read


def _fake_run_writing(contigs_text=None, ace_text=None, returncode=0, stderr=""):
    """Build a subprocess.run stand-in that writes CAP3's output into cwd."""

    def run(command, cwd=None, **_kwargs):
        if contigs_text is not None:
            (Path(cwd) / "reads.fasta.cap.contigs").write_text(contigs_text, encoding="utf-8")
        if ace_text is not None:
            (Path(cwd) / "reads.fasta.cap.ace").write_text(ace_text, encoding="utf-8")
        return SimpleNamespace(returncode=returncode, stdout="", stderr=stderr)

    return run


def test_assemble_parses_cap3_output(monkeypatch):
    monkeypatch.setattr(assembler, "cap3_exe", lambda: Path("cap3"))
    monkeypatch.setattr(
        assembler.subprocess, "run", _fake_run_writing(SAMPLE_CONTIGS, SAMPLE_ACE)
    )
    contigs = Cap3Assembler().assemble(
        {"r1": "ACGTACGT", "r2": "ACGTACGT", "r3": "TTTTGGGG"}
    )
    assert [contig.id for contig in contigs] == ["Contig1", "Contig2"]
    assert contigs[0].num_reads == 12  # read support comes from the ACE file


def test_assemble_returns_empty_when_no_contigs_file(monkeypatch):
    monkeypatch.setattr(assembler, "cap3_exe", lambda: Path("cap3"))
    monkeypatch.setattr(assembler.subprocess, "run", _fake_run_writing(contigs_text=None))
    assert Cap3Assembler().assemble({"r1": "ACGTACGT", "r2": "TGCATGCA"}) == []


def test_assemble_raises_on_nonzero_returncode(monkeypatch):
    monkeypatch.setattr(assembler, "cap3_exe", lambda: Path("cap3"))
    monkeypatch.setattr(
        assembler.subprocess, "run", _fake_run_writing(returncode=1, stderr="boom")
    )
    with pytest.raises(RuntimeError, match="CAP3 assembly failed: boom"):
        Cap3Assembler().assemble({"r1": "ACGTACGT", "r2": "TGCATGCA"})


def test_assemble_passes_overlap_options(monkeypatch):
    captured = {}

    def run(command, cwd=None, **_kwargs):
        captured["command"] = command
        (Path(cwd) / "reads.fasta.cap.contigs").write_text(SAMPLE_CONTIGS, encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(assembler, "cap3_exe", lambda: Path("cap3"))
    monkeypatch.setattr(assembler.subprocess, "run", run)
    Cap3Assembler(overlap_length="50", overlap_identity_pct="95").assemble(
        {"r1": "ACGTACGT", "r2": "TGCATGCA"}
    )
    command = captured["command"]
    assert command[command.index("-o") + 1] == "50"
    assert command[command.index("-p") + 1] == "95"


# --- config.cap3_exe resolver -------------------------------------------------

def test_cap3_exe_prefers_env_dir(tmp_path, monkeypatch):
    exe = tmp_path / config.tool_name("cap3")
    exe.write_text("", encoding="utf-8")
    monkeypatch.setenv("CAP3_BIN", str(tmp_path))
    assert config.cap3_exe() == exe


def test_cap3_exe_falls_back_to_path(monkeypatch):
    monkeypatch.delenv("CAP3_BIN", raising=False)
    monkeypatch.setattr(config.shutil, "which", lambda _name: "/usr/local/bin/cap3")
    assert config.cap3_exe() == Path("/usr/local/bin/cap3")


def test_cap3_exe_raises_when_missing(monkeypatch):
    monkeypatch.delenv("CAP3_BIN", raising=False)
    monkeypatch.setattr(config.shutil, "which", lambda _name: None)
    with pytest.raises(FileNotFoundError, match="Could not find the CAP3 assembler"):
        config.cap3_exe()


# --- per-taxon read grouping --------------------------------------------------

def test_group_read_ids_by_taxon_dedups_and_preserves_order():
    hits = [
        {"qseqid": "B0_Tmaritima_16S_3", "sseqid": "read-1"},
        {"qseqid": "B0_Tmaritima_16S_7", "sseqid": "read-2"},
        {"qseqid": "B0_Tmaritima_16S_3", "sseqid": "read-1"},  # duplicate read id
        {"qseqid": "A_Mjannaschii_16S_1", "sseqid": "read-9"},
        {"qseqid": "B0_Tmaritima_16S_3", "sseqid": ""},         # no read id -> skipped
    ]
    grouped = group_read_ids_by_taxon(hits)
    assert grouped["B0_Tmaritima_16S"] == ["read-1", "read-2"]
    assert grouped["A_Mjannaschii_16S"] == ["read-9"]
    assert list(grouped) == ["B0_Tmaritima_16S", "A_Mjannaschii_16S"]
