"""Unit tests for the single-walk SRA project helpers (no filesystem needed)."""

from pathlib import Path

from sra_workflow import find_blast_prefixes, find_fasta_files, find_sra_files


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
