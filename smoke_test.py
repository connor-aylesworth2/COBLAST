"""Minimal backend smoke test for COBLAST+.

The launcher uses this to prove that BLAST+, FASTA validation, and parsing work
before opening the web interface.
"""

from pathlib import Path
import subprocess
import tempfile

from apoe_summary import build_apoe_probe_summary
from blast_runner import run_blast, validate_fasta_input
from config import blast_exe


DEMO_ROOT = Path(tempfile.gettempdir()) / "blast_flask_demo"
NT_SAMPLE_FASTA = DEMO_ROOT / "toy_nt.fasta"
PROTEIN_SAMPLE_FASTA = DEMO_ROOT / "toy_protein.fasta"
DB_PREFIX = DEMO_ROOT / "db" / "toy_nt"
PROTEIN_DB_PREFIX = DEMO_ROOT / "db" / "toy_protein"
TOY_SEQUENCE = "ATGCGTACGTAGCTAGCTAGCTAGCTA" * 4
TOY_CODING_SEQUENCE = "ATGGCTATGGCTCCTCGTACTGAAATTAATTCTACTCGTATTAATGGT"
TOY_PROTEIN_SEQUENCE = "MAMAPRTEINSTRING"


def build_database(fasta_path: Path, db_prefix: Path, dbtype: str) -> None:
    """Create a tiny BLAST database from a generated FASTA file."""
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(blast_exe("makeblastdb")),
        "-in",
        str(fasta_path),
        "-dbtype",
        dbtype,
        "-out",
        str(db_prefix),
    ]
    subprocess.run(cmd, check=True)


def ensure_toy_db() -> None:
    """Write toy nucleotide/protein FASTA files and index them with makeblastdb."""
    DB_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    NT_SAMPLE_FASTA.write_text(
        f">toy_sequence_1\n{TOY_SEQUENCE}\n"
        f">toy_coding_sequence\n{TOY_CODING_SEQUENCE}\n"
        ">toy_sequence_2\nTTTTTTTTCCCCCCCCAAAAGGGG\n",
        encoding="utf-8",
    )
    PROTEIN_SAMPLE_FASTA.write_text(
        f">toy_protein_1\n{TOY_PROTEIN_SEQUENCE}\n",
        encoding="utf-8",
    )

    build_database(NT_SAMPLE_FASTA, DB_PREFIX, "nucl")
    build_database(PROTEIN_SAMPLE_FASTA, PROTEIN_DB_PREFIX, "prot")


def exercise_validation() -> None:
    """Check accepted nucleotide/protein input and one rejected bad nucleotide."""
    nucleotide = validate_fasta_input(f">nt\nAUGCRYSWKMBDHVN", "nucleotide")
    protein = validate_fasta_input(f">protein\n{TOY_PROTEIN_SEQUENCE}", "protein")
    print(f"validated_nucleotide={nucleotide.total_length}")
    print(f"validated_protein={protein.total_length}")

    try:
        validate_fasta_input(">bad_nt\nATGQ", "nucleotide")
    except ValueError as exc:
        print(f"invalid_nucleotide_rejected={exc}")


def exercise_apoe_summary() -> None:
    """Check APOE exact-probe counts and per-site T percentages."""
    summary = build_apoe_probe_summary(
        [
            {
                "database_id": 1,
                "display_name": "SRX123456 APOE pilot",
                "db_prefix_path": r"C:\COBLAST_data\sra\SRX123456\reads",
                "hit_count": 5,
                "hits": [
                    {"qseqid": "AE4_E4=C"},
                    {"qseqid": "AE4_E23=T"},
                    {"qseqid": "AE4_E23=T"},
                    {"qseqid": "AE2_E34=C"},
                    {"qseqid": "AE2_E2=T"},
                ],
                "error": "",
            }
        ]
    )
    row = summary[0]
    assert row["sample"] == "SRX123456"
    assert row["ae4_c_hits"] == 1
    assert row["ae4_t_hits"] == 2
    assert row["total_exact_probe_hits"] == 5
    assert row["c_to_t_percent"] == "60.00"
    print(
        "apoe_summary="
        f"{row['sample']} c_to_t_percent={row['c_to_t_percent']}"
    )


def main() -> None:
    """Run all supported BLAST programs against the toy databases."""
    ensure_toy_db()
    exercise_validation()
    exercise_apoe_summary()

    searches = [
        ("blastn", f">query\n{TOY_SEQUENCE}", DB_PREFIX),
        ("blastp", f">query\n{TOY_PROTEIN_SEQUENCE}", PROTEIN_DB_PREFIX),
        ("blastx", f">query\n{TOY_CODING_SEQUENCE}", PROTEIN_DB_PREFIX),
        ("tblastn", f">query\n{TOY_PROTEIN_SEQUENCE}", DB_PREFIX),
    ]
    for program, sequence, database in searches:
        result = run_blast(
            sequence=sequence,
            database=database,
            program=program,
            output_format="tabular",
        )
        print(f"program={program}")
        print(f"returncode={result.returncode}")
        print(f"runtime_seconds={result.runtime_seconds:.3f}")
        print(f"hits={len(result.hits)}")
        for hit in result.hits[:3]:
            print(hit)
        if result.stderr:
            print("stderr:")
            print(result.stderr)


if __name__ == "__main__":
    main()
