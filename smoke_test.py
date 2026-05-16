from pathlib import Path
import subprocess
import tempfile

from blast_runner import run_blastn
from config import blast_exe


DEMO_ROOT = Path(tempfile.gettempdir()) / "blast_flask_demo"
SAMPLE_FASTA = DEMO_ROOT / "toy_nt.fasta"
DB_PREFIX = DEMO_ROOT / "db" / "toy_nt"
TOY_SEQUENCE = "ATGCGTACGTAGCTAGCTAGCTAGCTA" * 4


def ensure_toy_db() -> None:
    DB_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    SAMPLE_FASTA.write_text(
        f">toy_sequence_1\n{TOY_SEQUENCE}\n"
        ">toy_sequence_2\nTTTTTTTTCCCCCCCCAAAAGGGG\n",
        encoding="utf-8",
    )

    cmd = [
        str(blast_exe("makeblastdb")),
        "-in",
        str(SAMPLE_FASTA),
        "-dbtype",
        "nucl",
        "-out",
        str(DB_PREFIX),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    ensure_toy_db()
    for output_format in ("tabular", "xml"):
        result = run_blastn(
            sequence=f">query\n{TOY_SEQUENCE}",
            database=DB_PREFIX,
            output_format=output_format,
        )
        print(f"output_format={output_format}")
        print(f"returncode={result.returncode}")
        print(f"hits={len(result.hits)}")
        for hit in result.hits:
            print(hit)
        if result.stderr:
            print("stderr:")
            print(result.stderr)


if __name__ == "__main__":
    main()
