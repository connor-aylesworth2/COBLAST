from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import subprocess
import tempfile
from collections.abc import Iterable
from typing import Any

from config import blast_exe

try:
    from Bio import SearchIO
except ImportError as exc:  # pragma: no cover - exercised only when dependency is absent
    SearchIO = None
    BIOPYTHON_IMPORT_ERROR = exc
else:
    BIOPYTHON_IMPORT_ERROR = None


OUTFMT6_FIELDS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
ALLOWED_BLASTN_TASKS = {"blastn", "blastn-short", "dc-megablast", "megablast"}
BLAST_OUTPUT_FORMATS = {
    "tabular": "6 " + " ".join(OUTFMT6_FIELDS),
    "xml": "5",
}


@dataclass(frozen=True)
class BlastResult:
    returncode: int
    hits: list[dict[str, str]]
    stdout: str
    stderr: str
    command: list[str]
    output_format: str


def validate_fasta(sequence: str) -> str:
    cleaned = sequence.strip()
    if not cleaned:
        raise ValueError("Enter a FASTA sequence.")
    if not cleaned.startswith(">"):
        cleaned = ">query\n" + cleaned
    return cleaned + "\n"


def require_searchio() -> None:
    if SearchIO is None:
        raise RuntimeError(
            "Biopython is required for BLAST result parsing. "
            "Install dependencies with: python -m pip install -r requirements.txt"
        ) from BIOPYTHON_IMPORT_ERROR


def format_float(value: Any, decimals: int) -> str:
    if value is None:
        return ""
    return f"{float(value):.{decimals}f}"


def format_evalue(value: Any) -> str:
    if value is None:
        return ""
    number = float(value)
    if number == 0:
        return "0.0"
    return f"{number:.2e}"


def percent_identity(hsp: Any) -> float | None:
    if getattr(hsp, "ident_pct", None) is not None:
        return float(hsp.ident_pct)

    ident_num = getattr(hsp, "ident_num", None)
    aln_span = getattr(hsp, "aln_span", None)
    if ident_num is None or not aln_span:
        return None
    return (float(ident_num) / float(aln_span)) * 100


def searchio_results_to_hits(qresults: Iterable[Any]) -> list[dict[str, str]]:
    hits: list[dict[str, str]] = []
    for qresult in qresults:
        for hit in qresult:
            for hsp in hit:
                hits.append(
                    {
                        "qseqid": getattr(hsp, "query_id", None) or qresult.id,
                        "sseqid": getattr(hsp, "hit_id", None) or hit.id,
                        "pident": format_float(percent_identity(hsp), 3),
                        "length": str(getattr(hsp, "aln_span", "")),
                        "evalue": format_evalue(getattr(hsp, "evalue", None)),
                        "bitscore": format_float(getattr(hsp, "bitscore", None), 1),
                    }
                )
    return hits


def parse_blast_tabular(stdout: str) -> list[dict[str, str]]:
    require_searchio()
    if not stdout.strip():
        return []
    qresults = SearchIO.parse(StringIO(stdout), "blast-tab", fields=OUTFMT6_FIELDS)
    return searchio_results_to_hits(qresults)


def parse_blast_xml(stdout: str) -> list[dict[str, str]]:
    require_searchio()
    if not stdout.strip():
        return []
    qresults = SearchIO.parse(StringIO(stdout), "blast-xml")
    return searchio_results_to_hits(qresults)


def parse_blast_output(stdout: str, output_format: str = "tabular") -> list[dict[str, str]]:
    if output_format == "tabular":
        return parse_blast_tabular(stdout)
    if output_format == "xml":
        return parse_blast_xml(stdout)
    raise ValueError(f"Unsupported BLAST output format: {output_format}")


def run_blastn(
    sequence: str,
    database: str | Path,
    timeout_seconds: int = 60,
    task: str = "blastn-short",
    output_format: str = "tabular",
) -> BlastResult:
    fasta = validate_fasta(sequence)
    if task not in ALLOWED_BLASTN_TASKS:
        raise ValueError(f"Unsupported blastn task: {task}")
    if output_format not in BLAST_OUTPUT_FORMATS:
        raise ValueError(f"Unsupported BLAST output format: {output_format}")
    db_path = str(database)

    with tempfile.TemporaryDirectory(prefix="blast_flask_") as tmpdir:
        query_path = Path(tmpdir) / "query.fasta"
        query_path.write_text(fasta, encoding="utf-8")

        cmd = [
            str(blast_exe("blastn")),
            "-query",
            str(query_path),
            "-db",
            db_path,
            "-task",
            task,
            "-outfmt",
            BLAST_OUTPUT_FORMATS[output_format],
        ]

        completed = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )

    return BlastResult(
        returncode=completed.returncode,
        hits=parse_blast_output(completed.stdout, output_format)
        if completed.returncode == 0
        else [],
        stdout=completed.stdout,
        stderr=completed.stderr,
        command=cmd,
        output_format=output_format,
    )
