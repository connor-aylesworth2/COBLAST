"""Small JSON result store used for CSV/TSV downloads.

Search results are rendered immediately, but saving a copy lets the results page
serve export links without rerunning BLAST.
"""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import csv
import io
import json
from pathlib import Path
from uuid import UUID, uuid4

from apoe_summary import APOE_SUMMARY_EXPORT_COLUMNS, build_apoe_probe_summary
from etol_summary import (
    ETOL_PROBE_EXPORT_COLUMNS,
    ETOL_SPECIES_EXPORT_COLUMNS,
    etol_preset_records,
    etol_probe_count_rows,
    etol_species_count_rows,
)
from blast_runner import BlastResult, wrap_sequence
from config import runtime_data_dir


RESULTS_DIR = runtime_data_dir() / "results"
BATCH_RESULTS_DIR = runtime_data_dir() / "batch_results"

RESULT_COLUMNS = [
    # The first value is the internal hit key; the second is the export header.
    ("qseqid", "Query"),
    ("sseqid", "Subject"),
    ("stitle", "Subject title"),
    ("pident", "Percent identity"),
    ("length", "Alignment length"),
    ("qcovs", "Query coverage"),
    ("evalue", "E-value"),
    ("bitscore", "Bit score"),
]


def result_path(run_id: str) -> Path:
    """Resolve a saved result path after validating the UUID-like run id."""
    try:
        safe_id = str(UUID(run_id))
    except ValueError as exc:
        raise FileNotFoundError("Invalid result identifier.") from exc
    return RESULTS_DIR / f"{safe_id}.json"


def batch_result_path(batch_id: str) -> Path:
    """Resolve a saved batch-result path after validating the UUID-like id."""
    try:
        safe_id = str(UUID(batch_id))
    except ValueError as exc:
        raise FileNotFoundError("Invalid batch result identifier.") from exc
    return BATCH_RESULTS_DIR / f"{safe_id}.json"


def save_result(result: BlastResult) -> str:
    """Serialize one BlastResult to JSON and return its generated run id."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    run_id = str(uuid4())
    payload = asdict(result)
    payload["run_id"] = run_id
    payload["saved_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

    result_path(run_id).write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    return run_id


def save_batch_result(payload: dict) -> str:
    """Serialize one batch BLAST payload and return its generated batch id."""
    BATCH_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    batch_id = str(uuid4())
    payload = dict(payload)
    payload["batch_id"] = batch_id
    payload["saved_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

    batch_result_path(batch_id).write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    return batch_id


def load_result(run_id: str) -> dict:
    """Load a previously saved result payload."""
    path = result_path(run_id)
    if not path.exists():
        raise FileNotFoundError("Result not found.")
    return json.loads(path.read_text(encoding="utf-8"))


def load_batch_result(batch_id: str) -> dict:
    """Load a previously saved batch result payload."""
    path = batch_result_path(batch_id)
    if not path.exists():
        raise FileNotFoundError("Batch result not found.")
    return json.loads(path.read_text(encoding="utf-8"))


def result_rows_as_delimited(result_data: dict, delimiter: str) -> str:
    """Render saved hits as CSV or TSV text."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in RESULT_COLUMNS])
    for hit in result_data.get("hits", []):
        writer.writerow([hit.get(key, "") for key, _ in RESULT_COLUMNS])
    return buffer.getvalue()


def batch_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render batch hits and per-database errors as CSV or TSV text."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow(
        [
            "Database",
            "Program",
            "Query",
            "Subject",
            "Subject title",
            "Percent identity",
            "Alignment length",
            "Query coverage",
            "E-value",
            "Bit score",
            "Runtime seconds",
            "Return code",
            "Error",
        ]
    )
    for database_result in batch_data.get("database_results", []):
        hits = database_result.get("hits", [])
        if hits:
            for hit in hits:
                writer.writerow(
                    [
                        database_result.get("display_name", ""),
                        batch_data.get("program", ""),
                        hit.get("qseqid", ""),
                        hit.get("sseqid", ""),
                        hit.get("stitle", ""),
                        hit.get("pident", ""),
                        hit.get("length", ""),
                        hit.get("qcovs", ""),
                        hit.get("evalue", ""),
                        hit.get("bitscore", ""),
                        database_result.get("runtime_seconds", ""),
                        database_result.get("returncode", ""),
                        database_result.get("error", ""),
                    ]
                )
        else:
            writer.writerow(
                [
                    database_result.get("display_name", ""),
                    batch_data.get("program", ""),
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    database_result.get("runtime_seconds", ""),
                    database_result.get("returncode", ""),
                    database_result.get("error", ""),
                ]
            )
    return buffer.getvalue()


def apoe_summary_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render APOE per-sample probe summaries as CSV or TSV text."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in APOE_SUMMARY_EXPORT_COLUMNS])

    summary_rows = batch_data.get("apoe_probe_summary")
    if summary_rows is None:
        summary_rows = build_apoe_probe_summary(batch_data.get("database_results", []))

    for row in summary_rows:
        writer.writerow([apoe_summary_export_value(row, key) for key, _ in APOE_SUMMARY_EXPORT_COLUMNS])
    return buffer.getvalue()


def _etol_records_for_batch(batch_data: dict) -> tuple:
    """Resolve the probe panel used by a saved eToL batch (defaults to full)."""
    return etol_preset_records(batch_data.get("etol_preset_key") or "etol_full")


def etol_summary_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render eToL per-species exact-hit counts (every taxon) as CSV or TSV text."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in ETOL_SPECIES_EXPORT_COLUMNS])
    records = _etol_records_for_batch(batch_data)
    for row in etol_species_count_rows(batch_data.get("database_results", []), records):
        writer.writerow([row.get(key, "") for key, _ in ETOL_SPECIES_EXPORT_COLUMNS])
    return buffer.getvalue()


def etol_probe_counts_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render full eToL per-probe exact-hit counts (every probe) as CSV or TSV text."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in ETOL_PROBE_EXPORT_COLUMNS])
    records = _etol_records_for_batch(batch_data)
    for row in etol_probe_count_rows(batch_data.get("database_results", []), records):
        writer.writerow([row.get(key, "") for key, _ in ETOL_PROBE_EXPORT_COLUMNS])
    return buffer.getvalue()


def etol_contigs_as_fasta(batch_data: dict) -> str:
    """Render every assembled eToL contig across a batch as multi-FASTA text.

    Each header encodes the sample, species/taxon, contig id, read support, and
    length, so the file can be taken straight to the next re-probing or
    species-identification step (e.g. BLAST against a curated rRNA database).
    """
    blocks: list[str] = []
    for database_result in batch_data.get("database_results", []):
        sample = str(database_result.get("display_name", "") or "sample")
        contigs_by_species = database_result.get("contigs") or {}
        for taxon, contigs in contigs_by_species.items():
            for contig in contigs:
                sequence = str(contig.get("sequence", "") or "")
                if not sequence:
                    continue
                header = (
                    f"{sample}|{taxon}|{contig.get('id', '')}"
                    f"|reads={contig.get('num_reads', 0)}|len={len(sequence)}"
                )
                blocks.append(f">{header}\n{wrap_sequence(sequence)}")
    return ("\n".join(blocks) + "\n") if blocks else ""


def apoe_summary_export_value(row: dict, key: str) -> object:
    """Return APOE export values with compatibility for older saved summaries."""
    if key == "sample_database":
        return (
            row.get("sample_database")
            or row.get("database_sample")
            or row.get("sample")
            or row.get("database")
            or ""
        )
    if key == "c_to_t_percent" and row.get(key, "") == "":
        try:
            ae4_t_hits = int(row.get("ae4_t_hits", 0))
            ae2_t_hits = int(row.get("ae2_t_hits", 0))
            total_hits = int(row.get("total_exact_probe_hits", 0))
        except (TypeError, ValueError):
            return ""
        if total_hits <= 0:
            return ""
        return f"{((ae4_t_hits + ae2_t_hits) / total_hits) * 100:.2f}"
    return row.get(key, "")
