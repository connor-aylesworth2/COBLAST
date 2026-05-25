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

from blast_runner import BlastResult
from config import resource_root, runtime_data_dir


PROJECT_ROOT = resource_root()
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
