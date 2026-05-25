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


def load_result(run_id: str) -> dict:
    """Load a previously saved result payload."""
    path = result_path(run_id)
    if not path.exists():
        raise FileNotFoundError("Result not found.")
    return json.loads(path.read_text(encoding="utf-8"))


def result_rows_as_delimited(result_data: dict, delimiter: str) -> str:
    """Render saved hits as CSV or TSV text."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in RESULT_COLUMNS])
    for hit in result_data.get("hits", []):
        writer.writerow([hit.get(key, "") for key, _ in RESULT_COLUMNS])
    return buffer.getvalue()
