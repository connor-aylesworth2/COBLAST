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
    build_etol_matrix,
    etol_preset_records,
    etol_probe_count_rows,
    etol_species_count_rows,
)
from etol_validation import compute_confusion
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


def _columns_delimited(columns, rows, delimiter: str) -> str:
    """Render ``rows`` as CSV/TSV using ``(key, header)`` column pairs."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in columns])
    for row in rows:
        writer.writerow([row.get(key, "") for key, _ in columns])
    return buffer.getvalue()


def result_rows_as_delimited(result_data: dict, delimiter: str) -> str:
    """Render saved hits as CSV or TSV text."""
    return _columns_delimited(RESULT_COLUMNS, result_data.get("hits", []), delimiter)


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


def batch_summary_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render the top-of-page Batch Summary panel as a two-column CSV/TSV.

    One ``Statistic, Value`` row per metric shown in the results-page summary
    panel, including the same rows the panel only renders conditionally (probe
    preset, host normalization, contig assembly, human filter, etc.). This is the
    "export what's on screen" companion to the richer per-hit / per-probe
    downloads, so the headline numbers can be pasted into a methods table.
    """
    rows: list[tuple[str, object]] = [
        ("Program", batch_data.get("program", "")),
        ("Databases", len(batch_data.get("database_results", []))),
        (
            "Total runtime (BLAST search, summed across databases, seconds)",
            f"{batch_data.get('total_runtime_seconds', 0.0):.3f}",
        ),
    ]
    if batch_data.get("wall_clock_seconds") is not None:
        rows.append(
            (
                "Wall-clock elapsed time (seconds)",
                f"{batch_data.get('wall_clock_seconds', 0.0):.3f}",
            )
        )
    if batch_data.get("batch_workers") is not None:
        rows.append(("Concurrency (databases at a time)", batch_data.get("batch_workers")))
    rows.append(("Query records", batch_data.get("query_count", 0)))
    rows.append(("Query total bases/residues", batch_data.get("query_total_length", 0)))
    rows.append(("Total hits", batch_data.get("total_hits", 0)))

    if batch_data.get("apoe_probe_preset"):
        rows.append(("Probe preset", "APOE exact-match probes"))
        rows.append(("Hit filter", batch_data.get("hit_filter", "")))
    elif batch_data.get("etol_probe_preset"):
        rows.append(
            ("Probe preset", batch_data.get("etol_preset_label") or "eToL net probes")
        )
        rows.append(("Hit filter", batch_data.get("hit_filter", "")))
        if batch_data.get("etol_dedup_removed"):
            rows.append(
                ("Reads de-duplicated (reallocated to best probe)", batch_data["etol_dedup_removed"])
            )
        if batch_data.get("etol_normalized"):
            rows.append(("Host normalization", "reads per host cell (PGK1/hNSE control probes)"))
        if batch_data.get("assemble_contigs"):
            rows.append(("Contigs assembled (CAP3)", batch_data.get("contig_count", 0)))
        if batch_data.get("reprobe_contigs"):
            rows.append(("Re-probing extra reads recovered", batch_data.get("reprobe_new_reads", 0)))
        if batch_data.get("identify_contigs"):
            rows.append(("Contigs identified", batch_data.get("contigs_identified", 0)))
            rows.append(("Contig species-ID database", batch_data.get("species_id_db", "")))
        if batch_data.get("contig_assembly_unavailable"):
            rows.append(("Contig assembly", "CAP3 not found - assembly skipped"))

    if batch_data.get("human_filter_enabled"):
        rows.append(("Human filter hits removed", batch_data.get("human_filter_hits_removed", 0)))
        rows.append(("Human filter database", batch_data.get("human_filter_db", "")))

    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow(["Statistic", "Value"])
    for label, value in rows:
        writer.writerow([label, value])
    return buffer.getvalue()


def apoe_summary_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render APOE per-sample probe summaries as CSV or TSV text."""
    summary_rows = batch_data.get("apoe_probe_summary")
    if summary_rows is None:
        summary_rows = build_apoe_probe_summary(batch_data.get("database_results", []))
    return _columns_delimited(APOE_SUMMARY_EXPORT_COLUMNS, summary_rows, delimiter)


def _etol_records_for_batch(batch_data: dict) -> tuple:
    """Resolve the probe panel used by a saved eToL batch (defaults to full)."""
    return etol_preset_records(batch_data.get("etol_preset_key") or "etol_full")


def etol_summary_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render eToL per-species exact-hit counts (every taxon) as CSV or TSV text."""
    records = _etol_records_for_batch(batch_data)
    rows = etol_species_count_rows(batch_data.get("database_results", []), records)
    return _columns_delimited(ETOL_SPECIES_EXPORT_COLUMNS, rows, delimiter)


def etol_probe_counts_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render full eToL per-probe exact-hit counts (every probe) as CSV or TSV text."""
    records = _etol_records_for_batch(batch_data)
    rows = etol_probe_count_rows(batch_data.get("database_results", []), records)
    return _columns_delimited(ETOL_PROBE_EXPORT_COLUMNS, rows, delimiter)


ETOL_CONFUSION_COLUMNS = [
    ("result", "Result"),
    ("virus", "Virus (WGS)"),
    ("sample", "Sample"),
    ("srx", "SRX"),
    ("wgs_count", "WGS count"),
    ("actual", "WGS present"),
    ("predicted", "eToL-V predicted"),
    ("raw_hits", "Raw net hits"),
    ("confirmed_hits", "Validated hits"),
]
# Group the per-cell rows so true/false positives and negatives read in order.
_CONFUSION_RESULT_ORDER = {"TP": 0, "FN": 1, "FP": 2, "TN": 3}


def etol_confusion_rows_as_delimited(batch_data: dict, delimiter: str) -> str:
    """Render the eToL-V confusion matrix as a per-cell CSV/TSV.

    One row per scored (virus, sample) cell, carrying both the raw net hit count
    and the validated (contig-confirmed) count so a false negative shows where it
    was lost (raw 0 = the net E-value gate; raw > 0 but validated 0 = contig
    assembly/identification). The 2x2 totals and metrics are derivable by pivoting
    on the Result column.
    """
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow([label for _, label in ETOL_CONFUSION_COLUMNS])

    records = _etol_records_for_batch(batch_data)
    matrix = build_etol_matrix(batch_data.get("database_results", []), records)
    confusion = compute_confusion(matrix)
    cells = sorted(
        confusion.get("cells", []),
        key=lambda c: (_CONFUSION_RESULT_ORDER.get(c["result"], 9), c["virus"], c["sample"]),
    )
    for cell in cells:
        writer.writerow(
            [
                cell["result"],
                cell["virus"],
                cell["sample"],
                cell["srx"],
                cell["wgs_count"],
                "yes" if cell["actual"] else "no",
                "yes" if cell["predicted"] else "no",
                cell["raw_hits"],
                cell["confirmed_hits"],
            ]
        )
    return buffer.getvalue()


def etol_matrix_payload(batch_data: dict, level: str = "species") -> dict:
    """Build the plot-ready eToL hit matrix for a saved batch (heatmap source).

    Reshapes the same per-probe/per-species counts the CSV exports use into a
    dense ``rows x samples`` matrix and tags it with the preset so the client can
    pick paper-faithful defaults (raw hit counts for the viral panel, log2 reads
    per host cell for the cellular panels).
    """
    records = _etol_records_for_batch(batch_data)
    design = batch_data.get("design_matrix") or None
    matrix = build_etol_matrix(
        batch_data.get("database_results", []),
        records,
        level=level,
        condition_index=design,
    )
    preset = batch_data.get("etol_preset_key") or "etol_full"
    matrix["preset"] = preset
    matrix["preset_label"] = batch_data.get("etol_preset_label") or ""
    matrix["is_viral"] = preset == "etol_v"
    # Surface the applied design matrix (source, label set, and any samples it
    # did not cover) so the heatmap can show a banner + a label-driven legend.
    if design:
        matrix["design_matrix"] = {
            "source": design.get("source", ""),
            "conditions": design.get("conditions", []),
            "row_count": design.get("row_count", 0),
            "unmatched_samples": matrix.get("unmatched_samples") or [],
            "warnings": design.get("warnings", []),
        }
    return matrix


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
                # Append contig identification when the run produced it.
                homolog = str(contig.get("closest_homolog", "") or "")
                if homolog:
                    header += (
                        f"|confirmed={contig.get('confirmed_reads', 0)}"
                        f"|homolog={homolog}"
                    )
                blocks.append(f">{header}\n{wrap_sequence(sequence)}")
    return ("\n".join(blocks) + "\n") if blocks else ""
