"""Design-matrix parsing for eToL heatmap condition labels.

The eToL Result Heatmap annotates each sample column with a condition swatch
(AD vs control, etc.). Historically that label was guessed by regex from the
database name/path (:func:`etol_summary._sample_condition`), which silently
mislabels samples -- e.g. the auto-generated ``"SRA <acc> reads"`` name makes
every sample match the bare ``AD`` alternative ("re**ad**s").

This module lets the user upload an explicit *design matrix* on the batch BLAST
page that maps each sample to a condition. When supplied it is authoritative:
:func:`etol_summary.build_etol_matrix` looks the condition up here instead of
running the regex.

Format (strict): CSV or TSV with a header row containing (case-insensitively,
in any order) a ``sample`` column and a ``condition`` column. ``sample`` is an
SRA accession (SRR/SRX/ERR/...) or a database display name; ``condition`` is the
free-text label to show. One row per sample. Extra columns are ignored, leaving
room for a future multi-factor mode.
"""

from __future__ import annotations

import csv
import io

# Reuse the same accession grammar the heatmap uses to label columns, so a row
# keyed by an accession matches whichever accession the column resolves to.
from etol_summary import ETOL_ACCESSION_PATTERN

REQUIRED_COLUMNS = ("sample", "condition")


class DesignMatrixError(ValueError):
    """A fatal problem with an uploaded design matrix (rejected before BLAST)."""


def _detect_delimiter(filename: str, header_line: str) -> str:
    """Pick the delimiter from the file extension, falling back to a sniff."""
    name = (filename or "").lower()
    if name.endswith((".tsv", ".tab")):
        return "\t"
    if name.endswith(".csv"):
        return ","
    # No decisive extension (e.g. ``.txt``): prefer tab only when the header row
    # is clearly tab-delimited, otherwise default to comma.
    if "\t" in header_line and "," not in header_line:
        return "\t"
    return ","


def parse_design_matrix(text: str, *, filename: str = "") -> dict:
    """Parse design-matrix text into a JSON-serializable lookup index.

    Returns a dict with ``by_accession`` (ACCESSION -> label, accession upper-
    cased), ``by_name`` (lower-cased sample id -> label), the distinct
    ``conditions`` in first-seen order, the matched ``row_count``, the parsed
    ``source`` filename, and non-fatal ``warnings``. Raises
    :class:`DesignMatrixError` on a fatal format problem so the batch route can
    reject the upload up front (before the long BLAST run).
    """
    if text is None:
        raise DesignMatrixError("The design matrix file is empty.")
    # Tolerate a UTF-8 BOM from Windows editors (matches the FASTA upload path).
    text = text.lstrip("﻿")
    if not text.strip():
        raise DesignMatrixError("The design matrix file is empty.")

    first_line = text.splitlines()[0] if text.splitlines() else ""
    delimiter = _detect_delimiter(filename, first_line)
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))

    # The header is the first row that has any content (skip blank leading lines).
    header = None
    data_start = 0
    for index, row in enumerate(rows):
        if any(cell.strip() for cell in row):
            header = row
            data_start = index + 1
            break
    if header is None:
        raise DesignMatrixError("The design matrix file has no header row.")

    normalized = [cell.strip().lower() for cell in header]
    col_index: dict[str, int] = {}
    for name in REQUIRED_COLUMNS:
        if name not in normalized:
            present = ", ".join(cell for cell in header if cell.strip()) or "(none)"
            raise DesignMatrixError(
                f"The design matrix must have a '{name}' column. Found: {present}."
            )
        col_index[name] = normalized.index(name)

    by_accession: dict[str, str] = {}
    by_name: dict[str, str] = {}
    conditions: list[str] = []
    warnings: list[str] = []
    seen: set[str] = set()
    row_count = 0

    for row in rows[data_start:]:
        if not any(cell.strip() for cell in row):
            continue
        sample = row[col_index["sample"]].strip() if col_index["sample"] < len(row) else ""
        condition = (
            row[col_index["condition"]].strip()
            if col_index["condition"] < len(row)
            else ""
        )
        if not sample:
            warnings.append("Skipped a row with no sample identifier.")
            continue
        key = sample.lower()
        if key in seen:
            raise DesignMatrixError(
                f"Duplicate sample '{sample}' in the design matrix; "
                "each sample must appear once."
            )
        seen.add(key)
        row_count += 1
        by_name[key] = condition
        accession_match = ETOL_ACCESSION_PATTERN.search(sample)
        if accession_match:
            by_accession[accession_match.group(0).upper()] = condition
        if not condition:
            warnings.append(f"Sample '{sample}' has no condition; it will render unlabeled.")
        elif condition not in conditions:
            conditions.append(condition)

    if row_count == 0:
        raise DesignMatrixError("The design matrix has a header but no data rows.")

    return {
        "source": filename or "design matrix",
        "by_accession": by_accession,
        "by_name": by_name,
        "conditions": conditions,
        "row_count": row_count,
        "warnings": warnings,
    }
