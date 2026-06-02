"""APOE exact-match probe summary helpers."""

from __future__ import annotations

import re
from typing import Any


APOE_PROBE_DEFINITIONS = [
    {
        "query_id": "AE4_E4=C",
        "site": "ae4",
        "base": "c",
        "label": "AE4=C",
        "genotype": "APOE4",
    },
    {
        "query_id": "AE4_E23=T",
        "site": "ae4",
        "base": "t",
        "label": "AE4=T",
        "genotype": "APOE2+E3",
    },
    {
        "query_id": "AE2_E34=C",
        "site": "ae2",
        "base": "c",
        "label": "AE2=C",
        "genotype": "APOE3+E4",
    },
    {
        "query_id": "AE2_E2=T",
        "site": "ae2",
        "base": "t",
        "label": "AE2=T",
        "genotype": "APOE2",
    },
]

APOE_QUERY_IDS = {probe["query_id"] for probe in APOE_PROBE_DEFINITIONS}
APOE_ACCESSION_PATTERN = re.compile(r"\b(?:SRX|ERX|DRX)\d+\b", re.IGNORECASE)

APOE_SUMMARY_EXPORT_COLUMNS = [
    ("sample_database", "Sample/Database"),
    ("ae4_c_hits", "AE4=C hits"),
    ("ae4_t_hits", "AE4=T hits"),
    ("ae2_c_hits", "AE2=C hits"),
    ("ae2_t_hits", "AE2=T hits"),
    ("total_exact_probe_hits", "Total exact probe hits"),
    ("c_to_t_percent", "% C<->T"),
]


def apoe_probe_query_ids() -> set[str]:
    """Return the query IDs that belong to the stored APOE probe preset."""
    return set(APOE_QUERY_IDS)


def _percentage(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return ""
    return f"{(numerator / denominator) * 100:.2f}"


def _sample_label(database_result: dict[str, Any]) -> str:
    search_text = " ".join(
        str(database_result.get(key, ""))
        for key in ("display_name", "db_prefix_path")
    )
    accession_match = APOE_ACCESSION_PATTERN.search(search_text)
    if accession_match:
        return accession_match.group(0).upper()
    return (
        str(database_result.get("display_name") or "").strip()
        or f"Database {database_result.get('database_id', '')}".strip()
        or "Unknown sample"
    )


def build_apoe_probe_summary(database_results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarize exact APOE probe hits for each searched database/sample."""
    rows = []
    for database_result in database_results:
        counts = {query_id: 0 for query_id in APOE_QUERY_IDS}
        for hit in database_result.get("hits", []):
            query_id = hit.get("qseqid", "")
            if query_id in counts:
                counts[query_id] += 1

        ae4_c_hits = counts["AE4_E4=C"]
        ae4_t_hits = counts["AE4_E23=T"]
        ae2_c_hits = counts["AE2_E34=C"]
        ae2_t_hits = counts["AE2_E2=T"]
        ae4_total_hits = ae4_c_hits + ae4_t_hits
        ae2_total_hits = ae2_c_hits + ae2_t_hits
        total_exact_probe_hits = ae4_total_hits + ae2_total_hits
        error = str(database_result.get("error") or "")

        sample = _sample_label(database_result)

        rows.append(
            {
                "sample": sample,
                "sample_database": sample,
                "database_sample": sample,
                "database": database_result.get("display_name", ""),
                "db_prefix_path": database_result.get("db_prefix_path", ""),
                "ae4_c_hits": ae4_c_hits,
                "ae4_t_hits": ae4_t_hits,
                "ae4_total_hits": ae4_total_hits,
                "ae4_t_percent": _percentage(ae4_t_hits, ae4_total_hits),
                "ae2_c_hits": ae2_c_hits,
                "ae2_t_hits": ae2_t_hits,
                "ae2_total_hits": ae2_total_hits,
                "ae2_t_percent": _percentage(ae2_t_hits, ae2_total_hits),
                "total_exact_probe_hits": total_exact_probe_hits,
                "c_to_t_percent": _percentage(ae4_t_hits + ae2_t_hits, total_exact_probe_hits),
                "status": error or f"{database_result.get('hit_count', 0)} exact hit(s)",
                "error": error,
            }
        )
    return rows
