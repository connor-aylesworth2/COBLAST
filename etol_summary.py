"""eToL exact-match probe summary helpers.

The electronic Tree of Life (eToL) panel is a large set of 64-mer probes drawn
from 16S/18S rRNA sequences across the Tree of Life (Hu, Haas & Lathe, BMC
Microbiology 2022). This module mirrors :mod:`apoe_summary`, but instead of four
APOE probes it works with the full bundled eToL panel and aggregates exact
probe hits per sample into probe-, species-, and domain-level counts that feed
the batch results page, the CSV/TSV exports, and downstream species plots.
"""

from __future__ import annotations

from collections import OrderedDict
from functools import lru_cache
import re
from typing import Any

from config import resource_path


ETOL_PROBE_FASTA_PATH = resource_path("data", "eToL_probes.fasta")
ETOL_EXACT_MATCH_FILTER = "100% identity and 100% query coverage"

# Reused to label a sample by its SRA accession when one is present in the
# database name/path; otherwise we fall back to the display name.
ETOL_ACCESSION_PATTERN = re.compile(r"\b(?:SRX|ERX|DRX)\d+\b", re.IGNORECASE)

# Probe headers look like ``A_Hsalinarum_16S_3`` (Archaeon) or ``PGK1_2`` (human
# housekeeping control). The leading token before the first underscore is the
# eToL class code (A, B0..B6, C1..C4, D, E0, F0..F6, H0..H3); the header minus
# the trailing ``_<n>`` index is the species/probe group ("taxon"). Domain
# labels are a convenience grouping derived from the class-code first letter per
# the eToL paper's scheme and can be refined without touching probe counts.
ETOL_DOMAIN_BY_LETTER = {
    "A": "Archaea",
    "B": "Bacteria",
    "C": "Chloroplastida",
    "D": "Basal Eukaryota",
    "E": "Eukaryota",
    "F": "Fungi",
    "H": "Holozoa/Metazoa",
}
# Human housekeeping/normalization probes are not microbial taxa.
ETOL_CONTROL_GROUPS = {"PGK1", "hNSE"}

ETOL_SPECIES_EXPORT_COLUMNS = [
    ("sample_database", "Sample/Database"),
    ("domain", "Domain"),
    ("group", "Class"),
    ("taxon", "Species/Taxon"),
    ("probes_in_panel", "Probes in panel"),
    ("probes_detected", "Probes detected"),
    ("exact_hits", "Total exact probe hits"),
]

ETOL_PROBE_EXPORT_COLUMNS = [
    ("sample_database", "Sample/Database"),
    ("probe", "Probe"),
    ("taxon", "Species/Taxon"),
    ("group", "Class"),
    ("domain", "Domain"),
    ("exact_hits", "Exact hits"),
]


def _class_code(header: str) -> str:
    """Return the eToL class code (token before the first underscore)."""
    return header.split("_", 1)[0]


def _taxon(header: str) -> str:
    """Return the species/probe group (header minus the trailing _<index>)."""
    return re.sub(r"_\d+$", "", header)


def _domain(group: str) -> str:
    """Map a class code to a coarse domain label for grouping/plots."""
    if group in ETOL_CONTROL_GROUPS:
        return "Human control"
    return ETOL_DOMAIN_BY_LETTER.get(group[:1], "Other")


@lru_cache(maxsize=1)
def load_etol_probe_fasta() -> str:
    """Return the bundled eToL probe FASTA text used as the batch query."""
    return ETOL_PROBE_FASTA_PATH.read_text(encoding="utf-8")


@lru_cache(maxsize=1)
def etol_probe_records() -> tuple[dict[str, str], ...]:
    """Parse the bundled FASTA into ordered probe metadata dictionaries."""
    records: list[dict[str, str]] = []
    for line in load_etol_probe_fasta().splitlines():
        line = line.strip()
        if not line.startswith(">"):
            continue
        probe_id = line[1:].strip()
        if not probe_id:
            continue
        group = _class_code(probe_id)
        records.append(
            {
                "probe": probe_id,
                "taxon": _taxon(probe_id),
                "group": group,
                "domain": _domain(group),
            }
        )
    return tuple(records)


@lru_cache(maxsize=1)
def etol_probe_query_ids() -> frozenset[str]:
    """Return the query IDs that belong to the bundled eToL probe panel."""
    return frozenset(record["probe"] for record in etol_probe_records())


def etol_probe_count() -> int:
    """Return how many probes are in the bundled eToL panel."""
    return len(etol_probe_records())


@lru_cache(maxsize=1)
def etol_taxa() -> tuple[dict[str, str], ...]:
    """Return ordered, de-duplicated species/taxon metadata for the panel."""
    taxa: "OrderedDict[str, dict[str, str]]" = OrderedDict()
    for record in etol_probe_records():
        taxon = record["taxon"]
        if taxon not in taxa:
            taxa[taxon] = {
                "taxon": taxon,
                "group": record["group"],
                "domain": record["domain"],
                "probes_in_panel": 0,
            }
        taxa[taxon]["probes_in_panel"] += 1
    return tuple(taxa.values())


def _sample_label(database_result: dict[str, Any]) -> str:
    """Pick a human-friendly sample label, preferring an SRA accession."""
    search_text = " ".join(
        str(database_result.get(key, ""))
        for key in ("display_name", "db_prefix_path")
    )
    accession_match = ETOL_ACCESSION_PATTERN.search(search_text)
    if accession_match:
        return accession_match.group(0).upper()
    return (
        str(database_result.get("display_name") or "").strip()
        or f"Database {database_result.get('database_id', '')}".strip()
        or "Unknown sample"
    )


def _probe_counts(database_result: dict[str, Any]) -> dict[str, int]:
    """Count exact probe hits per query ID for one searched database."""
    query_ids = etol_probe_query_ids()
    counts: dict[str, int] = {}
    for hit in database_result.get("hits", []):
        query_id = hit.get("qseqid", "")
        if query_id in query_ids:
            counts[query_id] = counts.get(query_id, 0) + 1
    return counts


def build_etol_probe_summary(database_results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarize exact eToL probe hits per sample for the results page.

    Each row aggregates the panel down to the species/taxa that were actually
    detected (at least one exact probe hit), which is the "species found in the
    sample" view the workflow ultimately plots.
    """
    rows = []
    for database_result in database_results:
        counts = _probe_counts(database_result)

        detected_species = []
        for taxon in etol_taxa():
            taxon_probes = [
                record["probe"]
                for record in etol_probe_records()
                if record["taxon"] == taxon["taxon"]
            ]
            exact_hits = sum(counts.get(probe, 0) for probe in taxon_probes)
            probes_detected = sum(1 for probe in taxon_probes if counts.get(probe, 0) > 0)
            if exact_hits > 0:
                detected_species.append(
                    {
                        "taxon": taxon["taxon"],
                        "group": taxon["group"],
                        "domain": taxon["domain"],
                        "probes_in_panel": taxon["probes_in_panel"],
                        "probes_detected": probes_detected,
                        "exact_hits": exact_hits,
                    }
                )

        detected_species.sort(key=lambda item: (-item["exact_hits"], item["taxon"]))
        total_exact_probe_hits = sum(counts.values())
        sample = _sample_label(database_result)
        error = str(database_result.get("error") or "")

        rows.append(
            {
                "sample": sample,
                "sample_database": sample,
                "database": database_result.get("display_name", ""),
                "db_prefix_path": database_result.get("db_prefix_path", ""),
                "total_probes": etol_probe_count(),
                "probes_detected": sum(1 for value in counts.values() if value > 0),
                "species_total": len(etol_taxa()),
                "species_detected": len(detected_species),
                "total_exact_probe_hits": total_exact_probe_hits,
                "detected_species": detected_species,
                "status": error or f"{total_exact_probe_hits} exact hit(s)",
                "error": error,
            }
        )
    return rows


def etol_probe_count_rows(database_results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Flat per-probe count rows (every probe, including zeros) for export."""
    rows = []
    for database_result in database_results:
        counts = _probe_counts(database_result)
        sample = _sample_label(database_result)
        for record in etol_probe_records():
            rows.append(
                {
                    "sample_database": sample,
                    "probe": record["probe"],
                    "taxon": record["taxon"],
                    "group": record["group"],
                    "domain": record["domain"],
                    "exact_hits": counts.get(record["probe"], 0),
                }
            )
    return rows


def etol_species_count_rows(database_results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Flat per-species count rows (every taxon, including zeros) for export."""
    rows = []
    for database_result in database_results:
        counts = _probe_counts(database_result)
        sample = _sample_label(database_result)
        for taxon in etol_taxa():
            taxon_probes = [
                record["probe"]
                for record in etol_probe_records()
                if record["taxon"] == taxon["taxon"]
            ]
            exact_hits = sum(counts.get(probe, 0) for probe in taxon_probes)
            probes_detected = sum(1 for probe in taxon_probes if counts.get(probe, 0) > 0)
            rows.append(
                {
                    "sample_database": sample,
                    "domain": taxon["domain"],
                    "group": taxon["group"],
                    "taxon": taxon["taxon"],
                    "probes_in_panel": taxon["probes_in_panel"],
                    "probes_detected": probes_detected,
                    "exact_hits": exact_hits,
                }
            )
    return rows
