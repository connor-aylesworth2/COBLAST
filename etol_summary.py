"""eToL exact-match probe summary helpers.

The electronic Tree of Life (eToL) panel is a set of 64-mer probes drawn from
16S/18S rRNA sequences across the Tree of Life (Hu, Haas & Lathe, BMC
Microbiology 2022;22:317). This module mirrors :mod:`apoe_summary`, but works
with the bundled eToL panel and aggregates exact probe hits per sample into
probe-, species-, and domain-level counts that feed the batch results page, the
CSV/TSV exports, and downstream species plots.

Three eToL batch presets share this machinery, differing only in which probes
are used as the query (and therefore counted):

* ``etol_full``    - the full microbial panel (no human controls).
* ``etol_control`` - the human housekeeping control probes only (PGK1, hNSE).
* ``etol_quick``   - one probe per species (the first probe of each species),
  a slim panel for fast test runs.
"""

from __future__ import annotations

from collections import OrderedDict
from functools import lru_cache
import re
from typing import Any

from blast_runner import _parse_panel_fasta, _pairs_to_fasta
from config import resource_path


ETOL_FULL_FASTA_PATH = resource_path("data", "eToL_probes.fasta")
ETOL_CONTROL_FASTA_PATH = resource_path("data", "eToL_control_probes.fasta")
# The eToL panels keep the paper's permissive net rather than exact matches:
# default megablast with no identity or coverage filter, gated on E-value < 0.01
# (Hu, Haas & Lathe 2022), so partial and mismatched rRNA matches are retained
# but statistically insignificant ones are not.
ETOL_NET_FILTER = "default megablast net, E-value < 0.01 (no identity or coverage filter)"

# Reused to label a sample by its SRA accession when one is present in the
# database name/path; otherwise we fall back to the display name.
ETOL_ACCESSION_PATTERN = re.compile(r"\b(?:SRX|ERX|DRX)\d+\b", re.IGNORECASE)

# Probe headers look like ``B0_Tmaritima_16S_3`` (Bacterium) or ``PGK1_2`` (human
# housekeeping control). The leading token before the first underscore is the
# eToL class code (A, B0..B6, C1..C4, D, E0, F0..F6, H0..H3); the header minus
# the trailing ``_<n>`` index is the species/probe group ("taxon"). Domain
# labels follow the class-code scheme defined in Hu, Haas & Lathe 2022:
# A Archaea; B Bacteria; C Chloroplastida; D Amoebozoa; E0 basal Eukaryota;
# F Fungi; H Holozoa/Metazoa.
ETOL_DOMAIN_BY_LETTER = {
    "A": "Archaea",
    "B": "Bacteria",
    "C": "Chloroplastida",
    "D": "Amoebozoa",
    "E": "Basal Eukaryota",
    "F": "Fungi",
    "H": "Holozoa/Metazoa",
}
# Human housekeeping/normalization probes are not microbial taxa.
ETOL_CONTROL_GROUPS = {"PGK1", "hNSE"}

# Species display name: drop the leading class code (e.g. "B0_") and the trailing
# rRNA subunit suffix (e.g. "_16S"/"_18S") so only the species label remains.
_CLASS_PREFIX_RE = re.compile(r"^[^_]+_")
_RRNA_SUFFIX_RE = re.compile(r"_\d+S$")

ETOL_SPECIES_EXPORT_COLUMNS = [
    ("sample_database", "Sample/Database"),
    ("domain", "Domain"),
    ("group", "Class"),
    ("species", "Species/Taxon"),
    ("probes_in_panel", "Probes in panel"),
    ("probes_detected", "Probes detected"),
    ("exact_hits", "Total exact probe hits"),
    ("host_cells", "Est. host cells"),
    ("normalized_abundance", "Reads per host cell"),
]

ETOL_PROBE_EXPORT_COLUMNS = [
    ("sample_database", "Sample/Database"),
    ("probe", "Probe"),
    ("species", "Species/Taxon"),
    ("group", "Class"),
    ("domain", "Domain"),
    ("exact_hits", "Exact hits"),
    ("normalized_abundance", "Reads per host cell"),
]


def _class_code(header: str) -> str:
    """Return the eToL class code (token before the first underscore)."""
    return header.split("_", 1)[0]


def _taxon(header: str) -> str:
    """Return the species/probe group (header minus the trailing _<index>)."""
    return re.sub(r"_\d+$", "", header)


def _domain(group: str) -> str:
    """Map a class code to its domain label per the eToL paper's scheme."""
    if group in ETOL_CONTROL_GROUPS:
        return "Human control"
    return ETOL_DOMAIN_BY_LETTER.get(group[:1], "Other")


def _species(taxon: str) -> str:
    """Return just the species label (no class prefix, no rRNA-unit suffix)."""
    name = _CLASS_PREFIX_RE.sub("", taxon)
    name = _RRNA_SUFFIX_RE.sub("", name)
    return name or taxon


def _record_meta(probe: str) -> dict[str, str]:
    """Build the metadata dictionary used throughout summaries/exports."""
    group = _class_code(probe)
    taxon = _taxon(probe)
    return {
        "probe": probe,
        "taxon": taxon,
        "group": group,
        "domain": _domain(group),
        "species": _species(taxon),
    }


@lru_cache(maxsize=1)
def _full_pairs() -> tuple[tuple[str, str], ...]:
    return tuple(_parse_panel_fasta(ETOL_FULL_FASTA_PATH.read_text(encoding="utf-8")))


@lru_cache(maxsize=1)
def _control_pairs() -> tuple[tuple[str, str], ...]:
    return tuple(_parse_panel_fasta(ETOL_CONTROL_FASTA_PATH.read_text(encoding="utf-8")))


@lru_cache(maxsize=1)
def _quick_pairs() -> tuple[tuple[str, str], ...]:
    """Return the first probe of each species in the full microbial panel."""
    seen: set[str] = set()
    chosen: list[tuple[str, str]] = []
    for header, sequence in _full_pairs():
        taxon = _taxon(header)
        if taxon not in seen:
            seen.add(taxon)
            chosen.append((header, sequence))
    return tuple(chosen)


# Each preset: form field name, UI labels/description, and the probe source.
ETOL_PRESETS: "OrderedDict[str, dict[str, Any]]" = OrderedDict(
    [
        (
            "etol_full",
            {
                "form_field": "etol_probe_preset",
                "label": "eToL Full probe batch",
                "short_label": "eToL Full",
                "panel_label": "eToL Full",
                "microbial": True,
                "description": (
                    "Microbial electronic Tree of Life panel across Archaea, "
                    "Bacteria, Chloroplastida, Amoebozoa, basal Eukaryota, Fungi, "
                    "and Holozoa/Metazoa; BLASTN default megablast net (E-value "
                    "< 0.01, no identity or coverage filter). Counts matching reads per probe "
                    "and species, de-duplicates reads across probes, and reports "
                    "reads per host cell using the PGK1/hNSE control probes."
                ),
                "pairs": _full_pairs,
            },
        ),
        (
            "etol_control",
            {
                "form_field": "etol_control_probe_preset",
                "label": "eToL Control probe batch",
                "short_label": "eToL Control (human control)",
                "panel_label": "eToL Control",
                "microbial": False,
                "description": (
                    "Same as the eToL Full preset, but only uses human sequences "
                    "(the PGK1 and hNSE housekeeping probes) for a control."
                ),
                "pairs": _control_pairs,
            },
        ),
        (
            "etol_quick",
            {
                "form_field": "etol_quick_probe_preset",
                "label": "eToL Quick probe batch",
                "short_label": "eToL Quick (one probe per species)",
                "panel_label": "eToL Quick",
                "microbial": True,
                "description": (
                    "Same as the eToL Full preset, but uses only the first probe "
                    "of each species (one probe per species) for fast test runs."
                ),
                "pairs": _quick_pairs,
            },
        ),
    ]
)


def etol_preset_keys() -> list[str]:
    """Return the eToL preset keys in display order."""
    return list(ETOL_PRESETS)


def etol_preset_form_field(key: str) -> str:
    """Return the form field name for a preset."""
    return ETOL_PRESETS[key]["form_field"]


def etol_preset_label(key: str) -> str:
    """Return the short label shown on the results page."""
    return ETOL_PRESETS[key]["short_label"]


def etol_preset_is_microbial(key: str) -> bool:
    """Return True for microbial panels (eligible for human-read filtering)."""
    return bool(ETOL_PRESETS[key]["microbial"])


def etol_preset_fasta(key: str) -> str:
    """Return the query FASTA text for a preset."""
    return _pairs_to_fasta(ETOL_PRESETS[key]["pairs"]())


@lru_cache(maxsize=None)
def etol_preset_records(key: str) -> tuple[dict[str, str], ...]:
    """Return ordered probe metadata for a preset's panel."""
    return tuple(_record_meta(header) for header, _ in ETOL_PRESETS[key]["pairs"]())


def etol_preset_query_ids(key: str) -> frozenset[str]:
    """Return the query IDs that belong to a preset's panel."""
    return frozenset(record["probe"] for record in etol_preset_records(key))


def etol_preset_probe_count(key: str) -> int:
    """Return how many probes are in a preset's panel."""
    return len(etol_preset_records(key))


# --- Host-cell normalization (Hu, Haas & Lathe 2022) -----------------------
# Each microbial panel is searched together with the housekeeping control probes
# (PGK1, hNSE). Dividing microbial read counts by the estimated host-cell count
# -- the mean per-gene control readcount over ~50 transcripts per cell -- yields
# the paper's "reads per host cell" abundance, normalizing for how much host
# material each library represents.
HOST_TRANSCRIPTS_PER_CELL = 50.0


@lru_cache(maxsize=1)
def etol_control_records() -> tuple[dict[str, str], ...]:
    """Return metadata for the housekeeping control probes (PGK1, hNSE)."""
    return tuple(_record_meta(header) for header, _ in _control_pairs())


def etol_control_query_ids() -> frozenset[str]:
    """Return the control-probe query ids used for host normalization."""
    return frozenset(record["probe"] for record in etol_control_records())


def etol_search_pairs(key: str) -> tuple[tuple[str, str], ...]:
    """Return the (header, sequence) probes actually BLASTed for a preset.

    Microbial presets append the housekeeping control probes so a single search
    yields both the microbial net and the host-normalization counts (the control
    counts are the host-cell normalization denominator). The standalone control
    preset is searched as-is.
    """
    pairs = tuple(ETOL_PRESETS[key]["pairs"]())
    if etol_preset_is_microbial(key):
        pairs = pairs + _control_pairs()
    return pairs


def etol_search_fasta(key: str) -> str:
    """Return the FASTA actually sent to BLAST (microbial panel + controls)."""
    return _pairs_to_fasta(etol_search_pairs(key))


def etol_search_query_ids(key: str) -> frozenset[str]:
    """Return every query id in the searched FASTA (microbial + control)."""
    return frozenset(header for header, _ in etol_search_pairs(key))


def control_gene_means(control_counts: dict[str, int]) -> "OrderedDict[str, float]":
    """Mean read count per housekeeping gene (PGK1, hNSE) across its probes."""
    by_gene: "OrderedDict[str, list[int]]" = OrderedDict()
    for probe, count in (control_counts or {}).items():
        by_gene.setdefault(_taxon(probe), []).append(int(count))
    return OrderedDict(
        (gene, sum(values) / len(values)) for gene, values in by_gene.items() if values
    )


def compute_host_cells(control_counts: dict[str, int]) -> float:
    """Estimate host-cell count from housekeeping control-probe read counts.

    Per Hu, Haas & Lathe 2022: average the read counts of each gene's probes,
    take the mean across genes, then divide by ~50 transcripts per host cell.
    Returns 0.0 when no control reads were counted (normalization undefined).
    """
    gene_means = list(control_gene_means(control_counts).values())
    if not gene_means:
        return 0.0
    return (sum(gene_means) / len(gene_means)) / HOST_TRANSCRIPTS_PER_CELL


def normalized_abundance(raw_count: int, host_cells: float) -> float | None:
    """Reads per host cell, or None when the host-cell estimate is unavailable."""
    if not host_cells or host_cells <= 0:
        return None
    return raw_count / host_cells


def _format_host_cells(host_cells: float) -> str:
    """Render the host-cell estimate for export/display (blank when undefined)."""
    return "" if not host_cells or host_cells <= 0 else f"{host_cells:.4g}"


def _format_normalized(value: float | None) -> str:
    """Render a normalized abundance for export/display (blank when undefined)."""
    return "" if value is None else f"{value:.4g}"


def etol_preset_options() -> list[dict[str, Any]]:
    """Return per-preset metadata (incl. probe counts) for the batch form."""
    options = []
    for key, preset in ETOL_PRESETS.items():
        options.append(
            {
                "key": key,
                "form_field": preset["form_field"],
                "label": preset["label"],
                "short_label": preset["short_label"],
                "panel_label": preset["panel_label"],
                "description": preset["description"],
                "microbial": preset["microbial"],
                "probe_count": etol_preset_probe_count(key),
            }
        )
    return options


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


def _taxa_from(records: tuple[dict[str, str], ...]) -> list[dict[str, Any]]:
    """Collapse probe records into ordered, de-duplicated species/taxon rows."""
    taxa: "OrderedDict[str, dict[str, Any]]" = OrderedDict()
    for record in records:
        taxon = record["taxon"]
        if taxon not in taxa:
            taxa[taxon] = {
                "taxon": taxon,
                "species": record["species"],
                "group": record["group"],
                "domain": record["domain"],
                "probes_in_panel": 0,
            }
        taxa[taxon]["probes_in_panel"] += 1
    return list(taxa.values())


def _probe_counts(database_result: dict[str, Any], query_ids: frozenset[str]) -> dict[str, int]:
    """Count exact probe hits per query ID for one searched database."""
    counts: dict[str, int] = {}
    for hit in database_result.get("hits", []):
        query_id = hit.get("qseqid", "")
        if query_id in query_ids:
            counts[query_id] = counts.get(query_id, 0) + 1
    return counts


def _probes_by_taxon(records: tuple[dict[str, str], ...]) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    for record in records:
        grouped.setdefault(record["taxon"], []).append(record["probe"])
    return grouped


def group_read_ids_by_taxon(
    hits: list[dict[str, Any]]
) -> "OrderedDict[str, list[str]]":
    """Group matched-read ids (``sseqid``) by the species/taxon of their probe.

    The contig-assembly step assembles the reads each species' probes recovered
    together (per taxon), mirroring the paper's per-group assembly for species
    identification. Read ids are de-duplicated within a taxon and first-seen
    order is preserved so the assembler input is stable.
    """
    grouped: "OrderedDict[str, list[str]]" = OrderedDict()
    seen: dict[str, set[str]] = {}
    for hit in hits:
        read_id = str(hit.get("sseqid", "") or "")
        if not read_id:
            continue
        taxon = _taxon(str(hit.get("qseqid", "") or ""))
        if taxon not in grouped:
            grouped[taxon] = []
            seen[taxon] = set()
        if read_id not in seen[taxon]:
            seen[taxon].add(read_id)
            grouped[taxon].append(read_id)
    return grouped


def build_etol_probe_summary(
    database_results: list[dict[str, Any]], records: tuple[dict[str, str], ...]
) -> list[dict[str, Any]]:
    """Summarize exact eToL probe hits per sample for the results page.

    Each row aggregates the panel down to the species/taxa that were actually
    detected (at least one exact probe hit), which is the "species found in the
    sample" view the workflow ultimately plots.
    """
    query_ids = frozenset(record["probe"] for record in records)
    taxa = _taxa_from(records)
    probes_by_taxon = _probes_by_taxon(records)

    rows = []
    for database_result in database_results:
        counts = _probe_counts(database_result, query_ids)
        # Host-cell normalization uses the housekeeping control-probe counts that
        # the microbial search recorded separately (never human-filtered).
        control_counts = database_result.get("etol_control_counts") or {}
        host_cells = compute_host_cells(control_counts)
        gene_means = control_gene_means(control_counts)

        detected_species = []
        for taxon in taxa:
            taxon_probes = probes_by_taxon[taxon["taxon"]]
            exact_hits = sum(counts.get(probe, 0) for probe in taxon_probes)
            probes_detected = sum(1 for probe in taxon_probes if counts.get(probe, 0) > 0)
            if exact_hits > 0:
                detected_species.append(
                    {
                        "taxon": taxon["taxon"],
                        "species": taxon["species"],
                        "group": taxon["group"],
                        "domain": taxon["domain"],
                        "probes_in_panel": taxon["probes_in_panel"],
                        "probes_detected": probes_detected,
                        "exact_hits": exact_hits,
                        "normalized_abundance": normalized_abundance(exact_hits, host_cells),
                        "normalized_label": _format_normalized(
                            normalized_abundance(exact_hits, host_cells)
                        ),
                    }
                )

        detected_species.sort(key=lambda item: (-item["exact_hits"], item["species"]))
        total_exact_probe_hits = sum(counts.values())
        sample = _sample_label(database_result)
        error = str(database_result.get("error") or "")

        rows.append(
            {
                "sample": sample,
                "sample_database": sample,
                "database": database_result.get("display_name", ""),
                "db_prefix_path": database_result.get("db_prefix_path", ""),
                "total_probes": len(records),
                "probes_detected": sum(1 for value in counts.values() if value > 0),
                "species_total": len(taxa),
                "species_detected": len(detected_species),
                "total_exact_probe_hits": total_exact_probe_hits,
                "detected_species": detected_species,
                "host_cells": host_cells,
                "host_cells_label": _format_host_cells(host_cells),
                "control_gene_means": dict(gene_means),
                "normalized": host_cells > 0,
                "status": error or f"{total_exact_probe_hits} matched read(s)",
                "error": error,
            }
        )
    return rows


def etol_probe_count_rows(
    database_results: list[dict[str, Any]], records: tuple[dict[str, str], ...]
) -> list[dict[str, Any]]:
    """Flat per-probe count rows (every probe, including zeros) for export."""
    query_ids = frozenset(record["probe"] for record in records)
    rows = []
    for database_result in database_results:
        counts = _probe_counts(database_result, query_ids)
        host_cells = compute_host_cells(database_result.get("etol_control_counts") or {})
        sample = _sample_label(database_result)
        for record in records:
            raw = counts.get(record["probe"], 0)
            rows.append(
                {
                    "sample_database": sample,
                    "probe": record["probe"],
                    "species": record["species"],
                    "group": record["group"],
                    "domain": record["domain"],
                    "exact_hits": raw,
                    "normalized_abundance": _format_normalized(
                        normalized_abundance(raw, host_cells)
                    ),
                }
            )
    return rows


def etol_species_count_rows(
    database_results: list[dict[str, Any]], records: tuple[dict[str, str], ...]
) -> list[dict[str, Any]]:
    """Flat per-species count rows (every taxon, including zeros) for export."""
    query_ids = frozenset(record["probe"] for record in records)
    taxa = _taxa_from(records)
    probes_by_taxon = _probes_by_taxon(records)
    rows = []
    for database_result in database_results:
        counts = _probe_counts(database_result, query_ids)
        host_cells = compute_host_cells(database_result.get("etol_control_counts") or {})
        sample = _sample_label(database_result)
        for taxon in taxa:
            taxon_probes = probes_by_taxon[taxon["taxon"]]
            exact_hits = sum(counts.get(probe, 0) for probe in taxon_probes)
            probes_detected = sum(1 for probe in taxon_probes if counts.get(probe, 0) > 0)
            rows.append(
                {
                    "sample_database": sample,
                    "domain": taxon["domain"],
                    "group": taxon["group"],
                    "species": taxon["species"],
                    "probes_in_panel": taxon["probes_in_panel"],
                    "probes_detected": probes_detected,
                    "exact_hits": exact_hits,
                    "host_cells": _format_host_cells(host_cells),
                    "normalized_abundance": _format_normalized(
                        normalized_abundance(exact_hits, host_cells)
                    ),
                }
            )
    return rows
