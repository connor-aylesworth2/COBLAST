"""Unit tests for the APOE and eToL probe-count summaries.

These exercise the pure aggregation logic over synthetic hit lists, so they do
not require BLAST+ or any patient data. The eToL tests read the bundled probe
panels under data/, so run pytest from the repository root.
"""

from apoe_summary import build_apoe_probe_summary
from etol_summary import (
    build_etol_probe_summary,
    etol_preset_probe_count,
    etol_preset_query_ids,
    etol_preset_records,
    etol_probe_count_rows,
    etol_species_count_rows,
)


# --- APOE ------------------------------------------------------------------

def _apoe_sample():
    return {
        "database_id": 1,
        "display_name": "SRX123456 APOE pilot",
        "db_prefix_path": r"C:\COBLAST_data\sra\SRX123456\reads",
        "hit_count": 5,
        "hits": [
            {"qseqid": "AE4_E4=C"},
            {"qseqid": "AE4_E23=T"},
            {"qseqid": "AE4_E23=T"},
            {"qseqid": "AE2_E34=C"},
            {"qseqid": "AE2_E2=T"},
        ],
        "error": "",
    }


def test_apoe_summary_counts_and_ratio():
    row = build_apoe_probe_summary([_apoe_sample()])[0]
    assert row["sample"] == "SRX123456"  # accession pulled from the path
    assert row["ae4_c_hits"] == 1
    assert row["ae4_t_hits"] == 2
    assert row["ae2_c_hits"] == 1
    assert row["ae2_t_hits"] == 1
    assert row["total_exact_probe_hits"] == 5
    # (AE4=T + AE2=T) / total = (2 + 1) / 5 = 60%
    assert row["c_to_t_percent"] == "60.00"


def test_apoe_summary_handles_zero_hits():
    row = build_apoe_probe_summary([{"display_name": "empty", "hits": []}])[0]
    assert row["total_exact_probe_hits"] == 0
    assert row["c_to_t_percent"] == ""  # no division by zero


# --- eToL ------------------------------------------------------------------

def test_etol_panel_sizes():
    assert etol_preset_probe_count("etol_full") == 1017
    assert etol_preset_probe_count("etol_quick") == 120
    assert etol_preset_probe_count("etol_control") == 4


def test_etol_quick_is_one_probe_per_species():
    # Quick keeps exactly one probe per species, so its size equals the species
    # count of the full panel.
    full_records = etol_preset_records("etol_full")
    species = {record["taxon"] for record in full_records}
    assert etol_preset_probe_count("etol_quick") == len(species)


def test_etol_control_probes_excluded_from_full():
    full_ids = etol_preset_query_ids("etol_full")
    control_ids = etol_preset_query_ids("etol_control")
    assert "PGK1_2" in control_ids
    assert "PGK1_2" not in full_ids


def test_etol_summary_aggregates_species():
    records = etol_preset_records("etol_full")
    sample = {
        "database_id": 1,
        "display_name": "SRX123456 brain pilot",
        "db_prefix_path": r"C:\COBLAST_data\sra\SRX123456\reads",
        "hits": [
            {"qseqid": "A_Hsalinarum_16S_1"},
            {"qseqid": "A_Hsalinarum_16S_1"},
            {"qseqid": "A_Hsalinarum_16S_3"},
            {"qseqid": "A_MethanocaldococcusSG1_16S_1"},
        ],
        "error": "",
    }
    row = build_etol_probe_summary([sample], records)[0]
    assert row["sample"] == "SRX123456"
    assert row["total_probes"] == 1017
    assert row["total_exact_probe_hits"] == 4
    assert row["probes_detected"] == 3  # two probes of one taxon + one of another
    assert row["species_detected"] == 2
    top = row["detected_species"][0]
    assert top["species"] == "Hsalinarum"
    assert top["exact_hits"] == 3


def test_etol_count_rows_include_every_probe_and_species():
    records = etol_preset_records("etol_quick")
    results = [{"display_name": "s", "hits": [{"qseqid": records[0]["probe"]}]}]
    probe_rows = etol_probe_count_rows(results, records)
    species_rows = etol_species_count_rows(results, records)
    # One row per probe and per species, including zero-count entries.
    assert len(probe_rows) == len(records)
    assert sum(row["exact_hits"] for row in probe_rows) == 1
    assert any(row["exact_hits"] == 1 for row in species_rows)
