"""Unit tests for the APOE and eToL probe-count summaries.

These exercise the pure aggregation logic over synthetic hit lists, so they do
not require BLAST+ or any patient data. The eToL tests read the bundled probe
panels under data/, so run pytest from the repository root.
"""

from apoe_summary import build_apoe_probe_summary
from etol_summary import (
    build_etol_probe_summary,
    compute_host_cells,
    etol_control_query_ids,
    etol_preset_fasta,
    etol_preset_probe_count,
    etol_preset_query_ids,
    etol_preset_records,
    etol_probe_count_rows,
    etol_search_query_ids,
    etol_species_count_rows,
    normalized_abundance,
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
    assert row["sample_database"] == "SRX123456"  # accession pulled from the path
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


# --- host-cell normalization ----------------------------------------------

def test_etol_microbial_search_appends_control_probes():
    # The microbial panels are searched together with the 4 control probes so a
    # single run yields both the net and the host-normalization counts.
    full_ids = etol_search_query_ids("etol_full")
    control_ids = etol_control_query_ids()
    assert control_ids <= full_ids
    assert len(full_ids) == 1017 + len(control_ids)
    # The standalone control preset is searched as-is (no microbial appended).
    assert etol_search_query_ids("etol_control") == control_ids


def test_compute_host_cells_means_genes_over_fifty():
    # mean(PGK1)=102, mean(hNSE)=98 -> combined 100 -> /50 = 2 host cells.
    host_cells = compute_host_cells(
        {"PGK1_2": 100, "PGK1_3": 104, "hNSE_2": 96, "hNSE_3": 100}
    )
    assert host_cells == 2.0
    assert compute_host_cells({}) == 0.0  # no control reads -> undefined
    assert normalized_abundance(10, host_cells) == 5.0
    assert normalized_abundance(10, 0.0) is None


def test_etol_summary_reports_normalized_abundance():
    records = etol_preset_records("etol_full")
    sample = {
        "display_name": "SRX999 brain",
        "hits": [{"qseqid": "A_Hsalinarum_16S_1"}] * 20,  # 20 microbial reads
        "etol_control_counts": {"PGK1_2": 100, "PGK1_3": 104, "hNSE_2": 96, "hNSE_3": 100},
        "error": "",
    }
    row = build_etol_probe_summary([sample], records)[0]
    assert row["host_cells"] == 2.0
    assert row["normalized"] is True
    top = row["detected_species"][0]
    assert top["exact_hits"] == 20
    assert top["normalized_abundance"] == 10.0  # 20 reads / 2 host cells


def test_etol_summary_without_control_counts_is_unnormalized():
    records = etol_preset_records("etol_full")
    sample = {"display_name": "s", "hits": [{"qseqid": "A_Hsalinarum_16S_1"}], "error": ""}
    row = build_etol_probe_summary([sample], records)[0]
    assert row["host_cells"] == 0.0
    assert row["normalized"] is False
    assert row["detected_species"][0]["normalized_abundance"] is None


# --- cross-probe read de-duplication ---------------------------------------

def test_deduplicate_reads_allocates_each_read_to_best_probe():
    from app import deduplicate_reads_to_best_probe

    hits = [
        {"qseqid": "A_x_16S_1", "sseqid": "read1", "bitscore": "50.0", "pident": "95", "qcovs": "90"},
        {"qseqid": "B_y_16S_1", "sseqid": "read1", "bitscore": "120.0", "pident": "100", "qcovs": "100"},
        {"qseqid": "A_x_16S_1", "sseqid": "read2", "bitscore": "80.0", "pident": "98", "qcovs": "100"},
    ]
    kept, removed = deduplicate_reads_to_best_probe(hits)
    assert removed == 1
    kept_pairs = {(hit["qseqid"], hit["sseqid"]) for hit in kept}
    # read1 goes to its higher-scoring probe; read2 is unique.
    assert ("B_y_16S_1", "read1") in kept_pairs
    assert ("A_x_16S_1", "read2") in kept_pairs
    assert ("A_x_16S_1", "read1") not in kept_pairs


def test_deduplicate_reads_is_order_independent():
    from app import deduplicate_reads_to_best_probe

    forward = [
        {"qseqid": "A_x_16S_1", "sseqid": "r", "bitscore": "120.0", "pident": "100", "qcovs": "100"},
        {"qseqid": "B_y_16S_1", "sseqid": "r", "bitscore": "60.0", "pident": "95", "qcovs": "90"},
    ]
    reverse = list(reversed(forward))
    assert deduplicate_reads_to_best_probe(forward)[0][0]["qseqid"] == "A_x_16S_1"
    assert deduplicate_reads_to_best_probe(reverse)[0][0]["qseqid"] == "A_x_16S_1"


# --- megablast / blastn-short partition of the eToL panel ------------------

def test_etol_full_panel_has_exactly_one_non_megablast_probe():
    # Only F3_Gpolymorpha_18S_7 lacks a 28-base unambiguous window, so it is the
    # single probe that falls back to blastn-short while the other 1,016 use
    # megablast. Guards the speedup split against future probe-panel edits.
    from blast_runner import _parse_panel_fasta, has_megablast_seed

    pairs = _parse_panel_fasta(etol_preset_fasta("etol_full"))
    short = [header for header, seq in pairs if not has_megablast_seed(seq)]
    assert short == ["F3_Gpolymorpha_18S_7"]


def test_etol_quick_panel_is_all_megablast_safe():
    from blast_runner import _parse_panel_fasta, has_megablast_seed

    pairs = _parse_panel_fasta(etol_preset_fasta("etol_quick"))
    assert all(has_megablast_seed(seq) for _, seq in pairs)
