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
    assert len(etol_control_query_ids()) == 4


def test_etol_quick_is_one_probe_per_species():
    # Quick keeps exactly one probe per species, so its size equals the species
    # count of the full panel.
    full_records = etol_preset_records("etol_full")
    species = {record["taxon"] for record in full_records}
    assert etol_preset_probe_count("etol_quick") == len(species)


def test_etol_control_probes_excluded_from_full():
    full_ids = etol_preset_query_ids("etol_full")
    control_ids = etol_control_query_ids()
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


# --- net E-value gate and control-read counting ----------------------------

def test_filter_net_probe_hits_gates_on_evalue_and_panel():
    from app import filter_net_probe_hits

    panel = {"A_x_16S_1", "B_y_16S_1"}
    hits = [
        {"qseqid": "A_x_16S_1", "sseqid": "r1", "evalue": "1e-30"},    # kept
        {"qseqid": "A_x_16S_1", "sseqid": "r2", "evalue": "0.5"},      # E >= 0.01
        {"qseqid": "A_x_16S_1", "sseqid": "r3", "evalue": "0.01"},     # not < 0.01
        {"qseqid": "Z_off_16S_1", "sseqid": "r4", "evalue": "1e-50"},  # off-panel
        {"qseqid": "B_y_16S_1", "sseqid": "r5"},                       # unparseable E
    ]
    assert [hit["sseqid"] for hit in filter_net_probe_hits(hits, panel)] == ["r1"]


def test_count_control_reads_dedups_to_best_control_probe():
    from app import count_control_reads

    control_ids = frozenset({"PGK1_2", "PGK1_3", "hNSE_2"})
    # read1 hits two PGK1 probes -> counted once, against the higher-bitscore one.
    hits = [
        {"qseqid": "PGK1_2", "sseqid": "read1", "bitscore": "50", "pident": "95", "qcovs": "90"},
        {"qseqid": "PGK1_3", "sseqid": "read1", "bitscore": "120", "pident": "100", "qcovs": "100"},
        {"qseqid": "hNSE_2", "sseqid": "read2", "bitscore": "80", "pident": "98", "qcovs": "100"},
    ]
    # Every control probe is present (incl. zeros); read1 is not double-counted.
    assert count_control_reads(hits, control_ids) == {"PGK1_2": 0, "PGK1_3": 1, "hNSE_2": 1}


