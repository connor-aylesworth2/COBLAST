"""Unit tests for the APOE and eToL probe-count summaries.

These exercise the pure aggregation logic over synthetic hit lists, so they do
not require BLAST+ or any patient data. The eToL tests read the bundled probe
panels under data/, so run pytest from the repository root.
"""

from apoe_summary import build_apoe_probe_summary
from design_matrix import parse_design_matrix
from etol_summary import (
    _sample_condition,
    build_etol_matrix,
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
    sort_results_by_condition,
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


# --- eToL-V (viral) preset -------------------------------------------------

def test_etol_v_panel_size_and_own_controls():
    # 115 viral probes + this panel's own 2 PGK controls (not the cellular 4).
    assert etol_preset_probe_count("etol_v") == 115
    v_controls = etol_control_query_ids("etol_v")
    assert v_controls == {"PGK1_1", "PGK1_2"}
    # The viral panel is searched together with its own controls (117 total) and
    # must NOT pull in the cellular hNSE controls.
    search_ids = etol_search_query_ids("etol_v")
    assert v_controls <= search_ids
    assert len(search_ids) == 115 + 2
    assert not any(qid.startswith("hNSE") for qid in search_ids)
    # Controls are not searched as a viral taxon.
    assert "PGK1_1" not in etol_preset_query_ids("etol_v")


def test_etol_v_records_map_to_viral_domain_and_families():
    records = etol_preset_records("etol_v")
    assert {r["domain"] for r in records} == {"Viruses"}
    groups = {r["group"] for r in records}
    assert groups == {"V-HHV", "V-HAdV", "V-HPV", "V-HCoV"}
    # Subunit is retained in the taxon so AdC penton is its own row (matches the
    # eToL-V dissertation's reporting granularity).
    by_probe = {r["probe"]: r for r in records}
    penton = by_probe["V-HAdV_AdC_penton_1"]
    assert penton["taxon"] == "V-HAdV_AdC_penton"
    assert penton["species"] == "AdC_penton"


def test_etol_v_summary_normalizes_on_own_pgk():
    records = etol_preset_records("etol_v")
    sample = {
        "display_name": "SRX555 brain virome",
        "hits": [{"qseqid": "V-HAdV_AdC_penton_1"}] * 12,
        # mean(PGK1)=300 -> /50 = 6 host cells; 12 reads / 6 = 2 per host cell.
        "etol_control_counts": {"PGK1_1": 280, "PGK1_2": 320},
        "error": "",
    }
    row = build_etol_probe_summary([sample], records)[0]
    assert row["host_cells"] == 6.0
    top = row["detected_species"][0]
    assert top["species"] == "AdC_penton"
    assert top["exact_hits"] == 12
    assert top["normalized_abundance"] == 2.0


# --- heatmap matrix --------------------------------------------------------

def test_build_etol_matrix_shape_and_counts():
    records = etol_preset_records("etol_quick")
    results = [
        {
            "display_name": "SRX111_AD",
            "db_prefix_path": r"C:\COBLAST_data\sra\SRX111\reads",
            "hits": [{"qseqid": records[0]["probe"]}] * 3,
            "etol_control_counts": {"PGK1_2": 100, "hNSE_2": 100},
        },
        {
            "display_name": "SRX222_CTRL",
            "db_prefix_path": r"C:\COBLAST_data\sra\SRX222\reads",
            "hits": [],
        },
    ]
    matrix = build_etol_matrix(results, records, level="species")
    # One row per species/taxon, one column per sample, hits aligned row x col.
    assert len(matrix["rows"]) == len({r["taxon"] for r in records})
    assert [c["sample"] for c in matrix["cols"]] == ["SRX111", "SRX222"]
    assert matrix["cols"][0]["condition"] == "AD"
    assert matrix["cols"][1]["condition"] == "CTRL"
    # mean(PGK1_2=100, hNSE_2=100)=100 -> /50 = 2 host cells.
    assert matrix["cols"][0]["host_cells"] == 2.0
    total = sum(value for row in matrix["hits"] for value in row)
    assert total == 3
    # No contig identification ran, so there is no validated layer.
    assert matrix["confirmed"] is None


def test_build_etol_matrix_labels_run_accessions():
    # SRR/ERR/DRR run accessions (not just SRX experiment accessions) should be
    # picked up so SRR-built patient databases get a clean column label.
    records = etol_preset_records("etol_quick")
    results = [
        {
            "display_name": "SRR9999999 AD brain",
            "db_prefix_path": r"C:\COBLAST_data\sra\SRR9999999\reads",
            "hits": [],
        }
    ]
    matrix = build_etol_matrix(results, records, level="species")
    assert matrix["cols"][0]["sample"] == "SRR9999999"
    assert matrix["cols"][0]["condition"] == "AD"


def test_build_etol_matrix_probe_level_and_confirmed_layer():
    records = etol_preset_records("etol_quick")
    taxon = records[0]["taxon"]
    results = [
        {
            "display_name": "SRX999_AD/LBD",
            "hits": [{"qseqid": records[0]["probe"]}] * 5,
            "contig_identification": {taxon: {"confirmed_reads": 4}},
        }
    ]
    probe_matrix = build_etol_matrix(results, records, level="probe")
    assert len(probe_matrix["rows"]) == len(records)
    assert probe_matrix["rows"][0]["key"] == records[0]["probe"]
    # Combined diagnosis suffix is preferred over a bare "AD".
    assert probe_matrix["cols"][0]["condition"] == "AD/LBD"
    # Contig identification populates the validated (confirmed) layer per taxon.
    assert probe_matrix["confirmed"] is not None
    assert probe_matrix["confirmed"][0][0] == 4


def test_build_etol_matrix_condition_index_overrides_regex():
    # An uploaded design matrix is authoritative for the condition strip, fixing
    # the name-guessing regex's failure on auto-generated "SRA <acc> reads" names
    # (every one of which the regex mislabels "AD" via the "ad" in "reads").
    records = etol_preset_records("etol_quick")
    results = [
        {
            "display_name": "SRA SRR21676105 reads",
            "db_prefix_path": r"C:\Users\x\Downloads\db\srr21676105",
            "hits": [],
        },
        {"display_name": "SRA SRR21676101 reads", "hits": []},
        {"display_name": "SRA SRR99999999 reads", "hits": []},  # absent from matrix
    ]
    index = parse_design_matrix(
        "sample,condition\nSRR21676105,CONTROL\nSRR21676101,AD/LBD\n",
        filename="design.csv",
    )
    matrix = build_etol_matrix(results, records, condition_index=index)
    assert [c["condition"] for c in matrix["cols"]] == ["CONTROL", "AD/LBD", ""]
    # The regex path would have mislabeled the control as "AD"; the matrix doesn't.
    assert _sample_condition(results[0]) == "AD"
    # Unmatched columns are reported (by their resolved accession label) not guessed.
    assert matrix["unmatched_samples"] == ["SRR99999999"]


def test_build_etol_matrix_condition_index_matches_by_display_name():
    # Custom (non-SRA) databases have no accession, so the matrix matches on the
    # database display name instead.
    records = etol_preset_records("etol_quick")
    results = [{"display_name": "Patient cortex pool", "hits": []}]
    index = parse_design_matrix(
        "sample,condition\nPatient cortex pool,AD\n", filename="design.csv"
    )
    matrix = build_etol_matrix(results, records, condition_index=index)
    assert matrix["cols"][0]["condition"] == "AD"
    assert matrix["unmatched_samples"] == []


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


# --- Batch summary export --------------------------------------------------

def _delimited_to_dict(text, delimiter):
    """Parse a two-column Statistic/Value export back into a dict."""
    import csv
    import io

    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    assert rows[0] == ["Statistic", "Value"]
    return {label: value for label, value in rows[1:]}


def test_batch_summary_export_mirrors_panel():
    from result_store import batch_summary_rows_as_delimited

    batch = {
        "program": "blastn",
        "database_results": [{"display_name": "SRR1"}, {"display_name": "SRR2"}],
        "total_runtime_seconds": 1032.318,
        "wall_clock_seconds": 5037.698,
        "batch_workers": 30,
        "query_count": 1021,
        "query_total_length": 65344,
        "total_hits": 23822,
        "etol_probe_preset": True,
        "etol_preset_label": "eToL Full",
        "hit_filter": "default megablast net, E-value < 0.01",
        "etol_dedup_removed": 253909,
        "etol_normalized": True,
        "assemble_contigs": True,
        "contig_count": 3041,
        "identify_contigs": True,
        "contigs_identified": 3019,
        "species_id_db": "ToL_rRNA",
        "human_filter_enabled": True,
        "human_filter_hits_removed": 166339,
        "human_filter_db": "nt_human_9606",
    }
    parsed = _delimited_to_dict(batch_summary_rows_as_delimited(batch, delimiter=","), ",")
    assert parsed["Program"] == "blastn"
    assert parsed["Databases"] == "2"
    # The two timers stay distinct and both accurate in the export.
    assert parsed["Total runtime (BLAST search, summed across databases, seconds)"] == "1032.318"
    assert parsed["Wall-clock elapsed time (seconds)"] == "5037.698"
    assert parsed["Concurrency (databases at a time)"] == "30"
    assert parsed["Total hits"] == "23822"
    assert parsed["Probe preset"] == "eToL Full"
    assert parsed["Contigs assembled (CAP3)"] == "3041"
    assert parsed["Human filter hits removed"] == "166339"
    # TSV round-trips identically with a tab delimiter.
    parsed_tsv = _delimited_to_dict(batch_summary_rows_as_delimited(batch, delimiter="\t"), "\t")
    assert parsed_tsv["Wall-clock elapsed time (seconds)"] == "5037.698"


def test_batch_summary_export_omits_absent_rows():
    from result_store import batch_summary_rows_as_delimited

    # A plain (non-preset) batch with no wall clock or concurrency recorded:
    # those rows are dropped, mirroring the panel's conditionals.
    batch = {
        "program": "blastn",
        "database_results": [{"display_name": "db1"}],
        "total_runtime_seconds": 4.0,
        "query_count": 1,
        "query_total_length": 300,
        "total_hits": 5,
    }
    parsed = _delimited_to_dict(batch_summary_rows_as_delimited(batch, delimiter=","), ",")
    assert "Wall-clock elapsed time (seconds)" not in parsed
    assert "Concurrency (databases at a time)" not in parsed
    assert "Probe preset" not in parsed
    assert "Human filter hits removed" not in parsed
    assert parsed["Total hits"] == "5"


# --- sample ordering by design-matrix condition ----------------------------

def test_sort_results_by_condition_groups_by_design_matrix_order():
    # Condition order follows the matrix's first-seen order (Control before AD),
    # samples sort alphabetically within a condition, and a sample with no matrix
    # row falls to the end rather than being scattered.
    results = [
        {"display_name": "SRA SRR3 reads", "db_prefix_path": ""},
        {"display_name": "SRA SRR1 reads", "db_prefix_path": ""},
        {"display_name": "SRA SRR2 reads", "db_prefix_path": ""},
        {"display_name": "SRA SRR9 reads", "db_prefix_path": ""},  # absent from matrix
    ]
    index = parse_design_matrix(
        "sample,condition\nSRR1,Control\nSRR3,Control\nSRR2,AD\n", filename="d.csv"
    )
    ordered = [r["display_name"] for r in sort_results_by_condition(results, index)]
    assert ordered == [
        "SRA SRR1 reads",
        "SRA SRR3 reads",
        "SRA SRR2 reads",
        "SRA SRR9 reads",
    ]
    # No matrix -> original order untouched (never sort on the unreliable regex).
    assert sort_results_by_condition(results, None) is results
