"""Unit tests for the design-matrix parser (eToL heatmap condition labels).

These exercise the strict CSV/TSV format and the accession/display-name index it
produces. No BLAST+ or patient data required; run pytest from the repo root.
"""

import pytest

from design_matrix import DesignMatrixError, parse_design_matrix


def test_parses_csv_with_accession_and_name_keys():
    text = (
        "sample,condition\n"
        "SRR21676105,CONTROL\n"
        "SRA SRX17674465 reads,AD/LBD\n"
    )
    index = parse_design_matrix(text, filename="design.csv")
    # Accession-bearing ids land in by_accession (upper-cased); every row is also
    # indexed by its lower-cased verbatim id.
    assert index["by_accession"] == {
        "SRR21676105": "CONTROL",
        "SRX17674465": "AD/LBD",
    }
    assert index["by_name"]["srr21676105"] == "CONTROL"
    assert index["by_name"]["sra srx17674465 reads"] == "AD/LBD"
    assert index["conditions"] == ["CONTROL", "AD/LBD"]
    assert index["row_count"] == 2
    assert index["warnings"] == []


def test_header_columns_are_case_and_order_insensitive():
    text = "Condition,Sample\nAD,SRR1\n"
    index = parse_design_matrix(text, filename="d.csv")
    assert index["by_accession"] == {"SRR1": "AD"}


def test_tab_separated_by_extension():
    text = "sample\tcondition\nSRR1\tAD\n"
    index = parse_design_matrix(text, filename="d.tsv")
    assert index["by_accession"] == {"SRR1": "AD"}


def test_tolerates_utf8_bom():
    text = "﻿sample,condition\nSRR1,AD\n"
    index = parse_design_matrix(text, filename="d.csv")
    assert index["row_count"] == 1
    assert index["by_accession"] == {"SRR1": "AD"}


def test_extra_columns_ignored():
    text = "sample,condition,sex,region\nSRR1,AD,female,HPC\n"
    index = parse_design_matrix(text, filename="d.csv")
    assert index["by_accession"] == {"SRR1": "AD"}


def test_blank_condition_warns_but_is_not_fatal():
    text = "sample,condition\nSRR1,\n"
    index = parse_design_matrix(text, filename="d.csv")
    assert index["row_count"] == 1
    assert index["by_name"]["srr1"] == ""
    assert index["conditions"] == []
    assert any("no condition" in w for w in index["warnings"])


def test_duplicate_sample_is_fatal():
    text = "sample,condition\nSRR1,AD\nSRR1,CONTROL\n"
    with pytest.raises(DesignMatrixError, match="Duplicate sample"):
        parse_design_matrix(text, filename="d.csv")


def test_missing_required_column_is_fatal():
    text = "foo,bar\n1,2\n"
    with pytest.raises(DesignMatrixError, match="must have a 'sample' column"):
        parse_design_matrix(text, filename="d.csv")


def test_empty_file_is_fatal():
    with pytest.raises(DesignMatrixError, match="empty"):
        parse_design_matrix("", filename="d.csv")


def test_header_only_is_fatal():
    with pytest.raises(DesignMatrixError, match="no data rows"):
        parse_design_matrix("sample,condition\n", filename="d.csv")
