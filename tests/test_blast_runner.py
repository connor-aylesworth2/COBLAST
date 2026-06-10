"""Unit tests for the pure (BLAST-binary-free) logic in blast_runner.

These cover FASTA validation, tabular parsing, parameter building, the
local-only command guard, and the exact-match probe overrides. None of them
invoke a BLAST+ executable, so they run in CI without NCBI BLAST+ installed.
"""

import pytest

from blast_runner import (
    BLAST_PROGRAMS,
    EXACT_MATCH_MAX_TARGET_SEQS,
    build_blast_parameters,
    coerce_to_fasta_text,
    enforce_local_blast_only,
    format_evalue,
    format_float,
    normalize_fasta_lines,
    parse_blast_tabular,
    validate_fasta_input,
)


try:  # FASTA validation needs Biopython; everything else does not.
    import Bio  # noqa: F401

    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

needs_biopython = pytest.mark.skipif(
    not HAS_BIOPYTHON, reason="Biopython is required for FASTA validation"
)


# --- FASTA coercion / normalization ---------------------------------------

def test_coerce_bare_sequence_gets_query_header():
    assert coerce_to_fasta_text("ACGT").splitlines()[0] == ">query"


def test_coerce_empty_raises():
    with pytest.raises(ValueError):
        coerce_to_fasta_text("   ")


def test_normalize_rejects_sequence_before_header():
    with pytest.raises(ValueError):
        normalize_fasta_lines("ACGT\n>missing\nACGT")


# --- validate_fasta_input --------------------------------------------------

@needs_biopython
def test_validate_nucleotide_converts_u_to_t():
    result = validate_fasta_input(">rna\nACGU", "nucleotide")
    assert "U" not in result.fasta
    assert "ACGT" in result.fasta
    assert result.total_length == 4


@needs_biopython
def test_validate_rejects_protein_letters_in_nucleotide():
    with pytest.raises(ValueError):
        validate_fasta_input(">x\nACGTE", "nucleotide")  # E is not a base


@needs_biopython
def test_validate_rejects_gap_characters():
    with pytest.raises(ValueError):
        validate_fasta_input(">x\nACG-T", "nucleotide")


@needs_biopython
def test_validate_rejects_duplicate_ids():
    with pytest.raises(ValueError):
        validate_fasta_input(">dup\nACGT\n>dup\nTTTT", "nucleotide")


@needs_biopython
def test_validate_protein_ok():
    result = validate_fasta_input(">p\nMAMAPRTEINSTRING", "protein")
    assert result.sequence_type == "protein"
    assert result.records[0].length == len("MAMAPRTEINSTRING")


# --- tabular parsing -------------------------------------------------------

def test_parse_tabular_row():
    line = "q1\ts1\tSubject one\t100.000\t64\t100\t1e-30\t120.5"
    hits = parse_blast_tabular(line)
    assert len(hits) == 1
    assert hits[0]["qseqid"] == "q1"
    assert hits[0]["sseqid"] == "s1"
    assert hits[0]["pident"] == "100.000"
    assert hits[0]["qcovs"] == "100.0"
    assert hits[0]["evalue"] == "1.00e-30"


def test_parse_tabular_wrong_column_count_raises():
    with pytest.raises(ValueError):
        parse_blast_tabular("too\tfew\tcolumns")


def test_parse_tabular_empty_is_empty():
    assert parse_blast_tabular("   \n") == []


# --- numeric formatting ----------------------------------------------------

def test_format_evalue_zero():
    assert format_evalue(0) == "0.0"


def test_format_evalue_scientific():
    assert format_evalue(1e-30) == "1.00e-30"


def test_format_float_blank():
    assert format_float("", 3) == ""


# --- local-only guard ------------------------------------------------------

def test_enforce_local_blast_only_blocks_remote():
    with pytest.raises(RuntimeError):
        enforce_local_blast_only(["blastn", "-remote", "-db", "nt"])


def test_enforce_local_blast_only_allows_clean_command():
    # Should not raise for an ordinary local command.
    enforce_local_blast_only(["blastn", "-db", "toy", "-query", "q.fasta"])


# --- parameter building ----------------------------------------------------

def test_blastn_default_task_is_megablast():
    # Regression guard: matches command-line BLAST+ so general searches are
    # functionally equivalent to running `blastn` on the CLI.
    assert BLAST_PROGRAMS["blastn"]["default_task"] == "megablast"


def test_parameters_use_preset_defaults():
    params = build_blast_parameters(
        program="blastn",
        sensitivity_preset="fast",
        evalue=None,
        max_target_seqs=None,
        word_size=None,
        perc_identity=None,
    )
    assert params["evalue"] == "10"
    assert params["max_target_seqs"] == "10"
    assert "qcov_hsp_perc" not in params
    assert "perc_identity" not in params


def test_parameters_apply_overrides():
    params = build_blast_parameters(
        program="blastn",
        sensitivity_preset="standard",
        evalue="1e-5",
        max_target_seqs="50",
        word_size="11",
        perc_identity="95",
        qcov_hsp_perc="90",
    )
    assert params["evalue"] == "1e-5"
    assert params["max_target_seqs"] == "50"
    assert params["word_size"] == "11"
    assert params["perc_identity"] == "95"
    assert params["qcov_hsp_perc"] == "90"


def test_exact_match_probe_overrides_preset_cap():
    # The 'fast' preset would cap targets at 10 and set no coverage filter; the
    # exact-match path must lift that cap and enforce full-length exact matches
    # so per-probe read counts are not silently truncated.
    params = build_blast_parameters(
        program="blastn",
        sensitivity_preset="fast",
        evalue=None,
        max_target_seqs=None,
        word_size=None,
        perc_identity=None,
        exact_match_probe=True,
    )
    assert params["perc_identity"] == "100"
    assert params["qcov_hsp_perc"] == "100"
    assert params["max_target_seqs"] == EXACT_MATCH_MAX_TARGET_SEQS
    assert int(params["max_target_seqs"]) >= 1_000_000


def test_exact_match_probe_ignores_user_max_target_seqs():
    params = build_blast_parameters(
        program="blastn",
        sensitivity_preset="standard",
        evalue=None,
        max_target_seqs="25",  # a user value must not truncate exact-match counts
        word_size=None,
        perc_identity=None,
        exact_match_probe=True,
    )
    assert params["max_target_seqs"] == EXACT_MATCH_MAX_TARGET_SEQS


def test_perc_identity_rejected_for_non_blastn():
    with pytest.raises(ValueError):
        build_blast_parameters(
            program="blastp",
            sensitivity_preset="standard",
            evalue=None,
            max_target_seqs=None,
            word_size=None,
            perc_identity="90",
        )


def test_invalid_qcov_hsp_perc_rejected():
    with pytest.raises(ValueError):
        build_blast_parameters(
            program="blastn",
            sensitivity_preset="standard",
            evalue=None,
            max_target_seqs=None,
            word_size=None,
            perc_identity=None,
            qcov_hsp_perc="150",
        )
