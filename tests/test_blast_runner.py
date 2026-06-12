"""Unit tests for the pure (BLAST-binary-free) logic in blast_runner.

These cover FASTA validation, tabular parsing, parameter building, the
local-only command guard, and the exact-match probe overrides. None of them
invoke a BLAST+ executable, so they run in CI without NCBI BLAST+ installed.
"""

import pytest

from blast_runner import (
    BLAST_PROGRAMS,
    COBLAST_NUM_THREADS_ENV,
    EXACT_MATCH_MAX_TARGET_SEQS,
    EXACT_MATCH_MT_MODE,
    NUM_THREADS_LIMIT,
    build_blast_parameters,
    coerce_to_fasta_text,
    enforce_local_blast_only,
    format_evalue,
    format_float,
    normalize_fasta_lines,
    parse_blast_tabular,
    parse_mt_mode,
    resolve_num_threads,
    run_jobs_concurrently,
    validate_fasta_input,
)
from config import allocate_batch_resources, available_cpu_count, default_thread_count


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


def test_parameters_omit_unset_fields_for_blast_defaults():
    # Sensitivity presets were removed: blank advanced fields are omitted so
    # BLAST+ applies its own defaults (e-value 10, max_target_seqs 500).
    params = build_blast_parameters(
        program="blastn",
        evalue=None,
        max_target_seqs=None,
        word_size=None,
        perc_identity=None,
    )
    assert "evalue" not in params
    assert "max_target_seqs" not in params
    assert "word_size" not in params
    assert "qcov_hsp_perc" not in params
    assert "perc_identity" not in params


def test_parameters_apply_overrides():
    params = build_blast_parameters(
        program="blastn",
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


def test_exact_match_probe_lifts_target_cap_and_enforces_full_coverage():
    # The exact-match path must enforce full-length exact matches and lift the
    # target cap so per-probe read counts are not silently truncated.
    params = build_blast_parameters(
        program="blastn",
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
            evalue=None,
            max_target_seqs=None,
            word_size=None,
            perc_identity="90",
        )


def test_invalid_qcov_hsp_perc_rejected():
    with pytest.raises(ValueError):
        build_blast_parameters(
            program="blastn",
            evalue=None,
            max_target_seqs=None,
            word_size=None,
            perc_identity=None,
            qcov_hsp_perc="150",
        )


# --- CPU threads / mt_mode -------------------------------------------------

def _params(**overrides):
    base = dict(
        program="blastn",
        evalue=None,
        max_target_seqs=None,
        word_size=None,
        perc_identity=None,
    )
    base.update(overrides)
    return build_blast_parameters(**base)


def test_config_default_thread_count_within_bounds():
    total = available_cpu_count()
    assert total >= 1
    assert 1 <= default_thread_count() <= total


def test_num_threads_added_when_requested():
    params = _params(num_threads=4)
    assert params["num_threads"] == "4"
    assert "mt_mode" not in params  # no mode requested, not an exact-match run


def test_num_threads_omitted_by_default():
    # The adaptive default is resolved in run_blast, not here.
    assert "num_threads" not in _params()


def test_exact_match_uses_auto_split_when_multithreaded():
    params = _params(num_threads=8, exact_match_probe=True)
    assert params["mt_mode"] == EXACT_MATCH_MT_MODE
    assert params["mt_mode"] == "0"  # BLAST chooses query or database splitting.


def test_mt_mode_omitted_when_single_threaded():
    # mt_mode is meaningless with a single thread, so it should not be emitted.
    assert "mt_mode" not in _params(num_threads=1, exact_match_probe=True)


def test_parse_mt_mode_rejects_bad_value():
    with pytest.raises(ValueError):
        parse_mt_mode("9")


def test_resolve_num_threads_explicit_request_wins(monkeypatch):
    monkeypatch.delenv(COBLAST_NUM_THREADS_ENV, raising=False)
    assert resolve_num_threads(3) == 3


def test_resolve_num_threads_env_override(monkeypatch):
    monkeypatch.setenv(COBLAST_NUM_THREADS_ENV, "2")
    assert resolve_num_threads(None) == 2


def test_resolve_num_threads_falls_back_to_default(monkeypatch):
    monkeypatch.delenv(COBLAST_NUM_THREADS_ENV, raising=False)
    assert resolve_num_threads(None) == default_thread_count()


def test_resolve_num_threads_rejects_out_of_range(monkeypatch):
    monkeypatch.delenv(COBLAST_NUM_THREADS_ENV, raising=False)
    with pytest.raises(ValueError):
        resolve_num_threads(NUM_THREADS_LIMIT + 1)


def test_run_jobs_concurrently_preserves_order_and_captures_errors():
    def work(x):
        if x == 2:
            raise ValueError("boom")
        return x * 10

    jobs = [{"x": 0}, {"x": 1}, {"x": 2}, {"x": 3}]
    results = run_jobs_concurrently(work, jobs, max_workers=2)
    assert results[0] == 0
    assert results[1] == 10
    assert isinstance(results[2], ValueError)  # captured in place, not raised
    assert results[3] == 30


# --- batch core-budget allocation -----------------------------------------

def test_allocate_batch_resources_never_oversubscribes(monkeypatch):
    monkeypatch.delenv("COBLAST_BATCH_WORKERS", raising=False)
    budget = default_thread_count()
    for jobs in (1, 2, 5, budget, budget + 7):
        workers, threads = allocate_batch_resources(jobs)
        assert workers >= 1 and threads >= 1
        assert workers <= jobs
        assert workers * threads <= budget  # no CPU oversubscription


def test_allocate_batch_resources_prefers_concurrency_when_many_jobs(monkeypatch):
    monkeypatch.delenv("COBLAST_BATCH_WORKERS", raising=False)
    budget = default_thread_count()
    workers, threads = allocate_batch_resources(budget + 10)
    assert workers == budget  # spend the budget on workers
    assert threads == 1


def test_allocate_batch_resources_single_job_uses_threads(monkeypatch):
    monkeypatch.delenv("COBLAST_BATCH_WORKERS", raising=False)
    budget = default_thread_count()
    workers, threads = allocate_batch_resources(1)
    assert workers == 1
    assert threads == budget  # one database: hand the spare cores to threads


def test_allocate_batch_resources_env_override(monkeypatch):
    monkeypatch.setenv("COBLAST_BATCH_WORKERS", "1")
    workers, _ = allocate_batch_resources(8)
    assert workers == 1


def test_allocate_batch_resources_requested_workers_wins(monkeypatch):
    monkeypatch.delenv("COBLAST_BATCH_WORKERS", raising=False)
    workers, _ = allocate_batch_resources(8, requested_workers=2)
    assert workers == 2


def test_allocate_batch_resources_requested_workers_clamped_to_jobs(monkeypatch):
    monkeypatch.delenv("COBLAST_BATCH_WORKERS", raising=False)
    workers, _ = allocate_batch_resources(3, requested_workers=99)
    assert workers == 3


def test_allocate_batch_resources_request_beats_env(monkeypatch):
    monkeypatch.setenv("COBLAST_BATCH_WORKERS", "1")
    workers, _ = allocate_batch_resources(8, requested_workers=4)
    assert workers == 4


def test_allocate_batch_resources_invalid_request_falls_back(monkeypatch):
    monkeypatch.delenv("COBLAST_BATCH_WORKERS", raising=False)
    workers, _ = allocate_batch_resources(8, requested_workers="abc")
    assert workers == min(default_thread_count(), 8)
