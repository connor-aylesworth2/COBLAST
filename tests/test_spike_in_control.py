"""Synthetic spike-in positive control for the COBLAST+ eToL net pipeline.

WHAT THIS IS (for the dissertation's verification section)
----------------------------------------------------------
A spike-in control is the microbiology idea of a *positive control* applied to
software: you construct a sample whose true composition you already know, run it
through the instrument untouched, and check the instrument returns the known
answer. Here the "instrument" is the real COBLAST+ eToL net path -- the same
functions the web app runs, not a reimplementation:

    run_blast_probe_panel  (megablast net over the whole probe panel)
      -> filter_net_probe_hits      (E-value < 0.01 net gate)
      -> split control vs microbial hits
      -> count_control_reads        (host-normalization denominator)
      -> deduplicate_reads_to_best_probe  (one read -> its best probe)
      -> build_etol_probe_summary   (species + /50 host-cell normalization)

The synthetic sample is a FASTA of reads with a KNOWN composition:
  * three microbial species from the panel, spiked at known read counts that
    straddle the paper's 3-5 reads/host-cell cellular cutoff (20, 10, 3);
  * the four housekeeping control probes (PGK1, hNSE) spiked so the host-cell
    denominator is an exact round number (100 reads each -> 2 host cells);
  * one panel species deliberately left ABSENT (0 reads);
  * a few hundred random-sequence reads as background noise.

Each spiked read embeds a real 64 bp panel probe inside random flanks, so the
net must actually recover it by BLAST -- the marker is really in the read, the
abundance is exactly what we put there.

WHAT IT VERIFIES (three claims a reader can check)
--------------------------------------------------
  1. Recovery / sensitivity - every spiked species is detected and its net read
     count equals what was spiked (abundance is preserved end to end).
  2. Specificity - the absent species is NOT reported, and the noise reads
     produce no spurious species; exactly the three spiked species come back.
     (Conserved-region cross-hits between probes are absorbed by the
     de-duplication step, which is itself part of what this exercises.)
  3. Normalization - reported normalized abundance == net reads / host cells,
     with host cells computed from the spiked control reads (here 2.0).

This is verification, not biological validation: it proves the pipeline reports
what is present at the abundance present and rejects what is absent. Validation
against real data (reproducing Lathe/Veso) is a separate, later step.

RUN IT
------
  * As a test:            pytest tests/test_spike_in_control.py
  * As a dissertation
    artifact (prints an
    expected-vs-observed
    table):               python tests/test_spike_in_control.py

Both need NCBI BLAST+ (blastn + makeblastdb) discoverable the same way the app
finds it (PATH, BLAST_BIN, or the bundled copy). Without it the pytest test
skips; the standalone run prints a setup hint and exits non-zero.
"""

from __future__ import annotations

import platform
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# Allow `python tests/test_spike_in_control.py` from the repo root as well as
# pytest (which puts the root on the path via pytest.ini).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import blast_exe, default_thread_count  # noqa: E402
from blast_runner import run_blast_probe_panel  # noqa: E402
from etol_summary import (  # noqa: E402
    build_etol_probe_summary,
    etol_control_query_ids,
    etol_preset_records,
    etol_search_fasta,
    etol_search_pairs,
    etol_search_query_ids,
)

PRESET = "etol_quick"  # one probe per species: fast and unambiguous to reason about
READ_LEN = 150         # a realistic Illumina-ish read; probe (64 bp) sits in the middle
SPIKE_LEVELS = (20, 10, 3)   # straddle the 3-5 reads/host-cell cellular cutoff
CONTROL_READS = 100    # per control probe -> host cells = mean(100,100)/50 = 2.0
NOISE_READS = 500      # random-sequence background: must produce zero detections
EXPECTED_HOST_CELLS = 2.0
RNG_SEED = 20240601    # fixed so the synthetic sample is byte-for-byte reproducible

_PURE_ACGT = re.compile(r"[ACGT]+")


def _probe_sequences() -> dict[str, str]:
    """Map probe id -> sequence for every probe actually BLASTed (panel + controls)."""
    return {header: seq for header, seq in etol_search_pairs(PRESET)}


def choose_spike_probes() -> tuple[list[dict[str, str]], dict[str, str]]:
    """Pick 3 'present' + 1 'absent' panel species, preferring distinct domains.

    Selection is deterministic (panel order is fixed) and restricted to pure-ACGT
    probes so megablast can always seed on the 64 bp marker. Preferring different
    domains keeps the spiked markers maximally divergent, so the specificity claim
    does not lean on the de-duplication step -- though that step would cover
    conserved-region cross-hits anyway. Returns (four probe records, id->seq).
    """
    seqs = _probe_sequences()
    controls = etol_control_query_ids(PRESET)
    chosen: list[dict[str, str]] = []
    seen_domains: set[str] = set()
    seen_taxa: set[str] = set()
    # First pass: one probe per new domain. Second pass: top up with new taxa.
    for prefer_new_domain in (True, False):
        for record in etol_preset_records(PRESET):
            probe = record["probe"]
            if probe in controls or record["taxon"] in seen_taxa:
                continue
            if not _PURE_ACGT.fullmatch(seqs.get(probe, "")):
                continue
            if prefer_new_domain and record["domain"] in seen_domains:
                continue
            chosen.append(record)
            seen_taxa.add(record["taxon"])
            seen_domains.add(record["domain"])
            if len(chosen) == 4:
                return chosen, seqs
    raise RuntimeError("Panel does not expose four seedable species to spike.")


def _random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _spiked_read(rng: random.Random, probe_seq: str) -> str:
    """Embed a probe inside random flanks so the read is a realistic length."""
    pad = READ_LEN - len(probe_seq)
    left = pad // 2
    return _random_seq(rng, left) + probe_seq + _random_seq(rng, pad - left)


def synthesize_reads(spike_probes: list[dict[str, str]], seqs: dict[str, str]) -> tuple[str, dict[str, int]]:
    """Build the spiked read FASTA and return (fasta_text, expected net counts).

    ``expected`` maps each present species' taxon to the number of reads spiked
    for it -- the ground truth the pipeline must recover.
    """
    rng = random.Random(RNG_SEED)
    present = spike_probes[:3]
    lines: list[str] = []
    expected: dict[str, int] = {}
    read_no = 0

    # Microbial spikes: known species at known abundance.
    for record, level in zip(present, SPIKE_LEVELS):
        expected[record["taxon"]] = level
        for _ in range(level):
            lines.append(f">read_{read_no}_{record['probe']}\n{_spiked_read(rng, seqs[record['probe']])}")
            read_no += 1

    # Host-cell controls: spike every control probe equally for a clean denominator.
    for probe in sorted(etol_control_query_ids(PRESET)):
        for _ in range(CONTROL_READS):
            lines.append(f">read_{read_no}_{probe}\n{_spiked_read(rng, seqs[probe])}")
            read_no += 1

    # Background noise: random reads that must match nothing.
    for _ in range(NOISE_READS):
        lines.append(f">read_{read_no}_noise\n{_random_seq(rng, READ_LEN)}")
        read_no += 1

    return "\n".join(lines) + "\n", expected


def _build_db(fasta_text: str, workdir: Path) -> str:
    """Write the reads and index them with makeblastdb (the reads are the DB)."""
    fasta = workdir / "spike_in_reads.fasta"
    fasta.write_text(fasta_text, encoding="utf-8")
    prefix = workdir / "spike_in_db"
    subprocess.run(
        [str(blast_exe("makeblastdb")), "-in", str(fasta), "-dbtype", "nucl",
         "-out", str(prefix)],
        check=True, capture_output=True, text=True,
    )
    return str(prefix)


def run_net_pipeline(db_prefix: str) -> dict[str, object]:
    """Drive the real eToL net path over the spiked DB and summarize per species.

    Mirrors app.py's batch route: net BLAST -> E-value gate -> split controls ->
    count controls -> de-duplicate microbial reads -> species/normalization
    summary. Imported from app inside the function to match the test-suite
    convention and avoid importing the Flask app at collection time.
    """
    from app import (
        count_control_reads,
        deduplicate_reads_to_best_probe,
        filter_net_probe_hits,
    )

    probe_query_ids = set(etol_search_query_ids(PRESET))
    control_query_ids = etol_control_query_ids(PRESET)

    result = run_blast_probe_panel(panel_fasta=etol_search_fasta(PRESET), database=db_prefix)
    panel_hits = filter_net_probe_hits(result.hits, probe_query_ids)
    control_hits = [h for h in panel_hits if h.get("qseqid", "") in control_query_ids]
    micro_hits = [h for h in panel_hits if h.get("qseqid", "") not in control_query_ids]
    control_counts = count_control_reads(control_hits, control_query_ids)
    micro_hits, _removed = deduplicate_reads_to_best_probe(micro_hits)

    db_result = {
        "display_name": "SPIKE_IN_CONTROL",
        "hits": micro_hits,
        "etol_control_counts": control_counts,
    }
    return build_etol_probe_summary([db_result], etol_preset_records(PRESET))[0]


def run_spike_in() -> dict[str, object]:
    """Build the synthetic sample, run the pipeline, and return everything needed
    to both assert and print the expected-vs-observed table."""
    spike_probes, seqs = choose_spike_probes()
    present = spike_probes[:3]
    absent = spike_probes[3]
    fasta_text, expected = synthesize_reads(spike_probes, seqs)

    with tempfile.TemporaryDirectory(prefix="coblast_spikein_") as tmp:
        db_prefix = _build_db(fasta_text, Path(tmp))
        summary = run_net_pipeline(db_prefix)

    detected = {d["taxon"]: d for d in summary["detected_species"]}
    return {
        "present": present,
        "absent": absent,
        "expected": expected,
        "summary": summary,
        "detected": detected,
    }


def check(outcome: dict[str, object]) -> None:
    """Assert the three claims: recovery, specificity, normalization."""
    present = outcome["present"]
    absent = outcome["absent"]
    expected: dict[str, int] = outcome["expected"]
    summary = outcome["summary"]
    detected: dict[str, dict] = outcome["detected"]

    host_cells = summary["host_cells"]

    # (3) Normalization denominator recovered exactly from the spiked controls.
    assert host_cells == EXPECTED_HOST_CELLS, f"host cells {host_cells} != {EXPECTED_HOST_CELLS}"

    # (2) Specificity: exactly the three spiked species, nothing else.
    assert summary["species_detected"] == 3, summary["species_detected"]
    assert absent["taxon"] not in detected, f"absent species {absent['taxon']} was reported"

    for record in present:
        taxon = record["taxon"]
        # (1) Recovery: detected, and net count equals what was spiked.
        assert taxon in detected, f"spiked species {taxon} was NOT detected"
        got = detected[taxon]["exact_hits"]
        assert got == expected[taxon], f"{taxon}: recovered {got} reads, spiked {expected[taxon]}"
        # (3) Normalization: reported abundance == net reads / host cells.
        expect_norm = expected[taxon] / host_cells
        assert detected[taxon]["normalized_abundance"] == expect_norm, (
            f"{taxon}: normalized {detected[taxon]['normalized_abundance']} != {expect_norm}"
        )


def format_table(outcome: dict[str, object]) -> str:
    """Render the expected-vs-observed table for the dissertation write-up."""
    present = outcome["present"]
    absent = outcome["absent"]
    expected: dict[str, int] = outcome["expected"]
    summary = outcome["summary"]
    detected: dict[str, dict] = outcome["detected"]
    host_cells = summary["host_cells"]

    lines = [
        "COBLAST+ eToL net -- synthetic spike-in positive control",
        f"panel = {PRESET}   host cells (from control reads) = {host_cells}   "
        f"(expected {EXPECTED_HOST_CELLS})",
        "",
        f"{'species (taxon)':<34}{'spiked':>7}{'recovered':>11}{'norm exp':>10}{'norm obs':>10}  result",
        "-" * 86,
    ]
    for record in present:
        taxon = record["taxon"]
        got = detected.get(taxon, {})
        recovered = got.get("exact_hits", 0)
        norm_obs = got.get("normalized_abundance")
        norm_exp = expected[taxon] / host_cells if host_cells else None
        ok = recovered == expected[taxon] and norm_obs == norm_exp
        lines.append(
            f"{taxon:<34}{expected[taxon]:>7}{recovered:>11}"
            f"{('%.2f' % norm_exp):>10}{(('%.2f' % norm_obs) if norm_obs is not None else '-'):>10}"
            f"  {'PASS' if ok else 'FAIL'}"
        )
    # Specificity rows.
    lines.append("-" * 86)
    absent_ok = absent["taxon"] not in detected
    lines.append(
        f"{absent['taxon']:<34}{0:>7}{'0' if absent_ok else 'DETECTED':>11}{'-':>10}{'-':>10}"
        f"  {'PASS (absent)' if absent_ok else 'FAIL'}"
    )
    noise_ok = summary["species_detected"] == 3
    lines.append(
        f"{'random noise reads (' + str(NOISE_READS) + ')':<34}{0:>7}"
        f"{('0' if noise_ok else 'SPURIOUS'):>11}{'-':>10}{'-':>10}"
        f"  {'PASS (no spurious)' if noise_ok else 'FAIL'}"
    )
    return "\n".join(lines)


# --- pytest entry point -----------------------------------------------------

def test_spike_in_recovers_abundance_and_rejects_absent():
    import pytest

    try:
        blast_exe("blastn")
        blast_exe("makeblastdb")
    except FileNotFoundError as exc:
        pytest.skip(f"BLAST+ not available: {exc}")

    check(run_spike_in())


# --- run provenance (so the artifact self-documents its conditions) ---------

def _tool_version(name: str) -> str:
    """First line of ``<tool> -version`` (e.g. 'blastn: 2.17.0+')."""
    out = subprocess.run(
        [str(blast_exe(name)), "-version"], capture_output=True, text=True, check=True
    )
    return out.stdout.strip().splitlines()[0]


def _git_provenance() -> str:
    """Commit hash + clean/dirty flag, or a note when this isn't a git checkout.

    A positive control verifies *committed* code; a dirty tree means the table
    describes code that lives nowhere, so the flag is worth surfacing.
    """
    here = Path(__file__).resolve().parent
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=here,
            capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"], cwd=here,
            capture_output=True, text=True, check=True).stdout.strip()
        return f"{head}{' (DIRTY - uncommitted changes)' if dirty else ' (clean)'}"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown (not a git checkout)"


def provenance() -> str:
    """Record the exact conditions this control ran under, above the table.

    The result is only citable if a reader can reproduce the setup: the
    committed code, the BLAST+ build (megablast tie-breaking shifts across
    releases), which binary resolved, the thread count (feeds tie-breaking into
    de-duplication), and the interpreter/OS. The seed is fixed, so the synthetic
    input is byte-identical across machines.
    """
    return "\n".join([
        "run provenance",
        "--------------",
        f"  git commit   : {_git_provenance()}",
        f"  blastn       : {_tool_version('blastn')}",
        f"  makeblastdb  : {_tool_version('makeblastdb')}",
        f"  blast binary : {blast_exe('blastn')}",
        f"  num_threads  : {default_thread_count()}",
        f"  python       : {platform.python_version()}",
        f"  platform     : {platform.platform()}",
        f"  rng seed     : {RNG_SEED}",
        "",
    ])


# --- standalone dissertation artifact ---------------------------------------

def main() -> int:
    try:
        blast_exe("blastn")
        blast_exe("makeblastdb")
    except FileNotFoundError as exc:
        print(f"BLAST+ not found: {exc}\nInstall BLAST+ or set BLAST_BIN, then retry.",
              file=sys.stderr)
        return 1

    print(provenance())
    outcome = run_spike_in()
    print(format_table(outcome))
    try:
        check(outcome)
    except AssertionError as exc:
        print(f"\nSPIKE-IN FAILED: {exc}", file=sys.stderr)
        return 1
    print("\nAll spike-in checks passed: the pipeline recovered every spiked "
          "abundance, rejected the absent species and the noise, and normalized "
          "correctly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
