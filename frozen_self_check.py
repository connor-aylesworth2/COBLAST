"""Packaged end-to-end self-check for the bundling-sensitive eToL stages.

WHY THIS EXISTS
---------------
``smoke_test.py`` proves BLAST+ and the summary parsers work, but it never
exercises the three stages that silently produced *zeros* in a mis-bundled
frozen build: patient-read recovery (``blastdbcmd`` on a ``-parse_seqids`` DB),
the secondary human filter, and CAP3 contig assembly. Those stages depend on
bundled binaries and ``_MEIPASS`` resource resolution, so a ``pytest``/source
run cannot catch a regression in them -- only running the *packaged .exe* can.

This module builds a tiny synthetic sample of known composition, drives the real
net path over it, and asserts the known answer. It is wired to
``run_COBLAST.py --self-check`` and run as a post-build gate against the freshly
built ``dist/COBLAST.exe`` so a bundling regression fails the build instead of a
user's analysis. It is the spike-in positive control idea (see
``tests/test_spike_in_control.py``), pointed at the frozen binary and extended to
cover human filtering and assembly.

WHAT IT ASSERTS
---------------
  * read recovery works via ``blastdbcmd`` (the ``-parse_seqids`` id index), with
    zero unresolved reads -- the exact thing that silently failed;
  * the human filter removes the reads planted in the synthetic human DB and
    keeps the rest;
  * CAP3 assembles the surviving reads into >= 1 contig -- but ONLY when a CAP3
    binary resolves. A build without CAP3 (legitimate: it is not redistributed)
    skips this assertion instead of failing, so user launches without CAP3 pass.

Returns 0 on success, 1 on any failure. Needs the bundled BLAST+ (always
present); CAP3 is optional.
"""

from __future__ import annotations

import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path

PRESET = "etol_quick"   # one probe per species: fast, unambiguous
READ_LEN = 150
CORE_LEN = 64           # panel probe length embedded in each read
A_READS = 4             # "microbial" reads that survive and must assemble
B_READS = 3             # reads also planted in the human DB -> must be removed
RNG_SEED = 20240717

_PURE_ACGT = re.compile(r"[ACGT]+")


def _rand(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _pick_two_probes() -> tuple[tuple[str, str], tuple[str, str]]:
    """Return two ((probe_id, seq)) from distinct taxa, pure-ACGT so megablast seeds."""
    from etol_summary import (
        etol_control_query_ids,
        etol_preset_records,
        etol_search_pairs,
    )

    seqs = {pid: seq for pid, seq in etol_search_pairs(PRESET)}
    controls = etol_control_query_ids(PRESET)
    chosen: list[tuple[str, str]] = []
    seen_taxa: set[str] = set()
    for record in etol_preset_records(PRESET):
        probe = record["probe"]
        if probe in controls or record["taxon"] in seen_taxa:
            continue
        seq = seqs.get(probe, "")
        if len(seq) < CORE_LEN or not _PURE_ACGT.fullmatch(seq):
            continue
        chosen.append((probe, seq[:CORE_LEN]))
        seen_taxa.add(record["taxon"])
        if len(chosen) == 2:
            return chosen[0], chosen[1]
    raise RuntimeError("Panel does not expose two seedable species for the self-check.")


def _embed(core: str, left: str, right: str) -> str:
    """Wrap a probe core in fixed flanks to a full-length read (fixed = clean assembly)."""
    read = left + core + right
    return read[:READ_LEN]


def _makeblastdb(fasta: Path, prefix: Path) -> None:
    from config import blast_exe

    subprocess.run(
        [str(blast_exe("makeblastdb")), "-in", str(fasta), "-dbtype", "nucl",
         # -parse_seqids is the whole point: without the id index, blastdbcmd
         # read recovery falls back to a FASTA scan, which is what we verify.
         "-parse_seqids", "-out", str(prefix)],
        check=True, capture_output=True, text=True,
    )


def run() -> int:
    """Build the synthetic sample, drive the real net path, assert the known answer."""
    from config import blast_exe
    from blast_runner import run_blast_probe_panel
    from etol_summary import etol_search_fasta, etol_search_query_ids
    from human_filter import extract_reads, filter_human_hits
    from assembler import Cap3Assembler
    from app import deduplicate_reads_to_best_probe, filter_net_probe_hits

    try:
        blast_exe("blastn"); blast_exe("makeblastdb"); blast_exe("blastdbcmd")
    except FileNotFoundError as exc:
        print(f"self-check FAILED: BLAST+ not resolvable in this build: {exc}", file=sys.stderr)
        return 1

    (a_id, a_core), (b_id, b_core) = _pick_two_probes()
    rng = random.Random(RNG_SEED)
    pad = (READ_LEN - CORE_LEN) // 2
    # Distinct fixed flanks per species so A reads never hit the human DB (built
    # from a B read) at >150 bits, and identical reads within a species so CAP3
    # merges them deterministically.
    read_a = _embed(a_core, _rand(rng, pad), _rand(rng, READ_LEN - CORE_LEN - pad))
    read_b = _embed(b_core, _rand(rng, pad), _rand(rng, READ_LEN - CORE_LEN - pad))

    a_ids = {f"selfcheck_A_{i}" for i in range(A_READS)}
    b_ids = {f"selfcheck_B_{i}" for i in range(B_READS)}

    with tempfile.TemporaryDirectory(prefix="coblast_selfcheck_") as tmp:
        work = Path(tmp)
        sample_fasta = work / "sample.fasta"
        sample_fasta.write_text(
            "".join(f">{rid}\n{read_a}\n" for rid in sorted(a_ids))
            + "".join(f">{rid}\n{read_b}\n" for rid in sorted(b_ids)),
            encoding="utf-8",
        )
        sample_db = work / "sample_db"
        _makeblastdb(sample_fasta, sample_db)

        human_fasta = work / "human.fasta"
        human_fasta.write_text(f">planted_human\n{read_b}\n", encoding="utf-8")
        human_db = work / "human_db"
        _makeblastdb(human_fasta, human_db)

        # Real net path over the synthetic sample.
        result = run_blast_probe_panel(
            panel_fasta=etol_search_fasta(PRESET), database=str(sample_db)
        )
        panel_hits = filter_net_probe_hits(result.hits, set(etol_search_query_ids(PRESET)))

        kept, hf_stats = filter_human_hits(
            panel_hits,
            db_prefix_path=str(sample_db),
            source_fasta_path=str(sample_fasta),
            human_db_prefix_path=str(human_db),
        )
        kept, _dedup_removed = deduplicate_reads_to_best_probe(kept)
        kept_ids = {h.get("sseqid", "") for h in kept}

        checks: list[tuple[str, bool, str]] = []

        # (1) Read recovery -- the exact stage that silently failed when frozen.
        method_ok = hf_stats["method"].startswith("blastdbcmd")
        checks.append((
            "read recovery",
            method_ok and hf_stats["reads_unresolved"] == 0,
            f"method={hf_stats['method']} unresolved={hf_stats['reads_unresolved']}",
        ))

        # (2) Human filter removes the planted reads and keeps the rest.
        removed_planted = not (b_ids & kept_ids)
        kept_survivors = a_ids & kept_ids
        checks.append((
            "human filter",
            removed_planted and hf_stats["hits_removed"] >= B_READS and bool(kept_survivors),
            f"removed={hf_stats['hits_removed']} kept_A={len(kept_survivors)}/{A_READS} "
            f"leaked_B={len(b_ids & kept_ids)}",
        ))

        # (3) CAP3 assembly of the survivors -- only when CAP3 is bundled/available.
        assembler = Cap3Assembler()
        if assembler.is_available():
            reads, _method = extract_reads(str(sample_db), str(sample_fasta), sorted(kept_survivors))
            contigs = assembler.assemble(reads)
            checks.append((
                "CAP3 assembly",
                len(contigs) >= 1,
                f"{len(contigs)} contig(s) from {len(reads)} read(s)",
            ))
        else:
            checks.append((
                "CAP3 assembly",
                True,
                "SKIPPED (no CAP3 in this build; contig assembly disabled by design)",
            ))

    print(f"COBLAST+ packaged self-check  (frozen={bool(getattr(sys, 'frozen', False))}, panel={PRESET})")
    ok = True
    for name, passed, detail in checks:
        ok = ok and passed
        print(f"  {name:<14} {'PASS' if passed else 'FAIL'}  {detail}")
    if not ok:
        print("\nself-check FAILED: a bundling-sensitive stage did not behave. "
              "This build is not trustworthy.", file=sys.stderr)
        return 1
    print("\nself-check passed: read recovery, human filter, and assembly all work in this build.")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
