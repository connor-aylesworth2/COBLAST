"""CPU-scaling benchmark for COBLAST+ local BLAST searches.

Measure how `-num_threads` (a single search) and process-level concurrency
(the batch / "many patients" path) affect wall-clock time ON YOUR hardware, so
the production defaults can be tuned to real numbers instead of guesses.

Requires NCBI BLAST+ discoverable the same way the app finds it (PATH or
BLAST_BIN). Run from the repository root. Examples:

    # How well does one search scale across cores?
    python benchmark.py --threads 1,2,4,8

    # Measure the eToL probe panel with BLAST's automatic threading mode:
    python benchmark.py --etol quick --threads 1,2,4,8

    # Compare explicit query/database splitting when diagnosing a real DB:
    python benchmark.py --etol full --db PATIENT_DB --threads 2,4,8 --mt-mode 2

    # The "100 patients" batch: how much does running patients concurrently help?
    python benchmark.py --mode batch --copies 8 --concurrency 1,2,4,8 --batch-threads 1

Nothing here touches patient data: it synthesizes a random nucleotide database
unless you pass --db / --query pointing at your own local files.
"""

from __future__ import annotations

import argparse
import random
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from time import perf_counter

from blast_runner import run_blast, run_jobs_concurrently
from config import available_cpu_count, blast_exe, default_thread_count

try:
    from etol_summary import etol_preset_fasta
except Exception:  # pragma: no cover - only needed for --etol
    etol_preset_fasta = None


def synthesize_database(workdir: Path, reads: int, read_len: int, seed: int = 0):
    """Write a random nucleotide FASTA and index it with makeblastdb."""
    random.seed(seed)
    sequences = [
        "".join(random.choice("ACGT") for _ in range(read_len)) for _ in range(reads)
    ]
    fasta = workdir / "synthetic_reads.fasta"
    with fasta.open("w", encoding="utf-8") as handle:
        for i, seq in enumerate(sequences):
            handle.write(f">read_{i}\n{seq}\n")

    prefix = workdir / "synthetic_db"
    subprocess.run(
        [str(blast_exe("makeblastdb")), "-in", str(fasta), "-dbtype", "nucl",
         "-out", str(prefix)],
        check=True, capture_output=True, text=True,
    )
    return sequences, str(prefix)


def sample_query(sequences: list[str], records: int, seed: int = 1) -> str:
    """Build a multi-FASTA query by sampling reads (so it actually has hits)."""
    random.seed(seed)
    chosen = random.sample(sequences, min(records, len(sequences)))
    return "".join(f">q_{i}\n{seq}\n" for i, seq in enumerate(chosen))


def median_time(thunk, repeat: int) -> float:
    times = []
    for _ in range(repeat):
        start = perf_counter()
        thunk()
        times.append(perf_counter() - start)
    return statistics.median(times)


def parse_int_list(text: str) -> list[int]:
    return [int(part) for part in text.split(",") if part.strip()]


def run_single(args, query, db):
    print("\n== Single-search scaling (-num_threads) ==")
    if args.mt_mode:
        print(f"   mt_mode = {args.mt_mode}")
    print(f"{'threads':>8} {'seconds':>10} {'speedup':>9} {'efficiency':>11}")
    baseline = None
    for threads in args.threads:
        seconds = median_time(
            lambda t=threads: run_blast(
                sequence=query, database=db, program=args.program,
                num_threads=t, mt_mode=args.mt_mode or None,
                exact_match_probe=args.exact,
            ),
            args.repeat,
        )
        if baseline is None:
            baseline = seconds
        speedup = baseline / seconds if seconds else float("nan")
        efficiency = speedup / threads
        print(f"{threads:>8} {seconds:>10.3f} {speedup:>8.2f}x {efficiency:>10.0%}")


def run_batch(args, query, db):
    dbs = [db] * args.copies  # each "patient" is an equal-sized DB
    print(f"\n== Batch scaling ({args.copies} patient DBs, {args.batch_threads} thread(s)/job) ==")
    print(f"{'concurrency':>11} {'oversub':>8} {'seconds':>10} {'speedup':>9}")
    baseline = None
    for concurrency in args.concurrency:
        jobs = [
            dict(sequence=query, database=patient_db, program=args.program,
                 num_threads=args.batch_threads, exact_match_probe=args.exact)
            for patient_db in dbs
        ]
        seconds = median_time(
            lambda j=jobs, c=concurrency: run_jobs_concurrently(run_blast, j, max_workers=c),
            args.repeat,
        )
        if baseline is None:
            baseline = seconds
        speedup = baseline / seconds if seconds else float("nan")
        oversub = concurrency * args.batch_threads / max(1, available_cpu_count())
        flag = "  <-- oversubscribed" if oversub > 1.0 else ""
        print(f"{concurrency:>11} {oversub:>7.2f}x {seconds:>10.3f} {speedup:>8.2f}x{flag}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["single", "batch"], default="single")
    parser.add_argument("--program", default="blastn")
    parser.add_argument("--db", help="Existing BLAST db prefix (default: synthesize one).")
    parser.add_argument("--query", help="Query FASTA path (default: sample from the db).")
    parser.add_argument("--etol", choices=["quick", "full", "control"],
                        help="Use an eToL probe panel as the query (sets exact-match mode).")
    parser.add_argument("--threads", type=parse_int_list, default=[1, 2, 4],
                        help="Comma-separated thread counts for single mode. Default: 1,2,4")
    parser.add_argument("--mt-mode", dest="mt_mode", choices=["0", "1", "2"], default=None)
    parser.add_argument("--concurrency", type=parse_int_list, default=[1, 2, 4],
                        help="Comma-separated concurrency levels for batch mode. Default: 1,2,4")
    parser.add_argument("--batch-threads", type=int, default=1,
                        help="Threads per job in batch mode. Default: 1")
    parser.add_argument("--copies", type=int, default=8,
                        help="Number of patient DBs to simulate in batch mode. Default: 8")
    parser.add_argument("--reads", type=int, default=20000, help="Synthetic DB read count.")
    parser.add_argument("--read-len", type=int, default=150, help="Synthetic read length.")
    parser.add_argument("--query-records", type=int, default=200,
                        help="How many reads to sample into the default query.")
    parser.add_argument("--repeat", type=int, default=1, help="Timed repeats (median reported).")
    args = parser.parse_args()
    args.exact = False

    try:
        blast_exe("blastn")
    except FileNotFoundError as exc:
        print(f"BLAST+ not found: {exc}\nInstall BLAST+ or set BLAST_BIN, then retry.", file=sys.stderr)
        return 1

    print(f"Machine: {available_cpu_count()} logical CPUs; adaptive default = "
          f"{default_thread_count()} thread(s)/job")

    with tempfile.TemporaryDirectory(prefix="coblast_bench_") as tmp:
        workdir = Path(tmp)

        if args.db:
            db = args.db
            sequences = []
        else:
            print(f"Synthesizing DB: {args.reads} reads x {args.read_len} bp ...")
            sequences, db = synthesize_database(workdir, args.reads, args.read_len)

        if args.etol:
            if etol_preset_fasta is None:
                print("Could not load eToL panels.", file=sys.stderr)
                return 1
            preset_key = {"quick": "etol_quick", "full": "etol_full", "control": "etol_control"}[args.etol]
            query = etol_preset_fasta(preset_key)
            args.program = "blastn"
            args.exact = True
        elif args.query:
            query = Path(args.query).read_text(encoding="utf-8")
        elif sequences:
            query = sample_query(sequences, args.query_records)
        else:
            print("Provide --query when using --db with a non-synthetic database.", file=sys.stderr)
            return 1

        # Keep all sweep values within what the machine actually has.
        ceiling = available_cpu_count()
        args.threads = [t for t in args.threads if 1 <= t <= ceiling] or [1]

        if args.mode == "single":
            run_single(args, query, db)
        else:
            run_batch(args, query, db)

    print("\nDone. Report the tables back and we'll tune the production defaults.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
