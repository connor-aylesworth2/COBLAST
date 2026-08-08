#!/usr/bin/env python3
"""C4 re-tested on a measure that never touches the host-cell estimate.

C4 (burden vs age, Kohen) is computed on reads/host-cell, and burden tracks its
own denominator (F11). This re-runs the same age contrast on microbial reads per
million library reads, which the host-cell estimate never enters. Per million
rather than raw counts because Kohen libraries span 10.6M-53.8M seqs.

Reads the SAME matrix F4 scores, so a disagreement here is the normalisation and
nothing else. All 29 samples are usable: dropping samples with no host-cell
estimate is a burden-only exclusion.

    python scripts/c4_rpm.py <kohen_batch_id> <age_metadata.csv>
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_thesis import (  # noqa: E402
    KOHEN_ELDERLY_MIN, KOHEN_YOUNG_MAX, load_matrix, spearman_r, welch_t,
)

SEQS_COLUMN = "# of Seqs"  # NOT '# of Bases'/'mean read length' -- those are misaligned.


def library_sizes(path: str) -> dict[str, tuple[float, float]]:
    """Accession -> (age, library seqs), keyed on Run and Experiment both."""
    index: dict[str, tuple[float, float]] = {}
    for row in csv.DictReader(open(path, encoding="utf-8-sig")):
        if SEQS_COLUMN not in row:
            raise SystemExit(f"{path} has no {SEQS_COLUMN!r} column: {list(row)}")
        entry = (float(row["age_at_death"]), float(row[SEQS_COLUMN].replace(",", "")))
        for key in (row["Run"], row["Experiment"]):
            index[key.strip().upper()] = entry
    return index


def rpm_rows(matrix: dict, meta: dict[str, tuple[float, float]]) -> list[dict]:
    rows = []
    for c, col in enumerate(matrix["cols"]):
        key = col["sample"].strip().upper()
        if key not in meta:
            raise SystemExit(f"Kohen sample {col['sample']} not in age metadata.")
        age, seqs = meta[key]
        reads = sum((matrix["hits"][r][c] or 0) for r in range(len(matrix["rows"])))
        rows.append({"sample": col["sample"], "age": age, "reads": reads,
                     "seqs": seqs, "rpm": reads * 1e6 / seqs})
    return rows


def main(batch_id: str, metadata: str) -> None:
    matrix = load_matrix(batch_id, "etol_full")["matrix"]
    rows = rpm_rows(matrix, library_sizes(metadata))

    print(f"n = {len(rows)} (all Kohen samples; no host-cell exclusion applies)")
    print(f"rho(age, RPM)  = {spearman_r([r['age'] for r in rows], [r['rpm'] for r in rows]):+.4f}")
    print(f"rho(age, reads) = {spearman_r([r['age'] for r in rows], [float(r['reads']) for r in rows]):+.4f}")

    elderly = [r["rpm"] for r in rows if r["age"] >= KOHEN_ELDERLY_MIN]
    young = [r["rpm"] for r in rows if r["age"] <= KOHEN_YOUNG_MAX]
    test = welch_t(elderly, young)
    if test is None:
        raise SystemExit("Too few samples per age group for Welch's t-test.")
    print(f"elderly n={len(elderly)} mean={test['mean_a']:.3f} RPM | "
          f"young n={len(young)} mean={test['mean_b']:.3f} RPM")
    print(f"Welch p = {test['p']:.4f}   direction = "
          f"{'INCREASES with age' if test['mean_a'] > test['mean_b'] else 'DECREASES with age'}")
    print("\n  age      RPM     reads       seqs  sample")
    for r in sorted(rows, key=lambda r: r["age"]):
        print(f"{r['age']:5.1f} {r['rpm']:8.3f} {r['reads']:9d} {r['seqs']:10.0f}  {r['sample']}")


def self_check() -> None:
    """The join and the sum, on a matrix small enough to check by hand."""
    matrix = {
        "rows": [{}, {}],
        "cols": [{"sample": "srr001"}, {"sample": "SRX002"}],
        "hits": [[3, 10], [None, 5]],  # sample 1: 3 reads, sample 2: 15 reads
    }
    meta = {"SRR001": (70.0, 1_000_000.0), "SRX002": (50.0, 3_000_000.0)}
    rows = {r["sample"]: r for r in rpm_rows(matrix, meta)}
    assert rows["srr001"]["reads"] == 3 and rows["srr001"]["rpm"] == 3.0
    assert rows["SRX002"]["reads"] == 15 and rows["SRX002"]["rpm"] == 5.0
    assert rows["srr001"]["age"] == 70.0, "age paired to the wrong library"
    print("self-check OK")


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-check"]:
        self_check()
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        raise SystemExit(__doc__)
