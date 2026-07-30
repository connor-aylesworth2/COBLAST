#!/usr/bin/env python3
"""Recovery figure for the synthetic spike-in positive control.

Renders the two panels a verification section needs, from the SAME run the
pytest control performs -- this script imports :func:`run_spike_in` and plots
what it returns, so the figure can never disagree with the table or the test:

  Panel A  Recovery: expected (spiked) vs observed (recovered) reads, one point
           per spiked organism, against the y = x identity line. The designed
           absent species sits at (0, 0), and any OFF-TARGET detection -- an
           organism reported that was never spiked -- is drawn as its own series
           at expected = 0, with the count annotated.
  Panel B  Host-cell normalisation: computed vs expected, as an observed/expected
           ratio against a 1.0 target line, for each housekeeping control gene
           and for the derived host-cell estimate. A ratio keeps reads (~100) and
           host cells (~2) on one comparable axis; the absolutes are printed on
           each bar.

Axes are LINEAR by design. The spike levels (20/10/3 reads) span well under one
decade, so log-log is unwarranted, and the expected = 0 series (absent species,
off-targets) cannot be drawn on a log axis at all.

Like :mod:`scripts.plot_etol` this is NOT part of the app runtime -- it is the
only other thing that pulls in matplotlib, so it stays out of the bundled exe.

Run from the repository root:

    python scripts/plot_spike_in.py --out docs/fig_spike_in.pdf

The output format follows the extension (.pdf/.svg vector, .png raster).
The figure footer stamps the commit and clean/dirty state it was generated
from, so it is citable: generate it from a CLEAN tree, and commit (or ignore)
the image so a re-run does not stamp itself DIRTY.

Requires ``matplotlib`` (``pip install matplotlib``); ``--print-only`` needs it not.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tests.test_spike_in_control import (  # noqa: E402
    CONTROL_READS,
    EXPECTED_HOST_CELLS,
    NOISE_READS,
    PRESET,
    RNG_SEED,
    _git_provenance,
    check,
    run_spike_in,
)


def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit(
            "matplotlib is required for plotting (pip install matplotlib); "
            "use --print-only for numbers without it."
        ) from exc


def collect(outcome: dict) -> dict:
    """Reduce a spike-in run to exactly the series the two panels plot."""
    expected: dict[str, int] = outcome["expected"]
    detected: dict[str, dict] = outcome["detected"]
    summary = outcome["summary"]

    spiked = [
        {
            "label": record["species"],
            "expected": expected[record["taxon"]],
            "observed": detected.get(record["taxon"], {}).get("exact_hits", 0),
        }
        for record in outcome["present"]
    ]
    # Off-target = reported but never spiked. Empty on a passing control; drawn
    # (and counted) anyway so a regression shows up as points instead of silence.
    off_target = [
        {"label": detected[taxon].get("species", taxon), "observed": detected[taxon]["exact_hits"]}
        for taxon in detected
        if taxon not in expected
    ]
    # Sorted by gene name: the summary's order follows the hit list, so bars
    # would otherwise shuffle between regenerations of the same figure.
    controls = [
        {"label": gene, "observed": mean, "expected": float(CONTROL_READS), "unit": "reads"}
        for gene, mean in sorted((summary.get("control_gene_means") or {}).items())
    ]
    controls.append(
        {
            "label": "Host cells",
            "observed": float(summary["host_cells"]),
            "expected": float(EXPECTED_HOST_CELLS),
            "unit": "cells",
        }
    )
    return {
        "spiked": spiked,
        "off_target": off_target,
        "controls": controls,
        "absent": outcome["absent"]["species"],
    }


def render(series: dict, out_path: str, passed: bool, failure: str = "") -> None:
    plt = _plt()

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.5, 4.9))

    # --- Panel A: recovery vs identity -------------------------------------
    exp = [s["expected"] for s in series["spiked"]]
    obs = [s["observed"] for s in series["spiked"]]
    top = max(exp + obs + [1]) * 1.18
    ax_a.plot([0, top], [0, top], linestyle="--", linewidth=1, color="#888",
              label="y = x (perfect recovery)", zorder=1)
    ax_a.scatter(exp, obs, s=70, color="#1f77b4", zorder=3, label="Spiked organism")
    for s in series["spiked"]:
        ax_a.annotate(s["label"], (s["expected"], s["observed"]),
                      textcoords="offset points", xytext=(8, -4), fontsize=8)

    # Designed negative: a panel species deliberately left out of the sample.
    ax_a.scatter([0], [0], s=70, marker="s", facecolors="none", edgecolors="#2ca02c",
                 zorder=3, label=f"Absent species ({series['absent']})")
    off = series["off_target"]
    ax_a.scatter([0] * len(off), [o["observed"] for o in off], s=70, marker="x",
                 color="#d62728", zorder=3,
                 label=f"Off-target detections (n={len(off)})")
    for o in off:
        ax_a.annotate(o["label"], (0, o["observed"]),
                      textcoords="offset points", xytext=(8, -4), fontsize=8, color="#d62728")

    ax_a.set_xlim(-top * 0.06, top)
    ax_a.set_ylim(-top * 0.06, top)
    ax_a.set_xlabel("Expected abundance (reads spiked)")
    ax_a.set_ylabel("Observed abundance (reads recovered)")
    ax_a.set_title("A. Recovery of spiked organisms")
    ax_a.legend(fontsize=8, loc="upper left")
    ax_a.grid(alpha=0.25)

    # --- Panel B: host-cell normalisation ----------------------------------
    labels = [c["label"] for c in series["controls"]]
    ratios = [(c["observed"] / c["expected"]) if c["expected"] else 0.0
              for c in series["controls"]]
    positions = range(len(labels))
    ax_b.bar(positions, ratios, width=0.55, color="#ff7f0e")
    ax_b.axhline(1.0, linestyle="--", linewidth=1, color="#888",
                 label="Expected (ratio = 1.0)")
    for x, (c, ratio) in enumerate(zip(series["controls"], ratios)):
        ax_b.annotate(f"{c['observed']:g} / {c['expected']:g}\n{c['unit']}",
                      (x, ratio), textcoords="offset points", xytext=(0, 6),
                      ha="center", fontsize=8)
    ax_b.set_xticks(list(positions))
    ax_b.set_xticklabels(labels)
    ax_b.set_ylim(0, max(ratios + [1.0]) * 1.45)
    ax_b.set_ylabel("Observed / expected")
    ax_b.set_title("B. Host-cell normalisation")
    ax_b.legend(fontsize=8, loc="lower right")
    ax_b.grid(alpha=0.25, axis="y")

    # Two lines: the 40-char commit hash plus a verdict overflows one line at
    # this width and silently clips at the figure edges.
    verdict = "ALL CHECKS PASSED" if passed else f"CHECKS FAILED: {failure}"
    fig.text(0.5, 0.045,
             f"COBLAST+ eToL net spike-in positive control  |  panel={PRESET}  "
             f"seed={RNG_SEED}  noise reads={NOISE_READS}  |  {verdict}",
             ha="center", fontsize=7.5, color="black" if passed else "#b00020")
    fig.text(0.5, 0.015, f"commit {_git_provenance()}",
             ha="center", fontsize=7.5, color="black" if passed else "#b00020")
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"wrote {out_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="fig_spike_in.pdf",
                        help="Figure path; format follows the extension.")
    parser.add_argument("--print-only", action="store_true",
                        help="Print the plotted series; no plotting libraries needed.")
    args = parser.parse_args(argv)

    outcome = run_spike_in()
    series = collect(outcome)

    # A figure drawn from a failed control would misrepresent the pipeline, so
    # the verdict is stamped on the figure and the exit code, never hidden.
    passed, failure = True, ""
    try:
        check(outcome)
    except AssertionError as exc:
        passed, failure = False, str(exc)

    for s in series["spiked"]:
        print(f"spiked   {s['label']:<28} expected={s['expected']:>4} observed={s['observed']:>4}")
    for o in series["off_target"]:
        print(f"OFF-TARGET {o['label']:<26} expected=   0 observed={o['observed']:>4}")
    print(f"off-target detections: {len(series['off_target'])}")
    for c in series["controls"]:
        print(f"control  {c['label']:<28} expected={c['expected']:>7g} observed={c['observed']:>7g} ({c['unit']})")
    if not passed:
        print(f"SPIKE-IN FAILED: {failure}", file=sys.stderr)

    if not args.print_only:
        render(series, args.out, passed, failure)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
