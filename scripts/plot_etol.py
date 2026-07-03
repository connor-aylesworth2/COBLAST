#!/usr/bin/env python3
"""Publication figures for an eToL-V batch.

Reproduces the eToL-V dissertation's confusion matrix (Figure 9: eToL-V calls vs
the eToL whole-genome-shotgun ground truth) and, optionally, the validated-hit
heatmap (Figure 10), as standalone image files for a thesis or paper.

This script is intentionally NOT part of the app runtime (it is the only thing
that pulls in matplotlib), so it stays out of the bundled COBLAST+ executable.
It reuses the exact same computation as the in-app panel
(:func:`etol_validation.compute_confusion`), so the figure matches the web view.

Run from the repository root:

    # from a completed eToL-V batch
    python scripts/plot_etol.py --batch-id <batch-uuid> --out fig9_confusion.png

    # from an exported matrix JSON (offline; the /etol-matrix.json download)
    python scripts/plot_etol.py --matrix-json matrix.json --out fig9.png --heatmap fig10.png

    # everything a dissertation needs in one shot: figures + data CSVs
    # (metrics vs Veso's published Fig 9, and the per-cell confusion table)
    python scripts/plot_etol.py --batch-id <uuid> --out fig9.png --heatmap fig10.png --data veso

    # numbers only, no plotting libraries required
    python scripts/plot_etol.py --matrix-json matrix.json --print-only

Requires ``matplotlib`` (``pip install matplotlib``) unless ``--print-only``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from etol_summary import build_etol_matrix, etol_preset_records  # noqa: E402
from etol_validation import _finalize, compute_confusion  # noqa: E402

# The dissertation's published Figure 9 counts (Edinburgh B270917) -- the target a
# COBLAST+ reproduction has to hit. Metrics are derived from these by the same
# _finalize the app uses, so the "published" row can't drift from the code.
VESO_FIG9 = {"tp": 9, "fp": 1, "fn": 35, "tn": 411}


def load_matrix(args: argparse.Namespace) -> tuple[dict, dict | None]:
    """Return (matrix, batch); ``batch`` is None for the --matrix-json path.

    The batch is needed only to emit the per-cell confusion CSV (--data), which
    reuses the app's ``etol_confusion_rows_as_delimited``.
    """
    if args.matrix_json:
        return json.loads(Path(args.matrix_json).read_text(encoding="utf-8")), None
    from result_store import load_batch_result

    batch = load_batch_result(args.batch_id)
    if batch.get("etol_preset_key") != "etol_v":
        raise SystemExit(f"Batch {args.batch_id} is not an eToL-V run.")
    matrix = build_etol_matrix(
        batch.get("database_results", []), etol_preset_records("etol_v")
    )
    return matrix, batch


def write_data(cm: dict, batch: dict | None, prefix: str) -> None:
    """Write the dissertation data files: metrics (reproduced vs published) + cells."""
    def r(value):
        return "" if value is None else round(value, 4)

    published = _finalize(dict(VESO_FIG9))
    metrics_path = f"{prefix}_metrics.csv"
    with open(metrics_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "TP", "FP", "FN", "TN", "N",
                    "accuracy", "precision", "recall", "F1"])
        for name, m in (("COBLAST+ (reproduced)", cm),
                        ("Veso dissertation Fig 9 (published)", published)):
            w.writerow([name, m["tp"], m["fp"], m["fn"], m["tn"], m["n"],
                        r(m["accuracy"]), r(m["precision"]), r(m["recall"]), r(m["f1"])])
    print(f"wrote {metrics_path}"
          + ("  [MATCHES published]" if (cm["tp"], cm["fp"], cm["fn"], cm["tn"])
             == (VESO_FIG9["tp"], VESO_FIG9["fp"], VESO_FIG9["fn"], VESO_FIG9["tn"])
             else "  [DIFFERS from published -- inspect cells CSV]"))

    if batch is not None:
        from result_store import etol_confusion_rows_as_delimited

        cells_path = f"{prefix}_cells.csv"
        Path(cells_path).write_text(
            etol_confusion_rows_as_delimited(batch, ","), encoding="utf-8"
        )
        print(f"wrote {cells_path}")


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


def _metrics_caption(cm: dict) -> str:
    def pct(value):
        return "n/a" if value is None else f"{value * 100:.0f}%"

    f1 = "n/a" if cm["f1"] is None else f"{cm['f1']:.2f}"
    return (
        f"Accuracy {pct(cm['accuracy'])}   Precision {pct(cm['precision'])}   "
        f"Recall {pct(cm['recall'])}   F1 {f1}   (N={cm['n']}, "
        f"{cm['scored_samples']} samples, {cm['stage']} stage)"
    )


def render_confusion(cm: dict, out_path: str, title: str) -> None:
    """Render the 2x2 confusion matrix in the dissertation's Figure 9 layout."""
    plt = _plt()
    import numpy as np

    data = np.array([[cm["tn"], cm["fp"]], [cm["fn"], cm["tp"]]])
    labels = [["TN", "FP"], ["FN", "TP"]]

    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    im = ax.imshow(data, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Predicted Negative", "Predicted Positive"])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Actual Negative", "Actual Positive"])
    ax.set_xlabel("Predicted Label")
    ax.set_ylabel("Actual Label")
    ax.set_title(title)

    threshold = data.max() / 2 if data.max() else 0
    for i in range(2):
        for j in range(2):
            ax.text(
                j, i, f"{labels[i][j]}: {data[i, j]}",
                ha="center", va="center",
                color="white" if data[i, j] > threshold else "black",
                fontsize=12,
            )
    fig.colorbar(im, ax=ax)
    fig.text(0.5, 0.02, _metrics_caption(cm), ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"wrote {out_path}")


def render_heatmap(matrix: dict, out_path: str, stage: str) -> None:
    """Render the rows-by-samples validated-hit heatmap (Figure 10 style)."""
    plt = _plt()
    import numpy as np

    counts = matrix.get("confirmed") if stage == "validated" else matrix.get("hits")
    if counts is None:
        counts = matrix.get("hits")
    data = np.array([[v or 0 for v in row] for row in counts], dtype=float)
    row_labels = [r.get("label") or r.get("key") for r in matrix["rows"]]
    col_labels = [c["sample"] for c in matrix["cols"]]

    height = max(4.0, len(row_labels) * 0.16)
    width = max(6.0, len(col_labels) * 0.28)
    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(data, cmap="coolwarm", aspect="auto")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=5)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=5, rotation=60, ha="right")
    ax.set_title(f"eToL-V probe hits ({stage})")
    fig.colorbar(im, ax=ax, label="Hit count")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"wrote {out_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--batch-id", help="A saved eToL-V batch id.")
    source.add_argument("--matrix-json", help="An exported etol-matrix.json file.")
    parser.add_argument("--out", default="etol_v_confusion.png", help="Confusion-matrix image path.")
    parser.add_argument("--heatmap", help="Also write the validated-hit heatmap here.")
    parser.add_argument("--stage", choices=("validated", "raw"), default="validated")
    parser.add_argument("--title", default="Confusion Matrix for eToL-V vs eToL Results")
    parser.add_argument("--data", metavar="PREFIX",
                        help="Also write <PREFIX>_metrics.csv (reproduced vs Veso's "
                             "published Fig 9) and, from a batch, <PREFIX>_cells.csv.")
    parser.add_argument("--print-only", action="store_true", help="Print numbers; no plotting.")
    args = parser.parse_args(argv)

    matrix, batch = load_matrix(args)
    cm = compute_confusion(matrix, stage=args.stage)

    print(
        f"TP={cm['tp']} FP={cm['fp']} FN={cm['fn']} TN={cm['tn']}  | "
        + _metrics_caption(cm)
    )
    if not cm["scored_samples"]:
        raise SystemExit(
            "No samples matched the bundled eToL WGS ground truth; nothing to plot. "
            "Run the eToL-V preset on the SRP398685 / EBB samples."
        )
    if args.data:
        write_data(cm, batch, args.data)
    if args.print_only:
        return 0

    render_confusion(cm, args.out, args.title)
    if args.heatmap:
        render_heatmap(matrix, args.heatmap, cm["stage"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
