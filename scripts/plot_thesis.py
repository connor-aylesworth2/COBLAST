#!/usr/bin/env python3
"""Dissertation figures and tables that are computed downstream of COBLAST+.

`plot_etol.py` renders the eToL-V confusion matrix and `plot_spike_in.py` the
verification control. This script covers everything else the write-up needs:
the analyses the tool deliberately does not perform, run on the tool's own
measurements.

  F4   Kohen age-burden           young vs elderly with every point overlaid,
                                  plus the continuous panel showing the age gap
  F5   Viral detection heatmap    individuals x regions, 5 detected viruses,
                                  ground truth beside COBLAST+ validated calls
  F6   Regional burden            paired lines, one per individual
  F8   Per-probe age association  volcano over the 120 taxa, BH threshold
  F9   AD-shortlist overlap       ranked lollipop, 9 genus + 1 domain (A3)
  F10  Region effects             forest plot, naive vs clustered (mixed model)
  T6   Validation scorecard       V2, C1-C5 against the registered criteria
  T7   Absolute burden            fold change scored, absolutes transcribed (A4)
  T8   FDR table                  per-taxon r, p, q

Governing registration: `docs/acceptance_criteria.md`, amendments A1-A5. The
decisions those encode are applied here and named at each use site, so a reader
can check the code against the registration line by line:

  A2  the cellular cutoff applies PER SPECIES BEFORE SUMMING (`CELLULAR_CUTOFF`)
  A3  the tenth AD-shortlist entry is matched at domain level, not genus
  A4  C5/T7 scores the max-to-control fold change, never the absolutes
  A5  V1 is a verification artefact and is absent here; V2 is a Miss

C1 is deliberately CUTOFF-FREE and computed from the RAW hits layer as a
pooled ratio of sums -- that is what `static/etol_pie.js` exports and what the
criterion is registered against. See :func:`domain_composition`.

Statistics are ported from `static/etol_pie.js` rather than taken from SciPy so
that a number in the write-up matches the number the tool itself would show;
`test_etol_pie_stats.js` pins the JS side against closed forms.

Run from the repository root:

    python3 scripts/plot_thesis.py \\
        --ebb-batch   1254c3b1-26ec-42f7-becd-303ccb1d4400 \\
        --viral-batch fdda627a-55cc-4ce0-9fb5-8da8ba705319 \\
        --kohen-batch c3a41385-27bc-4612-b809-47c6a64b8f64 \\
        --covariates  f10_covariates.csv \\
        --age-metadata /home/s2837739/COBLAST_2.0/age_metadata.csv \\
        --outdir figures

    # one item at a time while drafting
    python3 scripts/plot_thesis.py --only F6 ... --outdir figures

Requires matplotlib; F10 additionally requires statsmodels. Figures are written
as PDF (vector, per the figure spec) and every table as CSV alongside.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from etol_summary import build_etol_matrix, etol_preset_records  # noqa: E402
from etol_validation import compute_confusion, load_crosswalk  # noqa: E402
from result_store import load_batch_result  # noqa: E402
from tests.test_spike_in_control import _git_provenance  # noqa: E402

# --- registered constants -------------------------------------------------

# A2: a species contributes to a sample's burden only if its OWN reads per host
# cell reach this. Applied before summing, never to the total.
CELLULAR_CUTOFF = 4.0

# The four limbic regions in the fixed display order (F6's x-axis).
REGION_ORDER = ("AMYG", "BA24", "HPC", "HYPO")
# F10 contrasts every region against BA24, so BA24 is the model's reference
# level. That is a different thing from the display order above.
REGION_REFERENCE = "BA24"

# F5: the only five viruses with any positive cell in the corrected 24-virus
# ground truth. The other 19 are zero throughout and are stated, not drawn.
F5_VIRUSES = ("HAdV-C", "HSV1", "CMV", "EBV", "KSHV")

# C4's registered design. Asserted against the metadata before any Kohen number
# is computed -- a re-sorted spreadsheet that misaligns `age_at_death` would
# otherwise change every result silently and flag nothing.
KOHEN_N = 29
KOHEN_YOUNG_MAX = 59
KOHEN_ELDERLY_MIN = 68
KOHEN_EXPECTED = {"young_n": 16, "elderly_n": 13, "young_mean": 47.6, "elderly_mean": 85.2}
# The published statistic C4 scores against, annotated on F4 beside our own.
KOHEN_PUBLISHED_P = 0.0202

# A5's recorded V2 score. Recomputed here and asserted, so the scorecard cannot
# drift from the registration without failing loudly.
A5_V2_VALIDATED = {"tp": 2, "fp": 3, "fn": 43, "tn": 792}
A5_V2_RAW = {"tp": 9, "fp": 7, "fn": 36, "tn": 788}

# C5 published values, in ORGANISM units. A4: transcribed, never scored against.
C5_PUBLISHED = {"bacteria": 0.14, "fungi": 0.05, "combined": 0.19, "max_case": 1.8}
C5_PUBLISHED_FOLD = C5_PUBLISHED["max_case"] / C5_PUBLISHED["combined"]  # ~9.5x

# C3's published top-10. Nine named genera; the tenth is unnamed in the source
# and is matched at domain level under A3.
AD_SHORTLIST = (
    ("Cortinarius", "genus"),
    ("Tausonia", "genus"),
    ("Acrocalymma", "genus"),
    ("Aureobasidium", "genus"),
    ("Alternaria", "genus"),
    ("Komagataella", "genus"),
    ("Sphingomonas", "genus"),
    ("Streptococcus", "genus"),
    ("Staphylococcus", "genus"),
    ("uncharacterised Chloroplastida", "domain"),  # A3
)
# A3: the tenth entry counts as matched if ANY Chloroplastida taxon is detected
# and over-represented. Class codes C1-C4 carry that domain in the panel.
CHLOROPLASTIDA_CLASSES = ("C1", "C2", "C3", "C4")


# --- statistics, ported from static/etol_pie.js ---------------------------


def _lgamma(x: float) -> float:
    return math.lgamma(x)


def _betacf(a: float, b: float, x: float) -> float:
    """Modified Lentz continued fraction for the incomplete beta."""
    tiny = 1e-30
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-7:
            break
    return h


def ibeta(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta -- the one piece of numerics under every p."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        _lgamma(a + b) - _lgamma(a) - _lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - math.exp(
        _lgamma(a + b) - _lgamma(a) - _lgamma(b)
        + b * math.log(1.0 - x) + a * math.log(x)
    ) * _betacf(b, a, 1.0 - x) / b


def student_p(t: float, df: float) -> float:
    """Two-sided Student's t tail."""
    if df <= 0:
        return float("nan")
    return ibeta(df / 2.0, 0.5, df / (df + t * t))


def welch_t(a: list[float], b: list[float]) -> dict | None:
    """Welch's unequal-variance t-test -- what R's t.test() does by default.

    Chosen because it is the most likely thing behind the source's published
    't-test', which is what C4 scores against.
    """
    if len(a) < 2 or len(b) < 2:
        return None
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    sa, sb = va / len(a), vb / len(b)
    denom = sa + sb
    if denom <= 0:
        return None
    t = (ma - mb) / math.sqrt(denom)
    df = denom * denom / (
        sa * sa / (len(a) - 1) + sb * sb / (len(b) - 1)
    )
    return {"t": t, "df": df, "p": student_p(t, df), "mean_a": ma, "mean_b": mb}


def bh(ps: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, returned in input order.

    Descending step-up with a running minimum, so the output is monotone and
    capped at 1 -- the same construction as `bh()` in etol_pie.js.
    """
    n = len(ps)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: ps[i], reverse=True)
    out = [1.0] * n
    running = 1.0
    for k, i in enumerate(order):
        rank = n - k
        running = min(running, ps[i] * n / rank)
        out[i] = min(1.0, running)
    return out


def cohens_d(a: list[float], b: list[float]) -> float:
    """Standardised difference, pooled SD. F8's primary effect size.

    Primary rather than r because the Kohen ages are bimodal (16 at 29-59, 13
    at 68-95, nothing between): a correlation over that gap is a two-group
    difference wearing a correlation's clothes.
    """
    if len(a) < 2 or len(b) < 2:
        return 0.0
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    pooled = ((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2)
    return 0.0 if pooled <= 0 else (ma - mb) / math.sqrt(pooled)


def spearman_r(xs: list[float], ys: list[float]) -> float:
    """Spearman rho with midranks. F8's SECONDARY effect size only."""
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            mid = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k]] = mid
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den = math.sqrt(sum((v - mx) ** 2 for v in rx) * sum((v - my) ** 2 for v in ry))
    return 0.0 if den == 0 else num / den


# --- loading --------------------------------------------------------------


def load_matrix(batch_id: str, preset: str, level: str = "species") -> dict:
    """Batch id -> the rows x samples matrix the app itself renders from."""
    batch = load_batch_result(batch_id)
    key = batch.get("etol_preset_key")
    if key != preset:
        raise SystemExit(
            f"Batch {batch_id} is preset {key!r}, expected {preset!r}. "
            "Check the batch ids -- scoring the wrong preset silently produces "
            "a plausible, wrong figure."
        )
    matrix = build_etol_matrix(
        batch.get("database_results", []), etol_preset_records(preset), level=level
    )
    return {"batch": batch, "matrix": matrix}


def load_covariates(path: str) -> dict[str, dict]:
    """EBB covariates keyed by BOTH accessions.

    `ETOL_ACCESSION_PATTERN` matches SRX or SRR and takes whichever appears
    first in a database's name, so the sample label is not predictable from
    here. Keying on both makes the join work either way -- the same trick
    `etol_validation.load_crosswalk` uses.
    """
    index: dict[str, dict] = {}
    for row in csv.DictReader(open(path, encoding="utf-8-sig")):
        for key in (row["srx"], row["srr"], row["sample_name"]):
            index[key.strip().upper()] = row
    return index


def load_age_metadata(path: str) -> dict[str, dict]:
    """Kohen age/sex keyed by both accessions, with C4's design asserted.

    Only `Run`, `Experiment`, `age_at_death` and `gender` are read. The file's
    `# of Bases` / `mean read length` columns are NOT used: they were misaligned
    by a spreadsheet sort and disagree with `Bases` on every row. The human
    filter derives its threshold from the reads themselves (A1), so nothing in
    the pipeline depends on them either.
    """
    index: dict[str, dict] = {}
    ages: list[float] = []
    for row in csv.DictReader(open(path, encoding="utf-8-sig")):
        age = float(row["age_at_death"])
        entry = {
            "age": age,
            "sex": (row.get("gender") or "").strip(),
            "run": row["Run"].strip(),
            "srx": row["Experiment"].strip(),
        }
        for key in (entry["run"], entry["srx"]):
            index[key.upper()] = entry
        ages.append(age)

    young = [a for a in ages if a <= KOHEN_YOUNG_MAX]
    elderly = [a for a in ages if a >= KOHEN_ELDERLY_MIN]
    problems = []
    if len(ages) != KOHEN_N:
        problems.append(f"{len(ages)} samples, expected {KOHEN_N}")
    if len(young) != KOHEN_EXPECTED["young_n"]:
        problems.append(f"young n={len(young)}, expected {KOHEN_EXPECTED['young_n']}")
    if len(elderly) != KOHEN_EXPECTED["elderly_n"]:
        problems.append(f"elderly n={len(elderly)}, expected {KOHEN_EXPECTED['elderly_n']}")
    if len(young) + len(elderly) != len(ages):
        problems.append("a sample falls inside the registered 59-68 age gap")
    for label, got, want in (
        ("young mean", sum(young) / max(len(young), 1), KOHEN_EXPECTED["young_mean"]),
        ("elderly mean", sum(elderly) / max(len(elderly), 1), KOHEN_EXPECTED["elderly_mean"]),
    ):
        if abs(got - want) > 0.1:
            problems.append(f"{label} {got:.1f}, registered {want}")
    if problems:
        raise SystemExit(
            "Kohen metadata does not match C4's registered design: "
            + "; ".join(problems)
            + ".\nRefusing to compute -- a misaligned age column would change "
            "every result and flag nothing."
        )
    return index


# --- shared quantities ----------------------------------------------------


def per_cell(matrix: dict, row: int, col: int) -> float:
    """Reads per host cell for one taxon in one sample."""
    host_cells = matrix["cols"][col]["host_cells"] or 0.0
    if host_cells <= 0:
        return 0.0
    return (matrix["hits"][row][col] or 0) / host_cells


def unmeasurable_samples(matrix: dict) -> list[str]:
    """Samples with no host-cell estimate, which cannot be normalised at all.

    `per_cell` returns 0 when host_cells <= 0, so such a sample silently scores a
    burden of 0 and is then averaged in as though it were measured and empty.
    That is a structural zero, not an observation, and it biases every group
    mean it lands in. Callers exclude these and say so.
    """
    return [col["sample"] for col in matrix["cols"] if (col["host_cells"] or 0.0) <= 0]


def sample_burden(matrix: dict, col: int, cutoff: float = CELLULAR_CUTOFF) -> float:
    """A2: sum reads per host cell over species that individually reach cutoff.

    NOTE the denominator's leverage: burden is Sum(hits / host_cells), so when
    host_cells varies more across samples than microbial content does, this
    measure tracks 1/host_cells. Check `burden ~ host_cells` before reading any
    between-group difference as biology.
    """
    total = 0.0
    for row in range(len(matrix["rows"])):
        value = per_cell(matrix, row, col)
        if value >= cutoff:
            total += value
    return total


def domain_composition(matrix: dict, cols: list[int] | None = None) -> "OrderedDict[str, dict]":
    """C1's quantity, reproducing `exportCsv` in static/etol_pie.js exactly.

    Pooled ratio of sums -- total raw hits per domain divided by total host
    cells across the scope -- NOT a mean of per-sample ratios, and with NO
    cutoff applied. That is what the criterion is registered against; A2's
    cutoff names C2, C4 and T7, not C1.
    """
    cols = list(range(len(matrix["cols"]))) if cols is None else cols
    reads: dict[str, int] = defaultdict(int)
    detected_in: dict[str, int] = defaultdict(int)
    for r, row in enumerate(matrix["rows"]):
        for c in cols:
            hits = matrix["hits"][r][c] or 0
            if hits > 0:
                reads[row["domain"]] += hits
                detected_in[row["domain"]] += 1
    host_cells = sum(matrix["cols"][c]["host_cells"] or 0.0 for c in cols)
    grand = sum(reads.values())
    out: "OrderedDict[str, dict]" = OrderedDict()
    for domain in sorted(reads, key=lambda d: reads[d], reverse=True):
        out[domain] = {
            "reads": reads[domain],
            "percent": 100.0 * reads[domain] / grand if grand else 0.0,
            "per_host_cell": reads[domain] / host_cells if host_cells > 0 else 0.0,
            "samples_detected": detected_in[domain],
        }
    return out


def ebb_sample_index(matrix: dict, covariates: dict[str, dict]) -> list[dict]:
    """Attach individual/region/diagnosis to each EBB column, or fail loudly."""
    out = []
    missing = []
    for c, col in enumerate(matrix["cols"]):
        cov = covariates.get(col["sample"].strip().upper())
        if cov is None:
            missing.append(col["sample"])
            continue
        out.append({
            "col": c,
            "sample": col["sample"],
            "individual": cov["individual"],
            "region": cov["region"],
            "diagnosis": cov["diagnosis2"],
            "diagnosis_full": cov["diagnosis_full"],
            "burden": sample_burden(matrix, c),
        })
    if missing:
        raise SystemExit(
            f"{len(missing)} EBB sample(s) not in the covariates file: "
            f"{', '.join(missing[:5])}. The join key is the accession in the "
            "database name; check srx/srr coverage."
        )
    return out


def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit("matplotlib is required (pip install matplotlib)") from exc


def _finish(fig, note: str, path: Path) -> None:
    """Caption the figure with its provenance, reserve room for it, save, close.

    Every figure routes through here so the caption cannot silently clip at the
    figure edges. It did: at one line these notes plus a 40-char commit hash run
    off both sides, exactly the failure `plot_spike_in.py` hit and solved by
    splitting onto two lines. Wrapping to the figure's own width fixes it for
    any caption length instead of one.
    """
    import textwrap

    plt = _plt()
    lines = textwrap.wrap(f"{note}  |  commit {_git_provenance()}",
                          width=max(60, int(fig.get_size_inches()[0] * 18)))
    reserved = 0.021 * len(lines) + 0.022
    fig.tight_layout(rect=(0, reserved, 1, 1))
    for i, line in enumerate(reversed(lines)):
        fig.text(0.5, 0.008 + i * 0.019, line, ha="center", fontsize=6.5, color="#555")
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def _write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"wrote {path}")


# --- F5: viral detection heatmap -----------------------------------------


def load_viral_truth() -> dict[tuple[str, str], int]:
    """(individual, region) -> per-virus counts from the corrected ground truth."""
    path = REPO_ROOT / "data" / "etol_v_ground_truth.csv"
    rows = list(csv.reader(path.read_text(encoding="utf-8-sig").splitlines()))
    header = [h.strip() for h in rows[0][2:]]
    truth: dict[tuple[str, str], dict[str, int]] = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        number, _, year = row[0].strip().partition("^")
        individual = f"SD{int(number):03d}-{year}"
        truth[(individual, row[1].strip())] = {
            name: int(float(value or 0)) for name, value in zip(header, row[2:])
        }
    return truth


def figure_f5(viral: dict, covariates: dict[str, dict], outdir: Path) -> None:
    """Ground truth beside COBLAST+ validated calls, individuals x regions.

    Two panels because V2 missed: panel B is the honest second failure figure,
    showing where the calls are absent against a populated truth.
    """
    plt = _plt()
    import numpy as np

    matrix = viral["matrix"]
    truth = load_viral_truth()
    individuals = sorted({cov["individual"] for cov in covariates.values()})

    # COBLAST+ validated calls, summed per virus token over that virus's taxa.
    from etol_validation import universe_taxa, load_ground_truth

    _truth, universe = load_ground_truth()
    taxa_for = universe_taxa(universe)
    row_index = {row["key"]: i for i, row in enumerate(matrix["rows"])}
    confirmed = matrix.get("confirmed") or matrix["hits"]
    calls: dict[tuple[str, str], dict[str, int]] = {}
    for c, col in enumerate(matrix["cols"]):
        cov = covariates.get(col["sample"].strip().upper())
        if cov is None:
            continue
        cell = {}
        for virus in F5_VIRUSES:
            cell[virus] = sum(
                (confirmed[row_index[t]][c] or 0)
                for t in taxa_for.get(virus, frozenset())
                if t in row_index
            )
        calls[(cov["individual"], cov["region"])] = cell

    n_rows = len(F5_VIRUSES) * len(individuals)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, max(6.0, n_rows * 0.22)))
    for ax, (source, title) in zip(
        axes, ((truth, "A. Corrected ground truth"), (calls, "B. COBLAST+ validated calls"))
    ):
        data = np.full((n_rows, len(REGION_ORDER)), np.nan)
        labels = []
        for vi, virus in enumerate(F5_VIRUSES):
            for ii, individual in enumerate(individuals):
                r = vi * len(individuals) + ii
                labels.append(f"{virus}  {individual}")
                for ci, region in enumerate(REGION_ORDER):
                    cell = source.get((individual, region))
                    if cell is None:
                        continue  # not sampled -> stays NaN, drawn as absent
                    data[r, ci] = cell.get(virus, 0)

        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("#d9d9d9")  # SD042-18 BA24: absent, never zero
        shown = np.sqrt(np.ma.masked_invalid(data))
        im = ax.imshow(shown, cmap=cmap, aspect="auto")
        ax.set_xticks(range(len(REGION_ORDER)))
        ax.set_xticklabels(REGION_ORDER)
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(labels, fontsize=5.5)
        ax.set_title(title, fontsize=10)
        for r in range(n_rows):
            for ci in range(len(REGION_ORDER)):
                value = data[r, ci]
                text = "--" if np.isnan(value) else f"{int(value)}"
                ax.text(ci, r, text, ha="center", va="center", fontsize=5,
                        color="white" if (not np.isnan(value) and value > 0) else "#333")
        for vi in range(1, len(F5_VIRUSES)):
            ax.axhline(vi * len(individuals) - 0.5, color="black", linewidth=0.8)
        fig.colorbar(im, ax=ax, label="sqrt(reads)")

    _finish(fig, "F5  19 of the 24 panel viruses are zero in every cell and are not drawn; "
                "grey = not sampled, distinct from zero", outdir / "F5_viral_heatmap.pdf")


# --- F6: regional burden, paired lines ------------------------------------


def figure_f6(samples: list[dict], outdir: Path) -> None:
    """One line per individual across their own regions. The form is the argument."""
    plt = _plt()

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    colours = {"AD": "#d55e00", "CONTROL": "#0072b2"}
    by_individual: dict[str, dict[str, float]] = defaultdict(dict)
    diagnosis: dict[str, str] = {}
    for s in samples:
        by_individual[s["individual"]][s["region"]] = s["burden"]
        diagnosis[s["individual"]] = s["diagnosis"]

    for individual, burdens in sorted(by_individual.items()):
        # SD042-18 has no BA24: plot only the regions it contributed, so the
        # line carries a visible gap rather than an interpolation.
        xs = [i for i, region in enumerate(REGION_ORDER) if region in burdens]
        ys = [burdens[REGION_ORDER[i]] for i in xs]
        ax.plot(xs, ys, marker="o", markersize=4, linewidth=1.4,
                color=colours[diagnosis[individual]], alpha=0.85)
        ax.annotate(individual, (xs[-1], ys[-1]), textcoords="offset points",
                    xytext=(5, 0), fontsize=6, color=colours[diagnosis[individual]])

    for group, colour in colours.items():
        members = [i for i, d in diagnosis.items() if d == group]
        means = []
        for i, region in enumerate(REGION_ORDER):
            values = [by_individual[m][region] for m in members if region in by_individual[m]]
            means.append(sum(values) / len(values) if values else float("nan"))
        ax.plot(range(len(REGION_ORDER)), means, color=colour, linewidth=3,
                alpha=0.35, zorder=0, label=f"{group} mean (n={len(members)})")

    ax.set_xticks(range(len(REGION_ORDER)))
    ax.set_xticklabels(REGION_ORDER)
    ax.set_xlabel("Brain region")
    ax.set_ylabel(f"Microbial burden (reads per host cell, species >= {CELLULAR_CUTOFF:g})")
    ax.set_title("Regional microbial burden, paired within individual")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    _finish(fig, f"F6  n = {len(samples)} samples from {len(by_individual)} individuals; "
                "SD042-18 contributed 3 regions (no BA24) and its line is broken, not interpolated", outdir / "F6_regional_burden_paired.pdf")


# --- F4 / C4: Kohen age-burden --------------------------------------------


def kohen_burdens(matrix: dict, ages: dict[str, dict]) -> list[dict]:
    """Per-sample age, sex and total burden. One computation, two consumers.

    F4 draws these and the C4 scorecard row scores them, so the figure and the
    verdict can never disagree.
    """
    dropped = unmeasurable_samples(matrix)
    if dropped:
        print(f"WARNING: {len(dropped)} Kohen sample(s) have no host-cell estimate "
              f"and are excluded from burden: {', '.join(dropped)}. "
              "Report n analysed separately from n sequenced.")
    rows = []
    for c, col in enumerate(matrix["cols"]):
        meta = ages.get(col["sample"].strip().upper())
        if meta is None:
            raise SystemExit(
                f"Kohen sample {col['sample']} is not in the age metadata; "
                "check the Run/Experiment accessions."
            )
        if col["sample"] in dropped:
            continue
        rows.append({
            "col": c, "sample": col["sample"], "age": meta["age"],
            "sex": meta["sex"], "burden": sample_burden(matrix, c),
        })
    return rows


def kohen_age_test(rows: list[dict]) -> dict:
    """Welch on total burden, elderly vs young, with C4's verdict attached."""
    young = [r["burden"] for r in rows if r["age"] <= KOHEN_YOUNG_MAX]
    elderly = [r["burden"] for r in rows if r["age"] >= KOHEN_ELDERLY_MIN]
    test = welch_t(elderly, young)
    if test is None:
        raise SystemExit("Too few Kohen samples per group to run Welch's t-test.")
    increased = test["mean_a"] > test["mean_b"]
    return {
        **test, "young": young, "elderly": elderly, "increased": increased,
        "verdict": "PASS" if increased and test["p"] < 0.05
                   else "PARTIAL" if increased else "MISS",
    }


def figure_f4(rows: list[dict], outdir: Path) -> dict:
    """Two panels: the source's own grouping, and the continuous truth behind it."""
    plt = _plt()

    result = kohen_age_test(rows)
    young = [r for r in rows if r["age"] <= KOHEN_YOUNG_MAX]
    elderly = [r for r in rows if r["age"] >= KOHEN_ELDERLY_MIN]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.0, 5.0))

    # --- Panel A: the source's binarised presentation -----------------------
    groups = [[r["burden"] for r in young], [r["burden"] for r in elderly]]
    ax_a.boxplot(groups, positions=[0, 1], widths=0.5, showfliers=False,
                 medianprops={"color": "#333"})
    # At n=29 every point is shown: a box alone discards what a reader wants.
    rng = __import__("random").Random(20220317)  # fixed seed -> stable jitter
    for x, group in enumerate(groups):
        for value in group:
            ax_a.scatter(x + rng.uniform(-0.13, 0.13), value, s=26, zorder=3,
                         color="#0072b2" if x == 0 else "#d55e00", alpha=0.8)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels([
        f"Young\nn={len(young)}, mean {sum(r['age'] for r in young)/len(young):.1f} y",
        f"Elderly\nn={len(elderly)}, mean {sum(r['age'] for r in elderly)/len(elderly):.1f} y",
    ])
    ax_a.set_ylabel(f"Microbial burden (reads per host cell, species >= {CELLULAR_CUTOFF:g})")
    ax_a.set_title("A. Burden by age group, as the source grouped it")
    ax_a.annotate(
        f"Welch t = {result['t']:.2f}, p = {result['p']:.4f}\n"
        f"published P = {KOHEN_PUBLISHED_P}\nC4 verdict: {result['verdict']}",
        xy=(0.03, 0.97), xycoords="axes fraction", va="top", fontsize=8,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#bbb"},
    )
    ax_a.grid(alpha=0.25, axis="y")

    # --- Panel B: the continuous distribution, gap shown --------------------
    # No trend line. The 59-68 gap means a slope would be two clusters levering
    # a line, which is precisely the artifact this panel exists to expose.
    for group, colour, label in ((young, "#0072b2", "Young"), (elderly, "#d55e00", "Elderly")):
        ax_b.scatter([r["age"] for r in group], [r["burden"] for r in group],
                     s=34, color=colour, label=label, zorder=3)
        mean_burden = sum(r["burden"] for r in group) / len(group)
        ax_b.plot([min(r["age"] for r in group), max(r["age"] for r in group)],
                  [mean_burden, mean_burden], color=colour, linewidth=2, alpha=0.55)
    ax_b.axvspan(KOHEN_YOUNG_MAX, KOHEN_ELDERLY_MIN, color="#bbb", alpha=0.3, zorder=0)
    ax_b.annotate(f"no samples\n{KOHEN_YOUNG_MAX:.0f}-{KOHEN_ELDERLY_MIN:.0f} y",
                  xy=((KOHEN_YOUNG_MAX + KOHEN_ELDERLY_MIN) / 2, ax_b.get_ylim()[1]),
                  ha="center", va="top", fontsize=7.5, color="#555")
    ax_b.set_xlabel("Age at death (years)")
    ax_b.set_ylabel("Microbial burden (reads per host cell)")
    ax_b.set_title("B. Burden against age, continuous")
    ax_b.legend(fontsize=8)
    ax_b.grid(alpha=0.25)

    rho = spearman_r([r["age"] for r in rows], [r["burden"] for r in rows])
    _finish(fig, f"F4  n = {len(rows)} independent individuals, one sample each; horizontal bars "
                f"are group means, NOT a fitted slope -- ages are bimodal with no samples between "
                f"{KOHEN_YOUNG_MAX:.0f} and {KOHEN_ELDERLY_MIN:.0f} y, so Spearman rho = {rho:.2f} "
                "is reported as secondary only", outdir / "F4_kohen_age_burden.pdf")

    _write_csv(outdir / "F4_kohen_burdens.csv",
               ["sample", "age", "sex", "group", "burden"],
               [[r["sample"], r["age"], r["sex"],
                 "young" if r["age"] <= KOHEN_YOUNG_MAX else "elderly",
                 round(r["burden"], 4)] for r in sorted(rows, key=lambda r: r["age"])])
    return result


# --- F8 / T8: per-taxon age association -----------------------------------


def kohen_taxon_stats(matrix: dict, ages: dict[str, dict]) -> list[dict]:
    """Per-taxon young-vs-elderly test over the registered 120-taxon unit.

    Taxa with no reads in any sample carry no information and cannot be tested;
    they are excluded from the BH family and both counts are reported, so the
    correction is over what was actually tested rather than over 120 by fiat.
    """
    young_cols, elderly_cols = [], []
    for c, col in enumerate(matrix["cols"]):
        meta = ages.get(col["sample"].strip().upper())
        if meta is None:
            raise SystemExit(
                f"Kohen sample {col['sample']} is not in the age metadata; "
                "check the Run/Experiment accessions."
            )
        (young_cols if meta["age"] <= KOHEN_YOUNG_MAX else elderly_cols).append(c)

    age_by_col = {
        c: ages[col["sample"].strip().upper()]["age"]
        for c, col in enumerate(matrix["cols"])
    }
    results = []
    for r, row in enumerate(matrix["rows"]):
        elderly = [per_cell(matrix, r, c) for c in elderly_cols]
        young = [per_cell(matrix, r, c) for c in young_cols]
        if not any(elderly) and not any(young):
            continue  # all-zero taxon: nothing to test
        test = welch_t(elderly, young)
        if test is None:
            continue
        every = [per_cell(matrix, r, c) for c in sorted(age_by_col)]
        results.append({
            "taxon": row["key"],
            "species": row["species"],
            "domain": row["domain"],
            "group": row["group"],
            "d": cohens_d(elderly, young),
            "r": spearman_r([age_by_col[c] for c in sorted(age_by_col)], every),
            "p": test["p"],
            "mean_elderly": test["mean_a"],
            "mean_young": test["mean_b"],
        })
    qs = bh([res["p"] for res in results])
    for res, q in zip(results, qs):
        res["q"] = q
    return sorted(results, key=lambda res: res["p"])


def figure_f8_and_t8(results: list[dict], panel_taxa: int, outdir: Path) -> None:
    plt = _plt()

    survivors = [res for res in results if res["q"] < 0.05]
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    xs = [res["d"] for res in results]
    ys = [-math.log10(max(res["p"], 1e-300)) for res in results]
    ax.scatter(xs, ys, s=18, color="#999", label=f"tested (n={len(results)})")
    if survivors:
        ax.scatter([r["d"] for r in survivors],
                   [-math.log10(max(r["p"], 1e-300)) for r in survivors],
                   s=34, color="#d55e00", label=f"q < 0.05 (n={len(survivors)})")
        for res in survivors:
            ax.annotate(res["species"], (res["d"], -math.log10(max(res["p"], 1e-300))),
                        textcoords="offset points", xytext=(6, 2), fontsize=7)

    # BH threshold: the largest raw p that still clears q<0.05. Drawn only when
    # something clears it -- a line at an imaginary threshold would be a claim.
    cleared = [res["p"] for res in results if res["q"] < 0.05]
    if cleared:
        ax.axhline(-math.log10(max(cleared)), linestyle="--", color="#d55e00",
                   linewidth=1, label="BH threshold, q = 0.05")
    else:
        ax.text(0.5, 0.94, f"0 of {len(results)} taxa survive at q < 0.05",
                transform=ax.transAxes, ha="center", fontsize=10, color="#d55e00")

    ax.axvline(0, color="#555", linewidth=0.8)
    ax.set_xlabel("Standardised difference, elderly - young (Cohen's d)")
    ax.set_ylabel("-log10(p), Welch's t-test")
    ax.set_title("Per-taxon association with age, Kohen cohort")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)
    _finish(fig, f"F8  {len(results)} of {panel_taxa} panel taxa carried reads and were tested; "
                "BH over the tested set. n = 29 independent individuals, one sample each", outdir / "F8_age_volcano.pdf")

    reported = [res for res in results if res["q"] < 0.05] or results[:10]
    _write_csv(
        outdir / "T8_fdr_taxa.csv",
        ["taxon", "species", "domain", "cohens_d", "spearman_r", "p", "q",
         "mean_elderly", "mean_young", "survives"],
        [[res["taxon"], res["species"], res["domain"], round(res["d"], 4),
          round(res["r"], 4), f"{res['p']:.3g}", f"{res['q']:.3g}",
          round(res["mean_elderly"], 4), round(res["mean_young"], 4),
          "yes" if res["q"] < 0.05 else "no"] for res in reported],
    )


# --- F9: AD-shortlist overlap ---------------------------------------------


def homolog_genus(homolog: str) -> str:
    """Genus from a SILVA hit title.

    `closest_homolog` carries the FULL reference title: an accession, a space,
    then the whole ``;``-delimited lineage --

        LRBV01001970.46707.48505 Eukaryota;...;Quercus;Quercus lobata

    The genus is the SECOND-TO-LAST rank. The last rank is a species string that
    is routinely not a binomial (``uncultured fungus``, ``metagenome``,
    ``Sphingomonadaceae bacterium KVD-unk-19``), so it cannot be split on
    whitespace to recover a genus either.

    Splitting the whole title on its first space returns the ACCESSION, which
    matches no genus at all. That was the original implementation and it scored
    every named genus in the AD shortlist as absent -- a wrong result that looked
    exactly like a real one, which is why this has a test.
    """
    ranks = [part.strip() for part in homolog.split(";") if part.strip()]
    return ranks[-2] if len(ranks) >= 2 else ""


def figure_f9(ebb: dict, samples: list[dict], outdir: Path) -> list[dict]:
    """Ranked lollipop over the published top 10. Never a Venn."""
    plt = _plt()

    matrix = ebb["matrix"]
    batch = ebb["batch"]
    # Genus comes from contig identification, not probe labels: no shortlist
    # genus exists in the panel's taxon names (they are coded reference
    # organisms), which is why A3 exists at all.
    # A taxon's contig identifies independently in EVERY sample, and the calls
    # differ between them. Keeping only the first one seen (setdefault) scored
    # Sphingomonas as absent from EBB despite 27 hits across three SILVA
    # entries, so collect every genus each taxon ever resolved to.
    genera_by_taxon: dict[str, set[str]] = defaultdict(set)
    for result in batch.get("database_results", []):
        for taxon, ident in (result.get("contig_identification") or {}).items():
            genus = homolog_genus((ident.get("closest_homolog") or "").strip())
            if genus:
                genera_by_taxon[taxon].add(genus)

    ad_cols = [s["col"] for s in samples if s["diagnosis"] == "AD"]
    control_cols = [s["col"] for s in samples if s["diagnosis"] == "CONTROL"]

    def over_representation(rows: list[int]) -> tuple[float, float]:
        """(proportion of AD samples above the control mean, effect size)."""
        if not rows:
            return 0.0, 0.0
        control = [sum(per_cell(matrix, r, c) for r in rows) for c in control_cols]
        ad = [sum(per_cell(matrix, r, c) for r in rows) for c in ad_cols]
        control_mean = sum(control) / len(control) if control else 0.0
        above = sum(1 for v in ad if v > control_mean) / len(ad) if ad else 0.0
        return above, cohens_d(ad, control)

    entries = []
    for name, level in AD_SHORTLIST:
        if level == "genus":
            rows = [
                r for r, row in enumerate(matrix["rows"])
                if any(g.lower() == name.lower()
                       for g in genera_by_taxon.get(row["key"], ()))
            ]
            carrier = ", ".join(sorted({matrix["rows"][r]["species"] for r in rows})) or "-"
        else:  # A3: domain-level match for the unnamed chloroplastida
            rows = [
                r for r, row in enumerate(matrix["rows"])
                if row["group"] in CHLOROPLASTIDA_CLASSES
                and any(matrix["hits"][r][c] for c in ad_cols + control_cols)
            ]
            carrier = ", ".join(sorted({matrix["rows"][r]["species"] for r in rows})) or "-"
        proportion, effect = over_representation(rows)
        entries.append({
            "name": name, "level": level, "detected": bool(rows),
            "proportion": proportion, "effect": effect, "carrier": carrier,
            # Genera, not full SILVA titles: the titles are ~200 chars each and
            # made the CSV unreadable. What a reader needs is WHAT the carriers
            # identified as -- which for the Chloroplastida entry is the whole
            # argument against counting it.
            "homologs": ", ".join(sorted(
                {g for r in rows for g in genera_by_taxon.get(matrix["rows"][r]["key"], ())}
            )) or "-",
            # C3's rule: over-represented if the AD proportion clears half the
            # control proportion. Control proportion is 0.5 by construction of
            # the mean, so the bar is 0.25.
            "over": bool(rows) and proportion > 0.25,
        })

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ys = range(len(entries))
    for y, entry in zip(ys, entries):
        colour = "#d55e00" if entry["over"] else ("#888" if entry["detected"] else "#ccc")
        ax.plot([0, entry["effect"]], [y, y], color=colour, linewidth=1.5, zorder=1)
        ax.scatter([entry["effect"]], [y], s=70, color=colour, zorder=2,
                   marker="o" if entry["level"] == "genus" else "s")
    ax.axvline(0, color="#555", linewidth=0.9)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([
        f"{e['name']}{' (domain level)' if e['level']=='domain' else ''}" for e in entries
    ], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("AD-vs-control standardised difference (reads per host cell)")
    ax.set_title("Published AD shortlist, scored in EBB")
    matched = sum(1 for e in entries if e["over"])
    ax.grid(alpha=0.25, axis="x")
    _finish(fig, f"F9  {matched}/10 over-represented; square marker = matched at DOMAIN level "
                 "because the source named it via the 23S/28S LSU pass COBLAST+ does not "
                 "implement, so genus matching caps at 9/10 by construction (A3)",
            outdir / "F9_ad_shortlist_overlap.pdf")

    _write_csv(
        outdir / "F9_ad_shortlist_overlap.csv",
        ["entry", "match_level", "detected", "ad_proportion_above_control_mean",
         "cohens_d", "carrier_taxa", "closest_homologs", "over_represented"],
        [[e["name"], e["level"], "yes" if e["detected"] else "no",
          round(e["proportion"], 4), round(e["effect"], 4), e["carrier"],
          e["homologs"], "yes" if e["over"] else "no"] for e in entries],
    )
    return entries


# --- F10: region effects, naive vs clustered ------------------------------


def figure_f10(samples: list[dict], outdir: Path) -> dict:
    """Two series: naive sample-level OLS vs random-intercept mixed model.

    The comparison IS the figure. Showing only the clustered model would hide
    the size of the correction, which is the methodological point.
    """
    try:
        import pandas as pd
        import statsmodels.formula.api as smf
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit("F10 needs statsmodels and pandas (pip install statsmodels)") from exc

    plt = _plt()

    frame = pd.DataFrame(samples)
    others = [r for r in REGION_ORDER if r != REGION_REFERENCE]
    formula = f"burden ~ C(region, Treatment('{REGION_REFERENCE}'))"

    naive = smf.ols(formula, data=frame).fit()
    # Random intercept per individual: region is a WITHIN-individual contrast,
    # which is where this design's power actually lives.
    clustered = smf.mixedlm(formula, data=frame, groups=frame["individual"]).fit()

    def extract(fit):
        out = {}
        for region in others:
            term = f"C(region, Treatment('{REGION_REFERENCE}'))[T.{region}]"
            if term not in fit.params:
                continue
            estimate = float(fit.params[term])
            ci = fit.conf_int()
            out[region] = (estimate, float(ci.loc[term][0]), float(ci.loc[term][1]))
        return out

    series = {"Naive (35 samples as independent)": (extract(naive), "#888", -0.16),
              "Clustered (random intercept per individual)": (extract(clustered), "#d55e00", 0.16)}

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for label, (values, colour, offset) in series.items():
        first = True
        for i, region in enumerate(others):
            if region not in values:
                continue
            estimate, low, high = values[region]
            y = i + offset
            ax.plot([low, high], [y, y], color=colour, linewidth=2)
            ax.scatter([estimate], [y], color=colour, s=45, zorder=3,
                       label=label if first else None)
            first = False
    ax.axvline(0, color="#555", linewidth=0.9, linestyle="--")
    ax.set_yticks(range(len(others)))
    ax.set_yticklabels([f"{r} vs {REGION_REFERENCE}" for r in others])
    ax.invert_yaxis()
    ax.set_xlabel("Difference in burden (reads per host cell), 95% CI")
    ax.set_title("Region effects with and without individual-level clustering")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, axis="x")
    n_individuals = len({s["individual"] for s in samples})
    _finish(fig, f"F10  burden ~ region + (1|individual); n = {len(samples)} samples from "
                f"{n_individuals} individuals; SD042-18 contributed 3 regions and is retained "
                "(the mixed model does not require balance)", outdir / "F10_region_forest.pdf")

    # Diagnosis as a SECONDARY model, reported in the legend rather than plotted:
    # it is constant within individual, so it competes with the random intercept
    # and is estimated across 9 individuals (6 AD vs 3 control).
    secondary = smf.mixedlm(f"{formula} + diagnosis", data=frame,
                            groups=frame["individual"]).fit()
    rows = []
    for label, fit in (("naive_ols", naive), ("clustered_mixedlm", clustered),
                       ("clustered_plus_diagnosis", secondary)):
        for term in fit.params.index:
            if term in ("Intercept", "Group Var"):
                continue
            ci = fit.conf_int()
            rows.append([label, term, round(float(fit.params[term]), 4),
                         round(float(ci.loc[term][0]), 4), round(float(ci.loc[term][1]), 4),
                         f"{float(fit.pvalues[term]):.3g}"])
    _write_csv(outdir / "F10_region_models.csv",
               ["model", "term", "estimate", "ci_low", "ci_high", "p"], rows)
    return {"naive": naive, "clustered": clustered, "secondary": secondary}


# --- T7 / C5 and T6 -------------------------------------------------------


def table_t7(matrix: dict, samples: list[dict], outdir: Path) -> dict:
    """A4: the fold change is scored; the absolutes are transcription only."""
    controls = [s for s in samples if s["diagnosis"] == "CONTROL"]
    control_mean = sum(s["burden"] for s in controls) / len(controls)
    peak = max(samples, key=lambda s: s["burden"])
    fold = peak["burden"] / control_mean if control_mean > 0 else float("nan")

    # Domain strata for the control arm. Cutoff-free and pooled, so the
    # ordering is computed the same way C1 is -- C5 asks for consistency with
    # C1, and mixing definitions across the two would be the error A4 warns of.
    control_domains = domain_composition(matrix, [s["col"] for s in controls])
    bacteria = control_domains.get("Bacteria", {}).get("per_host_cell", 0.0)
    fungi = control_domains.get("Fungi", {}).get("per_host_cell", 0.0)

    verdict = (
        "PASS" if (5.0 <= fold <= 20.0 and peak["diagnosis"] == "AD" and fungi > bacteria)
        else "PARTIAL" if (2.0 <= fold <= 40.0 and peak["diagnosis"] == "AD")
        or (5.0 <= fold <= 20.0)
        else "MISS"
    )
    _write_csv(
        outdir / "T7_absolute_burden.csv",
        ["quantity", "published", "published_units", "coblast", "coblast_units", "scored"],
        [
            ["Bacteria, control mean", C5_PUBLISHED["bacteria"], "organisms per host cell",
             round(bacteria, 4), "reads per host cell", "NOT SCORED (A4: unit mismatch)"],
            ["Fungi, control mean", C5_PUBLISHED["fungi"], "organisms per host cell",
             round(fungi, 4), "reads per host cell", "NOT SCORED (A4: unit mismatch)"],
            ["Combined, control mean", C5_PUBLISHED["combined"], "organisms per host cell",
             round(control_mean, 4), "reads per host cell", "NOT SCORED (A4: unit mismatch)"],
            ["Maximum-burden case", C5_PUBLISHED["max_case"], "organisms per host cell",
             round(peak["burden"], 4), "reads per host cell", "NOT SCORED (A4: unit mismatch)"],
            ["Max-to-control fold change", round(C5_PUBLISHED_FOLD, 2), "ratio",
             round(fold, 2), "ratio", f"SCORED -> {verdict}"],
            ["Maximum-burden sample", "M66 (AD, male)", "donor",
             f"{peak['sample']} / {peak['individual']} ({peak['diagnosis_full']})", "donor",
             "donor identity not automated -- EBB age/sex is a manual lookup"],
        ],
    )
    return {"fold": fold, "verdict": verdict, "peak": peak, "control_mean": control_mean,
            "bacteria": bacteria, "fungi": fungi}


def table_t6(rows: list[list], outdir: Path) -> None:
    _write_csv(outdir / "T6_scorecard.csv",
               ["target", "published value", "your value", "verdict"], rows)


# --- main -----------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ebb-batch", required=True, help="eToL Full batch over the 35 EBB samples.")
    parser.add_argument("--kohen-batch", required=True, help="eToL Full batch over the 29 Kohen samples.")
    parser.add_argument("--viral-batch", required=True, help="eToL-V batch over the 35 EBB samples.")
    parser.add_argument("--covariates", required=True, help="f10_covariates.csv")
    parser.add_argument("--age-metadata", required=True, help="Kohen age_metadata.csv")
    parser.add_argument("--outdir", default="figures")
    parser.add_argument("--only", nargs="*", default=None,
                        help="Subset, e.g. --only F6 F10. Default: everything.")
    args = parser.parse_args(argv)

    wanted = {name.upper() for name in args.only} if args.only else None

    def run(name: str) -> bool:
        return wanted is None or name in wanted

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    covariates = load_covariates(args.covariates)
    ebb = load_matrix(args.ebb_batch, "etol_full")
    samples = ebb_sample_index(ebb["matrix"], covariates)
    print(f"EBB: {len(samples)} samples, "
          f"{len({s['individual'] for s in samples})} individuals")

    scorecard: list[list] = []

    # V2 -- recomputed and checked against A5 rather than transcribed, so the
    # scorecard cannot quietly diverge from the registration.
    if run("T6") or run("F5"):
        viral = load_matrix(args.viral_batch, "etol_v")
    if run("T6"):
        validated = compute_confusion(viral["matrix"], stage="validated")
        raw = compute_confusion(viral["matrix"], stage="raw")
        for label, got, want in (("validated", validated, A5_V2_VALIDATED),
                                 ("raw", raw, A5_V2_RAW)):
            actual = {k: got[k] for k in ("tp", "fp", "fn", "tn")}
            if actual != want:
                raise SystemExit(
                    f"V2 {label} arm is {actual}, but amendment A5 records {want}. "
                    "Either the batch changed or the registration is stale -- resolve "
                    "before reporting a scorecard."
                )
        scorecard.append([
            "V2 eToL-V end-to-end vs corrected ground truth",
            "TP 9 / FP 1 / FN 36 / TN 794; precision >= .90",
            f"TP {validated['tp']} / FP {validated['fp']} / FN {validated['fn']} / "
            f"TN {validated['tn']}; precision {validated['precision']:.2f}",
            "MISS (precision .40 < .90; A5)",
        ])

    if run("F5"):
        figure_f5(viral, covariates, outdir)

    if run("F6"):
        figure_f6(samples, outdir)

    # C1 -- cutoff-free, pooled, raw hits: the etol_pie.js quantity exactly.
    if run("T6") or run("C1"):
        composition = domain_composition(ebb["matrix"])
        top3 = list(composition)[:3]
        all_seven = len(composition) >= 7
        expected = ["Fungi", "Bacteria", "Chloroplastida"]
        c1 = ("PASS" if top3 == expected and all_seven
              else "PARTIAL" if sorted(top3) == sorted(expected) and all_seven
              else "MISS")
        _write_csv(outdir / "C1_domain_composition.csv",
                   ["domain", "reads", "percent", "reads_per_host_cell", "samples_detected"],
                   [[d, v["reads"], round(v["percent"], 2), round(v["per_host_cell"], 4),
                     v["samples_detected"]] for d, v in composition.items()])
        scorecard.append([
            "C1 domain composition",
            "top 3 = fungi, bacteria, chloroplastida; all 7 present",
            f"top 3 = {', '.join(top3)}; {len(composition)}/7 domains detected", c1,
        ])

    # C2 -- burden by region, with the A2 cutoff applied per species.
    if run("T6") or run("C2"):
        def median(values):
            s = sorted(values)
            n = len(s)
            return float("nan") if not n else (
                s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2)

        by_region = {r: median([s["burden"] for s in samples if s["region"] == r])
                     for r in REGION_ORDER}
        ba24_ad = median([s["burden"] for s in samples
                          if s["region"] == "BA24" and s["diagnosis"] == "AD"])
        ba24_ctrl = median([s["burden"] for s in samples
                            if s["region"] == "BA24" and s["diagnosis"] == "CONTROL"])
        highest = max(by_region, key=lambda r: by_region[r])
        conditions = [highest == "BA24", ba24_ad > ba24_ctrl]
        scorecard.append([
            "C2 regional burden",
            "BA24 highest of 4 regions; AD-BA24 > control-BA24",
            f"highest = {highest}; AD-BA24 {ba24_ad:.2f} vs control-BA24 {ba24_ctrl:.2f}",
            "PASS" if all(conditions) else "PARTIAL" if any(conditions) else "MISS",
        ])

    kohen = None
    if run("F4") or run("F8") or run("T8") or run("T6"):
        ages = load_age_metadata(args.age_metadata)
        kohen = load_matrix(args.kohen_batch, "etol_full")
        kohen_rows = kohen_burdens(kohen["matrix"], ages)

    if run("F4"):
        figure_f4(kohen_rows, outdir)

    if run("F8") or run("T8"):
        taxon_results = kohen_taxon_stats(kohen["matrix"], ages)
        figure_f8_and_t8(taxon_results, len(kohen["matrix"]["rows"]), outdir)

    # C4 -- the same Welch F4 annotates, so figure and verdict cannot disagree.
    if run("T6"):
        result = kohen_age_test(kohen_rows)
        scorecard.append([
            "C4 burden increases with age (Kohen)",
            f"significant increase, P = {KOHEN_PUBLISHED_P}, t-test",
            f"elderly {result['mean_a']:.2f} vs young {result['mean_b']:.2f}, "
            f"Welch p = {result['p']:.4f}",
            result["verdict"],
        ])

    entries = None
    if run("F9") or run("T6"):
        entries = figure_f9(ebb, samples, outdir)
        matched = sum(1 for e in entries if e["over"])
        scorecard.append([
            "C3 AD-shortlist over-representation",
            "at least 6 of 10 genera detected and over-represented in AD",
            f"{matched}/10 (tenth entry matched at domain level, A3)",
            "PASS" if matched >= 6 else "PARTIAL" if matched >= 4 else "MISS",
        ])

    if run("T7") or run("T6"):
        t7 = table_t7(ebb["matrix"], samples, outdir)
        scorecard.append([
            "C5 absolute burden (max-to-control fold change)",
            f"{C5_PUBLISHED_FOLD:.1f}x (1.8 / 0.19 organisms per host cell)",
            f"{t7['fold']:.1f}x; max case {t7['peak']['individual']} "
            f"({t7['peak']['diagnosis_full']})",
            t7["verdict"],
        ])

    if run("F10"):
        figure_f10(samples, outdir)

    if run("T6") and scorecard:
        order = {"V2": 0, "C1": 1, "C2": 2, "C3": 3, "C4": 4, "C5": 5}
        scorecard.sort(key=lambda row: order.get(row[0].split()[0], 99))
        table_t6(scorecard, outdir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
