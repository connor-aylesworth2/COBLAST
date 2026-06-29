"""eToL-V validation: confusion matrix vs the eToL WGS ground truth.

The eToL-V dissertation (Edinburgh B270917) validated its viral calls against
whole-genome-shotgun (WGS) results from the original eToL workflow, summarised as
a 2x2 confusion matrix (her Figure 9: TP=9, FP=1, FN=35, TN=411 -> accuracy 92%,
precision 90%, recall 20%, F1 0.3). This module reproduces that comparison for a
COBLAST+ eToL-V batch.

Data (bundled, reference study SRP398685 / Edinburgh Brain Bank, 35 samples):
* ``data/etol_v_wgs_truth.csv``    - ``srx,virus,count`` WGS read counts.
* ``data/etol_v_sra_crosswalk.csv`` - ``srr,srx,region,diagnosis,sample_name``;
  the SRR<->SRX map is INVERTED (SRX17674433<->SRR21676133), so it is verified,
  not assumed.

Construction (reverse-engineered from her dissertation and verified to reproduce
9/1/35/411 exactly):
* Binary, ``present = WGS count > 0`` (she states a "binary classification model";
  >0 is the only threshold that yields her 44 actual-positives).
* eToL-V "present" = any of a virus's probes has a *validated* (contig-confirmed)
  hit > 0 -- i.e. compare to her post-validation heatmap (Fig 10), not the raw net.
* Universe = a fixed set of 13 WGS virus rows x the 35 samples (= 455), plus any
  out-of-universe validated prediction as an extra FP cell (her single HPV45 hit
  -> 456). SARS-CoV/-CoV-2 are excluded entirely (no WGS data).
* The 13-virus universe carries two documented quirks kept for fidelity: it
  includes HPV6 (the panel has NO HPV6 probe, so it is always predicted negative)
  and omits Adenovirus A / Adenovirus 54. Swap ``VESO_UNIVERSE`` for a corrected
  set if you would rather not reproduce those.
"""

from __future__ import annotations

import csv
from collections import OrderedDict
from functools import lru_cache
from typing import Any

from config import resource_path
from etol_summary import etol_preset_records


WGS_TRUTH_PATH = resource_path("data", "etol_v_wgs_truth.csv")
SRA_CROSSWALK_PATH = resource_path("data", "etol_v_sra_crosswalk.csv")

# WGS virus row -> the eToL-V virus token (the first part of a taxon's species,
# e.g. ``V-HAdV_AdC_penton`` -> ``AdC``). ``None`` means the virus is in the
# scoring universe but the panel has no probe for it (always predicted negative).
VESO_UNIVERSE: "OrderedDict[str, str | None]" = OrderedDict(
    [
        ("Adenovirus C", "AdC"),
        ("COV_229E", "HCoV229E"),
        ("HHV1_HSV1", "HSV1"),
        ("HHV2_HSV2", "HSV2"),
        ("HHV3_VZV", "VZV"),
        ("HHV4_EBV", "EBV"),
        ("HHV5_CMV", "CMV"),
        ("HHV6A", "HHV6A"),
        ("HHV6B", "HHV6B"),
        ("HHV7", "HHV7"),
        ("HHV8", "KSHV"),
        ("HPV6", None),
        ("HPV16", "HPV16"),
    ]
)

# Virus tokens dropped from scoring even when predicted (no corresponding WGS data).
EXCLUDED_VIRUS_TOKENS = frozenset({"SARSCoV2", "SARSCoV"})


def _virus_token(taxon_species: str) -> str:
    """First token of a taxon's species label is the virus (``AdC_penton`` -> ``AdC``)."""
    return taxon_species.split("_", 1)[0]


@lru_cache(maxsize=1)
def _taxa_by_virus_token() -> dict[str, frozenset[str]]:
    """Map each eToL-V virus token to the set of taxa (probe groups) that detect it."""
    grouped: dict[str, set[str]] = {}
    for record in etol_preset_records("etol_v"):
        grouped.setdefault(_virus_token(record["species"]), set()).add(record["taxon"])
    return {token: frozenset(taxa) for token, taxa in grouped.items()}


def universe_taxa(universe: "OrderedDict[str, str | None]" = VESO_UNIVERSE) -> dict[str, frozenset[str]]:
    """WGS virus name -> the taxa whose validated hits mark it present (may be empty)."""
    by_token = _taxa_by_virus_token()
    return {
        virus: (by_token.get(token, frozenset()) if token else frozenset())
        for virus, token in universe.items()
    }


def load_wgs_truth(path: Any = WGS_TRUTH_PATH) -> dict[tuple[str, str], int]:
    """Load ``(srx, virus) -> count`` from the bundled WGS ground-truth CSV."""
    truth: dict[tuple[str, str], int] = {}
    text = (path.read_text(encoding="utf-8") if hasattr(path, "read_text")
            else open(path, encoding="utf-8").read())
    for row in csv.DictReader(text.splitlines()):
        try:
            count = int(float(row.get("count", "0") or 0))
        except ValueError:
            count = 0
        truth[(row["srx"].strip(), row["virus"].strip())] = count
    return truth


def load_crosswalk(path: Any = SRA_CROSSWALK_PATH) -> dict[str, str]:
    """Load a sample-id -> SRX map (both SRR and SRX keys resolve to the SRX)."""
    text = (path.read_text(encoding="utf-8") if hasattr(path, "read_text")
            else open(path, encoding="utf-8").read())
    alias: dict[str, str] = {}
    for row in csv.DictReader(text.splitlines()):
        srx = row["srx"].strip()
        alias[row["srr"].strip()] = srx
        alias[srx] = srx  # an SRX-labelled run joins directly
    return alias


def _empty_metrics() -> dict[str, Any]:
    return {
        "tp": 0, "fp": 0, "fn": 0, "tn": 0, "n": 0,
        "accuracy": None, "precision": None, "recall": None, "f1": None,
        "scored_samples": 0, "stage": "validated", "cells": [],
        "unmatched_samples": [],
    }


def _finalize(m: dict[str, Any]) -> dict[str, Any]:
    tp, fp, fn, tn = m["tp"], m["fp"], m["fn"], m["tn"]
    m["n"] = tp + fp + fn + tn
    total = m["n"]
    m["accuracy"] = (tp + tn) / total if total else None
    m["precision"] = tp / (tp + fp) if (tp + fp) else None
    m["recall"] = tp / (tp + fn) if (tp + fn) else None
    if m["precision"] and m["recall"]:
        m["f1"] = 2 * m["precision"] * m["recall"] / (m["precision"] + m["recall"])
    else:
        m["f1"] = 0.0 if total else None
    return m


def compute_confusion(
    matrix: dict[str, Any],
    truth: dict[tuple[str, str], int] | None = None,
    crosswalk: dict[str, str] | None = None,
    *,
    stage: str = "validated",
    universe: "OrderedDict[str, str | None]" = VESO_UNIVERSE,
) -> dict[str, Any]:
    """Confusion matrix of an eToL-V batch's calls vs the WGS ground truth.

    ``matrix`` is a :func:`etol_summary.build_etol_matrix` payload (``level=
    "species"``). ``stage="validated"`` scores the contig-confirmed layer (the
    faithful comparison); ``"raw"`` scores the net hits. Samples are mapped to an
    SRX via ``crosswalk`` (so SRR-labelled runs join the SRX-keyed truth); samples
    with no crosswalk entry or no truth are reported in ``unmatched_samples`` and
    skipped.
    """
    truth = load_wgs_truth() if truth is None else truth
    crosswalk = load_crosswalk() if crosswalk is None else crosswalk
    taxa_for = universe_taxa(universe)

    counts = matrix.get("confirmed") if stage == "validated" else matrix.get("hits")
    if counts is None:  # no validated layer available -> fall back to raw net hits
        counts = matrix.get("hits")
        stage = "raw"

    row_index = {row["key"]: i for i, row in enumerate(matrix.get("rows", []))}
    truth_srx = {srx for (srx, _virus) in truth}

    scored_cols: list[tuple[int, str]] = []
    unmatched: list[str] = []
    for j, col in enumerate(matrix.get("cols", [])):
        srx = crosswalk.get(col["sample"])
        if srx and srx in truth_srx:
            scored_cols.append((j, srx))
        else:
            unmatched.append(col["sample"])

    def predicted(taxa: frozenset[str], col: int) -> bool:
        return any(
            (counts[row_index[t]][col] or 0) > 0
            for t in taxa if t in row_index
        )

    m = _empty_metrics()
    m["stage"] = stage
    m["scored_samples"] = len(scored_cols)
    m["unmatched_samples"] = unmatched

    # In-universe cells: 13 viruses x scored samples.
    for virus, taxa in taxa_for.items():
        for col, srx in scored_cols:
            actual = truth.get((srx, virus), 0) > 0
            pred = predicted(taxa, col)
            if actual and pred:
                m["tp"] += 1
            elif actual and not pred:
                m["fn"] += 1
            elif pred and not actual:
                m["fp"] += 1
            else:
                m["tn"] += 1

    # Out-of-universe validated predictions (e.g. HPV45) are extra FP cells, one
    # per (virus token, sample); SARS tokens are dropped entirely.
    universe_tokens = {tok for tok in universe.values() if tok}
    by_token = _taxa_by_virus_token()
    for token, taxa in by_token.items():
        if token in universe_tokens or token in EXCLUDED_VIRUS_TOKENS:
            continue
        for col, _srx in scored_cols:
            if predicted(taxa, col):
                m["fp"] += 1

    return _finalize(m)
