"""eToL-V confusion-matrix validation tests.

The headline test reconstructs the eToL-V dissertation's Figure 9 from the bundled
WGS ground truth: with a validated-hit matrix that mirrors her surviving calls
(adenovirus C penton in nine WGS-positive samples + the lone HPV45 L1 hit), the
computation must reproduce her published TP=9 / FP=1 / FN=35 / TN=411.
"""

from etol_summary import etol_preset_records
from etol_validation import (
    VESO_UNIVERSE,
    compute_confusion,
    load_crosswalk,
    load_wgs_truth,
    universe_taxa,
)


def test_crosswalk_is_inverted_and_complete():
    cross = load_crosswalk()
    # 35 SRR aliases (+ their SRX self-aliases).
    srr_keys = {k for k in cross if k.startswith("SRR")}
    assert len(srr_keys) == 35
    # The SRR<->SRX mapping is inverted, not same-suffix.
    assert cross["SRR21676133"] == "SRX17674433"
    assert cross["SRR21676131"] == "SRX17674435"
    assert cross["SRX17674433"] == "SRX17674433"  # SRX joins directly


def test_wgs_truth_positive_cells_match_actual_positives():
    truth = load_wgs_truth()
    # Exactly the 13-virus universe carries 44 WGS-positive cells (= TP+FN).
    pos = sum(
        1
        for (srx, virus), count in truth.items()
        if virus in VESO_UNIVERSE and count > 0
    )
    assert pos == 44


def test_universe_taxa_quirks():
    taxa = universe_taxa()
    # HPV6 is in the universe but the panel has no probe for it.
    assert taxa["HPV6"] == frozenset()
    # Adenovirus C maps to its three structural-protein taxa.
    assert taxa["Adenovirus C"] == frozenset(
        {"V-HAdV_AdC_penton", "V-HAdV_AdC_hexon", "V-HAdV_AdC_fiber"}
    )


def _zero_matrix():
    """A build_etol_matrix-shaped payload: all eToL-V taxa x the 35 WGS samples."""
    records = etol_preset_records("etol_v")
    taxa = []
    seen = set()
    for r in records:
        if r["taxon"] not in seen:
            seen.add(r["taxon"])
            taxa.append(r["taxon"])
    cross = load_crosswalk()
    srx_to_srr = {v: k for k, v in cross.items() if k.startswith("SRR")}
    srx_order = sorted(srx_to_srr)
    cols = [{"sample": srx_to_srr[srx]} for srx in srx_order]
    confirmed = [[0 for _ in cols] for _ in taxa]
    rows = [{"key": t} for t in taxa]
    return {"rows": rows, "cols": cols, "confirmed": confirmed,
            "hits": [list(r) for r in confirmed]}, taxa, srx_order


def test_reproduces_veso_confusion_matrix():
    truth = load_wgs_truth()
    matrix, taxa, srx_order = _zero_matrix()
    row_of = {t: i for i, t in enumerate(taxa)}
    col_of = {srx: j for j, srx in enumerate(srx_order)}

    # eToL-V's surviving validated calls: adenovirus C penton in 9 of the 30
    # WGS adenovirus-C-positive samples (-> 9 TP, 35 FN) ...
    adc_pos = [srx for srx in srx_order if truth.get((srx, "Adenovirus C"), 0) > 0]
    for srx in adc_pos[:9]:
        matrix["confirmed"][row_of["V-HAdV_AdC_penton"]][col_of[srx]] = 1
    # ... and the single HPV45 L1 hit in SRX17674444 (-> the 1 FP).
    matrix["confirmed"][row_of["V-HPV_HPV45_L1"]][col_of["SRX17674444"]] = 1

    m = compute_confusion(matrix, truth=truth, stage="validated")
    assert (m["tp"], m["fp"], m["fn"], m["tn"]) == (9, 1, 35, 411)
    assert m["n"] == 456
    assert round(m["accuracy"], 4) == 0.9211
    assert round(m["precision"], 2) == 0.90
    assert round(m["recall"], 4) == 0.2045
    assert round(m["f1"], 2) == 0.33


def test_sars_predictions_are_excluded():
    truth = load_wgs_truth()
    matrix, taxa, srx_order = _zero_matrix()
    row_of = {t: i for i, t in enumerate(taxa)}
    # A SARS-CoV-2 validated hit must NOT count as a false positive (no WGS data).
    matrix["confirmed"][row_of["V-HCoV_SARSCoV2_S"]][0] = 5
    m = compute_confusion(matrix, truth=truth, stage="validated")
    assert m["fp"] == 0
    assert m["tp"] == 0
