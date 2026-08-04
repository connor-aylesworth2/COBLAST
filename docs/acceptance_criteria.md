---
editor_options: 
  markdown: 
    wrap: 72
---

# Pre-registered acceptance criteria per validation target

**Registered 2026-07-30 against commit `a7cf9a8`, before any scored
run.** Committing this file *is* the registration: the git timestamp is
the evidence that the criteria predate the results. Do not edit a
criterion after a run — add a dated amendment row instead.

Scope: **validation** only (does COBLAST+ reproduce a published
result?). Verification artefacts — the spike-in positive control, the
frozen self-check, the 150-test suite — belong in §3.4 and are
deliberately absent here.

Prior run disclosed: an eToL-V shakedown on 2026-06-29 scored TP 2 / FP
0 / FN 42 / TN 411 under the retired 13-virus universe. It is reported
as an uncontrolled pilot, not as a result against these criteria.

## Reference works

-   **Hu et al. 2022**, *BMC Microbiol* 22:317 — the eToL method paper.
    Cohorts: MSBB, Rockefeller, Miami. **Not** EBB.
-   **Hu, McKenzie, Smith, Haas & Lathe 2023**, bioRxiv
    10.1101/2023.02.06.527297, "The remarkable complexity of the brain
    microbiome in health and disease" — the **EBB paper**: four regions
    (AMYG, BA24, HPC, HYPO) from 3 control and 6 AD brains, one sample
    compromised, **35 samples**. This is the same cohort as
    `data/etol_v_sra_crosswalk.csv`, so C1–C3 below are same-data
    reproductions. AD/VaD and AD/LBD individuals were retained in the
    broad AD category → sample level **24 AD vs 11 control**.
-   **eToL-V dissertation** (Edinburgh B270917) — the viral
    port-fidelity target.

| \# | Target | Claim type | Reference value (frozen) | Scored from | Pass | Partial | Miss |
|----|----|----|----|----|----|----|----|
| V1 | Reference confusion matrix recomputed from its own ground truth (eToL-V, Fig. 9) | Value-level, exact | TP 9 / FP 1 / FN 35 / TN 411, N = 456; acc .9211, prec .90, rec .2045, F1 .33; 13-virus `VESO_UNIVERSE` × 35 | `compute_confusion(truth=load_wgs_truth(), universe=VESO_UNIVERSE)`; pinned by `tests/test_etol_validation.py` | All four cells match exactly | — | Any cell differs |
| V2 | End-to-end eToL-V run, 35 EBB samples, vs the current ground truth | Value-level, banded | 24 viruses × 35 = 840 cells, 45 positive (HAdV-C 31, HSV1 7, CMV 5, EBV 1, KSHV 1). Reference-equivalent score: **TP 9 / FP 1 / FN 36 / TN 794**, acc .956, prec .90, rec .20, F1 .33 | `/batch-results/<id>/etol-confusion.csv`, `stage="validated"`; `stage="raw"` as sensitivity arm | TN ≥ 794 **and** precision ≥ .90 **and** recall ≥ .20 **and** every TP an HAdV-C cell | Same but recall ≥ .10, **and** each FN attributed to a stage via `raw_hits`/`confirmed_hits` | Precision \< .90, or a TP outside HAdV-C, or FN causes unattributed |
| C1 | Domain composition of brain microbiome (2023, Fig. 1 / Fig. 4A) | Ordinal + set-level | Most abundant three domains: **fungi, bacteria, chloroplastida**; archaea, amoebozoa, basal eukaryota and holozoa/metazoa also detected (all 7 present) | Domain composition CSV (`static/etol_pie.js` export), reads per host cell, 35 EBB samples | Top 3 domains by reads per host cell are exactly {F, B, C} in that order, and all 7 domains have ≥1 detected sample | Top 3 set is {F, B, C} in a different order, all 7 present | A domain outside {F, B, C} enters the top 3, or a domain is wholly undetected |
| C2 | Regional burden across limbic brain (2023, Fig. 4B) | Ordinal + direction | Microbial burden **highest in cingulate cortex (BA24)** of the four regions, and **higher in BA24 of AD than control**. The paper reports **no statistics** for this figure ("statistical analysis was not undertaken") | Species Summary CSV, `Reads per host cell` summed per sample, grouped by region from the crosswalk | BA24 has the highest median burden of the 4 regions **and** AD-BA24 median \> control-BA24 median | One of the two holds | Neither holds |
| C3 | Species over-represented in AD brain, EBB arm (2023, Figs. 8C, 9A) | Set-level + direction | The published top-10 shortlist: fungi **Cortinarius** (highest absolute abundance; 2nd most differential), **Tausonia** (Cryptococcus group), **Acrocalymma**, **Aureobasidium**, **Alternaria** (Aspergillus-like Ascomycota), **Komagataella** (Candida group); bacteria **Sphingomonas** (or Ralstonia; highest-abundance bacterium), **Streptococcus**, **Staphylococcus** (or Bacillus); plus one uncharacterised chloroplastida. Over-representation reached up to 60 % of AD samples; the ranked top 10 spans 47 % → 33 % (mean of MSBB + EBB) | Probe Counts CSV → per-probe AD-sample-over-control-mean proportion (computed downstream); genus names from `Closest homolog (contig)` in the Species Summary CSV | ≥ 6 of the 10 genera detected in EBB **and** over-represented in AD (proportion of AD samples above the EBB control mean \> 0.5 × control proportion) | 4–5 of 10 | ≤ 3 of 10 |
| C4 *(conditional)* | Microbial burden increases with age (2023, Fig. 5) | Direction + significance | Human hippocampus, Kohen et al. 2014, **N = 29** (young N = 16, mean 47.6 y, range 29–59; elderly N = 13, mean 85.2 y, range 68–95): significant increase with age, **P = 0.0202**, *t*-test | Species Summary CSV, total reads per host cell per sample, two-group comparison | Direction is an increase **and** p \< 0.05 by the same *t*-test family (Welch, `parametricTest`) | Direction is an increase, p ≥ 0.05 | Direction reversed |

C4 runs only if the Kohen 2014 run accessions are obtained; nothing
about that dataset exists in this repository. If they are not obtained
before the registration date, strike the row rather than leaving it
open.

**Fixed run configuration for V2, C1–C3** (any deviation invalidates the
score): eToL Full preset (eToL-V for V2), net gate E \< 0.01, human
filter on against `nt_human_9606` at 150 bits, contig assembly on,
contig identification on (`ToL_rRNA` for C1–C3, `ToL_virus_val` for V2),
re-probing **off**, cellular cutoff **4 reads per host cell** (the
paper's range is 3–5, typically 4–5 — one value is chosen and fixed
here), uniform across all 35 samples.

**Design matrix must be collapsed to two conditions.** The paper keeps
AD/VaD and AD/LBD inside the broad AD category (24 AD vs 11 control).
The crosswalk's four-way labelling (AD 16, CONTROL 11, AD/VaD 4, AD/LBD
4) would test a contrast the paper never drew, and the per-domain tests
require ≥ 3 normalizable samples per condition. Ship a two-level design
matrix for the reproduction runs and keep the four-level one for the
application chapter only.

**Do not add statistics the source did not report.** C2's source figure
carries no test; reproduction is the rank order, not a p-value. Note
also that at the *individual* level EBB is 6 AD vs 3 control, and
`rankTestPlan`'s own floor says a 3-versus-3 design over 7 domains
cannot reach q \< 0.05 for any domain — so a clustered test here is
guaranteed null by construction. Clustered analysis belongs in §3.6 as
an application-leg limitation, not as a reproduction criterion.

## Two cells you must resolve before this is registrable

1.  **V2 — how SARS-CoV-2 predictions score.** The current ground truth
    has *zero* SARS-CoV-2 positive cells, but `SARS-CoV-2` is a truth
    column, so `compute_confusion` scores it in-universe: every
    validated SARS-CoV-2 call becomes an FP and pushes precision below
    the .90 bar. The reference's headline surviving signal included
    SARS-CoV-2. `EXCLUDED_VIRUS_TOKENS` drops those tokens only on the
    legacy path. Decide in writing: score them (and set the
    reference-equivalent FP count from the reference's own SARS-CoV-2
    sample list) or exclude them. The V2 value above assumes
    **excluded**.
2.  **C3 — inherited shortlist, not re-derived.** The paper's top-10 was
    selected from a *pooled* four-dataset ranking; you cannot re-derive
    that selection from EBB alone. The criterion therefore inherits the
    published shortlist and asks only whether those genera reproduce in
    EBB. If you want per-genus EBB percentages as reference values, they
    must be read off Figure 8C's EBB (yellow) bars — the text gives only
    the 47 %→33 % pooled range. Absent that, keep the criterion at
    set-overlap and direction, which the text does support.

## Excluded from validation, with reason

| Published result | Why it cannot be a criterion |
|----|----|
| Pooled differential abundance: top-100 mean overabundance **13.68** (range 2.23–181.97, SD 20.40, P = 0.01); cluster mean **12.59** (range 2.52–36.88, SD 7.83, P \< 0.001) | Computed over MSBB + EBB + Miami + Rockefeller (31 control, 48 AD). MSBB is access-gated via the AD Knowledge Portal; EBB alone cannot reproduce the pooled statistic |
| The 23S/28S (LSU) confirmation arm behind Figs. 8–9 | COBLAST+ has no LSU disambiguation/reprobe pass — the documented gap in §3.5.3. `ToL_rRNA` contains LSU, but nothing uses it this way. Future Work |
| Viral findings: AdC = **83 %** of viral transcripts, detected in **\~50 %** of samples; viral burden higher in AD **P = 0.04**; HHV↔AdC inverse relationship in AD only; cutoff **0.03** reads per host cell | Different instrument: the paper probes with *stripped whole genomes* of the top-20 Readhead viruses, not eToL-V's 120-mer structural-protein panel. A disagreement here is a method difference, not a port defect. Cite as context for V2's recall, never as a target |
| Endogenous retroelement / HERV abundances | The extended retroelement probe list is not in COBLAST+ |
| Brain microbiome ≈ **20 %** of gut microbiome diversity; gut \~4× more diverse | Needs gut/faecal RNA-seq; no host normalisation possible |
| Cross-species conservation from *Drosophila* to human (Fig. 2) | Needs non-human brain datasets and has no host-cell normalisation |
| Microbe/retroelement/virus correlations R = 0.394, 0.267, 0.139 | Depends on the retroelement and stripped-genome arms above. (The source labels these both "coefficients of determination" and "R2 values"; do not quote a value without resolving which it is) |
| Hu et al. 2022 cohort results (MSBB, Rockefeller, Miami) | Different cohorts, access-gated; the 2023 paper supersedes it as the same-data target |

## Not in this table, on purpose

Anything the tool cannot currently export. A criterion that needs code
written to score it is not pre-registered, it is a plan.
