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
frozen self-check, the 151-test suite — belong in §2.4 and are
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
|---------|---------|---------|---------|---------|---------|---------|---------|
| V1 | Reference confusion matrix recomputed from its own ground truth (eToL-V, Fig. 9) | Value-level, exact | TP 9 / FP 1 / FN 35 / TN 411, N = 456; acc .9211, prec .90, rec .2045, F1 .33; 13-virus `VESO_UNIVERSE` × 35 | `compute_confusion(truth=load_wgs_truth(), universe=VESO_UNIVERSE)`; pinned by `tests/test_etol_validation.py` | All four cells match exactly | — | Any cell differs |
| V2 | End-to-end eToL-V run, 35 EBB samples, vs the current ground truth | Value-level, banded | 24 viruses × 35 = 840 cells, 45 positive (HAdV-C 31, HSV1 7, CMV 5, EBV 1, KSHV 1). Reference-equivalent score: **TP 9 / FP 1 / FN 36 / TN 794**, acc .956, prec .90, rec .20, F1 .33 | `/batch-results/<id>/etol-confusion.csv`, `stage="validated"`; `stage="raw"` as sensitivity arm | TN ≥ 794 **and** precision ≥ .90 **and** recall ≥ .20 **and** every TP an HAdV-C cell | Same but recall ≥ .10, **and** each FN attributed to a stage via `raw_hits`/`confirmed_hits` | Precision \< .90, or a TP outside HAdV-C, or FN causes unattributed |
| C1 | Domain composition of brain microbiome (2023, Fig. 1 / Fig. 4A) | Ordinal + set-level | Most abundant three domains: **fungi, bacteria, chloroplastida**; archaea, amoebozoa, basal eukaryota and holozoa/metazoa also detected (all 7 present) | Domain composition CSV (`static/etol_pie.js` export), reads per host cell, 35 EBB samples | Top 3 domains by reads per host cell are exactly {F, B, C} in that order, and all 7 domains have ≥1 detected sample | Top 3 set is {F, B, C} in a different order, all 7 present | A domain outside {F, B, C} enters the top 3, or a domain is wholly undetected |
| C2 | Regional burden across limbic brain (2023, Fig. 4B) | Ordinal + direction | Microbial burden **highest in cingulate cortex (BA24)** of the four regions, and **higher in BA24 of AD than control**. The paper reports **no statistics** for this figure ("statistical analysis was not undertaken") | Species Summary CSV, `Reads per host cell` summed per sample, grouped by region from the crosswalk | BA24 has the highest median burden of the 4 regions **and** AD-BA24 median \> control-BA24 median | One of the two holds | Neither holds |
| C3 | Species over-represented in AD brain, EBB arm (2023, Figs. 8C, 9A) | Set-level + direction | The published top-10 shortlist: fungi **Cortinarius** (highest absolute abundance; 2nd most differential), **Tausonia** (Cryptococcus group), **Acrocalymma**, **Aureobasidium**, **Alternaria** (Aspergillus-like Ascomycota), **Komagataella** (Candida group); bacteria **Sphingomonas** (or Ralstonia; highest-abundance bacterium), **Streptococcus**, **Staphylococcus** (or Bacillus); plus one uncharacterised chloroplastida (identity resolved in **amendment A3**). Over-representation reached up to 60 % of AD samples; the ranked top 10 spans 47 % → 33 % (mean of MSBB + EBB) | Probe Counts CSV → per-probe AD-sample-over-control-mean proportion (computed downstream); genus names from `Closest homolog (contig)` in the Species Summary CSV | ≥ 6 of the 10 genera detected in EBB **and** over-represented in AD (proportion of AD samples above the EBB control mean \> 0.5 × control proportion) | 4–5 of 10 | ≤ 3 of 10 |
| C4 *(conditional)* | Microbial burden increases with age (2023, Fig. 5) | Direction + significance | Human hippocampus, Kohen et al. 2014, **N = 29** (young N = 16, mean 47.6 y, range 29–59; elderly N = 13, mean 85.2 y, range 68–95): significant increase with age, **P = 0.0202**, *t*-test | Species Summary CSV, total reads per host cell per sample, two-group comparison | Direction is an increase **and** p \< 0.05 by the same *t*-test family (Welch, `parametricTest`) | Direction is an increase, p ≥ 0.05 | Direction reversed |
| C5 | Absolute microbial burden per host cell (2023, Discussion; feeds **T7**) | Ratio + direction, banded | **Published values, as stated (organism units):** bacteria ≈ **0.14** per control host cell; fungi ≈ **0.05** per control host cell; the two together ≈ **0.19**; maximum-burden case **1.8 microbes per host cell** (AD, male **M66**, Fig. 4 — an EBB donor). The 0.14/0.05 pair is the 2022 method paper's estimate, *confirmed* in 2023, not an EBB-derived measurement. **Scored quantity: the max-to-control fold change, ≈ 9.5×** (1.8 ÷ 0.19); the absolutes are transcribed but not scored — see amendment **A4** | Species Summary CSV, `Reads per host cell` summed per sample with the A2 cutoff applied per species **first**; `Domain` supplies the bacteria and fungi strata; control stratum = the 11 CONTROL samples; max case = the single highest-burden sample of the 35 | Max-burden sample ÷ control-mean burden lies in **[5×, 20×]**, **and** the max-burden sample is an AD sample, **and** fungi exceed bacteria in reads per host cell across controls (consistent with C1) | Fold change in [2×, 40×] with the max case an AD sample, **or** fold change in [5×, 20×] with the max case a control sample | Fold change \< 2× or \> 40×, or the control domain ordering contradicts C1 |

C4 runs only if the Kohen 2014 run accessions are obtained; nothing
about that dataset exists in this repository. If they are not obtained
before the registration date, strike the row rather than leaving it
open.

**Fixed run configuration for V2, C1–C3, C5** (any deviation invalidates
the score): eToL Full preset (eToL-V for V2), net gate E \< 0.01, human
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
guaranteed null by construction. Clustered analysis belongs in §2.6 as
an application-leg limitation, not as a reproduction criterion.

## Amendments

The registration above is preserved as first written. Each amendment is
dated, states what changed and why, and was recorded **before** the run
it affects was scored. A1–A4 were authored on 2026-08-01 and committed on
2026-08-04; that gap is transcription into this file, not revision of the
reasoning, and no scored run of C3, C4 or C5 had been performed on either
date.

**A5 is the exception, and says so in its own heading.** It was recorded
*after* the eToL-V run it concerns. It changes no Pass/Partial/Miss band
and reverses no outcome; V2's Miss is recorded as the registered bands
produced it.

### A1 · 2026-08-01 — human-filter threshold scales with read length

**Affects C4 (Kohen) only. V2 and C1–C3 are unchanged.**

The fixed run configuration specifies the human filter at **150 bits**,
Hu, Haas & Lathe 2022's brain value. That value is correct for EBB and
unusable for Kohen, for a reason that is arithmetic rather than
judgement.

`find_human_read_ids` runs `blastn -task megablast`
(`human_filter.py:182`) under default scoring (reward 1, penalty −2),
yielding \~1.847 bits per perfectly matching base (λ = 1.28, K = 0.46,
ungapped). A read of length *L* therefore cannot exceed \~1.847·*L*
bits, and a 150-bit gate requires ≥ 81 perfectly matching bases:

| Dataset | Mean read length | Max attainable bitscore | 150-bit gate |
|----|----|----|----|
| EBB (SRP398685) | \~150 bp | 278 | fires; 81/150 = 54 % of the read |
| Kohen 2014 | \~50 bp | 93.5 | **cannot fire — 0 reads removed** |

At 50 bp the filter is a silent no-op: every probe-matched human read is
retained. Because the panel is rRNA and human rRNA is conserved against
microbial rRNA, those reads concentrate on the conserved probes rather
than distributing as noise, so burden inflates non-randomly. Scoring C4
under the unamended configuration would test the age hypothesis on
contaminated counts.

**Amendment: the threshold is set to the library's mean read length.**
Converting each cutoff the source reports into matched bases recovers
that rule from the source itself:

| Dataset     | Published cutoff | Matched bases | Implied read length           |
|----------------|----------------|----------------|--------------------------|
| MSBB        | \> 160           | 86.0          | \~160 bp                      |
| brain / EBB | \> 150           | 80.6          | \~150 bp (measured; confirms) |
| Rockefeller | \> 126           | 67.6          | \~126 bp                      |
| Miami       | \> 100           | 53.5          | \~100 bp                      |

The matched-base fraction is \~54 % across all four datasets, and EBB's
measured 150 bp independently confirms the brain value. **C4 is scored
with the human filter at 50 bits** — 50 bp × 1 bit per base of read
length, giving \~26 matched bases, 53 % of the read, the same fraction
the source applied to every one of its own cohorts.

This restores the source method's intent rather than departing from it,
but it is a deviation from the value registered above and is recorded as
one. It is also a fifth entry on F2's *deviating* track.

Registered before the Kohen run: no C4 number existed when this was
written.

### A2 · 2026-08-01 — cutoff application order made explicit

**Affects C2, C4, C5 and T7.**

The fixed run configuration sets a cellular cutoff of 4 reads per host
cell but does not say where it applies when burden is a sum. Filtering
species below 4 and then summing is not the same number as summing and
then filtering.

**Amendment: the cutoff applies per species, before summing.** A species
contributes to a sample's burden only if its own reads per host cell ≥
4. This is what a per-cell detection threshold means, and it is applied
identically wherever burden is computed.

Registered before any burden number was computed.

### A3 · 2026-08-01 — identity of the tenth AD-shortlist entry, and its matching level

**Affects C3 and F9.**

C3 inherits the paper's top-10 shortlist. Nine entries are named genera.
The tenth is not. The abstract calls it "an uncharacterized
chloroplastida (algae-related) species"; the Results section adds that
"An unknown Chloroplastida species was also detected that is possibly
related to the green alga, *Nephroselmis*". No genus is assigned to it
anywhere in the paper, and Figure 9A's identification column leaves it
unnamed. It is the tenth and lowest-ranked entry, at the 33 % end of the
47 % → 33 % over-representation range.

**Nephroselmis is not a member of the eToL panel.** All 29
Chloroplastida-group taxa were enumerated from `data/eToL_probes.fasta`
(class codes C1–C4, 258 of the 1,017 probes): C1 *Cmucronatum, Jlibera,
Pbrachykentron, Ppyriformis, Tvaginalis*; C2 *GlaucocystisARP2014,
Hglobosa, Igalbana, Pbrassicae, Pfalciparum, Pgallinaceum, Sminus*; C3
*Ccaeruleus, Ccaldarium, Cmerolae, Gsublittoralis, Gsulphuraria,
Paerugineum, Ppurpureum, Rviolacea*; C4 *Aechinata, Catmophyticus,
Cvulgaris, Ebilobata, Gavonlea, Glongispicula, Gtheta, Lmarina,
Pbrevispinosa*. None is *Nephroselmis* and none is a prasinophyte. The
paper reached that name through its **23S/28S LSU confirmation arm** —
the single pass COBLAST+ does not implement, and F2's load-bearing
absence. A probe-name or genus-name match can therefore never succeed,
and C3 would be capped at 9/10 by construction rather than by evidence.

**Amendment: the tenth entry is scored at domain level.** It counts as
matched if **any** Chloroplastida taxon (class code C1–C4) is detected in
EBB and over-represented in AD under the same rule applied to the nine
named genera. The matching level for this entry alone is domain; the
other nine remain genus/clade matches, per the source's own caution
against species-level identification.

The looser level is recorded, not hidden. Group C spans red algae,
glaucophytes and excavates as well as green algae, so a C-group match is
weaker evidence than a genus match, and **F9's legend must state that the
tenth entry was matched at domain level and why** — that the source
identified it by a method the port does not implement. Reporting *which*
C-group taxon carried the signal, and its `Closest homolog (contig)`, is
required alongside the score even though neither is scored.

Registered before any C3 or F9 score existed.

### A4 · 2026-08-01 — C5/T7 is scored as a ratio, not as an absolute

**Affects C5 and T7.**

The published burden values and the tool's export are **in different
units**, and comparing them directly would be an error of kind, not of
precision.

-   The source normalises microbial readcounts to host cells and reports
    **microbial transcripts (reads) per host cell**. Host cells are the
    mean PGK1/NSE control readcount ÷ \~50 transcripts per cell — the
    same rule COBLAST+ implements (`HOST_TRANSCRIPTS_PER_CELL`,
    `etol_summary.compute_host_cells`). This is what `Reads per host
    cell` contains.
-   The Discussion's 0.14 bacteria / 0.05 fungi / 1.8 microbes per host
    cell are **organisms** per host cell. Neither the 2023 paper nor the
    2022 method paper states a reads → organisms conversion; 2022 gives
    only an order-of-magnitude assumption (\~2,000 ribosomes per
    slow-growing bacterial cell) used for sensitivity, never applied to
    the reported abundances.

That the two units diverge is visible **inside the source itself**: by
reads per host cell the paper ranks **fungi above bacteria** (Fig. 1,
Fig. 4A — the basis of C1), while by organisms per host cell it reports
bacteria (0.14) above fungi (0.05). Both statements can hold only if
fungal per-organism read yield exceeds bacterial by more than \~2.8×,
which is biologically unsurprising and analytically fatal to a direct
comparison. It also means the bacteria:fungi ratio is **not** unit-free,
because the conversion factor differs between the two domains.

**Amendment: C5 scores the max-to-control fold change (≈ 9.5×), not the
absolutes.** That ratio holds the same domain mixture in numerator and
denominator, so it survives the units gap far better than any per-domain
figure. T7 still prints the bacteria, fungi and max-case absolutes in
reads per host cell — the transcription the table exists for — but the
comparison column marks them **not scored against 0.14 / 0.05 / 1.8**,
with the unit mismatch stated in one sentence in the caption.

Two consequences to carry into the write-up:

1.  **Do not write "bacteria exceed fungi" in T7.** In the units
    COBLAST+ exports, the source says the opposite, and C1 already
    pre-registers fungi first.
2.  **M66 is an EBB donor**, so the max-burden case is a same-data target
    rather than a borrowed one. Once T2 carries EBB age and sex, check
    whether COBLAST+'s highest-burden donor is the male aged 66; a hit is
    a strong, cheap validation point and a miss is worth one honest
    sentence.

Registered before any burden number was computed.

### A5 · 2026-08-04 — V1 reclassified as verification; V2 scored as a Miss; SARS-CoV-2 resolved

**Affects V1, V2 and unresolved cell 1. Recorded AFTER the run it
concerns — post-hoc by construction, unlike A1–A4.**

The 35-sample eToL-V run (batch `fdda627a`, per-cell export
`etol_v_cells.csv`) was scored before this text existed. Nothing below
moves a band or reverses an outcome. What changes is where V1 sits, and
what the eToL-V leg is permitted to claim.

**1. V1 is a verification artefact, not a validation criterion.**

`test_reproduces_veso_confusion_matrix` does not recompute the reference
matrix from the reference's data. It *constructs* a prediction matrix
stipulating the reference's surviving calls — nine `V-HAdV_AdC_penton`
cells plus the lone `V-HPV_HPV45_L1` hit in SRX17674444 — and then
asserts that the scoring arithmetic, the inverted SRR↔SRX join and the
universe construction yield TP 9 / FP 1 / FN 35 / TN 411. Ground truth
and universe are loaded from the bundled CSV; the predictions are
asserted. The registered "Scored from" cell omits `compute_confusion`'s
first argument, `matrix`, which is what made it read as an end-to-end
re-derivation.

That is a known-answer fixture — the same category as the spike-in
positive control and the frozen self-check, which this document's Scope
routes to §2.4.

**Amendment: V1 moves to §2.4 as a verification artefact**, described as
a scoring-logic and ground-truth-join known-answer check, cited as
`tests/test_etol_validation.py::test_reproduces_veso_confusion_matrix`.
It passes. Its pass is evidence that COBLAST+'s scoring layer is
faithful — not evidence about adenovirus. Its fixture encodes the
reference's ground truth, which Lathe's corrected ground truth
supersedes; that is exactly why it verifies machinery and not biology.

**2. V2 is scored as a Miss.** Trigger: precision .40 \< .90.

| Arm | TP | FP | FN | TN | Precision | Recall |
|---------|---------|---------|---------|---------|---------|---------|
| `stage="validated"` (scored) | 2 | 3 | 43 | 792 | .40 | .044 |
| `stage="raw"` (sensitivity) | 9 | 7 | 36 | 788 | .56 | .200 |

N = 840, accuracy .945. Of the registered Pass conditions only "every TP
an HAdV-C cell" held, on both arms. The Partial band's attribution
requirement was satisfied and is reported regardless: of 43 false
negatives, **36 carried no net read at all** (`raw = 0`) and **7 carried
net reads that were not contig-confirmed**, all seven HAdV-C. Every one
of the 14 non-adenovirus false negatives is in the `raw = 0` bucket, so
net-stage failure accounts for the entire non-adenovirus shortfall and
assembly depth accounts only for adenovirus C. Six of the seven carried
a single net read, below `Cap3Assembler.MIN_READS = 2`, so CAP3 was never
invoked and no contig could form under any parameterisation; both true
positives carried three. Contig *identification* is therefore not
implicated, and no paragraph should blame `ToL_virus_val`.

**3. SARS-CoV-2 predictions are scored, not excluded.** This closes
unresolved cell 1 against the assumption recorded there.

A truth column of zeros is a negative result, not absent data: the
corrected ground truth asserts these samples contain no SARS-CoV-2. All
seven raw-stage SARS-CoV-2 hits fall in SD030-18, SD042-18, SD014-17,
SD032-17 and SD025-19; if the two-digit suffix encodes collection year —
**to confirm with Lathe before this is asserted in text** — every one
predates the virus's emergence, making them definitionally false
positives and pointing at cross-reactivity of the single
`V-HCoV_SARSCoV2_S` spike probe. Three survived contig confirmation.
Excluding them would delete the most informative cells in the run, and
would remove precisely the cells that caused the precision failure.

**4. The corrected-ground-truth evaluation is Future Work, not a new
criterion.** A characterisation of eToL-V against Lathe's ground truth
has no failure condition, and every number in it is already known, so
registering it as a criterion would be theatre. The stage-resolved
matrix and the attribution above are reported descriptively in the
Results and Discussion; the case for a properly pre-registered
evaluation belongs in Future Work.

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
    **excluded**. *(Closed by amendment A5: they are **scored**, against
    the assumption recorded here.)*
2.  **C3 — inherited shortlist, not re-derived.** The paper's top-10 was
    selected from a *pooled* four-dataset ranking; you cannot re-derive
    that selection from EBB alone. The criterion therefore inherits the
    published shortlist and asks only whether those genera reproduce in
    EBB. If you want per-genus EBB percentages as reference values, they
    must be read off Figure 8C's EBB (yellow) bars — the text gives only
    the 47 %→33 % pooled range. Absent that, keep the criterion at
    set-overlap and direction, which the text does support. *(The
    sub-question of what the tenth, unnamed entry is and how it matches
    is now closed by amendment A3.)*

## Excluded from validation, with reason

| Published result | Why it cannot be a criterion |
|------------------------------------|------------------------------------|
| Pooled differential abundance: top-100 mean overabundance **13.68** (range 2.23–181.97, SD 20.40, P = 0.01); cluster mean **12.59** (range 2.52–36.88, SD 7.83, P \< 0.001) | Computed over MSBB + EBB + Miami + Rockefeller (31 control, 48 AD). MSBB is access-gated via the AD Knowledge Portal; EBB alone cannot reproduce the pooled statistic |
| The 23S/28S (LSU) confirmation arm behind Figs. 8–9 | COBLAST+ has no LSU disambiguation/reprobe pass — the documented gap in §2.5.3. `ToL_rRNA` contains LSU, but nothing uses it this way. Future Work. *(This is also why the tenth AD-shortlist entry cannot be matched by name — see A3)* |
| Viral findings: AdC = **83 %** of viral transcripts, detected in **\~50 %** of samples; viral burden higher in AD **P = 0.04**; HHV↔AdC inverse relationship in AD only; cutoff **0.03** reads per host cell | Different instrument: the paper probes with *stripped whole genomes* of the top-20 Readhead viruses, not eToL-V's 120-mer structural-protein panel. A disagreement here is a method difference, not a port defect. Cite as context for V2's recall, never as a target |
| Endogenous retroelement / HERV abundances | The extended retroelement probe list is not in COBLAST+ |
| Brain microbiome ≈ **20 %** of gut microbiome diversity; gut \~4× more diverse | Needs gut/faecal RNA-seq; no host normalisation possible |
| Cross-species conservation from *Drosophila* to human (Fig. 2) | Needs non-human brain datasets and has no host-cell normalisation |
| Microbe/retroelement/virus correlations R = 0.394, 0.267, 0.139 | Depends on the retroelement and stripped-genome arms above. (The source labels these both "coefficients of determination" and "R2 values"; do not quote a value without resolving which it is) |
| Absolute organism counts per host cell (0.14 bacteria, 0.05 fungi, 1.8 max) as *value-level* targets | Stated in organism units the tool does not export and neither paper converts to. Retained as **C5** in ratio form only — see amendment A4 |
| Hu et al. 2022 cohort results (MSBB, Rockefeller, Miami) | Different cohorts, access-gated; the 2023 paper supersedes it as the same-data target |

## Not in this table, on purpose

Anything the tool cannot currently export. A criterion that needs code
written to score it is not pre-registered, it is a plan.
