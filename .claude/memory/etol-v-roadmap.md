---
name: etol-v-roadmap
description: "How to fold the previous student's eToL-V viral-detection workflow into COBLAST+ (it's a preset, not a new pipeline)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 6a172d6a-d6d5-4e4c-8cb2-8db35f467042
---

eToL-V = "electronic Tree of Life for Viruses", the MSc dissertation (Edinburgh,
B270917, 2024-25, supervised by Lathe) the user inherited. Goal: screen human
brain RNA-seq for viruses using ~115 120-nt probes sliced from conserved viral
structural-protein CDS (herpes MCP/gB/gH, adeno hexon/penton/fibre, HPV L1/L2,
corona Spike) + 2 PGK housekeeping = 117 probes (`final_probes.fasta`).

KEY CONCLUSION (2026-06-25): eToL-V is mechanically the SAME shape as the eToL
preset already built in COBLAST+ — probe FASTA -> BLASTn net vs SRA reads ->
human-artifact removal -> PGK /50 normalization -> CAP3 contigs -> validate
contigs vs reference DB. So it is a NEW PRESET, not a new pipeline. The
dissertation's normalization (mean PGK / 50) is identical to
`HOST_TRANSCRIPTS_PER_CELL=50.0` in etol_summary.py. See [[etol-lathe-alignment]].

STATUS (2026-06-29): PRESET BUILT. All 3 pieces below are DONE. `etol_v`
registered in `ETOL_PRESETS` (etol_summary.py); 117 probes (115 viral + this
panel's OWN 2 PGK controls) loaded from `data/etol_v_probes.fasta` +
`data/etol_v_control_probes.fasta`. Headers renamed to the cellular grammar with
viral class codes V-HHV(77)/V-HAdV(7)/V-HPV(10)/V-HCoV(21); `ETOL_DOMAIN_BY_LETTER`
gained `"V"→Viruses`; controls were made per-preset (`controls` callable in each
preset dict + `_control_pairs_for(key)`; cellular presets unchanged, still 4
PGK1/hNSE controls). Subunit kept in the taxon so AdC penton is its own row.
101 tests pass (3 new etol_v tests in tests/test_summaries.py). Validation DB
build script shipped: `scripts/build_etol_v_validation_db.sh` (RefSeq viral +
GRCh38 incl. mito → `ToL_virus_val`; mixed DB so genuine viral contigs survive and
herpes→mito artifacts get a "Homo sapiens" homolog and are dropped by
contig_id.py). REMAINING before a real run: (a) run that script on the server to
build+register `ToL_virus_val`; (b) build patient BLAST DBs from the new RNAseq;
(c) shakedown run + verify net E<0.01 keeps Veso's 9 TPs. Original 3 pieces,
now done:
1. Register an `etol_v` preset in `ETOL_PRESETS` (etol_summary.py:155). For
   fidelity, normalize on the preset's OWN 2 PGK probes (shipped in
   `final_probes.fasta`), NOT the cellular preset's PGK1+hNSE. Norm = mean(PGK)/50
   (identical to `HOST_TRANSCRIPTS_PER_CELL=50`; Veso then ×10 for per-10-cell
   display).
2. Viral-aware header parsing: rename probes into the existing
   `Class_Taxon_Subunit_Index` grammar (e.g. `V-HHV_HSV1_gB_3`,
   `V-HCoV_SARSCoV2_S_1`) and extend `ETOL_DOMAIN_BY_LETTER` with viral class
   codes, so `_class_code/_taxon/_species/_domain` are reused, not forked.
3. Contig artifact rejection (FIDELITY-CRITICAL): `identify_contigs`
   (contig_id.py) already keeps the best contig homolog + flags Homo
   sapiens/mitochondrion. BUT the validation DB must be VIRUS-APPROPRIATE —
   `core_nt`, OR a local RefSeq-viral + human-genome + human-mito build. Do NOT
   reuse the cellular `ToL_rRNA` (SILVA) DB: viruses have no rRNA, so a herpes→
   human-mitochondria artifact contig won't match SILVA at all and would slip
   through as "no hit" instead of being caught. Herpes→mito was the DOMINANT
   false-positive source in Veso's results — this step is what kills it.

Prereq: `final_probes.fasta` is NOT in the repo — get it from the student's
GitHub (github.com/B270917-2024/MSc_Dissertation). Build order: get FASTA ->
rename headers -> register preset (net/human-filter/dedup/PGK/CAP3/CSV all reused)
-> add nt validation. Part I (t-SNE, NJ trees, probe design scripts) is
out of scope for a screening tool. Paper-reported perf: 90% precision, 20% recall.

VISUALIZATION (started 2026-06-29; supersedes the earlier "heatmap = YAGNI"
note): building viz for BOTH presets as 3 tiers. Both papers' core artifact is
the SAME rows×samples matrix COBLAST+ already exports (etol_probe_count_rows /
etol_species_count_rows); presets differ only in params (cellular = log2 reads/
host cell, cutoff 3-5, Hu/Lathe pheatmap/Morpheus; viral = RAW hit count + a
"Cell count" annotation row, two-stage raw→validated = Veso Fig 8→10 matplotlib;
plus Veso's confusion matrix Fig 9 vs eToL WGS ground truth, probe×sample,
SARS-CoV-2 excluded). TIER 1 DONE: in-app vanilla-JS SVG heatmap (no lib, exe
stays lean) — new `build_etol_matrix` + `_sample_condition` (etol_summary.py),
`etol_matrix_payload` (result_store.py), GET `/batch-results/<id>/etol-matrix.json`
(app.py), `static/etol_heatmap.js` + styles + panel in batch_results.html;
per-preset defaults, value/stage/cutoff controls, condition column swatches,
PNG/SVG export; 103 tests pass (+2 new). TIER 2 (confusion matrix): user supplied
Veso's actual WGS ground truth = `WHOLE GENOME SCAN OF EBB.xls` (Sheet1, 20
viruses × 35 samples, COUNT data). CONSTRUCTION FULLY REVERSE-ENGINEERED + VERIFIED
(2026-06-29) — reproduces her TP=9/FP=1/FN=35/TN=411, N=456 exactly:
 • BINARY, present = count>0 (paper p40 "binary classification model"; >0 is the
   ONLY threshold giving TP+FN=44 — confirmed by parsing the xls).
 • Sample key: WGS columns are shorthand "srxNN" = SRX176744NN (verified: yields
   exactly 35 cols, 24 AD / 11 CONTROL). USER RUNS ON SRR, and COBLAST+'s
   accession regex matches SRX/ERX/DRX NOT SRR — so a SRR↔SRX176744NN crosswalk is
   the one external blocker (SRA Run Selector / Veso Table S2/S3). Consider adding
   SRR|ERR|DRR to ETOL_ACCESSION_PATTERN (also fixes heatmap labels for SRR runs).
 • Universe = the 13 eToL-V-targetable viruses that HAVE a WGS row (AdenoC, COV_229E,
   HHV1_HSV1, HHV2_HSV2, HHV3_VZV, HHV4_EBV, HHV5_CMV, HHV6A, HHV6B, HHV7, HHV8,
   HPV6, HPV16) × 35 samples = 455, +1 for the single HPV45-L1 prediction lacking a
   WGS row (= her FP=1) = 456. WGS-positive cells: AdenoC=30, HSV1=7, CMV=5, EBV=1,
   HHV8=1 → 44 (=TP+FN). ✓
 • eToL-V "present" = any of that virus's probes has a VALIDATED (contig-confirmed)
   hit >0 (compare to Fig 10, not raw Fig 8). SARS-CoV-2 EXCLUDED (no WGS row);
   other panel viruses with no WGS row (HPV11/18/31/45, AdB, NL63/OC43/HKU1/MERS)
   score as WGS-absent → a validated prediction there is an FP.
 • Build: virus→COBLAST+ taxon crosswalk + loader consuming `sample,virus,count`
   (binarize >0) + TP/FP/FN/TN + Acc/Prec/Rec/F1 panel. Working xls parser already
   prototyped (Excel COM → CSV; no pandas/xlrd in env).
 CROSSWALK RESOLVED + VERIFIED (2026-06-29): user supplied `DATA UPLOADED TO NCBI
   INCLUDING SRAs (1).ods`; sheet ALL_NCBI_DATA is the full SRA run table for study
   SRP398685 with Run+Experiment+SampleName in one row. SRR↔SRX mapping is INVERTED,
   not same-suffix: SRX17674433↔SRR21676133, …34↔…132, …35↔…131 (verify, don't
   assume!). Closes the join: COBLAST+ label (SRR) → crosswalk → SRX176744NN → WGS
   truth. Verified 35/35 region agreement (diagnosis 27/35 — the 8 "misses" are
   AD/VaD+AD/LBD samples WGS coarsely called "AD"; same AD/control split). Crosswalk
   CSV + normalized truth parsed in scratchpad.
 UNIVERSE PINNED: Veso's 13-virus set = {Adenovirus C, COV_229E, HHV1_HSV1,
   HHV2_HSV2, HHV3_VZV, HHV4_EBV, HHV5_CMV, HHV6A, HHV6B, HHV7, HHV8, HPV6, HPV16}
   × 35 = 455, +1 HPV45 out-of-universe FP = 456. Reproduces 9/1/35/411 EXACTLY
   (44 positives: AdC 30, HSV1 7, CMV 5, EBV 1, HHV8 1). QUIRK (judgment call, worth
   Veso confirm): includes HPV6 though the panel has NO HPV6 probe (→always predicted
   negative), and EXCLUDES Adenovirus A / Adenovirus 54 though they're adenoviruses.
   SARS-CoV/SARS-CoV-2 fully excluded (no WGS row). Panel viruses w/ no WGS row +
   a validated hit = FP (only HPV45 fired). Map: AdC→V-HAdV_AdC_*, 229E→HCoV229E,
   HSV1/2,VZV,EBV,CMV→HHV1/2/3/4/5, HHV6A/6B/7, HHV8→KSHV, HPV16; HPV6→no taxon.
 TIER 2 CORE DONE (2026-06-29): shipped `data/etol_v_wgs_truth.csv` (srx,virus,count;
   700 rows) + `data/etol_v_sra_crosswalk.csv` (srr↔srx↔region↔diagnosis; 35) +
   `etol_validation.py` (loaders, `VESO_UNIVERSE` configurable, `compute_confusion`
   over a build_etol_matrix payload; validated stage, SRR→SRX join, out-of-universe
   FP, SARS excluded). `tests/test_etol_validation.py` REPRODUCES her Fig 9 exactly
   (9/1/35/411, acc .9211/prec .90/rec .2045/F1 .33) from the bundled truth — the
   fidelity proof. ALL VISUALIZATION TIERS NOW DONE:
 TIER 2 PANEL DONE (2026-06-29): confusion matrix computed in the /batch-blast POST
   route (eToL-V only, guarded), persisted in the batch payload, server-rendered as a
   2×2 + Acc/Prec/Rec/F1 panel in batch_results.html (.confusion CSS); no-overlap and
   error branches handled. Template render smoke-tested.
 TIER 3 DONE (2026-06-29): `scripts/plot_etol.py` (matplotlib-only, NOT in the exe,
   reuses compute_confusion) renders her Fig 9 confusion matrix (+ optional Fig 10
   heatmap). `--batch-id` or `--matrix-json`; `--print-only` needs no plotting libs and
   verified TP=9/FP=1/FN=35/TN=411. 110 tests pass.
 So eToL/eToL-V data-vis is COMPLETE: in-app heatmap (Tier 1) + in-app confusion-matrix
   panel (Tier 2) + publication script (Tier 3). Only runtime step left for the user: run
   the eToL-V preset on the 35 SRP398685/EBB SRR samples (with contig id on) to populate
   real figures; the SRR labels join the SRX-keyed truth via data/etol_v_sra_crosswalk.csv.
 CONFUSION CSV EXPORT DONE (2026-06-29): compute_confusion now emits per-cell `cells`
   (result, virus, sample, srx, wgs_count, actual, predicted, RAW_HITS, CONFIRMED_HITS);
   `etol_confusion_rows_as_delimited` (result_store) + GET /batch-results/<id>/etol-
   confusion.{csv,tsv} (app.py) + buttons in batch_results.html. 111 tests pass.
 FIRST REAL RUN (2026-06-29, user, 35 EBB samples, contig assembly on): TP2/FP0/FN42/
   TN411 (N=455) vs Veso TP9/FP1/FN35/TN411 (N=456). DIAGNOSIS: all 35 joined (N=455=
   13×35), TN identical. FP0 = HPV45 error not reproduced (user aware, the −1 cell). The
   7 extra FN are ALL adenovirus-C-penton (Veso validated AdC in 9 samples, user in 2).
   Lost at one of two stages, told apart by the new CSV's RAW vs VALIDATED columns:
   raw 0 = NET (COBLAST gates E<0.01 in filter_net_probe_hits; Veso web BLASTn had NO
   gate = default E≤10, so kept weak AdC reads COBLAST drops; max_target_seqs=5e6 so NOT
   truncation); raw>0 & validated 0 = contig stage (CAP3 MIN_READS=2 singleton floor +
   DEFAULT_CONFIRM_IDENTITY_PCT=99). Likely fix if raw=0 dominates: relax the eToL-V net
   E-gate toward Veso's permissiveness. NB even Veso only validated 9/30 AdC (recall ~20%)
   — contig validation is inherently lossy at low abundance. TIER 3:
optional `scripts/plot_etol.py` (matplotlib/seaborn/sklearn over the CSVs), NOT
bundled in the exe — paper-pixel-faithful clustermap + ConfusionMatrixDisplay.
User chose "both" (in-app + script). NOTE: user syncs in-repo git-tracked
.claude/memory over this global store — the Tier-1 code now lives in the repo.

CONFIRMED FROM VESO'S FULL DISSERTATION (2026-06-26):
- Validated brain virome is THIN: after contig validation only adenovirus-C penton
  + SARS-CoV-2 probes survive; ~all herpesvirus hits were artifacts matching the
  human MITOCHONDRIAL genome. Recall ~20% is probe-limited, NOT port-limited.
- FIDELITY TARGET: a faithful port should reproduce Veso's confusion matrix vs the
  eToL WGS ground truth on the 35 EBB samples number-for-number: TP=9, FP=1,
  FN=35, TN=411 (acc 92%, prec 90%, recall 20%, F1 0.3). This IS the proof the
  port matches her methods.
- Strand conversion (plusstrand.py) was probe PREP (Part I), not a runtime step —
  final_probes.fasta is already plus-strand; no per-run work.
- Net: Veso used web BLASTn, max_target_seqs 5000, counted all hits then validated
  via contigs (no explicit E gate). COBLAST nets at E<0.01 — keep it, but verify
  it doesn't drop the 9 TPs (AdC/SARS hits are strong, so safe).
- Probe set: advisors confirm it's been EXTENDED since Veso to adequate breadth;
  probe-generation automation NOT needed. Still pull the latest final_probes.fasta
  from the student GitHub.
- SARS-CoV-2 anomaly: pre-pandemic EBB samples (deaths 2017-2019) carry sequences
  that phylogenetically cluster with the pandemic SARS-CoV-2 clade — unresolved
  (contamination vs real). Optional high-interest, rabbit-hole-risk thread.
- COBLAST+ removes eToL-V's stated #1 limitation (can't BLASTn vs SRA locally) and
  beats ViromeScan (Veso: 50h + TB/batch, failed on 34/35 samples).
- See [[coblast-ad-application-plan]] for the staged AD-application plan this feeds:
  validation floor (reproduce the matrix above + Hu 2023 cellular shortlist) then
  the AdC-anchored cellular+viral integrated stretch.
