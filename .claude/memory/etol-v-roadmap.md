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
PNG/SVG export; 103 tests pass (+2 new). TIER 2 (confusion matrix, user said
"wire in now", HAS the WGS ground truth): pending — consume CSV `sample,taxon,
present` (omit rows = excluded, e.g. SARS-CoV-2), compute TP/FP/FN/TN at
probe×sample + Acc/Prec/Rec/F1; fidelity target TP=9/FP=1/FN=35/TN=411. TIER 3:
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
