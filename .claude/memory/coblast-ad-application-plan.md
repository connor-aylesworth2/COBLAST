---
name: coblast-ad-application-plan
description: Strategy + eToL-V porting spec for applying COBLAST+ full eToL(+V) to the brain-microbiome × AD dissertation
metadata: 
  node_type: memory
  type: project
  originSessionId: f4dae7d6-9eba-406e-aa16-863dc8db33bf
---

Strategy for the dissertation's AD application, refined 2026-06-26 after reading Veso's full dissertation (B270917, "eToL-V: development of a rapid method for detecting viruses", Edinburgh 2024-25) and re-syncing COBLAST+ state. Builds on [[etol-lathe-alignment]] and [[etol-v-roadmap]].

**Shape = staged:** validation floor (guaranteed) → ONE novel stretch. Data = public AD brain RNA-seq backbone + 35 EBB limbic samples if access holds.

**Validation floor has TWO precise targets:** (a) cellular = reproduce Hu 2023 AD-overabundant shortlist (Cortinarius, Aspergillus-group, Cryptococcus/Tausonia, Komagataella/Candida; Sphingomonas, Streptococcus, Staphylococcus) + AdC dominance; (b) viral = reproduce Veso's confusion matrix number-for-number (TP=9, FP=1, FN=35, TN=411; only AdC penton + SARS-CoV-2 survive validation; herpes→human-mitochondria artefacts eliminated). Target (b) IS the fidelity proof that the eToL-V port matches Veso.

**Novel stretch (reframed):** the brain virome is THIN (validated = AdC penton + SARS-CoV-2 anomaly; herpes artefactual; recall ~20%, probe-limited not port-limited). So NOT rich tri-kingdom co-occurrence. Instead: "Does adenovirus-C presence mark the high-cellular-burden (bacterial+fungal) AD subset?" — first time both layers run on the SAME local patient DB in one tool; engages Hu 2023's AdC-vs-HHV inverse relationship + "~half of AD microbe-positive / subtype" thread. Optional time-boxed branch: deepen the SARS-CoV-2 pre-pandemic (2017-2019 EBB deaths) anomaly that phylogenetically matches the pandemic clade — high interest, rabbit-hole risk.

**Headline narrative:** COBLAST+ IS the tool Veso explicitly called for — it removes eToL-V's stated #1 limitation (can't BLASTn vs SRA locally) by running locally vs patient-built DBs, and beats ViromeScan (her run: 50h + TB/batch, failed on 34/35 samples).

**eToL-V port = preset + viral headers + ONE new validation DB** (net/human-filter/dedup/PGK-norm/CAP3/CSV all reused):
1. Register `etol_v` preset; use the preset's own 2 PGK controls (not PGK1+hNSE). Norm = mean(PGK)/50 (identical to HOST_TRANSCRIPTS_PER_CELL=50; Veso ×10 for per-10-cell presentation).
2. Viral header parsing into Class_Taxon_Subunit_Index (e.g. `V-HHV_HSV1_gB_3`, `V-HAdV_AdC_penton_1`, `V-HCoV_SARSCoV2_S_1`); extend ETOL_DOMAIN_BY_LETTER with viral class codes.
3. **Contig validation DB must be virus-appropriate** (core_nt OR local RefSeq-viral + human-genome + human-mito) — NOT the rRNA ToL_rRNA DB (viruses have no rRNA; the mito artefact won't match SILVA and would slip through). Reuse contig_id.py Homo-sapiens/mito flagging — this is the fidelity-critical step.
Prereq: pull EXTENDED `final_probes.fasta` from github.com/B270917-2024/MSc_Dissertation (advisors confirm breadth now adequate; probe automation NOT needed). Fidelity decision: keep net E<0.01 gate but verify it doesn't drop Veso's 9 TPs.

**Why:** problem-choice framework — fix the floor (zero miracles; satisfies proposal's correctness criterion) before novel biology on a heterogeneous phenotype.

**How to apply:** Remaining cellular pipeline gaps are now small (contig species-ID/re-probing already DONE per [[etol-lathe-alignment]]): 23S/28S disambiguation hook + heatmap; read-length-adaptive cutoff optional. Offered to draft a formal research-strategy doc next. NOTE: user curates in-repo `.claude/memory` (git-tracked) and syncs it over this global store — ask before relying on this file persisting; consider updating repo `etol-v-roadmap.md` instead.
