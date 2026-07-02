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

**ADVISOR'S RESEARCH-QUESTION MENU (received 2026-06-29, "each could be a paper"):**
(1) Does Down syndrome brain contain elevated microbes? (2) Does microbial burden
increase with age? (3) 5xFAD transgenic mice die early — elevated microbes in their
brain (messed-up antimicrobial defence)? needs MOUSE host-subtraction. (4) How
frequent is APOE editing in higher primates? 1 probe × ~1000 SRAs (advisor has the
probes + a sequence-error-rate check). All answerable from public NCBI SRAs.

COMPARISON verdict: advisor's = BREADTH plays (apply same method to a NEW
cohort/axis/organism, each with a built-in comparator/gradient → de-risked,
paper-shaped, dodges the "not all AD positive" heterogeneity). My earlier A–D =
DEPTH plays within the human AD set (co-occurrence/spatial/subset/anomaly) — higher
novelty but data-hungry + fragile on small n. For the timeline, advisor's instinct
wins the SPINE. RECOMMENDATION: spine = **burden vs age (#2)** [or **Down syndrome
(#1)** as case-control], run through BOTH cellular + viral layers on one patient DB
(does cellular burden AND AdC rise with age?) — this folds in the AdC×high-burden
question AND exercises eToL-V (advisor's list has NO viral component otherwise).
Bonuses if time: APOE-in-primates (#4) is nearly free — the exact-match APOE preset
still exists in COBLAST+, advisor has probes. SHELVE: 5xFAD mouse (#3) needs a mouse
genome BLAST DB — but note `human_filter.py` is DB-agnostic so it's "build+register a
mouse DB" (adapt build_etol_v_validation_db.sh), main cost is model caveats; and the
SARS-CoV-2 anomaly (rabbit hole). Timeline anchor: as of 2026-06-29, ~7 weeks to full
dissertation (≈ mid-Aug 2026), writing NOT yet started.

**ADVISOR VERDICT (2026-06-29, Doug + Rich replied, both enthusiastic):**
- 2 more dev days APPROVED.
- Re-probing: make it OFF by default, keep as an easy switch (needed only for future
  less-specific probe sets). Doug wants a MANUAL check that 0-novel-reads isn't a
  coding bug: BLAST a key contig back at the patient DB (megablast, E<0.01,
  max_target_seqs 100000, outfmt "6 sseqid"), `comm -23` vs that taxon's net reads;
  empty = correct (contig is built FROM those reads), non-empty while COBLAST says 0
  = bug. Reuse scripts/run_manual_blast.sh.
- Human false-positive strategy (Rich's "any ideas?"): COBLAST+ ALREADY automates
  his manual confirm (read filter → contig-level BLAST-back-vs-human → drop "Homo
  sapiens" contigs), and does it at contig level over the whole batch (stronger than
  his 100-read spot check). Proposed additions: (a) make the human component of the
  validation DB comprehensive (genome + mito + human rRNA/transcriptome — can only
  reject what's in the DB); (b) emit a per-call provenance/audit line + COUNT the
  human contigs dropped per taxon instead of dropping silently (the reviewer/clinician
  audit trail); high-confidence tier = passes contig (+ LSU for headline) validation.
- LSU (26S/28S): Rich says no easy automation but easy to pull a probe bunch. POSITION:
  targeted confirmation of HEADLINE taxa only, not blanket; cheap to add as a 2nd
  preset IF he supplies probes (eToL-V proved the machinery reuses); else skip. YAGNI.
- **NOVEL STRETCH RESHAPED** (Rich: "point the telescope at several stars"): NOT one
  question — combine **microbes-vs-age + Down syndrome + 5xFAD mouse** into ONE
  "resolution-of-the-telescope" multi-cohort paper. 5xFAD is back ON (human_filter.py
  is DB-agnostic → mouse just needs a mouse-genome BLAST DB via the validation-DB
  script; no pipeline change). Pragmatic confirmation: spot-check high-abundance reads
  only, do NOT back-check every read — this is what makes 3 cohorts fit the timeline.
  Rich has age SRAs to resend. AdC×cellular-burden rides along as the cross-layer angle.
- Publish: BioRxiv preprint (after validation floor) + GitHub (already live). Rich
  wants tool access ("my hands on your telescope") — offer a build/walkthrough after the
  2-day push. NCBI intro: ask Rich. Edinburgh Innovation (commercialisation): after data.

**Headline narrative:** COBLAST+ IS the tool Veso explicitly called for — it removes eToL-V's stated #1 limitation (can't BLASTn vs SRA locally) by running locally vs patient-built DBs, and beats ViromeScan (her run: 50h + TB/batch, failed on 34/35 samples).

**eToL-V port = preset + viral headers + ONE new validation DB** (net/human-filter/dedup/PGK-norm/CAP3/CSV all reused):
1. Register `etol_v` preset; use the preset's own 2 PGK controls (not PGK1+hNSE). Norm = mean(PGK)/50 (identical to HOST_TRANSCRIPTS_PER_CELL=50; Veso ×10 for per-10-cell presentation).
2. Viral header parsing into Class_Taxon_Subunit_Index (e.g. `V-HHV_HSV1_gB_3`, `V-HAdV_AdC_penton_1`, `V-HCoV_SARSCoV2_S_1`); extend ETOL_DOMAIN_BY_LETTER with viral class codes.
3. **Contig validation DB must be virus-appropriate** (core_nt OR local RefSeq-viral + human-genome + human-mito) — NOT the rRNA ToL_rRNA DB (viruses have no rRNA; the mito artefact won't match SILVA and would slip through). Reuse contig_id.py Homo-sapiens/mito flagging — this is the fidelity-critical step.
Probe set (CORRECTED 2026-07-01): repo still has Veso's ORIGINAL 115 viral probes; extending the panel to what Dr. Lathe would require is a higher bar than "adequate breadth" and is now scoped as incomplete FUTURE WORK (eToL-V beyond the existing implementation is a contingent deliverable, not a results chapter — see [[etol-v-roadmap]]). Fidelity decision: keep net E<0.01 gate but verify it doesn't drop Veso's 9 TPs.

**Why:** problem-choice framework — fix the floor (zero miracles; satisfies proposal's correctness criterion) before novel biology on a heterogeneous phenotype.

**How to apply:** Remaining cellular pipeline gaps are now small (contig species-ID/re-probing already DONE per [[etol-lathe-alignment]]): 23S/28S disambiguation hook + heatmap; read-length-adaptive cutoff optional. Offered to draft a formal research-strategy doc next. NOTE: user curates in-repo `.claude/memory` (git-tracked) and syncs it over this global store — ask before relying on this file persisting; consider updating repo `etol-v-roadmap.md` instead.
