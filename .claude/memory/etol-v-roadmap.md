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

Only 3 genuinely-new pieces:
1. Register an `etol_v` preset in `ETOL_PRESETS` (etol_summary.py:155); reuse
   existing PGK1/hNSE control probes for normalization.
2. Viral-aware header parsing: rename probes into the existing
   `Class_Taxon_Subunit_Index` grammar (e.g. `V-HHV_HSV1_gB_3`,
   `V-HCoV_SARSCoV2_S_1`) and extend `ETOL_DOMAIN_BY_LETTER` with viral class
   codes, so `_class_code/_taxon/_species/_domain` are reused, not forked.
3. Contig artifact rejection: `identify_contigs` (contig_id.py) already keeps the
   best contig homolog; point it at nt/core_nt and flag contigs whose best hit is
   Homo sapiens / mitochondrion (the paper's herpes-mito false-positive problem).

Prereq: `final_probes.fasta` is NOT in the repo — get it from the student's
GitHub (github.com/B270917-2024/MSc_Dissertation). Build order: get FASTA ->
rename headers -> register preset (net/human-filter/dedup/PGK/CAP3/CSV all reused)
-> add nt validation. Heatmap + WGS precision/recall metrics are presentation/
analysis, NOT pipeline (YAGNI). Part I (t-SNE, NJ trees, probe design scripts) is
out of scope for a screening tool. Paper-reported perf: 90% precision, 20% recall.
