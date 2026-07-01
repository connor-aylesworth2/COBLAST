---
name: coblast-verification-anchor
description: "The spike-in positive control is COBLAST+'s dissertation \"verification\" artifact, distinct from result reproduction (\"validation\")"
metadata: 
  node_type: memory
  type: project
  originSessionId: f83d486d-845c-43f8-bd95-d7bdac6ad3f2
---

For the dissertation, COBLAST+'s credibility rests on two separate legs, and conflating them is the trap:

- **Verification** ("did I build the tool right?") — the existing unit/smoke tests prove the *code* is correct (arithmetic, parsers, and `tests/test_etol_validation.py` reproduces Veso's confusion matrix but by *hand-injecting* the surviving calls, so it only proves the scoring harness). The new `tests/test_spike_in_control.py` is the stronger verification anchor: a synthetic spike-in positive control that drives the *real* eToL net path (`run_blast_probe_panel` → `filter_net_probe_hits` → `count_control_reads` → `deduplicate_reads_to_best_probe` → `build_etol_probe_summary`) over a known-composition sample and asserts recovery, specificity, and /50 normalization. Run: `python tests/test_spike_in_control.py` (prints an expected-vs-observed table; needs BLAST+).

- **Validation** ("did I build a biologically right tool?") — reproducing Lathe 2022 / Veso 2023 on real data. Can't be proven off-data; a bug-free pipeline can still be biologically wrong. See [[coblast-ad-application-plan]] and [[etol-lathe-alignment]].

**Why:** a reader/examiner wants proof COBLAST+ is logically sound *before* it's pointed at real data, so that tool-vs-paper disagreements read as biology, not bugs. **How to apply:** cite the spike-in as the verification chapter (positive control), and structure the reproduction as a graded ladder whose first rung is a known-answer control (the eToL-V/Veso confusion matrix run end-to-end, currently TP2/FP0/FN42 vs 9/1/35 — the E-gate fix is the pending move). Don't answer the "is it sound?" worry by writing more unit tests — that's the trap.
