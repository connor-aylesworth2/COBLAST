---
name: coblast-frozen-stages-silent-noop
description: "Frozen COBLAST+ exe uses a different data dir than a source run, so resource-gated post-BLAST stages (human filter, CAP3 assembly, contig ID) can silently return zeros"
metadata:
  node_type: memory
  type: project
  originSessionId: f1519b98-d49a-4cee-b72f-15a9c22c2515
---

The PyInstaller exe and a source checkout do **not** share a data world:
`config.runtime_data_dir()` returns `<repo>/instance` from source but
`%LOCALAPPDATA%\COBLAST_data` when frozen — two separate `database_registry.sqlite`
files. The eToL post-BLAST stages are opt-in AND resource-gated: **human filter**
needs a registered human-genome DB, **contig ID** needs a reference rRNA DB, and
all of human-filter/assembly/contig-ID recover patient reads via `blastdbcmd` on
a `-parse_seqids` DB (else scan the stored source FASTA). If the frozen registry
lacks those DBs / a parse_seqids index / a readable source FASTA, every one of
those stages returns **zero** — no crash, just wrong-looking results.

**Diagnostic tell:** the raw BLAST is identical across builds; only the
resource-gated post-steps diverge. On 2026-07-02 a frozen-vs-source pair
reconciled exactly — 613 + 3614 (dedup) + 10162 (human) = 4408 + 9981 (dedup) + 0
= 14389 net hits both — proving the core pipeline was fine and only the human
filter (0 vs 10162 removed) had silently not run, cascading into dedup/assembly.

**Why:** `filter_human_hits` conservatively KEEPS unrecoverable reads and reports
0-removed + a `note`, and that note only rendered in the per-DB detail table, so
"couldn't run" looked like a clean "found none" in the summary card.

**How to apply:** the fix that landed this session — (1) `app.summarize_human_filter_warnings`
surfaces the note in the batch summary card; (2) `frozen_self_check.py` +
`run_COBLAST.py --self-check`, run as a post-build gate in `build_standalone_exe.py`,
drives the *frozen exe* through read-recovery + human-filter + CAP3 on a synthetic
sample and fails the build on regression (CAP3 assertion is conditional so
no-CAP3 machines still pass). Quick per-machine unblock for a user: point the exe
at the populated registry with `COBLAST_DATA_DIR=<repo>\instance`. Same bundling
class as [[pyinstaller-build-needs-biopython]]; the self-check is the packaged
sibling of [[coblast-verification-anchor]]'s spike-in.
