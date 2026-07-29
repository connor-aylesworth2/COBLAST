# Methodology evidence pack — extracted from the COBLAST+ codebase

Everything below is pulled from the repository at commit `e13e3b3` (branch `main`,
clean tree). Each fact carries its source file and line so it can be checked and
cited. Sections follow the dissertation outline's numbering.

**Read the "Gaps" box at the end of each section first.** Four of the nine
subsections describe work the codebase does not currently contain, and one
(§2.5.3) states a scope limitation that is factually wrong as written in the
outline.

---

## 2.1 Hardware / software environment

### Development and analysis machines

| Component | Value | Source |
|---|---|---|
| Workstation OS | Windows 11 Home, build 10.0.26200 (`Windows-11-10.0.26200-SP0`) | `docs/appendix_spikein.txt:9` |
| Workstation Python | 3.12.10 (source minimum 3.11) | `docs/appendix_spikein.txt:8`; `README.md:183` |
| Workstation threads used | 14 logical cores allocated to BLAST | `docs/appendix_spikein.txt:7` |
| Linux server | build host for the large reference databases, `/home/s2837738/COBLAST_2.0/blast_dbs`, 8 threads for `makeblastdb` | `scripts/build_etol_v_validation_db.sh` (recoverable at `git show 0db7356^:scripts/build_etol_v_validation_db.sh`) |
| CI runner | GitHub Actions `ubuntu-24.04`, Python 3.12 | `.github/workflows/ci.yml` |

### BLAST+ deployment

NCBI BLAST+ **2.17.0+** throughout. The Windows default search location is
`C:\Program Files\NCBI\blast-2.17.0+\bin` (`config.py:17`), overridable by the
`BLAST_BIN` environment variable or `--blast-bin`. CI downloads the pinned
`ncbi-blast-2.17.0+-x64-linux.tar.gz` from the NCBI FTP server, **verifies it
against the published MD5**, and asserts both `blastn -version` and
`makeblastdb -version` report `2.17.0+` before any test runs
(`.github/workflows/ci.yml`, steps "Install NCBI BLAST+" and "Verify BLAST+
installation"). Six executables are used: `blastn`, `blastp`, `blastx`,
`tblastn`, `makeblastdb`, `blastdbcmd`.

Remote BLAST is structurally disabled, not merely unused: `-remote` is in
`DISALLOWED_BLAST_OPTIONS` (`config.py:20`) and
`enforce_local_blast_only()` (`blast_runner.py:378`) rejects any generated
command containing it. No sequence data leaves the machine.

### Flask stack and dependency surface

The entire runtime dependency list is two packages (`requirements.txt`):

```
Flask>=3.0,<4.0
biopython>=1.83,<2.0
```

Biopython is used for FASTA validation and BLAST XML parsing only. **There is no
numpy, scipy, pandas, statsmodels or scikit-bio anywhere in the runtime.** Every
statistical procedure in §2.7 is implemented from first principles in
`static/etol_pie.js` — this is a deliberate constraint of the one-file
executable, and it is worth stating in the methodology because it means the
statistics are auditable line-by-line rather than delegated to a library version.

Development extras (`requirements-dev.txt`): `pytest>=8,<9`, `pyinstaller>=6,<7`.

The server binds `127.0.0.1` on port 5000 (`config.py:18-19`) and auto-increments
to the next free local port when 5000 is occupied.

### CAP3

CAP3 is **not redistributed** for licensing reasons. COBLAST+ auto-detects the
copy bundled with Unipro UGENE at `…\Unipro UGENE\tools\cap3`
(`config.py:353`, `_ugene_cap3_candidates`), or honours `CAP3_BIN`. Assembly is
optional throughout: a run that requests assembly with no CAP3 present continues
and reports that contigs were skipped (`docs/eToL.md:186-192`).

### Reference databases (all built locally, all searched locally)

| Database | Contents | Build |
|---|---|---|
| `ToL_rRNA` | SILVA SSU + LSU NR99, U→T converted, ~0.5 GB, `makeblastdb -parse_seqids` | user-built; registered under the `reference` category |
| `nt_human_9606` | *Homo sapiens* subset of NCBI `nt` (`blastdbcmd -db nt -taxids 9606`) — chosen over a bare GRCh38 assembly because curated mRNAs catch exon-junction-spanning brain RNA-seq reads a genome assembly misses | user-built; category `human` |
| `ToL_virus_val` | RefSeq viral release genomes **+** human GRCh38.p14 including the mitochondrion `NC_012920.1`; ~3–4 GB | `scripts/build_etol_v_validation_db.sh` |

The `ToL_virus_val` mixture is load-bearing, not incidental. The script's own
header explains why: the eToL-V dominant false positive was herpesvirus probes
recruiting reads from the human **mitochondrial** genome, which then assembled
into herpes-looking contigs. `contig_id.py` drops any contig whose closest
homolog contains `"Homo sapiens"` — but it can only do so if human and mito
sequence is *in* the reference database. A viruses-only reference would give
those artefact contigs no hit at all and pass them through as real virus. The
script hard-fails the build if the concatenated FASTA contains zero
`Homo sapiens` or zero `mitochondrion` headers.

> ⚠ **`scripts/build_etol_v_validation_db.sh` is currently a 0-byte file in the
> working tree.** Its 154 lines were emptied in commit `0db7356`
> ("Troubleshooting last PyInstaller build crashing issue"), which looks
> accidental — that commit is otherwise about PyInstaller. The content is intact
> in git history and recoverable with
> `git show 0db7356^:scripts/build_etol_v_validation_db.sh > scripts/build_etol_v_validation_db.sh`.
> Restore it before submitting; a methodology that cites a build script should
> not point at an empty file.

### SRA acquisition

SRA Toolkit (`prefetch`, `fasterq-dump --fasta`, `fastq-dump`) located via
`SRA_TOOLKIT_BIN` or a sibling `sratoolkit.<ver>` folder. The SRA-mart workbench
runs `prefetch` → `fasterq-dump --fasta` → `makeblastdb -parse_seqids` per
accession, landing each run BLAST-ready under `…\sra\<ACCESSION>\`. Only *run*
accessions (SRR/ERR/DRR) are accepted; study/experiment/sample accessions are
rejected so a whole multi-terabyte study cannot be queued by accident
(`sra_workflow.py:132`).

The `-parse_seqids` flag matters downstream: it is what lets the human filter
recover matched reads by id via `blastdbcmd` rather than falling back to a
linear scan of the source FASTA (`human_filter.py:63-91`).

### PyInstaller packaging

One-file Windows executable via `COBLAST.spec` / `build_standalone_exe.py`.
Bundled binaries: `blastn`, `blastp`, `blastx`, `tblastn`, `makeblastdb`,
`blastdbcmd`, plus `ncbi-vdb-md.dll` and `nghttp2.dll`, all into `blast/bin`.
Bundled data: `templates`, `static`, `sample_data`, `data`, `requirements.txt`.
`console=True`, `upx=False`, no `runtime_tmpdir` override.

Persistent data lives outside the executable at `%LOCALAPPDATA%\COBLAST_data`,
relocatable via a first-run folder picker, the Settings page, `--data-dir`, or
`COBLAST_DATA_DIR`.

Two filesystem constraints are enforced in code, both learned from failures:

- The data directory must be a **fixed NTFS volume**. exFAT/FAT (USB sticks)
  break `makeblastdb`, and removable drives remap letters after reboot.
  Guarded by `blast_incompatible_filesystem()` and
  `require_blast_capable_data_dir()` (`config.py:70`, `config.py:98`), with a
  live `makeblastdb` probe in `makeblastdb_probe_error()` (`config.py:133`).
- **No spaces in the path.** BLAST+ cannot build databases under `D:\My Data`.
  On Windows, BLAST-facing paths are stored using 8.3 short segments
  (`SCHOOL~1`) via `blast_safe_path()` while the registry keeps the true source
  path (`README.md:906-909`).

### Concurrency

`default_thread_count()` (`config.py:614`) reserves 1 core on machines with ≤4
logical cores and 2 on larger machines (`CPU_RESERVE_SMALL`/`CPU_RESERVE_LARGE`,
`config.py:605-606`). Batch databases are searched in parallel by a thread pool
(`run_jobs_concurrently`, `blast_runner.py:720`) — each BLAST is a separate OS
process, so threads give real parallelism despite the GIL.

One deliberate serialisation: the human-genome filter takes a **global lock**
(`_HUMAN_BLAST_LOCK`, `human_filter.py:50`). Thirty concurrent workers all
scanning the same multi-GB human database thrashed the page cache and pinned CPU
at 2–4%; serialised, the first scan warms the cache and the rest run warm with
all cores. This changes *when* a search runs, never its result — a safe thing to
state in the methodology.

---

## 2.2 COBLAST+ architecture as instrument

### The eToL cellular preset net path, in execution order

This is the figure you want for "Fig: COBLAST+ architecture". Every arrow below
is a named function, and `tests/test_spike_in_control.py:11-17` documents this
same chain as the canonical description of the path.

1. **Panel assembly** — `etol_search_fasta(preset)` (`etol_summary.py:341`)
   concatenates the microbial panel with its housekeeping control probes. The
   controls are appended automatically and are never separately selectable, so a
   single BLAST yields both the net and the normalisation denominator.
2. **The net** — `run_blast_probe_panel()` (`blast_runner.py:694`) runs
   `blastn -task megablast` with BLAST's default scoring, **no identity filter
   and no coverage filter**. The sole parameter override is lifting
   `max_target_seqs` to `5000000` (`EXACT_MATCH_MAX_TARGET_SEQS`,
   `blast_runner.py:112`) so a probe matching many reads in a deep library is
   counted in full rather than truncated at BLAST's default of 500.
3. **Net gate** — `filter_net_probe_hits()` (`app.py:187`) restricts hits to the
   active panel and drops any hit with **E ≥ 0.01**, or whose E-value cannot be
   parsed. `ETOL_EVALUE_THRESHOLD = 0.01` (`app.py:184`) is the paper's own net
   cutoff from `Abundance_ToL.py`. This matters: BLAST's default E-value is 10,
   so "default megablast" alone would have been far too permissive.
4. **Control split** — control-probe hits are separated **before** the human
   filter. The control reads are human by design; filtering them would zero the
   normalisation denominator.
5. **Control counting** — `count_control_reads()` (`app.py:252`) de-duplicates
   control hits to one read per best control probe, then counts per probe
   including zeros. Mirrors `Abundance_count.py` §2. Without this, a read
   matching several redundant control probes would inflate the denominator.
6. **Secondary human filter** (optional) — `filter_human_hits()`
   (`human_filter.py:227`). Matched reads are recovered by `sseqid`, then
   `blastn -task megablast -max_target_seqs 1 -outfmt "6 qseqid bitscore"`
   against the selected human database. A read is dropped when its best human
   alignment exceeds **150 bits** (`HUMAN_BITSCORE_THRESHOLD`,
   `human_filter.py:42`). Bitscore is the *sole* criterion — the E-value is set
   to `1e9` (`DEFAULT_HUMAN_EVALUE`, `human_filter.py:41`) specifically so it
   never pre-filters anything. Reads that cannot be recovered are **kept**, never
   dropped on a guess, and the failure is reported.
7. **Cross-probe de-duplication** — `deduplicate_reads_to_best_probe()`
   (`app.py:210`). ~38% of the panel's probes share sequence with at least one
   other, so the same read can be recovered several times. Each read is allocated
   to a single probe by bitscore → percent identity → query coverage → probe id
   (lexical, as a deterministic final tie-break, `_similarity_rank`,
   `app.py:236`). Order is human-filter **then** de-duplicate, matching the
   paper's text.
8. **Host-cell normalisation** — `compute_host_cells()` (`etol_summary.py:361`):
   average each control gene's probe counts, take the mean across genes, divide
   by `HOST_TRANSCRIPTS_PER_CELL = 50.0` (`etol_summary.py:305`). Abundance is
   then reported both raw and as **reads per host cell**. When no control reads
   are found, normalisation reports `n/a` rather than zero.
9. **Contig assembly** — reads grouped per taxon by `group_read_ids_by_taxon()`
   (`etol_summary.py:516`), assembled by `Cap3Assembler.assemble()`
   (`assembler.py:130`). CAP3's own `-o`/`-p` overlap defaults are used unless
   deliberately tuned, so results match the paper's web-server runs.
   `MIN_READS = 2` — a single read is its own singlet and produces no contig.
   Each assembly runs in its own temp directory because CAP3 is filename-driven
   and concurrent assemblies would otherwise collide on output names.
10. **Re-probing** (optional, **off by default**) — `reprobe_and_reassemble()`
    (`contig_id.py:424`). See §2.5.3; this is where the outline needs correcting.
11. **Species identification** — `identify_contigs()` (`contig_id.py:262`).
    (a) `name_contigs()` batches every contig into one
    `blastn -task megablast` against the registered reference database at
    E ≤ `1e-5` (`DEFAULT_NAME_EVALUE`, `contig_id.py:48`), keeping the
    best-bitscore hit's `stitle` as the closest homolog.
    (b) Any contig whose homolog string contains `"Homo sapiens"`
    (`HUMAN_HOMOLOG_MARKER`, `contig_id.py:73`) is dropped — this catches human
    rRNA reads that slipped the genome-level filter.
    (c) `confirm_contig_reads()` re-BLASTs each taxon's reads against its own
    contigs (`-subject`, no index needed) and counts distinct reads aligning at
    **≥99% identity** (`DEFAULT_CONFIRM_IDENTITY_PCT`, `contig_id.py:45`) — 99%
    rather than 100% tolerates one sequencing error across a short rRNA read.
12. **Summary and export** — `build_etol_probe_summary()` (`etol_summary.py:542`)
    and `build_etol_matrix()` (`etol_summary.py:688`) feed the results tables,
    the SVG heatmap, the pie/stacked-bar chart, and the CSV/TSV exports.

Synthetic ids `c0, c1, …` are used for contigs during naming because CAP3
restarts its `Contig1` numbering per taxon; they map back to `(taxon, contig)`.

### Interface

Flask, local-only, NCBI-inspired palette, clinician-facing controls separated
from advanced BLAST parameters. Pages: `/` (single BLAST), `/databases`,
`/batch-blast`, `/sra`, `/settings`. Default per-search wall-clock timeout 3,600 s
(`DEFAULT_TIMEOUT_SECONDS`, `blast_runner.py:93`), capped at the same value.
Query validation accepts IUPAC nucleotide ambiguity codes and converts U→T;
`MAX_FASTA_RECORDS = 1500` (`blast_runner.py:96`) accommodates the 1,021-probe
full panel with headroom.

### Database registry

SQLite (`instance/database_registry.sqlite` from source; under the data directory
when frozen). Records display name, type (`nucl`/`prot`), category, source FASTA
path, BLAST prefix path, description, and an availability status
(`available`/`missing`/`invalid`). Categories in use: `human`, `viral`,
`eToL-V`, `toy`, `custom`, `reference`. Program/database compatibility is enforced
by filtering the dropdown on database type after program selection.

### Panel composition (verified by running the loaders)

| Preset | Probes | Taxa | Class codes | Domains | Probe length | Controls appended |
|---|---|---|---|---|---|---|
| eToL Full | 1,017 | 120 | 25 | 7 | 64 nt | `PGK1_2`, `PGK1_3`, `hNSE_2`, `hNSE_3` |
| eToL Quick | 120 | 120 | 25 | 7 | 64 nt | same four |
| eToL-V | 115 | 50 | 4 | 1 (Viruses) | 120 nt | `PGK1_1`, `PGK1_2` |

Class-code → domain mapping (`ETOL_DOMAIN_BY_LETTER`, `etol_summary.py:72`),
following the paper's scheme: **A** Archaea (91 probes); **B0–B6** Bacteria (171);
**C1–C4** Chloroplastida (258); **D** Amoebozoa (64); **E0** basal Eukaryota (46);
**F0–F6** Fungi (241); **H0–H3** Holozoa/Metazoa (146).

eToL-V classes: **V-HHV** herpesviruses (77 probes), **V-HCoV** coronaviruses
(21), **V-HPV** papillomaviruses (10), **V-HAdV** adenoviruses (7), spanning 24
virus tokens resolved to 50 gene-level taxa (e.g. `V-HAdV_AdC_penton`).
Header grammar is `Class_Taxon_Subunit_Index`.

**A known, documented net limitation to declare:** megablast requires a 28-base
unambiguous window to seed. A probe whose ambiguous bases leave no such window
finds nothing and is silently dropped. In the bundled microbial panel this is
exactly one probe, `F3_Gpolymorpha_18S_7` (`docs/eToL.md:31-34`). The whole panel
is run at megablast rather than `blastn-short` because megablast is far faster on
whole-SRA databases.

### Design matrix (condition labels)

Name-based condition inference is unreliable — auto-generated SRA database names
(`SRA <acc> reads`) carry no diagnosis, and substring matching mislabels samples.
An uploaded design matrix is therefore **authoritative** and suppresses the guess
entirely (`design_matrix.py`; applied via
`etol_summary.build_etol_matrix(..., condition_index=...)`).

Format: CSV or TSV, UTF-8 (leading BOM tolerated), header row with `sample` and
`condition` columns, case-insensitive, any order. `sample` binds by SRA accession
first, then by exact database display name. One row per sample; a duplicate is an
error. A malformed file is rejected **on the form, before the BLAST run**.
Unmatched samples render as a neutral "unlabeled" swatch and are listed in a
warning under the heatmap. Template downloadable at `/design-matrix-template.csv`.

---

## 2.3 Datasets

### EBB — fully specified in the repository

`data/etol_v_sra_crosswalk.csv` is the authoritative structure file: 35 rows,
columns `srr, srx, region, diagnosis, sample_name`.

- **Study**: SRP398685, Edinburgh Brain Bank (`etol_validation.py:9`).
- **Runs**: `SRR21676099`–`SRR21676133` (35 consecutive run accessions).
- **Experiments**: `SRX17674433`–`SRX17674467`.
- **Individuals**: 9 (`SD001-17`, `SD005-19`, `SD014-17`, `SD025-19`,
  `SD030-18`, `SD032-17`, `SD035-15`, `SD037-18`, `SD042-18`).
- **Regions**: 4 — HYPO (hypothalamus, n=9), HPC (hippocampus, n=9),
  AMYG (amygdala, n=9), BA24 (anterior cingulate, **n=8**).

**The repeated-measures structure, derived from the file:**

| Individual | Diagnosis | n samples | Regions |
|---|---|---|---|
| SD001-17 | AD | 4 | AMYG, BA24, HPC, HYPO |
| SD005-19 | AD/VaD | 4 | AMYG, BA24, HPC, HYPO |
| SD014-17 | AD | 4 | AMYG, BA24, HPC, HYPO |
| SD025-19 | AD | 4 | AMYG, BA24, HPC, HYPO |
| SD030-18 | CONTROL | 4 | AMYG, BA24, HPC, HYPO |
| SD032-17 | AD/LBD | 4 | AMYG, BA24, HPC, HYPO |
| SD035-15 | CONTROL | 4 | AMYG, BA24, HPC, HYPO |
| SD037-18 | AD | 4 | AMYG, BA24, HPC, HYPO |
| **SD042-18** | **CONTROL** | **3** | **AMYG, HPC, HYPO — no BA24** |

So 9 × 4 = 36 minus one = **35**. The unbalanced individual your outline
anticipates is `SD042-18`, a **control**, missing **BA24**. This is worth naming
explicitly: the single missing cell is in the control arm, which is exactly where
a listwise-deletion approach would cost you the most power.

Sample-level diagnosis counts: AD 16, CONTROL 11, AD/VaD 4, AD/LBD 4.
Individual-level: 4 AD, 3 CONTROL, 1 AD/VaD, 1 AD/LBD.

The `ETOL_CONDITION_PATTERN` regex (`etol_summary.py:61`) recognises
`AD/LBD`, `AD/VaD`, `CTRL`, `CONTROL`, `LBD`, `VaD`, `AD`, ordering longer
combinations first so `AD/LBD` is preferred over a bare `AD`.

> **Gap — EBB clinical covariates.** The repository contains region, diagnosis
> and a de-identified individual code, and **nothing else**. There is no age,
> sex, post-mortem interval, Braak stage, Thal phase or APOE genotype anywhere in
> the codebase (verified by grep across `.py`, `.md`, `.js`, `.html`, `.csv`).
> Your outline's Table 3 ("EBB by individual — region, dx, age, sex, PMI, Braak,
> Thal, APOE") will need those columns sourced from the EBB metadata or the
> source publication and added by hand.

> **Gap — Kohen 2014 is entirely absent.** No accession, no crosswalk, no
> metadata file, no code path references it. Grep for "Kohen", "GSE", "PRJNA"
> returns nothing outside the SRP398685 references. The 29 hippocampus samples,
> the one-per-individual design, the 29–95 age range and the 60–67 age gap are
> all facts you will need to state from the source paper. Nothing about that
> dataset can be cited to this codebase.

---

## 2.4 VERIFICATION protocol

This section is the best-evidenced in the whole outline. The codebase contains a
three-layer verification harness, and the top layer is a **required CI gate**.

### Layer 1 — Synthetic spike-in positive control

`tests/test_spike_in_control.py` (411 lines, extensively self-documenting; its
module docstring at lines 1–56 is written specifically for this dissertation
section and can be quoted more or less directly).

**The idea, in its own words**: a spike-in control is the microbiology notion of a
positive control applied to software — construct a sample whose true composition
you already know, run it through the instrument untouched, and check the
instrument returns the known answer. The "instrument" is the *real* eToL net
path, imported from `app.py` and `etol_summary.py`, **not a reimplementation**
(`run_net_pipeline`, line 188, explicitly mirrors the batch route).

**Design of the synthetic sample** (constants at lines 83–89):

| Parameter | Value | Rationale (from the code) |
|---|---|---|
| `PRESET` | `etol_quick` | one probe per species: fast and unambiguous to reason about |
| `READ_LEN` | 150 nt | realistic Illumina-like read; the 64 bp probe sits in the middle inside random flanks |
| `SPIKE_LEVELS` | 20, 10, 3 reads | deliberately straddle the paper's 3–5 reads/host-cell cellular cutoff |
| `CONTROL_READS` | 100 per control probe | gives host cells = mean(100,100)/50 = **exactly 2.0**, a round denominator |
| `NOISE_READS` | 500 random-sequence reads | background that must produce zero detections |
| absent species | 1 panel species, 0 reads | the specificity target |
| `RNG_SEED` | 20240601 | fixed, so the synthetic input is byte-identical on every machine |

Species selection (`choose_spike_probes`, line 99) is deterministic — the panel
order is fixed — restricted to pure-ACGT probes so megablast can always seed on
the 64 bp marker, and prefers one species per domain so the specificity claim
does not lean on the de-duplication step (though de-duplication would absorb
conserved-region cross-hits anyway).

**The three claims asserted** (`check`, line 242):

1. **Recovery / sensitivity** — every spiked species is detected and its net read
   count *equals* what was spiked. Abundance is preserved end to end.
2. **Specificity** — the absent species is not reported, the 500 noise reads
   produce no spurious species, and exactly three species come back
   (`species_detected == 3`).
3. **Normalisation** — reported normalised abundance equals net reads ÷ host
   cells, with host cells recovered exactly (2.0) from the spiked controls.

**Correctness is independent of biology**: the docstring states this outright
(lines 41–43) — "it proves the pipeline reports what is present at the abundance
present and rejects what is absent. Validation against real data (reproducing
Lathe/Veso) is a separate, later step."

**Recorded run** — `docs/appendix_spikein.txt` (UTF-16, note the encoding if you
paste it). It carries a provenance block (`provenance()`, line 360) recording git
commit, `blastn`/`makeblastdb` versions, which binary resolved, thread count,
Python version, platform and RNG seed. The rationale for each field is given in
the docstring: megablast tie-breaking shifts across BLAST releases, and the
thread count feeds tie-breaking into de-duplication.

Results table from that run:

| Species | Spiked | Recovered | Norm. expected | Norm. observed | Result |
|---|---|---|---|---|---|
| `A_Hsalinarum_16S` | 20 | 20 | 10.00 | 10.00 | PASS |
| `B0_DesulfobacteraceaecLaKi_16S` | 10 | 10 | 5.00 | 5.00 | PASS |
| `C1_Pbrachykentron_18S` | 3 | 3 | 1.50 | 1.50 | PASS |
| `D_Dtestaceum_18S` | 0 | 0 | — | — | PASS (absent) |
| random noise (500 reads) | 0 | 0 | — | — | PASS (no spurious) |

Host cells from control reads = 2.0 (expected 2.0). Three domains represented
(Archaea, Bacteria, Chloroplastida), so the three spiked markers are maximally
divergent.

> ⚠ **Regenerate this artifact before submission.** The recorded run's provenance
> line reads `git commit: 6603714… (DIRTY - uncommitted changes)`. The code's own
> comment explains why that matters (lines 343–346): "a positive control verifies
> *committed* code; a dirty tree means the table describes code that lives
> nowhere". Re-run `python tests/test_spike_in_control.py > docs/appendix_spikein.txt`
> from the clean tree at your submission commit so the appendix documents code a
> reader can actually check out. (It will also fix the UTF-16 encoding if you
> redirect from a UTF-8 shell.)

### Layer 2 — Packaged-executable self-check

`frozen_self_check.py`. This exists because of a real defect class: three stages
depend on bundled binaries and PyInstaller `_MEIPASS` resource resolution, and in
a mis-bundled build they silently returned **zeros** rather than failing —
patient-read recovery via `blastdbcmd`, the secondary human filter, and CAP3
assembly. A `pytest` run from source cannot catch that; only running the packaged
`.exe` can.

It builds a tiny synthetic sample of known composition (`RNG_SEED = 20240717`,
4 surviving "microbial" reads, 3 reads also planted in a synthetic human database)
and asserts:

- read recovery works via the `-parse_seqids` id index with **zero** unresolved
  reads — the exact thing that silently failed;
- the human filter removes the planted human reads and keeps the rest;
- CAP3 assembles the survivors into ≥1 contig — but only when a CAP3 binary
  resolves, since a build without CAP3 is legitimate.

Wired to `run_COBLAST.py --self-check` and run as a **post-build gate** against
the freshly built `dist/COBLAST.exe`, so a bundling regression fails the build
instead of a user's analysis.

### Layer 3 — Smoke test and unit suite

`smoke_test.py` creates toy nucleotide and protein databases with `makeblastdb`,
then runs `blastn`, `blastp`, `blastx` and `tblastn` through the shared backend.
The launcher runs it before opening the web interface, proving BLAST+, FASTA
validation and output parsing work before a user sees the UI.

**Table 2.x — COBLAST+ automated test suite.** 150 tests across 14 pytest
modules, collected at commit `e13e3b3` (`python -m pytest --collect-only -q`).

| Module | What it verifies | Tests |
|---|---|---|
| `test_blast_runner` | FASTA validation, tabular parsing, parameter building, local-only command guard, exact-match probe overrides; no BLAST+ binary invoked | 43 |
| `test_summaries` | APOE and eToL probe-count aggregation: read deduplication, E-value net gate, host-cell normalisation, matrix assembly, condition sorting | 29 |
| `test_database_registry` | SQLite database registry: removal, path handling, route behaviour | 15 |
| `test_assembler` | CAP3 contig assembler and its binary resolver (`subprocess` monkeypatched) | 13 |
| `test_contig_id` | Contig species identification and confirmed abundance (blastn helpers monkeypatched) | 12 |
| `test_design_matrix` | Design-matrix parser: strict CSV/TSV format, accession and display-name index | 10 |
| `test_sra_workflow` | Single-walk SRA project helpers: fetch-script generation, file discovery | 8 |
| `test_etol_validation` | eToL-V confusion matrix; reproduces the published TP 9 / FP 1 / FN 35 / TN 411 | 7 |
| `test_human_filter` | Matched-read recovery for the secondary human filter | 4 |
| `test_settings_browse` | Settings folder picker: cancel and missing dialog leave the typed path usable | 3 |
| `test_launcher_port` | Readiness check gating browser launch on Flask actually listening | 3 |
| `test_batch_progress` | Batch progress endpoint mirrors the in-memory registry | 1 |
| `test_data_dir_prompt` | First-run prompt accepts a typed path without the native dialog | 1 |
| `test_spike_in_control` | Synthetic spike-in positive control through the real eToL net path; also run as a separate required CI step | 1 |
| **Total** | | **150** |

One further module, `tests/test_etol_pie_stats.js` (125 assertions, Node), is
run manually and is **not** part of the pytest total or of CI — see the gap note
below.

**CI** (`.github/workflows/ci.yml`) runs on every push and PR to `main`:
checkout → Python 3.12 → BLAST+ 2.17.0 with MD5 verification → version assertion
→ `pip install -r requirements-dev.txt` → `pip check` → `pip freeze` →
`pytest -vv -ra --tb=long --junitxml=…` → **and then a separate required step,
"Run required BLAST+ spike-in validation", executing
`python tests/test_spike_in_control.py`**. The spike-in is not merely a test that
may skip; it is a build gate that must pass with real BLAST+ present.

> **Gap — the statistics tests are not in CI.** `tests/test_etol_pie_stats.js`
> (125 assertions covering midranks, exact and permutation rank tests,
> Kruskal–Wallis, BH, PERMANOVA, Bray–Curtis, and the whole parametric arm
> against closed forms) runs only manually via `node tests/test_etol_pie_stats.js`.
> The CI workflow has no Node step. Either add one, or state in the methodology
> that the statistical routines are verified by a separately-invoked harness —
> do not imply CI covers them.

---

## 2.5 VALIDATION protocol

### 2.5.1 Reproduction targets and acceptance criteria

> **Gap — no pre-registered acceptance criteria exist in the codebase.** There is
> no criteria table, no claim-type taxonomy, no match/partial/miss rubric. The
> only thing resembling a stated target is the fidelity line rendered in the
> results panel: *"Fidelity target (Veso, Fig 9): TP 9 / FP 1 / FN 35 / TN 411"*
> (`templates/batch_results.html:432-433`), plus the assertions in
> `tests/test_etol_validation.py`. Your "Table: Pre-registered acceptance
> criteria per target" — which you correctly identify as the highest-leverage
> table in the section — must be authored from scratch. Write it before you run
> anything, and date it.

What the codebase *does* give you, which constrains how the criteria can be
phrased:

- **Grading level.** Genus/clade grading is directly supported. `_species()`
  (`etol_summary.py:145`) strips the class prefix and rRNA-unit suffix and
  expands abbreviated binomials (`B0_Tmaritima_16S` → `T. maritima`), leaving
  non-binomial labels untouched (`AcidobacteriaKBS96`, `HSV1_gB`). Domain and
  class-code columns are carried through every summary and export, so ordinal and
  set-membership claims can be scored at domain, class, or taxon level without
  new code. `species_from_homolog()` (`contig_id.py:76`) returns SILVA's *last*
  taxonomic rank from the homolog string — so contig-based calls arrive as SILVA
  taxonomy strings, not NCBI taxids.
- **Available magnitudes.** Per-probe counts, per-species counts, reads-per-host-
  cell, contig-confirmed read counts, and closest-homolog percent identity — all
  exported per sample including zeros (`ETOL_SPECIES_EXPORT_COLUMNS`,
  `ETOL_PROBE_EXPORT_COLUMNS`, `etol_summary.py:100-125`).
- **A stability caveat for direction+significance claims.** Re-probing rewrites
  species names and confirmed abundance in place, so on/off results are not
  directly comparable; the documentation is explicit that one setting must be
  applied uniformly across every sample in a study (`docs/eToL.md:198-202`).

### 2.5.2 eToL-V confusion matrix

`etol_validation.py` is a complete, tested specification. Its docstring
(lines 1–28) is the reconstructed specification your outline asks you to state.

**Bundled data:**
- `data/etol_v_wgs_truth.csv` — 700 rows of `srx,virus,count`, whole-genome-shotgun
  read counts from the original eToL workflow (20 viruses × 35 samples).
- `data/etol_v_sra_crosswalk.csv` — the 35-row SRR↔SRX map.

**A verified, not assumed, quirk:** the SRR↔SRX mapping is **inverted** —
`SRR21676133 ↔ SRX17674433`, `SRR21676131 ↔ SRX17674435`. The lowest SRR maps to
the highest SRX. `test_crosswalk_is_inverted_and_complete`
(`tests/test_etol_validation.py:19`) pins this explicitly, because a same-suffix
assumption would silently mis-join every sample and still produce a plausible-
looking matrix.

**The scoring specification** (`compute_confusion`, line 140):

| Element | Specification | Source |
|---|---|---|
| Ground-truth positive | WGS count **> 0** — binary classification; >0 is the only threshold yielding the reference's 44 actual-positives | `etol_validation.py:20-21`; asserted at `tests/test_etol_validation.py:30-38` |
| Prediction positive | any of a virus's probes has a **validated (contig-confirmed) hit > 0** — i.e. compared against the post-validation heatmap, not the raw net | `etol_validation.py:22`, `stage="validated"` |
| Universe | 13 WGS virus rows × 35 samples = **455** cells | `VESO_UNIVERSE`, line 48 |
| Out-of-universe | any validated prediction outside the universe becomes an extra FP cell, one per (virus token, sample) — the lone HPV45 hit → **456** | lines 214–223 |
| Exclusions | `SARSCoV2`, `SARSCoV` dropped entirely, no WGS data exists for them | `EXCLUDED_VIRUS_TOKENS`, line 67 |

**The 13-virus universe**, in order: Adenovirus C, COV_229E, HHV1_HSV1,
HHV2_HSV2, HHV3_VZV, HHV4_EBV, HHV5_CMV, HHV6A, HHV6B, HHV7, HHV8, HPV6, HPV16.

**Two documented quirks kept deliberately for fidelity** (lines 24–27):
1. **HPV6 is in the universe but the panel has no HPV6 probe**, so it is always
   predicted negative. `universe_taxa()["HPV6"]` returns an empty frozenset,
   pinned by `test_universe_taxa_quirks`.
2. **Adenovirus A and Adenovirus 54 are omitted** from the universe.

The module states plainly that these can be swapped out — "Swap `VESO_UNIVERSE`
for a corrected set if you would rather not reproduce those" — which is exactly
the configured-to-match-the-reference vs corrected-ground-truth distinction your
outline draws. **Reproduction is run with COBLAST+ configured to match the
reference implementation**, and the alternative is a one-constant change.

**The reproduction is verified, not asserted.** `test_reproduces_veso_confusion_matrix`
(`tests/test_etol_validation.py:70`) constructs a matrix mirroring the reference's
surviving calls — adenovirus C penton in 9 of the 30 WGS-adenovirus-C-positive
samples, plus the single HPV45 L1 hit in `SRX17674444` — and asserts:

```
(tp, fp, fn, tn) == (9, 1, 35, 411)     n == 456
accuracy 0.9211    precision 0.90    recall 0.2045    F1 0.33
```

It further asserts a 456-cell per-cell breakdown summing to those totals, and
that every true positive is an adenovirus-C cell with a validated hit. The module
docstring states the construction was "reverse-engineered from her dissertation
and verified to reproduce 9/1/35/411 exactly" — so the reference matrix **is**
recomputed from the ground-truth CSV rather than taken as published, which is the
claim your outline wants to make.

**On the threshold-sensitivity analysis.** The mechanism you need already exists,
but it is not a read-count threshold — it is a **stage** switch. `stage="validated"`
scores the contig-confirmed layer (the faithful comparison); `stage="raw"` scores
the raw net hits (line 148–155). When no validated layer is present, the function
degrades gracefully to `raw` and records the substitution in the returned dict.
Running both and reporting them side by side is the sensitivity analysis; phrase
it as *net vs contig-validated*, not as a numeric threshold sweep, because that is
what the code actually varies. Every cell also carries both `raw_hits` and
`confirmed_hits`, so a false negative reveals **where** it was lost: `raw = 0`
means the net missed it; `raw > 0, confirmed = 0` means contig
assembly/identification dropped it (lines 199–202).

Robustness paths are tested too: samples with no crosswalk entry score nothing,
return all-zero counts with `accuracy = None`, and are listed in
`unmatched_samples` rather than crashing (`test_no_ground_truth_overlap_is_graceful`).
SARS predictions are confirmed not to count as false positives
(`test_sars_predictions_are_excluded`).

**Outputs**: in-app confusion panel with TP/FP/FN/TN grid plus accuracy,
precision, recall, F1 and N (`templates/batch_results.html:412-467`); per-cell
CSV/TSV export (`result_store.etol_confusion_rows_as_delimited`, line 297); and
publication figures via `scripts/plot_etol.py`, which reuses
`compute_confusion()` so the figure cannot drift from the web view and is
deliberately kept out of the bundled executable (it is the only thing pulling in
matplotlib).

### 2.5.3 ★ SCOPE LIMITATION — the outline is factually wrong here

Your outline states: *"COBLAST+ builds contigs and identifies them against SILVA
SSU+LSU but does not generate new probes from LSU sequences and re-BLAST."*

**The first half of that is right; the second half is not.** COBLAST+ *does*
generate new probes from contigs and re-BLAST. `reprobe_and_reassemble()`
(`contig_id.py:424`) implements the paper's Box 3: it takes each taxon's
most-abundant assembled contig, BLASTs it back against the **same patient
database** as a fresh probe, pulls reads the original 64-mer net missed,
human-filters those new reads, and re-assembles the taxon from original + new
reads. Parameters: `REPROBE_TOP_CONTIGS = 1` (the paper used the two most
abundant; one is used because each extra contig is a full-length, highly-conserved
rRNA query against the whole patient library and is the dominant cost),
`REPROBE_EVALUE = 0.01` — the same gate as the net, so re-probing is never more
permissive than the search that found the taxon — and
`REPROBE_MAX_TARGET_SEQS = "100000"`. It is **off by default** and route order is
assemble → reprobe → identify.

**What is genuinely missing** — and what the scope limitation should actually
say — is the **23S/28S (LSU) confirmation/disambiguation step**. This is recorded
as an outstanding gap: *"23S/28S (+mtDNA) confirmation step to disambiguate
redundant probes (Bonferroni concern with >1000 probes); semi-manual in the
paper, no hook yet."* The reference database build *does* include LSU — SILVA
SSU+LSU NR99 — with the explicit note that "SSU alone spans all domains A–H
(16S pro + 18S euk); LSU adds 23S/28S for the validation step". So the LSU
sequence is present and searchable; what is absent is the pass that uses LSU to
disambiguate among redundant probes.

**Suggested replacement wording:**

> COBLAST+ assembles matched reads into contigs and identifies them against a
> local SILVA SSU+LSU NR99 reference, and it implements the paper's contig
> re-probing pass (contig-as-new-probe, re-BLASTed against the same patient
> library). It does **not** implement the paper's separate 23S/28S disambiguation
> step, in which redundant probe assignments are resolved by a targeted
> large-subunit confirmation. Because that step is the paper's own replication
> mechanism — the one that substituted for formal multiple-testing correction
> across a >1,000-probe panel — validation here is scoped to the 16S/18S layer,
> and the 23S/28S reprobe is named as Future Work.

Two further deviations belong in the same paragraph, since both are honest
divergences from the source method:

1. **Read grouping granularity.** The paper's fast path groups reads per
   phylogenetic group A–H; COBLAST+ groups **per taxon**, a finer and defensible
   deviation.
2. **Fixed human cutoff.** The paper varied its human-filter bitscore cutoff by
   each dataset's mean read length (>160 MSBB, >126 Rockefeller, >100 Miami);
   COBLAST+ fixes **150**, the value the paper applied to brain, and documents
   that libraries with very different read lengths may warrant a different value
   (`human_filter.py:33-40`).

Additionally: re-probing is **one round**. The paper hints at iterating until no
new reads are recovered; the implementation does a single round by choice.

This is the content for your "Fig: eToL net path as implemented, annotated where
COBLAST+ diverges". The divergence annotations are: E < 0.01 gate made explicit
(the paper's, but easy to lose); `max_target_seqs` lifted; megablast for the whole
panel with one probe silently unseedable; per-taxon rather than per-group
assembly; single re-probe round with the top-1 contig; fixed 150-bit human cutoff;
**no 23S/28S disambiguation pass**.

---

## 2.6 APPLICATION protocol

> **Gap — this is the largest one. None of §2.6 exists in the codebase.**
>
> - **2.6.1 Per-probe age analysis.** Not implemented. There is no age variable
>   anywhere in the codebase, no young/elderly grouping, no continuous-age model.
>   The Kohen dataset is absent entirely (§2.3). The existing significance tests
>   operate on the **7 Tree-of-Life domains**, not on the probe set — a
>   ~1,017-element multiple-comparison problem the current BH correction has never
>   been applied to. (The BH implementation itself, `bh()` in
>   `static/etol_pie.js:342`, is size-agnostic and would work; the per-probe
>   *values* are exported per sample in the eToL Probe Counts CSV, including
>   zeros, so the input matrix exists. What does not exist is anything that runs
>   a test per probe.)
> - **2.6.2 Overlap test.** Not implemented. There is no encoded shortlist of the
>   2023 paper's AD-overabundant taxa, and no set-intersection scoring anywhere.
> - **2.6.3 Regional burden with individual-level clustering.** Not implemented,
>   and this is the sharpest mismatch. Every test in `static/etol_pie.js` is an
>   **independent-samples** test. The per-domain path groups samples by
>   design-matrix condition and permutes labels freely; PERMANOVA permutes labels
>   across all samples. **Nothing accounts for the fact that the 35 EBB samples
>   come from 9 individuals.** There is no paired test, no Friedman, no mixed
>   model, no random intercept, and no restricted (within-individual) permutation
>   scheme. Run as-is on EBB, the per-domain tests would treat 35 samples as 35
>   independent replicates, which is precisely the inflation your §2.6.3 is
>   designed to avoid.
>
> The honest framing for the methodology: the application-leg statistics are
> **specified in this chapter and executed outside the tool**, using COBLAST+'s
> exported count matrices as input. The tool produces the measurements; the
> clustered/repeated-measures models are applied downstream.

**What the tool exports, i.e. what your application analyses actually consume:**

| Export | Shape | Route / source |
|---|---|---|
| eToL Probe Counts CSV/TSV | one row per probe per sample, **including zeros**: `Sample/Database, Probe, Species/Taxon, Class, Domain, Exact hits, Reads per host cell` | `ETOL_PROBE_EXPORT_COLUMNS`, `etol_summary.py:117` |
| eToL Species Summary CSV/TSV | one row per species per sample, including zeros, plus `Est. host cells`, `Reads per host cell`, `Closest homolog (contig)`, `Confirmed reads (contig)` | `ETOL_SPECIES_EXPORT_COLUMNS`, `etol_summary.py:100` |
| Matrix JSON | `/batch-results/<id>/etol-matrix.json?level=species` — rows × cols with `hits`, `confirmed`, per-column `sample`, `condition`, `host_cells` | `build_etol_matrix`, `etol_summary.py:688` |
| Domain composition CSV | per scope × domain: reads, percent, reads per host cell — plus a second block with p, q(BH), medians, n, and the parametric statistic per domain | `exportCsv`, `static/etol_pie.js:1033` |
| Contig multi-FASTA | `/batch-results/<id>/etol-contigs.fasta`, headers carry `\|confirmed=\|homolog=` | `app.py` |
| Confusion matrix CSV/TSV | per-cell: `Result, Virus (WGS), Sample, SRX, …, raw_hits, confirmed_hits` | `result_store.py:297` |

Column headers `Exact hits` and `Total exact probe hits` are deliberately kept
stable across versions to avoid breaking downstream plots — a small but real
reproducibility point worth a sentence.

---

## 2.7 Statistical approach

All of this lives in `static/etol_pie.js`, hand-implemented with no statistical
library, and is verified against closed forms in `tests/test_etol_pie_stats.js`.

### The rank test — one function, two tests

`rankStat()` (line 80) computes S = Σ<sub>g</sub> R<sub>g</sub>² / n<sub>g</sub>.
The comment above it explains why this single statistic covers both designs:

- **Two groups**: S is a parabola in the first group's rank sum whose vertex sits
  exactly at the null mean, so P(S ≥ S_obs) *is* the two-sided Wilcoxon rank-sum
  p — no separate two-sided handling.
- **More than two groups**: S is Kruskal–Wallis H up to constants that are fixed
  under permutation (the tie correction among them), so the permutation p is
  Kruskal–Wallis's.

**Exact vs sampled** (`rankTestP`, line 96): with two groups and
C(n₁+n₂, n₁) ≤ 100,000 the label space is **enumerated exhaustively** —
lexicographic combination generation, p = hits/trials, an exact p-value.
Otherwise 20,000 sampled permutations (`PERMUTATIONS`, line 65) via Fisher–Yates
with a **fixed-seed LCG** (seed `20220317`, line 124) so a p-value cannot change
between two renders of the same figure, and p = (hits+1)/(trials+1) — the
add-one correction that keeps a sampled p strictly positive. A 4×6 design has
~2.3×10¹² label assignments, which "neither needs nor admits enumeration".

Justification given in code (lines 92–95, 690–692): permuting is exact under
ties; eToL group sizes are far too small for the normal/chi-square approximations
the textbook formulae rely on; and read counts are skewed, zero-inflated and
small-n, so no distributional assumption is safe.

### Multiple-testing correction

`bh()` (line 342) — Benjamini–Hochberg adjusted p-values returned in input order,
computed by descending step-up with a running minimum (so the output is monotone
and capped at 1). Applied **across the domains in the legend**, i.e. the domains
actually tested, and the count is printed in the caption
("Benjamini–Hochberg FDR over N domains").

Stars: `***` q<0.001, `**` q<0.01, `*` q<0.05, unmarked = not significant
(`stars`, line 355).

### The parametric arm — the deviation-from-source machinery you need

This is directly relevant to your §2.7 framing, but note it runs **the opposite
way round** from how the outline describes it.

`parametricTest()` (line 270) runs **Welch's unequal-variance t-test** for two
groups and **one-way ANOVA** for k>2, on the same values in the same pass. Welch
is chosen explicitly because "it is what R's `t.test()` does by default, so it is
the most likely thing behind a published 't-test'" (line 236) — i.e. it is the
best available reconstruction of the source paper's test.

Tail probabilities come from a hand-rolled regularized incomplete beta
(`ibeta`, line 205) built on a Lanczos `lgamma` and a modified-Lentz continued
fraction, with `studentP` and `fUpperP` both reduced to it so there is one piece
of numerics rather than two.

**Crucially, the parametric result is reported but deliberately never wired to
the significance marks** (lines 156–162, 827–838). The stars always follow the
rank test. Both tests always run on every render, "so there is no opportunity to
pick the friendlier one after seeing it." The caption names any domain reaching
nominal p<0.05 under the parametric test and states in the same sentence that it
is *nominal, not corrected*, and that the marks follow the rank test. The full
per-domain table — n, median, p, q(BH), stars, parametric statistic, uncorrected
parametric p — goes to the CSV.

> **Reconcile this with your outline before writing.** Your §2.7 proposes: run
> the source's unadjusted test as the primary for the validation leg, report
> Wilcoxon+BH as sensitivity. The implementation does the reverse: rank+BH is
> primary and drives the marks, the uncorrected parametric test is the reported
> companion. Both numbers are produced and exported in every run, so either
> framing is defensible from the same output — but the *figures* the tool
> generates mark the rank test, and your prose must match the figures you paste
> in. The simplest fix is to describe the implementation as it is: two tests, one
> pass, pre-committed as to which one carries the marks.

### Whole-community test

`permanova()` (line 298) — PERMANOVA (Anderson 2001) on **Bray–Curtis**
distances between per-sample composition vectors. Pseudo-F is built from the
distance matrix directly (no ordination), p by permuting group labels with a
fixed seed (`19710407`), **9,999 permutations**, p = (hits+1)/(trials+1), and
R² = between/total = the share of total dissimilarity explained by the grouping.
It is computed on the same per-sample vectors over the same domains as the
per-domain tests, "so it cannot disagree with the per-domain results".

A limitation is flagged in the source itself and should be carried into the
methodology verbatim in substance (lines 286–288): there is **no betadisper
companion**. PERMANOVA cannot distinguish a shift in group centroid from a
difference in within-group spread, so a small p with visibly uneven scatter needs
that check before it is called a composition shift.

### Power floor — "a null FDR list is a reportable result", already implemented

`rankTestPlan()` (line 140) computes the smallest p the chosen test can return:
the exact floor 2/C(n, n₁) when enumerated, otherwise the larger of the sampling
resolution 1/(P+1) and the true combinatorial floor k!/space. `floorQ` is that
minimum floor multiplied by the number of tested domains — the smallest q a
single perfectly-separated domain could possibly reach.

When `floorQ > 0.05`, the figure prints, in warning colour:

> *"Underpowered: at these group sizes a perfectly separated domain reaches only
> q = X, so unmarked means undetectable, not absent. Needs about 6 samples per
> condition."*

The rationale (lines 758–762): "the caption has to say which, or a null result
gets read as a negative finding." The test file pins the arithmetic: 3v3 and 5v5
over 7 domains are unreachable at 0.05; **6v6 over 7 domains is reachable**; a
4×6 design is ample (`tests/test_etol_pie_stats.js:54-62`). That "about 6 samples
per condition" figure is a concrete, defensible power statement you can cite.

### What the tool refuses to compute, and why

Both refusals are explicit in code (lines 694–699) and are worth stating as
methodological choices rather than omissions:

1. **No test on raw matched reads.** Raw counts are not depth-normalized, so a
   "difference" can simply be a difference in library size. The UI returns an
   explanatory note instead of a p-value when the raw-count axis is selected.
2. **No chi-square / Fisher on the stack.** Reads within a sample are not
   independent draws, so that test returns p ≈ 0 for any two libraries and means
   nothing.

### Inclusion rules and units of replication

- The unit of replication is the **sample, not the read** — stated in the figure
  caption itself ("Replicates are samples, not reads").
- A sample that cannot be normalized (no control-probe reads, hence
  `host_cells = 0`) yields `null` and **drops out of its group** rather than
  counting as a zero it never measured (lines 716–730).
- A domain is tested only when **every group has ≥3 normalizable samples**;
  otherwise the figure reports "No significance test: needs at least 3
  normalizable samples per condition."
- At least two design-matrix conditions are required; otherwise "No significance
  test: needs at least two design-matrix conditions."
- Pooling a condition into one bar sums reads **and** host cells
  (depth-weighted), rather than averaging per-sample ratios, so a deeper library
  carries proportional weight — the same aggregation the pooled pie uses. The
  source notes that if per-sample spread matters more than the group total, the
  right move is to plot samples rather than add error bars.

> **Two accuracy notes for the write-up.**
> 1. `static/etol_pie.js:700` still carries the stale comment
>    `// ponytail: two groups only; add Kruskal-Wallis if a 3-condition design
>    shows up.` Kruskal–Wallis *was* added (commit `48327b7`, "Reworking
>    statistical analysis to account for more than 2 sample groups"). The comment
>    is wrong; don't quote it. Worth deleting in a tidy-up commit.
> 2. There is **no clustering support** — see the §2.6 gap box. Any claim in §2.7
>    about handling the EBB repeated-measures structure describes analysis done
>    outside the tool.

---

## 2.8 Development-time clinician input

> **Gap — no repository evidence.** Nothing in the tracked codebase records
> clinician or supervisor interaction: no acknowledgements file, no design notes,
> no issue trail. Grep for "clinician input", "supervisor", "co-supervisor"
> returns nothing.

The only record is in `.claude/memory/`, which is **git-ignored**
(`.gitignore:24`) and therefore local to your machine, not part of the
repository. It documents a dated exchange (2026-06-29) with two advisors who
replied on a set of development questions, covering: making re-probing off by
default while keeping it switchable; a proposed manual check that a zero
new-reads result is not a coding bug (BLAST a key contig back at the patient
database and `comm -23` against that taxon's net reads — empty means correct,
non-empty while the tool reports 0 means a bug; `scripts/run_manual_blast.sh`
exists for this); making the human component of the validation database
comprehensive since the filter can only reject what is in the database; emitting
a per-call provenance line and counting dropped human contigs rather than
dropping silently; and the position that large-subunit confirmation should be
targeted at headline taxa only rather than applied blanket.

That is genuinely **troubleshooting feedback during development** — which is
exactly how your outline says to characterise it. Describe it as collaboration
on implementation decisions. It is not clinician validation of outputs, and the
codebase contains no artefact that could support the stronger claim. Two of those
suggestions are visibly reflected in the shipped code (re-probing defaults off;
the human-filter note surfacing at batch level via
`summarize_human_filter_warnings`, `app.py:274`), which you can cite as evidence
the feedback was acted on.

---

## 2.9 GenAI and coding-agent disclosure

> **Gap — nothing in the repository.** There is no disclosure statement, no
> AI-contribution note, no `AI_USE.md`, and no tool attribution in commit
> trailers. The `.claude/` directory (agent configuration and memory) is
> git-ignored and ships with neither the repository nor the executable.

Factual material available for the disclosure, from git metadata:

- Two author identities appear in the history: `connor-aylesworth2
  <s2837739@ed.ac.uk>` and `ConnorAylesworth <connor.aylesworth@maine.edu>`.
  Worth a one-line note that both are you, so a reader counting contributors is
  not misled.
- Commit messages are descriptive of intent ("Reworking statistical analysis to
  account for more than 2 sample groups"; "Troubleshooting per-domain
  significance on stacked bars, PERMANOVA for whole community testing, adding
  parametric tests (Welch's t), & fixing how they're all reported"), which
  supports a claim of author-directed development.
- Several source files carry `ponytail:` comments marking deliberate
  simplifications with their known ceiling and upgrade path — e.g. the global
  human-BLAST lock (`human_filter.py:49`), one re-probe contig
  (`contig_id.py:56`), 100k re-probe target cap (`contig_id.py:63`), depth-
  weighted pooling (`static/etol_pie.js:634`), and the missing betadisper
  companion (`static/etol_pie.js:286`). These are an auditable record of where
  simplifications were made knowingly rather than by oversight, which is useful
  material for a disclosure that wants to demonstrate author understanding of the
  code.

The disclosure text itself must be authored and cross-referenced to your cover
sheet; nothing in the codebase can be cited for it.

---

## Consolidated action list before writing

**Fix in the repository:**
1. Restore `scripts/build_etol_v_validation_db.sh` from
   `git show 0db7356^:scripts/build_etol_v_validation_db.sh` — it is currently 0 bytes.
2. Re-run `python tests/test_spike_in_control.py` from a clean tree at your
   submission commit and regenerate `docs/appendix_spikein.txt` (also fixes the
   UTF-16 encoding).
3. Delete the stale "two groups only" comment at `static/etol_pie.js:700`.
4. Optionally add a Node step to CI so `tests/test_etol_pie_stats.js` is covered.

**Author from scratch (no codebase source exists):**
5. The pre-registered acceptance criteria table (§2.5.1) — write and date it
   before running anything.
6. Kohen 2014 dataset description in full (§2.3).
7. EBB clinical covariates: age, sex, PMI, Braak, Thal, APOE (§2.3 table).
8. The entire application protocol (§2.6) — and state plainly that these analyses
   run downstream of COBLAST+'s exported count matrices, not inside the tool.
9. Clinician-input paragraph (§2.8) and GenAI disclosure (§2.9).

**Correct before writing:**
10. §2.5.3 — COBLAST+ *does* generate contig-derived probes and re-BLAST. The
    real limitation is the absent 23S/28S disambiguation pass. Replacement
    wording is given in §2.5.3 above.
11. §2.7 — the implementation makes rank+BH primary and the parametric test the
    reported-but-unmarked companion, which is the reverse of the outline's
    framing. Match the prose to the figures.
12. §2.6.3 / §2.7 — no clustering, paired, or mixed-model support exists in the
    tool. Any repeated-measures handling is downstream analysis and must be
    described as such.
