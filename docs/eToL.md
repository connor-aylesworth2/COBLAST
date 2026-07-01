# COBLAST+ eToL Probe Presets

The eToL probe presets are an **optional** analysis built on COBLAST+'s general
batch-BLAST workflow (see the "SRA Pilot and Batch Workflow" section of the main
[README](../README.md)). They are not required for ordinary BLAST use. This
document covers the eToL presets in full: the net, the two-pass task split,
cross-probe de-duplication, host-cell normalization, the secondary human filter,
and the export formats.

## eToL probe presets (the net)

The batch page includes **eToL probe presets** built for the
electronic Tree of Life (eToL) workflow described in Hu, Haas & Lathe,
*BMC Microbiology* 2022;22:317. Each preset BLASTNs a stored probe panel against
the selected nucleotide databases. Unlike the APOE genotyper (which saves only
100% identity / 100% coverage hits), the eToL panels reproduce the paper's
permissive *net*: BLAST's **default megablast** scoring with **no identity or
coverage filter**, gated on **E-value < 0.01** (the paper's net cutoff, applied
in `Abundance_ToL.py`), so partial and mismatched rRNA matches are retained for
the secondary human filter and cross-probe de-duplication to adjudicate while
statistically insignificant matches (down to BLAST's default e-value of 10) are
dropped. The probe FASTA is supplied
automatically, so the query box is left empty (read-only) when a preset is on.

The only override the net applies is lifting the `max_target_seqs` cap (see the
`EXACT_MATCH_MAX_TARGET_SEQS` constant in `blast_runner.py`), so a probe that
matches many reads in a deep patient database is counted in full rather than
truncated at BLAST's default.

To keep whole-SRA runs fast, the entire eToL panel is searched with `megablast`
(much faster on large read databases). `megablast` needs a 28-base unambiguous
window to seed, so a probe whose ambiguous bases leave no such window finds
nothing and is silently dropped. In the bundled microbial panel that is a single
probe (`F3_Gpolymorpha_18S_7`).

**Cross-probe de-duplication.** Because the probe collection is partly redundant
(rRNA is conserved, so ~38% of probes share sequence with at least one other),
the same read can be recovered by several probes. After the human filter,
COBLAST allocates each matched read to the single probe with the highest
similarity (ranked by bitscore, then identity, then coverage), exactly as the
paper specifies, so a read is counted once rather than inflating several probes.

**Host-cell normalization.** The microbial presets are searched together with the
housekeeping control probes (PGK1, hNSE) in the same run. The control reads are
counted separately (never human-filtered — they are human by design); like the
microbial net they pass the same E-value < 0.01 gate and are de-duplicated so
each read is allocated to its single best control probe before counting
(`Abundance_count.py` section 2), so a read recovered by several redundant
control probes does not inflate the normalization denominator. They are used to
estimate host abundance: the host-cell count is the mean per-gene control
readcount divided by ~50 transcripts per cell (`HOST_TRANSCRIPTS_PER_CELL`).
Microbial counts are then reported both raw and as **reads per host cell**
(`raw / host cells`), normalizing for how much host material each library
represents, as in the paper. When no control reads are found, normalization is
reported as `n/a`.

There are two eToL presets, plus the APOE preset; **only one preset can be
active at a time** (selecting one clears the others). Both eToL presets append
the PGK1/hNSE control probes to their search for host-cell normalization:

- **eToL Full** — the full microbial panel, `data/eToL_probes.fasta` (1,017
  64-mer probes, 120 species across Archaea, Bacteria, Chloroplastida,
  Amoebozoa, basal Eukaryota, Fungi, and Holozoa/Metazoa).
- **eToL Quick** — one probe per species (the first probe of each of the 120
  species), a slim 120-probe panel for fast test runs.

This is intended for the patient-sample use case: register a patient's brain (or
other tissue) RNA-seq reads as a local nucleotide database, select it (or several
patients) in the batch picker, and run a preset to count probe matches (the net).

eToL batch results include an `eToL Probe Summary` section with one block per
selected sample/database. Each block reports the total probe-matched reads, how
many probes were detected, how many species/taxa were detected, and a table of
the detected species (grouped by domain and eToL class code, sorted by
matched-read count). Species are shown by species label only (the class prefix and rRNA-unit
suffix are stripped, e.g. `B0_Tmaritima_16S` is displayed as `Tmaritima`). Class
codes map to domains per the eToL paper: A Archaea; B Bacteria; C Chloroplastida;
D Amoebozoa; E0 basal Eukaryota; F Fungi; H Holozoa/Metazoa.

Two count exports are offered alongside the raw hit table:

- **eToL Probe Counts** (CSV/TSV) — one row per probe per sample for every probe
  in the active panel (including zeros): `Sample/Database, Probe, Species/Taxon,
  Class, Domain, Exact hits`. This is the full count matrix for species plots.
- **eToL Species Summary** (CSV/TSV) — one row per species/taxon per sample
  (including zeros): `Sample/Database, Domain, Class, Species/Taxon, Probes in
  panel, Probes detected, Total exact probe hits`.

Sample labels follow the same SRA-accession rule as the APOE preset. Because the
eToL panel contains far more than the previous 100-record limit, the maximum
number of FASTA query records per run is 1,500 (`MAX_FASTA_RECORDS` in
`blast_runner.py`), which accommodates the full panel with headroom.

## Heatmap condition labels (design matrix)

The eToL Result Heatmap annotates each sample column with a coloured **condition**
swatch (AD vs control, etc.). By default the label is inferred from the database
name, which is unreliable — auto-generated SRA database names (`SRA <acc> reads`)
have no diagnosis in them, and substring matching mislabels samples.

To set the labels explicitly, upload a **design matrix** in the batch form
(*Design matrix (condition labels)*; eToL presets only). When provided it is
authoritative and the name-based guess is not used.

Format (strict): CSV or TSV, UTF-8 (a leading BOM is tolerated), with a header
row containing two columns, case-insensitive and in any order:

- `sample` — an SRA accession (SRR/SRX/ERR/…) **or** a database display name.
- `condition` — the free-text label to show (e.g. `AD`, `CONTROL`, `AD/LBD`).

```
sample,condition
SRR21676099,AD
SRR21676105,CONTROL
SRR21676101,AD/LBD
SRR21676126,AD/VaD
```

Matching is flexible: a row binds to a column by SRA accession first, then by the
exact database display name. One row per sample (a duplicate `sample` is an error).
Extra columns are ignored (room for a future multi-factor mode). A malformed file
(missing header/column, duplicate sample, no data rows) is rejected on the form
*before* the BLAST run. Samples with no matching row render as a neutral
"unlabeled" swatch and are listed in a warning under the heatmap; a blank
`condition` is allowed and also renders unlabeled. A starter file is downloadable
from the form (`/design-matrix-template.csv`). Parsing lives in `design_matrix.py`;
`etol_summary.build_etol_matrix(..., condition_index=...)` applies it.

## Secondary human filter

The microbial eToL presets (eToL Full and eToL Quick) offer an optional
**secondary human filter** that removes matched patient reads that are actually
human-derived — the second-round host filtering required by Hu, Haas & Lathe
2022 because sequence similarity is non-transitive (a read can be human yet still
match a probe that itself has no human match). When enabled (checkbox on the
batch page, plus a human-genome database selector), COBLAST takes the net probe
hits for each patient database, recovers the full matched reads, BLASTs them
(`megablast`) against the selected human genome database, and drops every read
whose best human alignment scores **above 150 bits**; the bitscore is the sole
criterion, with no coverage or E-value cutoff (the E-value is set permissively so
it never filters). The >150-bit cutoff is the value the paper applied to brain
(and liver/skin) datasets, so it is COBLAST's default for brain samples — note
the paper *adjusted* this cutoff to each dataset's mean read length (>160 MSBB,
>126 Rockefeller, >100 Miami), so libraries with very different read lengths may
warrant a different threshold (`HUMAN_BITSCORE_THRESHOLD` in `human_filter.py`).
The results page reports how many hits were removed per sample, and all
summaries/exports reflect the filtered hit list.

Matched reads are recovered by their `sseqid` (which equals the read's FASTA
record id): first with `blastdbcmd` (when the patient database was built with
`-parse_seqids`), otherwise by scanning the database's stored source FASTA. If a
read cannot be recovered it is kept unfiltered (never dropped on a guess), and a
note is shown. The SRA workbench automatically associates a discovered BLAST
database with a same-named FASTA in that SRA project. Re-register an older SRA
database if its registry entry predates this association.

To set up the human genome database, build it once with `makeblastdb` straight
from the NCBI genome FASTA — `.fna` files are already FASTA, so no conversion is
needed — then register it on the Databases page with category `human`:

```text
makeblastdb -in GCF_000001405.40_GRCh38.p14_genomic.fna -dbtype nucl ^
  -title "Human GRCh38.p14" -out human_GRCh38
```

The filter is offered only for the microbial presets; the APOE panel is human by
design, so human-read filtering does not apply to it.

## Contig re-probing (optional second pass)

Re-probing (Hu, Haas & Lathe 2022, Box 3) takes each taxon's most-abundant
assembled contig, BLASTs it back against the *same* patient database as a fresh
probe, pulls any reads the original 64-mer net missed, human-filters those new
reads, and re-assembles the taxon from its original + new reads. Because a contig
is built from the sample's own reads, it is a more sensitive, sample-specific
probe than the fixed reference panel — it can reach reads from a divergent strain
or from gene regions that fall between the net's probes.

Requires **contig assembly** to be enabled, and applies only to the microbial
eToL presets.

Contig assembly (and therefore re-probing) uses the **CAP3** assembler, which
COBLAST+ does not ship for licensing reasons. Install
[Unipro UGENE](https://ugene.net/download-all.html) in its default location — it
bundles CAP3 and COBLAST+ auto-detects it (`…\Unipro UGENE\tools\cap3`) — or set
`CAP3_BIN` to a folder holding a CAP3 binary. Without CAP3, assembly-dependent
steps are skipped and the run reports it.

### When to leave it OFF (the default)

For routine microbiome profiling, **don't enable re-probing.** The net already
runs at a permissive E-value (< 0.01); on most libraries it has already captured
everything a contig would, so re-probing recovers **0 new reads** and changes
nothing while still costing a full-library BLAST per taxon. It is also off by
default for cohort work: re-probing rewrites species names and confirmed
abundance in place, so on/off results are not directly comparable — pick one
setting and apply it uniformly across every sample in a study.

### Characterise it once per data type, not once per sample

Whether re-probing finds anything is a property of the *data regime* — the probe
panel against a library's read length, rRNA-depletion protocol, and how far its
organisms diverge from the reference panel — not an independent roll per sample.
Validate it at that granularity: on a new kind of data, run a representative
subset **both** ways (net-only and net+re-probe) and compare. If
`reprobe_new_reads` is 0 across the subset, the two results are identical and you
can leave re-probing off for the rest of that regime, re-checking only when the
data type changes (new tissue, protocol, or read length). Running both and seeing
them agree is what reassures you; turning re-probing on "to be sure" assumes the
answer instead of demonstrating it — and is not automatically more accurate (see
below).

### When a run does recover reads

Read the **`reprobe_new_reads`** count in the results:

- **0** — the net had already saturated capture for that sample. Re-probing
  confirmed there was nothing to add; trust the first-pass result.
- **> 0** — re-probing extended one or more contigs. Inspect these: a longer,
  more conserved contig can also over-recruit reads from a *closely related*
  taxon, so confirm the recovered reads reflect a genuinely missed strain rather
  than a neighbour bleeding into the call.

A free gut-check needs no extra run: re-probing can only gain reads where a contig
extends *beyond* the probe footprint, so if a sample's contigs are barely longer
than the 64-mer probes that seeded them, re-probing has no territory to work in
and will be 0 by construction.

### Caveats

- Re-probing is a sensitivity check, not a substitute for the net, and it is
  lightly validated — treat a non-zero result as a prompt to investigate, not a
  finished answer.
- Enabling it currently re-runs the whole pipeline (net search included), not
  just the re-probe step.
