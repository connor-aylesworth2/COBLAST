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
coverage filter** (BLAST's default e-value is the only gate), so partial and
mismatched rRNA matches are all retained for the secondary human filter and
cross-probe de-duplication to adjudicate. The probe FASTA is supplied
automatically, so the query box is left empty (read-only) when a preset is on.

The only override the net applies is lifting the `max_target_seqs` cap (see the
`EXACT_MATCH_MAX_TARGET_SEQS` constant in `blast_runner.py`), so a probe that
matches many reads in a deep patient database is counted in full rather than
truncated at BLAST's default.

To keep whole-SRA runs fast, the eToL panel is searched in two passes (matching
the paper's two-task approach): probes that have a 28-base unambiguous window run
with `megablast` (much faster on large read databases), and the few whose
ambiguous bases leave no such window fall back to `blastn-short`. In the bundled
microbial panel that is a single probe (`F3_Gpolymorpha_18S_7`); the partition is
computed from each probe's sequence, so it self-adjusts if the panel changes.

**Cross-probe de-duplication.** Because the probe collection is partly redundant
(rRNA is conserved, so ~38% of probes share sequence with at least one other),
the same read can be recovered by several probes. After the human filter,
COBLAST allocates each matched read to the single probe with the highest
similarity (ranked by bitscore, then identity, then coverage), exactly as the
paper specifies, so a read is counted once rather than inflating several probes.

**Host-cell normalization.** The microbial presets are searched together with the
housekeeping control probes (PGK1, hNSE) in the same run. The control reads are
counted separately (never human-filtered — they are human by design) and used to
estimate host abundance: the host-cell count is the mean per-gene control
readcount divided by ~50 transcripts per cell (`HOST_TRANSCRIPTS_PER_CELL`).
Microbial counts are then reported both raw and as **reads per host cell**
(`raw / host cells`), normalizing for how much host material each library
represents, as in the paper. When no control reads are found, normalization is
reported as `n/a`. The standalone **eToL Control** preset still runs the control
probes on their own for QC.

There are three eToL presets, plus the APOE preset; **only one preset can be
active at a time** (selecting one clears the others):

- **eToL Full** — the full microbial panel, `data/eToL_probes.fasta` (1,017
  64-mer probes, 120 species across Archaea, Bacteria, Chloroplastida,
  Amoebozoa, basal Eukaryota, Fungi, and Holozoa/Metazoa).
- **eToL Control** — the human housekeeping control probes only,
  `data/eToL_control_probes.fasta` (PGK1, hNSE). Same workflow as eToL Full but
  uses human sequences as a control.
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

The filter is offered only for the microbial presets; the APOE and eToL Control
panels are human by design, so human-read filtering does not apply to them.
