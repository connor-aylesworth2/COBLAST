---
name: etol-lathe-alignment
description: "Progress aligning COBLAST+ eToL preset to Lathe et al. 2022; what's done and what remains"
metadata: 
  node_type: memory
  type: project
  originSessionId: e9745756-996a-4ea2-928b-cfdda637d8a6
---

Goal: make COBLAST+'s eToL batch analysis match the method in Hu, Haas & Lathe,
*BMC Microbiology* 2022;22:317 (the "electronic tree of life" net), so results
can be compared against standard CLI BLAST+ runs. The eToL preset was originally
built as a copy of the APOE exact-match genotyper, which was wrong for eToL.

Done (as of 2026-06-15):
1. **First search = the net.** Dropped the 100%-identity/100%-coverage exact-match
   enforcement; the eToL panels now run BLAST **default megablast** with no
   identity or coverage filter (default e-value is the only gate), `blastn-short`
   for the 1 ambiguous probe (`F3_Gpolymorpha_18S_7`). Only `max_target_seqs` is
   lifted so deep DBs aren't truncated. Code: `etol_net_probe` path in
   `blast_runner.py`; `filter_net_probe_hits` (panel-restriction only) in `app.py`.
2. **Secondary human filter = bitscore.** Changed from "100% read query coverage"
   to "best human alignment > 150 bits" (the paper's brain/liver/skin cutoff),
   bitscore as sole criterion, permissive e-value. `HUMAN_BITSCORE_THRESHOLD` in
   `human_filter.py`.
3. **Cross-probe read de-duplication.** After the human filter, each matched read
   is allocated to the single best probe (ranked bitscore → identity → coverage).
   `deduplicate_reads_to_best_probe` in `app.py`.
4. **Host-cell normalization.** Control probes (PGK1, hNSE) are appended to the
   microbial search (Full/Quick = 1021 probes). Control hits are split out
   **before** the human filter (control reads are human by design — must NOT be
   filtered, or the denominator zeroes). host cells = mean(per-gene control
   readcount)/50; reports reads-per-host-cell in UI + CSV/TSV. Functions in
   `etol_summary.py` (`etol_search_fasta`, `compute_host_cells`, etc.).
5. **E-value < 0.01 net gate + control dedup** (added 2026-06-16, matching the
   original `Abundance_ToL.py`/`Abundance_count.py`). `filter_net_probe_hits`
   (app.py) now drops hits with E ≥ 0.01 or unparseable E, for both microbial and
   control probes — BLAST default e-value is **10**, so "default e-value is the
   only gate" was wrong/too permissive. Controls are now de-duplicated to their
   single best control probe via new `count_control_reads` (reuses
   `deduplicate_reads_to_best_probe`) before counting, so a read hitting several
   redundant control probes no longer inflates the normalization denominator.
   Dedup ordering kept as human-filter→dedup (paper's text; the original code's
   comment intended this but ran dedup first).
- 84 tests pass. Verified at unit/template/export level only — NOT yet run on a
  real SRA + human genome DB.

Contig assembly DONE (as of 2026-06-22): CAP3 backend behind an `Assembler`
protocol (`assembler.py`, `Cap3Assembler`, CAP3 defaults), reads grouped
**per taxon** via `group_read_ids_by_taxon` (etol_summary.py), assembled in a
batch-wide thread pool, stored per row + downloadable as multi-FASTA
(`/batch-results/<id>/etol-contigs.fasta`, app.py). Note: paper's fast path
groups per phylogenetic group A-H; per-taxon is a finer, defensible deviation.

Contig species-ID + confirmed abundance DONE (as of 2026-06-23): new
`contig_id.py` (`identify_contigs`, `name_contigs`, `confirm_contig_reads`,
pure parsers `_best_homolog_per_query`/`_confirmed_reads_from_tabular`).
(a) **Species naming** = one batched `blastn -task megablast` of every contig
(synthetic ids `c0,c1,..` map back to (taxon,contig) since CAP3 restarts
`Contig1` numbering per taxon) vs a registered **reference** rRNA DB; keeps the
best-bitscore hit's `stitle` as the closest homolog. (b) **Confirmed abundance**
= per taxon, `blastn -subject` of that taxon's reads vs its own contigs, count
distinct reads at **≥99% identity** (`DEFAULT_CONFIRM_IDENTITY_PCT`, the user's
"near-100%" choice). Annotates each contig dict in place + returns per-taxon
`{closest_homolog, homolog_pident, confirmed_reads}` consumed by
`build_etol_probe_summary`/`etol_species_count_rows`. UI: new "Contig species
identification" checkbox + reference-DB picker (mirrors the human-filter picker),
gated on assembly being on; results add a "Closest homolog (contig)" + "Confirmed
reads" column to the species + contig tables and to the species-summary CSV/TSV;
contig FASTA header gains `|confirmed=|homolog=`. New `reference` registry
category. Stays LOCAL — reference DB is the user's SILVA SSU+LSU NR99 build
`ToL_rRNA` (see local-DB note below). NOT yet run on a real SRA + ToL_rRNA.

Re-probing DONE (as of 2026-06-24): `reprobe_and_reassemble` in `contig_id.py`
(+ pure parser `_reprobe_reads_by_query`, blastn helper `_reprobe_hits`). One
round (Box 3): each taxon's top `REPROBE_TOP_CONTIGS=2` contigs by read support
are batched into ONE `blastn megablast` vs the SAME patient DB already in the
batch (no extra DB picker), matches gated at E<`REPROBE_EVALUE=0.01` (the net's
gate), `-max_target_seqs 100000`; genuinely new reads (not already in the taxon)
are extracted, optionally human-filtered (reuses the run's human DB +
`find_human_read_ids`, as the paper filters re-probe matches), then each taxon is
re-assembled from original+new reads. Mutates `contigs_by_species`/
`reads_by_taxon`/`all_reads` in place so the later naming+confirm steps use the
extended contigs. Route order: assemble → reprobe → identify. New
`reprobe_contigs` checkbox (gated on assembly, microbial). Serial re-assembly
(ponytail-noted ceiling). 9 tests in `tests/test_contig_id.py`; 100 pass total.
NOT yet run on a real SRA.

Reference DB (decided 2026-06-23): contig species-ID searches a LOCAL curated
rRNA DB, not nt — SILVA **SSU+LSU NR99** (`tax_silva`, U→T converted, one
`makeblastdb -parse_seqids`), built by the user as `ToL_rRNA` (~0.5 GB; orders of
magnitude under nt). SSU alone spans all domains A–H (16S pro + 18S euk); LSU
adds 23S/28S for the validation step. Registered under the new `reference`
category. Homologs come back as SILVA taxonomy strings (no NCBI taxids).

Remaining to match the paper (largest first):
- The entire contig-based species-ID arc (closest-homolog naming +
  contig-confirmed abundance + re-probing) is now IMPLEMENTED and local — see the
  three DONE blocks above. Only multi-round/iterative re-probing is out: the
  implementation does 1 round (the user's chosen scope); the paper hints at
  iterating until no new reads. Add a capped loop only if a real run shows one
  round leaves contigs short.
- **23S/28S (+mtDNA) confirmation** step to disambiguate redundant probes
  (Bonferroni concern with >1000 probes); semi-manual in the paper, no hook yet.
- **Read-length-adaptive human cutoff** (paper varied >150/>160/>126/>100 by mean
  read size; COBLAST fixes 150 for brain).
- **Heatmap visualization** (paper uses Morpheus/pheatmap; COBLAST exports CSV).
- Viral "stripping" method is a separate paper method, not implemented (out of
  current scope).

Human-filter DB choice (decided 2026-06-16): the secondary human filter must use
the **Homo sapiens subset of NCBI nt** (`blastdbcmd -db nt -taxids 9606 -outfmt
"%f"` → `makeblastdb`), NOT a bare GRCh38 genome assembly. The paper ran
`blastn -db nt -taxids 9606`, and for brain RNA-seq the nt curated mRNAs catch
exon-junction-spanning human reads the genome assembly misses. User has built
`nt_human_9606` and registers it as category `human`. COBLAST's human_filter.py
is DB-agnostic (no code change). Full-nt + `-taxids` plumbing was offered but
declined for now.

Design notes worth keeping: eToL preset still forces blastn+tabular; APOE preset
is untouched (still 100/100 exact). CSV export column headers kept stable
(`Exact hits`, `Total exact probe hits`) to avoid breaking downstream plots.
See [[coblast-architecture]] for module layout.
