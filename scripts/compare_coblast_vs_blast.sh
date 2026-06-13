#!/usr/bin/env bash
#
# compare_coblast_vs_blast.sh
#
# Reproduce, with raw NCBI BLAST+ commands and standard Unix tools, exactly the
# analysis COBLAST+ runs for the "eToL Full" preset WITH the secondary human
# filter, across all of your SRA read sets. Then compare the per-probe counts
# this script produces against the counts COBLAST+ reports, to confirm the
# wrapper matches command-line BLAST+.
#
# Per SRA, identical to COBLAST+:
#   1. eToL exact-probe search, split by megablast seed eligibility:
#        - probes with a >=28-base unambiguous window  -> blastn -task megablast
#        - probes that cannot seed (1 in the bundled panel: F3_Gpolymorpha_18S_7)
#                                                       -> blastn -task blastn-short
#      both with -perc_identity 100 -qcov_hsp_perc 100 -max_target_seqs 5000000,
#      outfmt "6 qseqid sseqid stitle pident length qcovs evalue bitscore".
#      Keep only rows with pident==100 AND qcovs==100 (the exact-match filter).
#   2. Secondary human filter: BLAST the matched reads against the human genome
#      (blastn -task megablast -evalue 1e-6 -qcov_hsp_perc 100
#      -max_target_seqs 1) and drop eToL hits whose read has a full-query-
#      coverage human-genome HSP.
#
# makeblastdb is run WITHOUT -parse_seqids and reads are recovered from the
# source FASTA, exactly as COBLAST+ does.
#
# ASSUMPTIONS (set the CONFIG block below to match your machine):
#   * Each SRA has its own directory directly under DATA_DIR holding the combined
#     read FASTA. SRA_GLOB matches it while SRA_EXCLUDE_GLOB skips any paired-end
#     mate files (so SRRxxxx.fasta is kept and SRRxxxx_1.fasta / _2.fasta are
#     skipped). Each sample is labelled by its folder under DATA_DIR.
#   * The human genome is either a FASTA file or an existing BLAST DB prefix.
#   * If your SRAs are still .sra/FASTQ, convert them to FASTA first
#     (e.g. fastq-dump --fasta), or ask for a version that adds that step.
#
set -euo pipefail

# ============================== CONFIG ==============================
# Path to your COBLAST+ checkout (must contain data/eToL_probes.fasta).
COBLAST_DIR="/path/to/COBLAST-"

# Directory holding your SRA read sets, one subdirectory per SRA.
DATA_DIR="/path/to/sra_data"

# Glob matching each SRA's read FASTA (searched recursively under DATA_DIR).
# Pick a pattern for the ONE combined read FASTA per SRA.
SRA_GLOB="SRR*.fasta"

# Filenames matching this glob are skipped, so paired-end mate files
# (e.g. SRRxxxx_1.fasta / SRRxxxx_2.fasta) are excluded and only the combined
# SRRxxxx.fasta is searched. Set to "" to skip nothing.
SRA_EXCLUDE_GLOB="*_*.fasta"

# Human genome: an absolute path to a FASTA file OR to an existing BLAST
# nucleotide DB prefix. If it is a FASTA with no DB yet, the DB is built once.
HUMAN_GENOME="/path/to/human_genome.fasta"

# Where results go. Keep this OUTSIDE DATA_DIR so generated files are not
# rescanned as if they were SRAs (a relative path lands beside where you run it).
OUT_DIR="./coblast_vs_blast_results"

# CPU threads per BLAST call (threads do not change the results, only the speed).
THREADS=8
# ====================================================================

PROBES="$COBLAST_DIR/data/eToL_probes.fasta"
OUTFMT="6 qseqid sseqid stitle pident length qcovs evalue bitscore"
HUMAN_EVALUE="1e-6"
HUMAN_QCOV="100"
MAX_TARGET_SEQS="5000000"

# --- sanity checks ---
for tool in blastn makeblastdb python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' is not on PATH." >&2; exit 1; }
done
[[ -f "$PROBES" ]] || { echo "ERROR: eToL panel not found: $PROBES" >&2; exit 1; }
[[ -d "$DATA_DIR" ]] || { echo "ERROR: DATA_DIR not found: $DATA_DIR" >&2; exit 1; }

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"        # absolute, so the find-prune below is reliable
WORK="$OUT_DIR/_work"
mkdir -p "$WORK"

# --- 1. Split the eToL panel into megablast-safe and blastn-short subsets ---
#        (28-base unambiguous-window rule, identical to COBLAST has_megablast_seed)
MB_PROBES="$WORK/etol_megablast.fasta"
SHORT_PROBES="$WORK/etol_blastn_short.fasta"
python3 - "$PROBES" "$MB_PROBES" "$SHORT_PROBES" <<'PY'
import sys
src, mb_out, short_out = sys.argv[1:4]

def longest_clean(seq):
    best = cur = 0
    for base in seq.upper():
        if base in "ACGT":
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best

recs, header, seq = [], None, []
for line in open(src):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if header is not None:
            recs.append((header, "".join(seq)))
        header, seq = line, []
    elif line.strip():
        seq.append(line.strip())
if header is not None:
    recs.append((header, "".join(seq)))

with open(mb_out, "w") as mb, open(short_out, "w") as sh:
    n_mb = n_sh = 0
    for header, seq in recs:
        if longest_clean(seq) >= 28:
            mb.write(header + "\n" + seq + "\n"); n_mb += 1
        else:
            sh.write(header + "\n" + seq + "\n"); n_sh += 1
print(f"  eToL panel split: {n_mb} megablast probes, {n_sh} blastn-short probes",
      file=sys.stderr)
PY

# --- 2. Build the human-genome DB if HUMAN_GENOME is a FASTA without a DB ---
if [[ -f "${HUMAN_GENOME}.nsq" || -f "${HUMAN_GENOME}.00.nsq" || -f "${HUMAN_GENOME}.nal" ]]; then
  HUMAN_DB="$HUMAN_GENOME"                          # already a BLAST DB prefix
elif [[ -f "$HUMAN_GENOME" ]]; then
  HUMAN_DB="$WORK/human_db"
  if [[ ! -f "${HUMAN_DB}.nsq" && ! -f "${HUMAN_DB}.00.nsq" ]]; then
    echo "Building human-genome BLAST DB (one time)..." >&2
    makeblastdb -in "$HUMAN_GENOME" -dbtype nucl -out "$HUMAN_DB" >/dev/null
  fi
else
  echo "ERROR: human genome not found as FASTA or DB prefix: $HUMAN_GENOME" >&2
  exit 1
fi

# --- 3. Discover the per-SRA read FASTAs (recursively; excluding outputs/human) ---
shopt -s nullglob
find_expr=( "$DATA_DIR" -path "$OUT_DIR" -prune -o -type f -name "$SRA_GLOB" )
[[ -n "$SRA_EXCLUDE_GLOB" ]] && find_expr+=( ! -name "$SRA_EXCLUDE_GLOB" )
find_expr+=( -print )
mapfile -t SRA_FILES < <( find "${find_expr[@]}" | sort )
KEEP=()
for f in "${SRA_FILES[@]:-}"; do
  [[ "$f" == "$HUMAN_GENOME" ]] && continue
  KEEP+=("$f")
done
SRA_FILES=("${KEEP[@]:-}")
(( ${#SRA_FILES[@]} > 0 )) || { echo "ERROR: no SRA FASTAs matched '$SRA_GLOB' under $DATA_DIR" >&2; exit 1; }
echo "Found ${#SRA_FILES[@]} SRA FASTA file(s)." >&2

# --- 4. Per-SRA pipeline ---
SUMMARY="$OUT_DIR/summary_probe_counts.tsv"
printf 'sample\tprobe\texact_hits\n' > "$SUMMARY"

for SRA in "${SRA_FILES[@]}"; do
  # Sample label = the SRA's folder directly under DATA_DIR.
  rel="${SRA#"$DATA_DIR"/}"
  name="${rel%%/*}"
  [[ "$name" == "$rel" ]] && name="$(basename "${SRA%.*}")"   # fallback: FASTA sits in DATA_DIR
  echo "=== $name ===" >&2
  sdir="$OUT_DIR/$name"; mkdir -p "$sdir"

  # 4a. Build a nucleotide DB (no -parse_seqids, exactly like COBLAST).
  db="$WORK/${name}_db"
  makeblastdb -in "$SRA" -dbtype nucl -out "$db" >/dev/null

  # 4b. eToL exact-probe search: megablast subset + blastn-short subset.
  blastn -task megablast    -query "$MB_PROBES"    -db "$db" -num_threads "$THREADS" \
         -perc_identity 100 -qcov_hsp_perc 100 -max_target_seqs "$MAX_TARGET_SEQS" \
         -outfmt "$OUTFMT" > "$sdir/etol_megablast.tsv"
  blastn -task blastn-short -query "$SHORT_PROBES" -db "$db" -num_threads "$THREADS" \
         -perc_identity 100 -qcov_hsp_perc 100 -max_target_seqs "$MAX_TARGET_SEQS" \
         -outfmt "$OUTFMT" > "$sdir/etol_blastn_short.tsv"

  # 4c. Merge + exact-match filter (pident==100 AND qcovs==100). Columns:
  #     1 qseqid  2 sseqid  3 stitle  4 pident  5 length  6 qcovs  7 evalue  8 bitscore
  cat "$sdir/etol_megablast.tsv" "$sdir/etol_blastn_short.tsv" \
    | awk -F'\t' '$4==100 && $6==100' > "$sdir/etol_exact_hits.tsv"

  # 4d. Recover the matched reads (by sseqid) from the SRA FASTA.
  cut -f2 "$sdir/etol_exact_hits.tsv" | sort -u > "$sdir/matched_read_ids.txt"
  python3 - "$sdir/matched_read_ids.txt" "$SRA" > "$sdir/matched_reads.fasta" <<'PY'
import sys
ids = {x.strip() for x in open(sys.argv[1]) if x.strip()}
keep = False
for line in open(sys.argv[2]):
    if line.startswith(">"):
        keep = line[1:].split()[0] in ids
    if keep:
        sys.stdout.write(line)
PY

  awk '/^>/{sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1]}' \
    "$sdir/matched_reads.fasta" | sort -u > "$sdir/recovered_read_ids.txt"
  comm -23 "$sdir/matched_read_ids.txt" "$sdir/recovered_read_ids.txt" \
    > "$sdir/unresolved_read_ids.txt"
  matched_reads=$(wc -l < "$sdir/matched_read_ids.txt")
  recovered_reads=$(wc -l < "$sdir/recovered_read_ids.txt")
  unresolved_reads=$(wc -l < "$sdir/unresolved_read_ids.txt")
  echo "  read recovery: $recovered_reads/$matched_reads recovered; $unresolved_reads unresolved" >&2

  # 4e. Human filter: BLAST matched reads vs human genome; collect read IDs that hit.
  if [[ -s "$sdir/matched_reads.fasta" ]]; then
    blastn -task megablast -query "$sdir/matched_reads.fasta" -db "$HUMAN_DB" \
           -num_threads "$THREADS" -evalue "$HUMAN_EVALUE" \
           -qcov_hsp_perc "$HUMAN_QCOV" -max_target_seqs 1 \
           -outfmt "6 qseqid sseqid pident length qcovhsp evalue bitscore" \
           > "$sdir/human_hits.tsv"
    awk -F'\t' '$5==100 {print $1}' "$sdir/human_hits.tsv" \
      | sort -u > "$sdir/human_read_ids.txt"
  else
    : > "$sdir/human_hits.tsv"
    : > "$sdir/human_read_ids.txt"
  fi
  human_reads=$(wc -l < "$sdir/human_read_ids.txt")
  echo "  human matches: $human_reads/$recovered_reads recovered read(s)" >&2

  # 4f. Drop eToL hits whose read hit the human genome.
  awk -F'\t' 'NR==FNR{h[$1]; next} !($2 in h)' \
      "$sdir/human_read_ids.txt" "$sdir/etol_exact_hits.tsv" \
      > "$sdir/etol_human_filtered.tsv"

  # 4g. Per-probe exact-hit counts after the human filter (what COBLAST reports).
  cut -f1 "$sdir/etol_human_filtered.tsv" | sort | uniq -c \
    | awk -v s="$name" '{print s"\t"$2"\t"$1}' >> "$SUMMARY"

  before=$(wc -l < "$sdir/etol_exact_hits.tsv")
  after=$(wc -l < "$sdir/etol_human_filtered.tsv")
  echo "  exact hits: $before; human-removed: $(( before - after )); kept: $after" >&2
done

echo "" >&2
echo "Done." >&2
echo "  Per-probe counts (long format) : $SUMMARY" >&2
echo "  Per-sample intermediates       : $OUT_DIR/<sample>/" >&2
echo "" >&2
echo "Note: this lists only probes with >=1 exact hit; COBLAST's eToL Probe Counts" >&2
echo "CSV also includes zero-count probes. Compare on the nonzero counts." >&2
