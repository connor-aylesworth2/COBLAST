#!/usr/bin/env bash
#
# run_manual_blast.sh
#
# Run ONE ordinary BLAST+ search (one query vs one local database) straight from
# the command line and save the result, with NO COBLAST+ involvement. This is the
# manual baseline you compare COBLAST+ against by hand: it reproduces exactly what
# COBLAST+'s general "Run BLAST" page does for a plain search with every Advanced
# setting left blank.
#
# Why this matches COBLAST+'s general search (traced through blast_runner.run_blast):
#   * Program blastn, default task "megablast" -> COBLAST passes "-task megablast"
#     explicitly, which is also command-line blastn's own default, so results are
#     identical with or without the flag.
#   * Every Advanced field left blank is omitted, so BLAST+ applies its OWN
#     defaults: e-value 10 and max_target_seqs 500. This script likewise passes
#     NO -evalue / -max_target_seqs / -perc_identity, so both sides use the same
#     BLAST defaults. (This is the key difference from the eToL preset, which
#     would force -max_target_seqs 5000000 and the permissive net; none of that
#     applies here.)
#   * Tabular output with COBLAST's column set:
#       6 qseqid sseqid stitle pident length qcovs evalue bitscore
#   * -num_threads is the only always-present extra flag. Threads change the
#     SPEED of a search, never which hits are found, so any value is a fair
#     comparison. COBLAST omits -mt_mode for a general (non-preset) search, so
#     this script omits it too.
#
# To compare by hand:
#   1. Run this script; open the saved .tsv.
#   2. In COBLAST+, open "Run BLAST" (NOT Batch, NO preset), choose blastn, paste
#      the same query, select the same database, leave Advanced settings blank,
#      and run. Download the result as TSV.
#   3. The two tables should be identical (row order aside).
#
set -euo pipefail

# ============================== CONFIG ==============================
# BLAST DB prefix = the path to the database files WITHOUT any extension. For the
# multi-volume SRR21676099 DB (a .nal alias plus .00.* / .01.* volumes) the
# prefix is the shared stem below.
DB_PREFIX="/home/s2837738/COBLAST_2.0/SRA_data/SRR21676099/blastdb/SRR21676099"

# BLAST program. blastn for a nucleotide query vs a nucleotide database.
PROGRAM="blastn"

# CPU threads. Affects speed only, never the hit set.
THREADS=8

# Where the query FASTA and results are written.
OUT_DIR="./manual_blast_results"
# ====================================================================

# The test query. PGK1_2 is just a convenient 64 bp nucleotide sequence to drive
# the search; nothing here is eToL-specific. Swap PROBE_ID/PROBE_SEQ (or point
# the script at your own FASTA) to test any other query.
PROBE_ID="PGK1_2"
PROBE_SEQ="TGATGAAGAGGGAGCCAAGATTGTCAAAGACCTAATGTCCAAAGCTGAGAAGAATGGTGTGAAG"

# blastn's default task; harmless to state explicitly (it is the CLI default too).
TASK="megablast"
OUTFMT="6 qseqid sseqid stitle pident length qcovs evalue bitscore"

# --- sanity checks ---
command -v "$PROGRAM" >/dev/null 2>&1 || { echo "ERROR: '$PROGRAM' is not on PATH." >&2; exit 1; }
if ! blastdbcmd -db "$DB_PREFIX" -info >/dev/null 2>&1; then
  echo "ERROR: BLAST cannot open the database at:" >&2
  echo "         $DB_PREFIX" >&2
  echo "       Point DB_PREFIX at the file stem (no extension). For a multi-volume" >&2
  echo "       DB that is the prefix shared by the .nal alias and the .00.* files." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
QUERY="$OUT_DIR/${PROBE_ID}.fasta"
RESULTS="$OUT_DIR/${PROBE_ID}_vs_$(basename "$DB_PREFIX").tsv"
CMD_LOG="$OUT_DIR/${PROBE_ID}_vs_$(basename "$DB_PREFIX").cmd.txt"

printf '>%s\n%s\n' "$PROBE_ID" "$PROBE_SEQ" > "$QUERY"

# Build the command exactly as COBLAST+'s general search does for a blank-Advanced
# blastn run: program, query, db, outfmt, task, num_threads. No other flags, so
# BLAST+ defaults (e-value 10, max_target_seqs 500) apply.
cmd=( "$PROGRAM"
      -query "$QUERY"
      -db "$DB_PREFIX"
      -outfmt "$OUTFMT"
      -task "$TASK"
      -num_threads "$THREADS" )

# Record the exact invocation alongside the results for reproducibility.
printf '%q ' "${cmd[@]}" > "$CMD_LOG"; printf '\n' >> "$CMD_LOG"

echo "Program  : $PROGRAM" >&2
echo "Query    : $QUERY ($PROBE_ID, ${#PROBE_SEQ} bp)" >&2
echo "Database : $DB_PREFIX" >&2
echo "Command  : ${cmd[*]}" >&2
echo "" >&2

# outfmt 6 prints one tab-delimited line per HSP, no header.
"${cmd[@]}" > "$RESULTS"

HITS=$(grep -c . "$RESULTS" || true)
echo "Hits (rows) : $HITS   (capped at BLAST's default max_target_seqs=500)" >&2
echo "Results     : $RESULTS" >&2
echo "Command log : $CMD_LOG" >&2
echo "" >&2
echo "Columns: qseqid  sseqid  stitle  pident  length  qcovs  evalue  bitscore" >&2
if (( HITS > 0 )); then
  echo "First rows:" >&2
  head -n 5 "$RESULTS" >&2
fi
