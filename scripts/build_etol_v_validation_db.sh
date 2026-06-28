#!/usr/bin/env bash
#
# build_etol_v_validation_db.sh
#
# Build the LOCAL contig-validation database for the eToL-V (viral) preset, the
# same way you built ToL_rRNA and nt_human_9606. Run it in a screen/tmux session
# on the server (it downloads ~1 GB and builds a ~3-4 GB BLAST DB).
#
# WHY THIS DB EXISTS (read before swapping in a different one):
#   The eToL-V net recruits reads with viral structural-protein probes, assembles
#   them into contigs, then COBLAST+ (contig_id.py) BLASTs each contig against a
#   reference DB to (a) name the closest homolog and (b) DROP any contig whose
#   closest homolog is "Homo sapiens". In the eToL-V dissertation the dominant
#   false positive was herpesvirus probes recruiting reads that actually came
#   from the human MITOCHONDRIAL genome and assembled into herpes-looking contigs.
#   This validation step is what kills them -- but ONLY if the reference DB
#   contains human (incl. mito) sequences so those artifact contigs get a
#   "Homo sapiens" top hit. A viruses-only DB would give them NO hit and let them
#   through as real virus. So this DB deliberately mixes:
#       1. RefSeq viral genomes  -> genuine viral contigs get a viral homolog and
#                                   survive (AdC penton, SARS-CoV-2 S, etc.).
#       2. Human genome GRCh38   -> artifact contigs (esp. herpes->mito) get a
#          (includes the mito)     "Homo sapiens ..." closest homolog and are
#                                   dropped by contig_id.py's human-contig filter.
#   This is exactly the memo's "RefSeq-viral + human-genome + human-mito" build.
#   Do NOT reuse ToL_rRNA here: viruses have no rRNA, so a mito artifact would not
#   match SILVA at all and would slip through.
#
# WHAT contig_id.py NEEDS FROM THE DEFLINES:
#   It only substring-matches "Homo sapiens" in the BLAST stitle. NCBI RefSeq
#   human deflines contain it (e.g. ">NC_012920.1 Homo sapiens mitochondrion,
#   complete genome"), so no defline munging is needed.
#
set -euo pipefail

# ============================== CONFIG ==============================
# Where the finished BLAST DB should live. Point this at the SAME folder where
# ToL_rRNA and nt_human_9606 already sit so COBLAST+ can find them together.
DEST_DIR="/home/s2837738/COBLAST_2.0/blast_dbs"

# DB name (the prefix you will register in COBLAST+). Mirrors your ToL_rRNA name.
DB_NAME="ToL_virus_val"

# CPU threads for makeblastdb's input parsing (modest effect; set to taste).
THREADS=8

# Scratch dir for downloads + the concatenated FASTA. Needs ~6 GB free.
BUILD_DIR="${DEST_DIR}/_build_${DB_NAME}"

# Delete the multi-GB intermediate FASTA + downloads after a successful build?
# 1 = clean up (recommended), 0 = keep for inspection / re-runs.
CLEANUP=1

# NCBI sources -------------------------------------------------------
# RefSeq viral release (all viral genomic FASTA shards: viral.1.*, viral.2.*, ...)
REFSEQ_VIRAL_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
# Human reference genome, RefSeq GRCh38.p14 (one file; includes every chromosome,
# unplaced scaffolds, and the mitochondrion NC_012920.1).
HUMAN_GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
# ====================================================================

# --- sanity checks ---
for tool in wget zcat makeblastdb blastdbcmd; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' is not on PATH." >&2; exit 1; }
done

mkdir -p "$DEST_DIR" "$BUILD_DIR"
VIRAL_DIR="$BUILD_DIR/refseq_viral"
HUMAN_GZ="$BUILD_DIR/GRCh38.p14_genomic.fna.gz"
COMBINED="$BUILD_DIR/etol_v_validation.fna"
DB_PREFIX="$DEST_DIR/$DB_NAME"
mkdir -p "$VIRAL_DIR"

echo "==> eToL-V contig-validation DB build" >&2
echo "    Destination : $DB_PREFIX" >&2
echo "    Build dir   : $BUILD_DIR" >&2
echo "" >&2

# --- 1. RefSeq viral genomes ---
# Recursive, flattened (-nd), no-parent (-np) fetch of every viral genomic shard.
echo "==> [1/4] Downloading RefSeq viral genomes ..." >&2
wget -q --show-progress -r -np -nd \
     -A 'viral.*.genomic.fna.gz' \
     -P "$VIRAL_DIR" \
     "$REFSEQ_VIRAL_URL"
VIRAL_COUNT=$(find "$VIRAL_DIR" -name 'viral.*.genomic.fna.gz' | wc -l)
[ "$VIRAL_COUNT" -gt 0 ] || { echo "ERROR: no RefSeq viral shards downloaded." >&2; exit 1; }
echo "    got $VIRAL_COUNT viral shard(s)." >&2

# --- 2. Human genome GRCh38 (with mitochondrion) ---
echo "==> [2/4] Downloading human genome GRCh38.p14 (~1 GB) ..." >&2
wget -q --show-progress -O "$HUMAN_GZ" "$HUMAN_GENOME_URL"

# --- 3. Decompress + concatenate into one FASTA ---
# zcat streams each gz's decompressed bytes; concatenating them yields one valid
# multi-FASTA. makeblastdb cannot read .gz directly, hence this step.
echo "==> [3/4] Decompressing + concatenating into one FASTA ..." >&2
zcat "$VIRAL_DIR"/viral.*.genomic.fna.gz "$HUMAN_GZ" > "$COMBINED"

# Sanity: confirm the two ingredients that make the validation work are present.
SEQ_TOTAL=$(grep -c '^>' "$COMBINED" || true)
HUMAN_HDRS=$(grep -c '^>.*Homo sapiens' "$COMBINED" || true)
MITO_HDRS=$(grep -c '^>.*[Mm]itochondrion' "$COMBINED" || true)
echo "    sequences total        : $SEQ_TOTAL" >&2
echo "    'Homo sapiens' headers : $HUMAN_HDRS   (must be > 0, or artifacts won't be dropped)" >&2
echo "    'mitochondrion' headers: $MITO_HDRS    (must be > 0; this is what catches herpes->mito)" >&2
[ "$HUMAN_HDRS" -gt 0 ] || { echo "ERROR: no human headers -> human/mito artifacts would slip through. Aborting." >&2; exit 1; }
[ "$MITO_HDRS" -gt 0 ]  || { echo "ERROR: no mitochondrion header found. Aborting." >&2; exit 1; }

# --- 4. Build the BLAST nucleotide DB ---
# No -parse_seqids: contig_id.py uses this only as a blastn '-db' target and reads
# 'stitle', so seqid indexing is not needed (and -parse_seqids on the whole human
# genome is slow and can choke on long deflines). Add it back only if you want
# blastdbcmd retrieval-by-accession from this DB.
echo "==> [4/4] Running makeblastdb ..." >&2
makeblastdb \
  -in "$COMBINED" \
  -dbtype nucl \
  -title "eToL-V contig validation: RefSeq viral + human GRCh38 (incl. mito)" \
  -out "$DB_PREFIX"

# --- verify + report ---
echo "" >&2
echo "==> Verifying DB ..." >&2
blastdbcmd -db "$DB_PREFIX" -info >&2

if [ "$CLEANUP" -eq 1 ]; then
  echo "==> Cleaning up intermediates (CLEANUP=1) ..." >&2
  rm -rf "$BUILD_DIR"
fi

cat >&2 <<EOF

==================================================================
 DONE. Validation DB prefix:
   $DB_PREFIX
------------------------------------------------------------------
 Register it in COBLAST+ (Databases page):
   - Prefix path : $DB_PREFIX
   - Type        : nucleotide
   - Category    : reference   (any nucl category works; "reference"
                   keeps it alongside ToL_rRNA in the picker)
 Then for an eToL-V run, on the Batch page select:
   - Preset                : eToL-V (viruses)
   - Secondary human filter: nt_human_9606
   - Assemble contigs      : on
   - Identify contigs      : on  -> reference DB = $DB_NAME  (NOT ToL_rRNA)
==================================================================
EOF
