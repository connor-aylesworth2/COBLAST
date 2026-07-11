#!/usr/bin/env bash
# Decompress + convert fastq.gz -> fasta in the current directory.
# Run from the folder containing the .fastq.gz files.
set -euo pipefail

read -rp "Enter the 8 fastq.gz filenames, separated by commas: " reply

IFS=',' read -ra files <<< "$reply"

for f in "${files[@]}"; do
    f="$(echo "$f" | xargs)"          # trim surrounding whitespace
    [ -z "$f" ] && continue
    if [ ! -f "$f" ]; then
        echo "SKIP: '$f' not found" >&2
        continue
    fi
    out="${f%.gz}"; out="${out%.fastq}"; out="${out%.fq}.fasta"
    # zcat = decompress; sed maps 4-line fastq records to 2-line fasta.
    zcat -- "$f" | sed -n '1~4s/^@/>/p;2~4p' > "$out"
    echo "OK:   $f -> $out"
done
