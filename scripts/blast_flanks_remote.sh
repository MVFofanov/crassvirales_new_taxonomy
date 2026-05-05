#!/usr/bin/env bash
set -euo pipefail

# Root with prophage subdirs (override by passing as $1)
ROOT="${1:-/mnt/c/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/integrated_prophages/integrated_prophages_flanking_regions}"

# BLAST knobs (kept identical to your example; customize via env if needed)
DB="${DB:-nt}"                             # e.g., nt or refseq_genomic
TASK="${TASK:-dc-megablast}"               # dc-megablast or megablast
ENTREZQ="${ENTREZQ:-txid2[ORGN]}"          # bacteria
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-50}"
MAX_HSPS="${MAX_HSPS:-1}"
OUTFMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle'

# Overwrite behavior: set OVERWRITE=1 to force rerun even if output exists
OVERWRITE="${OVERWRITE:-1}"

# Retry settings (handles common NCBI -remote hiccups)
RETRIES="${RETRIES:-3}"
BASE_SLEEP="${BASE_SLEEP:-5}"

command -v blastn >/dev/null 2>&1 || { echo "ERROR: blastn not found in PATH."; exit 1; }

echo "Root: $ROOT"
echo "DB: $DB | TASK: $TASK | ENTREZ: $ENTREZQ"
echo "Will process 10kb left/right flank FASTAs in each subdir."

# Find only files that are *individual* 10kb flanks (exclude *_all_*.fna)
# Structure assumed: $ROOT/<prophage_id>/*.fna
# Depth limited to one level of subdirectories.
while IFS= read -r -d '' fasta; do
  base="$(basename "$fasta" .fna)"
  dir="$(dirname "$fasta")"

  # Keep only _flanking_10kb_left_* or _flanking_10kb_right_* (exclude _all_*)
  case "$base" in
    *_flanking_10kb_left_*|*_flanking_10kb_right_*) ;;
    *) continue ;;
  esac

  outdir="$dir/blast"
  outfile="$outdir/${base}_blast.txt"
  tmpout="$outfile.tmp"

  mkdir -p "$outdir"

  if [[ "$OVERWRITE" != "1" && -s "$outfile" ]]; then
    echo "SKIP (exists): $outfile"
    continue
  fi

  echo "BLAST: $fasta"
  ok=0
  for try in $(seq 1 "$RETRIES"); do
    # Clean temp before attempt
    : > "$tmpout" || true

    if blastn -remote \
      -db "$DB" \
      -query "$fasta" \
      -entrez_query "$ENTREZQ" \
      -task "$TASK" \
      -max_target_seqs "$MAX_TARGET_SEQS" \
      -max_hsps "$MAX_HSPS" \
      -outfmt "$OUTFMT" \
      -out "$tmpout"; then
      # Success (remote finished without RPC error)
      if [[ -s "$tmpout" ]]; then
        mv -f "$tmpout" "$outfile"
        ok=1
        echo "OK -> $outfile"
        break
      fi
    fi

    # Retry with backoff
    sleep_time=$((BASE_SLEEP * try))
    echo "WARN: attempt $try failed for $fasta, retrying in ${sleep_time}s..."
    sleep "$sleep_time"
  done

  # If still not ok, keep temp (for debugging) and move on
  if [[ "$ok" -ne 1 ]]; then
    echo "FAIL: $fasta (see $tmpout if present)"
  fi

done < <(find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name '*.fna' -print0)

echo "Done."
