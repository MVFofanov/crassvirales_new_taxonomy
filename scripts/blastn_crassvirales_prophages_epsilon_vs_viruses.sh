#!/usr/bin/env bash
set -euo pipefail

export NCBI_API_KEY
export EMAIL

crassvirales_prophages_dir="/mnt/c/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages"

IN_DIR="${crassvirales_prophages_dir}/crassvirales_prophages_ncbi_epsilon"
OUT_DIR="${crassvirales_prophages_dir}/blast_prophages_vs_ncbi_and_gtdb/blast_remote_results/crassvirales_prophages_epsilon_vs_viruses"
mkdir -p "$OUT_DIR"

# Optional: API key helps with rate limits
# export NCBI_API_KEY=YOUR_KEY

run_one() {
  local fa="$1"
  # first header token as ID
  local header
  header=$(awk 'NR==1{gsub(/^>/,""); split($0,a,/[ \t]/); print a[1]; exit}' "$fa")
  local safe=${header//[^A-Za-z0-9._-]/_}
  local out="$OUT_DIR/${safe}.blastn.tsv"
  local err="$OUT_DIR/${safe}.stderr.log"

  [[ -s "$out" ]] && { echo "Skip $fa (exists)"; return 0; }

  blastn \
    -query "$fa" \
    -db nt \
    -remote \
    -task dc-megablast \
    -evalue 1e-5 \
    -max_target_seqs 50 \
    -max_hsps 5 \
    -entrez_query 'txid10239[ORGN]' \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen staxids stitle' \
    -out "$out" \
    2> "$err"
}

run_with_retry() {
  local fa="$1"
  local tries=0 max=3
  until run_one "$fa"; do
    ((tries++)) || true
    if (( tries >= max )); then
      echo "FAILED after $tries tries: $fa" | tee -a "$OUT_DIR/failed.log"
      return 1
    fi
    echo "Retry $tries for $fa..."
    sleep $((5 * tries))   # backoff
  done
}

shopt -s nullglob
for fa in "$IN_DIR"/*.fna; do
  echo "==> BLAST: $fa"
  run_with_retry "$fa" || true
  sleep 1   # be nice to NCBI
done
