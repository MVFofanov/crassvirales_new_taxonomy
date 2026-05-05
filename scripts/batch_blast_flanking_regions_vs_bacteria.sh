#!/usr/bin/env bash
# batch_blast_bacteria_10kb.sh
# Remote BLASTN vs nt for all {prophage_id}/*flanking_10kb_*.fna under a base directory.
# Output: {prophage_id}/blast_results/{file}.tsv
set -euo pipefail

# ---- optional: conda ----
if [[ -f "${HOME}/mambaforge/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1091
  source "${HOME}/mambaforge/etc/profile.d/conda.sh" && conda activate blast_env || true
fi

BASE_DIR="${1:-}"
if [[ -z "$BASE_DIR" || ! -d "$BASE_DIR" ]]; then
  echo "Usage: $0 /path/to/base_dir" >&2
  exit 1
fi

# ---- Config (override via env) ----
BLASTN_BIN="${BLASTN:-blastn}"   # or absolute path
DB="${DB:-nt}"

# Bacteria only; excludes environmental/metagenome entries
ENTREZ_QUERY="${ENTREZ_QUERY:-'(txid2[ORGN]) NOT metagenomes[filter] NOT \"environmental samples\"[organism]'}"

# Default to megablast for near-identical flanks (faster, cleaner).
# Switch to dc-megablast by: TASK=dc-megablast EVALUE=1e-10
TASK="${TASK:-megablast}"
EVALUE="${EVALUE:-1e-25}"

# Megablast tuning (ignored if TASK=dc-megablast)
WORD_SIZE="${WORD_SIZE:-28}"
PERC_ID="${PERC_ID:-60}"

# Downrank short/duplicated HSPs; keep only strong, longer matches
# QCOV="${QCOV:-80}"               # qcov_hsp_perc
CULL="${CULL:-1}"                # culling_limit
SUBJECT_BESTHIT="${SUBJECT_BESTHIT:-1}"

MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-50}"
MAX_HSPS="${MAX_HSPS:-1}"

# NOTE: omit staxids/sscinames to avoid local taxdb warnings in remote mode
OUTFMT="${OUTFMT:-6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle}"

# Retries & pacing to avoid Entrez aborts/rate limits
ATTEMPTS="${ATTEMPTS:-3}"              # total attempts per file
BACKOFF_BASE="${BACKOFF_BASE:-10}"     # seconds; grows per attempt
PAUSE_BETWEEN="${PAUSE_BETWEEN:-6}"    # pause after success
JITTER_MAX="${JITTER_MAX:-5}"          # add 0..JITTER_MAX random seconds
OVERWRITE="${OVERWRITE:-1}"

if ! command -v "$BLASTN_BIN" >/dev/null 2>&1; then
  echo "Error: blastn not found at '$BLASTN_BIN'." >&2
  exit 1
fi

run_blast() {
  local query="$1" out="$2"
  "$BLASTN_BIN" -remote -db "$DB" \
    -entrez_query "$ENTREZ_QUERY" \
    -query "$query" \
    -task megablast \
    -evalue 0.05 \
    -word_size 28 \
    -reward 1 -penalty -2 \
    -gapopen 0 -gapextend 2 \
    -soft_masking true -dust yes \
    -max_target_seqs 100 \
    -max_hsps 0 \
    -outfmt "$OUTFMT" -out "$out"
    # If you want species-specific repeats (e.g., human), add:
    # -window_masker_taxid 9606
}


should_long_sleep() {
  # Return 0 if stderr shows known remote issues
  local errfile="$1"
  grep -Eq 'search aborted by Entrez|Connection stream is in bad state|CPU usage limit was exceeded' "$errfile"
}

# --- ONLY 10 kb files ---
while IFS= read -r -d '' fna; do
  prophage_dir="$(dirname "$fna")"
  fasta_file="$(basename "$fna")"
  base_noext="${fasta_file%.fna}"

  outdir="${prophage_dir}/blast_results"
  outfile="${outdir}/${base_noext}.tsv"
  mkdir -p "$outdir"

  if [[ "$OVERWRITE" -ne 1 && -s "$outfile" ]]; then
    echo "[SKIP] $outfile exists (use OVERWRITE=1 to redo)"
    continue
  fi

  # sanitized temp copy to avoid issues with '&' etc.
  safe_name="${fasta_file//&/_}"
  tmp_query="${prophage_dir}/.${safe_name}.tmp.fna"
  cp -f -- "$fna" "$tmp_query"

  rel_in="$(realpath --relative-to="$BASE_DIR" "$fna" 2>/dev/null || echo "$fna")"
  rel_out="$(realpath --relative-to="$BASE_DIR" "$outfile" 2>/dev/null || echo "$outfile")"
  echo "[RUN ] ${rel_in} -> ${rel_out}"

  success=0
  errlog="$(mktemp)"
  for i in $(seq 1 "$ATTEMPTS"); do
    if run_blast "$tmp_query" "$outfile" 2> "$errlog"; then
      success=1
      # polite pause after success
      sleep $(( PAUSE_BETWEEN + RANDOM % (JITTER_MAX+1) ))
      break
    else
      if should_long_sleep "$errlog"; then
        wait_s=$(( BACKOFF_BASE * i + 20 + RANDOM % (JITTER_MAX+1) ))
        echo "[WARN] Remote refused/aborted; sleeping ${wait_s}s (attempt $i/$ATTEMPTS)" >&2
        sleep "$wait_s"
      else
        wait_s=$(( BACKOFF_BASE * i + RANDOM % (JITTER_MAX+1) ))
        echo "[WARN] BLAST failed; sleeping ${wait_s}s (attempt $i/$ATTEMPTS)" >&2
        sleep "$wait_s"
      fi
    fi
  done
  rm -f -- "$tmp_query" "$errlog"

  if [[ "$success" -ne 1 ]]; then
    echo "[FAIL] BLAST ultimately failed for: $fna" >&2
    echo "$fna" >> "${BASE_DIR%/}/blast_failed.list"
  fi
done < <(find "$BASE_DIR" -type f -name "*flanking_10kb_*.fna" -print0)

echo "Done."
