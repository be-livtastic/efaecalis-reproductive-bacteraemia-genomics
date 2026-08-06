#!/usr/bin/env bash
set -Eeuo pipefail
threads=2; input_root=""; output_root=""; expected_samples=72
while (($#)); do
  case "$1" in
    --threads) threads=$2; shift 2 ;;
    --input-root) input_root=$2; shift 2 ;;
    --output-root) output_root=$2; shift 2 ;;
    --expected-samples) expected_samples=$2; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$input_root" && -n "$output_root" ]] || { echo "--input-root and --output-root are required" >&2; exit 2; }
command -v mafft >/dev/null || { echo "mafft is required" >&2; exit 127; }
mkdir -p "$output_root"
for gene in gdh gyd pstS gki xpt yqiL pyrC groEL recA; do
  source_fasta="$input_root/${gene}_${expected_samples}_sequences.fasta"
  destination="$output_root/${gene}.aligned.fasta"
  stderr_log="$output_root/${gene}.mafft.stderr.txt"
  checksum="$output_root/${gene}.aligned.fasta.sha256"
  marker="$output_root/${gene}.complete"
  [[ -s "$source_fasta" ]] || { echo "Missing $source_fasta" >&2; exit 3; }
  [[ ! -e "$destination" && ! -e "$stderr_log" && ! -e "$checksum" && ! -e "$marker" ]] || {
    echo "Refusing to overwrite existing output for $gene" >&2; exit 4;
  }
  temporary_alignment=$(mktemp "$output_root/.${gene}.alignment.XXXXXX")
  temporary_log=$(mktemp "$output_root/.${gene}.mafft.XXXXXX")
  cleanup() { rm -f "$temporary_alignment" "$temporary_log"; }
  trap cleanup EXIT INT TERM
  echo "+ mafft --auto --thread $threads $source_fasta"
  start=$SECONDS
  if mafft --auto --thread "$threads" "$source_fasta" > "$temporary_alignment" 2> "$temporary_log"; then status=0; else status=$?; fi
  sequence_count=$(grep -c '^>' "$temporary_alignment" || true)
  echo "$gene: exit_status=$status runtime_seconds=$((SECONDS-start)) sequences=$sequence_count"
  if [[ "$status" -ne 0 ]]; then exit "$status"; fi
  [[ "$sequence_count" -eq "$expected_samples" ]] || { echo "$gene: expected $expected_samples aligned records" >&2; exit 3; }
  mv "$temporary_alignment" "$destination"
  mv "$temporary_log" "$stderr_log"
  sha256sum "$destination" > "$checksum"
  printf 'gene=%s\nsequences=%s\nstrategy=--auto\nthreads=%s\n' "$gene" "$expected_samples" "$threads" > "$marker"
  trap - EXIT INT TERM
done
