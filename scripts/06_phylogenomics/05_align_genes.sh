#!/usr/bin/env bash
set -Eeuo pipefail
threads=1; input_root=""; output_root=""
while (($#)); do
  case "$1" in
    --threads) threads=$2; shift 2 ;;
    --input-root) input_root=$2; shift 2 ;;
    --output-root) output_root=$2; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$input_root" && -n "$output_root" ]] || { echo "--input-root and --output-root are required" >&2; exit 2; }
command -v mafft >/dev/null || { echo "mafft is required" >&2; exit 127; }
mkdir -p "$output_root"
for gene in gdh gyd pstS gki aroE xpt yqiL pyrC groEL recA; do
  source_fasta="$input_root/${gene}_72_sequences.fasta"
  destination="$output_root/${gene}.aligned.fasta"
  stderr_log="$output_root/${gene}.mafft.stderr.txt"
  [[ -s "$source_fasta" ]] || { echo "Missing $source_fasta" >&2; exit 3; }
  [[ ! -e "$destination" && ! -e "$stderr_log" ]] || { echo "Refusing to overwrite $destination or $stderr_log" >&2; exit 4; }
  echo "+ mafft --auto --thread $threads $source_fasta"
  start=$SECONDS
  if mafft --auto --thread "$threads" "$source_fasta" > "$destination" 2> "$stderr_log"; then status=0; else status=$?; fi
  echo "$gene: exit_status=$status runtime_seconds=$((SECONDS-start)) sequences=$(grep -c '^>' "$destination")"
  [[ "$status" -eq 0 ]] || exit "$status"
done
