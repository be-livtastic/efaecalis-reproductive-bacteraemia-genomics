#!/usr/bin/env bash
set -Eeuo pipefail

# Annotate every selected genome with one accession-stable Prokka directory.
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
genome_root="${1:-$repo_root/data/raw/ncbi_genomes}"
output_root="${2:-$repo_root/data/processed/annotations}"
threads="${THREADS:-2}"

# --- Validate the complete input set before annotation ---
command -v prokka >/dev/null 2>&1 || { echo "ERROR: Prokka is not on PATH." >&2; exit 127; }
[[ -d "$genome_root" ]] || { echo "ERROR: genome directory not found: $genome_root" >&2; exit 2; }
[[ ! -e "$output_root" ]] || { echo "ERROR: refusing to overwrite: $output_root" >&2; exit 3; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: THREADS must be a positive integer." >&2; exit 2; }

mapfile -d '' genomes < <(find "$genome_root" -type f -name '*.fna' -print0 | sort -z)
[[ ${#genomes[@]} -eq 72 ]] || {
  echo "ERROR: expected 72 .fna genomes; found ${#genomes[@]} under $genome_root" >&2
  exit 2
}

mkdir -p "$output_root"
# --- Run one accession-stable annotation directory per genome ---
for genome in "${genomes[@]}"; do
  accession=$(basename "$genome" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' | head -n 1)
  [[ -n "$accession" ]] || { echo "ERROR: no assembly accession in $genome" >&2; exit 2; }
  prokka --outdir "$output_root/$accession" --prefix "$accession" \
    --locustag "${accession//[._]/}" --genus Enterococcus --species faecalis \
    --cpus "$threads" "$genome"
done

echo "Prokka annotations completed under: $output_root"
