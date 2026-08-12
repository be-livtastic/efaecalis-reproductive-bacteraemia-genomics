#!/usr/bin/env bash
set -euo pipefail

# Portable NCBI retrieval script. It uses repository-relative paths and can be
# run from any working directory on Linux or WSL. NCBI Datasets must be on PATH.
# Existing archives or extraction directories are never overwritten.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
ACCESSIONS="${1:-$PROJECT_ROOT/data/accession_lists/selected_72_accessions.tsv}"
OUTPUT_ROOT="${2:-$PROJECT_ROOT/data/raw/ncbi_genomes}"

# --- Validate tools and inputs before creating output ---
if ! command -v datasets >/dev/null 2>&1; then
    echo "ERROR: NCBI Datasets CLI is not available on PATH." >&2
    exit 1
fi
if [[ ! -f "$ACCESSIONS" ]]; then
    echo "ERROR: accession manifest not found: $ACCESSIONS" >&2
    exit 1
fi
if [[ -e "$OUTPUT_ROOT" ]]; then
    echo "ERROR: refusing to overwrite existing output: $OUTPUT_ROOT" >&2
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"
# --- Download each comparison group independently ---
for category in Reproductive Bacteraemia; do
    category_lower="${category,,}"
    accession_file="$OUTPUT_ROOT/${category_lower}_accessions.txt"
    archive="$OUTPUT_ROOT/${category_lower}_genomes.zip"
    extract_dir="$OUTPUT_ROOT/${category_lower}"

    # Derive the download list from the canonical two-column manifest.
    awk -F $'\t' -v category="$category" \
        'NR > 1 && $2 == category {print $1}' "$ACCESSIONS" > "$accession_file"
    expected=14
    [[ "$category" == "Bacteraemia" ]] && expected=58
    observed="$(awk 'NF {n++} END {print n+0}' "$accession_file")"
    if [[ "$observed" -ne "$expected" ]]; then
        echo "ERROR: expected $expected $category accessions; found $observed." >&2
        exit 1
    fi

    # Retain sequence and annotation inputs needed by downstream stages.
    datasets download genome accession \
        --inputfile "$accession_file" \
        --include genome,protein,gff3,seq-report \
        --filename "$archive"
    mkdir "$extract_dir"
    unzip -q "$archive" -d "$extract_dir"
done

echo "NCBI genome retrieval completed under: $OUTPUT_ROOT"
