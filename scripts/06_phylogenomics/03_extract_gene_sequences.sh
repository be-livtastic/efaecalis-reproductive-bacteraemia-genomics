#!/usr/bin/env bash
set -Eeuo pipefail

# --- Parse and validate command-line inputs ---
usage() { echo "Usage: $0 --coordinates FILE --output-root DIR [--expected-samples N]"; }
coordinates=""; output_root=""; expected_samples=72
while (($#)); do
  case "$1" in
    --coordinates) coordinates=$2; shift 2 ;;
    --output-root) output_root=$2; shift 2 ;;
    --expected-samples) expected_samples=$2; shift 2 ;;
    --force-run-id) echo "Versioned run identifier: $2"; shift 2 ;;
    *) usage >&2; exit 2 ;;
  esac
done
[[ -n "$coordinates" && -n "$output_root" ]] || { usage >&2; exit 2; }
[[ -s "$coordinates" ]] || { echo "Missing or empty coordinates: $coordinates" >&2; exit 2; }
command -v samtools >/dev/null || { echo "samtools is required" >&2; exit 127; }
command -v seqkit >/dev/null || { echo "seqkit is required" >&2; exit 127; }

by_sample="$output_root/by_sample"
by_gene="$output_root/by_gene"
mkdir -p "$by_sample" "$by_gene"

# --- Extract each reviewed interval in strand-aware orientation ---
tail -n +2 "$coordinates" |
while IFS=$'\t' read -r sample accession group gene matched_attribute matched_value feature_type contig start end strand phase locus_tag feature_id gene_name product expected rank reason qc gff fna source; do
  [[ "$qc" == "SELECTED" ]] || { echo "Unvalidated coordinate: $sample/$gene" >&2; exit 3; }
  sample_dir="$by_sample/$accession"
  mkdir -p "$sample_dir"
  destination="$sample_dir/${gene}.fasta"
  [[ ! -e "$destination" ]] || { echo "Refusing to overwrite $destination" >&2; exit 4; }
  [[ -s "$fna" ]] || { echo "Missing FNA: $fna" >&2; exit 4; }
  if [[ ! -s "${fna}.fai" || "$fna" -nt "${fna}.fai" ]]; then
    echo "+ samtools faidx $fna"
    samtools faidx "$fna"
  fi
  region="${contig}:${start}-${end}"
  echo "+ samtools faidx $fna $region"
  if [[ "$strand" == "-" ]]; then
    samtools faidx "$fna" "$region" | seqkit seq --quiet --seq-type dna -r -p |
      awk -v id="$accession" 'BEGIN{print ">" id} !/^>/{printf "%s",$0} END{print ""}' > "$destination"
  elif [[ "$strand" == "+" ]]; then
    samtools faidx "$fna" "$region" |
      awk -v id="$accession" 'BEGIN{print ">" id} !/^>/{printf "%s",$0} END{print ""}' > "$destination"
  else
    echo "Invalid strand for $sample/$gene: $strand" >&2; exit 4
  fi
done

# --- Confirm every locus contains the expected sample set ---
for gene in gdh gyd pstS gki xpt yqiL pyrC groEL recA; do
  destination="$by_gene/${gene}_${expected_samples}_sequences.fasta"
  [[ ! -e "$destination" ]] || { echo "Refusing to overwrite $destination" >&2; exit 4; }
  mapfile -t files < <(find "$by_sample" -mindepth 2 -maxdepth 2 -type f -name "${gene}.fasta" | sort)
  [[ ${#files[@]} -eq "$expected_samples" ]] || { echo "$gene: expected $expected_samples individual sequences; found ${#files[@]}" >&2; exit 5; }
  for file in "${files[@]}"; do
    sed -n '1,2p' "$file"
  done > "$destination"
  count=$(grep -c '^>' "$destination")
  unique=$(grep '^>' "$destination" | sort -u | wc -l)
  [[ "$count" -eq "$expected_samples" && "$unique" -eq "$expected_samples" ]] || { echo "$gene identifiers failed validation" >&2; exit 5; }
  echo "$gene: sequences=$count unique_ids=$unique"
done
