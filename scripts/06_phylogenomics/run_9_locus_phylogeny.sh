#!/usr/bin/env bash
set -Eeuo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"
# --- Resolve tools and parse run options ---
if command -v python >/dev/null; then python_cmd=python
elif command -v python3 >/dev/null; then python_cmd=python3
else echo "Python is required" >&2; exit 127
fi
threads=2; dry_run=false; force=false; requested_stage="all"
expected_samples=72; expected_loci=9; expected_records=648; output_prefix="efaecalis_72_genomes_9_locus_observed_indels"
prokka_root="data/processed/annotations"
while (($#)); do
  case "$1" in
    --threads) threads=$2; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    --force) force=true; shift ;;
    --stage) requested_stage=$2; shift 2 ;;
    --prokka-root) prokka_root=$2; shift 2 ;;
    -h|--help) echo "Usage: $0 [--threads N] [--dry-run] [--force] [--stage NAME] [--prokka-root DIR]"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "--threads must be a positive integer" >&2; exit 2; }
# --- Define versioned analysis, table and figure locations ---
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
base="analysis/phylogenomics/72_genomes_9_locus_observed_indels"
tables_dir="results/tables/phylogeny_72_genomes_9_locus_observed_indels"
figures_dir="results/figures/phylogeny_72_genomes_9_locus_observed_indels"
if "$force"; then work_root="$base/runs/$timestamp"; else work_root="$base"; fi
log_root="$work_root/logs/$timestamp"
mkdir -p "$log_root"
summary="$log_root/run_summary.tsv"
printf 'stage\tcommand\tstart_utc\tfinish_utc\truntime_seconds\texit_status\tinputs\toutputs\n' > "$summary"

# --- Execute one logged, auditable pipeline stage ---
run_stage() {
  local stage=$1 inputs=$2 outputs=$3; shift 3
  local command_string start finish elapsed status
  printf -v command_string '%q ' "$@"
  echo "+ $command_string"
  if "$dry_run"; then
    printf '%s\t%s\t%s\t%s\t0\tDRY_RUN\t%s\t%s\n' "$stage" "$command_string" "$timestamp" "$timestamp" "$inputs" "$outputs" >> "$summary"
    return 0
  fi
  start=$(date -u +%Y-%m-%dT%H:%M:%SZ); SECONDS=0
  set +e
  "$@" > >(tee "$log_root/${stage}.stdout.log") 2> >(tee "$log_root/${stage}.stderr.log" >&2)
  status=$?
  set -e
  elapsed=$SECONDS; finish=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$stage" "$command_string" "$start" "$finish" "$elapsed" "$status" "$inputs" "$outputs" >> "$summary"
  printf '%s\n' "$command_string" > "$log_root/${stage}.command.txt"
  printf '%s\n' "$status" > "$log_root/${stage}.exit_status.txt"
  [[ "$status" -eq 0 ]] || return "$status"
  printf '%s' "$command_string" | sha256sum > "$log_root/${stage}.command.sha256"
  printf 'stage=%s\nfinish_utc=%s\nexit_status=0\n' "$stage" "$finish" > "$log_root/${stage}.complete"
}
want() { [[ "$requested_stage" == all || "$requested_stage" == "$1" ]]; }
ensure_new() { [[ ! -e "$1" ]] || { echo "Refusing to overwrite $1 (use --force for a versioned run)" >&2; exit 4; }; }

echo "Nine-locus concatenated housekeeping-gene phylogeny (72 genomes; observed disrupted spans; no imputation)"
echo "Start: $timestamp; work root: $work_root; threads: $threads; dry-run: $dry_run"

# --- Policy and software provenance ---
if want policy || [[ "$requested_stage" == all ]]; then
  run_stage 00_policy_consistency "active nine-locus scripts and configuration" "$log_root/00_policy_consistency.complete" \
    "$python_cmd" scripts/06_phylogenomics/00_validate_9_locus_policy.py
fi

if want versions; then
  versions_target="$work_root/software_versions.txt"; ensure_new "$versions_target"
  run_stage versions "environment/environment.yml" "$versions_target" bash -c \
    '{ "$2" --version; samtools --version; mafft --version; (iqtree2 --version || iqtree --version); seqkit version; R --version; } > "$1"' _ "$versions_target" "$python_cmd"
fi
# --- Annotation discovery and reviewed coordinate selection ---
manifest="$work_root/manifests/prokka_annotation_manifest.tsv"
if want discovery || [[ "$requested_stage" == all ]]; then
  ensure_new "$manifest"
  run_stage 01_input_discovery "$prokka_root;data/accession_lists/selected_72_accessions.tsv" "$manifest" \
    "$python_cmd" scripts/06_phylogenomics/01_find_annotation_files.py --prokka-root "$prokka_root" \
      --accessions data/accession_lists/selected_72_accessions.tsv \
      --expected "$expected_samples" --output "$manifest"
fi
coordinates="$work_root/coordinates/selected_gene_coordinates.tsv"
if want candidates || [[ "$requested_stage" == all ]]; then
  [[ -s "$manifest" ]] || { echo "Manifest is required before candidate search" >&2; exit 3; }
  candidates="$work_root/coordinates/all_gene_candidates.tsv"
  raw_coordinates="$work_root/coordinates/selected_gene_coordinates.prokka.tsv"
  ensure_new "$candidates"; ensure_new "$raw_coordinates"; ensure_new "$coordinates"
  override_args=()
  [[ -s config/phylogeny_9_loci_overrides.tsv ]] && override_args=(--overrides config/phylogeny_9_loci_overrides.tsv)
  run_stage 02_gene_candidate_review "$manifest;config/phylogeny_9_loci_genes.tsv" "$candidates;$coordinates" \
    "$python_cmd" scripts/06_phylogenomics/02_extract_gene_coordinates.py --manifest "$manifest" \
      --genes config/phylogeny_9_loci_genes.tsv --all-candidates "$candidates" \
      --selected "$raw_coordinates" --expected-samples "$expected_samples" \
      --expected-loci "$expected_loci" "${override_args[@]}"
  run_stage 02b_observed_span_policy "$raw_coordinates;config/phylogeny_9_loci_observed_span_overrides.tsv" "$coordinates" \
    "$python_cmd" scripts/06_phylogenomics/02b_apply_observed_span_overrides.py \
      --coordinates "$raw_coordinates" --overrides config/phylogeny_9_loci_observed_span_overrides.tsv \
      --output "$coordinates"
fi
# --- Sequence extraction and biological QC ---
sequence_root="$work_root/extracted_sequences"
if want extraction || [[ "$requested_stage" == all ]]; then
  [[ $(($(wc -l < "$coordinates") - 1)) -eq "$expected_records" ]] || { echo "Extraction requires exactly $expected_records reviewed coordinates" >&2; exit 3; }
  run_stage 05_sequence_extraction "$coordinates" "$sequence_root" \
    bash scripts/06_phylogenomics/03_extract_gene_sequences.sh --coordinates "$coordinates" \
      --output-root "$sequence_root" --expected-samples "$expected_samples"
fi
if want sequence_qc || [[ "$requested_stage" == all ]]; then
  run_stage 05b_neighbourhood_qc "$coordinates" "$work_root/qc/gene_neighbourhoods.tsv" \
    "$python_cmd" scripts/06_phylogenomics/04_report_gene_neighbourhoods.py \
      --coordinates "$coordinates" --output "$work_root/qc/gene_neighbourhoods.tsv" --flank-count 2
  run_stage 06_sequence_qc "$sequence_root;$coordinates" "$work_root/qc/sequence_qc.tsv" \
    "$python_cmd" scripts/06_phylogenomics/04_validate_sequences.py --sequence-root "$sequence_root" \
      --coordinates "$coordinates" --qc-dir "$work_root/qc" --tables-dir "$tables_dir" \
      --thresholds config/phylogeny_9_loci_qc_thresholds.tsv \
      --reviewed-findings config/phylogeny_9_loci_primary_findings.tsv \
      --expected-samples "$expected_samples"
fi
if want reference_qc || [[ "$requested_stage" == all ]]; then
  run_stage 06b_pubmlst_qc "$sequence_root;references/phylogeny/pubmlst/efaecalis_2026-08-06" "$tables_dir/pubmlst_locus_qc.tsv" \
    "$python_cmd" scripts/06_phylogenomics/04c_validate_pubmlst.py \
      --sequence-root "$sequence_root" \
      --pubmlst-dir references/phylogeny/pubmlst/efaecalis_2026-08-06 \
      --output "$tables_dir/pubmlst_locus_qc.tsv" --expected-samples "$expected_samples"
  run_stage 06c_refseq_protein_qc "$sequence_root;config/phylogeny_9_loci_protein_references.tsv" "$tables_dir/protein_reference_qc.tsv" \
    "$python_cmd" scripts/06_phylogenomics/04b_validate_proteins.py \
      --sequence-root "$sequence_root" \
      --references config/phylogeny_9_loci_protein_references.tsv \
      --reference-dir references/phylogeny/proteins/efaecalis_v583_refseq_2026-08-06 \
      --findings config/phylogeny_9_loci_primary_findings.tsv \
      --output "$tables_dir/protein_reference_qc.tsv" --expected-samples "$expected_samples"
fi
# --- Alignment, concatenation and phylogenetic inference ---
if want alignment || [[ "$requested_stage" == all ]]; then
  [[ -s "$work_root/qc/sequence_qc.tsv" && -s "$work_root/qc/sequence_qc_failures.tsv" ]] || {
    echo "Completed sequence QC outputs are required before alignment" >&2; exit 3;
  }
  [[ $(($(wc -l < "$work_root/qc/sequence_qc_failures.tsv") - 1)) -eq 0 ]] || {
    echo "Alignment blocked: sequence QC failures remain" >&2; exit 3;
  }
  run_stage 07_mafft_alignment "$sequence_root/by_gene" "$work_root/alignments" \
    bash scripts/06_phylogenomics/05_align_genes.sh --threads "$threads" \
      --input-root "$sequence_root/by_gene" --output-root "$work_root/alignments" \
      --expected-samples "$expected_samples"
fi
if want concatenation || [[ "$requested_stage" == all ]]; then
  run_stage 08_alignment_qc_and_concatenation "$work_root/alignments" "$work_root/concatenated;$work_root/qc" \
    "$python_cmd" scripts/06_phylogenomics/06_concatenate_alignments.py --alignment-dir "$work_root/alignments" \
      --output-dir "$work_root/concatenated" --qc-dir "$work_root/qc" \
      --thresholds config/phylogeny_9_loci_qc_thresholds.tsv --expected-samples "$expected_samples" \
      --output-prefix "$output_prefix" --minimum-concatenated-coverage 0.95
fi
if want iqtree; then
  run_stage 10_iqtree "$work_root/concatenated" "$work_root/iqtree" \
    bash scripts/06_phylogenomics/07_run_iqtree.sh \
      --alignment "$work_root/concatenated/${output_prefix}_alignment.fasta" \
      --partitions "$work_root/concatenated/${output_prefix}_partitions.nex" \
      --output-dir "$work_root/iqtree" --threads "$threads" --output-prefix "$output_prefix"
fi
# --- Metadata-aware publication figures ---
if want visualisation; then
  tree="$work_root/iqtree/${output_prefix}.treefile"
  run_stage 11_metadata_join_and_visualisation "$tree;data/metadata/curated_metadata_72_genomes.csv" "$figures_dir" \
    Rscript scripts/06_phylogenomics/08_visualise_9_locus_tree.R --tree "$tree" \
      --metadata data/metadata/curated_metadata_72_genomes.csv \
      --output-dir "$figures_dir" --tables-dir "$tables_dir" --expected-tips "$expected_samples"
fi
echo "Finish: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Run manifest: $summary"
