#!/usr/bin/env bash
set -Eeuo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"
if command -v python >/dev/null; then python_cmd=python
elif command -v python3 >/dev/null; then python_cmd=python3
else echo "Python is required" >&2; exit 127
fi
threads=1; dry_run=false; force=false; requested_stage="all"
prokka_root="local_archive/large_outputs/Phylogeny_project/phylogeny_work/01_annotations"
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
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
base="analysis/phylogenomics/10_locus"
if "$force"; then work_root="$base/runs/$timestamp"; else work_root="$base"; fi
log_root="$work_root/logs/$timestamp"
mkdir -p "$log_root"
summary="$log_root/run_summary.tsv"
printf 'stage\tcommand\tstart_utc\tfinish_utc\truntime_seconds\texit_status\tinputs\toutputs\n' > "$summary"

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
}
want() { [[ "$requested_stage" == all || "$requested_stage" == "$1" ]]; }
ensure_new() { [[ ! -e "$1" ]] || { echo "Refusing to overwrite $1 (use --force for a versioned run)" >&2; exit 4; }; }

echo "Ten-locus concatenated housekeeping-gene phylogeny"
echo "Start: $timestamp; work root: $work_root; threads: $threads; dry-run: $dry_run"

if want versions; then
  versions_target="$work_root/software_versions.txt"; ensure_new "$versions_target"
  run_stage versions "environment/environment.yml" "$versions_target" bash -c \
    '{ "$2" --version; samtools --version; mafft --version; (iqtree2 --version || iqtree --version); seqkit version; R --version; } > "$1"' _ "$versions_target" "$python_cmd"
fi
manifest="$work_root/manifests/prokka_annotation_manifest.tsv"
if want discovery || [[ "$requested_stage" == all ]]; then
  ensure_new "$manifest"
  run_stage 01_input_discovery "$prokka_root;data/accession_lists/selected_72_accessions.tsv" "$manifest" \
    "$python_cmd" scripts/06_phylogenomics/01_find_annotation_files.py --prokka-root "$prokka_root" --output "$manifest"
fi
coordinates="$work_root/coordinates/selected_gene_coordinates.tsv"
if want candidates || [[ "$requested_stage" == all ]]; then
  [[ -s "$manifest" ]] || { echo "Manifest is required before candidate search" >&2; exit 3; }
  candidates="$work_root/coordinates/all_gene_candidates.tsv"
  ensure_new "$candidates"; ensure_new "$coordinates"
  override_args=()
  [[ -s config/phylogeny_10_loci_overrides.tsv ]] && override_args=(--overrides config/phylogeny_10_loci_overrides.tsv)
  run_stage 02_gene_candidate_review "$manifest;config/phylogeny_10_loci_genes.tsv" "$candidates;$coordinates" \
    "$python_cmd" scripts/06_phylogenomics/02_extract_gene_coordinates.py --manifest "$manifest" \
      --all-candidates "$candidates" --selected "$coordinates" "${override_args[@]}"
fi
sequence_root="$work_root/extracted_sequences"
if want extraction || [[ "$requested_stage" == all ]]; then
  [[ $(($(wc -l < "$coordinates") - 1)) -eq 720 ]] || { echo "Extraction requires exactly 720 reviewed coordinates" >&2; exit 3; }
  run_stage 05_sequence_extraction "$coordinates" "$sequence_root" \
    bash scripts/06_phylogenomics/03_extract_gene_sequences.sh --coordinates "$coordinates" --output-root "$sequence_root"
fi
if want sequence_qc || [[ "$requested_stage" == all ]]; then
  run_stage 06_sequence_qc "$sequence_root;$coordinates" "$work_root/qc/sequence_qc.tsv" \
    "$python_cmd" scripts/06_phylogenomics/04_validate_sequences.py --sequence-root "$sequence_root" \
      --coordinates "$coordinates" --qc-dir "$work_root/qc" --tables-dir results/tables/phylogeny_10_locus
fi
if want alignment || [[ "$requested_stage" == all ]]; then
  run_stage 07_mafft_alignment "$sequence_root/by_gene" "$work_root/alignments" \
    bash scripts/06_phylogenomics/05_align_genes.sh --threads "$threads" \
      --input-root "$sequence_root/by_gene" --output-root "$work_root/alignments"
fi
if want concatenation || [[ "$requested_stage" == all ]]; then
  run_stage 08_alignment_qc_and_concatenation "$work_root/alignments" "$work_root/concatenated;$work_root/qc" \
    "$python_cmd" scripts/06_phylogenomics/06_concatenate_alignments.py --alignment-dir "$work_root/alignments" \
      --output-dir "$work_root/concatenated" --qc-dir "$work_root/qc"
fi
if want iqtree; then
  run_stage 10_iqtree "$work_root/concatenated" "$work_root/iqtree" \
    bash scripts/06_phylogenomics/07_run_iqtree.sh \
      --alignment "$work_root/concatenated/efaecalis_72_10_locus_alignment.fasta" \
      --partitions "$work_root/concatenated/efaecalis_72_10_locus_partitions.nex" \
      --output-dir "$work_root/iqtree" --threads "$threads"
fi
if want visualisation; then
  tree="$work_root/iqtree/efaecalis_72_10_locus.treefile"
  run_stage 11_metadata_join_and_visualisation "$tree;data/metadata/curated_metadata_72_genomes.csv" "results/figures/phylogeny_10_locus" \
    Rscript scripts/06_phylogenomics/08_visualise_10_locus_tree.R --tree "$tree" \
      --metadata data/metadata/curated_metadata_72_genomes.csv \
      --output-dir results/figures/phylogeny_10_locus --tables-dir results/tables/phylogeny_10_locus
fi
echo "Finish: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Run manifest: $summary"
