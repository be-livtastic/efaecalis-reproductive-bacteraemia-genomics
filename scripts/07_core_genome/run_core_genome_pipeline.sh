#!/usr/bin/env bash
set -euo pipefail

# Resolve the repository and configurable Conda environment before dispatching a stage.
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=${EFAECALIS_PROJECT_ROOT:-$(cd "$script_dir/../.." && pwd)}
core_env=${EFAECALIS_CORE_ENV:-efaecalis_phylogeny}
stage=""
overwrite=0

usage() {
  echo "Usage: $0 --stage audit|panaroo|panaroo-qc|alignment-qc|iqtree|distances|pangenome-summary|comparison [--overwrite]" >&2
}

# Parse only explicit stages so the mandatory checkpoints cannot be bypassed by an all-in-one command.
while (($#)); do
  case "$1" in
    --stage) stage=${2:-}; shift 2 ;;
    --overwrite) overwrite=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -n "$stage" ]] || { usage; exit 2; }

mamba_bin=${MAMBA_EXE:-/home/belivtastic/miniforge3/bin/mamba}
[[ -x "$mamba_bin" ]] || { echo "ERROR: Mamba executable not found: $mamba_bin" >&2; exit 127; }
env_prefix=${EFAECALIS_CORE_ENV_PREFIX:-/home/belivtastic/miniforge3/envs/$core_env}
main_prefix=${EFAECALIS_MAIN_ENV_PREFIX:-/home/belivtastic/miniforge3/envs/efaecalis_phylogeny}
py=("$env_prefix/bin/python" "$script_dir/core_genome_analysis.py")
rscript=("$main_prefix/bin/Rscript")
[[ -x "${py[0]}" ]] || { echo "ERROR: Core environment Python is unavailable: ${py[0]}" >&2; exit 127; }
export PATH="$env_prefix/bin:$PATH"
overwrite_arg=()
((overwrite)) && overwrite_arg=(--overwrite)
log_dir="$project_root/analysis/logs/core_genome"
mkdir -p "$log_dir"

# Record a stage-specific command and elapsed time in a durable log.
run_logged() {
  local name=$1; shift
  local log="$log_dir/${name}.log"
  if [[ -e "$log" && $overwrite -eq 0 ]]; then
    echo "ERROR: Refusing to overwrite log: $log" >&2; exit 3
  fi
  local started ended status
  started=$(date +%s)
  {
    echo "Started_UTC: $(date -u +%FT%TZ)"
    printf 'Command:'; printf ' %q' "$@"; printf '\n'
    set +e
    "$@"
    status=$?
    set -e
    ended=$(date +%s)
    echo "Exit_status: $status"
    echo "Elapsed_seconds: $((ended-started))"
    echo "Finished_UTC: $(date -u +%FT%TZ)"
    exit "$status"
  } 2>&1 | tee "$log"
}

cd "$project_root"
case "$stage" in
  audit)
    run_logged A0_input_audit "${py[@]}" audit "${overwrite_arg[@]}"
    ;;
  panaroo)
    input_dir="$project_root/data/processed/core_genome/input_gffs_72"
    output_dir="$project_root/analysis/core_genome/panaroo_strict_core95"
    [[ $(find "$input_dir" -maxdepth 1 -type f -name '*.gff' | wc -l) -eq 72 ]] || { echo "ERROR: Run and pass audit first" >&2; exit 4; }
    # Create only Panaroo's parent because Panaroo itself requires its output directory not to exist.
    mkdir -p "$(dirname "$output_dir")"
    if [[ -e "$output_dir" && $overwrite -eq 0 ]]; then echo "ERROR: Refusing to overwrite $output_dir" >&2; exit 3; fi
    if [[ -e "$output_dir" ]]; then mv "$output_dir" "${output_dir}.backup.$(date -u +%Y%m%dT%H%M%SZ)"; fi
    mapfile -t gffs < <(find "$input_dir" -maxdepth 1 -type f -name '*.gff' | sort)
    [[ -x "$env_prefix/bin/panaroo" ]] || { echo "ERROR: Panaroo is unavailable in $env_prefix" >&2; exit 127; }
    run_logged A1_panaroo "$env_prefix/bin/panaroo" -i "${gffs[@]}" -o "$output_dir" --clean-mode strict -a core --aligner mafft --core_threshold 0.95 -t 4
    echo "CHECKPOINT: Panaroo finished. Run --stage panaroo-qc and inspect its outputs before alignment-qc or IQ-TREE." 
    ;;
  panaroo-qc)
    run_logged A1_panaroo_qc "${py[@]}" panaroo-qc "${overwrite_arg[@]}"
    echo "CHECKPOINT: Inspect panaroo_checkpoint_72.csv and panaroo_genome_representation_72.csv before proceeding."
    ;;
  alignment-qc)
    run_logged A2_alignment_qc "${py[@]}" alignment-qc "${overwrite_arg[@]}"
    ;;
  iqtree)
    alignment="$project_root/analysis/core_genome/panaroo_strict_core95/core_gene_alignment.aln"
    prefix="$project_root/analysis/core_genome/iqtree/core_gene_alignment"
    [[ -s "$project_root/results/tables/core_genome/core_alignment_qc_72.csv" ]] || { echo "ERROR: Alignment QC output is missing" >&2; exit 4; }
    mkdir -p "$(dirname "$prefix")"
    if compgen -G "${prefix}.*" >/dev/null && [[ $overwrite -eq 0 ]]; then echo "ERROR: Refusing to overwrite IQ-TREE outputs" >&2; exit 3; fi
    iq_args=(-s "$alignment" -m MFP -B 1000 --alrt 1000 -T 4 --prefix "$prefix")
    ((overwrite)) && iq_args+=(-redo)
    [[ -x "$env_prefix/bin/iqtree" ]] || { echo "ERROR: IQ-TREE is unavailable in $env_prefix" >&2; exit 127; }
    run_logged A3_iqtree "$env_prefix/bin/iqtree" "${iq_args[@]}"
    ;;
  distances)
    alignment="$project_root/analysis/core_genome/panaroo_strict_core95/core_gene_alignment.aln"
    output_dir="$project_root/analysis/core_genome/distances"
    matrix="$output_dir/core_alignment_snp_distances_72.tsv"
    mkdir -p "$output_dir"
    if [[ -e "$matrix" && $overwrite -eq 0 ]]; then echo "ERROR: Refusing to overwrite $matrix" >&2; exit 3; fi
    [[ -x "$env_prefix/bin/snp-dists" ]] || { echo "ERROR: snp-dists is unavailable in $env_prefix" >&2; exit 127; }
    run_logged A4_snp_dists bash -c '"$1" -j 4 -b "$2" > "$3"' _ "$env_prefix/bin/snp-dists" "$alignment" "$matrix"
    run_logged A4_nearest_neighbours "${py[@]}" distance-results "${overwrite_arg[@]}"
    ;;
  pangenome-summary)
    run_logged A5_pangenome_summary "${py[@]}" pangenome-summary "${overwrite_arg[@]}"
    ;;
  comparison)
    run_logged A6_tree_comparison "${rscript[@]}" "$script_dir/compare_core_phylogenies.R" "${overwrite_arg[@]}"
    ;;
  *) echo "ERROR: Unknown stage: $stage" >&2; usage; exit 2 ;;
esac
