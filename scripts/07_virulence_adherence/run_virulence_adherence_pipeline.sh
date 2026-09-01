#!/usr/bin/env bash
set -euo pipefail

# Resolve repository paths and dispatch one independently reviewable virulence stage.
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=${EFAECALIS_PROJECT_ROOT:-$(cd "$script_dir/../.." && pwd)}
analysis_env=${EFAECALIS_VIRULENCE_ENV:-efaecalis_phylogeny}
stage=""
overwrite=0

usage() {
  echo "Usage: $0 --stage reference-audit|pilot|screen|aggregation-review|resolve|analyse|core-integration [--overwrite]" >&2
}

# Require an explicit stage so reference and pilot checkpoints remain enforceable.
while (($#)); do
  case "$1" in
    --stage) stage=${2:-}; shift 2 ;;
    --overwrite) overwrite=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -n "$stage" ]] || { usage; exit 2; }

mamba_bin=${MAMBA_EXE:-}
if [[ -z "$mamba_bin" ]]; then
  mamba_bin=$(command -v mamba || command -v micromamba || true)
fi
[[ -n "$mamba_bin" && -x "$mamba_bin" ]] || {
  echo "ERROR: Mamba was not found. Activate Miniforge or set MAMBA_EXE." >&2
  exit 127
}
mamba_root=${MAMBA_ROOT_PREFIX:-$(cd "$(dirname "$mamba_bin")/.." && pwd)}
env_prefix=${EFAECALIS_VIRULENCE_ENV_PREFIX:-$mamba_root/envs/$analysis_env}
main_prefix=${EFAECALIS_MAIN_ENV_PREFIX:-$mamba_root/envs/efaecalis_phylogeny}
py=("$env_prefix/bin/python" "$script_dir/virulence_screen.py")
rscript=("$main_prefix/bin/Rscript")
[[ -x "${py[0]}" ]] || { echo "ERROR: Virulence environment Python is unavailable: ${py[0]}" >&2; exit 127; }
export PATH="$env_prefix/bin:$PATH"
overwrite_arg=()
((overwrite)) && overwrite_arg=(--overwrite)
log_dir="$project_root/analysis/logs/virulence_adherence"
mkdir -p "$log_dir"

# Record command, status and runtime for each independent stage.
run_logged() {
  local name=$1; shift
  local log="$log_dir/${name}.log"
  if [[ -e "$log" && $overwrite -eq 0 ]]; then echo "ERROR: Refusing to overwrite log: $log" >&2; exit 3; fi
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
  reference-audit)
    run_logged B1_reference_audit "${py[@]}" reference-audit "${overwrite_arg[@]}"
    echo "CHECKPOINT: Inspect virulence_reference_audit.csv and the pinned mapping before running the pilot."
    ;;
  pilot)
    [[ -s "$project_root/results/tables/virulence_adherence/virulence_reference_audit.csv" ]] || { echo "ERROR: Run and pass reference-audit first" >&2; exit 4; }
    run_logged B2_pilot "${py[@]}" pilot "${overwrite_arg[@]}"
    echo "CHECKPOINT: Inspect pilot HSPs, consolidated candidates and per-target status QC before screening all 72 genomes."
    ;;
  screen)
    [[ -s "$project_root/results/tables/virulence_adherence/virulence_pilot_status_qc.csv" ]] || { echo "ERROR: Run and approve the pilot first" >&2; exit 4; }
    run_logged B3_screen_72 "${py[@]}" screen "${overwrite_arg[@]}"
    ;;
  aggregation-review)
    [[ -s "$project_root/data/processed/virulence_adherence/virulence_target_status_by_genome_72.csv" ]] || { echo "ERROR: Full 72-genome status table is missing" >&2; exit 4; }
    run_logged B3_aggregation_substance_review "${py[@]}" aggregation-review "${overwrite_arg[@]}"
    echo "CHECKPOINT: Review competitive asa1/family assignments before resolving statuses or running prevalence/Fisher analysis."
    ;;
  resolve)
    [[ -s "$project_root/results/tables/virulence_adherence/aggregation_substance_genome_review_recommendations_72.csv" ]] || { echo "ERROR: Aggregation-substance review is missing" >&2; exit 4; }
    run_logged B4_resolved_statuses "${py[@]}" resolve "${overwrite_arg[@]}"
    ;;
  analyse)
    [[ -s "$project_root/data/processed/virulence_adherence/virulence_resolved_target_status_by_genome_72.csv" ]] || { echo "ERROR: Approved resolved 72-genome status table is missing" >&2; exit 4; }
    run_logged B4_B6_analysis "${rscript[@]}" "$script_dir/analyse_virulence_adherence.R" "${overwrite_arg[@]}"
    ;;
  core-integration)
    echo "ERROR: Core integration is intentionally unavailable until Pipeline A passes and its 72-tip tree is reviewed." >&2
    exit 4
    ;;
  *) echo "ERROR: Unknown stage: $stage" >&2; usage; exit 2 ;;
esac
