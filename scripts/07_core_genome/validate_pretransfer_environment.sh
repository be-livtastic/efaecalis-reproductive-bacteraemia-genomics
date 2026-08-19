#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
core_env=${EFAECALIS_CORE_ENV:-efaecalis_core_genome}
destination_dir="$project_root/analysis/core_genome/panaroo_strict_core95"
alignment_path="$destination_dir/core_gene_alignment.aln"

expected_python="3.11.15"
expected_iqtree="3.1.2"
expected_mafft="7.526"
expected_panaroo="1.8.0"
expected_snp_dists="1.2.0"
expected_biopython="1.87"
expected_pandas="3.0.5"
expected_numpy="2.3.5"

declare -a blockers=()

add_blocker() {
  blockers+=("$1")
}

if ! command -v git >/dev/null 2>&1; then
  add_blocker "git is not installed or not in PATH"
fi

if ! git -C "$project_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  add_blocker "workspace is not a readable Git repository: $project_root"
else
  current_branch=$(git -C "$project_root" branch --show-current 2>/dev/null || true)
  if [[ "$current_branch" != "main" ]]; then
    add_blocker "current branch is '$current_branch' (expected 'main')"
  fi

  if ! git -C "$project_root" status --short >/dev/null 2>&1; then
    add_blocker "git status is not readable for this workspace"
  fi

  if ! git -C "$project_root" rev-parse --verify origin/main^{commit} >/dev/null 2>&1; then
    add_blocker "origin/main is not resolvable to a commit"
  fi
fi

cpu_count=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 0)
if ! [[ "$cpu_count" =~ ^[0-9]+$ ]] || (( cpu_count < 4 )); then
  add_blocker "CPU check failed: found '$cpu_count' CPUs (requires >= 4)"
fi

mem_total_kb=$(awk '/MemTotal:/ {print $2}' /proc/meminfo 2>/dev/null || echo 0)
mem_available_kb=$(awk '/MemAvailable:/ {print $2}' /proc/meminfo 2>/dev/null || echo 0)
if ! [[ "$mem_total_kb" =~ ^[0-9]+$ ]] || (( mem_total_kb < 14680064 )); then
  add_blocker "RAM check failed: total RAM is below ~14 GiB (expected approximately 16 GiB)"
fi
if ! [[ "$mem_available_kb" =~ ^[0-9]+$ ]] || (( mem_available_kb <= 0 )); then
  add_blocker "RAM check failed: available RAM could not be determined"
fi

if ! command -v micromamba >/dev/null 2>&1; then
  add_blocker "micromamba is not installed or not in PATH"
else
  if ! micromamba env list | awk 'NR > 1 {print $1}' | grep -Fxq "$core_env"; then
    add_blocker "Conda environment '$core_env' does not exist"
  else
    python_version=$(micromamba run -n "$core_env" python --version | awk '{print $2}')
    [[ "$python_version" == "$expected_python" ]] || add_blocker "Python version is '$python_version' (expected '$expected_python')"

    iqtree_version=$(micromamba run -n "$core_env" iqtree --version | awk '/IQ-TREE version/ {print $3; exit}')
    [[ "$iqtree_version" == "$expected_iqtree" ]] || add_blocker "IQ-TREE version is '$iqtree_version' (expected '$expected_iqtree')"

    mafft_version=$(micromamba run -n "$core_env" mafft --version 2>&1 | awk 'NR==1 {gsub(/^v/, "", $1); print $1}')
    [[ "$mafft_version" == "$expected_mafft" ]] || add_blocker "MAFFT version is '$mafft_version' (expected '$expected_mafft')"

    panaroo_version=$(micromamba run -n "$core_env" panaroo --version | awk '{print $2}')
    [[ "$panaroo_version" == "$expected_panaroo" ]] || add_blocker "Panaroo version is '$panaroo_version' (expected '$expected_panaroo')"

    snp_dists_version=$(micromamba run -n "$core_env" snp-dists -v | awk '{print $2}')
    [[ "$snp_dists_version" == "$expected_snp_dists" ]] || add_blocker "snp-dists version is '$snp_dists_version' (expected '$expected_snp_dists')"

    readarray -t py_lib_versions < <(micromamba run -n "$core_env" python - <<'PY'
import Bio
import pandas
import numpy
print(f"Bio={Bio.__version__}")
print(f"pandas={pandas.__version__}")
print(f"numpy={numpy.__version__}")
PY
)
    biopython_version=${py_lib_versions[0]#Bio=}
    pandas_version=${py_lib_versions[1]#pandas=}
    numpy_version=${py_lib_versions[2]#numpy=}
    [[ "$biopython_version" == "$expected_biopython" ]] || add_blocker "Biopython version is '$biopython_version' (expected '$expected_biopython')"
    [[ "$pandas_version" == "$expected_pandas" ]] || add_blocker "pandas version is '$pandas_version' (expected '$expected_pandas')"
    [[ "$numpy_version" == "$expected_numpy" ]] || add_blocker "numpy version is '$numpy_version' (expected '$expected_numpy')"
  fi
fi

if [[ ! -d "$destination_dir" ]]; then
  add_blocker "destination directory is missing: $destination_dir"
fi

if [[ -e "$alignment_path" ]]; then
  add_blocker "alignment already exists and must not be transferred: $alignment_path"
fi

if ((${#blockers[@]} > 0)); then
  printf 'Pre-transfer validation failed:\n'
  for blocker in "${blockers[@]}"; do
    printf 'BLOCKER: %s\n' "$blocker"
  done
  printf 'STOPPED_PRETRANSFER\n'
  exit 1
fi

printf 'READY_FOR_ALIGNMENT_TRANSFER\n'
