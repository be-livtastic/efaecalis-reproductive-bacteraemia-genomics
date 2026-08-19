#!/usr/bin/env bash
set -uo pipefail

# Default prefix tracks the current handoff commit and can be overridden.
expected_prefix="${1:-ffa0b93}"
alignment_path="analysis/core_genome/panaroo_strict_core95/core_gene_alignment.aln"
declare -a blockers=()

record_blocker() {
  blockers+=("$1")
}

run_cmd() {
  local output
  if output="$("$@" 2>&1)"; then
    printf '%s\n' "$output"
    RUN_CMD_OUTPUT="$output"
    return 0
  fi
  printf '%s\n' "$output"
  RUN_CMD_OUTPUT="$output"
  return 1
}

echo "${CODESPACES:-}"
hostname
pwd
whoami
echo

if ! run_cmd git --version; then
  record_blocker "Git is not available in PATH."
fi
if ! run_cmd command -v git; then
  record_blocker "Git executable path could not be resolved."
fi
echo

head_sha=""
origin_sha=""
if run_cmd git rev-parse HEAD; then
  head_sha="$RUN_CMD_OUTPUT"
else
  record_blocker "Unable to resolve local HEAD commit."
fi
if run_cmd git rev-parse origin/main; then
  origin_sha="$RUN_CMD_OUTPUT"
else
  record_blocker "Unable to resolve origin/main commit."
fi
run_cmd git status --short --branch || true
run_cmd git remote -v || true
echo

if [[ -n "$head_sha" && -n "$origin_sha" && "$head_sha" != "$origin_sha" ]]; then
  record_blocker "HEAD (${head_sha}) does not match origin/main (${origin_sha})."
fi
if [[ -n "$head_sha" && "$head_sha" != "${expected_prefix}"* ]]; then
  record_blocker "HEAD (${head_sha}) does not start with expected prefix ${expected_prefix}."
fi

current_branch="$(git branch --show-current 2>/dev/null || true)"
if [[ "$current_branch" != "main" ]]; then
  record_blocker "Current branch is '${current_branch:-unknown}', expected 'main'."
fi

if run_cmd nproc; then
  if [[ "$RUN_CMD_OUTPUT" != "4" ]]; then
    record_blocker "CPU count is ${RUN_CMD_OUTPUT}, expected 4."
  fi
fi
run_cmd awk '/MemTotal|MemAvailable|SwapTotal|SwapFree/ {print}' /proc/meminfo || true
run_cmd df -h || true
echo

mem_total_kib="$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
if [[ -n "$mem_total_kib" ]]; then
  if (( mem_total_kib < 15000000 || mem_total_kib > 18000000 )); then
    record_blocker "MemTotal is ${mem_total_kib} KiB, expected approximately 16 GB."
  fi
fi

if ! run_cmd command -v micromamba; then
  record_blocker "micromamba is not available in PATH."
else
  if run_cmd micromamba env list; then
    if ! grep -Eq '^efaecalis_core_genome([[:space:]]|$)' <<<"$RUN_CMD_OUTPUT"; then
      record_blocker "Environment 'efaecalis_core_genome' was not found in micromamba env list."
    fi
  else
    record_blocker "Failed to list micromamba environments."
  fi
  echo

  if run_cmd micromamba run -n efaecalis_core_genome python --version; then
    [[ "$RUN_CMD_OUTPUT" == "Python 3.11.15" ]] || record_blocker "Python version is '${RUN_CMD_OUTPUT}', expected 'Python 3.11.15'."
  else
    record_blocker "Unable to run Python in efaecalis_core_genome."
  fi

  if run_cmd micromamba run -n efaecalis_core_genome iqtree --version; then
    [[ "$RUN_CMD_OUTPUT" == *"3.1.2"* ]] || record_blocker "IQ-TREE version output does not include 3.1.2."
  else
    record_blocker "Unable to run IQ-TREE in efaecalis_core_genome."
  fi

  if run_cmd micromamba run -n efaecalis_core_genome mafft --version; then
    [[ "$RUN_CMD_OUTPUT" == *"v7.526"* || "$RUN_CMD_OUTPUT" == *"7.526"* ]] || record_blocker "MAFFT version output does not include 7.526."
  else
    record_blocker "Unable to run MAFFT in efaecalis_core_genome."
  fi

  if run_cmd micromamba run -n efaecalis_core_genome panaroo --version; then
    [[ "$RUN_CMD_OUTPUT" == "1.8.0" ]] || record_blocker "Panaroo version is '${RUN_CMD_OUTPUT}', expected 1.8.0."
  else
    record_blocker "Unable to run Panaroo in efaecalis_core_genome."
  fi

  if run_cmd micromamba run -n efaecalis_core_genome snp-dists -v; then
    [[ "$RUN_CMD_OUTPUT" == *"1.2.0"* ]] || record_blocker "snp-dists version output does not include 1.2.0."
  else
    record_blocker "Unable to run snp-dists in efaecalis_core_genome."
  fi

  if run_cmd micromamba run -n efaecalis_core_genome python - <<'PY'
import Bio, pandas, numpy
print("Biopython", Bio.__version__)
print("pandas", pandas.__version__)
print("numpy", numpy.__version__)
PY
  then
    grep -q '^Biopython 1.87$' <<<"$RUN_CMD_OUTPUT" || record_blocker "Biopython version output does not match 1.87."
    grep -q '^pandas 3.0.5$' <<<"$RUN_CMD_OUTPUT" || record_blocker "pandas version output does not match 3.0.5."
    grep -q '^numpy 2.3.5$' <<<"$RUN_CMD_OUTPUT" || record_blocker "numpy version output does not match 2.3.5."
  else
    record_blocker "Unable to verify Biopython/pandas/numpy versions in efaecalis_core_genome."
  fi
fi
echo

if [[ -f "$alignment_path" ]]; then
  echo "ALIGNMENT_PRESENT"
  record_blocker "Alignment already exists at ${alignment_path}; expected ALIGNMENT_ABSENT before transfer."
else
  echo "ALIGNMENT_ABSENT"
fi

if ((${#blockers[@]} == 0)); then
  echo "READY_FOR_ALIGNMENT_TRANSFER"
else
  echo "STOPPED_PRETRANSFER"
  echo "${blockers[0]}"
fi
