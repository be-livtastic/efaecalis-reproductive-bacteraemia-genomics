#!/usr/bin/env bash
set -euo pipefail

environment_name="efaecalis_core_genome"
repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
environment_file="${repository_root}/environment/core_genome.yml"
mamba_root_prefix="${MAMBA_ROOT_PREFIX:-/opt/conda}"

export MAMBA_ROOT_PREFIX="${mamba_root_prefix}"

current_home_mode="$(stat -c '%a' "${HOME}")"
if [ "${current_home_mode}" != "755" ]; then
  chmod 755 "${HOME}"
fi

if ! command -v micromamba >/dev/null 2>&1; then
  printf 'micromamba is required but was not found in PATH.\n' >&2
  exit 1
fi

if ! micromamba env list | awk 'NR > 1 {print $1}' | grep -Fxq "${environment_name}"; then
  micromamba create --yes --file "${environment_file}"
else
  printf 'Environment %s already exists; leaving it unchanged.\n' "${environment_name}"
fi

shell_init_file="${HOME}/.bashrc"
shell_hook='eval "$(micromamba shell hook --shell=bash --root-prefix="${MAMBA_ROOT_PREFIX:-/opt/conda}")"'
activation_block="$(printf '%s\n%s\n%s' '# efaecalis core-genome environment' "${shell_hook}" 'micromamba activate efaecalis_core_genome')"

if ! grep -Fq '# efaecalis core-genome environment' "${shell_init_file}" 2>/dev/null; then
  printf '\n%s\n' "${activation_block}" >> "${shell_init_file}"
fi

micromamba run --name "${environment_name}" iqtree --version >/dev/null
micromamba run --name "${environment_name}" mafft --version >/dev/null
micromamba run --name "${environment_name}" panaroo --version >/dev/null
micromamba run --name "${environment_name}" snp-dists --version >/dev/null

printf 'Validated tools in %s:\n' "${environment_name}"
micromamba run --name "${environment_name}" bash -c 'command -v iqtree mafft panaroo snp-dists'