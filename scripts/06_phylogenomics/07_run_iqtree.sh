#!/usr/bin/env bash
set -Eeuo pipefail
alignment=""; partitions=""; output_dir=""; threads="AUTO"; output_prefix="efaecalis_72_genomes_9_locus_observed_indels"
while (($#)); do
  case "$1" in
    --alignment) alignment=$2; shift 2 ;;
    --partitions) partitions=$2; shift 2 ;;
    --output-dir) output_dir=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    --output-prefix) output_prefix=$2; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -s "$alignment" && -s "$partitions" && -n "$output_dir" ]] || { echo "Valid --alignment, --partitions and --output-dir are required" >&2; exit 2; }
if command -v iqtree2 >/dev/null; then iqtree_cmd=iqtree2
elif command -v iqtree >/dev/null; then iqtree_cmd=iqtree
else echo "IQ-TREE 2 is required" >&2; exit 127
fi
version=$("$iqtree_cmd" --version 2>&1 | head -n 2)
grep -Eq 'IQ-TREE.*version 2|IQ-TREE multicore version 2' <<<"$version" || { echo "Detected command is not IQ-TREE 2: $version" >&2; exit 3; }
help=$("$iqtree_cmd" -h 2>&1)
for option in '-p' '-m' '-B' '-alrt' '-T'; do
  grep -q -- "$option" <<<"$help" || { echo "Installed IQ-TREE does not advertise required option $option" >&2; exit 3; }
done
mkdir -p "$output_dir"
prefix="$output_dir/$output_prefix"
for suffix in treefile contree iqtree log best_scheme.nex mldist ckp.gz; do
  [[ ! -e "${prefix}.${suffix}" ]] || { echo "Refusing to overwrite ${prefix}.${suffix}" >&2; exit 4; }
done
echo "+ $iqtree_cmd -s $alignment -p $partitions -m MFP+MERGE -B 1000 -alrt 1000 -T $threads --prefix $prefix"
"$iqtree_cmd" -s "$alignment" -p "$partitions" -m MFP+MERGE -B 1000 -alrt 1000 -T "$threads" --prefix "$prefix"
[[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]] || { echo "IQ-TREE did not produce required final files" >&2; exit 5; }
echo "IQ-TREE completed: ${prefix}.treefile"
grep -E 'Best-fit model|Log-likelihood|Number of sequences|Alignment has|UFBoot|SH-aLRT' "${prefix}.iqtree" "${prefix}.log" 2>/dev/null | tail -n 20 || true
