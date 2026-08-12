#!/usr/bin/env bash
set -Eeuo pipefail
# --- Resolve project tools and canonical inputs ---
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd); cd "$repo_root"
python_cmd=${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}; [[ -x "${python_cmd:-}" ]] || python_cmd=python3
rscript_cmd=${CONDA_PREFIX:+$CONDA_PREFIX/bin/Rscript}; [[ -x "${rscript_cmd:-}" ]] || rscript_cmd=Rscript
alignment=analysis/phylogenomics/72_genomes_9_locus_observed_indels/concatenated/efaecalis_72_genomes_9_locus_observed_indels_alignment.fasta
partitions=analysis/phylogenomics/72_genomes_9_locus_observed_indels/concatenated/partition_coordinates.tsv
tables=results/tables/phylogeny_72_genomes_9_locus_observed_indels/pca
figures=results/figures/phylogeny_72_genomes_9_locus_observed_indels/pca
# --- Calculate, plot and fingerprint the completed PCA ---
"$python_cmd" scripts/06_phylogenomics/10_run_alignment_pca.py --alignment "$alignment" --partitions "$partitions" --metadata data/metadata/curated_metadata_72_genomes.csv --output-dir "$tables" --min-allele-count 2 --components 10
"$rscript_cmd" scripts/06_phylogenomics/10_plot_alignment_pca.R --scores "$tables/pca_scores.tsv" --variance "$tables/pca_explained_variance.tsv" --output-dir "$figures"
sha256sum "$alignment" > "$tables/input_alignment.sha256"
printf 'status\tCOMPLETE\nmethod\tSVD of standardised observed-allele features\n' > "$tables/PCA_COMPLETE.tsv"
