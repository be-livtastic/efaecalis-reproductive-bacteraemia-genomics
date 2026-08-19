# Core-genome pipeline

This staged workflow complements, and never regenerates, the validated nine-locus phylogeny. It uses exactly 72 canonical Prokka GFF files, preserves accession version suffixes, runs Panaroo with strict cleaning and a 95% core threshold, validates the resulting alignment, and only then permits IQ-TREE and distance analysis.

## Environment

The dependency preflight found Panaroo 1.8.0 incompatible with the validated environment's Python 3.12 pin. The isolated environment is specified reproducibly as:

```bash
mamba env create -f environment/core_genome.yml
export EFAECALIS_CORE_ENV=efaecalis_core_genome
```

The existing environment is not modified. Four threads are used because the current machine has eight CPUs but only about 3.7 GiB RAM.

The 2026-08-17 implementation run created this environment and confirmed the exact installed versions recorded in `environment/core_genome_software_versions.tsv`. The solver decision is documented in `environment/core_genome_environment_decision.md`.

## Stages and checkpoints

```bash
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage audit
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage panaroo
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage panaroo-qc
```

Stop after `panaroo-qc`. Inspect `results/tables/core_genome/panaroo_checkpoint_72.csv` and `panaroo_genome_representation_72.csv`. Continue only if all 72 accessions are represented and the alignment is valid:

```bash
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage alignment-qc
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage iqtree
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage distances
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage pangenome-summary
bash scripts/07_core_genome/run_core_genome_pipeline.sh --stage comparison
```

Every stage refuses existing outputs unless `--overwrite` is supplied. Computational outputs are stored under `analysis/core_genome`, final tables under `results/tables/core_genome`, figures under `results/figures/core_genome`, and runtime logs under `analysis/logs/core_genome`.

Core families are present in at least 69 of 72 genomes; families present in fewer genomes are accessory. Additional frequency categories are reported only when Panaroo supplies them directly. Distances are called **pairwise SNP distances derived from the concatenated core-gene alignment**, never whole-genome SNP distances. Default `snp-dists` ignores gaps and ambiguous bases.

All analytical trees remain unrooted. Midpoint rooting is display-only. No project-defined clades are assigned, no transmission is inferred, and Gubbins is not run automatically.

## Rebuilt-Codespace pre-transfer gate

After pushing commit `ffa0b93` (or its full SHA) to `main` and rebuilding Codespace `expert-guide-5gx777p954pq37p6v`, run:

```bash
bash scripts/07_core_genome/validate_pretransfer_environment.sh
```

This validates the exact pre-transfer constraints in one pass: `main` branch state, `HEAD`/`origin/main` agreement, expected core tool versions, expected package versions, and that `analysis/core_genome/panaroo_strict_core95/core_gene_alignment.aln` is absent before transfer.

Required terminal token:

- `READY_FOR_ALIGNMENT_TRANSFER` when every gate passes.
- `STOPPED_PRETRANSFER` plus the first blocker when any gate fails.

Do not transfer the alignment, run Panaroo, or start IQ-TREE until this gate returns `READY_FOR_ALIGNMENT_TRANSFER`.

Helper tests cover versioned accession parsing, missing-genome rejection, 1.5-IQR flag retention, tied nearest neighbours and SNP-matrix symmetry:

```bash
python -m unittest scripts/07_core_genome/test_core_genome_analysis.py
```
