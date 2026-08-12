# Script stages

Scripts are grouped by intended execution order. A stage directory exists only
when it contains an implemented script.

1. `00_setup`: create the local writable directory skeleton.
2. `01_genome_retrieval`: public genome download and accession-safe renaming.
3. `02_metadata_curation`: validate identifiers and generate analysis metadata
   plus manuscript Supplementary Table S1.
4. `03_amr_analysis`: retained AMRFinderPlus analysis and isolate summary.
5. `04_annotation`: portable Prokka runner; completed run summaries are under
   `results/tables/annotation_qc`.
6. `06_phylogenomics`: validated nine-locus extraction, QC, inference, PCA and
   tree visualisation.
7. `08_visualisation`: descriptive dataset overview plotting.

Stages `05_pangenome` and `07_statistics` are planned but do not have approved
methods, so empty placeholder directories are not retained.

Scripts must use project-relative/configurable paths and refuse silent
overwriting. Linux and WSL filesystems are case-sensitive.
