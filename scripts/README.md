# Script stages

Scripts are grouped by intended execution order. Empty stages are deliberate:
their scientific method has not yet been implemented or approved.

1. `00_setup`: setup helpers, if later required.
2. `01_genome_retrieval`: public genome download and accession-safe renaming.
3. `02_metadata_curation`: future reproducible metadata derivation.
4. `03_amr_analysis`: retained AMRFinderPlus analysis and isolate summary.
5. `04_annotation`: future portable Prokka runner; completed run summaries are
   under `results/tables/annotation_qc`.
6. `05_pangenome`: future Panaroo stage.
7. `06_phylogenomics`: deliberately empty while the workflow is revised.
8. `07_statistics`: future inferential analysis.
9. `08_visualisation`: exploratory dataset overview and retained plotting code.

Scripts must use project-relative/configurable paths and refuse silent
overwriting. Linux and WSL filesystems are case-sensitive.

