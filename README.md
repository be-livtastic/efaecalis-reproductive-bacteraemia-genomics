# Comparative Genomics of Reproductive- and Bacteraemia-Associated *Enterococcus faecalis*

> **Work in progress:** this MSc project is under active development. Results,
> scripts and interpretations may change following quality control and review.

## Project overview

This repository supports a comparative genomics investigation of publicly
available *Enterococcus faecalis* genomes associated with reproductive samples
and bacteraemia. It contains reproducible accession lists, curated non-sensitive
metadata, analysis scripts, compact quality-control records and selected
outputs. Public genome assemblies and large generated intermediates are
retrieved or regenerated locally and are not stored in Git.

## Research rationale and aim

*E. faecalis* occupies diverse ecological and clinical contexts. Comparative
genomics can help describe similarities and differences between genome sets,
but genomic associations alone do not establish pathogenicity, transmission or
phenotype.

The project aims to compare genomic characteristics of the two selected public
datasets, with particular attention to antimicrobial-resistance determinants
and, in later work, core-genome population structure.

## Dataset

- 14 reproductive-associated genomes
- 58 bacteraemia-associated genomes
- 72 genomes in the combined dataset

The groups are unequal in size. Accessions and category assignments are defined
in `data/accession_lists/selected_72_accessions.tsv`.

## Repository structure

- `config/`: portable example configuration.
- `data/accession_lists/`: public assembly accessions and dataset assignments.
- `data/metadata/`: curated metadata and its data dictionary.
- `environment/`: pinned software environment and recorded versions.
- `scripts/`: ordered retrieval, analysis and visualisation scripts.
- `results/`: selected compact tables and figures.
- `docs/`: workflow, decisions and troubleshooting documentation.
- `manuscript/`: reserved for an appropriate future manuscript draft.

`data/raw/`, `data/interim/`, current phylogeny work and `local_archive/` are
local-only.

## Analysis workflow

The intended workflow is:

NCBI retrieval → metadata curation → genome selection → AMRFinderPlus → Prokka
→ Panaroo → optional recombination filtering → IQ-TREE → statistical analysis
→ R visualisation.

AMRFinderPlus analysis and preliminary Prokka annotation have been run. Panaroo,
recombination filtering and final core-genome phylogenomics are not yet present.
Earlier single-gene `rpoB` work is intentionally excluded while the phylogeny
workflow is revised. See `docs/workflow/analysis_workflow.md`.

## Software

The reproducible environment is described in `environment/environment.yml`.
Recorded run versions are listed in `environment/software_versions.tsv`.
AMRFinderPlus database provenance is in
`environment/amrfinderplus_version.tsv`.

## Reproduction

1. Create the Conda/Mamba environment:
   `mamba env create -f environment/environment.yml`
2. Activate it: `mamba activate efaecalis_genomics`
3. Retrieve assemblies with
   `bash scripts/01_genome_retrieval/download_ncbi_genomes.sh`
4. Run scripts from their documented stage in numerical order.
5. Review input and output paths in `config/paths.example.yml`.

Scripts refuse to replace existing outputs by default. Personal absolute paths
are not required.

## Expected outputs

Expected repository-scale outputs include curated metadata, accession manifests,
AMRFinderPlus summary tables, Prokka annotation summaries, quality-control
tables and selected figures. Bulk genomes and pipeline intermediates remain
local and reproducible.

## Ethical and data considerations

The project uses public sequence accessions and public associated metadata.
Signed ethics forms, risk assessments, administrative records, signatures,
student numbers and personal contact information are not part of this
repository. The original working Excel workbook is retained locally and is
explicitly excluded.

## Limitations

Interpretation must account for unequal sample sizes, study and geographic
bias, incomplete public metadata, and the inability to infer antimicrobial
phenotype directly from genotype. Some accessions may represent repeated or
epidemiologically linked isolates. Reproductive-associated isolates should not
be assumed to be pathogenic, and dataset-level associations should not be
interpreted as causal.

## Citation status

This work has not yet reached a final citation or publication state. Software,
database and public accession citations will be completed as the dissertation
and analyses mature.

## Author

Olivia Williams  
MSc Biotechnology, University of Chester

