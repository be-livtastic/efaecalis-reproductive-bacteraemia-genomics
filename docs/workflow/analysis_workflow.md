# Analysis workflow

## Status key

- **Implemented:** represented by retained scripts and compact outputs.
- **Preliminary/archived:** previously explored but excluded pending revision.
- **Future:** planned and not yet implemented.

## 1. NCBI retrieval — implemented

The canonical 14- and 58-genome accession lists are consolidated into
`data/accession_lists/selected_72_accessions.tsv`. The retrieval script validates
group counts and downloads public assemblies with NCBI Datasets into ignored
local storage. Existing destinations cause a safe failure.

## 2. Metadata curation — implemented, under review

Public accession-linked fields are retained in
`data/metadata/curated_metadata_72_genomes.csv`. Personal paths, workbook row
numbers and the original Excel workbook are excluded. Field definitions and
missing-value treatment are documented in the data dictionary.

## 3. Genome selection — implemented

The current comparison contains 14 reproductive-associated and 58
bacteraemia-associated genomes. Accession/category assignments are the canonical
selection interface. Selection should be revalidated before each full rerun.

## 4. AMRFinderPlus — implemented

AMRFinderPlus was run using version 4.2.7 and database 2026-05-15.1. Retained
scripts preserve the existing scientific classification logic. They now use
portable paths, check required inputs and refuse to overwrite previous results
unless the user deliberately opts in.

## 5. Prokka — implemented annotation run; bulk outputs local

All 72 genomes were annotated with Prokka 1.15.6. Compact per-genome annotation
statistics and run provenance are retained. Bulk GFF, FASTA, GenBank, SQN and
log files remain in the ignored local archive and can be regenerated.

## 6. Panaroo — future

A Panaroo pangenome/core-genome workflow has not yet been added. Its parameters,
quality thresholds and outputs must be documented when the revised phylogeny
method is approved.

## 7. Optional recombination filtering — future

No recombination-filtering tool has been selected. Whether to apply filtering,
and with which tool and parameters, remains a scientific-method decision.

## 8. IQ-TREE — preliminary work archived; ten-locus workflow implemented

Earlier IQ-TREE 2.0.7 results were based on a single `rpoB` gene and are
archived. A replacement workflow now implements a defined ten-locus concatenated
housekeeping-gene analysis. Candidate annotations must pass a manual review gate
before sequence extraction, and IQ-TREE is not run automatically during initial
inspection.

### Ten-locus concatenated housekeeping-gene phylogeny

The analysis uses the seven established *E. faecalis* MLST loci (`gdh`, `gyd`,
`pstS`, `gki`, `aroE`, `xpt`, and `yqiL`) plus three conserved markers (`pyrC`,
`groEL`, and `recA`). Complete annotated CDSs are aligned independently and
concatenated in that fixed order. MLST primer and reference evidence supports
annotation disambiguation but does not define the extracted sequence span.

This multilocus analysis is more informative than the archived single-gene
`rpoB` exploration, but it is not equivalent to a pangenome-derived core-genome
phylogeny. It provides a defined, reproducible view that should be interpreted
alongside curated metadata and AMR findings. Phylogenetic proximity does not
demonstrate pathogenicity, transmission, or epidemiological linkage.

## 9. Statistical analysis — partial/future

Current AMR summaries are descriptive. Confirmatory comparative tests, covariate
handling and multiple-testing policy require explicit scientific approval.

## 10. R visualisation — implemented for AMR; exploratory overview retained

Selected AMR figures and an exploratory dataset-overview script are retained.
Figures should be regenerated from canonical tables and labelled according to
the completion status of their underlying analysis.
