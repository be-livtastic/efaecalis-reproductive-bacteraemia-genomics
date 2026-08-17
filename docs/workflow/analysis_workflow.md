# Analysis workflow

## Current dataset and provenance

The canonical comparison contains 72 complete human *Enterococcus faecalis*
genomes: 14 reproductive-associated and 58 bacteraemia-associated. Assembly
accessions, including version suffixes such as `.1`, and source assignments are
defined in `data/accession_lists/selected_72_accessions.tsv`. Curated public
metadata and strain names are in
`data/metadata/curated_metadata_72_genomes.csv`. Every analytical join must
match all 72 accessions exactly; genomes are never dropped silently.

Large assemblies, annotations and raw tool outputs are local/regenerable.
Compact audit records, scripts, decisions, tables and selected figures are
retained in the repository.

## 1. Retrieval, metadata and annotation — completed

The NCBI retrieval stage validates the 14/58 accession manifest before writing
assemblies to local storage. Metadata curation checks one-to-one
assembly–BioSample–strain mappings and writes the canonical 72-row metadata
table. Prokka 1.15.6 annotations exist for all 72 genomes and are the frozen
annotation input for the multilocus and Panaroo workflows.

## 2. Comparative AMRFinderPlus analysis — completed

The authoritative pipeline is
`scripts/03_amr_analysis/amrfinder_comparative_analysis_72_genomes.R`. It reads
72 per-genome AMRFinderPlus TSVs, requires the validated 23-column schema, and
preserves all 567 raw records (523 AMR and 44 STRESS) in an audit table.

Functional AMR presence is defined exactly as:

```text
Type == "AMR" and Method in {EXACTX, BLASTX, POINTX}
```

`PARTIALX`, `INTERNAL_STOP`, STRESS and other excluded records remain visible
with exclusion reasons but never create a functional 1. Coverage below 80% is
an audit flag, not an additional exclusion rule. The 72-row binary matrix is
the basis of prevalence, cautious shared/dataset-only descriptions, drug-class
summaries, Fisher tests and BH correction. Genomic detection is not phenotypic
resistance.

The unclustered and clustered AMR presence/absence heatmaps both use strain
names from the curated metadata and an explicit blue/red binary palette. The
clustered variant applies binary Jaccard distance to unscaled 0/1 calls and
average-linkage clustering to both strains and genes. Its full row/column
distance audit and clustering policy are written beside the heatmaps. A
`--heatmaps-only` mode regenerates only these heatmap artifacts from the frozen
validated AMR tables.

The older `amrfinder_analysis_72_genomes.R` is retained for legacy tabular
outputs and selected descriptive plots. Obsolete per-genome AMR-gene-count,
confidence-group and reference-coverage figures are not part of the current
figure set.

## 3. Nine-locus phylogeny — completed and frozen

The validated unrooted phylogeny uses the fixed order `gdh`, `gyd`, `pstS`,
`gki`, `xpt`, `yqiL`, `pyrC`, `groEL`, and `recA`. `aroE` is excluded from this
whole-span analysis under the documented observed-sequence policy; this does
not affect its use in formal seven-locus MLST. The 11,341-nt alignment contains
all 72 genomes, and IQ-TREE selected `TN+F+I+G4` with 1,000 ultrafast bootstrap
and 1,000 SH-aLRT replicates. The original analytical tree is unrooted;
midpoint rooting is display-only.

The matched 71-genome sensitivity analysis excluding `GCA_029011395.1` differs
materially after pruning to common tips, so it accompanies interpretation.
Detailed methods and limitations are in `docs/workflow/phylogeny_9_locus.md`.

## 4. Formal MLST and HLGR-associated genotype proxy — completed

The R workflow under `scripts/05_mlst/` assigns the official seven-locus order
`gdh-gyd-pstS-gki-aroE-xpt-yqiL` from pinned PubMLST resources. Known alleles
require 100% identity and 100% allele-length coverage; an ST requires an exact
seven-allele profile match. Sixty-nine genomes have formal ST assignments and
three remain explicitly unassigned. No nearest allele, ST or clonal complex is
invented.

The frozen AMR matrix supplies the primary
`aac(6')-Ie/aph(2'')-Ia` HLGR-associated genotype proxy. The validated state is
29 detected carriers (5/14 reproductive and 24/58 bacteraemia). This is a
genomic proxy, not a confirmed HLGR phenotype. ST prevalence is descriptive;
the ST-concentration, patristic-distance and AMR/tree concordance tests are
labelled exploratory. All analytical distance calculations use the frozen
unrooted nine-locus tree.

## 5. Targeted virulence/adherence analysis — completed through B4

The independent workflow under `scripts/07_virulence_adherence/` uses a pinned
official VFDB core Set A snapshot and BLASTN against the canonical assemblies.
Every candidate HSP and consolidated locus remains auditable. Accepted,
partial, review-required, ambiguous-multiple-hit and not-detected statuses are
mutually exclusive and reconcile to 14/58 for every target.

The approved resolved analysis applies competitive best-reference assignment.
It reports `aggregation_substance_family_detected` as the primary broad binary
feature and `asa1_specific` as a secondary lower-confidence breakdown. The
`GCA_029011395.1` `cylA` call is retained as accepted with copy number two.
Unresolved counts and maximum-possible uncertainty bounds are shown before
Fisher results. Optional integration with a core-genome tree remains pending
because no accepted core tree exists yet.

## 6. Panaroo core-genome analysis — A0/A1 complete; A2 review hold

The staged pipeline under `scripts/07_core_genome/` complements rather than
replaces the nine-locus tree. Panaroo 1.8.0 ran with strict cleaning, MAFFT and
a 0.95 core threshold on all 72 accession-labelled Prokka GFFs. A1 completed
successfully on 2026-08-17:

- 7,290 total pangenome families;
- 2,192 families present in at least 69/72 genomes;
- 5,098 accessory families below that threshold;
- 72 taxa in `core_gene_alignment.aln`;
- concatenated alignment length 2,118,553 nt.

Panaroo QC passed exact accession and 14/58 reconciliation. Alignment QC found
81,093 variable and 59,907 parsimony-informative sites, no ambiguous bases, and
one prespecified gate failure: `GCA_029011745.1` has 91.9047% non-missing
sequence (171,503 gaps), below the 95% minimum. The genome remains in the
alignment; it has not been silently excluded. IQ-TREE, SNP distances and
core/nine-locus tree comparison have not been run pending an explicit review
decision.

When that gate is resolved, the intended computationally bounded tree setting
is a fixed `GTR+G4` model with 1,000 ultrafast bootstrap replicates and four
threads, without ModelFinder or SH-aLRT. This setting must be confirmed before
execution; a fast fixed-model tree without bootstrap is the documented
time-constrained fallback.

## 7. Environments and reproducibility

The main environment is pinned in `environment/environment.yml`, with direct
package versions in `environment/software_versions.tsv` and historical run
versions in the recorded provenance files. The isolated Panaroo environment is
defined by `environment/core_genome.yml`; its solver decision and installed
versions are recorded separately. Scripts use project-relative paths, preserve
accession version suffixes and refuse to overwrite outputs unless an explicit
overwrite option is supplied.

## 8. Interpretation boundaries

The source groups are small and unbalanced, so prevalence and effect sizes are
reported alongside adjusted p-values. “Detected only in one dataset” is not a
specificity claim. Public-study structure, geography, incomplete metadata and
related isolates can confound source comparisons. Genomic AMR or virulence
detection does not establish expression, phenotype, pathogenicity,
transmission or clinical risk.
