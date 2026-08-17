# Formal MLST and HLGR-proxy integration

This workflow assigns formal seven-locus *Enterococcus faecalis* sequence types and integrates them with the validated AMRFinderPlus analysis and the existing 72-tip nine-locus tree. It does not modify or rerun either upstream pipeline.

## Inputs and frozen boundaries

- Canonical accessions: `data/accession_lists/selected_72_accessions.tsv`
- Canonical genomes: the 14- and 58-genome `02_genomes_fna` directories in `local_archive/large_outputs`
- Study metadata: `data/metadata/sample_metadata_source_72.tsv`; `Study_ID` is the BioProject accession
- Frozen AMR matrix: `results/tables/amr/amr_gene_presence_absence.csv`
- Frozen accepted hits: `data/processed/amr/amr_functional_hits_with_contigs_72.csv`
- Frozen unrooted tree: the primary 72-genome IQ-TREE `.treefile`

The integration scripts never reopen raw AMRFinderPlus TSVs or reinterpret `Type` and `Method`. The known-invalid legacy HLGR columns in `curated_metadata_72_genomes.csv` are not analytical inputs.

## PubMLST reference and assignment policy

`01_prepare_pubmlst_reference.R` retrieves official resources from the PubMLST REST API and pins the response, endpoints, retrieval time, schema, record counts and SHA-256 checksums. The formal locus order is:

`gdh-gyd-pstS-gki-aroE-xpt-yqiL`

An official allele number requires 100% nucleotide identity and 100% allele-length coverage. Formal ST assignment requires seven accepted exact alleles and an exact match to one pinned official profile. Nearest alleles, nearest profiles and clonal complexes are never invented. A curator grouping is imported only if it is an explicit profile-response field.

`aroE` is valid here because conventional MLST uses an internal fragment. This does not change the separate whole-span decision made for the nine-locus phylogeny.

Locus statuses are `Multiple_exact_alleles`, `Multiple_locus_copies`, `Exact_known_allele`, `Novel_allele_candidate`, `Missing_locus`, `Non_exact_best_match`, and `Unresolved`. A duplicated locus contributes only when all exact copies carry the same official allele. Unambiguous novel sequences receive stable `NOVEL_<hash>` tokens for the dataset-profile sensitivity analysis but never receive official allele numbers.

## HLGR-associated genotype proxy

The primary proxy is accepted upstream detection of `aac(6')-Ie/aph(2'')-Ia`. Outputs use only “HLGR-associated genotype proxy detected/not detected.” Genomic detection is not confirmed high-level gentamicin resistance; phenotypic susceptibility testing would be required.

`config/hlgr_proxy_determinants.tsv` records the approved proxy and review watchlist. An accepted non-primary aminoglycoside determinant on the watchlist, or one whose frozen AMRFinderPlus subclass indicates gentamicin, produces `HLGR_proxy_review_required` and stops final classification.

The integration gate expects the validated upstream state of 29 carriers: 5/14 reproductive-associated and 24/58 bacteraemia-associated. A changed count stops the run for review.

## Statistics and interpretation

ST prevalence is descriptive and primary. No per-ST Fisher tests are run.

The exploratory ST concentration statistic is `T = sum_s choose(k_s, 2)`, where `k_s` is the carrier count in ST `s`. The one-sided 100,000-permutation test holds ST assignments and the assigned carrier count fixed; larger values are more concentrated. It is repeated for complete dataset-specific allele profiles.

Phylogenetic clustering uses mean pairwise patristic distance among the 29 carriers on the original unrooted tree. Unrestricted and source-stratified 100,000-permutation nulls use the prespecified lower tail because smaller distances indicate clustering. Midpoint rooting is presentation-only.

The AMR Jaccard/patristic Mantel-style analysis and pairwise similarity strata are exploratory. Pairwise summaries are descriptive because genome pairs are non-independent and source, ST and study origin overlap.

## Outputs

Processed MLST outputs in `data/processed/mlst`:

- `mlst_allele_hit_audit_72.tsv`
- `mlst_locus_assignment_qc_72.csv`
- `formal_mlst_assignments_72.csv`
- `mlst_assignment_qc_summary_72.csv`

Integrated tables in `results/tables/mlst_amr_phylogeny`:

- `integrated_genome_mlst_amr_72.csv`
- `hlgr_proxy_supporting_hits_72.csv`
- `hlgr_proxy_review_audit_72.csv`
- `st_hlgr_proxy_prevalence_72.csv`
- `st_concentration_permutation_exploratory_72.csv`
- `phylogenetic_clustering_permutations_72.csv`
- `amr_phylogeny_mantel_exploratory_72.csv`
- `amr_pairwise_relationship_audit_72.csv`
- `amr_pairwise_similarity_summary_72.csv`
- `mlst_amr_phylogeny_qc_summary_72.csv`

Figures in `results/figures/mlst_amr_phylogeny` include a presentation-only midpoint tree and metadata tracks. No automatic tree-group classification is produced.

`analysis/logs/mlst_amr_phylogeny/mlst_amr_phylogeny_session_info_72.txt` records the command sequence and R session/software versions.

## Running

Use the project environment, then run once without overwrite:

```bash
Rscript scripts/05_mlst/01_prepare_pubmlst_reference.R --date=2026-08-16
Rscript scripts/05_mlst/02_assign_mlst_72_genomes.R \
  --reference-dir references/mlst/pubmlst/efaecalis_2026-08-16
Rscript scripts/05_mlst/03_integrate_mlst_amr_phylogeny.R
Rscript scripts/05_mlst/04_visualise_mlst_amr_phylogeny.R
```

Every stage refuses existing outputs unless `--overwrite` is supplied. Review the pinned-reference provenance, locus QC, assignment QC and integration QC before interpreting statistical tables or figures.
