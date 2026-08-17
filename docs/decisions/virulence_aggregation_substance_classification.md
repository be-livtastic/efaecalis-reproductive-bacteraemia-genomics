# Aggregation-substance classification and multiple-locus decision

## Decision status

The active targeted virulence/adherence policy reports two separate features. `aggregation_substance_family_detected` is the primary comparison and records accepted genomic detection of any retained VFDB aggregation-substance family reference after coordinate-local locus consolidation. `asa1_specific` is a secondary, lower-confidence narrow breakdown based on competitive reference support.

The original B3 calls and all BLAST HSPs remain immutable audit evidence. This decision is applied only in `virulence_resolved_target_status_by_genome_72.csv`.

## Evidence prompting review

The B3 screen found 67 accepted-threshold matches to the explicit VFDB `asa1` reference across 54 genomes. Competitive comparison at the same genomic coordinates showed that 15 loci were best supported by the explicit `asa1` reference, whereas 52 were better supported by retained `prgB/asc10`, `EF0149`, or `EF0485` references. Treating every accepted-threshold `asa1` query match as gene-specific presence would therefore conflate related aggregation-substance family members.

Multiple matches from different references at overlapping coordinates are one genomic locus. Only coordinate-separated loci remain distinct. All reference IDs, coordinates, strands, identity, coverage, HSP counts and contig identifiers remain in the locus-level audits.

## Approved analytical policy

- `aggregation_substance_family_detected` is accepted when at least one coordinate-consolidated family locus meets at least 80% nucleotide identity and 80% reference coverage against any retained aggregation-substance family reference.
- Multiple accepted family loci contribute one binary presence call and are retained as copy-number/context information.
- `asa1_specific` is accepted only when the explicit `asa1` reference has the highest BLAST bit score among retained family references overlapping that locus.
- A genome with accepted-threshold `asa1` query matches but no locus competitively best supported by `asa1` remains `review_required` for the narrow feature and contributes no narrow binary presence.
- This yields 13 `asa1_specific` genomes. The count is an implementation assertion.
- Prevalence and Fisher testing use `aggregation_substance_family_detected` as the primary aggregation-substance feature. `asa1_specific` is reported alongside as a secondary, lower-confidence breakdown.

## Separate `cylA` resolution

`GCA_029011395.1` contains two full-length `cylA` matches at 99.516% identity and 100% coverage: one on chromosome `CP118057.1` and one on plasmid `CP118058.1`. Its resolved genome-level call is accepted presence with copy number 2. Both loci remain in the candidate audit; binary presence remains one.

## Interpretation limitation

These are reference-based genomic classifications, not evidence of expression, functional aggregation substance, cytolytic phenotype, or virulence in a host. The narrow `asa1_specific` category remains lower confidence because homologous aggregation-substance genes are highly similar and definitive naming may require broader locus and functional review.

## Machine-readable evidence

- `results/tables/virulence_adherence/aggregation_substance_locus_identity_review_72.csv`
- `results/tables/virulence_adherence/aggregation_substance_genome_review_recommendations_72.csv`
- `results/tables/virulence_adherence/asa1_ambiguous_loci_72.csv`
- `data/processed/virulence_adherence/aggregation_substance_family_accepted_loci_72.csv`
- `data/processed/virulence_adherence/virulence_resolved_target_status_by_genome_72.csv`
