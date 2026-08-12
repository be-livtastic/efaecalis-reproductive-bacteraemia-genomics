# Nine-locus concatenated housekeeping-gene phylogeny

## Status and biological rationale

This completed pipeline implements an extended MLST-style, whole-locus
analysis of 14 reproductive-associated and 58 bacteraemia-associated
*Enterococcus faecalis* assemblies. It is not a core-genome phylogeny.

The fixed locus order is `gdh`, `gyd`, `pstS`, `gki`, `xpt`, `yqiL`, `pyrC`,
`groEL`, and `recA`. The first six are retained from the seven-locus
*E. faecalis* MLST scheme; `aroE` is excluded consistently from all genomes
because one assembly contains a confirmed one-base deletion and no defensible
intact whole CDS. The additional three markers are conserved loci selected for
a defined multilocus analysis. Complete Prokka CDS coordinates are normally
used; two confirmed disrupted loci use exact homologous genomic spans from the
underlying FNA. Conventional internal MLST fragments are used only as
independent identity evidence, not as substitutes in the phylogenetic matrix.

## Reproducible environment

Create or update the dedicated Miniforge environment:

```bash
mamba env create -f environment/environment.yml
# If it already exists:
mamba env update -f environment/environment.yml --prune
mamba activate efaecalis_phylogeny
```

For shells without the activation hook:

```bash
source "$(mamba info --base)/etc/profile.d/conda.sh"
conda activate efaecalis_phylogeny
```

The container uses Miniforge and mamba. The ignored Prokka root must be mounted
read-only; it is never copied into the image.

## Annotation discovery and review

The discovery stage requires exactly 72 unique accessions, one GFF and one FNA
per accession, matching basenames, and complete agreement between GFF and FNA
contig identifiers. Final sample identifiers are bare assembly accessions.

GFF3 attributes are parsed structurally. Exact normalized `gene`, `Name`, and
`product` evidence is configured in `config/phylogeny_9_loci_genes.tsv`.
Published MLST primers and pinned reference information provide supporting
evidence for ambiguous annotations. Vague substrings are never sufficient.

All candidates are written to a review table. A unique high-confidence CDS may
be selected automatically; multiple, split, partial, pseudogene, or otherwise
ambiguous candidates block the run. Reviewed exceptions use
`config/phylogeny_9_loci_overrides.tsv`, based on the committed example, and
must name the feature, reviewer, date, and rationale. The pipeline never joins
split features or excludes a genome silently.

The revised policy requires 648 coordinate records (72 genomes by nine loci),
comprising 505 unique exact annotations and 143 reviewed `gyd`/`pstS`
overrides. Those paralog choices are supported independently by MLST primer
placement, length, conserved gene neighbourhood and alignment behaviour. A
reference-guided check of `GCA_029011745.1 aroE` found a one-base deletion at
homologous alignment column 194: all 71 intact genomes contain `C`, while the
underlying FNA does not. Joining the split features yields 866 nt and multiple
translation stops; correction would require inventing a nucleotide. `aroE` is
therefore excluded from every genome. The complete evidence is documented in
`docs/decisions/phylogeny_9_locus_annotation_review.md`.

## Coordinate extraction and strand handling

Coordinates are validated against contig lengths. `samtools faidx` indexes an
FNA only if its index is absent or older than the FNA. It extracts the annotated
interval in reference orientation. SeqKit is explicitly given
`--seq-type dna` and reverse-complements minus-strand features with
`seqkit seq --seq-type dna -r -p`. Final FASTA headers contain only the assembly
accession, with exactly 72 unique sequences required for every locus.

## Sequence and alignment QC

Biopython reports nucleotide length, modulo three, start and terminal codons,
internal stops, ambiguous bases and residues, translation length, median-length
deviation, possible truncation, and possible frameshift. Alternative bacterial
start codons are reported rather than rejected solely for not being `ATG`.
`config/phylogeny_9_loci_qc_thresholds.tsv` defines locus-specific expected,
PASS and REVIEW length ranges plus alignment gap and divergence limits. Values
outside the REVIEW envelope FAIL. `config/phylogeny_9_loci_reviewed_findings.tsv`
adds evidence-backed codes such as `SPLIT_ANNOTATION`, `ASSEMBLY_INDEL`, and
`POSSIBLE_FRAMESHIFT`.

The validator also writes an assembly-level summary and flags genomes with more
than one suspicious locus. Formal strand-aware neighbourhood output records the
two upstream and two downstream CDSs for every selected locus.

Internal stops, frameshifts, major truncation, duplicates, or unresolved
annotations block MAFFT unless a named record has an explicit reviewed
observed-disruption policy. Each locus is aligned independently with
`mafft --auto --thread N`. Alignment QC requires the same 72 identifiers,
uniform aligned length, and reports gaps, ambiguity, variable sites, and
parsimony-informative sites.

The project default is two threads to limit resource pressure. The diagnostic
reference alignments were also rerun with `mafft --auto --thread 2`; they
confirmed the yqiL insertion at alignment column 1226 and the gdh deletions at
columns 1427-1428. On the review machine these 72-sequence checks took 1:42
(about 169 MiB peak resident memory) for yqiL and 0:39 (about 74 MiB) for gdh.
Production runtimes will vary by locus and computer.

Post-MAFFT QC reports gap and ambiguity proportions, variable and
parsimony-informative sites, all pairwise divergences, and each sequence's
median divergence from its peers. Excessive-gap or divergence outliers block
concatenation before FASTA or partition outputs are written.

The final extracted set contains 648 sequences and all nine locus FASTAs have
72 unique identifiers. `aroE` is no longer an input. The explicit policy keeps
the exact 2413-nt yqiL and 1046-nt recA genomic spans from
`GCA_029011395.1`; the one-base insertion/deletion and resulting disruption
are genuine observed characters, and no nucleotide is imputed. The shorter
1488-nt `GCA_029011535.1 gdh` allele is retained as an approved terminal
variant. Regenerated sequence QC reports 648 PASS, 0 REVIEW, and 0 FAIL.

Reference-guided follow-up confirms that the remaining failures are not simple
Prokka boundary errors. `GCA_029011395.1` has a one-base insertion that splits
the normal 2412 nt acetyl-CoA acetyltransferase/HMG-CoA reductase fusion into
separate `thlA` and `mvaA` ORFs, plus a different one-base deletion that splits
the terminal portion of `recA`. The short `GCA_029011535.1 gdh` also has a
near-terminal two-base deletion relative to the cohort. None can be restored to
a cohort-like whole CDS without changing underlying assembly nucleotides.

## Concatenation, partitioning, and inference

Only aligned loci are concatenated, in the fixed order above. The pipeline
validates that partitions cover the concatenated alignment exactly without
gaps or overlaps. IQ-TREE 2 syntax is checked before running a partitioned
ModelFinder/merge analysis with 1,000 ultrafast bootstrap and 1,000 SH-aLRT
replicates. Branch lengths and support values are distinct quantities.

## Metadata integration and visualisation

Tree tips join to `data/metadata/curated_metadata_72_genomes.csv` using
`assembly_accession`. Acceptance requires exactly 72 tips, 72 unique metadata
keys, and no unmatched records. The original inferred topology is retained
unrooted. Midpoint rooting is a display transformation only. Reproductive-only
and bacteraemia-only panels are pruned from the combined tree rather than
inferred separately.

## Running safely

```bash
scripts/06_phylogenomics/run_9_locus_phylogeny.sh --threads 2 --dry-run
scripts/06_phylogenomics/run_9_locus_phylogeny.sh --threads 2
```

The first unresolved annotation returns a non-zero review-gate status.
Destinations are never overwritten. `--force` creates a timestamped versioned
workspace. IQ-TREE and visualisation require explicit `--stage iqtree` and
`--stage visualisation` invocations after all upstream review gates pass.

Every invoked stage records the command, standard output, standard error,
runtime, exit status, inputs, outputs, and timestamp. Successful stages receive
a completion marker and command hash. MAFFT writes temporary files first; only
a successful 72-record alignment is promoted, hashed, and marked complete.
The policy preflight checks locus order/count, the 648-record invariant, active
prefixes, and figure naming. Bulk intermediates remain local and ignored.

## Trusted-locus validation and sensitivity result

The pre-alignment workflow performs two small, versioned checks. PubMLST
alleles dated 2026-08-06 validate the six active official MLST loci: 432/432
sample-locus comparisons pass. Nine pinned RefSeq V583 proteins provide BLASTP
identity plus bidirectional coverage checks: 646/648 pass and the yqiL/recA
records are the two expected disruptions, with no unexpected failures. Full
reference coverage also acts as a protein-length/domain-architecture screen;
no claim of a separate Pfam/InterPro domain search is made.

The primary 72-genome alignment is 11341 nt and every genome has at least
99.6738% non-missing characters (policy minimum 95%). IQ-TREE selected
`TN+F+I+G4` and completed 1000 ultrafast bootstrap and 1000 SH-aLRT replicates.
The matched nine-locus sensitivity tree excluding `GCA_029011395.1` differs
materially after pruning to 71 shared tips (normalized RF 0.352941), so this
sensitivity must accompany biological interpretation.

## Limitations

The result is a **nine-locus concatenated housekeeping-gene phylogeny**. It is
not a core-genome phylogeny and, because whole CDSs and two additional markers
are used, it is not a standard seven-locus MLST phylogeny.

Nine loci represent only a small fraction of the genome and may not reproduce
relationships from a recombination-aware core-genome analysis. Individual
loci may have different evolutionary histories. Unequal sample sizes,
study/geographic bias, incomplete public metadata, genotype–phenotype
limitations, and repeated or epidemiologically linked isolates constrain
interpretation. Genomic proximity to a bacteraemia-associated isolate does not
prove pathogenicity or clinical risk.

## Downstream topology comparison and PCA

`09_compare_completed_trees.R` compares jointly supported bipartitions,
reproductive-isolate neighbourhoods, and terminal branch lengths across the
primary 72-genome nine-locus tree, the 71-genome exclusion sensitivity tree,
and the 72-genome seven-locus sensitivity tree.

`run_9_locus_pca.sh` runs a reproducible PCA of the primary concatenated
alignment. `10_run_alignment_pca.py` validates identifiers, encodes observed
alternative alleles occurring in at least two genomes, standardises features,
and performs SVD. `10_plot_alignment_pca.R` creates labelled PC1-PC2 and
PC1-PC3 plots. Missing numeric values are assigned the feature mean after
centring solely for matrix decomposition; no sequence file is changed.
