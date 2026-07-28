# Ten-locus concatenated housekeeping-gene phylogeny

## Status and biological rationale

This work-in-progress pipeline implements an extended MLST-style, whole-CDS
analysis of 14 reproductive-associated and 58 bacteraemia-associated
*Enterococcus faecalis* assemblies. It is not a core-genome phylogeny.

The fixed locus order is `gdh`, `gyd`, `pstS`, `gki`, `aroE`, `xpt`, `yqiL`,
`pyrC`, `groEL`, and `recA`. The first seven constitute the established
*E. faecalis* MLST scheme. The additional markers are conserved loci selected
for a defined multilocus analysis. Complete Prokka CDS coordinates are used;
the conventional internal MLST fragments are not substituted for them.

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
`product` evidence is configured in `config/phylogeny_10_loci_genes.tsv`.
Published MLST primers and pinned reference information provide supporting
evidence for ambiguous annotations. Vague substrings are never sufficient.

All candidates are written to a review table. A unique high-confidence CDS may
be selected automatically; multiple, split, partial, pseudogene, or otherwise
ambiguous candidates block the run. Reviewed exceptions use
`config/phylogeny_10_loci_overrides.tsv`, based on the committed example, and
must name the feature, reviewer, date, and rationale. The pipeline never joins
split features or excludes a genome silently.

## Coordinate extraction and strand handling

Coordinates are validated against contig lengths. `samtools faidx` indexes an
FNA only if its index is absent or older than the FNA. It extracts the annotated
interval in reference orientation. SeqKit reverse-complements minus-strand
features. Final FASTA headers contain only the assembly accession, with exactly
72 unique sequences required for every locus.

## Sequence and alignment QC

Biopython reports nucleotide length, modulo three, start and terminal codons,
internal stops, ambiguous bases and residues, translation length, median-length
deviation, possible truncation, and possible frameshift. Alternative bacterial
start codons are reported rather than rejected solely for not being `ATG`.
The default length warning threshold is ±10% from the locus median.

Internal stops, frameshifts, major truncation, duplicates, or unresolved
annotations block MAFFT. Each locus is aligned independently with
`mafft --auto --thread N`. Alignment QC requires the same 72 identifiers,
uniform aligned length, and reports gaps, ambiguity, variable sites, and
parsimony-informative sites.

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
scripts/06_phylogenomics/run_10_locus_phylogeny.sh --threads 4 --dry-run
scripts/06_phylogenomics/run_10_locus_phylogeny.sh --threads 4
```

The first unresolved annotation returns a non-zero review-gate status.
Destinations are never overwritten. `--force` creates a timestamped versioned
workspace. IQ-TREE and visualisation require explicit `--stage iqtree` and
`--stage visualisation` invocations after all upstream review gates pass.

Every invoked stage records the command, standard output, standard error,
runtime, exit status, inputs, outputs, and timestamp. Bulk intermediates and
checkpoints remain local and ignored.

## Limitations

Ten loci represent only a small fraction of the genome and may not reproduce
relationships from a recombination-aware core-genome analysis. Individual
loci may have different evolutionary histories. Unequal sample sizes,
study/geographic bias, incomplete public metadata, genotype–phenotype
limitations, and repeated or epidemiologically linked isolates constrain
interpretation. Genomic proximity to a bacteraemia-associated isolate does not
prove pathogenicity or clinical risk.
