# Ten-locus concatenated housekeeping-gene phylogeny run summary

Status: **Not yet run beyond input/candidate review.**

This file is a reproducibility template. Values must be populated from a successfully
validated, timestamped pipeline run and must not be inferred or invented.

## Inputs and locus validation

- Genomes expected: 72 (14 reproductive-associated; 58 bacteraemia-associated)
- Genomes found: pending
- Validated GFF/FNA pairs: pending
- Loci expected: 10 per genome; 720 total
- Loci successfully extracted: pending
- Missing or duplicate loci: pending reviewed report

## Sequence and alignment QC

- Sequence QC: pending
- Per-locus alignment lengths: pending
- Variable and parsimony-informative sites: pending
- Concatenated length: pending
- Partition coordinates: pending

## Phylogenetic inference and metadata

- IQ-TREE selected models/partition scheme: not run
- Ultrafast bootstrap and SH-aLRT completion: not run
- Metadata join: pending; acceptance requires 72 tree tips and 72 exact matches
- Generated figures: none

## Warnings and limitations

This is a defined ten-locus whole-CDS analysis based on seven established MLST
loci and three additional conserved markers. It is not a pangenome-derived
core-genome phylogeny. Unequal group sizes, sampling and geographic bias,
incomplete public metadata, genotype–phenotype limitations, and potentially
repeated or epidemiologically linked isolates constrain interpretation.
Phylogenetic proximity does not establish pathogenicity or transmission.

## Reproduction

```bash
mamba env create -f environment/environment.yml
mamba activate efaecalis_phylogeny
scripts/06_phylogenomics/run_10_locus_phylogeny.sh --threads 4 --dry-run
scripts/06_phylogenomics/run_10_locus_phylogeny.sh --threads 4
```

The second command intentionally stops at the reviewed-locus gate. Exact downstream
commands will be appended only after candidate selections pass scientific review.
