# Nine-locus concatenated housekeeping-gene phylogeny run summary

Status: **Sequence QC blocked; inference not run.**

## Inputs and locus validation

- Genomes: 72 (14 reproductive-associated; 58 bacteraemia-associated)
- Loci: 9 per genome; 648 total
- Fixed order: `gdh`, `gyd`, `pstS`, `gki`, `xpt`, `yqiL`, `pyrC`, `groEL`, `recA`
- `aroE`: excluded consistently after reference-guided review found a one-base
  deletion in `GCA_029011745.1` and no defensible intact whole CDS
- Reviewed overrides: 143 (`gyd` = 71; `pstS` = 72)

## Sequence and alignment QC

- Locus-specific sequence QC: 648 records; 645 PASS; 1 REVIEW; 2 FAIL
- Outstanding failures: assembly-level one-base frameshifts affecting `yqiL`
  and `recA` in `GCA_029011395.1`
- Formal REVIEW: `gdh` in `GCA_029011535.1`; cohort alignment shows a
  near-terminal two-base deletion and early termination
- Assembly QC: `GCA_029011395.1` has two failed loci and is flagged
  `multiple_suspicious_loci=true`
- Alignment: blocked by sequence-QC exit status 3 until the remaining failures
  are resolved
- Concatenation and partitions: not generated

## Phylogenetic inference and metadata

- IQ-TREE selected models/partition scheme: not run
- Ultrafast bootstrap and SH-aLRT completion: not run
- Metadata join: pending; acceptance requires 72 tree tips and 72 exact matches
- Generated figures: none

## Warnings and limitations

This is a **nine-locus concatenated housekeeping-gene phylogeny** containing six
established MLST loci and three additional conserved markers. It is neither a
pangenome-derived core-genome phylogeny nor standard seven-locus MLST. The
excluded `aroE` locus and all remaining sequence caveats must be reported with
any eventual tree.

## Reproduction

```bash
mamba env create -f environment/environment.yml
mamba activate efaecalis_phylogeny
scripts/06_phylogenomics/run_9_locus_phylogeny.sh --threads 2 --dry-run
scripts/06_phylogenomics/run_9_locus_phylogeny.sh --threads 2
```

The second command intentionally stops if the remaining sequence-QC failures
have not been resolved. IQ-TREE and visualisation remain explicit later stages.
