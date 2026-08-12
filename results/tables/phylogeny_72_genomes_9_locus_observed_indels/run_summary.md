# 72-genome nine-locus observed-indel phylogeny

Status: **Primary inference and requested sensitivity inference completed.**

## Explicit analysis policy

All 72 genomes and all nine loci are retained. For confirmed assembly-level
disruptions, the analysis uses only nucleotides present in the FNA across the
homologous genomic locus. No base is inserted, deleted, corrected, or imputed.

- `GCA_029011395.1 yqiL`: exact 2413-nt `mvaA-thlA` genomic union; observed
  one-base insertion retained.
- `GCA_029011395.1 recA`: exact 1046-nt recA/adjacent-fragment union; observed
  one-base deletion retained.
- `GCA_029011535.1 gdh`: observed 1488-nt terminal variant retained; its 2.36%
  gap burden is below the configured 3% locus limit.

Frameshift and internal-stop signals remain recorded as information codes for
the two explicitly reviewed disrupted loci. This exception is not available to
unreviewed records.

## Regenerated validation results

- Annotation inputs: 72/72 valid; 14 reproductive and 58 bacteraemia.
- Coordinate selection: 648/648; no unresolved combinations.
- Observed-span overrides: 2; imputed nucleotides: 0.
- Extracted FASTAs: 72 unique identifiers for every locus.
- Sequence QC: 648 PASS, 0 REVIEW, 0 FAIL under the explicit disruption policy.
- Neighbourhood QC: 648 selected loci and 3240 selected/flanking CDS rows.
- PubMLST fragment QC: 432/432 PASS for `gdh`, `gyd`, `pstS`, `gki`, `xpt`, and
  `yqiL` against the 2026-08-06 allele snapshot.
- RefSeq V583 BLASTP QC: 646 PASS, 2 EXPECTED_DISRUPTION, 0 unexpected FAIL.
- Alignment: nine independent `mafft --auto --thread 2` runs, each with 72 taxa.
- Concatenated matrix: 72 taxa, 11341 nt, 320 variable sites, 253
  parsimony-informative sites.
- Minimum per-genome non-missing proportion: 0.996738 (required: 0.95).
- IQ-TREE composition test: 0 failures.

## Primary inference

- IQ-TREE 2: `MFP+MERGE`, 1000 ultrafast bootstrap replicates, and 1000 SH-aLRT
  replicates using two threads.
- Selected merged model: `TN+F+I+G4` across all nine loci.
- Best log-likelihood: -19854.5835.
- Tree and figures retain all 72 exact metadata-matched tips.

## Sensitivity analysis

An otherwise identical nine-locus inference excluding only
`GCA_029011395.1` selected the same `TN+F+I+G4` model. Comparing its 71-taxon
tree with the primary tree after pruning that tip gave Robinson-Foulds distance
48/136 (normalized 0.352941). The refitted topology is therefore not identical.
This sensitivity is material and must accompany interpretation, especially
because the alignment contains many identical concatenated sequences.

This result is a **nine-locus concatenated housekeeping-gene phylogeny**, not a
core-genome phylogeny.
