# Nine-locus annotation review and sequence-QC decision

## Decision status

The active policy excludes `aroE` consistently from all 72 genomes. Script 02
must select exactly 648 records: 72 genomes by nine loci, with 505 unique exact
annotations and 143 reviewed overrides. Script 03 must produce 648 individual
FASTAs and nine locus FASTAs, each containing 72 unique assembly-accession
identifiers.

Tree inference must not proceed yet. Excluding `aroE` resolves its inconsistent
whole-CDS representation, but the short `yqiL` and `recA` still fail sequence
QC.

## Independent evidence used for overrides

Primer matching was evaluated against the nucleotide sequences in strand-aware
CDS orientation. It therefore provides a positional and sequence-identity
check independent of the Prokka product label. Primer evidence was considered
together with copy number, length, conserved neighbouring genes, coding
integrity and temporary MAFFT alignment behaviour.

### `gyd`

- One 1002 nt copy was selected per genome.
- Both pinned MLST primers matched exactly in all 72 selected copies.
- The 71 rejected 1011 nt paralogs had four mismatches to each primer.
- All selected copies shared the conserved `cggR-gyd-pgk-tpiA` context.
- Rejected copies instead shared an `obg-gyd-fur-cvfB` context.
- The 72-sequence alignment had 1002 columns, no gaps, six variable sites and
  a maximum consensus distance of 0.30%.

The 71 ambiguous `gyd` overrides are retained.

### `pstS`

- One 891 nt copy was selected per genome.
- Sixty-seven selected copies matched both primers exactly; five had one
  mismatch at the reverse-primer site, consistent with allelic variation.
- Rejected 855 nt paralogs had five or six forward-primer mismatches and nine
  reverse-primer mismatches.
- Selected copies shared the `dapX-pstS-phoR-phoP` phosphate-regulation
  neighbourhood.
- Rejected copies instead occurred beside `ftsX` and `pstC1`.
- The 72-sequence alignment had 891 columns, no gaps, 19 variable sites and a
  maximum consensus distance of 0.56%.

The 72 ambiguous `pstS` overrides are retained.

### Reference-guided `aroE` investigation and exclusion decision

`GCA_029011745.1` contains overlapping Prokka features `aroE_1` (663 nt) and
`aroE_2` (228 nt), compared with an intact cohort length of 867 nt. The 663 nt
segment contains both MLST primers; the 228 nt segment contains neither. The
broader locus remains beside `aroF`, consistent with the expected genomic
region.

This evidence confirms the locus and the internal MLST target, but not an intact
whole CDS. In the 72-sequence `aroE` alignment, this record introduces 204 gaps
and has the highest consensus distance (1.51%).

Reference-guided inspection of the existing FNA then established why Prokka
split the feature:

- The union of coordinates `CP118085.1:1624528-1625393` on the minus strand is
  866 nt, rather than the intact 867 nt cohort length.
- MAFFT places one deletion at homologous alignment column 194. All 71 intact
  genomes carry `C` at that column; the nucleotide is absent from this FNA.
- Translation of the full 866 nt union contains multiple stops after the
  deletion, and scanning the surrounding sequence finds no intact 800-900 nt
  bacterial ORF.
- The only long clean ORF is the already annotated downstream 663 nt fragment.

Producing an intact 867 nt sequence would require inserting a nucleotide not
present in the assembly. That would be imputation, not coordinate correction.
Reference-guided reconstruction was therefore rejected. The approved fallback
excludes `aroE` from all 72 genomes, yielding a nine-locus dataset while
retaining the complete 72-genome cohort. The prior ten-locus files remain only
as review provenance and must not be used for inference.

## Length caveats

### Short `gdh`: `GCA_029011535.1`

The selected `gdh` is 1488 nt compared with the 1524 nt cohort mode, a 36 nt
difference. The sequence is divisible by three, has no internal stop, and is
within the configured 10% length tolerance. The difference corresponds to 12
codons and could reflect an alternative start site, terminal deletion,
annotation truncation or genuine allele variation rather than a frameshift.
Direct cohort alignment shows that this is not solely an alternative start-site
call. The sequence lacks two bases at homologous alignment columns 1427-1428
and then has 34 terminal gap positions relative to the 1524 nt cohort sequence.
The selected Prokka feature remains in-frame because it terminates after the
resulting altered C-terminal segment, but the evidence is consistent with a
near-C-terminal frameshift and early stop in the underlying assembly. Restoring
the cohort-like CDS would require inserting bases that are not present. It
passes the automated ±10% rule but remains biologically unresolved.

### Short `yqiL`: `GCA_029011395.1`

The selected `yqiL` is 1245 nt compared with the 2412 nt cohort mode: a 1167 nt
difference and approximately half the usual length. It is divisible by three,
starts with `TTG`, ends with `TGA`, and has no internal stop or ambiguous base,
but official QC identifies it as outside length tolerance and possibly
truncated. Plausible explanations include a split annotation, incorrect
paralog selection, gene disruption, annotation-boundary error or a frameshift
outside the selected feature. This is a critical pre-tree record.

Reference-guided investigation shows that the cohort's 2412 nt feature is a
real acetyl-CoA acetyltransferase/HMG-CoA reductase fusion rather than a simple
803-aa thiolase. Its first approximately 409 amino acids align to `thlA`, while
its final approximately 394 amino acids align to `mvaA`. In
`GCA_029011395.1`, the homologous genomic union is 2413 nt and contains one
extra `A` at alignment column 1226. That insertion breaks the fusion into the
overlapping, individually clean 1245 nt `thlA` and 1191 nt `mvaA` predictions.
The pinned `yqiL` primers occur only within the `thlA` portion (positions 56 and
557 in CDS orientation). The original MLST method assays an internal fragment,
so primer recovery does not validate the complete fusion CDS. Producing a
2412 nt whole CDS for this assembly would require deleting a nucleotide that is
present in its FNA and is therefore not defensible coordinate correction.

This interpretation is consistent with the original *E. faecalis* MLST method,
which explicitly amplified internal fragments, and with the V583 reference
protein annotated as an acetyl-CoA acetyltransferase/HMG-CoA reductase fusion:

- Ruiz-Garbajosa et al. (2006), DOI
  [10.1128/JCM.02596-05](https://doi.org/10.1128/JCM.02596-05)
- NCBI protein AAG02439.1,
  [acetyl-CoA acetyltransferase/HMG-CoA reductase](https://www.ncbi.nlm.nih.gov/protein/AAG02439.1)

### Short `recA`: `GCA_029011395.1`

The selected `recA` is 942 nt compared with the 1047 nt cohort mode, a 105 nt
or 35-codon difference. It is divisible by three, starts with `TTG`, ends with
`TGA`, and has no internal stop or ambiguous base. Its 10.03% deviation narrowly
exceeds the configured 10% tolerance. Because `recA` is normally highly
conserved, its annotation boundary and alignment require direct inspection.

The adjacent 108 nt hypothetical feature overlaps the short `recA` by four
bases and lies in the same orientation. Their genomic union is 1046 nt, one
base shorter than the 1047 nt cohort CDS. Alignment places a missing `C` at
homologous column 894. Translation in the original frame then stops after the
314-aa Prokka `recA`, while the residual terminal region is predicted
separately. This supports a genuine assembly-level frameshift rather than a
recoverable boundary mistake. Restoring a complete 1047 nt `recA` would require
inserting a base absent from the FNA.

The short `yqiL` and `recA` occur in the same assembly, `GCA_029011395.1`.
Their co-occurrence raises an assembly- or annotation-quality concern rather
than supporting an assumption of two independent biological truncations. The
direct sequence evidence shows two distinct one-base indels in that assembly:
an insertion in the `yqiL` fusion and a deletion in `recA`.

## Official sequence-QC result

Biopython 1.87 is installed in the `efaecalis_phylogeny` environment. The
historical ten-locus validation reported 720 records and three failures:

| Assembly | Locus | Length | Cohort median/mode | QC result |
|---|---:|---:|---:|---|
| `GCA_029011745.1` | `aroE` | 663 | 867 | Outside tolerance; possible truncation |
| `GCA_029011395.1` | `yqiL` | 1245 | 2412 | Outside tolerance; possible truncation |
| `GCA_029011395.1` | `recA` | 942 | 1047 | Outside tolerance |

None of these three contains an internal stop, ambiguous nucleotide or
out-of-frame length. Those properties do not overrule the comparative length
and alignment evidence.

Under the nine-locus policy, the `aroE` row is absent by design. The former
universal 10% rule reported 646 passes and two failures. Locus-specific QC now
reports 645 PASS, one REVIEW (`GCA_029011535.1 gdh`), and two FAIL
(`GCA_029011395.1 yqiL` and `recA`). The failures carry `SPLIT_ANNOTATION` and
`ASSEMBLY_INDEL`; gdh carries `LENGTH_REVIEW` and `POSSIBLE_FRAMESHIFT`. The
assembly summary flags `GCA_029011395.1` with multiple suspicious loci.

## Final primary analysis policy

The manual checkpoint now approves retention of all 72 genomes and all nine
loci using an observed-disrupted-locus policy. For `GCA_029011395.1`, yqiL is
the exact 2413-nt `mvaA-thlA` FNA union and recA is the exact 1046-nt
recA/adjacent-fragment FNA union. Their observed one-base insertion/deletion is
retained. Nothing is imputed. The 1488-nt gdh terminal variant is retained with
an explicit caveat because its 2.36% gap burden is below the configured 3%
locus threshold.

Regenerated QC reports 648 PASS, 0 REVIEW, and 0 FAIL under this narrowly scoped
policy. The affected genome has 99.6738% non-missing characters in the 11341-nt
concatenation, exceeding the required 95%. PubMLST fragment checks pass 432/432
applicable comparisons; pinned V583 RefSeq BLASTP checks pass 646 records and
classify only the two known disrupted translations as EXPECTED_DISRUPTION.

Primary IQ-TREE inference completed with all 72 genomes. A same-locus
sensitivity inference excluding `GCA_029011395.1` gave normalized
Robinson-Foulds distance 0.352941 after matching the 71 shared tips, so topology
sensitivity is explicitly retained as a limitation.

