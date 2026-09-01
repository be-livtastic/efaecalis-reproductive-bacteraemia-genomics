# Targeted virulence/adherence pipeline

This workflow screens only `ace`, `efaA`, `ebpA`, `ebpB`, `ebpC`, explicit `asa1`, `gelE`, `sprE`, `esp`, and `cylA`. It reads canonical FNA assemblies and frozen MLST/AMR/HLGR outputs. It does not reopen raw AMRFinderPlus files or infer phenotype.

The launcher discovers `mamba` or `micromamba` from `MAMBA_EXE` or `PATH` and
derives environments from `MAMBA_ROOT_PREFIX`. Non-standard installations can
set `EFAECALIS_VIRULENCE_ENV_PREFIX` and `EFAECALIS_MAIN_ENV_PREFIX` explicitly.

## Reference and detection policy

The reference is the official VFDB core Set A nucleotide release dated 2026-02-06. The raw archive, retrieval date, URL, SHA-256 checksums, selected sequences and symbol mapping are pinned under `references/virulence`. Explicit `asa1` is accepted; other aggregation-substance family members remain contextual/review-required. Additional cytolysin components remain visible without implying a complete operon or cytolytic phenotype.

Final statuses are mutually exclusive:

- `accepted_present`: at least 80% identity and 80% reference coverage;
- `flagged_partial`: at least 80% identity and 50–<80% coverage;
- `review_required`: 70–<80% identity with at least 80% coverage, or unresolved mapping;
- `ambiguous_multiple_hit`: multiple distinct accepted loci requiring review;
- `not_detected`: no retained candidate.

Only `accepted_present` contributes a 1. `Ebp_operon_complete` requires accepted `ebpA`, `ebpB`, and `ebpC`.

Overlapping matches from multiple VFDB references are consolidated to the same genomic locus before status assignment. Only multiple genuinely distinct accepted coordinate loci produce `ambiguous_multiple_hit`. The B3 screen writes `asa1_ambiguous_loci_72.csv` with genome/source/ST/study metadata and one coordinate-rich row per accepted ambiguous locus; these genomes remain excluded from binary presence and subsequent inference until reviewed.

## Stages and checkpoints

```bash
bash scripts/07_virulence_adherence/run_virulence_adherence_pipeline.sh --stage reference-audit
bash scripts/07_virulence_adherence/run_virulence_adherence_pipeline.sh --stage pilot
```

Stop after the four-genome pilot. Inspect raw HSPs, consolidated candidates and `virulence_pilot_status_qc.csv` before continuing:

```bash
bash scripts/07_virulence_adherence/run_virulence_adherence_pipeline.sh --stage screen
bash scripts/07_virulence_adherence/run_virulence_adherence_pipeline.sh --stage aggregation-review
bash scripts/07_virulence_adherence/run_virulence_adherence_pipeline.sh --stage resolve
bash scripts/07_virulence_adherence/run_virulence_adherence_pipeline.sh --stage analyse
```

The post-screen aggregation review compares every accepted explicit-`asa1` query locus competitively with the retained `prgB/asc10`, `EF0149` and `EF0485` VFDB references at the same coordinates. It writes locus-level evidence and genome-level recommendations but deliberately does not rewrite B3 statuses. Review and approve those recommendations before creating a resolved input for prevalence or Fisher testing.

The approved resolved layer reports `aggregation_substance_family_detected` as the primary binary comparison and `asa1_specific` as a secondary, lower-confidence breakdown. The rationale is retained in the staged analysis notes and the project’s QC records rather than in a standalone docs folder.

Every target/source prevalence row reports accepted, partial, review-required, ambiguous and not-detected counts. It also reports a maximum possible percentage if all unresolved calls were confirmed; this is an uncertainty bound, not another prevalence estimate. Figures label unresolved counts directly, and a Results-ready summary precedes statistical interpretation.

All 72 genomes must be accounted for, with 14 reproductive and 58 bacteraemia genomes for every target. Two-sided Fisher tests are run only after building the accepted-detection matrix, with BH correction across tested features. ST-level patterns are descriptive.

The Python tests cover interval consolidation, distant same-contig loci, all detection thresholds and multiple accepted loci. The R helper tests cover conservative Ebp completeness and BH correction:

```bash
python -m unittest scripts/07_virulence_adherence/test_virulence_screen.py
Rscript scripts/07_virulence_adherence/test_analysis_helpers.R
```
