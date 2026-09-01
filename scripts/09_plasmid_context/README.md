# Plasmid-context module

This downstream-only module maps each accepted hit in the frozen AMRFinderPlus
long table to its canonical NCBI Datasets `sequence_report.jsonl` replicon
record. It joins canonical metadata and formal MLST assignments, and it does not
rerun AMRFinderPlus, MLST, metadata curation, or genome retrieval.

From the repository root:

```bash
python3 -m pip install -r scripts/downstream_requirements.txt
python3 scripts/09_plasmid_context/plasmid_context.py
```

The default sequence-report location is the existing ignored
`local_archive/large_outputs/` tree. A relocated canonical download can be used
with `--sequence-reports-root PATH`. Use `--overwrite` to replace this module's
deterministic outputs. Tables are written to `results/tables/plasmid_context/`
and plots to `results/figures/plasmid_context/`.

The plasmid AMR-cargo heatmap colours detected genes by a documented display
class and groups class-adjacent columns with separators. The display classes are
aminoglycoside, macrolide/lincosamide, tetracycline,
vancomycin/glycopeptide, and other; the canonical detailed AMR class remains in
the hit-context table. A companion grouped bar chart reports accepted AMR-hit
burden for chromosome versus plasmid replicons by the same display classes, with
its underlying counts in
`amr_burden_by_gene_class_and_localisation_72.tsv`.

Localisation is evidence-limited: the exact AMRFinderPlus contig accession is
matched to the NCBI sequence report and classified from
`assignedMoleculeLocationType`. It is not inferred from plasmid prediction.
The optional `optional_similar_plasmids/` table is created only when ST6 proxy
plasmids have similar names, lengths (within 10%), and AMR cargo (Jaccard at
least 0.5). Such rows are screening candidates, not evidence of sequence-level
plasmid identity.
