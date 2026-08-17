# Core-genome environment decision

Date: 2026-08-17

The existing `efaecalis_phylogeny` environment was inspected first. It already contained compatible MAFFT and IQ-TREE installations but did not contain Panaroo or `snp-dists`. A non-mutating Mamba solver dry run found that adding Panaroo 1.8.0 would require changing the validated Python 3.12 dependency set because the available `intbitset` build requires Python 3.11 or earlier. The existing environment was therefore left unchanged.

The separate `efaecalis_core_genome` environment was created from `core_genome.yml`. Direct version checks after creation confirmed Python 3.11.15, Panaroo 1.8.0, `snp-dists` 1.2.0, MAFFT 7.526, IQ-TREE 3.1.2, Biopython 1.87, pandas 3.0.5 and NumPy 2.3.5. This environment is used only for the core-genome stages; downstream R summaries use the established project environment.
