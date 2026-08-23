# Repository organisation manifest (2026-08-19)

## Scope

Local organisation/reconciliation only. No biological analyses were run.

## Archive discovery and verification

Two checksum-verified archives with the same base name were found:

1. `C:/Users/bless/Downloads/e_faecalis_disseration_project/efaecalis_core_genome_results_20260819.tar.gz`  
   - Size: `3,337,308` bytes  
   - SHA-256 (calculated): `36db6b806c0aadc7ec1428567a81da46677deef061f21e03f96ea433ff2df110`  
   - SHA-256 (recorded): `36db6b806c0aadc7ec1428567a81da46677deef061f21e03f96ea433ff2df110`  
   - Verification: `PASS`
2. `C:/Users/bless/Downloads/e_faecalis_disseration_project/backups/core_genome_20260819/efaecalis_core_genome_results_20260819.tar.gz`  
   - Size: `3,380,212` bytes  
   - SHA-256 (calculated): `758b3b52dcc779ccd1579340e06ee6ee8e1e15686619a5bb7a73b29f491e304d`  
   - SHA-256 (recorded): `758b3b52dcc779ccd1579340e06ee6ee8e1e15686619a5bb7a73b29f491e304d`  
   - Verification: `PASS`

Both were extracted only to staging:

- `backups/core_genome_20260819/extracted_staging/root_archive/`
- `backups/core_genome_20260819/extracted_staging/backup_archive/`

Reconciliation inventories:

- `backups/core_genome_20260819/extracted_staging/reconciliation_root_archive.tsv`
- `backups/core_genome_20260819/extracted_staging/reconciliation_backup_archive.tsv`

## Reconciliation outcome (hash-based)

Primary (`root_archive`) classification counts:

- A already present + identical: 20
- B already present + archive newer: 1
- C already present + different: 25
- D genuinely new: 26

Secondary (`backup_archive`) classification counts:

- A already present + identical: 28
- B already present + archive newer: 1
- C already present + different: 41
- D genuinely new: 26

No C/B file was silently overwritten.

## Files moved/organised

### Added from staging to analysis working areas

- `analysis/core_genome/iqtree_primary_runs/full_core_72_gtrg4_20260819T042202Z.*` (complete working run set, including `.ckp.gz`, `.pid`, `.launcher.log`, `.bionj`)
- `analysis/core_genome/qc/output/*` (detailed QC/audit outputs)

### Superseded developmental QC attempts archived (not deleted)

- `analysis/core_genome/qc/archive/core_tree_duplicate_sequence_pairs_20260819T044833Z.tsv`
- `analysis/core_genome/qc/archive/core_tree_terminal_branch_lengths_20260819T044833Z.tsv`

### Download/archive hygiene

- Root download copy moved to:  
  `backups/core_genome_20260819/source_archives/root_download/efaecalis_core_genome_results_20260819_rootcopy.tar.gz`  
  `backups/core_genome_20260819/source_archives/root_download/efaecalis_core_genome_results_20260819_rootcopy.tar.gz.sha256`
- Zero-byte transfer fragment moved to:  
  `backups/repository_cleanup_20260819/GCF_053117695.1.zero_byte_transfer_fragment`

## Duplicate/variant files identified

- Two different, checksum-valid archives share the same name (`efaecalis_core_genome_results_20260819.tar.gz`) but have different content hashes and file counts (80 vs 106 archive entries).
- Multiple files with matching intended roles but differing hashes were detected (category C/B), including:
  - `results/core_genome/iqtree_primary/full_core_72_gtrg4_20260819T042202Z.*`
  - several `results/core_genome/qc/*` timestamped files
  - several `results/figures/core_genome_phylogeny/*.svg`
  - `scripts/07_core_genome/08_visualise_core_genome_tree.R`
  - `scripts/07_core_genome/README_core_genome.md`
  - `.devcontainer/devcontainer.json`
  - `environment/core_genome.yml`

These were preserved for manual scientific/provenance review; no overwrite applied.

## Git-tracked vs local-only separation

`.gitignore` was updated so local-only working artefacts remain untracked:

- `analysis/core_genome/iqtree_primary_runs/**`
- `analysis/core_genome/qc/output/**`
- `analysis/core_genome/qc/archive/**`
- backup reconciliation/staging folders under `backups/core_genome_20260819/`
- cleanup quarantine folder `backups/repository_cleanup_20260819/`

Compact reportable outputs remain under tracked `results/core_genome/`, `results/figures/core_genome_phylogeny/`, and `results/tables/core_genome_phylogeny/`.

## Documentation/path updates applied

- `scripts/07_core_genome/README_core_genome.md` (core-genome output path wording updated)
- `docs/workflow/analysis_workflow.md` (core-genome section updated to current organised locations)
- `README.md` (core-genome workflow status/path summary updated)
- `scripts/07_core_genome/compare_core_phylogenies.R` (core-tree input and figure output paths aligned)

## Final core-genome structure (concise)

```text
analysis/core_genome/
  panaroo_strict_core95/
  iqtree_primary_runs/
  qc/
    output/
    archive/

results/core_genome/
  iqtree_primary/
  qc/
  trees/

results/figures/core_genome_phylogeny/
results/tables/core_genome/
results/tables/core_genome_phylogeny/

scripts/07_core_genome/
environment/
backups/core_genome_20260819/
backups/repository_cleanup_20260819/
```
