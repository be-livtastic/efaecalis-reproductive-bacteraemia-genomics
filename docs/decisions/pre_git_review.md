# Pre-Git review

## Stop point

The project has been reorganised locally. Nothing has been staged, committed,
renamed on GitHub or pushed. The existing repository and its history remain in
place.

The requested next repository name is
`efaecalis-reproductive-bacteraemia-genomics`. Renaming the current repository
will preserve its existing history, including historical raw genome objects.
Remote privacy has not yet been verified.

## Proposed commit

The exact proposed additions and tracked removals are recorded in
`proposed_commit_files.tsv`.

- Proposed additions: 55 files after including this report.
- Proposed tracked removals: 97 legacy paths.
- Approximate added content: 3.2 MB.
- Largest proposed file: approximately 0.55 MB.
- Files larger than 5 MB: none.
- Files approaching GitHub's 100 MB per-file limit: none.

The tracked removals remove raw genomes and legacy generated data from the next
tree, but do not erase them from existing Git history.

## Deliberately excluded

- `local_archive/`: complete legacy material, current phylogeny work, original
  scripts/configuration, R session state and original workbook.
- Original Excel workbook: local-only and explicitly ignored.
- Raw NCBI genome downloads and duplicate FASTA collections: public,
  regenerable and large.
- Bulk Prokka outputs: locally retained; only compact summary/manifests proposed.
- Panaroo and IQ-TREE intermediates: regenerable and/or not yet final.
- Current phylogeny scripts and outputs: held locally pending revision.
- Caches, environments, logs, temporary documents and private configuration.

The complete exclusion rules are in the repository-root `.gitignore`.

## Secret and privacy scan

- Repository-facing working tree: no secret indicators and no personal absolute
  paths detected.
- Original Excel workbook: excluded without publication.
- Existing history: no actionable credential was identified by the pattern
  scan. AWS-style `AKIA` matches occurred inside public protein FASTA sequences
  because amino-acid sequences use the same uppercase alphabet; these are
  biological-sequence false positives.
- Existing history contains a personal installation path in the historical
  AMRFinder version file. The proposed replacement version record omits it.

Before staging, the scan should be repeated over the exact proposed set. After
staging, it should be repeated over the index. Repository privacy must be
verified before the next push.

## Validation completed

- Archive move checksums verified against original-path aggregates.
- Accession manifest validated: 14 reproductive, 58 bacteraemia, 72 unique.
- Curated metadata validated: 72 rows and 72 unique accessions.
- Curated metadata contains no path, workbook-row, email, student-number or
  signature columns.
- Bash, Python and R syntax checks passed.
- Retrieval and helper scripts passed missing-input/existing-output safe-failure
  tests.
- No workbook, raw genome, bulk Prokka file or current phylogeny material occurs
  in the proposed additions.

