# Repository cleanup audit

This directory records the cleanup without modifying scientific results.

- `pre_cleanup_manifest.tsv`: file size, SHA-256, Git status, classification
  and content source for the pre-cleanup working set.
- `post_cleanup_manifest.tsv`: the same inventory after approved cleanup.
- `cleanup_delta.tsv`: every added, removed, moved or modified file found by
  comparing the two manifests.
- `cleanup_summary.tsv`: counts and byte totals from the delta.
- `cleanup_path_mappings.tsv`: explicit old-to-new paths used to recognise
  intentional relocations.
- `pre_cleanup_git_status.txt` and `pre_cleanup_git_diff_stat.txt`: the dirty
  worktree state that existed before cleanup changes.
- `post_cleanup_git_status.txt` and `post_cleanup_git_diff_stat.txt`: the
  corresponding state immediately before staging and commit.
- `retention_decisions.tsv`: material deliberately retained despite apparent
  redundancy, plus the disclosed local-session exception.

The manifest builder excludes `.git/`, this audit-output directory and the
builder itself. The pre-cleanup inventory was formalised after PCA plotting
work began, so its five affected PCA paths were reconstructed from clean Git
`HEAD` blobs. The captured initial Git status showed those paths were clean;
the manifest identifies reconstructed entries as `content_source=git_HEAD`.
The plotting run also briefly normalised the companion metadata table to LF;
its original CRLF representation was restored and verified by the captured
frozen-bundle SHA-256. That manifest row records
`content_source=pre_cleanup_CRLF_restored`.

Directories do not have content SHA-256 values and therefore are represented
by their file records. Empty-directory removal is reported separately in the
final cleanup report.
