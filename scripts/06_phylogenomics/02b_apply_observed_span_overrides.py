#!/usr/bin/env python3
"""Apply reviewed FNA span overrides without changing any nucleotide."""
from __future__ import annotations

import argparse
from pathlib import Path

from phylogeny_common import fasta_lengths, read_tsv, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--overrides", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows, overrides = read_tsv(args.coordinates), read_tsv(args.overrides)
    keyed = {(r["sample_id"], r["canonical_gene"]): r for r in overrides}
    if len(keyed) != len(overrides):
        raise SystemExit("Duplicate observed-span override")
    applied = set()
    for row in rows:
        key = (row["sample_id"], row["canonical_gene"])
        if key not in keyed:
            continue
        override = keyed[key]
        start, end = int(override["start"]), int(override["end"])
        if end - start + 1 != int(override["observed_length_nt"]):
            raise SystemExit(f"Override length mismatch for {key}")
        lengths = fasta_lengths(Path(row["fna_path"]))
        if override["contig"] not in lengths or start < 1 or end > lengths[override["contig"]]:
            raise SystemExit(f"Override outside FNA bounds for {key}")
        row.update(contig=override["contig"], start=str(start), end=str(end), strand=override["strand"],
                   sequence_length_expected=override["observed_length_nt"],
                   selection_source="reviewed_observed_fna_span")
        applied.add(key)
    if applied != set(keyed):
        raise SystemExit(f"Unapplied observed-span overrides: {sorted(set(keyed) - applied)}")
    write_tsv(args.output, rows, list(rows[0]))
    print(f"Observed-span policy: records={len(rows)} overrides_applied={len(applied)} imputed_nucleotides=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
