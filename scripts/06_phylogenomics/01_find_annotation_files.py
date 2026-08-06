#!/usr/bin/env python3
"""Discover and validate a policy-defined set of local Prokka annotations."""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

from phylogeny_common import accession_from_text, fasta_lengths, read_tsv, write_tsv

FIELDS = ["sample_id", "assembly_accession", "dataset_group", "prokka_directory",
          "gff_path", "fna_path", "ffn_path", "faa_path", "gff_exists",
          "fna_exists", "basename_match", "notes"]


def gff_contigs(path: Path) -> set[str]:
    contigs: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 9:
                contigs.add(parts[0])
    return contigs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prokka-root", type=Path, required=True)
    parser.add_argument("--accessions", type=Path, default=Path("data/accession_lists/selected_72_accessions.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=72)
    parser.add_argument("--exclude-accessions", type=Path)
    args = parser.parse_args()

    expected_rows = read_tsv(args.accessions)
    excluded = set()
    if args.exclude_accessions:
        excluded = {r["assembly_accession"].strip() for r in read_tsv(args.exclude_accessions)}
        expected_rows = [r for r in expected_rows if r["assembly_accession"].strip() not in excluded]
    expected = {r["assembly_accession"].strip(): r["dataset_category"].strip() for r in expected_rows}
    if len(expected) != args.expected or len(expected_rows) != args.expected:
        raise SystemExit(f"Expected {args.expected} unique accessions; found {len(expected)} unique in {len(expected_rows)} rows")

    discovered: dict[str, list[Path]] = defaultdict(list)
    for directory in sorted(p for p in args.prokka_root.rglob("*") if p.is_dir()):
        accession = accession_from_text(directory.name)
        if accession and list(directory.glob("*.gff")):
            discovered[accession].append(directory)

    rows: list[dict] = []
    unmatched: list[str] = []
    valid = 0
    for accession in sorted(set(expected) | (set(discovered) - excluded)):
        directories = discovered.get(accession, [])
        if accession not in expected:
            unmatched.append(accession)
        if not directories:
            rows.append({"sample_id": accession, "assembly_accession": accession,
                         "dataset_group": expected.get(accession, ""), "notes": "MISSING_ANNOTATION_DIRECTORY"})
            continue
        for directory in directories:
            gffs, fnas = sorted(directory.glob("*.gff")), sorted(directory.glob("*.fna"))
            ffns, faas = sorted(directory.glob("*.ffn")), sorted(directory.glob("*.faa"))
            notes: list[str] = []
            if len(directories) != 1: notes.append("DUPLICATE_ANNOTATION_DIRECTORY")
            if len(gffs) != 1: notes.append(f"GFF_COUNT={len(gffs)}")
            if len(fnas) != 1: notes.append(f"FNA_COUNT={len(fnas)}")
            basename_match = bool(len(gffs) == len(fnas) == 1 and gffs[0].stem == fnas[0].stem)
            if len(gffs) == len(fnas) == 1:
                absent = sorted(gff_contigs(gffs[0]) - set(fasta_lengths(fnas[0])))
                if absent: notes.append(f"GFF_CONTIGS_ABSENT_FROM_FNA={len(absent)}")
            if accession not in expected: notes.append("UNEXPECTED_ACCESSION")
            if not notes and basename_match:
                valid += 1
            rows.append({
                "sample_id": accession, "assembly_accession": accession,
                "dataset_group": expected.get(accession, ""),
                "prokka_directory": str(directory.resolve()),
                "gff_path": str(gffs[0].resolve()) if len(gffs) == 1 else "",
                "fna_path": str(fnas[0].resolve()) if len(fnas) == 1 else "",
                "ffn_path": str(ffns[0].resolve()) if len(ffns) == 1 else "",
                "faa_path": str(faas[0].resolve()) if len(faas) == 1 else "",
                "gff_exists": str(len(gffs) == 1).lower(),
                "fna_exists": str(len(fnas) == 1).lower(),
                "basename_match": str(basename_match).lower(),
                "notes": ";".join(notes)
            })

    write_tsv(args.output, rows, FIELDS)
    group_counts = Counter(r["dataset_group"].casefold() for r in rows if not r["notes"])
    missing = sorted(set(expected) - set(discovered))
    duplicates = sorted(k for k, v in discovered.items() if len(v) != 1)
    print(f"Annotation directories found: {sum(map(len, discovered.values()))}")
    print(f"Valid matching GFF/FNA pairs: {valid}/{args.expected}")
    print(f"Missing accessions: {','.join(missing) or 'none'}")
    print(f"Duplicate accessions: {','.join(duplicates) or 'none'}")
    print(f"Unexpected accessions: {','.join(unmatched) or 'none'}")
    print(f"Counts by group: reproductive={group_counts['reproductive']}; bacteraemia={group_counts['bacteraemia']}")
    if valid != args.expected or missing or duplicates or unmatched:
        print(f"ERROR: input discovery did not satisfy the {args.expected}-genome invariant", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
