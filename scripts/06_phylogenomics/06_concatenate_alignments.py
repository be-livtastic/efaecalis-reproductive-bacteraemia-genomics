#!/usr/bin/env python3
"""Validate locus alignments, calculate QC, and concatenate in fixed order."""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

from phylogeny_common import atomic_text, read_fasta, write_tsv

LOCI = ["gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL", "pyrC", "groEL", "recA"]


def site_counts(records: dict[str, str]) -> tuple[int, int]:
    variable = informative = 0
    for column in zip(*records.values()):
        states = {}
        for base in column:
            if base in "ACGT":
                states[base] = states.get(base, 0) + 1
        if len(states) > 1: variable += 1
        if sum(count >= 2 for count in states.values()) >= 2: informative += 1
    return variable, informative


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--qc-dir", type=Path, required=True)
    args = parser.parse_args()
    alignments, sample_set = {}, None
    summaries, missingness = [], []
    for locus in LOCI:
        records = read_fasta(args.alignment_dir / f"{locus}.aligned.fasta")
        if len(records) != 72: raise SystemExit(f"{locus}: expected 72 records, found {len(records)}")
        lengths = {len(v) for v in records.values()}
        if len(lengths) != 1: raise SystemExit(f"{locus}: sequences have unequal aligned lengths")
        if sample_set is None: sample_set = set(records)
        if set(records) != sample_set: raise SystemExit(f"{locus}: identifiers differ from the first locus")
        length = lengths.pop()
        variable, informative = site_counts(records)
        summaries.append({"gene": locus, "sequence_count": 72, "alignment_length": length,
                          "variable_sites": variable, "parsimony_informative_sites": informative})
        for sample, seq in records.items():
            missingness.append({"sample_id": sample, "gene": locus,
                "gap_proportion": seq.count("-") / length,
                "ambiguous_proportion": sum(x not in "ACGT-" for x in seq) / length})
        alignments[locus] = records
        print(f"{locus}: sequences=72 length={length} variable={variable} informative={informative}")

    coordinates, start = [], 1
    for summary in summaries:
        end = start + summary["alignment_length"] - 1
        coordinates.append({"gene": summary["gene"], "order": LOCI.index(summary["gene"]) + 1,
                            "start": start, "end": end, "aligned_length": summary["alignment_length"],
                            "cumulative_length": end})
        start = end + 1
    concatenated = {sample: "".join(alignments[locus][sample] for locus in LOCI) for sample in sorted(sample_set)}
    lengths = {len(seq) for seq in concatenated.values()}
    if len(lengths) != 1 or start - 1 != next(iter(lengths)):
        raise SystemExit("Concatenated alignment or partition coverage validation failed")
    fasta = "".join(f">{sample}\n{sequence}\n" for sample, sequence in concatenated.items())
    partitions = "#nexus\nbegin sets;\n" + "".join(
        f"  charset {r['gene']} = {r['start']}-{r['end']};\n" for r in coordinates) + "end;\n"
    atomic_text(args.output_dir / "efaecalis_72_10_locus_alignment.fasta", fasta)
    atomic_text(args.output_dir / "efaecalis_72_10_locus_partitions.nex", partitions)
    write_tsv(args.output_dir / "partition_coordinates.tsv", coordinates, list(coordinates[0]))
    write_tsv(args.qc_dir / "alignment_summary.tsv", summaries, list(summaries[0]))
    write_tsv(args.qc_dir / "alignment_missingness.tsv", missingness, list(missingness[0]))
    print(f"Concatenated alignment: taxa=72 total_length={start - 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
