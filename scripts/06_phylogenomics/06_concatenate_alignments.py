#!/usr/bin/env python3
"""Validate locus alignments, calculate QC, and concatenate in fixed order."""
from __future__ import annotations

import argparse
import itertools
import statistics
from pathlib import Path

from phylogeny_common import atomic_text, read_fasta, read_tsv, write_tsv

LOCI = ["gdh", "gyd", "pstS", "gki", "xpt", "yqiL", "pyrC", "groEL", "recA"]


# --- Count variable and parsimony-informative sites ---
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


# --- Calculate pairwise divergence while ignoring gaps ---
def divergence(first: str, second: str) -> float:
    comparable = [(a, b) for a, b in zip(first, second) if a in "ACGT" and b in "ACGT"]
    return sum(a != b for a, b in comparable) / max(1, len(comparable))


# --- Validate, concatenate and partition locus alignments ---
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path,
                        default=Path("config/phylogeny_9_loci_qc_thresholds.tsv"))
    parser.add_argument("--expected-samples", type=int, default=72)
    parser.add_argument("--output-prefix", default="efaecalis_72_9_locus")
    parser.add_argument("--minimum-concatenated-coverage", type=float, default=0.95)
    args = parser.parse_args()
    thresholds = {r["canonical_gene"]: r for r in read_tsv(args.thresholds)}
    alignments, sample_set = {}, None
    summaries, missingness, pairwise_rows, outliers = [], [], [], []
    for locus in LOCI:
        records = read_fasta(args.alignment_dir / f"{locus}.aligned.fasta")
        if len(records) != args.expected_samples:
            raise SystemExit(f"{locus}: expected {args.expected_samples} records, found {len(records)}")
        lengths = {len(v) for v in records.values()}
        if len(lengths) != 1: raise SystemExit(f"{locus}: sequences have unequal aligned lengths")
        if sample_set is None: sample_set = set(records)
        if set(records) != sample_set: raise SystemExit(f"{locus}: identifiers differ from the first locus")
        length = lengths.pop()
        variable, informative = site_counts(records)
        pairwise = []
        for (sample_a, seq_a), (sample_b, seq_b) in itertools.combinations(records.items(), 2):
            value = divergence(seq_a, seq_b)
            pairwise.append(value)
            pairwise_rows.append({"gene": locus, "sample_a": sample_a, "sample_b": sample_b,
                                  "pairwise_divergence": value})
        maximum_pairwise = max(pairwise, default=0.0)
        median_pairwise = statistics.median(pairwise) if pairwise else 0.0
        summaries.append({"gene": locus, "sequence_count": args.expected_samples, "alignment_length": length,
                          "variable_sites": variable, "parsimony_informative_sites": informative,
                          "median_pairwise_divergence": median_pairwise,
                          "maximum_pairwise_divergence": maximum_pairwise})
        for sample, seq in records.items():
            gap = seq.count("-") / length
            ambiguous = sum(x not in "ACGT-" for x in seq) / length
            peers = [divergence(seq, other) for other_sample, other in records.items() if other_sample != sample]
            consensus_distance = statistics.median(peers)
            codes = []
            if gap > float(thresholds[locus]["maximum_gap_proportion"]): codes.append("EXCESSIVE_GAPS")
            if consensus_distance > float(thresholds[locus]["maximum_pairwise_divergence"]):
                codes.append("DIVERGENCE_OUTLIER")
            row = {"sample_id": sample, "gene": locus, "gap_proportion": gap,
                   "ambiguous_proportion": ambiguous, "median_divergence_from_peers": consensus_distance,
                   "qc_status": "FAIL" if codes else "PASS", "failure_codes": ";".join(codes)}
            missingness.append(row)
            if codes: outliers.append(row)
        alignments[locus] = records
        print(f"{locus}: sequences={args.expected_samples} length={length} variable={variable} informative={informative}")

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
    coverage_rows = []
    for sample, sequence in concatenated.items():
        covered = sum(base in "ACGT" for base in sequence)
        coverage = covered / len(sequence)
        coverage_rows.append({"sample_id": sample, "alignment_length": len(sequence),
                              "non_missing_characters": covered, "non_missing_proportion": coverage,
                              "minimum_required": args.minimum_concatenated_coverage,
                              "qc_status": "PASS" if coverage >= args.minimum_concatenated_coverage else "FAIL"})
    write_tsv(args.qc_dir / "alignment_summary.tsv", summaries, list(summaries[0]))
    write_tsv(args.qc_dir / "alignment_missingness.tsv", missingness, list(missingness[0]))
    write_tsv(args.qc_dir / "alignment_pairwise_divergence.tsv", pairwise_rows, list(pairwise_rows[0]))
    write_tsv(args.qc_dir / "alignment_qc_failures.tsv", outliers, list(missingness[0]))
    if outliers:
        raise SystemExit(f"Alignment QC failed for {len(outliers)} sample/locus records; concatenation is blocked")
    write_tsv(args.qc_dir / "concatenated_coverage_qc.tsv", coverage_rows, list(coverage_rows[0]))
    coverage_failures = [row for row in coverage_rows if row["qc_status"] == "FAIL"]
    if coverage_failures:
        raise SystemExit(f"Concatenated coverage failed for {len(coverage_failures)} genomes; concatenation is blocked")
    atomic_text(args.output_dir / f"{args.output_prefix}_alignment.fasta", fasta)
    atomic_text(args.output_dir / f"{args.output_prefix}_partitions.nex", partitions)
    write_tsv(args.output_dir / "partition_coordinates.tsv", coordinates, list(coordinates[0]))
    print(f"Concatenated alignment: taxa={args.expected_samples} total_length={start - 1}")
    print(f"Minimum per-genome non-missing proportion: {min(r['non_missing_proportion'] for r in coverage_rows):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
