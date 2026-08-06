#!/usr/bin/env python3
"""Validate whole-CDS sequences using locus-specific PASS/REVIEW/FAIL rules."""
from __future__ import annotations

import argparse
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq

from phylogeny_common import read_fasta, read_tsv, write_tsv

LOCI = ["gdh", "gyd", "pstS", "gki", "xpt", "yqiL", "pyrC", "groEL", "recA"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-root", type=Path, required=True)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=72)
    parser.add_argument("--thresholds", type=Path,
                        default=Path("config/phylogeny_9_loci_qc_thresholds.tsv"))
    parser.add_argument("--reviewed-findings", type=Path,
                        default=Path("config/phylogeny_9_loci_reviewed_findings.tsv"))
    args = parser.parse_args()
    coordinates = read_tsv(args.coordinates)
    threshold_rows = read_tsv(args.thresholds)
    thresholds = {row["canonical_gene"]: row for row in threshold_rows}
    if set(thresholds) != set(LOCI):
        raise SystemExit(f"QC threshold loci differ from active loci: {sorted(thresholds)}")
    findings: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    if args.reviewed_findings.exists():
        for finding in read_tsv(args.reviewed_findings):
            findings[(finding["sample_id"], finding["canonical_gene"])].append(finding)
    strand = {(r["assembly_accession"], r["canonical_gene"]): r["strand"] for r in coordinates}
    by_locus: dict[str, dict[str, str]] = {}
    for locus in LOCI:
        path = args.sequence_root / "by_gene" / f"{locus}_{args.expected_samples}_sequences.fasta"
        by_locus[locus] = read_fasta(path)
        if len(by_locus[locus]) != args.expected_samples:
            raise SystemExit(f"{locus}: expected {args.expected_samples} sequences, found {len(by_locus[locus])}")

    medians = {locus: statistics.median(map(len, records.values())) for locus, records in by_locus.items()}
    rows, failures, length_rows = [], [], []
    all_samples = sorted(set().union(*(set(records) for records in by_locus.values())))
    presence_rows = []
    for sample in all_samples:
        presence_rows.append({"sample_id": sample, **{locus: int(sample in by_locus[locus]) for locus in LOCI}})
    for locus, records in by_locus.items():
        lengths = [len(seq) for seq in records.values()]
        length_rows.append({"gene": locus, "sequence_count": len(records), "minimum_length": min(lengths),
                            "median_length": medians[locus], "maximum_length": max(lengths)})
        for sample, sequence in records.items():
            record_findings = findings.get((sample, locus), [])
            observed_disruption = any(f["code"] == "OBSERVED_DISRUPTED_LOCUS" for f in record_findings)
            translation_input = sequence + "N" * ((-len(sequence)) % 3)
            translated = str(Seq(translation_input).translate(table=11))
            internal = translated[:-1].count("*")
            ambiguous_nt = sum(base not in "ACGT" for base in sequence)
            ambiguous_aa = translated.count("X")
            deviation = abs(len(sequence) - medians[locus]) / medians[locus]
            rule = thresholds[locus]
            length = len(sequence)
            pass_min, pass_max = int(rule["pass_min_nt"]), int(rule["pass_max_nt"])
            review_min, review_max = int(rule["review_min_nt"]), int(rule["review_max_nt"])
            failure_codes, review_codes, information_codes = [], [], []
            if length % 3:
                target = information_codes if observed_disruption else failure_codes
                target.extend(["LENGTH_NOT_DIVISIBLE_BY_3", "OBSERVED_FRAMESHIFT"])
            if internal:
                (information_codes if observed_disruption else failure_codes).append("INTERNAL_STOP")
            if ambiguous_nt: failure_codes.append("AMBIGUOUS_BASES")
            if not (review_min <= length <= review_max):
                failure_codes.append("LENGTH_OUTLIER")
                if length < review_min: failure_codes.append("POSSIBLE_TRUNCATION")
            elif not (pass_min <= length <= pass_max):
                review_codes.append("LENGTH_REVIEW")
            if sequence[-3:] not in {"TAA", "TAG", "TGA"}:
                (information_codes if observed_disruption else failure_codes).append("MISSING_TERMINAL_STOP")
            for finding in record_findings:
                target = (failure_codes if finding["severity"] == "FAIL" else
                          review_codes if finding["severity"] == "REVIEW" else information_codes)
                target.append(finding["code"])
            failure_codes = list(dict.fromkeys(failure_codes))
            review_codes = list(dict.fromkeys(review_codes))
            status = "FAIL" if failure_codes else ("REVIEW" if review_codes else "PASS")
            row = {"sample_id": sample, "canonical_gene": locus, "nucleotide_length": len(sequence),
                   "length_modulo_3": len(sequence) % 3, "start_codon": sequence[:3],
                   "terminal_stop_codon": sequence[-3:], "internal_stop_codons": internal,
                   "ambiguous_bases": ambiguous_nt, "translated_amino_acid_length": len(translated),
                   "ambiguous_residue_proportion": ambiguous_aa / max(1, len(translated)),
                   "median_locus_length": medians[locus], "length_deviation_proportion": deviation,
                   "expected_length_nt": rule["expected_length_nt"],
                   "length_pass_range_nt": f"{pass_min}-{pass_max}",
                   "length_review_range_nt": f"{review_min}-{review_max}",
                   "possible_truncation": str("POSSIBLE_TRUNCATION" in failure_codes).lower(),
                   "possible_frameshift": str("POSSIBLE_FRAMESHIFT" in failure_codes).lower(),
                   "reverse_complemented": str(strand.get((sample, locus)) == "-").lower(),
                   "translation_policy": "OBSERVED_DISRUPTED_LOCUS" if observed_disruption else "COMPLETE_CDS",
                   "qc_status": status, "failure_codes": ";".join(failure_codes),
                   "review_codes": ";".join(review_codes),
                   "information_codes": ";".join(information_codes)}
            rows.append(row)
            if status != "PASS": failures.append(row)

    qc_fields = list(rows[0])
    write_tsv(args.qc_dir / "sequence_qc.tsv", rows, qc_fields)
    write_tsv(args.qc_dir / "locus_length_summary.tsv", length_rows, list(length_rows[0]))
    write_tsv(args.qc_dir / "sequence_qc_failures.tsv", failures, qc_fields)
    write_tsv(args.tables_dir / "gene_presence_qc.tsv", presence_rows, ["sample_id", *LOCI])
    hit_rows = [{"canonical_gene": locus, "expected": args.expected_samples, "observed": len(by_locus[locus]),
                 "qc_reviews": sum(r["canonical_gene"] == locus and r["qc_status"] == "REVIEW" for r in rows),
                 "qc_failures": sum(r["canonical_gene"] == locus and r["qc_status"] == "FAIL" for r in rows)}
                for locus in LOCI]
    write_tsv(args.tables_dir / "gene_hit_summary.tsv", hit_rows, list(hit_rows[0]))
    assembly_rows = []
    for sample in all_samples:
        sample_rows = [r for r in rows if r["sample_id"] == sample]
        reviews = [r["canonical_gene"] for r in sample_rows if r["qc_status"] == "REVIEW"]
        failed = [r["canonical_gene"] for r in sample_rows if r["qc_status"] == "FAIL"]
        assembly_rows.append({"sample_id": sample, "loci_evaluated": len(sample_rows),
                              "pass_count": sum(r["qc_status"] == "PASS" for r in sample_rows),
                              "review_count": len(reviews), "fail_count": len(failed),
                              "review_loci": ";".join(reviews), "failed_loci": ";".join(failed),
                              "assembly_qc_status": "FAIL" if failed else ("REVIEW" if reviews else "PASS"),
                              "multiple_suspicious_loci": str(len(reviews) + len(failed) > 1).lower()})
    write_tsv(args.qc_dir / "assembly_qc_summary.tsv", assembly_rows, list(assembly_rows[0]))
    write_tsv(args.tables_dir / "assembly_qc_summary.tsv", assembly_rows, list(assembly_rows[0]))
    review_count = sum(r["qc_status"] == "REVIEW" for r in rows)
    fail_count = sum(r["qc_status"] == "FAIL" for r in rows)
    print(f"Sequence QC: records={len(rows)} pass={len(rows)-len(failures)} review={review_count} fail={fail_count}")
    for row in length_rows:
        print(f"{row['gene']}: n={row['sequence_count']} min={row['minimum_length']} median={row['median_length']} max={row['maximum_length']}")
    if failures:
        print("REVIEW GATE: sequence QC reviews or failures remain; alignment is blocked.", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
