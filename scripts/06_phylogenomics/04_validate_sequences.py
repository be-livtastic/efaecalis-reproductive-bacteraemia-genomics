#!/usr/bin/env python3
"""Validate extracted whole-CDS sequences and emit human-readable QC tables."""
from __future__ import annotations

import argparse
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq

from phylogeny_common import read_fasta, read_tsv, write_tsv

LOCI = ["gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL", "pyrC", "groEL", "recA"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-root", type=Path, required=True)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--length-tolerance", type=float, default=0.10)
    args = parser.parse_args()
    coordinates = read_tsv(args.coordinates)
    strand = {(r["assembly_accession"], r["canonical_gene"]): r["strand"] for r in coordinates}
    by_locus: dict[str, dict[str, str]] = {}
    for locus in LOCI:
        path = args.sequence_root / "by_gene" / f"{locus}_72_sequences.fasta"
        by_locus[locus] = read_fasta(path)
        if len(by_locus[locus]) != 72:
            raise SystemExit(f"{locus}: expected 72 sequences, found {len(by_locus[locus])}")

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
            translated = str(Seq(sequence).translate(table=11))
            internal = translated[:-1].count("*")
            ambiguous_nt = sum(base not in "ACGT" for base in sequence)
            ambiguous_aa = translated.count("X")
            deviation = abs(len(sequence) - medians[locus]) / medians[locus]
            flags = []
            if len(sequence) % 3: flags.append("LENGTH_NOT_DIVISIBLE_BY_3")
            if internal: flags.append("INTERNAL_STOP")
            if ambiguous_nt: flags.append("AMBIGUOUS_BASES")
            if deviation > args.length_tolerance: flags.append("LENGTH_OUTSIDE_TOLERANCE")
            if len(sequence) < medians[locus] * 0.8: flags.append("POSSIBLE_TRUNCATION")
            if len(sequence) % 3: flags.append("POSSIBLE_FRAMESHIFT")
            row = {"sample_id": sample, "canonical_gene": locus, "nucleotide_length": len(sequence),
                   "length_modulo_3": len(sequence) % 3, "start_codon": sequence[:3],
                   "terminal_stop_codon": sequence[-3:], "internal_stop_codons": internal,
                   "ambiguous_bases": ambiguous_nt, "translated_amino_acid_length": len(translated),
                   "ambiguous_residue_proportion": ambiguous_aa / max(1, len(translated)),
                   "median_locus_length": medians[locus], "length_deviation_proportion": deviation,
                   "possible_truncation": str("POSSIBLE_TRUNCATION" in flags).lower(),
                   "possible_frameshift": str("POSSIBLE_FRAMESHIFT" in flags).lower(),
                   "reverse_complemented": str(strand.get((sample, locus)) == "-").lower(),
                   "qc_status": "PASS" if not flags else ";".join(flags)}
            rows.append(row)
            if flags: failures.append(row)

    qc_fields = list(rows[0])
    write_tsv(args.qc_dir / "sequence_qc.tsv", rows, qc_fields)
    write_tsv(args.qc_dir / "locus_length_summary.tsv", length_rows, list(length_rows[0]))
    write_tsv(args.qc_dir / "sequence_qc_failures.tsv", failures, qc_fields)
    write_tsv(args.tables_dir / "gene_presence_qc.tsv", presence_rows, ["sample_id", *LOCI])
    hit_rows = [{"canonical_gene": locus, "expected": 72, "observed": len(by_locus[locus]),
                 "qc_failures": sum(r["canonical_gene"] == locus for r in failures)} for locus in LOCI]
    write_tsv(args.tables_dir / "gene_hit_summary.tsv", hit_rows, list(hit_rows[0]))
    print(f"Sequence QC: records={len(rows)} failures={len(failures)}")
    for row in length_rows:
        print(f"{row['gene']}: n={row['sequence_count']} min={row['minimum_length']} median={row['median_length']} max={row['maximum_length']}")
    if failures:
        print("REVIEW GATE: sequence QC failures remain; alignment is blocked.", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
