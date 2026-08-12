#!/usr/bin/env python3
"""Validate official MLST fragment identity against a pinned PubMLST snapshot."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from phylogeny_common import read_fasta, write_tsv

LOCI = ["gdh", "gyd", "pstS", "gki", "xpt", "yqiL"]


# --- Check MLST loci against the frozen PubMLST snapshot ---
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-root", type=Path, required=True)
    parser.add_argument("--pubmlst-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=72)
    args = parser.parse_args()
    fields = "qseqid sseqid pident length mismatch gaps qlen slen qstart qend sstart send evalue bitscore".split()
    rows = []
    for locus in LOCI:
        query = args.sequence_root / "by_gene" / f"{locus}_{args.expected_samples}_sequences.fasta"
        subject = args.pubmlst_dir / f"{locus}.alleles.fasta"
        command = ["blastn", "-task", "blastn", "-query", str(query), "-subject", str(subject),
                   "-max_target_seqs", "20", "-outfmt", "6 " + " ".join(fields)]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        best = {}
        for line in completed.stdout.splitlines():
            hit = dict(zip(fields, line.split("\t")))
            if hit["qseqid"] not in best or float(hit["bitscore"]) > float(best[hit["qseqid"]]["bitscore"]):
                best[hit["qseqid"]] = hit
        for sample in read_fasta(query):
            hit = best.get(sample)
            if not hit:
                rows.append({"sample_id": sample, "canonical_gene": locus, "qc_status": "FAIL",
                             "failure_code": "NO_PUBMLST_ALLELE_HIT"})
                continue
            aligned, qlen, slen = int(hit["length"]), int(hit["qlen"]), int(hit["slen"])
            identity = float(hit["pident"]) / 100
            qcov, allele_cov = aligned / qlen, aligned / slen
            passed = identity >= 0.95 and allele_cov >= 0.95
            rows.append({"sample_id": sample, "canonical_gene": locus, "best_allele": hit["sseqid"],
                         "percent_identity": identity, "whole_locus_query_coverage": qcov,
                         "pubmlst_allele_coverage": allele_cov, "alignment_length_nt": aligned,
                         "gaps": hit["gaps"], "mismatches": hit["mismatch"], "query_length_nt": qlen,
                         "allele_length_nt": slen, "query_start": hit["qstart"], "query_end": hit["qend"],
                         "allele_start": hit["sstart"], "allele_end": hit["send"], "evalue": hit["evalue"],
                         "bitscore": hit["bitscore"], "qc_status": "PASS" if passed else "FAIL",
                         "failure_code": "" if passed else "PUBMLST_IDENTITY_OR_ALLELE_COVERAGE"})
    write_tsv(args.output, rows, list(rows[0]))
    failures = sum(r["qc_status"] == "FAIL" for r in rows)
    print(f"PubMLST QC: records={len(rows)} pass={len(rows)-failures} fail={failures}")
    return 3 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
