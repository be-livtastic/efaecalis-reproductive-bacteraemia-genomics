#!/usr/bin/env python3
"""BLAST translated loci against a small version-pinned trusted protein set."""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO

from phylogeny_common import read_fasta, read_tsv, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-root", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--findings", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=72)
    args = parser.parse_args()
    refs = read_tsv(args.references)
    disrupted = {(r["sample_id"], r["canonical_gene"]) for r in read_tsv(args.findings)
                 if r["code"] == "OBSERVED_DISRUPTED_LOCUS"}
    with tempfile.TemporaryDirectory(prefix="protein_qc.") as temporary_name:
        temporary = Path(temporary_name)
        query_path, subject_path = temporary / "queries.faa", temporary / "references.faa"
        with subject_path.open("w", encoding="utf-8") as handle:
            for ref in refs:
                source = args.reference_dir / f"{ref['accession']}.fasta"
                record = next(SeqIO.parse(source, "fasta"))
                handle.write(f">{ref['canonical_gene']}|{ref['accession']}\n{str(record.seq)}\n")
        with query_path.open("w", encoding="utf-8") as handle:
            for ref in refs:
                locus = ref["canonical_gene"]
                records = read_fasta(args.sequence_root / "by_gene" / f"{locus}_{args.expected_samples}_sequences.fasta")
                for sample, nt in records.items():
                    padded = nt + "N" * ((-len(nt)) % 3)
                    aa = str(Seq(padded).translate(table=11)).replace("*", "X")
                    handle.write(f">{sample}|{locus}\n{aa}\n")
        fields = "qseqid sseqid pident length mismatch gaps qlen slen qstart qend sstart send evalue bitscore".split()
        command = ["blastp", "-query", str(query_path), "-subject", str(subject_path), "-max_target_seqs", "20",
                   "-seg", "no", "-outfmt", "6 " + " ".join(fields)]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
    best = {}
    for line in completed.stdout.splitlines():
        values = dict(zip(fields, line.split("\t")))
        sample, locus = values["qseqid"].split("|", 1)
        subject_locus, accession = values["sseqid"].split("|", 1)
        if locus != subject_locus:
            continue
        key = (sample, locus)
        if key not in best or float(values["bitscore"]) > float(best[key]["bitscore"]):
            best[key] = values | {"sample_id": sample, "canonical_gene": locus,
                                  "reference_accession": accession}
    rows = []
    for ref in refs:
        locus = ref["canonical_gene"]
        for sample in read_fasta(args.sequence_root / "by_gene" / f"{locus}_{args.expected_samples}_sequences.fasta"):
            hit = best.get((sample, locus))
            if not hit:
                rows.append({"sample_id": sample, "canonical_gene": locus, "qc_status": "FAIL",
                             "failure_code": "NO_TRUSTED_PROTEIN_HIT"})
                continue
            aligned = int(hit["length"]); qlen = int(hit["qlen"]); slen = int(hit["slen"])
            qcov, scov, identity = aligned / qlen, aligned / slen, float(hit["pident"]) / 100
            expected_disruption = (sample, locus) in disrupted
            passed = identity >= 0.70 and qcov >= 0.90 and scov >= 0.90
            status = "EXPECTED_DISRUPTION" if expected_disruption else ("PASS" if passed else "FAIL")
            rows.append({"sample_id": sample, "canonical_gene": locus,
                         "reference_accession": hit["reference_accession"], "percent_identity": identity,
                         "query_coverage": qcov, "reference_coverage": scov,
                         "alignment_length_aa": aligned, "gaps": hit["gaps"], "mismatches": hit["mismatch"],
                         "query_length_aa": qlen, "reference_length_aa": slen,
                         "query_start": hit["qstart"], "query_end": hit["qend"],
                         "reference_start": hit["sstart"], "reference_end": hit["send"],
                         "evalue": hit["evalue"], "bitscore": hit["bitscore"], "qc_status": status,
                         "failure_code": "" if passed or expected_disruption else "PROTEIN_COVERAGE_OR_IDENTITY"})
    write_tsv(args.output, rows, list(rows[0]))
    counts = {status: sum(r["qc_status"] == status for r in rows) for status in {r["qc_status"] for r in rows}}
    print("Protein QC: " + " ".join(f"{key}={counts[key]}" for key in sorted(counts)))
    return 3 if counts.get("FAIL", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
