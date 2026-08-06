#!/usr/bin/env python3
"""Enumerate GFF candidates and select only unambiguous, high-confidence CDS hits."""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

from phylogeny_common import fasta_lengths, normalize, parse_attributes, read_tsv, write_tsv

CANDIDATE_FIELDS = ["sample_id", "assembly_accession", "dataset_group", "canonical_gene",
    "matched_attribute", "matched_value", "feature_type", "contig", "start", "end", "strand",
    "phase", "locus_tag", "feature_id", "gene_name", "product", "sequence_length_expected",
    "candidate_rank", "match_reason", "qc_status"]
SELECTED_FIELDS = CANDIDATE_FIELDS + ["gff_path", "fna_path", "selection_source"]


def exact_evidence(attrs: dict[str, str], config: dict[str, str]) -> list[tuple[int, str, str, str]]:
    aliases = {normalize(x) for x in config["aliases"].split(";") if x}
    products = {normalize(x) for x in config["search_terms"].split(";") if x}
    evidence: list[tuple[int, str, str, str]] = []
    for key in ("gene", "Name"):
        value = attrs.get(key, "")
        base = normalize(value).removesuffix(" 1").removesuffix(" 2")
        if base in aliases:
            evidence.append((1, key, value, "exact_symbol_or_alias"))
    product = attrs.get("product", "")
    if normalize(product) in products:
        evidence.append((2, "product", product, "exact_product"))
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--genes", type=Path, default=Path("config/phylogeny_9_loci_genes.tsv"))
    parser.add_argument("--overrides", type=Path)
    parser.add_argument("--all-candidates", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=72)
    parser.add_argument("--expected-loci", type=int, default=9)
    args = parser.parse_args()
    manifest, genes = read_tsv(args.manifest), read_tsv(args.genes)
    expected_records = len(manifest) * len(genes)
    if len(manifest) != args.expected_samples or len(genes) != args.expected_loci:
        raise SystemExit(f"Expected {args.expected_samples} manifest rows and {args.expected_loci} genes; found {len(manifest)} and {len(genes)}")

    overrides = {}
    if args.overrides and args.overrides.exists():
        for row in read_tsv(args.overrides):
            key = (row["sample_id"].strip(), row["canonical_gene"].strip())
            if key in overrides: raise SystemExit(f"Duplicate override: {key}")
            overrides[key] = row

    all_rows: list[dict] = []
    candidates: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for sample in manifest:
        fna_lengths = fasta_lengths(Path(sample["fna_path"]))
        with Path(sample["gff_path"]).open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"): continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9: continue
                contig, _, feature_type, start, end, _, strand, phase, attribute_text = parts
                if feature_type != "CDS": continue
                attrs = parse_attributes(attribute_text)
                for gene in genes:
                    evidence = exact_evidence(attrs, gene)
                    if not evidence: continue
                    rank, matched_attribute, matched_value, reason = sorted(evidence)[0]
                    flags: list[str] = []
                    if contig not in fna_lengths: flags.append("CONTIG_ABSENT")
                    if int(start) < 1 or int(end) > fna_lengths.get(contig, 0): flags.append("COORDINATE_OUT_OF_RANGE")
                    text = normalize(attribute_text)
                    if "pseudo true" in text or "pseudogene" in text: flags.append("PSEUDOGENE")
                    row = {
                        "sample_id": sample["sample_id"], "assembly_accession": sample["assembly_accession"],
                        "dataset_group": sample["dataset_group"], "canonical_gene": gene["canonical_gene"],
                        "matched_attribute": matched_attribute, "matched_value": matched_value,
                        "feature_type": feature_type, "contig": contig, "start": start, "end": end,
                        "strand": strand, "phase": phase, "locus_tag": attrs.get("locus_tag", ""),
                        "feature_id": attrs.get("ID", ""), "gene_name": attrs.get("gene", attrs.get("Name", "")),
                        "product": attrs.get("product", ""), "sequence_length_expected": int(end)-int(start)+1,
                        "candidate_rank": rank, "match_reason": reason,
                        "qc_status": ";".join(flags) if flags else "CANDIDATE"
                    }
                    all_rows.append(row)
                    candidates[(sample["sample_id"], gene["canonical_gene"])].append(row)

    selected: list[dict] = []
    unresolved: list[tuple[str, str, int]] = []
    for sample in manifest:
        for gene in genes:
            key = (sample["sample_id"], gene["canonical_gene"])
            options = candidates.get(key, [])
            override = overrides.get(key)
            chosen = []
            source = ""
            if override:
                chosen = [r for r in options if r["feature_id"] == override["selected_feature_id"]]
                source = "reviewed_override"
                if len(chosen) != 1:
                    raise SystemExit(f"Override for {key} did not select exactly one candidate")
            elif len(options) == 1 and not options[0]["qc_status"].replace("CANDIDATE", ""):
                chosen = options
                source = "unique_exact_annotation"
            if len(chosen) == 1:
                row = dict(chosen[0], gff_path=sample["gff_path"], fna_path=sample["fna_path"],
                           selection_source=source, qc_status="SELECTED")
                selected.append(row)
            else:
                unresolved.append((key[0], key[1], len(options)))

    write_tsv(args.all_candidates, all_rows, CANDIDATE_FIELDS)
    write_tsv(args.selected, selected, SELECTED_FIELDS)
    counts = Counter(r["canonical_gene"] for r in all_rows)
    print("Candidate counts: " + "; ".join(f"{g['canonical_gene']}={counts[g['canonical_gene']]}" for g in genes))
    print(f"Selected records: {len(selected)}/{expected_records}")
    print(f"Unresolved sample/locus combinations: {len(unresolved)}")
    for sample, gene, count in unresolved[:30]:
        print(f"REVIEW_REQUIRED\t{sample}\t{gene}\tcandidates={count}")
    if unresolved or len(selected) != expected_records:
        print("REVIEW GATE: downstream extraction is blocked. Copy the override example, document decisions, and rerun.", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
