#!/usr/bin/env python3
"""PCA of observed alleles in an aligned multilocus nucleotide matrix."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from Bio import SeqIO


# --- Command-line interface ---
def arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--alignment", required=True, type=Path)
    p.add_argument("--partitions", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--min-allele-count", type=int, default=2)
    p.add_argument("--components", type=int, default=10)
    return p.parse_args()


# --- Small deterministic TSV writer ---
def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


# --- Encode observed alleles and run the alignment PCA ---
def main() -> int:
    a = arguments()
    records = list(SeqIO.parse(a.alignment, "fasta"))
    if not records or len({len(x.seq) for x in records}) != 1:
        raise SystemExit("Alignment is empty or not rectangular")
    ids = [x.id for x in records]
    if len(ids) != len(set(ids)):
        raise SystemExit("Duplicate FASTA identifiers")
    seq = np.array([list(str(x.seq).upper()) for x in records], dtype="U1")
    with a.partitions.open(newline="", encoding="utf-8-sig") as handle:
        partitions = list(csv.DictReader(handle, delimiter="\t"))
    def locus_for(position: int) -> str:
        hits = [r["gene"] for r in partitions if int(r["start"]) <= position <= int(r["end"])]
        if len(hits) != 1:
            raise SystemExit(f"Partition coverage error at alignment position {position}")
        return hits[0]
    valid = set("ACGT")
    features, info = [], []
    variable_sites = singleton_sites = excluded_missing = 0
    for j in range(seq.shape[1]):
        calls = [x for x in seq[:, j] if x in valid]
        counts = {x: calls.count(x) for x in sorted(set(calls))}
        if len(counts) < 2:
            continue
        variable_sites += 1
        ordered = sorted(counts, key=lambda x: (-counts[x], x))
        alternatives = [x for x in ordered[1:] if counts[x] >= a.min_allele_count]
        if not alternatives:
            singleton_sites += 1
            continue
        if not calls:
            excluded_missing += 1
            continue
        for allele in alternatives:
            col = np.array([1.0 if x == allele else 0.0 if x in valid else np.nan for x in seq[:, j]])
            mean = np.nanmean(col)
            centred = col - mean
            centred[np.isnan(centred)] = 0.0  # mean value after centring; analytical only
            sd = np.sqrt(np.mean(centred**2))
            if sd == 0:
                continue
            features.append(centred / sd)
            info.append({"locus": locus_for(j + 1), "alignment_position_1based": j + 1, "major_allele": ordered[0],
                         "encoded_allele": allele, "encoded_allele_count": counts[allele],
                         "called_genomes": len(calls), "missing_genomes": len(records) - len(calls)})
    if not features:
        raise SystemExit("No PCA features passed filtering")
    matrix = np.column_stack(features)
    u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
    eigen = singular**2 / (len(records) - 1)
    explained = eigen / eigen.sum()
    ncomp = min(a.components, len(singular))
    scores = u[:, :ncomp] * singular[:ncomp]
    loadings = vt[:ncomp, :].T

    with a.metadata.open(newline="", encoding="utf-8-sig") as handle:
        metadata = {r["assembly_accession"]: r for r in csv.DictReader(handle)}
    if set(ids) != set(metadata):
        raise SystemExit(f"Metadata mismatch: tree-only={len(set(ids)-set(metadata))}, metadata-only={len(set(metadata)-set(ids))}")
    a.output_dir.mkdir(parents=True, exist_ok=True)
    score_fields = [
        "assembly_accession", "biosample_accession", "strain",
        "reproductive_bacteraemia_category"
    ] + [f"PC{i}" for i in range(1, ncomp + 1)]
    score_rows = []
    for i, sample in enumerate(ids):
        strain = metadata[sample].get("strain", "").strip() or sample
        row = {
            "assembly_accession": sample,
            "biosample_accession": metadata[sample].get("biosample_accession", ""),
            "strain": strain,
            "reproductive_bacteraemia_category": metadata[sample]["reproductive_bacteraemia_category"],
        }
        row.update({f"PC{k+1}": scores[i, k] for k in range(ncomp)})
        score_rows.append(row)
    write_tsv(a.output_dir / "pca_scores.tsv", score_fields, score_rows)
    variance_rows = [{"component": f"PC{i+1}", "eigenvalue": eigen[i],
                      "explained_variance_proportion": explained[i],
                      "cumulative_variance_proportion": explained[: i + 1].sum()} for i in range(ncomp)]
    write_tsv(a.output_dir / "pca_explained_variance.tsv", list(variance_rows[0]), variance_rows)
    for i, row in enumerate(info):
        for k in range(ncomp): row[f"PC{k+1}_loading"] = loadings[i, k]
    write_tsv(a.output_dir / "pca_feature_loadings.tsv", list(info[0]), info)
    qc = [{"genomes": len(records), "alignment_length_nt": seq.shape[1], "variable_sites": variable_sites,
           "sites_excluded_below_minimum_allele_count": singleton_sites,
           "encoded_features": matrix.shape[1], "minimum_allele_count": a.min_allele_count,
           "ambiguous_or_gap_calls": int(np.sum(~np.isin(seq, list(valid)))),
           "missing_value_numeric_policy": "feature mean (zero after centring); source alignment unchanged"}]
    write_tsv(a.output_dir / "pca_input_qc.tsv", list(qc[0]), qc)
    print(f"PCA complete: {len(records)} genomes, {seq.shape[1]} nt, {variable_sites} variable sites, {matrix.shape[1]} encoded features")
    print("Explained variance:", ", ".join(f"PC{i+1}={explained[i]*100:.2f}%" for i in range(min(5, ncomp))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
