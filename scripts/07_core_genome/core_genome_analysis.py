#!/usr/bin/env python3
"""QC and summarisation helpers for the staged 72-genome core analysis."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO

ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+)")
EXPECTED_GROUPS = {"Reproductive": 14, "Bacteraemia": 58}
EXPECTED_GENOMES = 72
CORE_MIN = 69


def project_root() -> Path:
    override = __import__("os").environ.get("EFAECALIS_PROJECT_ROOT")
    return Path(override).resolve() if override else Path(__file__).resolve().parents[2]


def atomic_csv(frame: pd.DataFrame, path: Path, overwrite: bool, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, index=False, sep=sep)
    temporary.replace(path)


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def extract_accession(text: str) -> str:
    match = ACCESSION_RE.search(text)
    if not match:
        raise ValueError(f"Could not parse a versioned accession from: {text}")
    return match.group(1)


def read_manifest(root: Path) -> pd.DataFrame:
    path = root / "data/accession_lists/selected_72_accessions.tsv"
    frame = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(frame, ["assembly_accession", "dataset_category"], str(path))
    if len(frame) != EXPECTED_GENOMES or frame["assembly_accession"].nunique() != EXPECTED_GENOMES:
        raise ValueError("Canonical manifest must contain exactly 72 unique accessions")
    counts = frame["dataset_category"].value_counts().to_dict()
    if counts != EXPECTED_GROUPS:
        raise ValueError(f"Canonical source counts must be {EXPECTED_GROUPS}; observed {counts}")
    frame = frame.rename(columns={"assembly_accession": "Genome", "dataset_category": "Source"})
    return frame.sort_values("Genome").reset_index(drop=True)


def gff_metrics(path: Path) -> dict[str, object]:
    cds = 0
    fasta = False
    sequence_headers = 0
    genome_length = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            stripped = line.rstrip("\n\r")
            if stripped == "##FASTA":
                fasta = True
                continue
            if not fasta:
                if not line.startswith("#"):
                    fields = stripped.split("\t")
                    if len(fields) >= 3 and fields[2] == "CDS":
                        cds += 1
            elif line.startswith(">"):
                sequence_headers += 1
            else:
                genome_length += sum(base in "ACGTNacgtn" for base in stripped)
    return {
        "CDS_count": cds,
        "Genome_length_bp": genome_length,
        "GFF_size_bytes": path.stat().st_size,
        "Embedded_FASTA": fasta,
        "Embedded_sequence_count": sequence_headers,
    }


def add_iqr_flags(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    flags: list[list[str]] = [[] for _ in range(len(frame))]
    for column in columns:
        q1, q3 = frame[column].quantile([0.25, 0.75])
        iqr = q3 - q1
        low, high = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        for index, value in enumerate(frame[column]):
            if value < low or value > high:
                flags[index].append(f"{column}_outside_1.5_IQR")
    frame["Outlier_flags"] = [";".join(items) if items else "none" for items in flags]
    return frame


def audit_inputs(args: argparse.Namespace) -> None:
    root = project_root()
    manifest = read_manifest(root)
    gff_root = root / "local_archive/large_outputs/Phylogeny_project/phylogeny_work/01_annotations"
    gffs = sorted(gff_root.rglob("*.gff"))
    if len(gffs) != EXPECTED_GENOMES:
        raise ValueError(f"Expected exactly 72 canonical Prokka GFFs; found {len(gffs)}")
    rows: list[dict[str, object]] = []
    for path in gffs:
        accession = extract_accession(path.name)
        row = {"Genome": accession, "Original_GFF": str(path.relative_to(root))}
        row.update(gff_metrics(path))
        rows.append(row)
    qc = pd.DataFrame(rows)
    if qc["Genome"].nunique() != EXPECTED_GENOMES:
        duplicates = qc.loc[qc["Genome"].duplicated(False), "Genome"].tolist()
        raise ValueError(f"Duplicate GFF accessions detected: {duplicates}")
    expected, observed = set(manifest["Genome"]), set(qc["Genome"])
    if expected != observed:
        raise ValueError(f"GFF/manifest mismatch; missing={sorted(expected-observed)}, unexpected={sorted(observed-expected)}")
    if not qc["Embedded_FASTA"].all() or (qc["Embedded_sequence_count"] < 1).any():
        bad = qc.loc[(~qc["Embedded_FASTA"]) | (qc["Embedded_sequence_count"] < 1), "Genome"].tolist()
        raise ValueError(f"GFF files lack embedded sequence required by Panaroo: {bad}")
    qc = manifest.merge(qc, on="Genome", how="left", validate="one_to_one")
    qc = add_iqr_flags(qc, ["CDS_count", "Genome_length_bp", "GFF_size_bytes"])
    stage_dir = root / "data/processed/core_genome/input_gffs_72"
    if stage_dir.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite staged GFF directory: {stage_dir}")
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True)
    for row in qc.itertuples(index=False):
        source = root / row.Original_GFF
        destination = stage_dir / f"{row.Genome}.gff"
        shutil.copy2(source, destination)
        if extract_accession(destination.name) != row.Genome:
            raise AssertionError("Staged filename accession changed unexpectedly")
    qc["Staged_GFF"] = qc["Genome"].map(lambda value: f"data/processed/core_genome/input_gffs_72/{value}.gff")
    atomic_csv(qc, root / "results/tables/core_genome/core_genome_input_qc_72.csv", args.overwrite)
    atomic_csv(qc[["Genome", "Source", "Original_GFF", "Staged_GFF"]], root / "data/processed/core_genome/core_genome_input_manifest_72.tsv", args.overwrite, "\t")
    summary = pd.DataFrame([
        {"Metric": "genomes", "Value": len(qc)},
        {"Metric": "reproductive_genomes", "Value": int((qc.Source == "Reproductive").sum())},
        {"Metric": "bacteraemia_genomes", "Value": int((qc.Source == "Bacteraemia").sum())},
        {"Metric": "gffs_with_embedded_fasta", "Value": int(qc.Embedded_FASTA.sum())},
        {"Metric": "annotation_outliers_flagged", "Value": int((qc.Outlier_flags != "none").sum())},
    ])
    atomic_csv(summary, root / "results/tables/core_genome/core_genome_input_qc_summary_72.csv", args.overwrite)


def panaroo_table(root: Path) -> tuple[pd.DataFrame, list[str]]:
    path = root / "analysis/core_genome/panaroo_strict_core95/gene_presence_absence.csv"
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    manifest = read_manifest(root)
    genomes = manifest["Genome"].tolist()
    missing = [genome for genome in genomes if genome not in frame.columns]
    if missing:
        raise ValueError(f"Panaroo gene table lacks canonical sample columns: {missing}")
    unexpected_accession_columns = [column for column in frame.columns if ACCESSION_RE.fullmatch(column) and column not in genomes]
    if unexpected_accession_columns:
        raise ValueError(f"Unexpected accession columns in Panaroo table: {unexpected_accession_columns}")
    return frame, genomes


def alignment_records(path: Path) -> tuple[list[str], list[str]]:
    alignment = AlignIO.read(path, "fasta")
    names = [record.id for record in alignment]
    sequences = [str(record.seq).upper() for record in alignment]
    if len(set(map(len, sequences))) != 1:
        raise ValueError("Core alignment sequences do not have equal length")
    return names, sequences


def validate_distance_matrix(matrix: pd.DataFrame, expected: set[str]) -> None:
    """Require the canonical square, symmetric distance matrix with a zero diagonal."""
    if matrix.shape != (len(expected), len(expected)) or set(matrix.index) != expected or set(matrix.columns) != expected:
        raise ValueError("SNP-distance matrix must be square with the canonical accession set")
    values = matrix.to_numpy(dtype=float)
    if not np.allclose(values, values.T) or not np.allclose(np.diag(values), 0):
        raise ValueError("SNP-distance matrix must be symmetric with a zero diagonal")


def tied_minimum_neighbours(distances: pd.Series, genome: str) -> tuple[float, list[str]]:
    """Return every accession tied at the minimum non-self distance."""
    candidates = distances.drop(genome)
    minimum = float(candidates.min())
    tied = sorted(candidates.index[np.isclose(candidates.to_numpy(dtype=float), minimum)].tolist())
    return minimum, tied


def panaroo_qc(args: argparse.Namespace) -> None:
    root = project_root()
    frame, genomes = panaroo_table(root)
    presence = frame[genomes].ne("")
    family_counts = presence.sum(axis=1)
    alignment_path = root / "analysis/core_genome/panaroo_strict_core95/core_gene_alignment.aln"
    names, sequences = alignment_records(alignment_path)
    expected = set(genomes)
    if len(names) != EXPECTED_GENOMES or len(set(names)) != EXPECTED_GENOMES or set(names) != expected:
        raise ValueError(f"Core alignment must contain the exact 72 canonical accessions; observed {len(names)} records")
    supplied = []
    for candidate in ["summary_statistics.txt", "summary_statistics.csv", "summary_statistics.tsv"]:
        path = root / "analysis/core_genome/panaroo_strict_core95" / candidate
        if path.exists():
            supplied.append(candidate)
    metrics = [
        ("status", "pass"),
        ("genomes_in_gene_table", len(genomes)),
        ("alignment_taxa", len(names)),
        ("alignment_length", len(sequences[0])),
        ("total_gene_families", len(frame)),
        ("core_gene_families_at_least_69_of_72", int((family_counts >= CORE_MIN).sum())),
        ("accessory_gene_families_below_69_of_72", int((family_counts < CORE_MIN).sum())),
        ("panaroo_supplied_summary_files", ";".join(supplied) if supplied else "none"),
    ]
    atomic_csv(pd.DataFrame(metrics, columns=["Metric", "Value"]), root / "results/tables/core_genome/panaroo_checkpoint_72.csv", args.overwrite)
    per_genome = pd.DataFrame({"Genome": genomes, "Panaroo_family_count": [int(presence[column].sum()) for column in genomes]})
    per_genome = read_manifest(root).merge(per_genome, on="Genome", validate="one_to_one")
    atomic_csv(per_genome, root / "results/tables/core_genome/panaroo_genome_representation_72.csv", args.overwrite)


def alignment_qc(args: argparse.Namespace) -> None:
    root = project_root()
    path = root / "analysis/core_genome/panaroo_strict_core95/core_gene_alignment.aln"
    names, sequences = alignment_records(path)
    expected = set(read_manifest(root).Genome)
    if set(names) != expected or len(names) != EXPECTED_GENOMES:
        raise ValueError("Alignment taxon set does not exactly match the canonical 72 genomes")
    matrix = np.vstack([np.frombuffer(sequence.encode("ascii"), dtype=np.uint8) for sequence in sequences])
    base_codes = np.array([ord(base) for base in "ACGT"], dtype=np.uint8)
    base_counts = np.vstack([(matrix == code).sum(axis=0) for code in base_codes])
    allele_count = (base_counts > 0).sum(axis=0)
    informative_alleles = (base_counts >= 2).sum(axis=0)
    nonmissing = np.isin(matrix, base_codes)
    gaps = matrix == ord("-")
    ambiguous = ~(nonmissing | gaps)
    per_genome = read_manifest(root).set_index("Genome").loc[names].reset_index()
    per_genome["Alignment_length"] = matrix.shape[1]
    per_genome["Non_missing_bases"] = nonmissing.sum(axis=1)
    per_genome["Non_missing_proportion"] = per_genome["Non_missing_bases"] / matrix.shape[1]
    per_genome["Gap_bases"] = gaps.sum(axis=1)
    per_genome["Ambiguous_bases"] = ambiguous.sum(axis=1)
    per_genome["Below_95_percent_non_missing"] = per_genome["Non_missing_proportion"] < 0.95
    summary = pd.DataFrame([
        ("taxa", matrix.shape[0]),
        ("alignment_length", matrix.shape[1]),
        ("variable_sites", int((allele_count >= 2).sum())),
        ("parsimony_informative_sites", int((informative_alleles >= 2).sum())),
        ("invariant_sites_among_observed_bases", int((allele_count == 1).sum())),
        ("all_missing_sites", int((allele_count == 0).sum())),
        ("total_gap_characters", int(gaps.sum())),
        ("total_ambiguous_characters", int(ambiguous.sum())),
        ("minimum_genome_non_missing_proportion", float(per_genome.Non_missing_proportion.min())),
        ("genomes_below_95_percent_non_missing", int(per_genome.Below_95_percent_non_missing.sum())),
    ], columns=["Metric", "Value"])
    atomic_csv(per_genome, root / "results/tables/core_genome/core_alignment_per_genome_qc_72.csv", args.overwrite)
    atomic_csv(summary, root / "results/tables/core_genome/core_alignment_qc_72.csv", args.overwrite)
    if per_genome.Below_95_percent_non_missing.any():
        bad = per_genome.loc[per_genome.Below_95_percent_non_missing, "Genome"].tolist()
        raise ValueError(f"Alignment QC failed: genomes below 95% non-missing: {bad}")


def distance_results(args: argparse.Namespace) -> None:
    root = project_root()
    path = root / "analysis/core_genome/distances/core_alignment_snp_distances_72.tsv"
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    expected = set(read_manifest(root).Genome)
    validate_distance_matrix(matrix, expected)
    metadata = pd.read_csv(root / "results/tables/mlst_amr_phylogeny/integrated_genome_mlst_amr_72.csv", dtype=str)
    require_columns(metadata, ["Genome", "Source", "ST", "HLGR_proxy"], "integrated MLST/AMR table")
    lookup = metadata.set_index("Genome")
    rows = []
    reproductive = sorted(lookup.index[lookup.Source == "Reproductive"])
    for genome in reproductive:
        minimum, tied = tied_minimum_neighbours(matrix.loc[genome], genome)
        for neighbour in tied:
            rows.append({
                "Genome": genome,
                "ST": lookup.at[genome, "ST"],
                "Closest_core_alignment_neighbour": neighbour,
                "Neighbour_source": lookup.at[neighbour, "Source"],
                "Neighbour_ST": lookup.at[neighbour, "ST"],
                "Pairwise_core_alignment_SNP_distance": int(minimum),
                "Same_source": lookup.at[genome, "Source"] == lookup.at[neighbour, "Source"],
                "Same_ST": lookup.at[genome, "ST"] == lookup.at[neighbour, "ST"] and pd.notna(lookup.at[genome, "ST"]),
                "HLGR_proxy_concordance": lookup.at[genome, "HLGR_proxy"] == lookup.at[neighbour, "HLGR_proxy"],
            })
    result = pd.DataFrame(rows)
    atomic_csv(result, root / "results/tables/core_genome/reproductive_core_nearest_neighbours_72.csv", args.overwrite)
    grouped = result.groupby("Genome")
    summary = pd.DataFrame([
        ("reproductive_genomes", len(reproductive)),
        ("nearest_set_includes_reproductive", int(grouped.Neighbour_source.apply(lambda x: (x == "Reproductive").any()).sum())),
        ("nearest_set_includes_bacteraemia", int(grouped.Neighbour_source.apply(lambda x: (x == "Bacteraemia").any()).sum())),
        ("ties_spanning_both_sources", int(grouped.Neighbour_source.apply(lambda x: x.nunique() > 1).sum())),
    ], columns=["Metric", "Value"])
    atomic_csv(summary, root / "results/tables/core_genome/reproductive_core_nearest_neighbour_summary_72.csv", args.overwrite)


def pangenome_summary(args: argparse.Namespace) -> None:
    root = project_root()
    frame, genomes = panaroo_table(root)
    metadata = read_manifest(root).set_index("Genome")
    presence = frame[genomes].ne("")
    counts = presence.sum(axis=1)
    reproductive = [genome for genome in genomes if metadata.at[genome, "Source"] == "Reproductive"]
    bacteraemia = [genome for genome in genomes if metadata.at[genome, "Source"] == "Bacteraemia"]
    rep_counts = presence[reproductive].sum(axis=1)
    bac_counts = presence[bacteraemia].sum(axis=1)
    gene_column = "Gene" if "Gene" in frame.columns else frame.columns[0]
    annotation_column = "Annotation" if "Annotation" in frame.columns else None
    families = pd.DataFrame({
        "Gene_family": frame[gene_column],
        "Annotation": frame[annotation_column] if annotation_column else "unavailable",
        "Reproductive_count": rep_counts,
        "Bacteraemia_count": bac_counts,
        "Total_count": counts,
        "Core_status": np.where(counts >= CORE_MIN, "Core (>=69/72)", "Accessory (<69/72)"),
    })
    families["Dataset_distribution"] = np.select(
        [(families.Reproductive_count > 0) & (families.Bacteraemia_count > 0),
         (families.Reproductive_count > 0) & (families.Bacteraemia_count == 0),
         (families.Reproductive_count == 0) & (families.Bacteraemia_count > 0)],
        ["Shared between datasets", "Detected only in the reproductive-associated dataset (n=14)",
         "Detected only in the bacteraemia-associated dataset (n=58)"], default="Not detected")
    accessory = presence.loc[counts < CORE_MIN]
    per_genome = pd.DataFrame({"Genome": genomes, "Accessory_gene_family_count": accessory.sum(axis=0).astype(int).values})
    per_genome = read_manifest(root).merge(per_genome, on="Genome", validate="one_to_one")
    source_summary = per_genome.groupby("Source").Accessory_gene_family_count.agg(["count", "mean", "median", "min", "max"]).reset_index()
    totals = pd.DataFrame([
        ("total_gene_families", len(frame)),
        ("core_gene_families_at_least_69_of_72", int((counts >= CORE_MIN).sum())),
        ("accessory_gene_families_below_69_of_72", int((counts < CORE_MIN).sum())),
    ], columns=["Metric", "Value"])
    atomic_csv(families.sort_values(["Total_count", "Gene_family"], ascending=[False, True]), root / "results/tables/core_genome/panaroo_gene_family_dataset_distribution_72.csv", args.overwrite)
    atomic_csv(per_genome, root / "results/tables/core_genome/panaroo_accessory_by_genome_72.csv", args.overwrite)
    atomic_csv(source_summary, root / "results/tables/core_genome/panaroo_accessory_by_source_72.csv", args.overwrite)
    atomic_csv(totals, root / "results/tables/core_genome/panaroo_pangenome_summary_72.csv", args.overwrite)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["audit", "panaroo-qc", "alignment-qc", "distance-results", "pangenome-summary"])
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    functions = {
        "audit": audit_inputs,
        "panaroo-qc": panaroo_qc,
        "alignment-qc": alignment_qc,
        "distance-results": distance_results,
        "pangenome-summary": pangenome_summary,
    }
    try:
        functions[args.stage](args)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
