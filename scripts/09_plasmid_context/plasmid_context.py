#!/usr/bin/env python3
"""Map frozen accepted AMRFinderPlus hits to canonical NCBI replicon records."""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "efaecalis-matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch


PROXY_GENE = "aac(6')-Ie/aph(2'')-Ia"
ACCEPTED_STATUS = "included_in_functional_presence"
DISPLAY_CLASS_ORDER = (
    "Aminoglycoside",
    "Macrolide/lincosamide",
    "Tetracycline",
    "Vancomycin/glycopeptide",
    "Other",
)
DISPLAY_CLASS_COLORS = {
    "Aminoglycoside": "#4477AA",
    "Macrolide/lincosamide": "#AA3377",
    "Tetracycline": "#EE7733",
    "Vancomycin/glycopeptide": "#228833",
    "Other": "#999933",
}
TABLE_NAMES = (
    "accepted_amr_hit_replicon_context_72.tsv",
    "amr_gene_by_replicon_matrix_72.tsv",
    "amr_burden_by_gene_class_and_localisation_72.tsv",
    "st6_hlgr_replicon_context_13.tsv",
    "plasmid_context_qc_summary.tsv",
)
FIGURE_NAMES = (
    "amr_hit_localisation_by_source_72.png",
    "plasmid_amr_gene_by_replicon_heatmap_72.png",
    "hlgr_proxy_localisation_by_source_72.png",
    "amr_burden_by_gene_class_and_localisation_72.png",
)


def project_root() -> Path:
    override = os.environ.get("EFAECALIS_PROJECT_ROOT")
    return Path(override).resolve() if override else Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="Replace this module's deterministic outputs.")
    parser.add_argument(
        "--sequence-reports-root",
        type=Path,
        help="Optional directory beneath which the 72 canonical sequence_report.jsonl files are located.",
    )
    return parser.parse_args()


def preflight(paths: list[Path], overwrite: bool) -> None:
    if len(paths) != len(set(paths)):
        raise ValueError("Output registry contains duplicate paths")
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError("Refusing to overwrite existing outputs; rerun with --overwrite: " + ", ".join(existing))
    for directory in {path.parent for path in paths}:
        directory.mkdir(parents=True, exist_ok=True)


def sequence_report_paths(root: Path, override: Path | None) -> list[Path]:
    search_root = override.resolve() if override else root / "local_archive/large_outputs"
    paths = sorted(search_root.glob("**/sequence_report.jsonl"))
    # Limit discovery to the canonical AMRFinderPlus NCBI dataset layout.
    paths = [p for p in paths if "AMRFinder" in str(p) and "ncbi_dataset/data" in p.as_posix()]
    if len(paths) != 72:
        raise ValueError(f"Expected 72 canonical AMRFinderPlus sequence reports below {search_root}; observed {len(paths)}")
    return paths


def load_replicons(paths: list[Path]) -> tuple[pd.DataFrame, dict[tuple[str, str], dict]]:
    records: list[dict] = []
    aliases: dict[tuple[str, str], dict] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                record = json.loads(line)
                record["sequence_report_path"] = str(path)
                records.append(record)
                for field in ("genbankAccession", "refseqAccession", "sequenceName"):
                    alias = record.get(field)
                    if alias:
                        key = (str(record.get("assemblyAccession")), str(alias))
                        if key in aliases and aliases[key] != record:
                            raise ValueError(f"Replicon alias is not unique within assembly: {key}")
                        aliases[key] = record
    assemblies = {record.get("assemblyAccession") for record in records}
    if len(assemblies) != 72:
        raise ValueError(f"Sequence reports contain {len(assemblies)} unique assemblies rather than 72")
    return pd.DataFrame(records), aliases


def localisation(record: dict) -> str:
    value = str(record.get("assignedMoleculeLocationType", "")).strip().lower()
    if value == "chromosome":
        return "Chromosome"
    if value == "plasmid":
        return "Plasmid"
    return "Unplaced/other"


def accession(record: dict, contig: str) -> str:
    for field in ("genbankAccession", "refseqAccession", "sequenceName"):
        if record.get(field) == contig:
            return str(record[field])
    return str(record.get("genbankAccession") or record.get("refseqAccession") or record.get("sequenceName") or contig)


def display_amr_class(amr_class: str) -> str:
    value = str(amr_class).upper()
    if value == "AMINOGLYCOSIDE":
        return "Aminoglycoside"
    if "LINCOSAMIDE" in value or "MACROLIDE" in value or "STREPTOGRAMIN" in value:
        return "Macrolide/lincosamide"
    if value == "TETRACYCLINE":
        return "Tetracycline"
    if value == "GLYCOPEPTIDE":
        return "Vancomycin/glycopeptide"
    return "Other"


def build_context(amr: pd.DataFrame, metadata: pd.DataFrame, mlst: pd.DataFrame, aliases: dict[tuple[str, str], dict], root: Path) -> pd.DataFrame:
    accepted = amr[amr["Functional_matrix_status"].eq(ACCEPTED_STATUS)].copy()
    if accepted.empty:
        raise ValueError("Frozen AMR table has no accepted functional hits")
    missing = sorted(
        f"{row.Genome}:{row.Contig_ID}" for row in accepted.itertuples(index=False)
        if (row.Genome, row.Contig_ID) not in aliases
    )
    if missing:
        raise ValueError("Accepted AMR contigs absent from sequence reports: " + ", ".join(missing))
    meta = metadata[["assembly_accession", "dataset_category", "bioproject_accession"]]
    sts = mlst[["Genome", "ST"]]
    joined = accepted.merge(meta, left_on="Genome", right_on="assembly_accession", validate="many_to_one")
    joined = joined.merge(sts, on="Genome", validate="many_to_one")
    if len(joined) != len(accepted) or not (joined.Source == joined.dataset_category).all():
        raise ValueError("AMR, metadata, and MLST joins are incomplete or source-inconsistent")
    rows = []
    for row in joined.itertuples(index=False):
        record = aliases[(row.Genome, row.Contig_ID)]
        if record.get("assemblyAccession") != row.Genome:
            raise ValueError(f"Contig {row.Contig_ID} resolves to the wrong assembly")
        loc = localisation(record)
        report_path = Path(record["sequence_report_path"])
        rows.append({
            "genome_accession": row.Genome,
            "source": row.Source,
            "ST": row.ST,
            "BioProject": row.bioproject_accession,
            "AMR_gene": row.Gene,
            "AMR_class": row.Class_for_analysis,
            "AMR_plot_class": display_amr_class(row.Class_for_analysis),
            "AMRFinderPlus_contig_ID": row.Contig_ID,
            "replicon_accession": accession(record, row.Contig_ID),
            "replicon_description": record.get("chrName") or record.get("sequenceName") or "",
            "localisation_classification": loc,
            "evidence_used": (
                "Accepted AMRFinderPlus hit (Functional_matrix_status=" + ACCEPTED_STATUS + "); "
                "exact contig-accession match to NCBI Datasets sequence_report.jsonl; "
                "assignedMoleculeLocationType=" + str(record.get("assignedMoleculeLocationType", "not supplied"))
            ),
            "replicon_length_bp": record.get("length", ""),
            "replicon_role": record.get("role", ""),
            "sequence_report_path": report_path.relative_to(root).as_posix(),
            "AMRFinderPlus_raw_record_ID": row.Raw_record_ID,
        })
    return pd.DataFrame(rows).sort_values(
        ["genome_accession", "replicon_accession", "AMR_gene", "AMRFinderPlus_raw_record_ID"]
    ).reset_index(drop=True)


def gene_matrix(context: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "genome_accession", "source", "ST", "BioProject", "replicon_accession",
        "replicon_description", "localisation_classification", "replicon_length_bp",
    ]
    matrix = context.assign(value=1).pivot_table(
        index=keys, columns="AMR_gene", values="value", aggfunc="sum", fill_value=0
    ).reset_index()
    matrix.columns.name = None
    gene_columns = sorted(set(context.AMR_gene))
    return matrix[keys + gene_columns].sort_values(["genome_accession", "replicon_accession"]).reset_index(drop=True)


def burden_by_class(context: pd.DataFrame) -> pd.DataFrame:
    index = pd.MultiIndex.from_product(
        [DISPLAY_CLASS_ORDER, ("Chromosome", "Plasmid")],
        names=["AMR_plot_class", "localisation_classification"],
    )
    burden = (context[context.localisation_classification.isin(["Chromosome", "Plasmid"])]
              .groupby(["AMR_plot_class", "localisation_classification"]).size()
              .reindex(index, fill_value=0).rename("accepted_AMR_hit_count").reset_index())
    burden["AMR_plot_class"] = pd.Categorical(
        burden.AMR_plot_class, categories=DISPLAY_CLASS_ORDER, ordered=True
    )
    return burden.sort_values(["AMR_plot_class", "localisation_classification"]).reset_index(drop=True)


def st6_table(context: pd.DataFrame, mlst: pd.DataFrame) -> pd.DataFrame:
    st6_genomes = set(mlst.loc[mlst.ST.astype(str).eq("6"), "Genome"])
    if len(st6_genomes) != 13:
        raise ValueError(f"Expected 13 ST6 genomes; observed {len(st6_genomes)}")
    proxy = context[context.genome_accession.isin(st6_genomes) & context.AMR_gene.eq(PROXY_GENE)].copy()
    if len(proxy) != 13 or proxy.genome_accession.nunique() != 13:
        raise ValueError("Each of the 13 ST6 genomes must have exactly one accepted HLGR proxy hit")
    cargo = (context.groupby(["genome_accession", "replicon_accession"])["AMR_gene"]
             .agg(lambda values: ";".join(sorted(set(values) - {PROXY_GENE}))).rename("other_AMR_genes_on_replicon"))
    proxy = proxy.join(cargo, on=["genome_accession", "replicon_accession"])
    proxy["plasmid_name"] = np.where(proxy.localisation_classification.eq("Plasmid"), proxy.replicon_description, "")
    proxy["plasmid_length_bp"] = np.where(proxy.localisation_classification.eq("Plasmid"), proxy.replicon_length_bp, "")
    proxy["other_AMR_genes_on_plasmid"] = np.where(
        proxy.localisation_classification.eq("Plasmid"), proxy.other_AMR_genes_on_replicon, ""
    )
    return proxy.rename(columns={
        "replicon_accession": "carrying_replicon",
        "localisation_classification": "chromosome_plasmid_status",
        "replicon_length_bp": "carrying_replicon_length_bp",
    })[[
        "genome_accession", "source", "ST", "BioProject", "AMR_gene", "carrying_replicon",
        "chromosome_plasmid_status", "replicon_description", "carrying_replicon_length_bp",
        "plasmid_name", "plasmid_length_bp", "other_AMR_genes_on_replicon",
        "other_AMR_genes_on_plasmid", "evidence_used",
    ]].sort_values("genome_accession").reset_index(drop=True)


def optional_similar_plasmids(st6: pd.DataFrame, output_dir: Path) -> int:
    plasmids = st6[st6.chromosome_plasmid_status.eq("Plasmid")]
    candidates = []
    for (_, left), (_, right) in combinations(plasmids.iterrows(), 2):
        lname = re.sub(r"[^a-z0-9]", "", str(left.plasmid_name).lower())
        rname = re.sub(r"[^a-z0-9]", "", str(right.plasmid_name).lower())
        name_similar = bool(lname and rname and (lname == rname or lname in rname or rname in lname))
        llen, rlen = float(left.plasmid_length_bp), float(right.plasmid_length_bp)
        length_ratio = min(llen, rlen) / max(llen, rlen)
        lcargo = set(filter(None, str(left.other_AMR_genes_on_plasmid).split(";")))
        rcargo = set(filter(None, str(right.other_AMR_genes_on_plasmid).split(";")))
        union = lcargo | rcargo
        cargo_jaccard = len(lcargo & rcargo) / len(union) if union else 1.0
        if name_similar and length_ratio >= 0.90 and cargo_jaccard >= 0.50:
            candidates.append({
                "genome_1": left.genome_accession, "plasmid_1": left.carrying_replicon,
                "genome_2": right.genome_accession, "plasmid_2": right.carrying_replicon,
                "normalised_name_similar": True, "length_ratio": length_ratio,
                "AMR_cargo_jaccard": cargo_jaccard,
                "interpretation": "screening candidate only; sequence similarity was not assessed",
            })
    if candidates:
        optional_dir = output_dir / "optional_similar_plasmids"
        optional_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(candidates).to_csv(optional_dir / "similar_plasmid_candidates.tsv", sep="\t", index=False)
    return len(candidates)


def localisation_plot(context: pd.DataFrame, path: Path) -> None:
    counts = context.groupby(["source", "localisation_classification"]).size().unstack(fill_value=0)
    counts = counts.reindex(index=["Reproductive", "Bacteraemia"], columns=["Chromosome", "Plasmid", "Unplaced/other"], fill_value=0)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    bottom = np.zeros(len(counts))
    colors = {"Chromosome": "#4477AA", "Plasmid": "#EE6677", "Unplaced/other": "#BBBBBB"}
    for column in counts.columns:
        ax.bar(counts.index, counts[column], bottom=bottom, label=column, color=colors[column])
        bottom += counts[column].to_numpy()
    ax.set_ylabel("Accepted AMRFinderPlus hits")
    ax.set_title("AMR-hit replicon localisation by source")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plasmid_heatmap(context: pd.DataFrame, path: Path) -> None:
    plasmid = context[context.localisation_classification.eq("Plasmid")]
    matrix = plasmid.assign(value=1).pivot_table(
        index=["genome_accession", "replicon_accession"], columns="AMR_gene", values="value", aggfunc="sum", fill_value=0
    )
    gene_classes = (plasmid[["AMR_gene", "AMR_plot_class"]].drop_duplicates()
                    .set_index("AMR_gene")["AMR_plot_class"].to_dict())
    columns = sorted(matrix.columns, key=lambda gene: (DISPLAY_CLASS_ORDER.index(gene_classes[gene]), gene))
    matrix = matrix.reindex(columns, axis=1)
    fig, ax = plt.subplots(figsize=(max(8, 0.42 * len(matrix.columns)), max(6, 0.22 * len(matrix))))
    present = matrix.to_numpy() > 0
    pixels = np.empty((len(matrix), len(matrix.columns), 4), dtype=float)
    pixels[:] = to_rgba("#F1F3F5")
    for column_index, gene in enumerate(matrix.columns):
        pixels[present[:, column_index], column_index, :] = to_rgba(DISPLAY_CLASS_COLORS[gene_classes[gene]])
    ax.imshow(pixels, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(range(len(matrix)), [f"{g} | {r}" for g, r in matrix.index], fontsize=6)
    ax.set_xlabel("AMR gene (presence/absence)")
    ax.set_ylabel("Genome | plasmid replicon")
    ax.set_title("AMR cargo detected on plasmid replicons")
    class_sequence = [gene_classes[gene] for gene in matrix.columns]
    for boundary in range(1, len(class_sequence)):
        if class_sequence[boundary] != class_sequence[boundary - 1]:
            ax.axvline(boundary - 0.5, color="white", linewidth=2.0)
            ax.axvline(boundary - 0.5, color="#555555", linewidth=0.6)
    present_classes = [name for name in DISPLAY_CLASS_ORDER if name in class_sequence]
    handles = [Patch(facecolor=DISPLAY_CLASS_COLORS[name], label=name) for name in present_classes]
    handles.append(Patch(facecolor="#F1F3F5", edgecolor="#CCCCCC", label="Not detected"))
    ax.legend(handles=handles, title="AMR class", frameon=False, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    fig.tight_layout(rect=(0, 0, 0.83, 1))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def burden_plot(burden: pd.DataFrame, path: Path) -> None:
    x = np.arange(len(DISPLAY_CLASS_ORDER))
    width = 0.39
    colors = {"Chromosome": "#4477AA", "Plasmid": "#EE6677"}
    fig, ax = plt.subplots(figsize=(10.2, 5.6))
    for offset, location in zip((-width / 2, width / 2), colors):
        values = (burden[burden.localisation_classification.eq(location)]
                  .set_index("AMR_plot_class").reindex(DISPLAY_CLASS_ORDER))
        bars = ax.bar(x + offset, values.accepted_AMR_hit_count, width,
                      color=colors[location], label=location)
        ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_xticks(x, DISPLAY_CLASS_ORDER, rotation=25, ha="right")
    ax.set_xlabel("AMR class")
    ax.set_ylabel("Count of accepted AMR hits")
    ax.set_title("Chromosome versus plasmid AMR burden by gene class")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def proxy_plot(context: pd.DataFrame, path: Path) -> None:
    proxy = context[context.AMR_gene.eq(PROXY_GENE)]
    counts = proxy.groupby(["source", "localisation_classification"]).size().unstack(fill_value=0)
    counts = counts.reindex(index=["Reproductive", "Bacteraemia"], columns=["Chromosome", "Plasmid", "Unplaced/other"], fill_value=0)
    x = np.arange(len(counts.index)); width = 0.25
    colors = {"Chromosome": "#4477AA", "Plasmid": "#EE6677", "Unplaced/other": "#BBBBBB"}
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for i, column in enumerate(counts.columns):
        bars = ax.bar(x + (i - 1) * width, counts[column], width, label=column, color=colors[column])
        ax.bar_label(bars, padding=2, fontsize=9)
    ax.set_xticks(x, counts.index)
    ax.set_ylabel("Accepted HLGR proxy hits")
    ax.set_title("HLGR proxy replicon localisation by source")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = project_root()
    paths = {
        "amr": root / "data/processed/amr/amr_standardised_long_72.csv",
        "metadata": root / "data/metadata/sample_metadata_source_72.tsv",
        "mlst": root / "data/processed/mlst/formal_mlst_assignments_72.csv",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen canonical input(s): " + ", ".join(missing))
    report_paths = sequence_report_paths(root, args.sequence_reports_root)
    table_dir = root / "results/tables/plasmid_context"
    figure_dir = root / "results/figures/plasmid_context"
    outputs = [table_dir / name for name in TABLE_NAMES] + [figure_dir / name for name in FIGURE_NAMES]
    preflight(outputs, args.overwrite)

    amr = pd.read_csv(paths["amr"], dtype=str, keep_default_na=False)
    metadata = pd.read_csv(paths["metadata"], sep="\t", dtype=str, keep_default_na=False)
    mlst = pd.read_csv(paths["mlst"], dtype=str, keep_default_na=False)
    if len(metadata) != 72 or metadata.assembly_accession.nunique() != 72:
        raise ValueError("Canonical metadata must contain 72 unique genomes")
    if len(mlst) != 72 or mlst.Genome.nunique() != 72:
        raise ValueError("Frozen MLST table must contain 72 unique genomes")
    replicons, aliases = load_replicons(report_paths)
    context = build_context(amr, metadata, mlst, aliases, root)
    matrix = gene_matrix(context)
    burden = burden_by_class(context)
    st6 = st6_table(context, mlst)
    similar_count = optional_similar_plasmids(st6, table_dir)
    qc = pd.DataFrame([
        {"metric": "canonical_genomes", "value": metadata.assembly_accession.nunique()},
        {"metric": "sequence_reports", "value": len(report_paths)},
        {"metric": "sequence_report_replicons", "value": len(replicons)},
        {"metric": "accepted_AMR_hits", "value": len(context)},
        {"metric": "accepted_hits_with_resolved_replicon", "value": context.replicon_accession.notna().sum()},
        {"metric": "chromosome_localised_hits", "value": context.localisation_classification.eq("Chromosome").sum()},
        {"metric": "plasmid_localised_hits", "value": context.localisation_classification.eq("Plasmid").sum()},
        {"metric": "ST6_genomes", "value": st6.genome_accession.nunique()},
        {"metric": "ST6_proxy_hits_on_plasmids", "value": st6.chromosome_plasmid_status.eq("Plasmid").sum()},
        {"metric": "optional_similar_plasmid_candidate_pairs", "value": similar_count},
    ])

    context.to_csv(table_dir / TABLE_NAMES[0], sep="\t", index=False)
    matrix.to_csv(table_dir / TABLE_NAMES[1], sep="\t", index=False)
    burden.to_csv(table_dir / TABLE_NAMES[2], sep="\t", index=False)
    st6.to_csv(table_dir / TABLE_NAMES[3], sep="\t", index=False)
    qc.to_csv(table_dir / TABLE_NAMES[4], sep="\t", index=False)
    localisation_plot(context, figure_dir / FIGURE_NAMES[0])
    plasmid_heatmap(context, figure_dir / FIGURE_NAMES[1])
    proxy_plot(context, figure_dir / FIGURE_NAMES[2])
    burden_plot(burden, figure_dir / FIGURE_NAMES[3])
    print(f"Wrote 5 tables to {table_dir} and 4 figures to {figure_dir}")


if __name__ == "__main__":
    main()
