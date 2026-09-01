#!/usr/bin/env python3
"""Downstream ST6-exclusion sensitivity analysis for the frozen HLGR proxy calls."""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "efaecalis-matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from plot_st_proxy_dotplot import generate_plot as save_st_dotplot
from plot_st_proxy_grouped_bar import generate_plot as save_st_grouped_bar


EXPECTED_SOURCE_COUNTS = {"Reproductive": 14, "Bacteraemia": 58}
POSITIVE_LABEL = "HLGR-associated genotype proxy detected"
NEGATIVE_LABEL = "HLGR-associated genotype proxy not detected"
OUTPUT_NAMES = (
    "hlgr_source_prevalence_sensitivity.tsv",
    "hlgr_source_comparison_sensitivity.tsv",
    "hlgr_st_proxy_prevalence.tsv",
)
FIGURE_NAMES = (
    "hlgr_prevalence_before_after_st6_exclusion.png",
    "hlgr_odds_ratio_forest_st6_sensitivity.png",
    "hlgr_st_proxy_dotplot.png",
    "hlgr_st_proxy_grouped_bar.png",
)


def project_root() -> Path:
    override = os.environ.get("EFAECALIS_PROJECT_ROOT")
    return Path(override).resolve() if override else Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="Replace this module's deterministic outputs.")
    return parser.parse_args()


def preflight(paths: list[Path], overwrite: bool) -> None:
    if len(paths) != len(set(paths)):
        raise ValueError("Output registry contains duplicate paths")
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError("Refusing to overwrite existing outputs; rerun with --overwrite: " + ", ".join(existing))
    for directory in {path.parent for path in paths}:
        directory.mkdir(parents=True, exist_ok=True)


def validate_inputs(metadata: pd.DataFrame, proxy: pd.DataFrame) -> pd.DataFrame:
    if len(metadata) != 72 or metadata["assembly_accession"].nunique() != 72:
        raise ValueError("Canonical metadata must contain 72 unique assembly accessions")
    observed = metadata["dataset_category"].value_counts().to_dict()
    if observed != EXPECTED_SOURCE_COUNTS:
        raise ValueError(f"Unexpected canonical source denominators: {observed}")
    if len(proxy) != 72 or proxy["Genome"].nunique() != 72:
        raise ValueError("Canonical proxy table must contain 72 unique genomes")
    if set(proxy["HLGR_proxy"]) != {POSITIVE_LABEL, NEGATIVE_LABEL}:
        raise ValueError("Canonical proxy table contains an unexpected proxy label")
    joined = metadata[["assembly_accession", "dataset_category"]].merge(
        proxy[["Genome", "Source", "ST", "HLGR_proxy"]],
        left_on="assembly_accession",
        right_on="Genome",
        validate="one_to_one",
    )
    if len(joined) != 72 or not (joined["dataset_category"] == joined["Source"]).all():
        raise ValueError("Canonical metadata and proxy table do not agree one-to-one on source")
    joined["proxy_positive"] = joined["HLGR_proxy"].eq(POSITIVE_LABEL)
    joined["ST"] = joined["ST"].astype("string").fillna("Unassigned")
    return joined


def comparison(data: pd.DataFrame, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for source in EXPECTED_SOURCE_COUNTS:
        subset = data[data["Source"].eq(source)]
        positive = int(subset["proxy_positive"].sum())
        denominator = len(subset)
        rows.append(
            {
                "analysis": label,
                "source": source,
                "denominator": denominator,
                "proxy_positive": positive,
                "proxy_negative": denominator - positive,
                "prevalence": positive / denominator,
                "prevalence_percent": 100 * positive / denominator,
            }
        )
    prevalence = pd.DataFrame(rows)
    repro = prevalence.set_index("source").loc["Reproductive"]
    bact = prevalence.set_index("source").loc["Bacteraemia"]
    table = np.array(
        [[repro.proxy_positive, repro.proxy_negative], [bact.proxy_positive, bact.proxy_negative]],
        dtype=float,
    )
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    if np.any(table == 0):
        ci_table = table + 0.5
        ci_method = "Woolf log-odds interval with Haldane-Anscombe 0.5 correction"
    else:
        ci_table = table
        ci_method = "Woolf log-odds interval"
    log_or = math.log((ci_table[0, 0] * ci_table[1, 1]) / (ci_table[0, 1] * ci_table[1, 0]))
    se = math.sqrt(float(np.sum(1 / ci_table)))
    result = pd.DataFrame(
        [{
            "analysis": label,
            "excluded_ST": "None" if label == "Full dataset" else "ST6",
            "total_denominator": len(data),
            "reproductive_denominator": int(repro.denominator),
            "reproductive_proxy_positive": int(repro.proxy_positive),
            "reproductive_proxy_negative": int(repro.proxy_negative),
            "reproductive_prevalence": repro.prevalence,
            "bacteraemia_denominator": int(bact.denominator),
            "bacteraemia_proxy_positive": int(bact.proxy_positive),
            "bacteraemia_proxy_negative": int(bact.proxy_negative),
            "bacteraemia_prevalence": bact.prevalence,
            "prevalence_difference_reproductive_minus_bacteraemia": repro.prevalence - bact.prevalence,
            "odds_ratio_reproductive_vs_bacteraemia": odds_ratio,
            "odds_ratio_95ci_lower": math.exp(log_or - 1.96 * se),
            "odds_ratio_95ci_upper": math.exp(log_or + 1.96 * se),
            "odds_ratio_ci_method": ci_method,
            "fisher_exact_two_sided_p": p_value,
        }]
    )
    return prevalence, result


def save_prevalence_plot(prevalence: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    analyses = ["Full dataset", "ST6 excluded"]
    x = np.arange(len(analyses))
    width = 0.34
    colors = {"Reproductive": "#4477AA", "Bacteraemia": "#CC6677"}
    for offset, source in zip((-width / 2, width / 2), colors):
        values = prevalence[prevalence.source.eq(source)].set_index("analysis").loc[analyses]
        bars = ax.bar(x + offset, values.prevalence, width, label=source, color=colors[source])
        for bar, (_, row) in zip(bars, values.iterrows()):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                    f"{row.proxy_positive:.0f}/{row.denominator:.0f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x, analyses)
    ax.set_ylim(0, max(0.65, prevalence.prevalence.max() + 0.12))
    ax.set_ylabel("HLGR proxy prevalence")
    ax.set_title("HLGR proxy prevalence before and after ST6 exclusion")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_forest_plot(results: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    plot = results.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot))
    estimates = plot["odds_ratio_reproductive_vs_bacteraemia"].to_numpy()
    lower = plot["odds_ratio_95ci_lower"].to_numpy()
    upper = plot["odds_ratio_95ci_upper"].to_numpy()
    ax.errorbar(estimates, y, xerr=[estimates - lower, upper - estimates], fmt="o", color="#332288", capsize=4)
    ax.axvline(1, color="#777777", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_yticks(y, plot["analysis"])
    ax.set_xlabel("Odds ratio (Reproductive vs Bacteraemia), 95% CI")
    ax.set_title("HLGR proxy source comparison")
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = project_root()
    metadata_path = root / "data/metadata/sample_metadata_source_72.tsv"
    proxy_path = root / "results/tables/mlst_amr_phylogeny/integrated_genome_mlst_amr_72.csv"
    for path in (metadata_path, proxy_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen canonical input: {path}")
    table_dir = root / "results/tables/hlgr_sensitivity"
    figure_dir = root / "results/figures/hlgr_sensitivity"
    outputs = [table_dir / name for name in OUTPUT_NAMES] + [figure_dir / name for name in FIGURE_NAMES]
    preflight(outputs, args.overwrite)
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str, keep_default_na=False)
    proxy = pd.read_csv(proxy_path, dtype=str, keep_default_na=False)
    data = validate_inputs(metadata, proxy)
    st6 = data.ST.eq("6")
    if int(st6.sum()) != 13:
        raise ValueError(f"Expected exactly 13 ST6 genomes, observed {int(st6.sum())}")

    full_prevalence, full_result = comparison(data, "Full dataset")
    excluded_prevalence, excluded_result = comparison(data.loc[~st6], "ST6 excluded")
    prevalence = pd.concat([full_prevalence, excluded_prevalence], ignore_index=True)
    results = pd.concat([full_result, excluded_result], ignore_index=True)
    st = (data.groupby(["ST", "Source"], dropna=False)["proxy_positive"]
          .agg(denominator="size", proxy_positive="sum").reset_index())
    st["proxy_positive"] = st.proxy_positive.astype(int)
    st["proxy_negative"] = st.denominator - st.proxy_positive
    st["prevalence"] = st.proxy_positive / st.denominator
    st = st.rename(columns={"Source": "source"}).sort_values(["ST", "source"])

    prevalence.to_csv(table_dir / OUTPUT_NAMES[0], sep="\t", index=False, float_format="%.10g")
    results.to_csv(table_dir / OUTPUT_NAMES[1], sep="\t", index=False, float_format="%.10g")
    st.to_csv(table_dir / OUTPUT_NAMES[2], sep="\t", index=False, float_format="%.10g")
    save_prevalence_plot(prevalence, figure_dir / FIGURE_NAMES[0])
    save_forest_plot(results, figure_dir / FIGURE_NAMES[1])
    save_st_dotplot(st, figure_dir / FIGURE_NAMES[2])
    save_st_grouped_bar(st, figure_dir / FIGURE_NAMES[3])
    print(f"Wrote 3 tables to {table_dir} and 4 figures to {figure_dir}")


if __name__ == "__main__":
    main()
