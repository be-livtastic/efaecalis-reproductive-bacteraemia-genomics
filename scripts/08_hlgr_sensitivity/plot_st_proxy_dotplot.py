#!/usr/bin/env python3
"""Plot source-specific HLGR proxy prevalence as ST-level dots."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "efaecalis-matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SOURCES = ("Reproductive", "Bacteraemia")
COLORS = {"Reproductive": "#4477AA", "Bacteraemia": "#CC6677"}
OUTPUT_NAME = "hlgr_st_proxy_dotplot.png"


def project_root() -> Path:
    override = os.environ.get("EFAECALIS_PROJECT_ROOT")
    return Path(override).resolve() if override else Path(__file__).resolve().parents[2]


def st_sort_key(value: str) -> tuple[int, int | str]:
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def generate_plot(st: pd.DataFrame, path: Path) -> None:
    required = {"ST", "source", "denominator", "proxy_positive", "prevalence"}
    missing = required - set(st.columns)
    if missing:
        raise ValueError("ST prevalence table lacks columns: " + ", ".join(sorted(missing)))
    if not set(st["source"]).issubset(SOURCES):
        raise ValueError("ST prevalence table contains an unexpected source")

    plot = st.copy()
    plot["ST"] = plot["ST"].astype(str)
    order = sorted(plot["ST"].unique(), key=st_sort_key)
    positions = {value: index for index, value in enumerate(order)}
    offsets = {"Reproductive": -0.12, "Bacteraemia": 0.12}

    fig, ax = plt.subplots(figsize=(max(11.5, 0.62 * len(order)), 6.2))
    for source in SOURCES:
        subset = plot[plot.source.eq(source)].copy()
        x = np.array([positions[value] + offsets[source] for value in subset.ST])
        sizes = 55 + 36 * np.sqrt(subset.denominator.astype(float))
        ax.scatter(
            x,
            subset.prevalence,
            s=sizes,
            color=COLORS[source],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.9,
            label=source,
            zorder=3,
        )
        for x_value, row in zip(x, subset.itertuples(index=False)):
            ax.annotate(
                f"{int(row.proxy_positive)}/{int(row.denominator)}",
                (x_value, float(row.prevalence)),
                xytext=(0, 7 if source == "Reproductive" else 18),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
            )

    ax.set_xticks(np.arange(len(order)), order, rotation=45, ha="right")
    ax.set_ylim(-0.04, 1.14)
    ax.set_xlabel("Sequence type (ST)")
    ax.set_ylabel("HLGR proxy prevalence")
    ax.set_title("HLGR proxy prevalence by sequence type and source")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.legend(frameon=False, ncol=2, title="Source")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="Replace the deterministic dot plot.")
    args = parser.parse_args()
    root = project_root()
    input_path = root / "results/tables/hlgr_sensitivity/hlgr_st_proxy_prevalence.tsv"
    output_path = root / "results/figures/hlgr_sensitivity" / OUTPUT_NAME
    if not input_path.is_file():
        raise FileNotFoundError(f"Missing downstream ST prevalence table: {input_path}")
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {output_path}; rerun with --overwrite")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(input_path, sep="\t", dtype={"ST": str})
    generate_plot(st, output_path)
    print(f"Wrote ST proxy dot plot to {output_path}")


if __name__ == "__main__":
    main()
