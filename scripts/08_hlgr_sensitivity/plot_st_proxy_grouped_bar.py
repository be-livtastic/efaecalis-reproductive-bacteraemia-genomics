#!/usr/bin/env python3
"""Plot grouped source-specific HLGR proxy prevalence bars for each ST."""

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
OUTPUT_NAME = "hlgr_st_proxy_grouped_bar.png"


def project_root() -> Path:
    override = os.environ.get("EFAECALIS_PROJECT_ROOT")
    return Path(override).resolve() if override else Path(__file__).resolve().parents[2]


def st_sort_key(value: str) -> tuple[int, int | str]:
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def generate_plot(st: pd.DataFrame, path: Path) -> None:
    required = {"ST", "source", "denominator", "prevalence"}
    missing = required - set(st.columns)
    if missing:
        raise ValueError("ST prevalence table lacks columns: " + ", ".join(sorted(missing)))
    if not set(st["source"]).issubset(SOURCES):
        raise ValueError("ST prevalence table contains an unexpected source")

    order = sorted(st["ST"].astype(str).unique(), key=st_sort_key)
    complete = pd.MultiIndex.from_product([order, SOURCES], names=["ST", "source"])
    plot = (st.assign(ST=st.ST.astype(str)).set_index(["ST", "source"])
            .reindex(complete).reset_index())
    plot["denominator"] = plot["denominator"].fillna(0).astype(int)
    plot["prevalence"] = plot["prevalence"].fillna(0.0)

    x = np.arange(len(order))
    width = 0.39
    fig, ax = plt.subplots(figsize=(max(11.5, 0.62 * len(order)), 6.2))
    for offset, source in zip((-width / 2, width / 2), SOURCES):
        values = plot[plot.source.eq(source)].set_index("ST").loc[order]
        bars = ax.bar(x + offset, values.prevalence, width, color=COLORS[source], label=source)
        for bar, denominator in zip(bars, values.denominator):
            if denominator == 0:
                continue
            y = max(bar.get_height() + 0.018, 0.018)
            ax.text(bar.get_x() + bar.get_width() / 2, y, f"n={denominator}",
                    ha="center", va="bottom", rotation=90, fontsize=7)

    ax.set_xticks(x, order, rotation=45, ha="right")
    ax.set_ylim(0, 1.13)
    ax.set_xlabel("Sequence type (ST)")
    ax.set_ylabel("HLGR proxy prevalence")
    ax.set_title("HLGR proxy prevalence by sequence type and source")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="Replace the deterministic grouped bar chart.")
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
    print(f"Wrote grouped ST prevalence chart to {output_path}")


if __name__ == "__main__":
    main()
