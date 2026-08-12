"""Summarise one combined AMRFinderPlus TSV without silent overwriting.

The script is working-directory independent. Paths are supplied explicitly so
it behaves consistently on Windows, WSL and Linux.
"""

import argparse
import sys

import pandas as pd
from pathlib import Path

# --- Command-line interface ---
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of an existing output after manual review.",
    )
    return parser.parse_args()


args = parse_args()
input_file = args.input_tsv.resolve()
output_file = args.output_csv.resolve()

# --- Input and overwrite safeguards ---
if not input_file.is_file():
    sys.exit(f"Input TSV does not exist: {input_file}")
if output_file.exists() and not args.overwrite:
    sys.exit(f"Refusing to overwrite existing output: {output_file}")
output_file.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

print("Columns in AMRFinderPlus output:")
print(df.columns.tolist())

name_col = "Name"
gene_col = "Element symbol"
type_col = "Type"
class_col = "Class"

required_columns = {name_col, gene_col, type_col, class_col}
missing_columns = sorted(required_columns.difference(df.columns))
if missing_columns:
    sys.exit("Missing required AMRFinderPlus columns: " + ", ".join(missing_columns))

# --- Reusable unique-value formatter ---
def collapse_unique(series):
    vals = sorted(set(
        str(x) for x in series.dropna()
        if str(x).strip() not in ["", "nan", "None"]
    ))
    return "; ".join(vals)

summary_rows = []

# --- Summarise resistance, virulence and stress hits per genome ---
for genome, sub in df.groupby(name_col):
    amr_sub = sub[sub[type_col].astype(str).str.upper().str.contains("AMR", na=False)]
    vir_sub = sub[sub[type_col].astype(str).str.upper().str.contains("VIRULENCE", na=False)]
    stress_sub = sub[sub[type_col].astype(str).str.upper().str.contains("STRESS", na=False)]

    all_genes = sub[gene_col].dropna().astype(str).tolist()
    amr_genes = amr_sub[gene_col].dropna().astype(str).tolist()

    aminoglycoside_sub = sub[
        sub[gene_col].astype(str).str.contains("aac|aph|ant", case=False, na=False) |
        sub[class_col].astype(str).str.contains("aminoglycoside", case=False, na=False)
    ]

    hlgr_detected = any(
        "aac(6')-Ie-aph(2'')" in gene or
        "aac(6')-Ie-aph(2" in gene
        for gene in all_genes
    )

    summary_rows.append({
        "assembly_accession": genome,
        "total_amrfinder_hits": len(sub),
        "amr_gene_count": len(set(amr_genes)),
        "amr_genes": collapse_unique(amr_sub[gene_col]),
        "amr_classes": collapse_unique(amr_sub[class_col]),
        "aminoglycoside_genes": collapse_unique(aminoglycoside_sub[gene_col]),
        "hlgr_detected": hlgr_detected,
        "virulence_gene_count": len(set(vir_sub[gene_col].dropna().astype(str))),
        "virulence_genes": collapse_unique(vir_sub[gene_col]),
        "stress_gene_count": len(set(stress_sub[gene_col].dropna().astype(str))),
        "stress_genes": collapse_unique(stress_sub[gene_col])
    })

# --- Write the compact per-genome table ---
summary = pd.DataFrame(summary_rows)
summary.to_csv(output_file, index=False)

print(f"\nSaved summary to: {output_file}")
print(summary)
