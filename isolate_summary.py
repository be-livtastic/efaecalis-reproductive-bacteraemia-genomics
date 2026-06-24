import pandas as pd
from pathlib import Path

input_file = Path("06_summary/amrfinder_all_14.tsv")
output_file = Path("06_summary/amrfinder_summary_by_genome.csv")

df = pd.read_csv(input_file, sep="\t")

print("Columns in AMRFinderPlus output:")
print(df.columns.tolist())

name_col = "Name"
gene_col = "Element symbol"
type_col = "Type"
class_col = "Class"

def collapse_unique(series):
    vals = sorted(set(
        str(x) for x in series.dropna()
        if str(x).strip() not in ["", "nan", "None"]
    ))
    return "; ".join(vals)

summary_rows = []

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

summary = pd.DataFrame(summary_rows)
summary.to_csv(output_file, index=False)

print(f"\nSaved summary to: {output_file}")
print(summary)
