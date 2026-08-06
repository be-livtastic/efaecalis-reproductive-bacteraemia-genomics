#!/usr/bin/env python3
"""Fail when active configuration and code drift from the 72-genome policy."""
from __future__ import annotations

import ast
from pathlib import Path

from phylogeny_common import read_tsv

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = ["gdh", "gyd", "pstS", "gki", "xpt", "yqiL", "pyrC", "groEL", "recA"]


def python_loci(path: Path) -> list[str] | None:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "LOCI" for t in node.targets):
            return ast.literal_eval(node.value)
    return None


def main() -> int:
    errors = []
    genes = read_tsv(ROOT / "config/phylogeny_9_loci_genes.tsv")
    configured = [r["canonical_gene"] for r in sorted(genes, key=lambda r: int(r["gene_order"]))]
    if configured != EXPECTED: errors.append(f"gene configuration order is {configured}")
    thresholds = read_tsv(ROOT / "config/phylogeny_9_loci_qc_thresholds.tsv")
    if [r["canonical_gene"] for r in thresholds] != EXPECTED: errors.append("QC threshold locus order differs")
    observed = read_tsv(ROOT / "config/phylogeny_9_loci_observed_span_overrides.tsv")
    if len(observed) != 2: errors.append(f"expected 2 observed-span overrides, found {len(observed)}")
    if any(r.get("imputed_bases", "") not in ("", "0") for r in observed):
        errors.append("an observed-span override records imputed bases")
    for name in ("04_validate_sequences.py", "06_concatenate_alignments.py"):
        value = python_loci(ROOT / "scripts/06_phylogenomics" / name)
        if value != EXPECTED: errors.append(f"{name} LOCI differs: {value}")
    runner = (ROOT / "scripts/06_phylogenomics/run_9_locus_phylogeny.sh").read_text(encoding="utf-8")
    for required in ("648", "efaecalis_72_genomes_9_locus_observed_indels", "phylogeny_72_genomes_9_locus_observed_indels", "04b_validate_proteins.py", "04c_validate_pubmlst.py"):
        if required not in runner: errors.append(f"runner lacks required policy token {required}")
    figure = (ROOT / "scripts/06_phylogenomics/08_visualise_9_locus_tree.R").read_text(encoding="utf-8")
    if "Nine-locus" not in figure: errors.append("figure script lacks nine-locus caption")
    if errors:
        raise SystemExit("POLICY CONSISTENCY FAILURE:\n- " + "\n- ".join(errors))
    print("Policy consistency: 9 loci; 72 genomes; 648 records; 2 observed-span overrides; no imputation; active prefixes agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
