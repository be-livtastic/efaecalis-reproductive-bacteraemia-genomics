#!/usr/bin/env python3
"""Record strand-aware neighbouring CDS annotations for every selected locus."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from phylogeny_common import parse_attributes, read_tsv, write_tsv


# --- Report flanking CDS context around every selected locus ---
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--flank-count", type=int, default=2)
    args = parser.parse_args()
    if args.flank_count < 1:
        raise SystemExit("--flank-count must be positive")

    selected = read_tsv(args.coordinates)
    by_gff: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        by_gff[row["gff_path"]].append(row)
    output = []
    for gff_name, selected_rows in by_gff.items():
        features: dict[str, list[dict[str, object]]] = defaultdict(list)
        with Path(gff_name).open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                attrs = parse_attributes(fields[8])
                features[fields[0]].append({"start": int(fields[3]), "end": int(fields[4]),
                    "strand": fields[6], "feature_id": attrs.get("ID", ""),
                    "gene": attrs.get("gene", attrs.get("Name", "")),
                    "product": attrs.get("product", "")})
        for contig in features:
            features[contig].sort(key=lambda x: (x["start"], x["end"]))
        for selected_row in selected_rows:
            contig_features = features[selected_row["contig"]]
            indices = [i for i, feature in enumerate(contig_features)
                       if feature["feature_id"] == selected_row["feature_id"]]
            if len(indices) != 1:
                raise SystemExit(f"Could not uniquely locate {selected_row['feature_id']} in {gff_name}")
            index = indices[0]
            for offset in range(-args.flank_count, args.flank_count + 1):
                neighbour_index = index + offset
                if neighbour_index < 0 or neighbour_index >= len(contig_features):
                    continue
                feature = contig_features[neighbour_index]
                # Relative transcriptional direction is more biologically useful than file order.
                relative = offset if selected_row["strand"] == "+" else -offset
                output.append({"sample_id": selected_row["sample_id"],
                    "canonical_gene": selected_row["canonical_gene"],
                    "selected_feature_id": selected_row["feature_id"], "relative_position": relative,
                    "relationship": "selected" if relative == 0 else ("upstream" if relative < 0 else "downstream"),
                    "contig": selected_row["contig"], "feature_id": feature["feature_id"],
                    "gene_name": feature["gene"], "product": feature["product"],
                    "start": feature["start"], "end": feature["end"], "strand": feature["strand"]})
    fields = ["sample_id", "canonical_gene", "selected_feature_id", "relative_position", "relationship",
              "contig", "feature_id", "gene_name", "product", "start", "end", "strand"]
    write_tsv(args.output, output, fields)
    print(f"Neighbourhood QC: selected_loci={len(selected)} rows={len(output)} flank_count={args.flank_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
