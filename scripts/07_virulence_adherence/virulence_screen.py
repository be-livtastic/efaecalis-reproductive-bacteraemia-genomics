#!/usr/bin/env python3
"""Pinned VFDB reference preparation and audited BLAST screening."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

VFDB_URL = "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz"
VFDB_RELEASE = "2026-02-06"
PRIMARY = {
    "ace": "ace", "efaa": "efaA", "ebpa": "ebpA", "ebpb": "ebpB", "ebpc": "ebpC",
    "asa1": "asa1_or_validated_equivalent", "gele": "gelE", "spre": "sprE", "esp": "esp", "cyla": "cylA",
}
TARGET_ORDER = ["ace", "efaA", "ebpA", "ebpB", "ebpC", "asa1_or_validated_equivalent", "gelE", "sprE", "esp", "cylA"]
EXPECTED_GROUPS = {"Reproductive": 14, "Bacteraemia": 58}


def root() -> Path:
    override = os.environ.get("EFAECALIS_PROJECT_ROOT")
    return Path(override).resolve() if override else Path(__file__).resolve().parents[2]


def atomic_frame(frame: pd.DataFrame, path: Path, overwrite: bool, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, index=False, sep=sep)
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_records(handle):
    header = None
    sequence: list[str] = []
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence)
            header, sequence = line[1:], []
        elif line:
            sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def parse_header(header: str) -> dict[str, str]:
    match = re.match(r"^(VFG\d+)\(gb\|([^)]*)\) \(([^)]*)\) (.*?) \[", header)
    if not match:
        raise ValueError(f"Unrecognised VFDB Set A header: {header}")
    vfg, accession, symbol, annotation = match.groups()
    return {"VFG_ID": vfg, "Reference_accession": accession, "Database_gene_symbol": symbol,
            "Database_annotation": annotation, "Original_header": header}


def reference_audit(args: argparse.Namespace) -> None:
    project = root()
    reference_dir = project / f"references/virulence/vfdb_setA_{VFDB_RELEASE}"
    reference_dir.mkdir(parents=True, exist_ok=True)
    archive = reference_dir / "VFDB_setA_nt.fas.gz"
    if not archive.exists() or args.overwrite:
        temporary = archive.with_suffix(archive.suffix + ".tmp")
        request = urllib.request.Request(VFDB_URL, headers={"User-Agent": "Mozilla/5.0 reproducible-academic-download"})
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
            shutil_copyfileobj = __import__("shutil").copyfileobj
            shutil_copyfileobj(response, output)
        temporary.replace(archive)
    selected: list[tuple[dict[str, str], str]] = []
    with gzip.open(archive, "rt", encoding="latin-1") as handle:
        for header, sequence in fasta_records(handle):
            if "Enterococcus faecalis" not in header:
                continue
            record = parse_header(header)
            symbol_key = record["Database_gene_symbol"].lower()
            annotation_lower = record["Database_annotation"].lower()
            header_lower = header.lower()
            if symbol_key in PRIMARY:
                record["Project_target"] = PRIMARY[symbol_key]
                record["Mapping_status"] = "approved_primary"
                record["Mapping_reason"] = "Explicit VFDB symbol matches the predefined project target"
            elif "aggregation substance" in annotation_lower or "[as (vf0352)" in header_lower:
                record["Project_target"] = "aggregation_substance_context"
                record["Mapping_status"] = "review_required_context"
                record["Mapping_reason"] = "Aggregation-substance family member is retained but not equated with asa1"
            elif "cytolysin" in header_lower:
                record["Project_target"] = "cytolysin_operon_context"
                record["Mapping_status"] = "context_only"
                record["Mapping_reason"] = "Recognised cytolysin component retained without inferring operon function or phenotype"
            else:
                continue
            record["Reference_length"] = str(len(sequence))
            record["Query_ID"] = f"{record['VFG_ID']}|{record['Reference_accession']}|{record['Database_gene_symbol']}|{record['Project_target']}"
            selected.append((record, sequence))
    mapping = pd.DataFrame([record for record, _ in selected])
    primary_counts = mapping.loc[mapping.Mapping_status == "approved_primary", "Project_target"].value_counts().to_dict()
    missing = [target for target in TARGET_ORDER if primary_counts.get(target, 0) < 1]
    if missing:
        raise ValueError(f"VFDB reference audit failed; required project targets absent: {missing}")
    fasta_path = reference_dir / "vfdb_enterococcus_target_and_context.fna"
    if fasta_path.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite extracted reference FASTA: {fasta_path}")
    with fasta_path.open("w") as handle:
        for record, sequence in selected:
            handle.write(f">{record['Query_ID']}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start:start + 70] + "\n")
    atomic_frame(mapping, reference_dir / "virulence_reference_mapping.tsv", args.overwrite, "\t")
    checksums = pd.DataFrame([
        {"File": archive.name, "SHA256": sha256(archive)},
        {"File": fasta_path.name, "SHA256": sha256(fasta_path)},
    ])
    atomic_frame(checksums, reference_dir / "SHA256SUMS.tsv", args.overwrite, "\t")
    provenance = pd.DataFrame([{
        "Resource": "VFDB core dataset Set A nucleotide sequences",
        "VFDB_release_date": VFDB_RELEASE,
        "Retrieval_date": __import__("datetime").date.today().isoformat(),
        "Source_URL": VFDB_URL,
        "Raw_record_count_selected": len(mapping),
        "Approved_primary_reference_count": int((mapping.Mapping_status == "approved_primary").sum()),
        "Primary_project_target_count": mapping.loc[mapping.Mapping_status == "approved_primary", "Project_target"].nunique(),
        "Licence": "CC BY-NC 4.0 for academic/non-commercial use; see VFDB terms",
    }])
    atomic_frame(provenance, reference_dir / "provenance.tsv", args.overwrite, "\t")
    qc = mapping.groupby(["Project_target", "Mapping_status"], dropna=False).size().reset_index(name="Reference_count")
    atomic_frame(qc, project / "results/tables/virulence_adherence/virulence_reference_audit.csv", args.overwrite)


def read_metadata(project: Path) -> pd.DataFrame:
    path = project / "results/tables/mlst_amr_phylogeny/integrated_genome_mlst_amr_72.csv"
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"Genome", "Source", "HLGR_proxy"}
    if not required.issubset(frame.columns) or len(frame) != 72 or frame.Genome.nunique() != 72:
        raise ValueError("Frozen integrated MLST/AMR table must contain 72 unique mapped genomes")
    if frame.Source.value_counts().to_dict() != EXPECTED_GROUPS:
        raise ValueError("Frozen integrated source counts are not 14 reproductive and 58 bacteraemia")
    return frame


def fna_manifest(project: Path) -> dict[str, Path]:
    roots = [project / "local_archive/large_outputs/Efaecalis_14_AMRFinder/02_genomes_fna",
             project / "local_archive/large_outputs/Efaecalis_58_AMRFinder/02_genomes_fna"]
    paths = [path for source in roots for path in source.glob("*.fna")]
    mapping = {path.stem: path for path in paths}
    expected = set(read_metadata(project).Genome)
    if len(paths) != 72 or len(mapping) != 72 or set(mapping) != expected:
        raise ValueError("Canonical FNA inputs do not exactly match the 72 frozen integrated genomes")
    return mapping


def run_one_blast(genome: str, fna: Path, query: Path, output: Path, overwrite: bool) -> tuple[str, float]:
    if output.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite BLAST output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    command = ["blastn", "-task", "blastn", "-dust", "no", "-soft_masking", "false", "-evalue", "1e-20",
               "-query", str(query), "-subject", str(fna), "-outfmt",
               "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore", "-out", str(output)]
    started = time.monotonic()
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"blastn failed for {genome}: {result.stderr.strip()}")
    return genome, time.monotonic() - started


def union_length(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, stop in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, stop])
        else:
            merged[-1][1] = max(merged[-1][1], stop)
    return sum(stop - start + 1 for start, stop in merged)


def split_coordinate_loci(group: pd.DataFrame, query_length: int) -> list[pd.DataFrame]:
    """Split HSPs on one contig/strand when their subject spans are not locally adjacent."""
    ordered = group.assign(
        Subject_low=group[["sstart", "send"]].min(axis=1),
        Subject_high=group[["sstart", "send"]].max(axis=1),
    ).sort_values(["Subject_low", "Subject_high"])
    maximum_gap = max(5000, query_length)
    clusters: list[list[int]] = []
    current: list[int] = []
    current_high = -1
    for index, row in ordered.iterrows():
        if current and int(row.Subject_low) > current_high + maximum_gap:
            clusters.append(current)
            current = []
            current_high = -1
        current.append(index)
        current_high = max(current_high, int(row.Subject_high))
    if current:
        clusters.append(current)
    return [group.loc[indexes].copy() for indexes in clusters]


def assign_distinct_target_loci(candidates: pd.DataFrame) -> pd.DataFrame:
    """Collapse overlapping reference matches into genomic loci without merging separated copies."""
    if candidates.empty:
        candidates = candidates.copy()
        candidates["Distinct_target_locus_number_on_contig"] = pd.Series(dtype="Int64")
        candidates["Distinct_target_locus_ID"] = pd.Series(dtype=str)
        return candidates
    annotated = candidates.copy()
    annotated["Distinct_target_locus_number_on_contig"] = 0
    annotated["Distinct_target_locus_ID"] = ""
    keys = ["Genome", "Target_gene", "Contig_ID", "Strand"]
    for (genome, target, contig, strand), group in annotated.groupby(keys, sort=True):
        ordered = group.sort_values(["Start", "Stop", "Reference_accession"])
        locus_number = 0
        current_stop = -1
        for index, row in ordered.iterrows():
            if int(row.Start) > current_stop:
                locus_number += 1
                current_stop = int(row.Stop)
            else:
                current_stop = max(current_stop, int(row.Stop))
            annotated.at[index, "Distinct_target_locus_number_on_contig"] = locus_number
            annotated.at[index, "Distinct_target_locus_ID"] = f"{contig}|{strand}|{locus_number}"
    annotated["Distinct_target_locus_number_on_contig"] = annotated["Distinct_target_locus_number_on_contig"].astype(int)
    return annotated


def classify_candidates(raw: pd.DataFrame, genomes: list[str], mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    query_map = mapping.set_index("Query_ID").to_dict("index")
    consolidated_rows = []
    if not raw.empty:
        raw["Strand"] = raw.apply(lambda row: "+" if row.sstart <= row.send else "-", axis=1)
        for keys, contig_group in raw.groupby(["Genome", "qseqid", "sseqid", "Strand"], sort=True):
            genome, query_id, contig, strand = keys
            info = query_map[query_id]
            query_length = int(contig_group.qlen.iloc[0])
            for locus_number, group in enumerate(split_coordinate_loci(contig_group, query_length), start=1):
                covered = union_length(list(zip(group.qstart.astype(int), group.qend.astype(int))))
                weighted_identity = float((group.pident * group.length).sum() / group.length.sum())
                consolidated_rows.append({
                    "Genome": genome, "Query_ID": query_id, "Target_gene": info["Project_target"],
                    "Database_gene_symbol": info["Database_gene_symbol"], "Reference_accession": info["Reference_accession"],
                    "Reference_name": info["Database_annotation"], "Mapping_status": info["Mapping_status"],
                    "Contig_ID": contig, "Locus_number_on_contig": locus_number,
                    "Start": int(group[["sstart", "send"]].min().min()),
                    "Stop": int(group[["sstart", "send"]].max().max()), "Strand": strand,
                    "Identity": weighted_identity, "Coverage": 100.0 * covered / query_length,
                    "HSP_count": len(group), "Best_bitscore": float(group.bitscore.max()),
                })
    candidates = assign_distinct_target_loci(pd.DataFrame(consolidated_rows))
    status_rows = []
    for genome in genomes:
        for target in TARGET_ORDER:
            subset = candidates[(candidates.Genome == genome) & (candidates.Target_gene == target)] if not candidates.empty else candidates
            accepted = subset[(subset.Identity >= 80) & (subset.Coverage >= 80)] if not subset.empty else subset
            review = subset[(subset.Identity >= 70) & (subset.Identity < 80) & (subset.Coverage >= 80)] if not subset.empty else subset
            partial = subset[(subset.Identity >= 80) & (subset.Coverage >= 50) & (subset.Coverage < 80)] if not subset.empty else subset
            accepted_locus_count = accepted.Distinct_target_locus_ID.nunique() if not accepted.empty else 0
            if accepted_locus_count > 1:
                status, reason, best = "ambiguous_multiple_hit", "Multiple distinct accepted loci require review", accepted.sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False).iloc[0]
            elif accepted_locus_count == 1:
                status, reason, best = "accepted_present", "One distinct genomic locus meets >=80% identity and >=80% reference coverage", accepted.sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False).iloc[0]
            elif len(review):
                status, reason, best = "review_required", "Coverage is sufficient but identity is 70-<80%", review.sort_values(["Best_bitscore", "Coverage"], ascending=False).iloc[0]
            elif len(partial):
                status, reason, best = "flagged_partial", "Identity is >=80% but reference coverage is 50-<80%", partial.sort_values(["Coverage", "Best_bitscore"], ascending=False).iloc[0]
            else:
                status, reason, best = "not_detected", "No candidate met a retained detection-status threshold", None
            status_rows.append({
                "Genome": genome, "Target_gene": target, "Detection_status": status, "Status_reason": reason,
                "Database_gene_symbol": "" if best is None else best.Database_gene_symbol,
                "Reference_accession": "" if best is None else best.Reference_accession,
                "Contig_ID": "" if best is None else best.Contig_ID,
                "Start": "" if best is None else best.Start, "Stop": "" if best is None else best.Stop,
                "Strand": "" if best is None else best.Strand,
                "Identity": "" if best is None else best.Identity, "Coverage": "" if best is None else best.Coverage,
                "Candidate_locus_count": subset.Distinct_target_locus_ID.nunique() if not subset.empty else 0,
            })
    return candidates, pd.DataFrame(status_rows)


def ambiguity_locus_table(candidates: pd.DataFrame, statuses: pd.DataFrame, metadata: pd.DataFrame,
                          target: str = "asa1_or_validated_equivalent") -> pd.DataFrame:
    """Return one fully annotated row per accepted locus in ambiguous genomes for a target."""
    columns = ["Genome", "Source", "ST", "Study_ID", "Contig_ID", "Locus_number", "Start", "Stop",
               "Strand", "Reference_ID", "Identity", "Coverage", "HSP_count", "Detection_status"]
    ambiguous = statuses.loc[
        (statuses.Target_gene == target) & (statuses.Detection_status == "ambiguous_multiple_hit"), "Genome"
    ].tolist()
    accepted = candidates[
        candidates.Genome.isin(ambiguous) & candidates.Target_gene.eq(target) &
        candidates.Identity.ge(80) & candidates.Coverage.ge(80)
    ].copy()
    if accepted.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    metadata_lookup = metadata.set_index("Genome")
    for genome, genome_rows in accepted.groupby("Genome", sort=True):
        locus_groups = list(genome_rows.groupby("Distinct_target_locus_ID", sort=True))
        locus_groups.sort(key=lambda item: (item[1].Contig_ID.iloc[0], int(item[1].Start.min())))
        for locus_number, (_, locus) in enumerate(locus_groups, start=1):
            best = locus.sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False).iloc[0]
            rows.append({
                "Genome": genome, "Source": metadata_lookup.at[genome, "Source"],
                "ST": metadata_lookup.at[genome, "ST"] if "ST" in metadata_lookup.columns else "",
                "Study_ID": metadata_lookup.at[genome, "Study_ID"] if "Study_ID" in metadata_lookup.columns else "",
                "Contig_ID": best.Contig_ID, "Locus_number": locus_number,
                "Start": int(locus.Start.min()), "Stop": int(locus.Stop.max()), "Strand": best.Strand,
                "Reference_ID": ";".join(sorted(locus.Reference_accession.astype(str).unique())),
                "Identity": best.Identity, "Coverage": best.Coverage,
                "HSP_count": int(locus.HSP_count.sum()), "Detection_status": "ambiguous_multiple_hit",
            })
    return pd.DataFrame(rows, columns=columns)


def assembly_replicon_map(fna: Path) -> dict[str, tuple[str, str]]:
    """Read contig descriptions and classify only explicitly labelled chromosome/plasmid replicons."""
    result = {}
    with fna.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            contig, _, description = header.partition(" ")
            lower = description.lower()
            kind = "plasmid" if "plasmid" in lower else "chromosome" if "chromosome" in lower else "other/unlabelled"
            result[contig] = (kind, description)
    return result


def aggregation_identity_review(project: Path, candidates: pd.DataFrame, statuses: pd.DataFrame,
                                metadata: pd.DataFrame, fnas: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Competitively compare accepted asa1-reference loci with retained aggregation-family references."""
    asa = candidates[
        candidates.Target_gene.eq("asa1_or_validated_equivalent") &
        candidates.Identity.ge(80) & candidates.Coverage.ge(80)
    ].copy()
    context = candidates[candidates.Target_gene.eq("aggregation_substance_context")].copy()
    status_lookup = statuses.set_index(["Genome", "Target_gene"])
    metadata_lookup = metadata.set_index("Genome")
    replicons = {genome: assembly_replicon_map(path) for genome, path in fnas.items()}
    rows = []
    for (genome, locus_id), locus in asa.groupby(["Genome", "Distinct_target_locus_ID"], sort=True):
        best_asa = locus.sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False).iloc[0]
        overlapping = context[
            context.Genome.eq(genome) & context.Contig_ID.eq(best_asa.Contig_ID) &
            context.Strand.eq(best_asa.Strand) & context.Stop.ge(best_asa.Start) & context.Start.le(best_asa.Stop)
        ].sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False)
        best_context = None if overlapping.empty else overlapping.iloc[0]
        context_score = float("nan") if best_context is None else float(best_context.Best_bitscore)
        asa_is_best = best_context is None or float(best_asa.Best_bitscore) > context_score
        replicon_type, description = replicons[genome].get(best_asa.Contig_ID, ("other/unlabelled", "unavailable"))
        rows.append({
            "Genome": genome, "Source": metadata_lookup.at[genome, "Source"],
            "ST": metadata_lookup.at[genome, "ST"] if "ST" in metadata_lookup.columns else "",
            "Study_ID": metadata_lookup.at[genome, "Study_ID"] if "Study_ID" in metadata_lookup.columns else "",
            "Original_genome_status": status_lookup.at[(genome, "asa1_or_validated_equivalent"), "Detection_status"],
            "Contig_ID": best_asa.Contig_ID, "Replicon_type": replicon_type,
            "Replicon_description": description, "Start": int(locus.Start.min()), "Stop": int(locus.Stop.max()),
            "Strand": best_asa.Strand, "Asa1_reference_ID": best_asa.Reference_accession,
            "Asa1_identity": best_asa.Identity, "Asa1_coverage": best_asa.Coverage,
            "Asa1_bitscore": best_asa.Best_bitscore, "Asa1_HSP_count": int(locus.HSP_count.sum()),
            "Best_context_symbol": "none" if best_context is None else best_context.Database_gene_symbol,
            "Best_context_reference_ID": "none" if best_context is None else best_context.Reference_accession,
            "Best_context_identity": "" if best_context is None else best_context.Identity,
            "Best_context_coverage": "" if best_context is None else best_context.Coverage,
            "Best_context_bitscore": "" if best_context is None else best_context.Best_bitscore,
            "Competitive_locus_interpretation": "asa1_reference_best_supported" if asa_is_best else "related_aggregation_substance_reference_best_supported",
        })
    locus_review = pd.DataFrame(rows)
    summaries = []
    for genome, group in locus_review.groupby("Genome", sort=True):
        asa_supported = int(group.Competitive_locus_interpretation.eq("asa1_reference_best_supported").sum())
        context_supported = int(group.Competitive_locus_interpretation.eq("related_aggregation_substance_reference_best_supported").sum())
        summaries.append({
            "Genome": genome, "Source": metadata_lookup.at[genome, "Source"],
            "ST": metadata_lookup.at[genome, "ST"] if "ST" in metadata_lookup.columns else "",
            "Study_ID": metadata_lookup.at[genome, "Study_ID"] if "Study_ID" in metadata_lookup.columns else "",
            "Original_detection_status": group.Original_genome_status.iloc[0],
            "Accepted_asa1_query_locus_count": len(group),
            "Asa1_reference_best_supported_locus_count": asa_supported,
            "Related_family_best_supported_locus_count": context_supported,
            "Review_recommendation": "accepted_presence_with_copy_number_and_context_annotation" if asa_supported else "review_required",
            "Recommendation_reason": "At least one locus is competitively best matched by the explicit asa1 reference" if asa_supported else "Every accepted asa1-query locus is better matched by a related aggregation-substance reference",
        })
    return locus_review, pd.DataFrame(summaries)


def run_aggregation_review(args: argparse.Namespace) -> None:
    """Generate the post-B3 review without changing any B3 status or downstream binary call."""
    project = root()
    candidates = pd.read_csv(project / "data/processed/virulence_adherence/virulence_consolidated_candidates_72.csv")
    statuses = pd.read_csv(project / "data/processed/virulence_adherence/virulence_target_status_by_genome_72.csv")
    metadata = read_metadata(project)
    locus_review, genome_review = aggregation_identity_review(project, candidates, statuses, metadata, fna_manifest(project))
    output = project / "results/tables/virulence_adherence"
    atomic_frame(locus_review, output / "aggregation_substance_locus_identity_review_72.csv", args.overwrite)
    atomic_frame(genome_review, output / "aggregation_substance_genome_review_recommendations_72.csv", args.overwrite)


def best_candidate_fields(subset: pd.DataFrame) -> dict[str, object]:
    """Return display fields for the best retained candidate without discarding the locus audit."""
    if subset.empty:
        return {"Database_gene_symbol": "", "Reference_accession": "", "Contig_ID": "", "Start": "",
                "Stop": "", "Strand": "", "Identity": "", "Coverage": ""}
    best = subset.sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False).iloc[0]
    return {column: best[column] for column in ["Database_gene_symbol", "Reference_accession", "Contig_ID",
                                                "Start", "Stop", "Strand", "Identity", "Coverage"]}


def resolve_for_analysis(args: argparse.Namespace) -> None:
    """Create an approved resolved layer while preserving all original B3 evidence unchanged."""
    project = root()
    candidates = pd.read_csv(project / "data/processed/virulence_adherence/virulence_consolidated_candidates_72.csv")
    original = pd.read_csv(project / "data/processed/virulence_adherence/virulence_target_status_by_genome_72.csv")
    recommendations = pd.read_csv(project / "results/tables/virulence_adherence/aggregation_substance_genome_review_recommendations_72.csv")
    metadata = read_metadata(project)
    genomes = sorted(metadata.Genome)
    recommendation_lookup = recommendations.set_index("Genome")
    rows = []

    # Retain non-aggregation calls, resolving only the approved two-copy cylA case.
    for row in original[original.Target_gene.ne("asa1_or_validated_equivalent")].itertuples(index=False):
        record = row._asdict()
        record["Original_detection_status"] = row.Detection_status
        subset = candidates[(candidates.Genome == row.Genome) & (candidates.Target_gene == row.Target_gene)]
        accepted = subset[subset.Identity.ge(80) & subset.Coverage.ge(80)]
        accepted_count = accepted.Distinct_target_locus_ID.nunique() if not accepted.empty else 0
        if row.Genome == "GCA_029011395.1" and row.Target_gene == "cylA":
            if accepted_count != 2:
                raise ValueError("Approved GCA_029011395.1 cylA resolution requires exactly two accepted loci")
            record["Detection_status"] = "accepted_present"
            record["Status_reason"] = "Reviewed two-copy cylA detection: one chromosomal and one plasmid locus"
        elif row.Detection_status == "ambiguous_multiple_hit":
            raise ValueError(f"Unresolved non-aggregation multiple-locus call remains: {row.Genome} {row.Target_gene}")
        record["Accepted_locus_count"] = accepted_count if record["Detection_status"] == "accepted_present" else 0
        record["Resolution_policy"] = "approved_cylA_copy_number_2" if row.Genome == "GCA_029011395.1" and row.Target_gene == "cylA" else "unchanged_B3_status"
        rows.append(record)

    # Resolve the narrow asa1-specific call using the approved competitive best-reference decision.
    for genome in genomes:
        source = metadata.set_index("Genome").at[genome, "Source"]
        original_row = original[(original.Genome == genome) & (original.Target_gene == "asa1_or_validated_equivalent")].iloc[0]
        subset = candidates[(candidates.Genome == genome) & (candidates.Target_gene == "asa1_or_validated_equivalent")]
        display = best_candidate_fields(subset)
        if genome in recommendation_lookup.index:
            recommendation = recommendation_lookup.loc[genome]
            supported = int(recommendation.Asa1_reference_best_supported_locus_count)
            if recommendation.Review_recommendation == "accepted_presence_with_copy_number_and_context_annotation":
                status = "accepted_present"
                reason = "Explicit asa1 reference is the competitive best match at at least one accepted locus"
            else:
                status = "review_required"
                reason = "Accepted-threshold asa1 query match is better supported by a related aggregation-substance reference"
        else:
            supported = 0
            status = original_row.Detection_status
            reason = original_row.Status_reason
        rows.append({
            "Genome": genome, "Source": source, "Target_gene": "asa1_specific", "Detection_status": status,
            "Status_reason": reason, **display, "Candidate_locus_count": subset.Distinct_target_locus_ID.nunique() if not subset.empty else 0,
            "Original_detection_status": original_row.Detection_status, "Accepted_locus_count": supported,
            "Resolution_policy": "approved_competitive_best_reference_rule",
        })

    # Build the broader family call after collapsing all reference matches at shared genomic coordinates.
    family = candidates[candidates.Target_gene.isin(["asa1_or_validated_equivalent", "aggregation_substance_context"])].copy()
    family["Target_gene"] = "aggregation_substance_family_detected"
    family = assign_distinct_target_loci(family)
    family_locus_rows = []
    for genome in genomes:
        subset = family[family.Genome == genome]
        accepted = subset[subset.Identity.ge(80) & subset.Coverage.ge(80)]
        review = subset[subset.Identity.ge(70) & subset.Identity.lt(80) & subset.Coverage.ge(80)]
        partial = subset[subset.Identity.ge(80) & subset.Coverage.ge(50) & subset.Coverage.lt(80)]
        accepted_count = accepted.Distinct_target_locus_ID.nunique() if not accepted.empty else 0
        if accepted_count:
            status, reason, display_subset = "accepted_present", "At least one accepted aggregation-substance-family locus detected", accepted
        elif not review.empty:
            status, reason, display_subset = "review_required", "Only 70-<80% identity full-coverage family candidates detected", review
        elif not partial.empty:
            status, reason, display_subset = "flagged_partial", "Only partial aggregation-substance-family candidates detected", partial
        else:
            status, reason, display_subset = "not_detected", "No retained aggregation-substance-family candidate", subset.iloc[0:0]
        rows.append({
            "Genome": genome, "Source": metadata.set_index("Genome").at[genome, "Source"],
            "Target_gene": "aggregation_substance_family_detected", "Detection_status": status,
            "Status_reason": reason, **best_candidate_fields(display_subset),
            "Candidate_locus_count": subset.Distinct_target_locus_ID.nunique() if not subset.empty else 0,
            "Original_detection_status": original[(original.Genome == genome) & (original.Target_gene == "asa1_or_validated_equivalent")].Detection_status.iloc[0],
            "Accepted_locus_count": accepted_count if status == "accepted_present" else 0,
            "Resolution_policy": "approved_broad_aggregation_substance_family_rule",
        })
        for locus_number, (_, locus) in enumerate(accepted.groupby("Distinct_target_locus_ID", sort=True), start=1):
            best = locus.sort_values(["Best_bitscore", "Coverage", "Identity"], ascending=False).iloc[0]
            family_locus_rows.append({
                "Genome": genome, "Source": metadata.set_index("Genome").at[genome, "Source"],
                "ST": metadata.set_index("Genome").at[genome, "ST"], "Study_ID": metadata.set_index("Genome").at[genome, "Study_ID"],
                "Contig_ID": best.Contig_ID, "Locus_number": locus_number, "Start": int(locus.Start.min()),
                "Stop": int(locus.Stop.max()), "Strand": best.Strand,
                "Best_reference_symbol": best.Database_gene_symbol, "Best_reference_ID": best.Reference_accession,
                "Identity": best.Identity, "Coverage": best.Coverage, "HSP_count": int(locus.HSP_count.sum()),
                "Matching_reference_IDs": ";".join(sorted(locus.Reference_accession.astype(str).unique())),
                "Detection_status": "accepted_present",
            })

    resolved = pd.DataFrame(rows)
    expected_targets = [target for target in TARGET_ORDER if target != "asa1_or_validated_equivalent"] + ["asa1_specific", "aggregation_substance_family_detected"]
    if len(resolved) != 72 * len(expected_targets) or resolved.Genome.nunique() != 72 or resolved.Target_gene.nunique() != len(expected_targets):
        raise ValueError("Resolved layer must contain exactly 72 genomes across 11 targets")
    if resolved.groupby(["Genome", "Target_gene"]).size().max() != 1:
        raise ValueError("Resolved layer contains duplicate genome-target rows")
    asa_count = resolved[(resolved.Target_gene == "asa1_specific") & (resolved.Detection_status == "accepted_present")].Genome.nunique()
    if asa_count != 13:
        raise ValueError(f"Approved narrow asa1-specific count is 13; resolved {asa_count}")
    cyl = resolved[(resolved.Genome == "GCA_029011395.1") & (resolved.Target_gene == "cylA")].iloc[0]
    if cyl.Detection_status != "accepted_present" or int(cyl.Accepted_locus_count) != 2:
        raise ValueError("Approved cylA two-copy resolution was not retained")
    output = project / "data/processed/virulence_adherence"
    atomic_frame(resolved.sort_values(["Genome", "Target_gene"]), output / "virulence_resolved_target_status_by_genome_72.csv", args.overwrite)
    atomic_frame(pd.DataFrame(family_locus_rows), output / "aggregation_substance_family_accepted_loci_72.csv", args.overwrite)
    qc = resolved.groupby(["Target_gene", "Detection_status"]).size().reset_index(name="Genome_count")
    atomic_frame(qc, project / "results/tables/virulence_adherence/virulence_resolution_qc_72.csv", args.overwrite)


def screen(args: argparse.Namespace) -> None:
    project = root()
    mapping_path = project / f"references/virulence/vfdb_setA_{VFDB_RELEASE}/virulence_reference_mapping.tsv"
    query = mapping_path.parent / "vfdb_enterococcus_target_and_context.fna"
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    metadata = read_metadata(project)
    fnas = fna_manifest(project)
    if args.stage == "pilot":
        metadata = metadata.assign(Proxy_group=metadata.HLGR_proxy.eq("HLGR-associated genotype proxy detected"))
        selected = metadata.sort_values("Genome").groupby(["Source", "Proxy_group"], as_index=False).first()
        if len(selected) != 4:
            raise ValueError("Pilot requires four non-empty Source x HLGR-proxy strata")
        genomes = sorted(selected.Genome)
        label = "pilot"
    else:
        genomes = sorted(metadata.Genome)
        label = "72"
    output_dir = project / f"analysis/virulence_adherence/blast_{label}"
    runtimes = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(run_one_blast, genome, fnas[genome], query, output_dir / f"{genome}.tsv", args.overwrite): genome for genome in genomes}
        for future in as_completed(futures):
            runtimes.append(future.result())
    columns = ["qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    raw_frames = []
    for genome in genomes:
        path = output_dir / f"{genome}.tsv"
        try:
            frame = pd.read_csv(path, sep="\t", names=columns)
        except pd.errors.EmptyDataError:
            frame = pd.DataFrame(columns=columns)
        frame.insert(0, "Genome", genome)
        frame.insert(1, "Source", metadata.set_index("Genome").at[genome, "Source"])
        frame["Original_tool_output_reference"] = str(path.relative_to(project))
        raw_frames.append(frame)
    raw = pd.concat(raw_frames, ignore_index=True)
    candidates, statuses = classify_candidates(raw, genomes, mapping)
    statuses = metadata[["Genome", "Source"]].merge(statuses, on="Genome", validate="one_to_many")
    expected_rows = len(genomes) * len(TARGET_ORDER)
    if len(statuses) != expected_rows or statuses.groupby(["Genome", "Target_gene"]).size().max() != 1:
        raise AssertionError("Every screened genome-target pair must have exactly one final status")
    suffix = "pilot" if args.stage == "pilot" else "72"
    atomic_frame(raw, project / f"data/processed/virulence_adherence/virulence_blast_hsp_audit_{suffix}.csv", args.overwrite)
    atomic_frame(candidates, project / f"data/processed/virulence_adherence/virulence_consolidated_candidates_{suffix}.csv", args.overwrite)
    atomic_frame(statuses, project / f"data/processed/virulence_adherence/virulence_target_status_by_genome_{suffix}.csv", args.overwrite)
    runtime = pd.DataFrame(runtimes, columns=["Genome", "Elapsed_seconds"])
    atomic_frame(runtime, project / f"results/tables/virulence_adherence/virulence_blast_runtime_{suffix}.csv", args.overwrite)
    qc = statuses.groupby(["Target_gene", "Detection_status"]).size().reset_index(name="Genome_count")
    atomic_frame(qc, project / f"results/tables/virulence_adherence/virulence_{suffix}_status_qc.csv", args.overwrite)
    # Surface every accepted asa1 locus from ambiguous genomes before prevalence or statistical analysis.
    asa1_ambiguity = ambiguity_locus_table(candidates, statuses, metadata)
    atomic_frame(asa1_ambiguity, project / f"results/tables/virulence_adherence/asa1_ambiguous_loci_{suffix}.csv", args.overwrite)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["reference-audit", "pilot", "screen", "aggregation-review", "resolve"])
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        if args.stage == "reference-audit":
            reference_audit(args)
        elif args.stage == "aggregation-review":
            run_aggregation_review(args)
        elif args.stage == "resolve":
            resolve_for_analysis(args)
        else:
            screen(args)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
