#!/usr/bin/env python3
"""Build a deterministic file-level repository cleanup manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path


PCA_PRE_CLEANUP_HEAD_PATHS = {
    "scripts/06_phylogenomics/10_plot_alignment_pca.R",
    "results/figures/phylogeny_72_genomes_9_locus_observed_indels/pca/pca_pc1_pc2.png",
    "results/figures/phylogeny_72_genomes_9_locus_observed_indels/pca/pca_pc1_pc2.pdf",
    "results/figures/phylogeny_72_genomes_9_locus_observed_indels/pca/pca_pc1_pc3.png",
    "results/figures/phylogeny_72_genomes_9_locus_observed_indels/pca/pca_pc1_pc3.pdf",
}


def git_lines(root: Path, *args: str) -> set[str]:
    result = subprocess.run(
        ["git", *args, "-z"], cwd=root, check=True, capture_output=True
    )
    return {
        item.decode("utf-8", errors="surrogateescape")
        for item in result.stdout.split(b"\0")
        if item
    }


def classify(path: str) -> str:
    parts = Path(path).parts
    name = Path(path).name

    if (
        "__pycache__" in parts
        or ".Rproj.user" in parts
        or path.startswith("local_archive/local_session/")
        or name.endswith((".pyc", ".pyo"))
        or name.endswith(".zero_byte_transfer_fragment")
    ):
        return "temporary"

    obsolete_prefixes = (
        "analysis/phylogenomics/9_locus/",
        "results/tables/phylogeny_9_locus/",
        "local_archive/large_outputs/AMR_analysis_outputs/",
        "local_archive/superseded_versions/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/02_rpob_fastas/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/03_clustalo/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/04_iqtree/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/05_quality_checks/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/06_logs/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/07_tree_figures/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/08_improved_tree_figures/",
    )
    obsolete_files = {
        "scripts/03_amr_analysis/amrfinder_analysis_72_genomes.R",
        "scripts/03_amr_analysis/isolate_summary.py",
        "local_archive/exploratory_scripts/setup script",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/phylogeny_run_summary.txt",
    }
    if path.startswith(obsolete_prefixes) or path in obsolete_files:
        return "obsolete"
    if "_assembly_accession_labels." in name:
        return "obsolete"

    duplicate_prefixes = (
        "local_archive/large_outputs/Phylogeny_preparation/",
        "local_archive/large_outputs/Phylogeny_project/combined_72_genomes/",
        "local_archive/large_outputs/Phylogeny_project/genomes/",
        "analysis/virulence_adherence/pilot_initial_failed_consolidation/",
    )
    if path.startswith(duplicate_prefixes):
        return "duplicate"

    provenance_prefixes = (
        "analysis/logs/",
        "backups/",
        "environment/",
        "local_archive/original_metadata/",
        "local_archive/large_outputs/Efaecalis_14_AMRFinder/01_ncbi_download/",
        "local_archive/large_outputs/Efaecalis_14_AMRFinder/05_amrfinder_results/",
        "local_archive/large_outputs/Efaecalis_58_AMRFinder/01_ncbi_download/",
        "local_archive/large_outputs/Efaecalis_58_AMRFinder/05_amrfinder_results/",
        "local_archive/large_outputs/Phylogeny_project/phylogeny_work/01_annotations/",
        "references/",
    )
    if path.startswith(provenance_prefixes):
        return "provenance"

    return "canonical"


def digest_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def digest_head(root: Path, path: str) -> tuple[int, str]:
    content = subprocess.run(
        ["git", "show", f"HEAD:{path}"], cwd=root, check=True, capture_output=True
    ).stdout
    return len(content), hashlib.sha256(content).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--git-status-output", type=Path)
    parser.add_argument("--git-diff-stat-output", type=Path)
    parser.add_argument(
        "--pre-cleanup-pca-from-head",
        action="store_true",
        help="Use HEAD content for the five PCA paths changed after baseline capture.",
    )
    return parser.parse_args()


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        handle.write(content)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    output = args.output if args.output.is_absolute() else root / args.output
    audit_prefix = "results/tables/repository_cleanup_audit/"
    tool_path = Path(__file__).resolve().relative_to(root).as_posix()

    tracked = git_lines(root, "ls-files")
    untracked = git_lines(root, "ls-files", "--others", "--exclude-standard")
    modified = git_lines(root, "diff", "--name-only") | git_lines(
        root, "diff", "--cached", "--name-only"
    )

    paths = []
    for candidate in root.rglob("*"):
        if not candidate.is_file() or candidate.is_symlink():
            continue
        relative = candidate.relative_to(root).as_posix()
        if relative.startswith(".git/") or relative.startswith(audit_prefix):
            continue
        # The manifest builder is audit machinery introduced after the baseline.
        if relative == tool_path:
            continue
        paths.append(relative)

    rows = []
    for relative in sorted(paths):
        use_head = args.pre_cleanup_pca_from_head and relative in PCA_PRE_CLEANUP_HEAD_PATHS
        size, sha256 = (
            digest_head(root, relative) if use_head else digest_file(root / relative)
        )
        if relative in tracked:
            git_status = "tracked_modified" if relative in modified else "tracked_clean"
        elif relative in untracked:
            git_status = "untracked"
        else:
            git_status = "ignored"
        rows.append(
            {
                "path": relative,
                "file_size_bytes": size,
                "sha256": sha256,
                "git_status": git_status,
                "classification": classify(relative),
                "content_source": "git_HEAD" if use_head else "working_tree",
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=output.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "path",
                "file_size_bytes",
                "sha256",
                "git_status",
                "classification",
                "content_source",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, output)
    if args.git_status_output:
        status_path = args.git_status_output if args.git_status_output.is_absolute() else root / args.git_status_output
        status = subprocess.run(
            ["git", "status", "--short"], cwd=root, check=True, capture_output=True, text=True
        ).stdout
        atomic_text(status_path, status)
    if args.git_diff_stat_output:
        stat_path = args.git_diff_stat_output if args.git_diff_stat_output.is_absolute() else root / args.git_diff_stat_output
        stat = subprocess.run(
            ["git", "diff", "--stat"], cwd=root, check=True, capture_output=True, text=True
        ).stdout
        atomic_text(stat_path, stat)
    print(f"Wrote {len(rows)} file records to {output.relative_to(root)}")


if __name__ == "__main__":
    main()
