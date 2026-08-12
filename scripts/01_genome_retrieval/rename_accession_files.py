"""Rename files inside an extracted NCBI dataset to accession-based names.

Pass ``data/raw/ncbi_genomes`` (or one category below it). The default is a
dry run; ``--apply`` is required to rename files. Existing files are never
overwritten.
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

# --- Command-line interface ---
parser = argparse.ArgumentParser()
parser.add_argument("dataset_dir", type=Path)
parser.add_argument(
    "--apply", action="store_true", help="Perform renames; default is dry run."
)
args = parser.parse_args()
dataset_dir = args.dataset_dir.resolve()

# --- Discover sequence and annotation files recursively ---
if not dataset_dir.is_dir():
    sys.exit(f"Dataset directory does not exist: {dataset_dir}")
extensions = {".fna", ".faa", ".gff"}
files = sorted(
    path for path in dataset_dir.rglob("*")
    if path.is_file() and path.suffix.lower() in extensions
)
if not files:
    sys.exit(f"No .fna, .faa or .gff files found under: {dataset_dir}")

# --- Preview or apply collision-safe renames ---
for file in files:
    match = re.search(r"(GC[AF]_\d+\.\d+)", str(file))
    if not match:
        print(f"Could not find accession in path: {file}")
        continue
    accession = match.group(1)
    destination = file.with_name(f"{accession}{file.suffix.lower()}")
    if file == destination:
        continue
    if destination.exists():
        sys.exit(f"Refusing to overwrite existing destination: {destination}")
    if args.apply:
        shutil.move(str(file), str(destination))
        print(f"Renamed {file.relative_to(dataset_dir)} -> {destination.name}")
    else:
        print(f"Would rename {file.relative_to(dataset_dir)} -> {destination.name}")
