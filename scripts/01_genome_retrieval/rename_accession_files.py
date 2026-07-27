"""Rename downloaded NCBI files to stable accession-based names safely.

Run this after retrieval and pass the dataset directory explicitly. The
default is a dry run; --apply is required to move files. Existing destination
files are never overwritten.
"""

import argparse
import re
import shutil
import sys
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("dataset_dir", type=Path)
parser.add_argument(
    "--apply", action="store_true", help="Perform renames; default is dry run."
)
args = parser.parse_args()
dataset_dir = args.dataset_dir.resolve()

folders = {
    "02_genomes_fna": ".fna",
    "03_proteins_faa": ".faa",
    "04_gff": ".gff"
}

for folder, new_ext in folders.items():
    folder_path = dataset_dir / folder
    if not folder_path.is_dir():
        sys.exit(f"Required input directory does not exist: {folder_path}")
    for file in sorted(folder_path.iterdir()):
        if not file.is_file():
            continue
        name = file.name

        match = re.search(r"(GC[AF]_\d+\.\d+)", name)
        if match:
            accession = match.group(1)
            new_name = folder_path / f"{accession}{new_ext}"
            if file != new_name:
                if new_name.exists():
                    sys.exit(
                        f"Refusing to overwrite existing destination: {new_name}"
                    )
                if args.apply:
                    shutil.move(str(file), str(new_name))
                    print(f"Renamed {file.name} -> {new_name.name}")
                else:
                    print(f"Would rename {file.name} -> {new_name.name}")
        else:
            print(f"Could not find accession in: {file.name}")
