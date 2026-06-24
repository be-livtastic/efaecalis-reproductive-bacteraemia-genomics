import os
import re
import shutil
from pathlib import Path

folders = {
    "02_genomes_fna": ".fna",
    "03_proteins_faa": ".faa",
    "04_gff": ".gff"
}

for folder, new_ext in folders.items():
    folder_path = Path(folder)
    for file in folder_path.iterdir():
        name = file.name

        match = re.search(r"(GC[AF]_\d+\.\d+)", name)
        if match:
            accession = match.group(1)
            new_name = folder_path / f"{accession}{new_ext}"
            if file != new_name:
                shutil.move(str(file), str(new_name))
                print(f"Renamed {file.name} -> {new_name.name}")
        else:
            print(f"Could not find accession in: {file.name}")
