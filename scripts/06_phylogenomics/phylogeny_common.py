#!/usr/bin/env python3
"""Shared, dependency-light helpers for the nine-locus phylogeny pipeline."""
from __future__ import annotations

import csv
import io
import re
import tempfile
from pathlib import Path
from urllib.parse import unquote

ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


# --- Identifier and tabular I/O helpers ---
def accession_from_text(value: str) -> str | None:
    match = ACCESSION_RE.search(value)
    return match.group(1) if match else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    # Some retained Windows-authored tables contain CR immediately before a tab
    # rather than as part of CRLF. Normalize in memory; never alter the source.
    # Decode bytes directly so Python's universal-newline handling cannot turn
    # the embedded CR into a row boundary before it is normalized.
    content = path.read_bytes().decode("utf-8-sig").replace("\r", "")
    return list(csv.DictReader(io.StringIO(content), delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing file: {path}")
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile("w", newline="", encoding="utf-8", dir=path.parent,
                                         prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        temporary.rename(path)
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


# --- GFF attribute parsing and normalization ---
def parse_attributes(text: str) -> dict[str, str]:
    """Parse GFF3 attributes without assuming attribute order."""
    result: dict[str, str] = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        if sep:
            result[unquote(key)] = unquote(value)
    return result


def normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


# --- FASTA readers shared by validation stages ---
def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in lengths:
                    raise ValueError(f"Duplicate FASTA identifier {current} in {path}")
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())
    return lengths


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA identifier {current} in {path}")
                records[current] = []
            elif current is not None:
                records[current].append(line.strip().upper())
    return {key: "".join(parts) for key, parts in records.items()}


# --- Atomic text output for interruption-safe results ---
def atomic_text(path: Path, content: str) -> None:
    """Create a new file atomically; never replace an existing destination."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing file: {path}")
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                         prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        temporary.rename(path)
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
