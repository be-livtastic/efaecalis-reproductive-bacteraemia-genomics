#!/usr/bin/env python3
"""Compare repository manifests and write a deterministic cleanup ledger."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


FIELDS = [
    "action",
    "old_path",
    "new_path",
    "old_size_bytes",
    "new_size_bytes",
    "old_sha256",
    "new_sha256",
    "classification",
]


def read_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["path"]: row for row in csv.DictReader(handle, delimiter="\t")}


def read_mappings(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["old_path"]: row["new_path"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--before", required=True, type=Path)
    parser.add_argument("--after", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--mappings", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    before = read_tsv(args.before)
    after = read_tsv(args.after)
    mappings = read_mappings(args.mappings)
    rows: list[dict[str, str]] = []
    handled_before: set[str] = set()
    handled_after: set[str] = set()

    for old_path, new_path in sorted(mappings.items()):
        if old_path not in before:
            raise ValueError(f"Mapped old path is absent from before manifest: {old_path}")
        if new_path not in after:
            raise ValueError(f"Mapped new path is absent from after manifest: {new_path}")
        old, new = before[old_path], after[new_path]
        action = "moved" if old["sha256"] == new["sha256"] else "moved_modified"
        rows.append(
            {
                "action": action,
                "old_path": old_path,
                "new_path": new_path,
                "old_size_bytes": old["file_size_bytes"],
                "new_size_bytes": new["file_size_bytes"],
                "old_sha256": old["sha256"],
                "new_sha256": new["sha256"],
                "classification": new["classification"],
            }
        )
        handled_before.add(old_path)
        handled_after.add(new_path)

    for path in sorted(before.keys() & after.keys()):
        if path in handled_before or path in handled_after:
            continue
        old, new = before[path], after[path]
        if old["sha256"] == new["sha256"]:
            continue
        rows.append(
            {
                "action": "modified",
                "old_path": path,
                "new_path": path,
                "old_size_bytes": old["file_size_bytes"],
                "new_size_bytes": new["file_size_bytes"],
                "old_sha256": old["sha256"],
                "new_sha256": new["sha256"],
                "classification": new["classification"],
            }
        )

    for path in sorted(before.keys() - after.keys() - handled_before):
        old = before[path]
        rows.append(
            {
                "action": "removed",
                "old_path": path,
                "new_path": "",
                "old_size_bytes": old["file_size_bytes"],
                "new_size_bytes": "",
                "old_sha256": old["sha256"],
                "new_sha256": "",
                "classification": old["classification"],
            }
        )

    for path in sorted(after.keys() - before.keys() - handled_after):
        new = after[path]
        rows.append(
            {
                "action": "added",
                "old_path": "",
                "new_path": path,
                "old_size_bytes": "",
                "new_size_bytes": new["file_size_bytes"],
                "old_sha256": "",
                "new_sha256": new["sha256"],
                "classification": new["classification"],
            }
        )

    rows.sort(key=lambda row: (row["action"], row["old_path"], row["new_path"]))
    write_tsv(args.output, rows, FIELDS)

    actions = sorted({row["action"] for row in rows})
    summary = [
        {"metric": "before_file_count", "value": str(len(before))},
        {"metric": "after_file_count", "value": str(len(after))},
        {"metric": "changed_file_records", "value": str(len(rows))},
    ]
    summary.extend(
        {
            "metric": f"action_{action}",
            "value": str(sum(row["action"] == action for row in rows)),
        }
        for action in actions
    )
    summary.extend(
        [
            {
                "metric": "removed_bytes",
                "value": str(
                    sum(
                        int(row["old_size_bytes"] or 0)
                        for row in rows
                        if row["action"] == "removed"
                    )
                ),
            },
            {
                "metric": "added_bytes",
                "value": str(
                    sum(
                        int(row["new_size_bytes"] or 0)
                        for row in rows
                        if row["action"] == "added"
                    )
                ),
            },
        ]
    )
    write_tsv(args.summary, summary, ["metric", "value"])
    print(f"Wrote {len(rows)} changed file records to {args.output}")


if __name__ == "__main__":
    main()
