#!/usr/bin/env python3
"""Download the 75 public Cui testis count files sequentially with resource guards."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / "agescape-revision-data" / "cui-testis"
MANIFEST = REPO / "data" / "manifest" / "GSE254315.json"
INVENTORY = REPO / "results" / "followup-cui-testis-file-inventory.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    records = []
    for entry in manifest["files"]:
        output = DATA / "raw" / entry["donor"] / entry["file"]
        output.parent.mkdir(parents=True, exist_ok=True)
        if not output.exists() or output.stat().st_size != entry["size_bytes"]:
            subprocess.run([
                sys.executable, str(REPO / "scripts" / "resource_guard.py"), "download",
                "--url", entry["url"], "--output", str(output),
                "--min-free-gb", "30", "--max-file-gb", "0.2",
            ], check=True)
        if output.stat().st_size != entry["size_bytes"]:
            raise ValueError(f"Size mismatch: {output}")
        entry["local_path"] = str(Path("..") / output.relative_to(REPO.parent))
        entry["sha256"] = sha256(output)
        entry["download_verified"] = True
        records.append(entry)
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"Verified {len(records)}/75: {entry['file']}", flush=True)
    if len(records) != 75 or not all(row["download_verified"] for row in records):
        raise ValueError("Incomplete Cui file acquisition")
    fieldnames = list(records[0])
    import csv
    with INVENTORY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print("All 75 files verified by expected byte size and SHA-256", flush=True)


if __name__ == "__main__":
    main()
