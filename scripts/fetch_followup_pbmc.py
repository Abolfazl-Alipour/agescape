#!/usr/bin/env python3
"""Fetch per-donor PBMC files sequentially through the resource guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parent.parent
DEFAULT_DATA = REPO.parent / "agescape-revision-data" / "pbmc"


def guarded_download(url: str, output: Path, expected_size: int) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.stat().st_size == expected_size:
        print(f"Already complete: {output}", flush=True)
    else:
        if output.exists():
            output.unlink()
        command = [
            sys.executable, str(REPO / "scripts" / "resource_guard.py"), "download",
            "--url", url, "--output", str(output), "--min-free-gb", "30",
            "--max-file-gb", "3",
        ]
        subprocess.run(command, check=True)
    actual = output.stat().st_size
    if actual != expected_size:
        raise ValueError(f"Size mismatch for {output}: {actual} != {expected_size}")
    return {"url": url, "path": str(output), "size_bytes": actual}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True, choices=["GSE213516", "GSE214546"])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    manifest = json.loads((REPO / "data" / "manifest" / "pbmc_geo_cohorts.json").read_text())
    cohort = manifest[args.cohort]
    records = []
    if args.cohort == "GSE213516":
        annotation = cohort["annotations"]
        records.append({
            "kind": "annotation", "donor": None,
            **guarded_download(annotation["url"], args.data_root / "raw" / args.cohort / annotation["file"], annotation["size_bytes"]),
        })
    for donor in cohort["donors"]:
        kinds = ["h5"] if args.cohort == "GSE214546" else ["matrix", "barcodes", "features"]
        for kind in kinds:
            source = donor[kind]
            filename = source["url"].rsplit("/", 1)[-1]
            output = args.data_root / "raw" / args.cohort / donor["donor"] / filename
            records.append({
                "kind": kind, "donor": donor["donor"], "gsm": donor["gsm"],
                **guarded_download(source["url"], output, source["size_bytes"]),
            })
    inventory = args.data_root / "interim" / args.cohort
    inventory.mkdir(parents=True, exist_ok=True)
    path = inventory / "download-inventory.json"
    path.write_text(json.dumps({"cohort": args.cohort, "files": records}, indent=2) + "\n")
    print(f"Wrote {path}; {len(records)} files verified by byte size", flush=True)


if __name__ == "__main__":
    main()
