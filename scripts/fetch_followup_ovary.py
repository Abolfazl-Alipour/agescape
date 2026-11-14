#!/usr/bin/env python3
"""Sequentially download the frozen Step 7 ovary expression inputs."""

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / f"{REPO.name}-data" / "ovary"
COHORTS = REPO / "data" / "manifest" / "ovary_cohorts.json"
RESULTS = REPO / "results"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def guarded_download(url, output, expected_size, cap_gb):
    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.exists() or output.stat().st_size != expected_size:
        subprocess.run([
            sys.executable, str(REPO / "scripts" / "resource_guard.py"), "download",
            "--url", url, "--output", str(output), "--min-free-gb", "30",
            "--max-file-gb", str(cap_gb),
        ], check=True)
    if output.stat().st_size != expected_size:
        raise ValueError(f"Downloaded size mismatch for {output}")
    return {
        "url": url,
        "local_path": str(output),
        "size_bytes": output.stat().st_size,
        "sha256": sha256(output),
        "guarded_download": True,
        "maximum_file_gb": cap_gb,
        "minimum_free_gb": 30,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True, choices=["GSE202601", "GSE329417"])
    args = parser.parse_args()
    started = time.time()
    manifest = json.loads(COHORTS.read_text())
    files = []
    if args.cohort == "GSE202601":
        for donor in manifest[args.cohort]["donors"]:
            entry = donor["rna_filtered_count_matrix"]
            filename = entry["url"].rsplit("/", 1)[-1]
            output = DATA / args.cohort / filename
            print(f"Jin {donor['donor']}: {filename}", flush=True)
            files.append({"donor": donor["donor"], **guarded_download(entry["url"], output, entry["size_bytes"], 1)})
    else:
        entry = next(row for row in manifest[args.cohort]["public_processed_objects"]
                     if row["url"].endswith("final_annotated_atlas.h5ad.gz"))
        output = DATA / args.cohort / entry["url"].rsplit("/", 1)[-1]
        files.append(guarded_download(entry["url"], output, entry["size_bytes"], 7))
    result = {
        "cohort": args.cohort,
        "execution_order": 1 if args.cohort == "GSE202601" else 2,
        "files": files,
        "download_bytes": sum(row["size_bytes"] for row in files),
        "raw_tar_downloaded": False,
        "fastq_downloaded": False,
        "sequential_resource_guarded_downloads": True,
        "elapsed_seconds": time.time() - started,
    }
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / f"followup-ovary-{args.cohort}-download.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
