#!/usr/bin/env python3
"""
fetch_data.py — Downloads missing datasets declared in data/manifest/*.json.
Downloads to data/raw/.partial/ using curl -C -, verifies SHA-256, and renames into place.
Supports optional --only <name> filter.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def compute_sha256(file_path: Path, chunk_size: int = 65536) -> str:
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Fetch missing data defined in manifests.")
    parser.add_argument("--only", type=str, default=None, help="Download only the specified filename")
    parser.add_argument("--include-optional", action="store_true", help="Also fetch datasets marked optional in the manifest")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    manifest_dir = repo_root / "data" / "manifest"
    raw_dir = repo_root / "data" / "raw"
    partial_dir = raw_dir / ".partial"

    raw_dir.mkdir(parents=True, exist_ok=True)
    partial_dir.mkdir(parents=True, exist_ok=True)

    manifest_files = sorted(manifest_dir.glob("*.json"))
    if not manifest_files:
        print(f"Error: No manifests found in {manifest_dir}")
        sys.exit(1)

    for m_path in manifest_files:
        with open(m_path, "r") as f:
            meta = json.load(f)

        filename = meta["filename"]
        if args.only and filename != args.only:
            continue

        if meta.get("optional") and not args.include_optional and not args.only:
            print(f"[SKIP] {filename} is optional (use --include-optional).")
            continue

        target_path = raw_dir / filename
        expected_sha = meta["sha256"]
        url = meta.get("url")

        if target_path.exists():
            actual_sha = compute_sha256(target_path)
            if actual_sha.lower() == expected_sha.lower():
                print(f"[OK] {filename} exists and SHA-256 verified.")
                continue
            else:
                print(f"[WARN] {filename} exists but SHA-256 mismatch. Re-downloading...")

        if not url:
            print(f"[ERROR] No URL provided in manifest for {filename} and file is missing/mismatched.")
            sys.exit(1)

        partial_path = partial_dir / filename
        print(f"[FETCH] Downloading {filename} from {url}...")
        cmd = ["curl", "-f", "-C", "-", "-L", "--retry", "5", "--retry-delay", "2", "-o", str(partial_path), url]
        res = subprocess.run(cmd)
        if res.returncode != 0:
            print(f"[ERROR] Download failed for {filename} (curl exit code {res.returncode})")
            partial_path.unlink(missing_ok=True)
            sys.exit(1)

        downloaded_sha = compute_sha256(partial_path)
        if downloaded_sha.lower() != expected_sha.lower():
            print(f"[ERROR] Checksum mismatch for downloaded file {filename}:")
            print(f"  Expected: {expected_sha}")
            print(f"  Got:      {downloaded_sha}")
            partial_path.unlink(missing_ok=True)
            sys.exit(1)

        partial_path.rename(target_path)
        print(f"[DONE] Successfully downloaded and verified {filename}")


if __name__ == "__main__":
    main()
