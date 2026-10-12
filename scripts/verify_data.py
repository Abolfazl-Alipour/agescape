#!/usr/bin/env python3
"""
verify_data.py — Verifies data files against manifests in data/manifest/*.json.
Exits 0 if all files match SHA-256, non-zero if missing or mismatched.
"""

import hashlib
import json
import os
import sys
from pathlib import Path


def compute_sha256(file_path: Path, chunk_size: int = 65536) -> str:
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def main():
    repo_root = Path(__file__).resolve().parent.parent
    manifest_dir = repo_root / "data" / "manifest"
    raw_dir = repo_root / "data" / "raw"

    manifest_files = sorted(manifest_dir.glob("*.json"))
    if not manifest_files:
        print(f"Error: No manifests found in {manifest_dir}")
        sys.exit(1)

    print(f"{'Filename':<35} {'Size (MB)':<12} {'Status':<10} {'SHA-256 match'}")
    print("-" * 80)

    all_passed = True
    for m_path in manifest_files:
        with open(m_path, "r") as f:
            meta = json.load(f)

        filename = meta["filename"]
        expected_sha = meta["sha256"]
        expected_size = meta["size_bytes"]
        target_path = raw_dir / filename

        if not target_path.exists():
            print(f"{filename:<35} {'MISSING':<12} {'FAIL':<10} File not found")
            all_passed = False
            continue

        actual_size = target_path.stat().st_size
        actual_sha = compute_sha256(target_path)
        size_mb = f"{actual_size / (1024 * 1024):.1f}"

        if actual_sha.lower() == expected_sha.lower():
            status = "PASS"
            match_str = f"OK ({actual_sha[:8]}...)"
        else:
            status = "FAIL"
            match_str = f"MISMATCH (expected {expected_sha[:8]}, got {actual_sha[:8]})"
            all_passed = False

        print(f"{filename:<35} {size_mb:<12} {status:<10} {match_str}")

    print("-" * 80)
    if all_passed:
        print("All raw datasets verified successfully.")
        sys.exit(0)
    else:
        print("Data verification FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
