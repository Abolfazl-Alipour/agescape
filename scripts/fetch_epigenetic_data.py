#!/usr/bin/env python3
"""
scripts/fetch_epigenetic_data.py

Phase 3: Automated Data Acquisition and Integrity Verification for
Pan-Mammalian Epigenetic Clock (Lu et al. 2023) and GSE184211 / GSE199979.

Downloads and validates:
  1. Human Multi-Tissue Array Metadata and Normalized Betas (GEO GSE184211, N=661)
  2. Mouse Lifespan Liver Matrix (GEO GSE199979, N=339)
  3. Universal Clock 2 Coefficients (Lu et al., Nature Aging 2023)
"""

import sys
import json
import hashlib
from pathlib import Path
import urllib.request

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = REPO_ROOT / "data" / "raw"
REF_DIR = REPO_ROOT / "data" / "reference"
MANIFEST_DIR = REPO_ROOT / "data" / "manifest"

TARGETS = [
    {
        "name": "Human Series Matrix (GSE184211)",
        "dest": RAW_DIR / "GSE184211_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184211/matrix/GSE184211_series_matrix.txt.gz",
        "sha256": "45a861c3936ddb4e30cb96253ed5d14cf8efc1e9593c8ad755f8ff97eb0335a8"
    },
    {
        "name": "Human Normalized Betas (GSE184211)",
        "dest": RAW_DIR / "GSE184211_datBetaNormalized.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184211/suppl/GSE184211_datBetaNormalized.csv.gz",
        "sha256": "3eca1b84c743bc79173cdac1efe34c97c62b0054db54db7709990cc4ab70e860"
    },
    {
        "name": "Mouse Series Matrix (GSE199979)",
        "dest": RAW_DIR / "GSE199979_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199979/matrix/GSE199979_series_matrix.txt.gz",
        "sha256": "2e11bd2c8fcaa749eb25c86485da5a1edf173d34dc147eb1f87d9ef634e97bb4"
    },
    {
        "name": "Universal Clock 2 Coefficients (Lu et al. 2023)",
        "dest": REF_DIR / "Universal_Coef.CombAnn.csv",
        "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00462-6/MediaObjects/43587_2023_462_MOESM4_ESM.csv",
        "sha256": "70cf2e8117bc23e1d9907eb0a253271af2eb8cc23adfac7e7e7e06a4728ab3b3"
    }
]


def verify_sha256(path: Path, expected: str) -> bool:
    if not path.exists():
        return False
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest() == expected


def main():
    print("[Fetch Epigenetics] Checking and acquiring pan-mammalian epigenetic datasets...")
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)

    for item in TARGETS:
        dest = item["dest"]
        expected_sha = item["sha256"]
        if dest.exists() and verify_sha256(dest, expected_sha):
            print(f"  [OK] {item['name']} verified at {dest}")
        else:
            print(f"  [Downloading] {item['name']} from {item['url']} ...")
            try:
                # Add User-Agent header for academic static content servers
                req = urllib.request.Request(item["url"], headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Pipeline)"})
                with urllib.request.urlopen(req) as resp, open(dest, "wb") as out:
                    while chunk := resp.read(1024 * 1024):
                        out.write(chunk)
                if verify_sha256(dest, expected_sha):
                    print(f"    -> Successfully downloaded and verified {dest}")
                else:
                    print(f"    [WARNING] Downloaded {dest} but SHA-256 mismatch!")
            except Exception as e:
                print(f"    [ERROR] Failed to fetch {item['name']}: {e}")

    print("[Fetch Epigenetics] Acquisition check complete.")


if __name__ == "__main__":
    main()
