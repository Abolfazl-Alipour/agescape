#!/usr/bin/env python3
"""
check_app_assets.py — verifies the static web app assets against their manifest.

Checks (exit 0 only if all pass, else exit 1 with diagnostics):
  1. app/data/manifest.json exists and parses.
  2. Every asset listed in the manifest exists under app/ and is <= 1 MB.
  3. Every source file listed for an asset exists (relative to the repo root)
     and its SHA-256 matches the digest recorded at generation time —
     this is what guarantees every displayed number traces to results/*.csv.
  4. Total size of files under app/data/ is <= 5 MB.

Run: python scripts/check_app_assets.py
"""

import hashlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
APP_DIR = REPO_ROOT / "app"
APP_DATA_DIR = APP_DIR / "data"
MANIFEST_PATH = APP_DATA_DIR / "manifest.json"

MAX_ASSET_BYTES = 1_048_576          # 1 MB per individual asset file
MAX_DATA_TOTAL_BYTES = 5_242_880     # 5 MB total under app/data


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    print("=== Checking Web Application Assets ===")
    errors = []

    # 1. manifest exists and parses
    if not MANIFEST_PATH.exists():
        print(f"ERROR: manifest not found at {MANIFEST_PATH}", file=sys.stderr)
        return 1
    try:
        with open(MANIFEST_PATH, "r") as f:
            manifest = json.load(f)
    except Exception as e:
        print(f"ERROR: failed to parse manifest.json: {e}", file=sys.stderr)
        return 1
    assets = manifest.get("assets")
    if not isinstance(assets, dict) or not assets:
        print("ERROR: manifest has no 'assets' mapping.", file=sys.stderr)
        return 1
    print(f"manifest: generated_by={manifest.get('generated_by')} "
          f"at {manifest.get('generated_at_utc')} ({len(assets)} assets)")

    # 2. + 3. every listed asset and source
    manifest_asset_paths = set()
    for name, meta in assets.items():
        rel = meta.get("path")
        if not rel:
            errors.append(f"asset '{name}' has no 'path'")
            continue
        manifest_asset_paths.add(rel)
        asset_file = APP_DIR / rel
        if not asset_file.exists():
            errors.append(f"asset '{name}': file missing at app/{rel}")
        else:
            size = asset_file.stat().st_size
            if size > MAX_ASSET_BYTES:
                errors.append(
                    f"asset '{name}' ({size} bytes) exceeds the 1 MB per-file cap")
            print(f"  [OK] {name:42s} app/{rel} ({size / 1024:.1f} KB)")

        sources = meta.get("sources", [])
        if not sources:
            errors.append(f"asset '{name}' lists no sources (provenance gap)")
        for s in sources:
            src_rel = s.get("path")
            expected = s.get("sha256")
            if not src_rel or not expected:
                errors.append(f"asset '{name}': source entry missing path/sha256")
                continue
            src_file = REPO_ROOT / src_rel
            if not src_file.exists():
                errors.append(f"asset '{name}': source file missing: {src_rel}")
                continue
            actual = sha256_of(src_file)
            if actual != expected:
                errors.append(
                    f"asset '{name}': SHA-256 mismatch for {src_rel}\n"
                    f"    manifest: {expected}\n"
                    f"    actual:   {actual}\n"
                    f"    -> re-run notebooks/07-app-assets.ipynb")
            else:
                print(f"       src {src_rel} sha256 OK")

    # unlisted files under app/data are not covered by provenance
    for p in sorted(APP_DATA_DIR.iterdir()):
        if p.is_file() and p.name != "manifest.json":
            rel = f"data/{p.name}"
            if rel not in manifest_asset_paths:
                errors.append(f"app/data/{p.name} is not listed in the manifest "
                              f"(stale or untracked asset — regenerate or remove)")

    # 4. total budget under app/data
    total = sum(p.stat().st_size for p in APP_DATA_DIR.iterdir() if p.is_file())
    print(f"\nTotal app/data size: {total / 1024:.1f} KB (cap 5 MB)")
    if total > MAX_DATA_TOTAL_BYTES:
        errors.append(f"app/data total {total} bytes exceeds 5 MB cap")

    if errors:
        print("\nERRORS ENCOUNTERED:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1
    print("ALL APP ASSET CHECKS PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
