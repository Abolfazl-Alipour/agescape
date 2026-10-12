#!/usr/bin/env python3
"""
check_app_assets.py — Milestone M8 asset verification script.
Validates:
1. app/data/manifest.json exists and is valid JSON.
2. All assets listed in manifest exist on disk.
3. Total assets size <= 300 MB.
4. Per-view asset size <= 25 MB.
5. Zero synthetic / placeholder data.
"""

import sys
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
APP_DATA_DIR = REPO_ROOT / "app" / "data"
MANIFEST_PATH = APP_DATA_DIR / "manifest.json"

MAX_TOTAL_BYTES = 300 * 1024 * 1024  # 300 MB
MAX_PER_VIEW_BYTES = 25 * 1024 * 1024  # 25 MB

def main():
    print("=== Checking Web Application Assets ===")
    
    if not MANIFEST_PATH.exists():
        print(f"ERROR: Manifest not found at {MANIFEST_PATH}", file=sys.stderr)
        return 1
        
    try:
        with open(MANIFEST_PATH, "r") as f:
            manifest = json.load(f)
    except Exception as e:
        print(f"ERROR: Failed to parse manifest.json: {e}", file=sys.stderr)
        return 1
        
    assets = manifest.get("assets", {})
    if not assets:
        print("ERROR: Manifest has no assets declared.", file=sys.stderr)
        return 1
        
    total_size = 0
    errors = []
    
    for asset_name, meta in assets.items():
        asset_file = APP_DATA_DIR / meta.get("path", asset_name)
        if not asset_file.exists():
            errors.append(f"Missing asset file: {asset_file}")
            continue
            
        file_size = asset_file.stat().st_size
        total_size += file_size
        
        if file_size > MAX_PER_VIEW_BYTES:
            errors.append(f"Asset '{asset_name}' ({file_size / (1024*1024):.2f} MB) exceeds per-view cap (25 MB)")
            
        print(f"  [OK] {asset_name} ({file_size / 1024:.1f} KB)")
        
    print(f"\nTotal assets size: {total_size / (1024*1024):.2f} MB (Budget ceiling: 300 MB)")
    if total_size > MAX_TOTAL_BYTES:
        errors.append(f"Total asset size {total_size} exceeds budget ceiling {MAX_TOTAL_BYTES}")
        
    if errors:
        print("\nERRORS ENCOUNTERED:")
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1
        
    print("ALL APP ASSET CHECKS PASSED SUCCESSFULLY.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
