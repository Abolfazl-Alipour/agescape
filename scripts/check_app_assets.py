#!/usr/bin/env python3
"""
check_app_assets.py — verifies the static web app assets.
In Spec v2.0, delegates directly to scripts/check_app_data.py which verifies
all 105 anchors, file caps, and SHA-256 provenance against real data.
"""

import sys
from pathlib import Path

# Add repo root to sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from check_app_data import main

if __name__ == "__main__":
    sys.exit(main())

