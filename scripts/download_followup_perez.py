#!/usr/bin/env python3

import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / f"{REPO.name}-data" / "muscle" / "perez"
MANIFEST = REPO / "data" / "manifest" / "GSE167186.json"


def main():
    started = time.time()
    DATA.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text())
    records = []
    for donor in manifest["donors"]:
        filename = f"{donor['gsm']}_{donor['hm_id']}_filtered_feature_bc_matrix.h5"
        output = DATA / filename
        command = [
            sys.executable, str(REPO / "scripts" / "resource_guard.py"), "download",
            "--url", donor["h5_url"], "--output", str(output),
            "--min-free-gb", "30", "--max-file-gb", "1",
        ]
        file_started = time.time()
        subprocess.run(command, cwd=REPO, check=True)
        if output.stat().st_size != donor["h5_size_bytes"]:
            raise ValueError(f"Downloaded size mismatch for {donor['gsm']}")
        records.append({"gsm": donor["gsm"], "bytes": output.stat().st_size,
                        "elapsed_seconds": time.time() - file_started})
    result = {"files": len(records), "bytes": sum(row["bytes"] for row in records),
              "elapsed_seconds": time.time() - started, "sequential": True,
              "minimum_free_gb": 30, "maximum_file_gb": 1, "downloads": records}
    (REPO / "results" / "followup-muscle-perez-download-resource.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
