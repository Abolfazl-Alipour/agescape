#!/usr/bin/env python3

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def decompress(source, target):
    if target.exists():
        print(f"[OK] {target}")
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rb") as handle_in, open(target, "wb") as handle_out:
        shutil.copyfileobj(handle_in, handle_out)
    print(f"[DONE] {target}")


def extract(source, target):
    if target.exists():
        print(f"[OK] {target}")
        return
    subprocess.run([
        sys.executable, str(REPO / "scripts" / "extract_seurat_h5ad.py"),
        str(source), str(target),
    ], check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", choices=["all", "afca", "worm", "zebrafish"], default="all")
    args = parser.parse_args()
    if args.only in ("all", "afca"):
        for organ in ("body", "head"):
            decompress(
                REPO / "data" / "raw" / "afca" / f"GSE218661_adata_{organ}_S_v1.0.h5ad.gz",
                REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad",
            )
    if args.only in ("all", "worm"):
        extract(
            REPO / "data" / "raw" / "worm-gao" / "GSE229022_integrated_seurat_revision.rds.gz",
            REPO / "data" / "interim" / "worm-gao-integrated.h5ad",
        )
    if args.only in ("all", "zebrafish"):
        extract(
            REPO / "data" / "raw" / "zebrafish-testis" / "seurat" / "testis_5_12_20_22.rds",
            REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad",
        )
        extract(
            REPO / "data" / "raw" / "zebrafish-testis" / "seurat" / "mo27testis.rds",
            REPO / "data" / "interim" / "zebrafish-testis-27mo.h5ad",
        )


if __name__ == "__main__":
    main()
