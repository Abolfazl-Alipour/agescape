#!/usr/bin/env python3
"""Stream-decompress and inspect the frozen GSE329417 atlas without outcomes."""

import argparse
import gzip
import hashlib
import json
import resource
import shutil
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / f"{REPO.name}-data" / "ovary" / "GSE329417"
COMPRESSED = DATA / "GSE329417_human_ovary_snucRNAseq_final_annotated_atlas.h5ad.gz"
ATLAS = DATA / "GSE329417_human_ovary_snucRNAseq_final_annotated_atlas.h5ad"
RESULTS = REPO / "results"
PREFIX = "followup-ovary-postmenopausal"
EXPECTED_COMPRESSED = 5746830392
MAX_DECOMPRESSED = 40 * 1024**3
MIN_FREE = 30 * 1024**3


def decompress():
    started = time.time()
    if not COMPRESSED.is_file() or COMPRESSED.stat().st_size != EXPECTED_COMPRESSED:
        raise ValueError("The complete guarded atlas download is absent or has an unexpected size")
    if ATLAS.exists() and ATLAS.stat().st_size > 0:
        print(f"Atlas already decompressed: {ATLAS}", flush=True)
    else:
        partial = ATLAS.with_suffix(".h5ad.partial")
        digest = hashlib.sha256()
        written = 0
        try:
            with gzip.open(COMPRESSED, "rb") as source, partial.open("wb") as target:
                while True:
                    block = source.read(8 * 1024 * 1024)
                    if not block:
                        break
                    written += len(block)
                    if written > MAX_DECOMPRESSED:
                        raise ValueError("Decompressed atlas exceeds the frozen 40 GiB cap")
                    if shutil.disk_usage(DATA).free < MIN_FREE:
                        raise ValueError("Decompression would violate the 30 GiB disk reserve")
                    target.write(block)
                    digest.update(block)
                    if written % (1024**3) < 8 * 1024 * 1024:
                        print(f"decompressed {written / 1024**3:.2f} GiB", flush=True)
            partial.replace(ATLAS)
        except Exception:
            partial.unlink(missing_ok=True)
            raise
    digest = hashlib.sha256()
    with ATLAS.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    result = {
        "compressed_path": str(COMPRESSED),
        "compressed_size_bytes": COMPRESSED.stat().st_size,
        "atlas_path": str(ATLAS),
        "decompressed_size_bytes": ATLAS.stat().st_size,
        "decompressed_sha256": digest.hexdigest(),
        "maximum_decompressed_gb": 40,
        "minimum_free_gb": 30,
        "memory_limit_gb": 8,
        "elapsed_seconds": time.time() - started,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
    }
    path = RESULTS / f"{PREFIX}-decompression.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"wrote {path}", flush=True)


def inspect():
    started = time.time()
    with h5py.File(ATLAS, "r") as handle:
        obs_group = handle["obs"]
        n_obs = len(handle["X/indptr"]) - 1
        n_vars = len(handle["var/_index"])
        columns = {}
        for column in obs_group:
            item = obs_group[column]
            if isinstance(item, h5py.Group) and "categories" in item:
                categories = item["categories"][:].astype(str)
                codes = item["codes"][:]
                unique = categories.tolist()
                n_unique = len(unique)
                dtype = "categorical"
            else:
                values = item[:]
                unique = np.unique(values.astype(str)).tolist()
                n_unique = len(unique)
                dtype = str(values.dtype)
            columns[column] = {
                "dtype": dtype,
                "n_unique": int(n_unique),
                "values": unique if len(unique) <= 100 else unique[:100],
            }
        obs_names = obs_group["_index"][:].astype(str)
        var_names = handle["var/_index"][:].astype(str)
        layers = {}
        for name, layer in handle["layers"].items():
            layers[name] = {
                "type": "CSR" if isinstance(layer, h5py.Group) else str(type(layer)),
                "dtype": str(layer["data"].dtype) if isinstance(layer, h5py.Group) else str(layer.dtype),
                "shape": [n_obs, n_vars],
            }
        schema = {
            "n_obs": int(n_obs),
            "n_vars": int(n_vars),
            "obs_names_unique": bool(pd.Index(obs_names).is_unique),
            "obs_names_examples": obs_names[:10].tolist(),
            "obs_columns": columns,
            "var_names_unique": bool(pd.Index(var_names).is_unique),
            "var_names_examples": var_names[:20].tolist(),
            "var_columns": list(handle["var"].keys()),
            "layers": layers,
            "x_type": "CSR",
            "x_dtype": str(handle["X/data"].dtype),
            "elapsed_seconds": time.time() - started,
            "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
            "outcomes_examined": False,
        }
    path = RESULTS / f"{PREFIX}-schema.json"
    path.write_text(json.dumps(schema, indent=2, sort_keys=True) + "\n")
    print(json.dumps(schema, indent=2, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["decompress", "inspect"])
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    if args.command == "decompress":
        decompress()
    else:
        inspect()


if __name__ == "__main__":
    main()
