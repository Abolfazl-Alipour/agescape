#!/usr/bin/env python3
"""Split the combined TMS QC-cleaned matrix into mechanical per-tissue partitions."""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DEFAULT_SOURCE = Path("/home/reza/projects/agescape/data/interim/m1-qc-cleaned.h5ad")
DEFAULT_OUTPUT = REPO / "data" / "interim"
EXPECTED_COUNTS = {
    "Marrow": 39471,
    "Spleen": 34989,
    "Limb_Muscle": 28760,
    "Lung": 24458,
    "Kidney": 20922,
    "Tongue": 20031,
    "Mammary_Gland": 12064,
    "Thymus": 9173,
    "Bladder": 8669,
    "Heart_and_Aorta": 8335,
    "Trachea": 7898,
    "Liver": 7084,
    "Fat": 6383,
    "Pancreas": 4212,
    "Skin": 3962,
    "Large_Intestine": 1775,
}
EXPECTED_CELLS = 238186
EXPECTED_GENES = 20138


def split(source: Path, output_dir: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(source, backed="r")
    if adata.shape != (EXPECTED_CELLS, EXPECTED_GENES):
        raise ValueError(f"Combined shape is {adata.shape}, expected {(EXPECTED_CELLS, EXPECTED_GENES)}")
    if "tissue" not in adata.obs:
        raise KeyError("Combined obs has no tissue column")
    tissues = adata.obs["tissue"].astype(str).to_numpy()
    observed_counts = pd.Series(tissues).value_counts().to_dict()
    if observed_counts != EXPECTED_COUNTS:
        raise ValueError(f"Tissue counts differ from expected: {observed_counts}")

    for tissue, expected in EXPECTED_COUNTS.items():
        positions = np.flatnonzero(tissues == tissue)
        if len(positions) != expected:
            raise ValueError(f"{tissue}: found {len(positions)} cells, expected {expected}")
        destination = output_dir / f"m1-qc-cleaned-{tissue}.h5ad"
        temporary = output_dir / f"m1-qc-cleaned-{tissue}.tmp.h5ad"
        print(f"Writing {tissue}: {expected} cells x {adata.n_vars} genes", flush=True)
        part = adata[positions].to_memory()
        if part.shape != (expected, EXPECTED_GENES):
            raise ValueError(f"{tissue}: partition shape is {part.shape}")
        if not part.obs["tissue"].astype(str).eq(tissue).all():
            raise ValueError(f"{tissue}: partition contains another tissue label")
        part.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, destination)
        del part
        gc.collect()
    adata.file.close()


def verify(source: Path, output_dir: Path) -> None:
    paths = sorted(output_dir.glob("m1-qc-cleaned-*.h5ad"))
    paths = [path for path in paths if ".tmp." not in path.name]
    expected_paths = {output_dir / f"m1-qc-cleaned-{tissue}.h5ad" for tissue in EXPECTED_COUNTS}
    if set(paths) != expected_paths or len(paths) != 16:
        raise ValueError(f"Expected exactly 16 tissue files, found {len(paths)}: {paths}")

    combined = ad.read_h5ad(source, backed="r")
    combined_index = pd.Index(combined.obs_names.astype(str))
    combined_var = pd.Index(combined.var_names.astype(str))
    if combined.shape != (EXPECTED_CELLS, EXPECTED_GENES):
        raise ValueError(f"Combined shape is {combined.shape}")
    if not combined_index.is_unique:
        raise ValueError("Combined obs indices are not unique")

    split_indices: list[pd.Index] = []
    total_cells = 0
    for path in paths:
        tissue = path.stem.removeprefix("m1-qc-cleaned-")
        part = ad.read_h5ad(path, backed="r")
        if part.shape != (EXPECTED_COUNTS[tissue], EXPECTED_GENES):
            raise ValueError(f"{path}: shape {part.shape} does not match expected")
        if not part.obs["tissue"].astype(str).eq(tissue).all():
            raise ValueError(f"{path}: cells do not all carry tissue label {tissue}")
        part_var = pd.Index(part.var_names.astype(str))
        if not part_var.equals(combined_var):
            raise ValueError(f"{path}: gene list differs from combined file")
        part_index = pd.Index(part.obs_names.astype(str))
        if not part_index.is_unique:
            raise ValueError(f"{path}: obs indices are not unique")
        split_indices.append(part_index)
        total_cells += part.n_obs
        print(f"Verified {path.name}: {part.n_obs} cells, {part.n_vars} genes", flush=True)
        part.file.close()

    union = split_indices[0].append(split_indices[1:])
    if total_cells != EXPECTED_CELLS or len(union) != EXPECTED_CELLS:
        raise ValueError(f"Split cell total is {total_cells}, index total is {len(union)}")
    if not union.is_unique:
        raise ValueError("Obs indices are duplicated across split files")
    if set(union) != set(combined_index):
        missing = set(combined_index) - set(union)
        extra = set(union) - set(combined_index)
        raise ValueError(f"Obs index union differs: {len(missing)} missing, {len(extra)} extra")
    combined.file.close()
    print(
        f"Split verification passed: 16 files, {total_cells} cells, "
        f"{EXPECTED_GENES} genes per file, unique complete obs index union",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if not args.verify_only:
        split(args.source, args.output_dir)
    verify(args.source, args.output_dir)


if __name__ == "__main__":
    main()
