#!/usr/bin/env python3
"""Validate raw Cui counts and compute donor-level descriptive QC without outcomes."""

from __future__ import annotations

import argparse
import gzip
import gc
import hashlib
import json
from pathlib import Path
import resource
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.followup_cui_testis import summarize_metric, validate_raw_counts

DATA = REPO.parent / "agescape-revision-data" / "cui-testis"
RESULTS = REPO / "results"
PREFIX = "followup-cui-testis"
CAVEAT = (
    "Changing germ-cell-stage composition can mimic broad expression displacement; "
    "testis hits/nulls cannot support general whole-body aging claims."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def donor_qc(donor: dict, entries: list[dict]) -> tuple[dict, pd.DataFrame, list[str]]:
    paths = {entry["kind"]: REPO / entry["local_path"] for entry in entries}
    for kind, path in paths.items():
        expected = next(entry for entry in entries if entry["kind"] == kind)
        if not path.exists() or path.stat().st_size != expected["size_bytes"]:
            raise ValueError(f"Missing or wrong-sized {kind} file for {donor['donor']}")
        if sha256(path) != expected["sha256"]:
            raise ValueError(f"Checksum mismatch for {path}")
    matrix = sp.csr_matrix(scipy.io.mmread(str(paths["matrix"])))
    barcodes = pd.read_csv(paths["barcodes"], sep="\t", header=None, dtype=str)[0].to_numpy()
    features = pd.read_csv(paths["features"], sep="\t", header=None, dtype=str)
    if features.shape[1] < 2:
        raise ValueError(f"Invalid features for {donor['donor']}")
    features = features.iloc[:, :2]
    features.columns = ["feature_id", "symbol"]
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(f"Matrix dimensions do not match metadata for {donor['donor']}")
    if not pd.Index(barcodes).is_unique or not features["feature_id"].is_unique:
        raise ValueError(f"Duplicate barcode or feature ID for {donor['donor']}")
    validation = validate_raw_counts(matrix)
    umi = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(float)
    mt_mask = features["symbol"].str.startswith("MT-").to_numpy()
    if not mt_mask.any():
        raise ValueError(f"No mitochondrial genes for {donor['donor']}")
    mt_umi = np.asarray(matrix[mt_mask, :].sum(axis=0)).ravel().astype(float)
    mt_fraction = np.divide(mt_umi, umi, out=np.full_like(umi, np.nan), where=umi > 0)
    if (umi <= 0).any() or not np.isfinite(mt_fraction).all():
        raise ValueError(f"Zero-library or invalid mitochondrial fraction for {donor['donor']}")
    qc_pass = (
        (detected >= 200) & (detected <= 5500) & (umi >= 800) & (umi <= 17000)
        & (mt_fraction <= 0.25)
    )
    row = {
        "donor": donor["donor"], "gsm": donor["gsm"], "library": donor["library"],
        "primary_age": donor["primary_age"], "primary_decade": donor["primary_decade"],
        "geo_age": donor["geo_age"], "geo_decade": donor["geo_decade"], "bmi": donor["bmi"],
        "genes": int(matrix.shape[0]), "raw_barcodes": int(matrix.shape[1]),
        **validation, **summarize_metric(umi, "umi"),
        **summarize_metric(detected, "detected_features"),
        **summarize_metric(mt_fraction, "mitochondrial_fraction"),
        "paper_qc_pass_cells": int(qc_pass.sum()),
        "paper_qc_fail_cells": int((~qc_pass).sum()),
        "paper_qc_pass_fraction": float(qc_pass.mean()),
        "mitochondrial_genes": int(mt_mask.sum()),
        "all_stored_counts_finite_nonnegative_integer": True,
        "annotation_crosswalk_available": False,
        "analysis": "descriptive_raw_count_qc_only_no_age_outcome_model",
        "caveat": CAVEAT,
    }
    del matrix, umi, detected, mt_umi, mt_fraction, qc_pass
    gc.collect()
    return row, features, barcodes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=DATA)
    args = parser.parse_args()
    manifest_path = REPO / "data" / "manifest" / "GSE254315.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["annotation_gate"]["public_barcode_to_cell_type_crosswalk"]:
        raise ValueError("This QC-only script must not be reused as an unblocked analysis")
    checkpoint = args.data_root / "interim" / "qc"
    checkpoint.mkdir(parents=True, exist_ok=True)
    by_donor = {}
    for entry in manifest["files"]:
        by_donor.setdefault(entry["donor"], []).append(entry)
    rows = []
    feature_reference = None
    feature_hashes = set()
    for entry in manifest["files"]:
        if entry["kind"] != "features":
            continue
        with gzip.open(REPO / entry["local_path"], "rb") as handle:
            feature_hashes.add(hashlib.sha256(handle.read()).hexdigest())
    if len(feature_hashes) != 1:
        raise ValueError("Decompressed feature content differs across donors")
    feature_hash = feature_hashes.pop()
    donor_validations = []
    for donor in manifest["donors"]:
        cache = checkpoint / f"{donor['donor']}.json"
        expected_hashes = {entry["kind"]: entry["sha256"] for entry in by_donor[donor["donor"]]}
        if cache.exists():
            saved = json.loads(cache.read_text())
            if saved["source_sha256"] == expected_hashes:
                rows.append(saved["qc"])
                donor_validations.append({"donor": donor["donor"], "checkpoint_verified": True})
                continue
        print(f"Validating {donor['donor']} ({len(rows) + 1}/25)", flush=True)
        row, features, barcodes = donor_qc(donor, by_donor[donor["donor"]])
        if feature_reference is None:
            feature_reference = features.copy()
        elif not features.equals(feature_reference):
            raise ValueError(f"Feature reference differs for {donor['donor']}")
        rows.append(row)
        donor_validations.append({
            "donor": donor["donor"], "barcodes_unique": True, "feature_ids_unique": True,
            "matrix_dimensions_match": True, "feature_reference_sha256": feature_hash,
        })
        cache.write_text(json.dumps({"source_sha256": expected_hashes, "qc": row}, indent=2) + "\n")
        del features, barcodes
        gc.collect()
    if len(rows) != 25:
        raise ValueError(f"Expected 25 QC rows, found {len(rows)}")
    output = RESULTS / f"{PREFIX}-qc.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    total_barcodes = sum(row["raw_barcodes"] for row in rows)
    total_qc_pass = sum(row["paper_qc_pass_cells"] for row in rows)
    validation = {
        "accession": "GSE254315", "donor_count": len(rows),
        "raw_barcodes_total": total_barcodes,
        "geo_summary_cell_count": 99130,
        "raw_barcode_count_matches_geo_summary": total_barcodes == 99130,
        "paper_qc_pass_cells": total_qc_pass,
        "geo_summary_minus_paper_qc_pass": 99130 - total_qc_pass,
        "published_cell_count_reconciled": total_qc_pass == 99130,
        "feature_reference_identical_decompressed_content": True,
        "feature_reference_sha256": feature_hash,
        "donors": donor_validations,
        "public_barcode_to_cell_type_crosswalk": False,
        "expression_outcomes_computed": False,
        "age_associations_computed": False,
        "cell_level_pvalues_computed": False,
        "caveat": CAVEAT,
    }
    (RESULTS / f"{PREFIX}-validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest["descriptive_qc"] = {
        "raw_barcodes": total_barcodes, "paper_qc_pass_cells": total_qc_pass,
        "geo_summary_cells": 99130, "geo_summary_minus_paper_qc_pass": 99130 - total_qc_pass,
        "published_cell_count_reconciled": total_qc_pass == 99130,
        "decompressed_feature_content_identical": True,
        "feature_content_sha256": feature_hash,
        "expression_outcomes_computed": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    resource_path = RESULTS / f"{PREFIX}-resource-use.json"
    previous_rss = json.loads(resource_path.read_text())["max_rss_kib"] if resource_path.exists() else 0
    resources = {
        "guard_memory_gb": 8, "guard_min_free_gb": 30, "parallel_heavy_jobs": 1,
        "maximum_single_file_gb": 0.2,
        "max_rss_kib": max(previous_rss, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "sequential_donor_processing": True,
        "external_checkpoints": str(Path("..") / checkpoint.relative_to(REPO.parent)),
    }
    (RESULTS / f"{PREFIX}-resource-use.json").write_text(json.dumps(resources, indent=2) + "\n")
    print(f"Wrote descriptive QC for {len(rows)} donors and {total_barcodes} raw barcodes", flush=True)


if __name__ == "__main__":
    main()
