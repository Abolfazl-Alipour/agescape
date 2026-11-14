#!/usr/bin/env python3
"""Validate Perez filtered 10x counts and fit whole-nucleus donor models."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.followup_muscle import (  # noqa: E402
    MIN_DONORS, MIN_RESIDUAL_DF, add_bh, assert_count_gate, build_design, fit_ols,
)

DATA = REPO.parent / f"{REPO.name}-data" / "muscle" / "perez"
MANIFEST = REPO / "data" / "manifest" / "GSE167186.json"
RESULTS = REPO / "results"
PREFIX = "followup-muscle-perez"
EXPECTED_NUCLEI = 134486


def write_csv(frame: pd.DataFrame, name: str) -> None:
    path = RESULTS / name
    frame.to_csv(path, index=False)
    print(f"wrote {path}: {len(frame)} rows", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def donor_path(donor: dict) -> Path:
    return DATA / f"{donor['gsm']}_{donor['hm_id']}_filtered_feature_bc_matrix.h5"


def read_donor(donor: dict) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    path = donor_path(donor)
    if not path.is_file() or path.stat().st_size != donor["h5_size_bytes"]:
        raise ValueError(f"Missing or unexpected Perez file size: {path}")
    started = time.time()
    with h5py.File(path, "r") as handle:
        root_keys = sorted(handle.keys())
        if root_keys != ["matrix"]:
            raise ValueError(f"Unexpected non-matrix groups in Perez H5 {donor['gsm']}: {root_keys}")
        matrix_group = handle["matrix"]
        shape = tuple(int(x) for x in matrix_group["shape"][:])
        if shape[1] != donor["filtered_nuclei"]:
            raise ValueError(f"Released nuclei mismatch for {donor['gsm']}: {shape}")
        data = matrix_group["data"]
        if data.dtype != np.dtype("int32"):
            raise ValueError(f"Unexpected Perez count dtype for {donor['gsm']}: {data.dtype}")
        n_stored = len(data)
        count_sum = 0
        for start in range(0, len(data), 5_000_000):
            values = data[start:start + 5_000_000]
            assert_count_gate(values)
            count_sum += int(values.sum(dtype=np.int64))
        matrix = sp.csc_matrix(
            (data[:], matrix_group["indices"][:], matrix_group["indptr"][:]), shape=shape
        )
        barcodes = matrix_group["barcodes"][:]
        if len(barcodes) != shape[1]:
            raise ValueError(f"Barcode mismatch for {donor['gsm']}")
        gene_ids = matrix_group["features"]["id"][:].astype(str)
        gene_names = matrix_group["features"]["name"][:].astype(str)
        if len(gene_ids) != shape[0] or len(set(gene_ids)) != len(gene_ids):
            raise ValueError(f"Feature gate failed for {donor['gsm']}")
        gene_sums = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
        depths = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
        if gene_sums.sum() != count_sum or depths.sum() != count_sum:
            raise ValueError(f"Count-sum reconciliation failed for {donor['gsm']}")
        if (depths <= 0).any():
            raise ValueError(f"Zero-UMI released nucleus in {donor['gsm']}")
    summary = {
        "gsm": donor["gsm"], "donor_id": donor["hm_id"], "sample_title": donor["sample_title"],
        "group": donor["group"], "old": int(donor["group"] == "Old"), "age_years": int(donor["age"]),
        "n_nuclei": shape[1], "n_genes": shape[0], "stored_count_values": n_stored,
        "total_raw_umi": count_sum, "mean_umi": float(depths.mean()),
        "median_umi": float(np.median(depths)), "mean_log_umi": float(np.log(depths).mean()),
        "minimum_umi": int(depths.min()), "maximum_umi": int(depths.max()),
        "all_counts_finite_nonnegative_integer_valued": True,
        "h5_root_groups": ";".join(root_keys), "released_nucleus_annotations_present": False,
        "elapsed_seconds": time.time() - started,
    }
    return summary, gene_sums, gene_ids, gene_names


def continuous_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    design = np.column_stack([
        np.ones(len(frame)), frame["age_years"].to_numpy(float) / 10,
        frame["mean_log_umi"].to_numpy(float),
    ])
    if len(frame) < MIN_DONORS or np.linalg.matrix_rank(design) != 3:
        raise ValueError("Perez exact-age design is not full rank")
    if len(frame) - design.shape[1] < MIN_RESIDUAL_DF or not np.isfinite(design).all():
        raise ValueError("Perez exact-age design lacks residual degrees of freedom")
    return design, ["intercept", "age_per_10y", "mean_log_umi"]


def run() -> None:
    started = time.time()
    manifest = json.loads(MANIFEST.read_text())
    summaries, count_rows = [], []
    expected_gene_ids = expected_gene_names = None
    for donor in manifest["donors"]:
        summary, gene_sums, gene_ids, gene_names = read_donor(donor)
        if expected_gene_ids is None:
            expected_gene_ids, expected_gene_names = gene_ids, gene_names
        elif not np.array_equal(gene_ids, expected_gene_ids) or not np.array_equal(gene_names, expected_gene_names):
            raise ValueError(f"Perez feature order differs for {donor['gsm']}")
        summaries.append(summary)
        count_rows.append(gene_sums)
        print(f"Perez {donor['gsm']} {donor['hm_id']}: {summary['n_nuclei']} nuclei validated", flush=True)
    donors = pd.DataFrame(summaries).sort_values("gsm").reset_index(drop=True)
    if len(donors) != 17 or donors["n_nuclei"].sum() != EXPECTED_NUCLEI:
        raise ValueError("Perez donor or released-nuclei total gate failed")
    if donors.groupby("group")["gsm"].nunique().to_dict() != {"Old": 11, "Young": 6}:
        raise ValueError("Perez released group sizes changed")
    counts = np.vstack(count_rows).astype(float)
    totals = counts.sum(axis=1)
    expression = np.log1p(counts * (1e6 / totals[:, None]))
    gene_filter = (counts.sum(axis=0) >= 10) & ((counts > 0).sum(axis=0) >= 2)
    rows = []
    for model in ["primary_group_depth_adjusted", "exact_age_10y_sensitivity"]:
        if model == "primary_group_depth_adjusted":
            design, names, omitted = build_design(donors, depth=True, covariates=())
        else:
            design, names = continuous_design(donors)
            omitted = []
        fitted = fit_ols(expression[:, gene_filter], design)
        testable = np.isfinite(fitted["pvalue"]) & np.isfinite(fitted["se"]) & (fitted["se"] > 0)
        gene_indices = np.flatnonzero(gene_filter)[testable]
        for position, gene_index in enumerate(gene_indices):
            rows.append({
                "assay": "snRNA", "unit": "whole_nucleus", "model": model,
                "gene_id": expected_gene_ids[gene_index], "gene": expected_gene_names[gene_index],
                "effect": fitted["effect"][testable][position],
                "effect_definition": "old_minus_young" if model.startswith("primary") else "per_10_exact_years",
                "se": fitted["se"][testable][position], "ci95_low": fitted["ci95_low"][testable][position],
                "ci95_high": fitted["ci95_high"][testable][position],
                "pvalue": fitted["pvalue"][testable][position], "r2": fitted["r2"][testable][position],
                "residual_df": fitted["residual_df"][testable][position],
                "n_donors": len(donors), "n_young": int((donors["old"] == 0).sum()),
                "n_old": int((donors["old"] == 1).sum()), "n_nuclei": int(donors["n_nuclei"].sum()),
                "design_columns": ";".join(names), "omitted_covariates": ";".join(omitted),
                "effect_scale": "log1p_CPM_donor_summed_raw_counts",
                "composition_confounded": True,
            })
    genes = add_bh(pd.DataFrame(rows), ["model"])
    genes["multiplicity_family"] = "all_whole_nucleus_genes_within_model"
    write_csv(donors, f"{PREFIX}-donor-summary.csv")
    write_csv(genes, f"{PREFIX}-gene-models.csv")
    endpoints = pd.DataFrame([
        {"endpoint": "whole_nucleus_descriptive", "status": "completed", "outcomes_examined": True,
         "reason": "17 donor summaries from released filtered nuclei."},
        {"endpoint": "whole_nucleus_donor_pseudobulk_gene_models", "status": "completed_composition_confounded", "outcomes_examined": True,
         "reason": "Group and exact-age sensitivity models; cell-composition changes cannot be separated from cell-intrinsic expression."},
        {"endpoint": "P1_transport", "status": "blocked_composition_confounded", "outcomes_examined": False,
         "reason": "No released per-nucleus annotation; no whole-nucleus transport claim."},
        {"endpoint": "cell_type_models", "status": "blocked_no_released_annotations", "outcomes_examined": False,
         "reason": "No annotations were inferred or created."},
        {"endpoint": "donor_CLR_composition", "status": "blocked_no_released_annotations", "outcomes_examined": False,
         "reason": "No cell-type counts are available."},
        {"endpoint": "Hellinger_displacement", "status": "blocked_no_released_annotations", "outcomes_examined": False,
         "reason": "No defensible cell-type-specific young reference."},
        {"endpoint": "bulk_72_man_arm", "status": "excluded", "outcomes_examined": False,
         "reason": "Bulk RNA-seq is outside single-nucleus inference and was not downloaded."},
    ])
    write_csv(endpoints, f"{PREFIX}-endpoint-status.csv")
    verification = {
        "files": len(donors), "donors": len(donors), "young_donors": 6, "old_donors": 11,
        "released_filtered_nuclei": int(donors["n_nuclei"].sum()),
        "expected_filtered_nuclei": EXPECTED_NUCLEI,
        "stored_count_values": int(donors["stored_count_values"].sum()),
        "total_raw_umi": int(donors["total_raw_umi"].sum()),
        "all_counts_finite_nonnegative_integer_valued": True,
        "all_files_int32": True, "feature_order_identical_across_donors": True,
        "genes_after_filter": int(gene_filter.sum()), "gene_models_tested": len(genes),
        "primary_group_model": "old + mean_log_umi", "exact_age_sensitivity": "age_per_10y + mean_log_umi",
        "whole_nucleus_models_composition_confounded": True,
        "bulk_arm_downloaded_or_analyzed": False, "raw_fastq_downloaded": False,
        "elapsed_seconds": time.time() - started,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "memory_limit_gb": 6, "minimum_free_gb": 30,
    }
    (RESULTS / f"{PREFIX}-verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")


def finalize() -> None:
    manifest = json.loads(MANIFEST.read_text())
    verification = json.loads((RESULTS / f"{PREFIX}-verification.json").read_text())
    for donor in manifest["donors"]:
        local = donor_path(donor)
        donor.update({
            "h5_downloaded": True, "h5_local_path": str(Path("..") / local.relative_to(REPO.parent)),
            "h5_sha256": sha256(local), "h5_count_values_validated": True,
        })
    manifest["scope"].update({"metadata_and_access_only": False, "expression_matrices_downloaded": True,
                              "h5_expression_values_read": True, "biological_outcomes_examined": True})
    manifest["snrna_arm"]["h5_schema"]["count_values_read"] = True
    manifest["snrna_arm"]["h5_schema"]["all_count_values_finite_nonnegative_integer_valued"] = True
    manifest["snrna_arm"]["released_per_nucleus_annotations_available"] = False
    manifest["analysis"] = {
        "completed": "2026-11-14", "frozen_plan": "data/manifest/muscle_analysis_plan.json",
        "script": "scripts/run_followup_perez.py", "report": "reports/followup-muscle-analysis.md",
        "inferential_unit": "donor", "primary_contrast": "released Old versus Young plus mean_log_umi",
        "exact_age_sensitivity": "age per 10 years plus mean_log_umi",
        "scope": "whole-nucleus descriptive and donor-pseudobulk gene models only",
        "composition_confounded": True, "cell_type_CLR_Hellinger_P1": "blocked_no_released_annotations",
        "annotation_search": {
            "completed": "2026-11-14", "result": "no_released_barcode_level_annotations_found",
            "sources": ["GEO series and per-donor supplementary inventories",
                        "all 17 downloaded filtered H5 root schemas (matrix only)",
                        "primary full text and supplementary inventory",
                        "public processed-object and repository searches"],
            "secondary_reannotations_not_substituted": True,
        },
        "bulk_72_man_arm": "excluded_not_downloaded", "raw_fastq_downloaded": False,
        "peak_rss_bytes_observed": verification["peak_rss_bytes"],
    }
    manifest["verdict"].update({"expression_analysis_performed": True,
                                "snrna_raw_counts": "validated_all_17_int32_h5_values_and_134486_nuclei",
                                "exact_age_analysis": "completed_predeclared_sensitivity",
                                "cell_type_analysis": "blocked_no_released_per_nucleus_annotations"})
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["run", "finalize"])
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    if args.command == "run":
        run()
    else:
        finalize()


if __name__ == "__main__":
    main()
