#!/usr/bin/env python3
"""Parallel Yazar mature A2 expression variance.

Reproduces scripts/run_mature_between_subject.py::run_expression_a2 for
Yazar exactly, one (contrast, group) unit per worker, then assembles the
joint stratified-global rows and the BH families in the parent with the
same code as the serial tail. Every per-unit computation (seeds, perm
orders, covariate pruning, detected-gene filter) is copied verbatim from
run_expression_a2; only the loop nesting changes (units run in parallel
instead of serially). The serial run was killed after ~95 minutes inside
A2 with nothing written: a single Yazar-scale batch-arm call benchmarks
at ~400 s and depth-only at ~55 s, so the serial A2 is ~4 to 6 hours.
Beast has 12 vCPUs and the serial run used one.
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
from scipy import sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_between_subject_gao_yazar import (  # noqa: E402
    A2_N_PERM,
    MASTER_SEED,
    MIN_CELLS_PER_DONOR,
    MIN_DONORS_A2_BIN,
    YAZAR_CACHE,
    YAZAR_H5AD,
    _prune_covariate_matrix,
    _residualize_against_depth,
    _summarize_a2_perm_stats,
    build_human_gene_projection,
    compute_a2_variance_metrics,
    reconcile_yazar_metadata,
)
from scripts.run_mature_between_subject import (  # noqa: E402
    CONTRASTS,
    RESULTS,
    prepare_donors,
    recompute_rarefied,
)

DATASET = "yazar"

# Set by the parent before forking; inherited copy-on-write by workers.
META: dict[str, Any] = {}
DONORS: pd.DataFrame | None = None
A2_INPUTS: dict[str, Any] = {}
MANIFEST: dict[str, str] = {}
CACHE_PROJECTION = None


def _load_counts(group: str, selected_donors: list[str]) -> np.ndarray:
    pooled_name = META["pooled_name"]
    donor_ids = META["donor_ids"]
    donor_to_index = {donor: index for index, donor in enumerate(donor_ids)}
    if group == pooled_name:
        return A2_INPUTS["a2_pool"][[donor_to_index[donor] for donor in selected_donors]]
    with (YAZAR_CACHE / MANIFEST[group]).open("rb") as handle:
        raw = pickle.load(handle)
    raw.columns = raw.columns.astype(str)
    dense = raw[selected_donors].to_numpy(dtype=float).T
    return (sp.csr_matrix(dense) @ CACHE_PROJECTION).toarray()


def run_a2_unit(payload: tuple[int, str]) -> dict[str, Any]:
    """One (contrast, group) unit, copied verbatim from run_expression_a2."""
    contrast_index, group_name = payload
    contrast, reference_bin, comparison_bin, reference_age, comparison_age = CONTRASTS[contrast_index]
    dataset = DATASET
    donors = DONORS
    cell_types = META["cell_types"]
    pooled_name = META["pooled_name"]
    rows: list[dict[str, Any]] = []
    ineligible: list[dict[str, Any]] = []
    stratified: dict[tuple[str, float], list[dict[str, Any]]] = {}
    unique_donors = sorted(donors.loc[donors["bin"].isin([reference_bin, comparison_bin]), "donor_id"].unique())
    donor_positions = {donor: index for index, donor in enumerate(unique_donors)}
    rng = np.random.default_rng(MASTER_SEED + 500 + contrast_index)
    global_ranks = np.array([rng.permutation(len(unique_donors)) for _ in range(A2_N_PERM)], dtype=np.int32)
    frame = donors[
        (donors["group"] == group_name)
        & donors["bin"].isin([reference_bin, comparison_bin])
        & (donors["n_cells_all"] >= MIN_CELLS_PER_DONOR)
    ].copy()
    if group_name != pooled_name and group_name not in MANIFEST:
        ineligible.append({"dataset": dataset, "contrast": contrast, "group": group_name, "reason": "not_in_yazar_cache"})
        return {"rows": rows, "ineligible": ineligible, "stratified": stratified}
    n_reference = int((frame["bin"] == reference_bin).sum())
    n_comparison = int((frame["bin"] == comparison_bin).sum())
    if min(n_reference, n_comparison) < MIN_DONORS_A2_BIN:
        ineligible.append({
            "dataset": dataset, "contrast": contrast, "group": group_name,
            "n_reference": n_reference, "n_comparison": n_comparison, "reason": "lt6_in_reference_or_comparison",
        })
        return {"rows": rows, "ineligible": ineligible, "stratified": stratified}
    order = frame["donor_id"].tolist()
    raw_counts = _load_counts(group_name, order)
    detected = (raw_counts > 0).sum(axis=0) >= 0.5 * len(order)
    detected_counts = raw_counts[:, detected]
    totals = detected_counts.sum(axis=1, keepdims=True)
    expression = np.log1p(detected_counts / np.maximum(totals, 1.0) * 1e4)
    is_comparison = (frame["bin"] == comparison_bin).to_numpy(dtype=bool)
    depth_expression = _residualize_against_depth(expression, frame["mean_log_umi_all"].to_numpy(float))
    frame["_is_comp"] = is_comparison.astype(float)
    batch_frame, batch_columns, _, batch_status = _prune_covariate_matrix(
        frame, ["sex_male", "mean_log_umi_all"], "study_or_pool", 2,
        10, None,
        require_is_comp_collinear_check="_is_comp", atomic_batch_block=True,
    )
    pool_status = "evaluated"
    n_excluded = 0
    if batch_status == "evaluated":
        keep = [order.index(donor) for donor in batch_frame["donor_id"].tolist()]
        batch_comparison = is_comparison[keep]
        if min((~batch_comparison).sum(), batch_comparison.sum()) < MIN_DONORS_A2_BIN:
            batch_status = "insufficient_after_singleton_drop"
        else:
            n_excluded = len(order) - len(keep)
            batch_raw = raw_counts[keep]
            batch_detected = (batch_raw > 0).sum(axis=0) >= 0.5 * len(keep)
            batch_counts = batch_raw[:, batch_detected]
            batch_totals = batch_counts.sum(axis=1, keepdims=True)
            batch_expression = np.log1p(batch_counts / np.maximum(batch_totals, 1.0) * 1e4)
            z_batch = np.column_stack([np.ones(len(keep)), batch_frame[batch_columns].to_numpy(float)])
            q_check, _ = np.linalg.qr(np.column_stack([z_batch, batch_comparison.astype(float)]), mode="reduced")
            if np.any(1 - np.sum(q_check ** 2, axis=1) <= 1e-4):
                batch_status = "extreme_leverage_guard"
    if batch_status != "evaluated":
        pool_status = "not_evaluable_within_pool_df_slice"
        n_excluded = 0
        batch_frame, batch_columns, _, _ = _prune_covariate_matrix(
            frame, ["sex_male", "mean_log_umi_all"], None, 1, 5, None,
            require_is_comp_collinear_check="_is_comp", atomic_batch_block=True,
        )
        keep = [order.index(donor) for donor in batch_frame["donor_id"].tolist()]
        batch_comparison = is_comparison
        batch_expression = expression
        z_batch = np.column_stack([np.ones(len(order)), batch_frame[batch_columns].to_numpy(float)])
    batch_order = [order[index] for index in keep]
    depth_orders = np.argsort(global_ranks[:, [donor_positions[donor] for donor in order]], axis=1)
    batch_orders = depth_orders if batch_order == order else np.argsort(global_ranks[:, [donor_positions[donor] for donor in batch_order]], axis=1)
    is_pool = group_name == pooled_name
    for epsilon in (1e-8, 1e-3):
        batch_stats, batch_permutations = compute_a2_variance_metrics(
            batch_expression, batch_comparison, z_batch, n_perm=A2_N_PERM,
            seed=MASTER_SEED + (0 if epsilon < 1e-5 else 1), eps=epsilon,
            custom_perm_orders=None if is_pool else batch_orders, return_perm=True,
        )
        depth_stats, depth_permutations = compute_a2_variance_metrics(
            depth_expression, is_comparison, None, n_perm=A2_N_PERM,
            seed=MASTER_SEED + (2 if epsilon < 1e-5 else 3), eps=epsilon,
            custom_perm_orders=None if is_pool else depth_orders, return_perm=True,
        )
        disagreement = batch_stats["mean_log_var_ratio"] * depth_stats["mean_log_var_ratio"] < 0
        common = {
            "dataset": dataset, "contrast": contrast,
            "reference_group": reference_bin, "comparison_group": comparison_bin,
            "reference_age_range": reference_age, "comparison_age_range": comparison_age,
            "eps_used": epsilon, "group": group_name, "is_pooled": is_pool,
            "is_stratified_global": False, "arm_direction_disagreement_flag": bool(disagreement),
        }
        batch_gene_count = batch_stats.pop("n_genes")
        depth_gene_count = depth_stats.pop("n_genes")
        row_batch = {
            **common, "arm": "batch_residualized", "pool_fe_status": pool_status,
            "n_donors_ref_young": int((~batch_comparison).sum()), "n_donors_comp": int(batch_comparison.sum()),
            "n_donors_excluded_singleton_batch": n_excluded,
            "n_cells_total": int(frame.iloc[keep]["n_cells_all"].sum()),
            **batch_stats, "gene_count": batch_gene_count,
        }
        row_depth = {
            **common, "arm": "depth_only_residualized", "pool_fe_status": "depth_only_all_donors",
            "n_donors_ref_young": n_reference, "n_donors_comp": n_comparison,
            "n_donors_excluded_singleton_batch": 0, "n_cells_total": int(frame["n_cells_all"].sum()),
            **depth_stats, "gene_count": depth_gene_count,
        }
        for row in (row_batch, row_depth):
            row.update({
                "q_bh_eta2": np.nan, "q_bh_lvr_one_sided": np.nan,
                "q_bh_lvr_two_sided": np.nan, "q_bh_tot_lvr_one_sided": np.nan,
                "family_size_m": np.nan,
            })
        rows.extend([row_batch, row_depth])
        if not is_pool:
            if pool_status == "evaluated":
                stratified.setdefault(("batch_residualized", epsilon), []).append({
                    "row": row_batch, "perm": batch_permutations,
                    "reference": {batch_order[i] for i, value in enumerate(batch_comparison) if not value},
                    "comparison": {batch_order[i] for i, value in enumerate(batch_comparison) if value},
                })
            stratified.setdefault(("depth_only_residualized", epsilon), []).append({
                "row": row_depth, "perm": depth_permutations,
                "reference": {order[i] for i, value in enumerate(is_comparison) if not value},
                "comparison": {order[i] for i, value in enumerate(is_comparison) if value},
            })
    return {"rows": rows, "ineligible": ineligible, "stratified": stratified}


def assemble(units: dict[tuple[int, str], dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Joint stratified rows and BH, copied verbatim from run_expression_a2 tail."""
    rows: list[dict[str, Any]] = []
    ineligible: list[dict[str, Any]] = []
    for contrast_index, (contrast, reference_bin, comparison_bin, reference_age, comparison_age) in enumerate(CONTRASTS):
        stratified: dict[tuple[str, float], list[dict[str, Any]]] = {}
        for group_name in [*META["cell_types"], META["pooled_name"]]:
            out = units[(contrast_index, group_name)]
            rows.extend(out["rows"])
            ineligible.extend(out["ineligible"])
            for key, items in out["stratified"].items():
                stratified.setdefault(key, []).extend(items)
        for epsilon in (1e-8, 1e-3):
            joint_rows = {}
            for arm in ("batch_residualized", "depth_only_residualized"):
                items = stratified.get((arm, epsilon), [])
                if len(items) < 2:
                    continue
                joint = _summarize_a2_perm_stats(
                    n_g=int(round(np.mean([item["row"]["gene_count"] for item in items]))),
                    resid_df_model=int(round(np.mean([item["row"]["residual_df_resid_model"] for item in items]))),
                    obs_mean_eta2=float(np.mean([item["row"]["mean_eta2"] for item in items])),
                    obs_med_eta2=float(np.median([item["row"]["mean_eta2"] for item in items])),
                    obs_mean_f=float(np.mean([item["row"]["mean_F"] for item in items])),
                    obs_mean_lvr=float(np.mean([item["row"]["mean_log_var_ratio"] for item in items])),
                    obs_med_lvr=float(np.median([item["row"]["mean_log_var_ratio"] for item in items])),
                    obs_tot_lvr=float(np.mean([item["row"]["log_total_var_ratio"] for item in items])),
                    obs_frac_gt=float(np.mean([item["row"]["frac_genes_comp_gt_young"] for item in items])),
                    perm_eta2=np.mean([item["perm"]["perm_eta2"] for item in items], axis=0),
                    perm_lvr=np.mean([item["perm"]["perm_lvr"] for item in items], axis=0),
                    perm_tot_lvr=np.mean([item["perm"]["perm_tot_lvr"] for item in items], axis=0),
                )
                gene_count = joint.pop("n_genes")
                row = {
                    "dataset": DATASET, "contrast": contrast, "reference_group": reference_bin,
                    "comparison_group": comparison_bin, "reference_age_range": reference_age,
                    "comparison_age_range": comparison_age, "arm": arm, "eps_used": epsilon,
                    "group": "Stratified_CellType_Global", "is_pooled": True, "is_stratified_global": True,
                    "pool_fe_status": f"stratified_{len(items)}_cell_types",
                    "n_donors_ref_young": len(set().union(*(item["reference"] for item in items))),
                    "n_donors_comp": len(set().union(*(item["comparison"] for item in items))),
                    "n_donors_excluded_singleton_batch": 0,
                    "n_cells_total": sum(item["row"]["n_cells_total"] for item in items),
                    **joint, "gene_count": gene_count, "q_bh_eta2": np.nan,
                    "q_bh_lvr_one_sided": np.nan, "q_bh_lvr_two_sided": np.nan,
                    "q_bh_tot_lvr_one_sided": np.nan, "family_size_m": np.nan,
                    "arm_direction_disagreement_flag": False,
                }
                joint_rows[arm] = row
                rows.append(row)
            if len(joint_rows) == 2:
                disagreement = joint_rows["batch_residualized"]["mean_log_var_ratio"] * joint_rows["depth_only_residualized"]["mean_log_var_ratio"] < 0
                for row in joint_rows.values():
                    row["arm_direction_disagreement_flag"] = bool(disagreement)
    result = pd.DataFrame(rows)
    non_stratified = (~result["is_pooled"]) & (~result["pool_fe_status"].str.startswith("not_evaluable"))
    p_q = (
        ("p_eta2_perm", "q_bh_eta2"),
        ("p_lvr_one_sided_comp_gt_young", "q_bh_lvr_one_sided"),
        ("p_lvr_two_sided", "q_bh_lvr_two_sided"),
        ("p_tot_lvr_one_sided", "q_bh_tot_lvr_one_sided"),
    )
    for _, indices in result[non_stratified].groupby(["contrast", "arm", "eps_used"]).groups.items():
        result.loc[indices, "family_size_m"] = len(indices)
        for p_column, q_column in p_q:
            result.loc[indices, q_column] = multipletests(result.loc[indices, p_column], method="fdr_bh")[1]
    pooled = result["is_pooled"] & ~result["is_stratified_global"]
    result.loc[pooled, "family_size_m"] = 1
    for p_column, q_column in p_q:
        result.loc[pooled, q_column] = result.loc[pooled, p_column]
    for _, indices in result[non_stratified].groupby(["contrast", "arm", "eps_used"]).groups.items():
        for p_column, q_column in p_q:
            expected = multipletests(result.loc[indices, p_column], method="fdr_bh")[1]
            if not np.allclose(result.loc[indices, q_column], expected, rtol=1e-12, atol=1e-15):
                raise ValueError("A2 BH verification failed")
    return result, pd.DataFrame(ineligible)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--units", nargs="*", default=None,
                        help="Optional subset of units as contrast_index:group, for testing.")
    parser.add_argument("--out-prefix", default="yazar_a2_parallel",
                        help="Prefix for intermediate pickle outputs under results/.")
    args = parser.parse_args()
    global META, DONORS, A2_INPUTS, MANIFEST, CACHE_PROJECTION
    t0 = time.monotonic()
    META = reconcile_yazar_metadata(YAZAR_H5AD, YAZAR_CACHE)
    DONORS = pd.read_csv(RESULTS / "followup-mature-between-subject-yazar-donors.csv")
    projection_info = build_human_gene_projection(META["gene_symbols"])
    original = pd.read_csv(RESULTS / "followup-between-subject-yazar-donors.csv")
    original = original[original["age_years"] >= 30].copy()
    rarefied, reference_audit, A2_INPUTS = recompute_rarefied(META, original, projection_info)
    # The pipeline writes prepare_donors(meta, rarefied), not rarefied itself.
    donors_check = prepare_donors(META, rarefied)
    # Cross-check the deterministic recompute against the files the serial run wrote.
    prev_donors = pd.read_csv(RESULTS / "followup-mature-between-subject-yazar-donors.csv")
    assert list(donors_check.columns) == list(prev_donors.columns) and len(donors_check) == len(prev_donors)
    for col in donors_check.columns:
        a, b = donors_check[col], prev_donors[col]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            if not np.allclose(a.to_numpy(dtype=float), b.to_numpy(dtype=float), rtol=1e-12, atol=1e-12, equal_nan=True):
                raise ValueError(f"prepare_donors differs from serial run files in column {col}")
        elif not a.astype(str).equals(b.astype(str)):
            raise ValueError(f"prepare_donors differs from serial run files in column {col}")
    prev_audit = pd.read_csv(RESULTS / "followup-mature-between-subject-yazar-reference-audit.csv")
    assert list(reference_audit.columns) == list(prev_audit.columns) and len(reference_audit) == len(prev_audit)
    for col in reference_audit.columns:
        a, b = reference_audit[col], prev_audit[col]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            if not np.allclose(a.to_numpy(dtype=float), b.to_numpy(dtype=float), rtol=1e-12, atol=1e-12, equal_nan=True):
                raise ValueError(f"recompute_rarefied reference audit differs in column {col}")
        elif not a.astype(str).equals(b.astype(str)):
            raise ValueError(f"recompute_rarefied reference audit differs in column {col}")
    print(f"setup done in {time.monotonic() - t0:.1f}s; a2_pool shape {A2_INPUTS['a2_pool'].shape}", flush=True)
    MANIFEST = json.loads((YAZAR_CACHE / "pb_manifest.json").read_text())
    sample = next(iter(MANIFEST.values()))
    with (YAZAR_CACHE / sample).open("rb") as handle:
        sample_frame = pickle.load(handle)
    cache_info = build_human_gene_projection(sample_frame.index.astype(str))
    if cache_info["audit"]["retained_symbols"] != A2_INPUTS["retained_symbols"]:
        raise ValueError("Yazar cache genes differ from h5ad genes")
    CACHE_PROJECTION = cache_info["proj"]
    del sample_frame
    groups = [*META["cell_types"], META["pooled_name"]]
    print(f"groups: {len(groups)}; contrasts: {len(CONTRASTS)}; units: {len(groups) * len(CONTRASTS)}", flush=True)
    if args.units:
        wanted = set()
        for spec in args.units:
            ci, _, grp = spec.partition(":")
            wanted.add((int(ci), grp))
        unit_list = sorted(wanted)
    else:
        unit_list = [(ci, g) for ci in range(len(CONTRASTS)) for g in groups]
    results: dict[tuple[int, str], dict[str, Any]] = {}
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx) as pool:
        futures = {pool.submit(run_a2_unit, u): u for u in unit_list}
        done = 0
        for fut in as_completed(futures):
            unit = futures[fut]
            results[unit] = fut.result()
            done += 1
            print(f"[{done}/{len(unit_list)}] {unit} rows={len(results[unit]['rows'])} "
                  f"ineligible={len(results[unit]['ineligible'])} elapsed={time.monotonic() - t0:.0f}s", flush=True)
    out_path = RESULTS / f"followup-mature-between-subject-{args.out_prefix}.pkl"
    with out_path.open("wb") as handle:
        pickle.dump(results, handle, protocol=4)
    print(f"wrote {out_path} ({len(results)} units) in {time.monotonic() - t0:.0f}s", flush=True)
    if not args.units:
        variance, ineligible = assemble(results)
        variance.to_csv(RESULTS / "followup-mature-between-subject-variance-yazar.csv", index=False)
        ineligible.to_csv(RESULTS / "followup-mature-between-subject-variance-yazar-ineligible.csv", index=False)
        print(f"wrote variance {len(variance)} rows, ineligible {len(ineligible)} rows", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
