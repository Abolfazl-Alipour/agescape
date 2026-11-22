#!/usr/bin/env python3
"""Rerun Gao and Yazar between-subject Parts A1, A2, and A3 at age 30-plus."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import pickle
import resource
import sys
import time
from typing import Any

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import h5py
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import reference_without_donor
from scripts.run_followup_de_man_entropy import multinomial_downsample
from scripts.run_followup_between_subject_gao_yazar import (
    A2_N_PERM,
    GAO_H5AD,
    MASTER_SEED,
    MIN_CELLS_PER_DONOR,
    MIN_DONORS_A2_BIN,
    MIN_LOO_REF_CELLS,
    YAZAR_CACHE,
    YAZAR_H5AD,
    _fit_stratified_global_row,
    _prune_covariate_matrix,
    _residualize_against_depth,
    _summarize_a2_perm_stats,
    apply_slice_bh,
    build_human_gene_projection,
    compute_a2_variance_metrics,
    fit_contrast_adj_diff,
    fit_ols_model_row,
    read_csr_chunk,
    reconcile_gao_metadata,
    reconcile_yazar_metadata,
)
from scripts.run_h3_expression_vector import hellinger_distance
from scripts.run_mature_nov17_donor_spread import build_reference_stats
from src.followup_muscle import count_gate

RESULTS = REPO / "results"
RUNTIME = RESULTS / "followup-mature-between-subject-runtime.csv"
ITEM8_DISTANCES = RESULTS / "followup-mature-nov17-donor-spread-distances.csv"
CHUNK_ROWS = 8000
MIN_REFERENCE_DONORS = 2
N_BOOT = 2000
BOOT_SEED = 7
DEPTHS = {"1500": 1500, "500": 500, "3000": 3000}
CONTRASTS = (
    ("OLD_vs_MATURE_30_39", "YOUNG", "OLD", "30-39", "60-plus"),
    ("MID_vs_MATURE_30_39", "YOUNG", "MID", "30-39", "40-59"),
    ("OLD_vs_MID", "MID", "OLD", "40-59", "60-plus"),
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def append_runtime(pipeline: str, start_utc: str, started: float, status: str, details: str = "") -> None:
    row = pd.DataFrame([{
        "pipeline": pipeline,
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "duration_seconds": time.monotonic() - started,
        "model_used": "Gemini 4 Neon XHigh",
        "status": status,
        "max_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        "details": details,
    }])
    if RUNTIME.exists():
        previous = pd.read_csv(RUNTIME)
        previous = previous[previous["pipeline"] != pipeline]
        row = pd.concat([previous, row], ignore_index=True)
    row.to_csv(RUNTIME, index=False)


def matrix_group(handle: h5py.File, slot: str) -> h5py.Group:
    if slot == "layers/counts":
        return handle["layers"]["counts"]
    if slot == "raw/X":
        return handle["raw"]["X"]
    return handle["X"]


def accumulate_pairs(target: np.ndarray, matrix: sp.csr_matrix, donors: np.ndarray,
                     cell_types: np.ndarray, n_groups: int, pooled: bool,
                     include_pooled: bool = True) -> None:
    if pooled:
        pairs = donors
        source = matrix
    elif include_pooled:
        pairs = np.concatenate([donors * n_groups + cell_types, donors * n_groups + n_groups - 1])
        source = sp.vstack([matrix, matrix], format="csr")
    else:
        pairs = donors * n_groups + cell_types
        source = matrix
    unique, inverse = np.unique(pairs, return_inverse=True)
    indicator = sp.csr_matrix(
        (np.ones(len(inverse)), (inverse, np.arange(len(inverse)))),
        shape=(len(unique), len(inverse)),
    )
    sums = (indicator @ source).tocsr()
    for index, pair in enumerate(unique):
        begin, end = int(sums.indptr[index]), int(sums.indptr[index + 1])
        target[int(pair), sums.indices[begin:end]] += sums.data[begin:end]


def accumulate_reference(target: np.ndarray, matrix: sp.csr_matrix, donors: np.ndarray,
                         cell_types: np.ndarray, ref_map: dict[int, int], n_groups: int) -> None:
    pairs = np.concatenate([donors * n_groups + cell_types, donors * n_groups + n_groups - 1])
    source = sp.vstack([matrix, matrix], format="csr")
    unique, inverse = np.unique(pairs, return_inverse=True)
    indicator = sp.csr_matrix(
        (np.ones(len(inverse)), (inverse, np.arange(len(inverse)))),
        shape=(len(unique), len(inverse)),
    )
    sums = (indicator @ source).tocsr()
    for index, pair in enumerate(unique):
        donor = int(pair // n_groups)
        group = int(pair % n_groups)
        begin, end = int(sums.indptr[index]), int(sums.indptr[index + 1])
        values = np.rint(sums.data[begin:end]).astype(np.int32)
        target[ref_map[donor], group, sums.indices[begin:end]] += values


def count_table(meta: dict[str, Any], original: pd.DataFrame) -> np.ndarray:
    donor_lookup = {donor: index for index, donor in enumerate(meta["donor_ids"])}
    groups = [*meta["cell_types"], meta["pooled_name"]]
    group_lookup = {group: index for index, group in enumerate(groups)}
    counts = np.zeros((len(meta["donor_ids"]), len(groups), 4), dtype=np.int64)
    columns = ["n_cells_all", "n_cells_ge500", "n_cells_ge1500", "n_cells_ge3000"]
    for row in original.itertuples():
        donor = donor_lookup[str(row.donor_id)]
        group = group_lookup[str(row.group)]
        counts[donor, group] = [getattr(row, column) for column in columns]
    return counts


def recompute_rarefied(meta: dict[str, Any], original: pd.DataFrame,
                       projection_info: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    dataset = meta["dataset"]
    ds_code = 1 if dataset == "gao" else 2
    donor_ids = meta["donor_ids"]
    ages = meta["donors_df"]["age_years"].to_numpy(dtype=float)
    mature = np.flatnonzero(ages >= 30)
    reference_donors = np.flatnonzero((ages >= 30) & (ages <= 39))
    mature_mask = np.isin(meta["cell_donor_idx"], mature)
    groups = [*meta["cell_types"], meta["pooled_name"]]
    n_donors = len(donor_ids)
    n_groups = len(groups)
    n_raw_genes = len(meta["gene_symbols"])
    projection = projection_info["proj"]
    n_genes = projection.shape[1]
    ref_map = {int(donor): index for index, donor in enumerate(reference_donors)}
    ref_counts = {
        depth: np.zeros((len(reference_donors), n_groups, n_genes), dtype=np.int32)
        for depth in DEPTHS
    }
    counts = count_table(meta, original)
    a2_ct = np.zeros((n_donors * (n_groups - 1), n_genes), dtype=np.float64) if dataset == "gao" else None
    a2_pool = np.zeros((n_donors, n_genes), dtype=np.float64)
    a2_cell_counts = np.zeros((n_donors, n_groups), dtype=np.int64) if dataset == "gao" else None

    with h5py.File(meta["h5_path"], "r") as handle:
        group = matrix_group(handle, meta["matrix_slot"])
        for start in range(0, meta["n_obs"], CHUNK_ROWS):
            stop = min(start + CHUNK_ROWS, meta["n_obs"])
            local = mature_mask[start:stop]
            if not local.any():
                continue
            chunk = read_csr_chunk(group, start, stop, n_raw_genes, meta["indptr"])[local]
            gate = count_gate(chunk.data)
            if gate["nonfinite"] or gate["negative"] or gate["noninteger"]:
                raise ValueError(f"Count gate failed: {gate}")
            clean = (chunk @ projection).tocsr()
            del chunk
            clean.eliminate_zeros()
            umi = np.rint(np.asarray(clean.sum(axis=1)).ravel()).astype(np.int64)
            donors = meta["cell_donor_idx"][start:stop][local]
            cell_types = meta["cell_ct_idx"][start:stop][local]
            if dataset == "gao":
                accumulate_pairs(a2_ct, clean, donors, cell_types, n_groups - 1, False, include_pooled=False)
                accumulate_pairs(a2_pool, clean, donors, cell_types, n_groups, True)
                pos_umi = umi >= 1
                pairs = np.concatenate([
                    donors[pos_umi] * n_groups + cell_types[pos_umi],
                    donors[pos_umi] * n_groups + n_groups - 1,
                ])
                a2_cell_counts += np.bincount(pairs, minlength=n_donors * n_groups).reshape(n_donors, n_groups)
            else:
                accumulate_pairs(a2_pool, clean, donors, cell_types, n_groups, True)
            ref_local = np.isin(donors, reference_donors)
            for depth_label, depth_value in DEPTHS.items():
                floor = 1500 if depth_label == "500" else depth_value
                eligible = umi >= floor
                if not eligible.any():
                    continue
                ref_in_eligible = ref_local[eligible]
                if not ref_in_eligible.any():
                    continue
                seed_depth = 5001500 if depth_label == "500" else depth_value
                rng = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, seed_depth, start // CHUNK_ROWS, 0]))
                rarefied = multinomial_downsample(clean[eligible], depth_value, rng)
                rarefied.eliminate_zeros()
                accumulate_reference(
                    ref_counts[depth_label],
                    rarefied[ref_in_eligible],
                    donors[eligible][ref_in_eligible],
                    cell_types[eligible][ref_in_eligible],
                    ref_map,
                    n_groups,
                )
                del rarefied
            del clean
            if start % (CHUNK_ROWS * 25) == 0:
                print(f"{dataset} rarefied reference pass at {stop}/{meta['n_obs']}", flush=True)

    if dataset == "gao":
        expected = counts[:, :, 0]
        if not np.array_equal(a2_cell_counts, expected):
            raise ValueError("Gao recomputed all-cell counts differ from the committed donor table")
    references: dict[str, dict[int, dict[str, Any]]] = {}
    audit_rows = []
    for depth_label in DEPTHS:
        depth_value = DEPTHS[depth_label]
        floor_index = 2 if depth_label in {"1500", "500"} else 3
        references[depth_label], audit = build_reference_stats(
            meta, counts[:, :, floor_index], ref_counts[depth_label], reference_donors, n_groups, n_genes
        )
        for record in audit:
            group_index = int(record["group_index"])
            audit_rows.append({
                "dataset": dataset,
                "depth": depth_label,
                "group": groups[group_index],
                "n_reference_donors": int(record["n_reference_donors"]),
                "min_loo_reference_cells": int(record.get("min_loo_reference_cells", 0)),
                "eligible": bool(record["eligible"]),
                "reference_age_min": 30.0,
                "reference_age_max": 39.0,
            })
    del ref_counts
    gc.collect()

    distance_sums = {depth: np.zeros((n_donors, n_groups)) for depth in DEPTHS}
    distance_counts = {depth: np.zeros((n_donors, n_groups), dtype=np.int64) for depth in DEPTHS}
    with h5py.File(meta["h5_path"], "r") as handle:
        group = matrix_group(handle, meta["matrix_slot"])
        for start in range(0, meta["n_obs"], CHUNK_ROWS):
            stop = min(start + CHUNK_ROWS, meta["n_obs"])
            local = mature_mask[start:stop]
            if not local.any():
                continue
            chunk = read_csr_chunk(group, start, stop, n_raw_genes, meta["indptr"])[local]
            clean = (chunk @ projection).tocsr()
            del chunk
            clean.eliminate_zeros()
            umi = np.rint(np.asarray(clean.sum(axis=1)).ravel()).astype(np.int64)
            donors = meta["cell_donor_idx"][start:stop][local]
            cell_types = meta["cell_ct_idx"][start:stop][local]
            for depth_label, depth_value in DEPTHS.items():
                floor = 1500 if depth_label == "500" else depth_value
                floor_index = 2 if depth_label in {"1500", "500"} else 3
                selected = umi >= floor
                if not selected.any():
                    continue
                seed_depth = 5001500 if depth_label == "500" else depth_value
                rng = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, seed_depth, start // CHUNK_ROWS, 0]))
                rarefied = multinomial_downsample(clean[selected], depth_value, rng)
                rarefied.eliminate_zeros()
                selected_donors = donors[selected]
                selected_types = cell_types[selected]
                pairs = np.concatenate([
                    selected_donors * n_groups + selected_types,
                    selected_donors * n_groups + n_groups - 1,
                ])
                matrix = sp.vstack([rarefied, rarefied], format="csr")
                eligible = counts[:, :, floor_index] >= MIN_CELLS_PER_DONOR
                for group_index in range(n_groups):
                    if group_index not in references[depth_label]:
                        eligible[:, group_index] = False
                unique_pairs = np.unique(pairs)
                unique_pairs = unique_pairs[eligible.flat[unique_pairs]]
                for pair in unique_pairs:
                    donor = int(pair // n_groups)
                    group_index = int(pair % n_groups)
                    reference = reference_without_donor(references[depth_label][group_index], donor_ids[donor])
                    if reference["n_cells"] < MIN_LOO_REF_CELLS or reference["counts"].sum() <= 0:
                        continue
                    chosen = pairs == pair
                    distances = hellinger_distance(matrix[chosen], reference["counts"])
                    distance_sums[depth_label][donor, group_index] += float(distances.sum())
                    distance_counts[depth_label][donor, group_index] += len(distances)
                del rarefied, matrix
            del clean
            if start % (CHUNK_ROWS * 25) == 0:
                print(f"{dataset} rarefied scoring pass at {stop}/{meta['n_obs']}", flush=True)

    rows = []
    for donor in mature:
        for group_index, group_name in enumerate(groups):
            row: dict[str, Any] = {"donor_id": donor_ids[donor], "group": group_name}
            for depth_label in DEPTHS:
                count = int(distance_counts[depth_label][donor, group_index])
                row[f"a1_dist_{depth_label}"] = (
                    float(distance_sums[depth_label][donor, group_index] / count) if count else np.nan
                )
                reference = references[depth_label].get(group_index)
                row[f"n_reference_donors_{depth_label}"] = int(reference["n_reference_donors"]) if reference else 0
                row[f"min_loo_reference_cells_{depth_label}"] = int(reference["min_loo_reference_cells"]) if reference else 0
            rows.append(row)
    a2_inputs = {"a2_ct_gao": a2_ct, "a2_pool": a2_pool, "retained_symbols": projection_info["audit"]["retained_symbols"]}
    return pd.DataFrame(rows), pd.DataFrame(audit_rows), a2_inputs


def prepare_donors(meta: dict[str, Any], rarefied: pd.DataFrame) -> pd.DataFrame:
    dataset = meta["dataset"]
    original = pd.read_csv(RESULTS / f"followup-between-subject-{dataset}-donors.csv")
    original = original[original["age_years"] >= 30].copy()
    drop_old_young_ref_cols = [
        "a1_dist_full_1500", "a1_dist_1500", "a1_dist_500", "a1_dist_3000",
        "a1_dist_full_500_sens", "a1_dist_mature_ref_1500", "a1_self_non_loo",
        "a1_split_half_h2_to_h1", "eligible_a1_young_ref", "eligible_a2_oy", "eligible_a2_my",
    ]
    original = original.drop(columns=[col for col in drop_old_young_ref_cols if col in original.columns])
    original["bin"] = np.select(
        [original["age_years"] <= 39, original["age_years"] <= 59], ["YOUNG", "MID"], default="OLD"
    )
    original["mature_bin"] = np.select(
        [original["age_years"] <= 39, original["age_years"] <= 59],
        ["MATURE_30_39", "MID_40_59"], default="OLD_60_PLUS",
    )
    item8 = pd.read_csv(ITEM8_DISTANCES)
    item8 = item8[item8["dataset"].str.lower().str.startswith(dataset)].copy()
    full_columns = [
        "donor_id", "group", "n_cells_ge1500", "distance_mature_30_39", "n_reference_donors",
        "min_loo_reference_cells", "reference_age_min", "reference_age_max",
    ]
    merged = original.merge(item8[full_columns], on=["donor_id", "group"], how="left", validate="one_to_one")
    if not np.array_equal(merged["n_cells_ge1500_x"], merged["n_cells_ge1500_y"]):
        raise ValueError(f"{dataset}: Item 8 cell counts do not match the original donor table")
    merged = merged.drop(columns=["n_cells_ge1500_y"]).rename(columns={"n_cells_ge1500_x": "n_cells_ge1500"})
    merged["a1_dist_full_1500"] = merged["distance_mature_30_39"]
    merged["a1_full_1500_source"] = "results/followup-mature-nov17-donor-spread-distances.csv"
    merged = merged.merge(rarefied, on=["donor_id", "group"], how="left", validate="one_to_one")
    merged["eligible_a1_1500"] = (
        (merged["n_cells_ge1500"] >= MIN_CELLS_PER_DONOR)
        & merged[["a1_dist_full_1500", "a1_dist_1500", "a1_dist_500"]].notna().all(axis=1)
    )
    merged["eligible_a1_3000"] = (
        (merged["n_cells_ge3000"] >= MIN_CELLS_PER_DONOR) & merged["a1_dist_3000"].notna()
    )
    merged["eligible_a3_1500"] = (
        (merged["n_cells_ge1500"] >= MIN_CELLS_PER_DONOR)
        & merged[["a3_disp_1500", "a3_disp_full_1500"]].notna().all(axis=1)
    )
    if merged["age_years"].min() < 30:
        raise ValueError("Under-30 donor in mature donor table")
    if merged.loc[merged["a1_dist_full_1500"].notna(), "reference_age_min"].min() < 30:
        raise ValueError("Under-30 donor in full-depth reference")
    if merged.loc[merged["a1_dist_full_1500"].notna(), "reference_age_max"].max() > 39:
        raise ValueError("Over-39 donor in mature baseline reference")
    return merged


def fit_models(dataset: str, donors: pd.DataFrame) -> pd.DataFrame:
    is_gao = dataset == "gao"
    pooled_name = "All_LV_Nuclei" if is_gao else "All_Cells"
    batch_col = "study_or_pool"
    ref_level = "ENCODE v4 (Snyder)" if is_gao else None
    arms = [
        ("A1", "hellinger_mature_30_39_ref", "full_1500", "mature_30_39_loo", "a1_dist_full_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_mature_30_39_ref", "1500", "mature_30_39_loo", "a1_dist_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_mature_30_39_ref", "500", "mature_30_39_loo", "a1_dist_500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_mature_30_39_ref", "3000", "mature_30_39_loo", "a1_dist_3000", "n_cells_ge3000", "mean_log_umi_ge3000"),
        ("A3", "within_donor_loo_disp", "1500", "none", "a3_disp_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A3", "within_donor_loo_disp", "full_1500", "none", "a3_disp_full_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
    ]
    rows: list[dict[str, Any]] = []

    def specs(base: pd.DataFrame) -> list[tuple]:
        if is_gao:
            return [
                ("primary", base, ["sex_male", "_depth"], batch_col, 2, 5, ref_level),
                ("batch_adj_no_depth", base, ["sex_male"], batch_col, 2, 5, ref_level),
                ("age_only", base, [], None, 1, 5, None),
                ("age_plus_sex", base, ["sex_male"], None, 1, 5, None),
                ("10x_only_no_penn", base[base["technology"] != "Dropseq"], ["sex_male", "_depth"], batch_col, 2, 5, ref_level),
                ("studies_ge8", base[~base["study_or_pool"].isin(["Kuppe 2022", "Sim 2021"])], ["sex_male", "_depth"], batch_col, 2, 5, ref_level),
                ("encode_v4_only", base[base["study_or_pool"] == "ENCODE v4 (Snyder)"], ["sex_male", "_depth"], None, 1, 5, None),
            ]
        return [
            ("primary", base, ["sex_male", "_depth"], batch_col, 2, 10, None),
            ("pool_fe_no_depth", base, ["sex_male"], batch_col, 2, 10, None),
            ("no_pool", base, ["sex_male", "_depth"], None, 1, 5, None),
            ("age_ge_30", base[base["age_years"] >= 30], ["sex_male", "_depth"], batch_col, 2, 10, None),
        ]

    groups = sorted(donors["group"].unique(), key=lambda group: (group == pooled_name, group))
    for part, metric, depth, reference_type, outcome, cell_col, depth_col in arms:
        for group_name in groups:
            frame = donors[donors["group"] == group_name].copy()
            frame["n_cells"] = frame[cell_col]
            frame["_depth"] = frame[depth_col]
            eligibility = "eligible_a1_3000" if depth == "3000" else ("eligible_a3_1500" if part == "A3" else "eligible_a1_1500")
            frame = frame[frame[eligibility]].copy()
            if frame.empty:
                continue
            is_pool = bool(frame["is_pooled"].iloc[0])
            fallback = "not_evaluable_within_pool_df" if (not is_gao and not bool(frame["yazar_pool_fe_all_evaluable"].iloc[0])) else "not_evaluable_within_pool_df_slice"
            for spec_name, subset, nonbatch, batch, min_batch, min_df, reference_level in specs(frame):
                fitted = fit_ols_model_row(
                    subset, outcome, "age_per_decade", nonbatch, batch, min_batch, min_df, reference_level,
                    {
                        "dataset": dataset, "part": part, "metric": metric, "depth": depth,
                        "reference_type": reference_type, "model_spec": spec_name, "group": group_name,
                        "is_pooled": is_pool, "is_stratified_global": False, "eligibility_status": "eligible",
                        "default_fallback_status": fallback,
                    }, fallback_to_nobatch=not is_gao,
                )
                if fitted is None:
                    continue
                fitted["mature_30_39_mean"] = float(subset.loc[subset["bin"] == "YOUNG", outcome].mean()) if part == "A1" else np.nan
                fitted["mid_mean"] = float(subset.loc[subset["bin"] == "MID", outcome].mean()) if part == "A1" else np.nan
                fitted["old_mean"] = float(subset.loc[subset["bin"] == "OLD", outcome].mean()) if part == "A1" else np.nan
                if part == "A1":
                    old_mature = fit_contrast_adj_diff(subset, outcome, "YOUNG", "OLD", nonbatch, batch, min_batch, min_df, reference_level)
                    old_mid = fit_contrast_adj_diff(subset, outcome, "MID", "OLD", nonbatch, batch, min_batch, min_df, reference_level)
                else:
                    old_mature = old_mid = {"diff": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
                for prefix, values in (("old_minus_mature_30_39", old_mature), ("old_minus_mid", old_mid)):
                    fitted[f"{prefix}_diff_adj"] = values["diff"]
                    fitted[f"{prefix}_se_adj"] = values["se"]
                    fitted[f"{prefix}_ci95_low_adj"] = values["ci_low"]
                    fitted[f"{prefix}_ci95_high_adj"] = values["ci_high"]
                    fitted[f"{prefix}_p_adj"] = values["p"]
                rows.append(fitted)
        primary = [
            row for row in rows
            if row["part"] == part and row["depth"] == depth and row["model_spec"] == "primary" and not row["is_pooled"]
        ]
        source = donors[donors["eligible_a1_3000" if depth == "3000" else ("eligible_a3_1500" if part == "A3" else "eligible_a1_1500")]]
        stratified = _fit_stratified_global_row(
            source, [row["group"] for row in primary], outcome, cell_col, depth_col, dataset,
            part, metric, depth, reference_type, include_part_a_cols=False,
        )
        if stratified is not None:
            for column in ("mature_30_39_mean", "mid_mean", "old_mean"):
                stratified[column] = np.nan
            for prefix in ("old_minus_mature_30_39", "old_minus_mid"):
                for suffix in ("diff_adj", "se_adj", "ci95_low_adj", "ci95_high_adj", "p_adj"):
                    stratified[f"{prefix}_{suffix}"] = np.nan
            rows.append(stratified)
    result = apply_mature_slice_bh(pd.DataFrame(rows), ["dataset", "part", "metric", "depth", "reference_type", "model_spec"])
    verify_model_bh(result)
    return result


def apply_mature_slice_bh(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    out["q_bh_within_slice"] = np.nan
    out["family_size_m"] = np.nan
    eligible_mask = (~out["is_pooled"].astype(bool)) & (~out["is_stratified_global"].astype(bool))
    pool_fe_arm = out["model_spec"].astype(str).isin(["primary", "pool_fe", "pool_fe_no_depth", "age_ge_20", "age_ge_30"])
    is_fallback = out["yazar_pool_fe_status"].astype(str).str.contains("not_evaluable")
    eligible_mask = eligible_mask & (~(pool_fe_arm & is_fallback))
    for _, idx in out.loc[eligible_mask].groupby(group_cols, sort=False).groups.items():
        pvals = out.loc[idx, "pvalue"].to_numpy(dtype=float)
        valid = np.isfinite(pvals)
        m = int(valid.sum())
        if m == 0:
            continue
        out.loc[idx, "family_size_m"] = m
        valid_idx = np.asarray(list(idx))[valid]
        out.loc[valid_idx, "q_bh_within_slice"] = multipletests(pvals[valid], method="fdr_bh")[1]
    pool_idx = out["is_pooled"].astype(bool) & (~out["is_stratified_global"].astype(bool))
    out.loc[pool_idx, "family_size_m"] = 1
    out.loc[pool_idx, "q_bh_within_slice"] = out.loc[pool_idx, "pvalue"]
    return out


def verify_model_bh(frame: pd.DataFrame) -> None:
    eligible = (~frame["is_pooled"].astype(bool)) & (~frame["is_stratified_global"].astype(bool))
    pool_arms = frame["model_spec"].isin(["primary", "pool_fe", "pool_fe_no_depth", "age_ge_20", "age_ge_30"])
    eligible &= ~(pool_arms & frame["yazar_pool_fe_status"].astype(str).str.contains("not_evaluable"))
    for _, indices in frame[eligible].groupby(["dataset", "part", "metric", "depth", "reference_type", "model_spec"]).groups.items():
        expected = multipletests(frame.loc[indices, "pvalue"], method="fdr_bh")[1]
        if not np.allclose(frame.loc[indices, "q_bh_within_slice"], expected, rtol=1e-12, atol=1e-15):
            raise ValueError("Model BH verification failed")
    pooled = frame["is_pooled"] & ~frame["is_stratified_global"]
    if not np.allclose(frame.loc[pooled, "q_bh_within_slice"], frame.loc[pooled, "pvalue"], rtol=0, atol=0):
        raise ValueError("Pooled singleton BH verification failed")


def load_yazar_pool(a2_inputs: dict[str, Any], meta: dict[str, Any]) -> None:
    path = YAZAR_CACHE / "pb_All_Cells.pkl"
    with path.open("rb") as handle:
        raw = pickle.load(handle)
    raw.columns = raw.columns.astype(str)
    projection_info = build_human_gene_projection(raw.index.astype(str))
    if projection_info["audit"]["retained_symbols"] != a2_inputs["retained_symbols"]:
        raise ValueError("Yazar pooled cache gene projection differs from h5ad")
    matrix = (sp.csr_matrix(raw[meta["donor_ids"]].to_numpy(dtype=float).T) @ projection_info["proj"]).toarray()
    a2_inputs["a2_pool"] = matrix
    del raw


def run_expression_a2(meta: dict[str, Any], donors: pd.DataFrame,
                      a2_inputs: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    dataset = meta["dataset"]
    is_gao = dataset == "gao"
    donor_ids = meta["donor_ids"]
    donor_to_index = {donor: index for index, donor in enumerate(donor_ids)}
    cell_types = meta["cell_types"]
    pooled_name = meta["pooled_name"]
    n_cell_types = len(cell_types)
    if not is_gao:
        manifest = json.loads((YAZAR_CACHE / "pb_manifest.json").read_text())
        # Keep doublet-free singlet a2_pool accumulated in recompute_rarefied (matching original Pass 1).
        sample = next(iter(manifest.values()))
        with (YAZAR_CACHE / sample).open("rb") as handle:
            sample_frame = pickle.load(handle)
        projection_info = build_human_gene_projection(sample_frame.index.astype(str))
        if projection_info["audit"]["retained_symbols"] != a2_inputs["retained_symbols"]:
            raise ValueError("Yazar cache genes differ from h5ad genes")
        cache_projection = projection_info["proj"]
        del sample_frame
    else:
        manifest = {}
        cache_projection = None

    def load_counts(group: str, selected_donors: list[str]) -> np.ndarray:
        if group == pooled_name:
            return a2_inputs["a2_pool"][[donor_to_index[donor] for donor in selected_donors]]
        cell_index = cell_types.index(group)
        if is_gao:
            indices = [donor_to_index[donor] * n_cell_types + cell_index for donor in selected_donors]
            return a2_inputs["a2_ct_gao"][indices]
        with (YAZAR_CACHE / manifest[group]).open("rb") as handle:
            raw = pickle.load(handle)
        raw.columns = raw.columns.astype(str)
        dense = raw[selected_donors].to_numpy(dtype=float).T
        return (sp.csr_matrix(dense) @ cache_projection).toarray()

    rows = []
    ineligible = []
    for contrast_index, (contrast, reference_bin, comparison_bin, reference_age, comparison_age) in enumerate(CONTRASTS):
        unique_donors = sorted(donors.loc[donors["bin"].isin([reference_bin, comparison_bin]), "donor_id"].unique())
        donor_positions = {donor: index for index, donor in enumerate(unique_donors)}
        rng = np.random.default_rng(MASTER_SEED + 500 + contrast_index)
        global_ranks = np.array([rng.permutation(len(unique_donors)) for _ in range(A2_N_PERM)], dtype=np.int32)
        stratified: dict[tuple[str, float], list[dict[str, Any]]] = {}
        for group_name in [*cell_types, pooled_name]:
            frame = donors[
                (donors["group"] == group_name)
                & donors["bin"].isin([reference_bin, comparison_bin])
                & (donors["n_cells_all"] >= MIN_CELLS_PER_DONOR)
            ].copy()
            if not is_gao and group_name != pooled_name and group_name not in manifest:
                ineligible.append({"dataset": dataset, "contrast": contrast, "group": group_name, "reason": "not_in_yazar_cache"})
                continue
            n_reference = int((frame["bin"] == reference_bin).sum())
            n_comparison = int((frame["bin"] == comparison_bin).sum())
            if min(n_reference, n_comparison) < MIN_DONORS_A2_BIN:
                ineligible.append({
                    "dataset": dataset, "contrast": contrast, "group": group_name,
                    "n_reference": n_reference, "n_comparison": n_comparison, "reason": "lt6_in_reference_or_comparison",
                })
                continue
            order = frame["donor_id"].tolist()
            raw_counts = load_counts(group_name, order)
            detected = (raw_counts > 0).sum(axis=0) >= 0.5 * len(order)
            detected_counts = raw_counts[:, detected]
            totals = detected_counts.sum(axis=1, keepdims=True)
            expression = np.log1p(detected_counts / np.maximum(totals, 1.0) * 1e4)
            is_comparison = (frame["bin"] == comparison_bin).to_numpy(dtype=bool)
            depth_expression = _residualize_against_depth(expression, frame["mean_log_umi_all"].to_numpy(float))
            frame["_is_comp"] = is_comparison.astype(float)
            batch_frame, batch_columns, _, batch_status = _prune_covariate_matrix(
                frame, ["sex_male", "mean_log_umi_all"], "study_or_pool", 2,
                5 if is_gao else 10, "ENCODE v4 (Snyder)" if is_gao else None,
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
                pool_status = f"not_evaluable_within_{'study' if is_gao else 'pool'}_df_slice"
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
                    "dataset": dataset, "contrast": contrast, "reference_group": reference_bin,
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


def residualize_batch(frame: pd.DataFrame, value: str) -> tuple[pd.DataFrame, int]:
    work = frame[frame[value].notna()].copy()
    sizes = work.groupby("study_or_pool")["donor_id"].transform("size")
    singleton = int((sizes == 1).sum())
    work = work[sizes > 1].copy()
    if work.empty:
        # No donor survives the non-null / non-singleton filters (e.g. a
        # cell type with zero eligible distances, such as Epicardial and
        # Mast at full_1500 in Gao). Keep the expected output column so
        # callers see an empty frame and record UNTESTABLE instead of
        # raising KeyError: residualized_distance.
        work["residualized_distance"] = pd.Series(dtype=float)
        return work, singleton
    grand = work[value].mean()
    means = work.groupby("study_or_pool")[value].transform("mean")
    work["residualized_distance"] = work[value] - means + grand
    return work, singleton


def distance_comparison(a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> dict[str, float]:
    statistic, pvalue = stats.levene(a, b, center="median")
    var_a = float(np.var(a, ddof=1))
    var_b = float(np.var(b, ddof=1))
    bootstraps = np.empty(N_BOOT)
    for index in range(N_BOOT):
        sample_a = rng.choice(a, len(a), replace=True)
        sample_b = rng.choice(b, len(b), replace=True)
        bootstraps[index] = np.var(sample_a, ddof=1) / np.var(sample_b, ddof=1)
    low, high = np.percentile(bootstraps, [2.5, 97.5])
    return {
        "levene_stat": float(statistic), "levene_p": float(pvalue),
        "var_comparison": var_a, "var_reference": var_b, "var_ratio": var_a / var_b,
        "bootstrap_ci95_low": float(low), "bootstrap_ci95_high": float(high),
    }


def run_distance_a2(dataset: str, donors: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(BOOT_SEED)
    rows = []
    pooled_name = "All_LV_Nuclei" if dataset == "gao" else "All_Cells"
    versions = (
        ("mature_30_39_full_1500", "a1_dist_full_1500", donors, CONTRASTS),
        ("retained_mid_40_59_full_1500", "a1_dist_mature_ref_mid_loo", donors[donors["age_years"] >= 40], CONTRASTS[2:]),
    )
    for version, value, source, contrasts in versions:
        for group_name, group_frame in source.groupby("group"):
            is_pool = group_name == pooled_name
            for adjustment in ("raw", "batch_mean_residualized"):
                if adjustment == "raw":
                    adjusted = group_frame.copy()
                    adjusted["residualized_distance"] = adjusted[value]
                    singleton = 0
                else:
                    adjusted, singleton = residualize_batch(group_frame, value)
                for contrast, reference_bin, comparison_bin, reference_age, comparison_age in contrasts:
                    reference = adjusted.loc[adjusted["bin"] == reference_bin, "residualized_distance"].dropna().to_numpy(float)
                    comparison = adjusted.loc[adjusted["bin"] == comparison_bin, "residualized_distance"].dropna().to_numpy(float)
                    floor = 2 if is_pool else 20
                    base = {
                        "dataset": dataset, "distance_version": version, "adjustment": adjustment,
                        "contrast": contrast, "group": group_name, "is_pooled": is_pool,
                        "reference_age_range": reference_age, "comparison_age_range": comparison_age,
                        "n_donors_reference": len(reference), "n_donors_comparison": len(comparison),
                        "n_singleton_batch_donors_excluded": singleton,
                        "distance_source": "Item 8 full-depth distances" if version.startswith("mature") else "Original 40-59 reference distances",
                        "n_bootstrap": N_BOOT, "bootstrap_seed": BOOT_SEED,
                    }
                    if min(len(reference), len(comparison)) < floor:
                        rows.append({**base, "status": "UNTESTABLE", "reason": f"fewer_than_{floor}_donors_per_group"})
                        continue
                    rows.append({
                        **base, "status": "TESTED", "reason": "",
                        "sd_reference": float(np.std(reference, ddof=1)),
                        "sd_comparison": float(np.std(comparison, ddof=1)),
                        **distance_comparison(comparison, reference, rng),
                    })
    result = pd.DataFrame(rows)
    result["q_bh_levene"] = np.nan
    result["family_size_m"] = np.nan
    eligible = (result["status"] == "TESTED") & ~result["is_pooled"]
    for _, indices in result[eligible].groupby(["distance_version", "adjustment", "contrast"]).groups.items():
        expected = multipletests(result.loc[indices, "levene_p"], method="fdr_bh")[1]
        result.loc[indices, "q_bh_levene"] = expected
        result.loc[indices, "family_size_m"] = len(indices)
        if not np.allclose(result.loc[indices, "q_bh_levene"], expected, rtol=1e-12, atol=1e-15):
            raise ValueError("Distance BH verification failed")
    pooled = (result["status"] == "TESTED") & result["is_pooled"]
    result.loc[pooled, "q_bh_levene"] = result.loc[pooled, "levene_p"]
    result.loc[pooled, "family_size_m"] = 1
    return result


def run_dataset(dataset: str) -> None:
    start_utc = utc_now()
    started = time.monotonic()
    try:
        meta = reconcile_gao_metadata(GAO_H5AD) if dataset == "gao" else reconcile_yazar_metadata(YAZAR_H5AD, YAZAR_CACHE)
        original = pd.read_csv(RESULTS / f"followup-between-subject-{dataset}-donors.csv")
        original = original[original["age_years"] >= 30].copy()
        projection_info = build_human_gene_projection(meta["gene_symbols"])
        rarefied, reference_audit, a2_inputs = recompute_rarefied(meta, original, projection_info)
        donors = prepare_donors(meta, rarefied)
        donors.to_csv(RESULTS / f"followup-mature-between-subject-{dataset}-donors.csv", index=False)
        reference_audit.to_csv(RESULTS / f"followup-mature-between-subject-{dataset}-reference-audit.csv", index=False)
        print(f"{dataset}: wrote {len(donors)} mature donor rows", flush=True)
        models = fit_models(dataset, donors)
        models.to_csv(RESULTS / f"followup-mature-between-subject-{dataset}-models.csv", index=False)
        print(f"{dataset}: wrote {len(models)} A1/A3 model rows", flush=True)
        variance, ineligible = run_expression_a2(meta, donors, a2_inputs)
        variance.to_csv(RESULTS / f"followup-mature-between-subject-variance-{dataset}.csv", index=False)
        ineligible.to_csv(RESULTS / f"followup-mature-between-subject-variance-{dataset}-ineligible.csv", index=False)
        print(f"{dataset}: wrote {len(variance)} A2 expression rows", flush=True)
        distance = run_distance_a2(dataset, donors)
        distance.to_csv(RESULTS / f"followup-mature-between-subject-distance-variance-{dataset}.csv", index=False)
        details = (
            f"donors={donors['donor_id'].nunique()}; baseline_donors={donors.loc[donors['bin'] == 'YOUNG', 'donor_id'].nunique()}; "
            f"models={len(models)}; expression_a2={len(variance)}; distance_a2={len(distance)}"
        )
        append_runtime(f"mature_between_subject_{dataset}", start_utc, started, "completed", details)
        print(details, flush=True)
    except Exception as exc:
        append_runtime(f"mature_between_subject_{dataset}", start_utc, started, "FAILED", f"{type(exc).__name__}: {exc}")
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["gao", "yazar", "all"], default="all")
    args = parser.parse_args()
    for dataset in (["gao", "yazar"] if args.dataset == "all" else [args.dataset]):
        run_dataset(dataset)


if __name__ == "__main__":
    main()
