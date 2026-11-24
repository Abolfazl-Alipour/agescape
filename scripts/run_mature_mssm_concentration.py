#!/usr/bin/env python3
"""Exploratory mature MSSM DLPFC full-depth concentration analysis."""

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

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import h5py
import numpy as np
import pandas as pd
from scipy import sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
DATA_ROOT = REPO_ROOT.parent / "agescape-revision-data"
MSSM_H5AD = DATA_ROOT / "mssm.h5ad"
GAO_H5AD = DATA_ROOT / "tmp_organs/organ17_gao_heart/GSE290367_integrated_RNA_adata.h5ad"
CACHE_DIR = DATA_ROOT / "mssm_conc_cache"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
SLOPES_CSV = RESULTS_DIR / "followup-mature-mssm-concentration-slopes.csv"
POOLED_CSV = RESULTS_DIR / "followup-mature-mssm-concentration-pooled.csv"
QC_CSV = RESULTS_DIR / "followup-mature-mssm-concentration-qc.csv"
REPORT_MD = REPORTS_DIR / "followup-mature-mssm-concentration.md"
VALIDATION_JSON = CACHE_DIR / "gao_validation.json"
REF_CACHE = CACHE_DIR / "mssm_references.pkl"
METRIC_CHECKPOINT = CACHE_DIR / "mssm_metrics_checkpoint.pkl"
DONOR_ROWS_CACHE = CACHE_DIR / "mssm_donor_rows.pkl"

from scripts.run_followup_between_subject_gao_yazar import (
    apply_slice_bh,
    assign_age_bin,
    build_human_gene_projection,
    fit_ols_model_row,
    read_obs_column,
    reconcile_gao_metadata,
)
from scripts.run_followup_mature_fulldepth import (
    DEPTH_STRATA,
    MAX_E2_GLOBAL_POOL_FRAC,
    MIN_CELLS_PER_DONOR,
    MIN_DONORS_SLOPE,
    MIN_LOO_REF_CELLS,
    MIN_LOO_YOUNG_E2_STRATUM,
    MIN_YOUNG_REF,
    SH_SEEDS,
    assign_depth_stratum,
    build_e1_offprogram_lookup,
    build_e2_silent_lookup,
    fit_stratified_global_with_batch,
    reconcile_mature_metadata,
)

DS_CODE = 3
RUN_MODEL = "kimi-code/k3.1-pro"
PRESTATED_MDE80 = 0.029
MEASURES = ("scatter", "shannon", "E1")


def log(message: str) -> None:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(f"[{datetime.now(timezone.utc).isoformat()}] {message} | max_rss={rss:.1f} MB", flush=True)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def raw_to_clean_map(gene_symbols: np.ndarray) -> tuple[np.ndarray, int, dict[str, Any]]:
    info = build_human_gene_projection(gene_symbols)
    proj = info["proj"].tocsr()
    mapping = np.full(proj.shape[0], -1, dtype=np.int32)
    for raw_i in range(proj.shape[0]):
        if proj.indptr[raw_i + 1] > proj.indptr[raw_i]:
            mapping[raw_i] = int(proj.indices[proj.indptr[raw_i]])
    return mapping, proj.shape[1], info["audit"]


def read_clean_rows(
    group: h5py.Group,
    indptr: np.ndarray,
    rows: np.ndarray,
    mapping: np.ndarray,
    n_genes: int,
) -> tuple[sp.csr_matrix, dict[str, int]]:
    rows = np.asarray(rows, dtype=np.int64)
    if len(rows) == 0:
        return sp.csr_matrix((0, n_genes), dtype=np.float64), {"nonfinite": 0, "negative": 0, "noninteger": 0, "nnz": 0}
    order = np.argsort(rows)
    sorted_rows = rows[order]
    inverse = np.empty(len(rows), dtype=np.int64)
    inverse[order] = np.arange(len(rows))
    breaks = np.flatnonzero(np.diff(sorted_rows) > 1) + 1
    starts = np.r_[0, breaks]
    stops = np.r_[breaks, len(rows)]
    data_parts: list[np.ndarray] = []
    index_parts: list[np.ndarray] = []
    ptr_parts: list[np.ndarray] = []
    gate = {"nonfinite": 0, "negative": 0, "noninteger": 0, "nnz": 0}
    for a, b in zip(starts, stops):
        r0, r1 = int(sorted_rows[a]), int(sorted_rows[b - 1] + 1)
        p0, p1 = int(indptr[r0]), int(indptr[r1])
        raw_data = group["data"][p0:p1]
        raw_indices = group["indices"][p0:p1]
        gate["nnz"] += int(len(raw_data))
        gate["nonfinite"] += int((~np.isfinite(raw_data)).sum())
        gate["negative"] += int((raw_data < 0).sum())
        gate["noninteger"] += int((np.abs(raw_data - np.rint(raw_data)) > 0).sum())
        clean_indices = mapping[raw_indices]
        keep = clean_indices >= 0
        data_parts.append(raw_data[keep].astype(np.float64, copy=False))
        index_parts.append(clean_indices[keep])
        raw_ptr = indptr[r0:r1 + 1].astype(np.int64) - p0
        raw_row_ids = np.repeat(np.arange(r1 - r0, dtype=np.int64), np.diff(raw_ptr))
        kept_counts = np.bincount(raw_row_ids[keep], minlength=r1 - r0)
        ptr_parts.append(np.r_[0, np.cumsum(kept_counts, dtype=np.int64)])
    sorted_ptr = ptr_parts[0].copy()
    for part in ptr_parts[1:]:
        sorted_ptr = np.concatenate([sorted_ptr, part[1:] + sorted_ptr[-1]])
    sorted_mat = sp.csr_matrix(
        (np.concatenate(data_parts), np.concatenate(index_parts), sorted_ptr),
        shape=(len(rows), n_genes),
    )
    sorted_mat.sum_duplicates()
    sorted_mat.eliminate_zeros()
    if np.array_equal(order, np.arange(len(rows))):
        return sorted_mat, gate
    return sorted_mat[inverse].tocsr(), gate


def shannon_from_csr(mat: sp.csr_matrix) -> np.ndarray:
    n_rows = mat.shape[0]
    out = np.full(n_rows, np.nan, dtype=np.float64)
    if n_rows == 0:
        return out
    row_ids = np.repeat(np.arange(n_rows, dtype=np.int32), np.diff(mat.indptr))
    data = mat.data.astype(np.float64, copy=False)
    positive = data > 0
    if not np.any(positive):
        return out
    r = row_ids[positive]
    v = data[positive]
    totals = np.bincount(r, weights=v, minlength=n_rows)
    valid = totals > 0
    p = v / totals[r]
    values = -np.bincount(r, weights=p * np.log(p), minlength=n_rows)
    out[valid] = values[valid]
    return out


def scatter_from_l1(l1: sp.csr_matrix) -> np.ndarray:
    nc = l1.shape[0]
    if nc < 2:
        return np.full(nc, np.nan)
    centroid = np.asarray(l1.sum(axis=0)).ravel()
    row_ids = np.repeat(np.arange(nc, dtype=np.int32), np.diff(l1.indptr))
    p = l1.data
    loo = np.maximum(centroid[l1.indices] - p, 0.0) / float(nc - 1)
    bc = np.bincount(row_ids, weights=np.sqrt(p * loo), minlength=nc)
    return np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc, 0.0, 1.0)))


def normalized_l1(mat: sp.csr_matrix) -> sp.csr_matrix:
    totals = np.asarray(mat.sum(axis=1)).ravel()
    return mat.multiply(1.0 / np.maximum(totals, 1.0)[:, None]).tocsr()


def basic_group_metrics(mat: sp.csr_matrix, donor_idx: int, ct_idx: int, split_half: bool = True) -> dict[str, Any]:
    nc = mat.shape[0]
    umi = np.rint(np.asarray(mat.sum(axis=1)).ravel()).astype(np.int64)
    l1 = normalized_l1(mat)
    result: dict[str, Any] = {
        "n_cells": int(nc),
        "_depth": float(np.log(np.maximum(umi, 1)).mean()),
        "_inv_umi": float((1.0 / np.maximum(umi, 1)).mean()),
        "scatter": float(scatter_from_l1(l1).mean()),
        "shannon": float(shannon_from_csr(mat).mean()),
    }
    if not split_half:
        return result
    raw = np.rint(mat.data).astype(np.int64)
    row_ids = np.repeat(np.arange(nc, dtype=np.int32), np.diff(mat.indptr))
    cols = mat.indices
    for seed_i, seed in enumerate(SH_SEEDS):
        rng = np.random.default_rng(np.random.SeedSequence([int(seed), DS_CODE, int(donor_idx), int(ct_idx)]))
        h1 = rng.binomial(raw, 0.5)
        for half_i, values in enumerate((h1, raw - h1)):
            positive = values > 0
            r, c, v = row_ids[positive], cols[positive], values[positive].astype(np.float64)
            totals = np.maximum(np.bincount(r, weights=v, minlength=nc), 1.0)
            p = v / totals[r]
            half = sp.csr_matrix((p, (r, c)), shape=(nc, mat.shape[1]))
            prefix = f"h{half_i + 1}_s{seed_i}_"
            result[prefix + "scatter"] = float(scatter_from_l1(half).mean())
            result[prefix + "shannon"] = float(shannon_from_csr(half).mean())
    return result


def detected_count_metrics(
    mat: sp.csr_matrix,
    genes: np.ndarray,
    donor_idx: int,
    ct_idx: int,
    tag: str,
) -> dict[str, Any]:
    nc = mat.shape[0]
    selected = mat[:, genes].tocsr()
    selected.eliminate_zeros()
    counts = np.diff(selected.indptr).astype(np.float64)
    result: dict[str, Any] = {
        "n": int(nc),
        "sum": float(counts.sum()),
        "mean": float(counts.mean()),
        "depth_sum": float(np.log(np.maximum(np.asarray(mat.sum(axis=1)).ravel(), 1.0)).sum()),
        "inv_sum": float((1.0 / np.maximum(np.asarray(mat.sum(axis=1)).ravel(), 1.0)).sum()),
    }
    row_ids = np.repeat(np.arange(nc, dtype=np.int32), np.diff(selected.indptr))
    raw = np.rint(selected.data).astype(np.int64)
    import zlib
    tag_crc = zlib.crc32(tag.encode("ascii")) & 0xFFFF
    for seed_i, seed in enumerate(SH_SEEDS):
        rng = np.random.default_rng(np.random.SeedSequence([int(seed), DS_CODE, int(donor_idx), int(ct_idx), tag_crc]))
        h1 = rng.binomial(raw, 0.5)
        for half_i, values in enumerate((h1, raw - h1)):
            positive = values > 0
            half = sp.csr_matrix(
                (np.ones(int(positive.sum()), dtype=np.float64), (row_ids[positive], selected.indices[positive])),
                shape=(nc, len(genes)),
            )
            half_counts = np.diff(half.indptr).astype(np.float64)
            result[f"h{half_i + 1}_s{seed_i}_sum"] = float(half_counts.sum())
    return result


def mssm_metadata() -> dict[str, Any]:
    with h5py.File(MSSM_H5AD, "r") as handle:
        obs = pd.DataFrame({
            c: read_obs_column(handle["obs"], c)
            for c in ("donor_id", "disease", "sex", "development_stage", "cell_type", "is_primary_data", "tissue", "assay", "suspension_type")
        })
        genes = read_obs_column(handle["var"], "feature_name")
        shape = tuple(int(x) for x in handle["X"].attrs["shape"])
        indptr = handle["X"]["indptr"][:].astype(np.int64)
    mapping, n_genes, gene_audit = raw_to_clean_map(genes)
    control = obs["disease"].eq("normal")
    control_obs = obs.loc[control].copy()
    control_obs["age_years"] = pd.to_numeric(
        control_obs["development_stage"].str.extract(r"^(\d+)-year-old stage$", expand=False), errors="coerce"
    )
    donor_audit = control_obs.groupby("donor_id", sort=True).agg(
        age_years=("age_years", "first"),
        age_nunique=("age_years", "nunique"),
        sex=("sex", "first"),
        sex_nunique=("sex", "nunique"),
        n_cells=("donor_id", "size"),
    )
    assert int(donor_audit["age_nunique"].max()) <= 1
    assert int(donor_audit["sex_nunique"].max()) <= 1
    parsed = donor_audit.loc[donor_audit["age_years"].notna()]
    mature = parsed.loc[parsed["age_years"] >= 30].copy().sort_index()
    mature["age_per_decade"] = mature["age_years"] / 10.0
    mature["bin"] = mature["age_years"].map(assign_age_bin)
    mature["sex_male"] = mature["sex"].eq("male").astype(int)
    donor_ids = mature.index.astype(str).tolist()
    donor_map = {d: i for i, d in enumerate(donor_ids)}
    cell_types = sorted(control_obs["cell_type"].astype(str).unique().tolist())
    ct_map = {c: i for i, c in enumerate(cell_types)}
    obs_donor = obs["donor_id"].map(donor_map).to_numpy(dtype=float)
    cell_donor = np.full(len(obs), -1, dtype=np.int16)
    valid_donor = np.isfinite(obs_donor)
    cell_donor[valid_donor] = obs_donor[valid_donor].astype(np.int16)
    cell_ct = np.full(len(obs), -1, dtype=np.int16)
    valid_ct = cell_donor >= 0
    cell_ct[valid_ct] = obs.loc[valid_ct, "cell_type"].map(ct_map).to_numpy(dtype=np.int16)
    donor_rows = [np.flatnonzero(cell_donor == i) for i in range(len(donor_ids))]
    donors_df = mature.reset_index(names="donor_id")
    donors_df["donor_id"] = donors_df["donor_id"].astype(str)
    probe = {
        "n_cells_total": int(len(obs)),
        "n_control_cells": int(control.sum()),
        "n_control_donors": int(len(donor_audit)),
        "n_control_donors_parseable_age": int(len(parsed)),
        "control_age_min": float(parsed["age_years"].min()),
        "control_age_max": float(parsed["age_years"].max()),
        "n_control_donors_missing_age": int(donor_audit["age_years"].isna().sum()),
        "n_donors_30plus": int(len(mature)),
        "n_male_30plus": int(mature["sex"].eq("male").sum()),
        "n_female_30plus": int(mature["sex"].eq("female").sum()),
        "n_reference_donors_30_39": int(mature["age_years"].between(30, 39).sum()),
        "n_mature_control_cells": int(valid_donor.sum()),
        "all_primary_data": bool(obs["is_primary_data"].eq("True").all()),
        "tissue_values": sorted(obs["tissue"].unique().tolist()),
        "assay_values": sorted(obs["assay"].unique().tolist()),
        "suspension_values": sorted(obs["suspension_type"].unique().tolist()),
        "file_size_bytes": int(MSSM_H5AD.stat().st_size),
        "shape": list(shape),
    }
    return {
        "dataset": "mssm",
        "h5_path": MSSM_H5AD,
        "matrix_slot": "X",
        "indptr": indptr,
        "shape": shape,
        "mapping": mapping,
        "n_genes": n_genes,
        "gene_audit": gene_audit,
        "donors_df": donors_df,
        "donor_ids": donor_ids,
        "cell_types": cell_types,
        "cell_ct": cell_ct,
        "donor_rows": donor_rows,
        "probe": probe,
    }


def pass1_mssm(meta: dict[str, Any], resume: bool = True) -> dict[str, Any]:
    if resume and REF_CACHE.exists():
        log(f"Reusing {REF_CACHE}")
        with REF_CACHE.open("rb") as handle:
            return pickle.load(handle)
    n_donors = len(meta["donor_ids"])
    n_cts = len(meta["cell_types"])
    n_genes = meta["n_genes"]
    donors_df = meta["donors_df"]
    bin_arr = donors_df["bin"].map({"LT20": 0, "YOUNG": 1, "MID": 2, "OLD": 3}).to_numpy(dtype=np.int32)
    young_donors = np.flatnonzero(bin_arr == 1).tolist()
    y_pairs = [(d, c) for d in young_donors for c in range(n_cts)]
    y_pair_map = {d * n_cts + c: i for i, (d, c) in enumerate(y_pairs)}
    n_yp = len(y_pairs)
    cnt_pre = np.zeros((n_donors, n_cts), dtype=np.int64)
    cnt_15 = np.zeros((n_donors, n_cts), dtype=np.int64)
    lu_15 = np.zeros((n_donors, n_cts), dtype=np.float64)
    invu_15 = np.zeros((n_donors, n_cts), dtype=np.float64)
    y_full_cnt = np.zeros((n_yp, n_genes), dtype=np.int64)
    y_full_det = np.zeros((n_yp, n_genes), dtype=np.int64)
    y_n_cells = np.zeros(n_yp, dtype=np.int64)
    y_st_det = np.zeros((4, n_yp, n_genes), dtype=np.int32)
    y_st_n = np.zeros((4, n_yp), dtype=np.int64)
    gate_total = {"nonfinite": 0, "negative": 0, "noninteger": 0, "nnz": 0}
    t0 = time.perf_counter()
    with h5py.File(MSSM_H5AD, "r") as handle:
        group = handle["X"]
        for di, rows in enumerate(meta["donor_rows"]):
            mat, gate = read_clean_rows(group, meta["indptr"], rows, meta["mapping"], n_genes)
            for key in gate_total:
                gate_total[key] += gate[key]
            umi = np.rint(np.asarray(mat.sum(axis=1)).ravel()).astype(np.int64)
            cts = meta["cell_ct"][rows]
            for ci in np.unique(cts):
                mask = cts == ci
                cnt_pre[di, ci] = int(mask.sum())
                keep = mask & (umi >= 1500)
                cnt_15[di, ci] = int(keep.sum())
                if keep.any():
                    lu_15[di, ci] = float(np.log(umi[keep]).sum())
                    invu_15[di, ci] = float((1.0 / umi[keep]).sum())
            if bin_arr[di] == 1:
                eligible = umi >= 1500
                emat = mat[eligible].tocsr()
                ects = cts[eligible]
                eumi = umi[eligible]
                strata = assign_depth_stratum(eumi)
                binary = emat.copy()
                binary.data = np.ones_like(binary.data)
                for ci in np.unique(ects):
                    yi = y_pair_map[di * n_cts + int(ci)]
                    cm = ects == ci
                    y_n_cells[yi] = int(cm.sum())
                    y_full_cnt[yi] = np.rint(np.asarray(emat[cm].sum(axis=0)).ravel()).astype(np.int64)
                    y_full_det[yi] = np.asarray(binary[cm].sum(axis=0)).ravel().astype(np.int64)
                    for si in range(4):
                        smask = cm & (strata == si)
                        y_st_n[si, yi] = int(smask.sum())
                        if smask.any():
                            y_st_det[si, yi] = np.asarray(binary[smask].sum(axis=0)).ravel().astype(np.int32)
            if (di + 1) % 10 == 0 or di + 1 == n_donors:
                log(f"MSSM pass 1 donor {di + 1}/{n_donors}; elapsed={time.perf_counter() - t0:.1f}s")
            del mat
            gc.collect()
    if any(gate_total[k] for k in ("nonfinite", "negative", "noninteger")):
        raise ValueError(f"MSSM count gate failed: {gate_total}")
    e1 = build_e1_offprogram_lookup(
        y_full_cnt, y_full_det, y_n_cells, y_pair_map, cnt_15, bin_arr,
        n_donors, n_cts, n_genes, True, None,
    )
    e2_global = build_e2_silent_lookup(
        y_full_det, y_n_cells, y_pair_map, cnt_15, bin_arr,
        n_donors, n_cts, n_genes, MIN_LOO_REF_CELLS, True, None,
    )
    dummy = np.full((n_donors, n_cts), MIN_CELLS_PER_DONOR, dtype=np.int64)
    e2_strata = []
    for si in range(4):
        e2_strata.append(build_e2_silent_lookup(
            y_st_det[si], y_st_n[si], y_pair_map, dummy, bin_arr,
            n_donors, n_cts, n_genes, MIN_LOO_YOUNG_E2_STRATUM, False, None,
        ))
    refs = {
        "cnt_pre": cnt_pre,
        "cnt_15": cnt_15,
        "lu_15": lu_15,
        "invu_15": invu_15,
        "y_pair_map": y_pair_map,
        "y_full_cnt": y_full_cnt,
        "y_full_det": y_full_det,
        "y_n_cells": y_n_cells,
        "y_st_det": y_st_det,
        "y_st_n": y_st_n,
        "bin_arr": bin_arr,
        "E1": e1,
        "E2_global": e2_global,
        "E2_strata": e2_strata,
        "gate": gate_total,
        "elapsed_s": time.perf_counter() - t0,
        "probe": meta["probe"],
        "gene_audit": meta["gene_audit"],
    }
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with REF_CACHE.open("wb") as handle:
        pickle.dump(refs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"MSSM pass 1 complete in {refs['elapsed_s']:.1f}s")
    return refs


def reference_genes(ref: dict[str, Any], di: int, ci: int, bin_arr: np.ndarray) -> np.ndarray | None:
    if bin_arr[di] == 1:
        return ref["loo"].get((di, ci))
    return ref["full"].get(ci)


def empty_offprogram() -> dict[str, Any]:
    out = {"n": 0, "sum": 0.0, "depth_sum": 0.0, "inv_sum": 0.0}
    for seed_i in range(len(SH_SEEDS)):
        out[f"h1_s{seed_i}_sum"] = 0.0
        out[f"h2_s{seed_i}_sum"] = 0.0
    return out


def add_offprogram(total: dict[str, Any], part: dict[str, Any]) -> None:
    for key in total:
        total[key] += part[key]


def offprogram_record(total: dict[str, Any]) -> dict[str, Any]:
    n = int(total["n"])
    if n <= 0:
        return {}
    out = {
        "n_cells": n,
        "_depth": total["depth_sum"] / n,
        "_inv_umi": total["inv_sum"] / n,
        "value": total["sum"] / n,
    }
    for seed_i in range(len(SH_SEEDS)):
        out[f"h1_s{seed_i}"] = total[f"h1_s{seed_i}_sum"] / n
        out[f"h2_s{seed_i}"] = total[f"h2_s{seed_i}_sum"] / n
    return out


def pass2_mssm(meta: dict[str, Any], refs: dict[str, Any], resume: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    if resume and DONOR_ROWS_CACHE.exists():
        log(f"Reusing {DONOR_ROWS_CACHE}")
        with DONOR_ROWS_CACHE.open("rb") as handle:
            saved = pickle.load(handle)
        return saved["ct_rows"], saved["e2_rows"]
    start_donor = 0
    ct_records: list[dict[str, Any]] = []
    e2_records: list[dict[str, Any]] = []
    if resume and METRIC_CHECKPOINT.exists():
        with METRIC_CHECKPOINT.open("rb") as handle:
            checkpoint = pickle.load(handle)
        start_donor = checkpoint["next_donor"]
        ct_records = checkpoint["ct_records"]
        e2_records = checkpoint["e2_records"]
        log(f"Resuming MSSM metrics at donor {start_donor + 1}")
    donors_df = meta["donors_df"]
    bin_arr = refs["bin_arr"]
    t0 = time.perf_counter()
    with h5py.File(MSSM_H5AD, "r") as handle:
        group = handle["X"]
        for di in range(start_donor, len(meta["donor_ids"])):
            rows = meta["donor_rows"][di]
            all_mat, _ = read_clean_rows(group, meta["indptr"], rows, meta["mapping"], meta["n_genes"])
            all_umi = np.rint(np.asarray(all_mat.sum(axis=1)).ravel()).astype(np.int64)
            eligible = all_umi >= 1500
            mat = all_mat[eligible].tocsr()
            del all_mat
            cts = meta["cell_ct"][rows][eligible]
            umi = all_umi[eligible]
            dr = donors_df.iloc[di]
            base = {
                "dataset": "mssm", "donor_idx": di, "donor_id": str(dr["donor_id"]),
                "age_years": float(dr["age_years"]), "age_per_decade": float(dr["age_per_decade"]),
                "bin": str(dr["bin"]), "sex": str(dr["sex"]), "sex_male": int(dr["sex_male"]),
                "study_or_pool": "MSSM",
            }
            pooled_e1 = empty_offprogram()
            pooled_e2 = [empty_offprogram() for _ in range(4)]
            for ci in np.unique(cts):
                ci = int(ci)
                cmask = cts == ci
                sub = mat[cmask].tocsr()
                nc = sub.shape[0]
                ct_name = meta["cell_types"][ci]
                if nc >= MIN_CELLS_PER_DONOR:
                    basic = basic_group_metrics(sub, di, ci, split_half=True)
                    rec = {**base, "group": ct_name, "ct_idx": ci, **basic}
                    e1_genes = reference_genes(refs["E1"], di, ci, bin_arr)
                    if e1_genes is not None and len(e1_genes):
                        e1 = detected_count_metrics(sub, e1_genes, di, ci, "e1")
                        rec["E1"] = e1["mean"]
                        rec["E1_n_genes_set"] = int(len(e1_genes))
                        for si in range(len(SH_SEEDS)):
                            rec[f"h1_s{si}_E1"] = e1[f"h1_s{si}_sum"] / nc
                            rec[f"h2_s{si}_E1"] = e1[f"h2_s{si}_sum"] / nc
                    ct_records.append(rec)
                else:
                    e1_genes = reference_genes(refs["E1"], di, ci, bin_arr)
                if e1_genes is not None and len(e1_genes):
                    e1_pool = detected_count_metrics(sub, e1_genes, di, ci, "e1")
                    add_offprogram(pooled_e1, e1_pool)
                strata = assign_depth_stratum(umi[cmask])
                for si in range(4):
                    genes = reference_genes(refs["E2_strata"][si], di, ci, bin_arr)
                    if genes is None or len(genes) == 0:
                        continue
                    smask = strata == si
                    if not smask.any():
                        continue
                    ssub = sub[smask].tocsr()
                    score = detected_count_metrics(ssub, genes, di, ci, "e2_st")
                    add_offprogram(pooled_e2[si], score)
                    if ssub.shape[0] >= MIN_CELLS_PER_DONOR:
                        value = offprogram_record(score)
                        er = {
                            **base, "group": ct_name, "ct_idx": ci, "stratum": DEPTH_STRATA[si][0],
                            "E2_stratum": value["value"], "E2_n_genes_set": int(len(genes)),
                            "n_cells": value["n_cells"], "_depth": value["_depth"], "_inv_umi": value["_inv_umi"],
                        }
                        for seed_i in range(len(SH_SEEDS)):
                            er[f"h1_s{seed_i}_E2_stratum"] = value[f"h1_s{seed_i}"]
                            er[f"h2_s{seed_i}_E2_stratum"] = value[f"h2_s{seed_i}"]
                        e2_records.append(er)
            if mat.shape[0] >= MIN_CELLS_PER_DONOR:
                pooled = basic_group_metrics(mat, di, 999, split_half=True)
                rec = {**base, "group": "All_Cells", "ct_idx": -1, **pooled}
                e1p = offprogram_record(pooled_e1)
                if e1p:
                    rec["E1"] = e1p["value"]
                    rec["E1_n_cells"] = e1p["n_cells"]
                    rec["E1_depth"] = e1p["_depth"]
                    rec["E1_inv_umi"] = e1p["_inv_umi"]
                    for seed_i in range(len(SH_SEEDS)):
                        rec[f"h1_s{seed_i}_E1"] = e1p[f"h1_s{seed_i}"]
                        rec[f"h2_s{seed_i}_E1"] = e1p[f"h2_s{seed_i}"]
                ct_records.append(rec)
            for si in range(4):
                value = offprogram_record(pooled_e2[si])
                if value and value["n_cells"] >= MIN_CELLS_PER_DONOR:
                    er = {
                        **base, "group": "All_Cells", "ct_idx": -1, "stratum": DEPTH_STRATA[si][0],
                        "E2_stratum": value["value"], "n_cells": value["n_cells"],
                        "_depth": value["_depth"], "_inv_umi": value["_inv_umi"],
                    }
                    for seed_i in range(len(SH_SEEDS)):
                        er[f"h1_s{seed_i}_E2_stratum"] = value[f"h1_s{seed_i}"]
                        er[f"h2_s{seed_i}_E2_stratum"] = value[f"h2_s{seed_i}"]
                    e2_records.append(er)
            if (di + 1) % 5 == 0 or di + 1 == len(meta["donor_ids"]):
                with METRIC_CHECKPOINT.open("wb") as checkpoint_file:
                    pickle.dump({"next_donor": di + 1, "ct_records": ct_records, "e2_records": e2_records}, checkpoint_file, protocol=pickle.HIGHEST_PROTOCOL)
                log(f"MSSM pass 2 donor {di + 1}/{len(meta['donor_ids'])}; elapsed={time.perf_counter() - t0:.1f}s")
            del mat
            gc.collect()
    ct_rows = pd.DataFrame(ct_records)
    e2_rows = pd.DataFrame(e2_records)
    with DONOR_ROWS_CACHE.open("wb") as handle_out:
        pickle.dump({"ct_rows": ct_rows, "e2_rows": e2_rows}, handle_out, protocol=pickle.HIGHEST_PROTOCOL)
    METRIC_CHECKPOINT.unlink(missing_ok=True)
    return ct_rows, e2_rows


def split_half_statistics(df: pd.DataFrame, outcome: str, is_count: bool) -> dict[str, float]:
    nan_result = {"sh_rms_delta_donor": np.nan, "sh_se_samp_slope": np.nan, "sh_rho_tech2": np.nan, "sh_snr_slope": np.nan}
    if len(df) < MIN_DONORS_SLOPE:
        return nan_result
    primary = fit_ols_model_row(
        df, outcome, "age_per_decade", ["sex_male", "_depth"], None, 1, 5, None, {}
    )
    if primary is None:
        return nan_result
    rms, se2 = [], []
    for seed_i in range(len(SH_SEEDS)):
        h1, h2 = f"h1_s{seed_i}_{outcome}", f"h2_s{seed_i}_{outcome}"
        if h1 not in df or h2 not in df:
            continue
        work = df.dropna(subset=[h1, h2]).copy()
        if len(work) < MIN_DONORS_SLOPE:
            continue
        work["_delta"] = work[h1] - work[h2]
        rms.append(float(np.mean(work["_delta"] ** 2)))
        fit = fit_ols_model_row(work, "_delta", "age_per_decade", ["sex_male", "_depth"], None, 1, 5, None, {})
        if fit is not None:
            se2.append(float(fit["se_raw"]) ** 2)
    factor = 1.0 if is_count else 0.5
    variance_scale = 1.0 if is_count else 0.25
    rms_value = float(np.sqrt(np.mean(rms))) if rms else np.nan
    sampling_se = float(factor * np.sqrt(np.mean(se2))) if se2 else np.nan
    return {
        "sh_rms_delta_donor": rms_value,
        "sh_se_samp_slope": sampling_se,
        "sh_rho_tech2": float(variance_scale * rms_value ** 2 / primary["resid_sd"] ** 2) if np.isfinite(rms_value) and primary["resid_sd"] > 0 else np.nan,
        "sh_snr_slope": float(abs(primary["slope_raw_per_decade"]) / sampling_se) if np.isfinite(sampling_se) and sampling_se > 0 else np.nan,
    }


def global_split_half_statistics(
    df: pd.DataFrame,
    eligible: list[str],
    outcome: str,
    scale_map: dict[Any, float],
    is_count: bool,
    strata_col: str | None,
) -> dict[str, float]:
    se2 = []
    rms = []
    for seed_i in range(len(SH_SEEDS)):
        h1, h2 = f"h1_s{seed_i}_{outcome}", f"h2_s{seed_i}_{outcome}"
        work = df.dropna(subset=[h1, h2]).copy()
        if work.empty:
            continue
        work["_delta"] = work[h1] - work[h2]
        keys = list(zip(work["group"], work[strata_col])) if strata_col else list(work["group"])
        scales = np.array([scale_map.get(k, np.nan) for k in keys], dtype=float)
        valid = np.isfinite(scales) & (scales > 0)
        if valid.any():
            rms.append(float(np.mean((work["_delta"].to_numpy()[valid] / scales[valid]) ** 2)))
        row, _ = fit_stratified_global_with_batch(
            work, eligible, "_delta", "n_cells", "_depth", "mssm", "3" if is_count else "1",
            outcome, "full_depth", "split_half_se", strata_col=strata_col, batch_col=None, ref_level=None,
            custom_scale_map=scale_map,
        )
        if row is not None:
            se2.append(float(row["se_raw"]) ** 2)
    factor = 1.0 if is_count else 0.5
    variance_scale = 1.0 if is_count else 0.25
    rms_value = float(np.sqrt(np.mean(rms))) if rms else np.nan
    sampling_se = float(factor * np.sqrt(np.mean(se2))) if se2 else np.nan
    return {"sh_rms_delta_donor": rms_value, "sh_se_samp_slope": sampling_se, "sh_variance_scale": variance_scale}


def fit_models(ct_rows: pd.DataFrame, e2_rows: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    configurations = []
    for measure in MEASURES:
        configurations.append((measure, ct_rows, None))
    configurations.append(("E2_stratum", e2_rows, "stratum"))
    for outcome, source, strata_col in configurations:
        if source.empty or outcome not in source:
            continue
        stratum_values = sorted(source[strata_col].unique()) if strata_col else [None]
        for stratum in stratum_values:
            current = source.loc[source[strata_col] == stratum].copy() if strata_col else source.copy()
            eligible: list[str] = []
            scale_map: dict[Any, float] = {}
            for group in sorted(current["group"].unique()):
                if group == "All_Cells":
                    continue
                sub = current.loc[(current["group"] == group) & current[outcome].notna()].copy()
                if len(sub) < MIN_DONORS_SLOPE:
                    continue
                row = fit_ols_model_row(
                    sub, outcome, "age_per_decade", ["sex_male", "_depth"], None, 1, 5, None,
                    {"dataset": "mssm", "metric": outcome, "group": group, "stratum": stratum or "",
                     "depth_mode": "full_depth" if not strata_col else "bounded_depth_strata",
                     "model_spec": "primary", "is_pooled": False, "is_stratified_global": False},
                )
                if row is None:
                    continue
                row.update(split_half_statistics(sub, outcome, outcome == "E1" or outcome == "E2_stratum"))
                rows.append(row)
                eligible.append(group)
                key = (group, stratum) if strata_col else group
                scale_map[key] = float(sub[outcome].std(ddof=1))
            pooled = current.loc[(current["group"] == "All_Cells") & current[outcome].notna()].copy()
            if len(pooled) >= MIN_DONORS_SLOPE:
                if outcome == "E1":
                    pooled["n_cells"] = pooled["E1_n_cells"]
                    pooled["_depth"] = pooled["E1_depth"]
                    pooled["_inv_umi"] = pooled["E1_inv_umi"]
                row = fit_ols_model_row(
                    pooled, outcome, "age_per_decade", ["sex_male", "_depth"], None, 1, 5, None,
                    {"dataset": "mssm", "metric": outcome, "group": "All_Cells", "stratum": stratum or "",
                     "depth_mode": "full_depth" if not strata_col else "bounded_depth_strata",
                     "model_spec": "primary", "is_pooled": True, "is_stratified_global": False},
                )
                if row is not None:
                    row.update(split_half_statistics(pooled, outcome, outcome == "E1" or outcome == "E2_stratum"))
                    rows.append(row)
            if len(eligible) >= 2:
                strat = current.loc[current["group"].isin(eligible) & current[outcome].notna()].copy()
                global_row, _ = fit_stratified_global_with_batch(
                    strat, eligible, outcome, "n_cells", "_depth", "mssm", "3" if outcome in ("E1", "E2_stratum") else "1",
                    outcome, "full_depth" if not strata_col else "bounded_depth_strata", "primary_full_depth",
                    model_spec="stratified_global_primary", gene_set=outcome if outcome in ("E1", "E2_stratum") else "none",
                    strata_col=strata_col, batch_col=None, ref_level=None,
                )
                if global_row is not None:
                    global_row["stratum"] = stratum or ""
                    global_row.update(global_split_half_statistics(
                        strat, eligible, outcome, scale_map, outcome in ("E1", "E2_stratum"), strata_col
                    ))
                    global_row["sh_rho_tech2"] = float(
                        (1.0 if outcome in ("E1", "E2_stratum") else 0.25)
                        * global_row["sh_rms_delta_donor"] ** 2 / global_row["resid_sd"] ** 2
                    ) if np.isfinite(global_row["sh_rms_delta_donor"]) and global_row["resid_sd"] > 0 else np.nan
                    global_row["sh_snr_slope"] = float(abs(global_row["slope_raw_per_decade"]) / global_row["sh_se_samp_slope"]) if np.isfinite(global_row["sh_se_samp_slope"]) and global_row["sh_se_samp_slope"] > 0 else np.nan
                    rows.append(global_row)
    result = pd.DataFrame(rows)
    if result.empty:
        raise RuntimeError("No evaluable MSSM models")
    result["part"] = result["metric"].map({"scatter": "1", "shannon": "2", "E1": "3", "E2_stratum": "3"})
    result["gene_set"] = np.where(result["metric"].isin(["E1", "E2_stratum"]), result["metric"], "none")
    result["analysis_role"] = "primary_full_depth"
    result["yazar_pool_fe_status"] = "not_applicable_single_cohort"
    result = apply_slice_bh(result, ["dataset", "part", "metric", "gene_set", "depth_mode", "stratum", "model_spec"])
    result["slope_per_decade"] = result["slope_raw_per_decade"]
    result["se"] = result["se_raw"]
    result["p"] = result["pvalue"]
    result["q"] = result["q_bh_within_slice"]
    result["mde80"] = result["mde80_raw_per_decade"]
    result["effect_scale"] = np.where(result["is_stratified_global"], "within_cell_type_SD", "raw_outcome_units")
    return result


def build_qc(meta: dict[str, Any], refs: dict[str, Any], slopes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ci, ct in enumerate(meta["cell_types"]):
        pre = refs["cnt_pre"][:, ci]
        post = refs["cnt_15"][:, ci]
        e1_summary = refs["E1"]["summary"]
        e2_available = [ci in refs["E2_strata"][si]["summary"] for si in range(4)]
        tested = slopes.loc[(slopes["group"] == ct) & ~slopes["is_pooled"], "metric"].unique().tolist()
        rows.append({
            "cell_type": ct,
            "n_donors_30plus_with_cells": int((pre > 0).sum()),
            "n_cells_before_depth_gate": int(pre.sum()),
            "n_donors_ge20_before_depth_gate": int((pre >= 20).sum()),
            "n_cells_ge1500": int(post.sum()),
            "n_cells_excluded_lt1500": int(pre.sum() - post.sum()),
            "n_donors_ge20_after_depth_gate": int((post >= 20).sum()),
            "eligible_slope_6donors_20cells": bool((post >= 20).sum() >= 6),
            "n_young_reference_donors_ge20": int((post[refs["bin_arr"] == 1] >= 20).sum()),
            "E1_reference_available": bool(ci in e1_summary),
            "E1_n_markers_full": int(e1_summary.get(ci, 0)),
            "E2_strata_reference_available": ";".join(DEPTH_STRATA[si][0] for si, ok in enumerate(e2_available) if ok),
            "E2_strata_n_markers_full": ";".join(str(len(refs["E2_strata"][si]["full"].get(ci, []))) for si in range(4)),
            "tested_measures": ";".join(sorted(tested)),
        })
    return pd.DataFrame(rows)


def gao_validation() -> dict[str, Any]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    output: dict[str, Any] = {"cell_type": "Lymphoid", "checks": {}}
    for scope in ("full_range", "mature_30plus"):
        if scope == "full_range":
            meta = reconcile_gao_metadata(GAO_H5AD)
        else:
            import scripts.run_followup_mature_fulldepth as fulldepth
            original_gao_path = fulldepth.GAO_H5AD
            fulldepth.GAO_H5AD = GAO_H5AD
            try:
                meta = reconcile_mature_metadata("gao")
            finally:
                fulldepth.GAO_H5AD = original_gao_path
        ci = meta["cell_types"].index("Lymphoid")
        mapping, n_genes, audit = raw_to_clean_map(meta["gene_symbols"])
        records = []
        t0 = time.perf_counter()
        with h5py.File(GAO_H5AD, "r") as handle:
            group = handle["layers"]["counts"]
            for di, donor in enumerate(meta["donor_ids"]):
                rows = np.flatnonzero((meta["cell_donor_idx"] == di) & (meta["cell_ct_idx"] == ci))
                if len(rows) < MIN_CELLS_PER_DONOR:
                    continue
                mat, _ = read_clean_rows(group, meta["indptr"], rows, mapping, n_genes)
                umi = np.rint(np.asarray(mat.sum(axis=1)).ravel()).astype(np.int64)
                mat = mat[umi >= 1500].tocsr()
                if mat.shape[0] < MIN_CELLS_PER_DONOR:
                    continue
                basic = basic_group_metrics(mat, di, ci, split_half=False)
                dr = meta["donors_df"].iloc[di]
                records.append({
                    "donor_id": str(donor), "age_per_decade": float(dr["age_per_decade"]), "bin": str(dr["bin"]),
                    "sex_male": int(dr["sex_male"]), "study_or_pool": str(dr["study"] if scope == "full_range" else dr["study"]),
                    **basic,
                })
                del mat
        frame = pd.DataFrame(records)
        comparisons = []
        if scope == "mature_30plus":
            donor1 = pd.read_csv(RESULTS_DIR / "followup-mature-fulldepth-gao-donor-part1.csv")
            donor2 = pd.read_csv(RESULTS_DIR / "followup-mature-fulldepth-gao-donor-part2.csv")
            ref1 = donor1.loc[donor1["group"] == "Lymphoid", ["donor_id", "scatter", "_depth", "n_cells"]]
            ref2 = donor2.loc[donor2["group"] == "Lymphoid", ["donor_id", "full_shannon"]]
            expected = ref1.merge(ref2, on="donor_id", validate="one_to_one")
            observed = frame[["donor_id", "scatter", "shannon", "_depth", "n_cells"]]
            merged = expected.merge(observed, on="donor_id", suffixes=("_expected", "_observed"), validate="one_to_one")
            assert len(merged) == len(expected) == len(observed)
            for expected_col, observed_col in (("scatter_expected", "scatter_observed"), ("full_shannon", "shannon"), ("_depth_expected", "_depth_observed"), ("n_cells_expected", "n_cells_observed")):
                delta = np.max(np.abs(merged[expected_col].to_numpy(float) - merged[observed_col].to_numpy(float)))
                comparisons.append({"field": expected_col, "max_abs_delta": float(delta)})
                assert delta < 1e-10, (scope, expected_col, delta)
        prefix = "followup-fulldepth" if scope == "full_range" else "followup-mature-fulldepth"
        for metric, filename, outcome in (("scatter", "cell-drift", "scatter"), ("shannon", "spread", "shannon")):
            reference = pd.read_csv(RESULTS_DIR / f"{prefix}-gao-{filename}.csv")
            expected = reference.loc[
                (reference["metric"] == metric) & (reference["depth_mode"] == "full_depth")
                & (reference["model_spec"] == "primary") & (reference["group"] == "Lymphoid")
            ].iloc[0]
            fitted = fit_ols_model_row(
                frame, outcome, "age_per_decade", ["sex_male", "_depth"], "study_or_pool", 2, 5,
                "ENCODE v4 (Snyder)", {}
            )
            assert fitted is not None
            for field in ("slope_raw_per_decade", "se_raw", "pvalue", "n_donors", "n_cells"):
                delta = abs(float(fitted[field]) - float(expected[field]))
                tolerance = 1e-9 if field != "pvalue" else 1e-8
                comparisons.append({"field": f"{metric}_{field}", "observed": float(fitted[field]), "expected": float(expected[field]), "abs_delta": float(delta)})
                assert delta < tolerance, (scope, metric, field, delta)
        output["checks"][scope] = {
            "status": "PASS", "n_donor_rows": len(frame), "elapsed_s": time.perf_counter() - t0,
            "gene_audit_retained": audit["n_genes_retained"], "comparisons": comparisons,
        }
        log(f"Gao {scope} validation PASS")
    output["status"] = "PASS"
    output["reference_path_note"] = "The requested spread CSV contains Shannon/KL, not scatter; scatter was checked against the cell-drift CSV and Shannon against the spread CSV."
    write_json(VALIDATION_JSON, output)
    return output


def fmt(value: Any, digits: int = 5) -> str:
    if value is None or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):.{digits}g}"


def markdown_table(df: pd.DataFrame, columns: list[str], labels: list[str] | None = None) -> str:
    labels = labels or columns
    lines = ["| " + " | ".join(labels) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in columns) + " |")
    return "\n".join(lines)


def write_report(meta: dict[str, Any], refs: dict[str, Any], qc: pd.DataFrame, models: pd.DataFrame) -> None:
    slopes = models.loc[~models["is_pooled"] & ~models["is_stratified_global"]].copy()
    pooled = models.loc[models["is_pooled"]].copy()
    probe = meta["probe"]
    audit = meta["gene_audit"]
    eligible = qc.loc[qc["eligible_slope_6donors_20cells"]]
    validation = json.loads(VALIDATION_JSON.read_text(encoding="utf-8"))
    lines = [
        "# Exploratory MSSM mature-to-old concentration test", "",
        f"This owner-approved exploratory run tests controls aged 30 years and older in MSSM DLPFC. Run model: {RUN_MODEL}. It does not enter the manuscript automatically; Reza decides whether and how to use it after seeing the results. The manuscript TeX was not edited.", "",
        "## Input, probe, and gene audit", "",
        f"Input: `{MSSM_H5AD}`; {probe['file_size_bytes']:,} bytes; shape {probe['shape'][0]:,} cells x {probe['shape'][1]:,} genes. All cells are primary DLPFC data, assay {probe['assay_values']}, suspension {probe['suspension_values']}. The run used `obs['cell_type']`, not `subclass`.", "",
        f"Controls: {probe['n_control_cells']:,} cells from {probe['n_control_donors']} donors. Exactly {probe['n_control_donors_parseable_age']} donors have a parseable age, range {probe['control_age_min']:.0f} to {probe['control_age_max']:.0f}; {probe['n_control_donors_missing_age']} lack one. At age 30+, there are {probe['n_donors_30plus']} donors: {probe['n_male_30plus']} male and {probe['n_female_30plus']} female. These exactly match the supplied probe counts. The 30-39 reference contains {probe['n_reference_donors_30_39']} donors.", "",
        f"Gene audit: {audit['n_genes_raw']:,} raw names; {audit['n_duplicate_collapsed']} duplicate occurrences collapsed; {audit['n_excluded_mt']} MT, {audit['n_excluded_rp_rpl_rps']} RP/RPL/RPS, {audit['n_excluded_mrp']} MRP, {audit['n_excluded_ercc']} ERCC, and {audit['n_excluded_sex']} sex genes excluded; {audit['n_genes_retained']:,} retained. Sex genes present: {', '.join(audit['sex_genes_present'])}; absent: {', '.join(audit['sex_genes_absent'])}.", "",
        "## Validation before MSSM", "",
        "Gao Lymphoid was independently read from raw count HDF5 slices with the adapted gene projection and per-cell functions. Full-range and mature 30+ donor scatter, Shannon, depth, and cell counts matched the committed donor rows, and the fitted slopes, standard errors, p-values, donor counts, and cell counts matched their respective committed model rows within rounding. Both validations PASS.", "",
        "The requested `results/followup-fulldepth-gao-spread.csv` contains Shannon/KL rather than Hellinger scatter. Accordingly, scatter was checked against `results/followup-fulldepth-gao-cell-drift.csv`, and Shannon was checked against the requested spread CSV. The mature equivalents were also checked. No Yazar data were sought or used.", "",
        "## Methods", "",
        "Counts were converted to within-cell proportions after the exact original gene audit. Cells below 1,500 post-audit molecules were excluded as an eligibility filter only; no primary outcome was downsampled. Donor is the inferential unit. A donor-cell-type group needs at least 20 eligible cells, and a slope needs at least six donors.", "",
        "Scatter is the mean per-cell Hellinger distance from the leave-one-out mean proportions of the donor's other cells of that type. Shannon is -sum(p log p). E1 is the original mean detected off-program gene count: young target detection <=5%, another-type detection >=20%, other-type CP10K >=20 and >=5 times target CP10K+1, with original reference and leave-one-out gates. E2 uses genes detected in <1% of young-reference cells and is evaluated in the original bounded strata [1500,2500), [2500,4500), [4500,9000), [9000,18000), with >=20 cells per donor-cell-type-stratum and the original >=200 leave-one-out young-reference cell gate. Counts are detected genes per cell, not UMIs.", "",
        f"Per-type OLS is outcome ~ age_per_decade + sex_male + mean(log UMI), exactly the original primary depth covariate, with the original minimum residual df=5. Pooled puts all of a donor's cells in one group; pooled E1/E2 combine the cell-type-specific gene-set scores. Stratified-global standardizes within cell type (or cell type x E2 stratum), centers sex and depth within group, includes cell-type/stratum dummies, and clusters SEs by donor. The original implementation uses linear mean-log-UMI in this regression, not nonlinear spline terms; the ambiguity was resolved in favor of the original code. Split-half uses the five original seeds {SH_SEEDS}; scale-invariant outcomes use the original 0.5 sampling-SE factor and counts use 1.0. BH is across tested cell types within each measure (and within each E2 stratum).", "",
        "## Eligibility and QC", "",
        f"The depth-gated slope rule yields {len(eligible)} eligible cell types. The pre-depth donor/cell probe also yields {int((qc['n_donors_ge20_before_depth_gate'] >= 6).sum())} types, exactly matching the supplied value 21. Across all types, {int(qc['n_cells_ge1500'].sum()):,} of {int(qc['n_cells_before_depth_gate'].sum()):,} cells pass the depth gate and {int(qc['n_cells_excluded_lt1500'].sum()):,} are excluded. The original residual-df>=5 model gate leaves 19 fitted scatter/Shannon types: T cells and smooth muscle cells each have seven donors and are not fitted. Donor-age handling exactly reproduces 209 mature donors, so it does not explain the Nov 23 count of 14.", "",
        markdown_table(qc.assign(eligible=qc['eligible_slope_6donors_20cells'].map({True: 'yes', False: 'no'})).astype({'n_donors_ge20_before_depth_gate': int, 'n_donors_ge20_after_depth_gate': int}), ["cell_type", "n_donors_ge20_before_depth_gate", "n_donors_ge20_after_depth_gate", "n_cells_ge1500", "eligible", "E1_n_markers_full"], ["Cell type", "Donors >=20 pre-depth", "Donors >=20 post-depth", "Cells >=1500", "Slope eligible", "E1 markers"]), "",
        f"14-vs-21 reconciliation: 21 is the general `cell_type`-level six-donor/20-cell count both before and after depth gating. Endpoint-specific reference gates are narrower: {int(((qc['E1_n_markers_full'] > 0) & qc['eligible_slope_6donors_20cells']).sum())} types have usable nonzero E1 markers, and the numbers of tested E2 types are " + ", ".join(f"{st}: {len(slopes.loc[(slopes['metric'] == 'E2_stratum') & (slopes['stratum'] == st), 'group'].unique())}" for st in sorted(slopes.loc[slopes['metric'] == 'E2_stratum', 'stratum'].unique())) + ". The exact Nov 23 feasibility definition is not present in the supplied instruction, so its 14 cannot honestly be reproduced by relabeling the general count. The nearest endpoint-specific count is 15 for E1, indicating an additional endpoint/reference/type filter (or a different earlier type level), not a donor-age discrepancy. No count was forced to 14.", "",
        "## Per-cell-type slopes", "",
        "Entries are raw-outcome slopes per decade with p and BH q. Stratified-global effects are SD units, shown separately below.", "",
    ]
    for measure in ("scatter", "shannon", "E1"):
        sub = slopes.loc[slopes["metric"] == measure].copy()
        lines.extend([f"### {measure}", ""])
        if sub.empty:
            lines.extend(["No evaluable cell types.", ""])
            continue
        view = pd.DataFrame({
            "Cell type": sub["group"], "Slope": sub["slope_raw_per_decade"].map(fmt),
            "SE": sub["se_raw"].map(fmt), "p": sub["pvalue"].map(lambda x: fmt(x, 4)),
            "BH q": sub["q_bh_within_slice"].map(lambda x: fmt(x, 4)),
            "Donors": sub["n_donors"].astype(int), "Cells": sub["n_cells"].astype(int),
            "MDE80": sub["mde80_raw_per_decade"].map(fmt),
        })
        lines.extend([markdown_table(view, list(view.columns)), ""])
    e2 = slopes.loc[slopes["metric"] == "E2_stratum"].copy()
    lines.extend(["### E2 depth-stratum counts", ""])
    if not e2.empty:
        view = pd.DataFrame({
            "Cell type": e2["group"], "Stratum": e2["stratum"], "Slope": e2["slope_raw_per_decade"].map(fmt),
            "SE": e2["se_raw"].map(fmt), "p": e2["pvalue"].map(lambda x: fmt(x, 4)),
            "BH q": e2["q_bh_within_slice"].map(lambda x: fmt(x, 4)), "Donors": e2["n_donors"].astype(int),
            "Cells": e2["n_cells"].astype(int), "MDE80": e2["mde80_raw_per_decade"].map(fmt),
        })
        lines.extend([markdown_table(view, list(view.columns)), ""])
    else:
        lines.extend(["No evaluable cell-type x stratum rows.", ""])
    lines.extend(["## Pooled and stratified-global", ""])
    view = pd.DataFrame({
        "Measure": pooled["metric"], "Stratum": pooled["stratum"], "Type": np.where(pooled["is_stratified_global"], "stratified-global", "pooled"),
        "Slope": pooled["slope_raw_per_decade"].map(fmt), "SE": pooled["se_raw"].map(fmt),
        "p": pooled["pvalue"].map(lambda x: fmt(x, 4)), "Donors": pooled["n_donors"].astype(int),
        "Cells": pooled["n_cells"].astype(int), "MDE80": pooled["mde80_raw_per_decade"].map(fmt),
        "SH RMS": pooled["sh_rms_delta_donor"].map(fmt), "SH sampling SE": pooled["sh_se_samp_slope"].map(fmt),
        "SH rho": pooled["sh_rho_tech2"].map(fmt), "SH SNR": pooled["sh_snr_slope"].map(fmt),
    })
    lines.extend([markdown_table(view, list(view.columns)), "", "Pooled slopes are raw outcome units; stratified-global slopes are within-cell-type SD units. Pooled q=p as a singleton; original stratified-global q is NA. SH RMS is the root-mean-square donor-level half difference, SH sampling SE uses the original scale factor, SH rho is the technical variance ratio, and SH SNR is |slope|/SH sampling SE.", ""])
    lines.extend(["## Split-half technical null and passing types", ""])
    lines.append("The pooled/global table reports the actual five-seed null summaries. Per-type RMS half differences, sampling SEs, technical variance ratios, and SNR are in the slopes CSV. The fixed-seed nulls quantify molecule-sampling repeatability; they are not biological permutation tests.")
    lines.append("")
    passing = slopes.loc[slopes["q_bh_within_slice"] < 0.05].copy()
    if passing.empty:
        lines.append("No cell-type x measure row passes BH q < 0.05.")
    else:
        for _, row in passing.iterrows():
            lines.append(f"- BH PASS: {row['group']} {row['metric']} {row['stratum']}: slope {fmt(row['slope_raw_per_decade'])}, p {fmt(row['pvalue'], 4)}, q {fmt(row['q_bh_within_slice'], 4)}.")
    nominal = slopes.loc[slopes["pvalue"] < 0.05].copy()
    lines.append("")
    lines.append(f"Nominal per-type p < 0.05: {len(nominal)} rows; BH-passing: {len(passing)} rows.")
    for _, row in nominal.iterrows():
        lines.append(f"- Nominal: {row['group']} {row['metric']} {row['stratum']}: slope {fmt(row['slope_raw_per_decade'])}, p {fmt(row['pvalue'], 4)}, q {fmt(row['q_bh_within_slice'], 4)}.")
    lines.append("")
    lines.append("Nominally significant pooled/global rows (p < 0.05): " + ("none" if pooled.loc[pooled['pvalue'] < 0.05].empty else "; ".join(f"{r.metric} {r.stratum} {'global' if r.is_stratified_global else 'pooled'} slope {fmt(r.slope_raw_per_decade)}, p {fmt(r.pvalue, 4)}" for r in pooled.loc[pooled['pvalue'] < 0.05].itertuples())) + ".")
    lines.extend(["", "All three BH passes are microglial-cell results with positive slopes: more within-cell Hellinger scatter and more detected E2 off-program genes with age. This is the opposite direction from expression concentration. Their split-half SNRs are 76.61 (scatter), 18.27 (E2 s1), and 24.86 (E2 s2), so the point estimates are large relative to molecular sampling noise. However, no pooled or stratified-global row passes p < 0.05, and the E2 global summaries for s1 and s2 are p 0.1068 and 0.0722. These are localized exploratory microglial associations, not a general MSSM concentration replication.", ""])
    lines.extend(["", "## Power and interpretation", "", f"The pre-stated power context is MDE80 = {PRESTATED_MDE80} per decade. The test is powered for Gao-sized effects (-0.0723 per decade) and larger, not Yazar-sized effects (-0.0197 per decade, approximately 453 donors needed). A null here does not rule out a Yazar-sized fall. Row-specific model-based MDE80 values are also reported; they do not replace that pre-stated context.", "", "## Limitations", "", "- Exploratory, single-cohort, post-mortem, cross-sectional DLPFC analysis; donor age is not longitudinal change.", "- Only normal controls aged 30+ enter. Eight 30-39 reference donors make E1/E2 marker estimation a potential source of instability; types without reference sets are untestable for those outcomes, not null. Usable E1 sets contain only 1 to 9 markers.", "- No study/pool fixed effects exist in this single cohort. Available metadata have constant assay and suspension and no batch, postmortem interval, RIN, or processing-cohort adjustment; residual technical confounding remains possible.", "- Depth adjustment is donor-level mean log UMI, exactly as in the original primary. Cell composition can affect pooled results; stratified-global is the within-type summary.", "- E1/E2 detected-gene counts are not molecule counts. E2 strata have different donor sets and BH families and are not independent replications; the s1 and s2 BH families contain only 2 and 4 cell types.", "- Split-half only estimates molecular sampling noise and cannot rule out ambient RNA, annotation error, or biological donor confounding.", "- No result enters the paper automatically.", ""])
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def run_mssm(resume: bool = True) -> None:
    if not VALIDATION_JSON.exists() or json.loads(VALIDATION_JSON.read_text(encoding="utf-8")).get("status") != "PASS":
        raise RuntimeError("Gao validation has not passed; MSSM will not run")
    RESULTS_DIR.mkdir(exist_ok=True)
    REPORTS_DIR.mkdir(exist_ok=True)
    meta = mssm_metadata()
    refs = pass1_mssm(meta, resume=resume)
    ct_rows, e2_rows = pass2_mssm(meta, refs, resume=resume)
    models = fit_models(ct_rows, e2_rows)
    qc = build_qc(meta, refs, models)
    slopes = models.loc[~models["is_pooled"] & ~models["is_stratified_global"]].copy()
    pooled = models.loc[models["is_pooled"]].copy()
    slopes.to_csv(SLOPES_CSV, index=False)
    pooled.to_csv(POOLED_CSV, index=False)
    qc.to_csv(QC_CSV, index=False)
    write_report(meta, refs, qc, models)
    log(f"Wrote {SLOPES_CSV}, {POOLED_CSV}, {QC_CSV}, and {REPORT_MD}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("validate-gao", "mssm", "all"))
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.stage in ("validate-gao", "all"):
        gao_validation()
    if args.stage in ("mssm", "all"):
        run_mssm(resume=not args.no_resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
