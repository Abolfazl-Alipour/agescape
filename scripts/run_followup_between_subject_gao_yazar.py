#!/usr/bin/env python3
"""Between-subject aging drift/variance/dispersion (Part A) and within-cell spread (Part B) in Gao Heart and Yazar Blood."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gc
import json
import os
from pathlib import Path
import pickle
import resource
import sys
import time
from typing import Any

for _env_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_env_key, "1")

import h5py
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats as scipy_stats
from scipy.special import digamma
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_cross_species_expression_vector import reference_without_donor, technical_gene_mask
from scripts.run_followup_de_man_entropy import SEEDS as ENTROPY_SEEDS
from scripts.run_followup_de_man_entropy import kl_uniform, multinomial_downsample, shannon
from scripts.run_followup_hcata_pool_sensitivity import add_yazar_pool_dummies
from scripts.run_followup_organ_17_heart_gao import STUDY_DUMMY_MAP, build_donor_audit
from scripts.run_h3_expression_vector import hellinger_distance
from scripts.run_rev3_randomness_variability_test import (
    _df_logvar_bias,
    _donor_pseudobulk_and_dispersion,
    _rarefy_csr,
    _residualize_against_depth,
    _stratum_stats,
)
from src.followup_muscle import count_gate
from src.followup_organs import ERCC_RE, assert_no_forbidden_dashes, write_clean_csv, write_clean_markdown

DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
GAO_H5AD = DATA_ROOT / "organ17_gao_heart" / "GSE290367_integrated_RNA_adata.h5ad"
YAZAR_H5AD = DATA_ROOT / "hcata" / "yazar_1e44db10.h5ad"
YAZAR_CACHE = DATA_ROOT / "hcata" / "yazar_cache"
CHECKPOINT_DIR = DATA_ROOT / "between_subject_cache"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"

MASTER_SEED = 42
A2_N_PERM = 2000
CHUNK_ROWS = 16000
SUBCHUNK_ROWS = 1000
MAX_WORKERS = 6
MIN_CELLS_PER_DONOR = 20
MIN_DONORS_SLOPE = 6
MIN_YOUNG_A1 = 2
MIN_LOO_REF_CELLS = 50
MIN_DONORS_A2_BIN = 6
DEPTH_LEVELS = (500, 1500, 3000)
YAZAR_DOUBLET_LABEL = "peripheral blood mononuclear cell"

GAO_EXPECTED_CELL_TYPES = (
    "Adipocyte",
    "Cardiomyocyte",
    "Endocardial",
    "Endothelial",
    "Epicardial",
    "Fibroblast",
    "LEC",
    "Lymphoid",
    "Mast",
    "Myeloid",
    "Neuronal",
    "Pericyte",
    "vSMC",
)

# Worker globals populated only in worker initializer before processing chunks
_W_GPROJ: sp.csc_matrix | None = None
_W_REFS: dict[str, Any] = {}


def _worker_init(g_proj: sp.csc_matrix | None = None, refs: dict[str, Any] | None = None) -> None:
    """Store read-only per-worker matrices inside spawned/forked worker globals and freeze CoW heap."""
    global _W_GPROJ, _W_REFS
    _W_GPROJ = g_proj
    _W_REFS = refs or {}
    gc.freeze()


def assign_age_bin(age_years: float) -> str:
    """Assign pre-specified cross-sectional age bins: <20=LT20, 20-39=YOUNG, 40-59=MID, 60+=OLD."""
    if age_years < 20.0:
        return "LT20"
    if age_years <= 39.0:
        return "YOUNG"
    if age_years <= 59.0:
        return "MID"
    return "OLD"


def _df_logvar_bias_from_df(df_val: float) -> float:
    """Exact Gaussian E[log(chi2_df / df)] = psi(df/2) - log(df/2) for arbitrary positive trace df."""
    if df_val <= 1e-6:
        raise ValueError(f"Degrees of freedom must be positive, got {df_val}")
    return float(digamma(0.5 * df_val) - np.log(0.5 * df_val))


def build_human_gene_projection(raw_symbols: np.ndarray | pd.Index | list[str]) -> dict[str, Any]:
    """Canonicalize symbols, collapse duplicates by column sum, and exclude human technical/sex/ERCC genes."""
    raw_arr = np.asarray([str(x).strip() for x in raw_symbols], dtype=object)
    n_raw = len(raw_arr)
    low = np.array([x.lower() for x in raw_arr], dtype=object)
    valid_mask = (raw_arr != "") & (low != "nan") & (low != "na")

    unique_ordered: list[str] = []
    sym_to_uidx: dict[str, int] = {}
    valid_raw_indices: list[int] = []
    valid_target_indices: list[int] = []

    for r_i in np.flatnonzero(valid_mask):
        sym = str(raw_arr[r_i])
        if sym not in sym_to_uidx:
            sym_to_uidx[sym] = len(unique_ordered)
            unique_ordered.append(sym)
        valid_raw_indices.append(int(r_i))
        valid_target_indices.append(sym_to_uidx[sym])

    unique_arr = np.array(unique_ordered, dtype=object)
    tech_base = technical_gene_mask(unique_arr)
    ercc_mask = np.array([bool(ERCC_RE.match(str(g))) for g in unique_arr], dtype=bool)
    excl_mask = tech_base | ercc_mask
    keep_unique_mask = ~excl_mask

    retained_symbols = unique_arr[keep_unique_mask]
    old_to_new = np.full(len(unique_arr), -1, dtype=np.int64)
    old_to_new[keep_unique_mask] = np.arange(int(keep_unique_mask.sum()), dtype=np.int64)

    rows_g: list[int] = []
    cols_g: list[int] = []
    for r_i, u_i in zip(valid_raw_indices, valid_target_indices):
        new_j = int(old_to_new[u_i])
        if new_j >= 0:
            rows_g.append(r_i)
            cols_g.append(new_j)

    proj = sp.csc_matrix(
        (np.ones(len(rows_g), dtype=np.float64), (rows_g, cols_g)),
        shape=(n_raw, int(len(retained_symbols))),
        dtype=np.float64,
    )

    upper_u = np.array([str(g).upper() for g in unique_arr], dtype=object)
    sex_targets = ("XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D")
    present_sex = [g for g in sex_targets if g in set(upper_u.tolist())]
    absent_sex = [g for g in sex_targets if g not in set(upper_u.tolist())]
    audit = {
        "n_genes_raw": int(n_raw),
        "n_invalid_names": int((~valid_mask).sum()),
        "n_duplicate_collapsed": int(valid_mask.sum() - len(unique_arr)),
        "n_genes_unique": int(len(unique_arr)),
        "n_excluded_mt": int(np.sum([str(g).startswith("MT-") for g in upper_u])),
        "n_excluded_rp_rpl_rps": int(np.sum([str(g).startswith(("RPL", "RPS")) for g in upper_u])),
        "n_excluded_mrp": int(np.sum([str(g).startswith(("MRPL", "MRPS")) for g in upper_u])),
        "n_excluded_ercc": int(ercc_mask.sum()),
        "n_excluded_sex": int(len(present_sex)),
        "sex_genes_present": present_sex,
        "sex_genes_absent": absent_sex,
        "n_technical_excluded": int(excl_mask.sum()),
        "n_genes_retained": int(len(retained_symbols)),
        "retained_symbols": [str(x) for x in retained_symbols.tolist()],
    }
    return {"proj": proj, "audit": audit}


def read_obs_column(obs_group: h5py.Group, column: str) -> np.ndarray:
    """Read an HDF5 anndata .obs column cleanly into string numpy array without full-obs load."""
    item = obs_group[column]
    if isinstance(item, h5py.Group) and "categories" in item:
        raw_cats = item["categories"][:]
        cats = np.array(
            [c.decode("utf-8", "replace") if isinstance(c, (bytes, np.bytes_)) else str(c) for c in raw_cats],
            dtype=object,
        )
        codes = item["codes"][:]
        vals = cats[np.maximum(codes, 0)]
        vals[codes < 0] = "nan"
        return vals
    raw_vals = item[:]
    return np.array(
        [v.decode("utf-8", "replace") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in raw_vals],
        dtype=object,
    )


def _decode_str_array(arr: np.ndarray) -> np.ndarray:
    return np.array(
        [v.decode("utf-8", "replace") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in arr],
        dtype=object,
    )


def read_gao_var_symbols(handle: h5py.File) -> np.ndarray:
    """Read gene symbols from Gao GSE290367 HDF5 var group."""
    var_grp = handle["var"]
    idx_key = var_grp.attrs.get("_index", "_index")
    if isinstance(idx_key, bytes):
        idx_key = idx_key.decode("utf-8")
    return _decode_str_array(var_grp[str(idx_key)][:])


def read_yazar_var_and_slot(handle: h5py.File) -> tuple[str, np.ndarray]:
    """Identify raw matrix group ('raw/X' if present else 'X') and feature_name/index gene symbols."""
    if "raw" in handle and "X" in handle["raw"] and "var" in handle["raw"]:
        slot = "raw/X"
        vgrp = handle["raw/var"]
    else:
        slot = "X"
        vgrp = handle["var"]
    if "feature_name" in vgrp:
        return slot, read_obs_column(vgrp, "feature_name")
    idx_key = vgrp.attrs.get("_index", "_index")
    if isinstance(idx_key, bytes):
        idx_key = idx_key.decode("utf-8")
    return slot, _decode_str_array(vgrp[str(idx_key)][:])


def read_csr_chunk(group: h5py.Group, start: int, stop: int, n_vars: int, indptr_all: np.ndarray) -> sp.csr_matrix:
    """Read CSR slice [start:stop] cleanly using pre-loaded indptr."""
    p0, p1 = int(indptr_all[start]), int(indptr_all[stop])
    data = group["data"][p0:p1]
    indices = group["indices"][p0:p1]
    return sp.csr_matrix((data, indices, indptr_all[start : stop + 1] - p0), shape=(stop - start, n_vars))


def build_stratified_split_half(donor_df: pd.DataFrame, stratify_cols: tuple[str, ...], seed: int = MASTER_SEED) -> dict[str, int]:
    """Deterministically assign YOUNG (20-39) donors 50/50 into Half 1 (ref=1) and Half 2 (eval=2)."""
    rng = np.random.default_rng(seed)
    young = donor_df.loc[donor_df["bin"] == "YOUNG"].reset_index(drop=True).sort_values(["donor_id"]).copy()
    assignment: dict[str, int] = {}
    if young.empty:
        return assignment

    # Hierarchical stratification fallback until >=2 donors per stratum where possible
    key_col = None
    for k in range(len(stratify_cols), 0, -1):
        cols = list(stratify_cols[:k])
        labels = young[cols].astype(str).agg("|".join, axis=1)
        if int((labels.value_counts() >= 2).sum()) >= 1:
            key_col = labels
            break
    if key_col is None:
        key_col = pd.Series(["global"] * len(young), index=young.index)

    young["_strat"] = key_col
    h1_total, h2_total = 0, 0
    for strat_val in sorted(young["_strat"].unique()):
        ids = young.loc[young["_strat"] == strat_val, "donor_id"].tolist()
        perm = list(rng.permutation(ids))
        for d in perm:
            if h1_total <= h2_total:
                assignment[str(d)] = 1
                h1_total += 1
            else:
                assignment[str(d)] = 2
                h2_total += 1
    return assignment


def reconcile_gao_metadata(h5_path: Path) -> dict[str, Any]:
    """Load Gao HDF5 obs, verify against committed donor CSV, and build donor/cell mappings."""
    with h5py.File(h5_path, "r") as f:
        obs_grp = f["obs"]
        cols = ["donor_id", "study", "age", "age_status", "disease", "sex", "technology", "region", "final_cell_type"]
        obs = pd.DataFrame({c: read_obs_column(obs_grp, c) for c in cols})
        gene_symbols = read_gao_var_symbols(f)
        indptr = f["layers"]["counts"]["indptr"][:].astype(np.int64)

    audit_df = build_donor_audit(obs)
    inc = audit_df[audit_df["included_primary"]].copy().sort_values(["age_num", "donor_id"])
    ref_csv = pd.read_csv(RESULTS_DIR / "followup-organ-17-heart-gao-donors.csv")
    ref_inc = ref_csv[ref_csv["included_primary"] == True].copy()  # noqa: E712
    assert set(inc["donor_id"]) == set(ref_inc["donor_id"]), "Gao included donors differ from committed CSV"

    inc["age_years"] = inc["age_num"].astype(float)
    inc["age_per_decade"] = inc["age_years"] / 10.0
    inc["bin"] = inc["age_years"].apply(assign_age_bin)
    inc["sex_male"] = (inc["sex"] == "male").astype(int)
    for st_name, col_name in STUDY_DUMMY_MAP.items():
        inc[col_name] = (inc["study"] == st_name).astype(float)

    bin_counts = inc["bin"].value_counts().to_dict()
    assert len(inc) == 117
    assert bin_counts.get("YOUNG", 0) == 20 and bin_counts.get("MID", 0) == 53 and bin_counts.get("OLD", 0) == 44 and bin_counts.get("LT20", 0) == 0

    split_half = build_stratified_split_half(inc, ("sex", "study"), seed=MASTER_SEED)
    inc["split_half_group"] = inc["donor_id"].map(lambda d: split_half.get(str(d), 0)).astype(int)

    donor_ids = sorted(inc["donor_id"].astype(str).tolist())
    donor_map = {d: i for i, d in enumerate(donor_ids)}
    cts = sorted(obs.loc[obs["donor_id"].isin(donor_map), "final_cell_type"].astype(str).unique().tolist())
    assert tuple(cts) == GAO_EXPECTED_CELL_TYPES, f"Unexpected Gao cell types: {cts}"
    ct_map = {c: i for i, c in enumerate(cts)}

    obs_donor_arr = obs["donor_id"].astype(str).to_numpy()
    obs_ct_arr = obs["final_cell_type"].astype(str).to_numpy()
    cell_donor_idx = np.array([donor_map.get(d, -1) for d in obs_donor_arr], dtype=np.int32)
    keep_mask = cell_donor_idx >= 0
    cell_ct_idx = np.array([ct_map.get(c, -1) if k else -1 for c, k in zip(obs_ct_arr, keep_mask)], dtype=np.int32)
    assert int(keep_mask.sum()) == 901946

    return {
        "dataset": "gao",
        "h5_path": h5_path,
        "matrix_slot": "layers/counts",
        "gene_symbols": gene_symbols,
        "indptr": indptr,
        "n_obs": len(obs),
        "keep_mask": keep_mask,
        "cell_donor_idx": cell_donor_idx,
        "cell_ct_idx": cell_ct_idx,
        "donors_df": inc.set_index(pd.Index(inc["donor_id"].astype(str), name=None)).loc[donor_ids].reset_index(drop=True),
        "donor_ids": donor_ids,
        "cell_types": cts,
        "pooled_name": "All_LV_Nuclei",
        "reconcile": {
            "n_donors": int(len(inc)),
            "n_cells_kept": int(keep_mask.sum()),
            "bins": {k: int(v) for k, v in bin_counts.items()},
            "age_min": float(inc["age_years"].min()),
            "age_max": float(inc["age_years"].max()),
        },
    }


def reconcile_yazar_metadata(h5_path: Path, cache_dir: Path) -> dict[str, Any]:
    """Load Yazar HDF5 obs and yazar_cache, verify against donor/eligibility CSVs, and exclude doublets."""
    audit_json = cache_dir / "obs_audit.json"
    pb_manifest_json = cache_dir / "pb_manifest.json"
    donor_meta_csv = cache_dir / "donor_meta.csv"
    cell_counts_csv = cache_dir / "cell_counts.csv"
    for req in (audit_json, pb_manifest_json, donor_meta_csv, cell_counts_csv):
        if not req.exists():
            raise FileNotFoundError(f"Missing reconstructed yazar_cache file: {req}")

    pb_manifest = json.loads(pb_manifest_json.read_text(encoding="utf-8"))
    for grp, fname in pb_manifest.items():
        if not (cache_dir / fname).exists():
            raise FileNotFoundError(f"Missing cache pickle {fname} for {grp}")

    with h5py.File(h5_path, "r") as f:
        obs_grp = f["obs"]
        cols = ["donor_id", "cell_type", "disease", "assay", "suspension_type", "development_stage", "sex", "pool_number"]
        obs = pd.DataFrame({c: read_obs_column(obs_grp, c) for c in cols})
        slot, gene_symbols = read_yazar_var_and_slot(f)
        indptr = (f["raw"]["X"]["indptr"][:] if slot == "raw/X" else f["X"]["indptr"][:]).astype(np.int64)

    base_mask = (
        (obs["disease"].astype(str) == "normal")
        & (obs["assay"].astype(str) == "10x 3' v2")
        & (obs["suspension_type"].astype(str) == "cell")
        & (obs["development_stage"].astype(str).str.endswith("-year-old stage"))
    )
    n_doublets = int((base_mask & (obs["cell_type"].astype(str) == YAZAR_DOUBLET_LABEL)).sum())
    keep_mask = (base_mask & (obs["cell_type"].astype(str) != YAZAR_DOUBLET_LABEL)).to_numpy(dtype=bool)

    yazar_csv = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-donors.csv")
    yazar_csv["donor_id"] = yazar_csv["donor_id"].astype(str)
    yazar_csv = yazar_csv.sort_values(["age_years", "donor_id"]).reset_index(drop=True)
    cache_donors = pd.read_csv(donor_meta_csv)
    cache_donors["donor_id"] = cache_donors["donor_id"].astype(str)
    assert set(cache_donors["donor_id"]) == set(yazar_csv["donor_id"]), "Cache donors mismatch CSV"

    yazar_csv["bin"] = yazar_csv["age_years"].astype(float).apply(assign_age_bin)
    yazar_csv["pool_number"] = yazar_csv["pool_number"].astype(int)
    if "age_per_decade" not in yazar_csv.columns:
        yazar_csv["age_per_decade"] = yazar_csv["age_years"].astype(float) / 10.0
    if "sex_male" not in yazar_csv.columns:
        yazar_csv["sex_male"] = (yazar_csv["sex"].astype(str).str.lower() == "male").astype(int)
    bin_counts = yazar_csv["bin"].value_counts().to_dict()
    assert len(yazar_csv) == 981
    assert (
        bin_counts.get("YOUNG", 0) == 103
        and bin_counts.get("MID", 0) == 193
        and bin_counts.get("OLD", 0) == 680
        and bin_counts.get("LT20", 0) == 5
    )

    yazar_csv["pool_bin"] = pd.qcut(yazar_csv["pool_number"], q=6, labels=False, duplicates="drop").astype(str)
    split_half = build_stratified_split_half(yazar_csv, ("sex", "pool_bin"), seed=MASTER_SEED)
    yazar_csv["split_half_group"] = yazar_csv["donor_id"].map(lambda d: split_half.get(str(d), 0)).astype(int)

    donor_ids = sorted(yazar_csv["donor_id"].tolist())
    donor_map = {d: i for i, d in enumerate(donor_ids)}
    cts = sorted(obs.loc[keep_mask, "cell_type"].astype(str).unique().tolist())
    ct_map = {c: i for i, c in enumerate(cts)}

    obs_d = obs["donor_id"].astype(str).to_numpy()
    obs_c = obs["cell_type"].astype(str).to_numpy()
    cell_donor_idx = np.array([donor_map.get(d, -1) if k else -1 for d, k in zip(obs_d, keep_mask)], dtype=np.int32)
    cell_ct_idx = np.array([ct_map.get(c, -1) if k else -1 for c, k in zip(obs_c, keep_mask)], dtype=np.int32)

    elig_csv = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-pool-fe-eligibility.csv")
    not_eval_set = set(elig_csv.loc[elig_csv["status"] == "not_evaluable_within_pool_df", "group"].astype(str).tolist())

    return {
        "dataset": "yazar",
        "h5_path": h5_path,
        "matrix_slot": slot,
        "gene_symbols": gene_symbols,
        "indptr": indptr,
        "n_obs": len(obs),
        "keep_mask": keep_mask,
        "cell_donor_idx": cell_donor_idx,
        "cell_ct_idx": cell_ct_idx,
        "donors_df": yazar_csv.set_index(pd.Index(yazar_csv["donor_id"].astype(str), name=None)).loc[donor_ids].reset_index(drop=True),
        "donor_ids": donor_ids,
        "cell_types": cts,
        "pooled_name": "All_Cells",
        "not_evaluable_all_umi_csv": sorted(not_eval_set),
        "reconcile": {
            "n_donors": int(len(yazar_csv)),
            "n_cells_base": int(base_mask.sum()),
            "n_doublets_excluded": int(n_doublets),
            "n_cells_kept": int(keep_mask.sum()),
            "bins": {k: int(v) for k, v in bin_counts.items()},
            "age_min": float(yazar_csv["age_years"].min()),
            "age_max": float(yazar_csv["age_years"].max()),
            "n_cell_types_non_doublet": int(len(cts)),
            "yazar_cache_groups": int(len(pb_manifest)),
            "yazar_cache_files_verified": sorted(pb_manifest.values()),
        },
    }


def _subchunk_pass0_worker(task: tuple[Any, ...]) -> dict[str, Any]:
    """Run count_gate, gene projection, and post-exclusion UMI sums on one 2,000-cell subchunk."""
    data, indices, indptr, n_raw_vars, d_idx, c_idx, n_donors, n_cts = task
    gate = count_gate(data)
    if gate["nonfinite"] or gate["negative"] or gate["noninteger"]:
        raise ValueError(f"Count gate failed in chunk: {gate}")

    sub = sp.csr_matrix((data.astype(np.float64, copy=False), indices, indptr), shape=(len(d_idx), n_raw_vars))
    sub_clean = (sub @ _W_GPROJ).tocsr()
    umi = np.rint(np.asarray(sub_clean.sum(axis=1)).ravel()).astype(np.int64)

    n_pairs = n_donors * n_cts
    p_idx = d_idx * n_cts + c_idx
    cnt_mat = np.zeros((n_pairs, 4), dtype=np.int64)
    pool_mat = np.zeros((n_donors, 4), dtype=np.int64)
    log_sum_ct = np.zeros((n_pairs, 4), dtype=np.float64)
    log_sum_pool = np.zeros((n_donors, 4), dtype=np.float64)

    thresholds = (1, 500, 1500, 3000)
    for k, thr in enumerate(thresholds):
        m = umi >= thr
        if not np.any(m):
            continue
        np.add.at(cnt_mat[:, k], p_idx[m], 1)
        np.add.at(pool_mat[:, k], d_idx[m], 1)
        lu = np.log(umi[m].astype(np.float64))
        np.add.at(log_sum_ct[:, k], p_idx[m], lu)
        np.add.at(log_sum_pool[:, k], d_idx[m], lu)

    del sub, sub_clean
    return {
        "gate": gate,
        "nnz": int(len(data)),
        "n_cells": int(len(d_idx)),
        "cnt_mat": cnt_mat,
        "pool_mat": pool_mat,
        "log_sum_ct": log_sum_ct,
        "log_sum_pool": log_sum_pool,
    }


def _subchunk_pass1_worker(task: tuple[Any, ...]) -> dict[str, Any]:
    """Accumulate Pass-1 Part B entropy/KL draws, A1 ref sums, and Gao/Yazar pseudobulks."""
    (
        data,
        indices,
        indptr,
        n_raw_vars,
        d_idx,
        c_idx,
        ds_id_code,
        subchunk_id,
        n_donors,
        n_cts,
        ref_pair_ym_lookup,
        ref_pair_y_lookup,
        split_half_arr,
        bin_arr,
        entropy_seeds,
        do_500_multiseed,
        gao_a3_pair_lookup,
        a2_pb_mode,
        partb_full_floor,
    ) = task

    sub = sp.csr_matrix((data.astype(np.float64, copy=False), indices, indptr), shape=(len(d_idx), n_raw_vars))
    sub_clean = (sub @ _W_GPROJ).tocsr()
    del sub
    sub_clean.eliminate_zeros()
    umi = np.rint(np.asarray(sub_clean.sum(axis=1)).ravel()).astype(np.int64)
    p_idx = d_idx.astype(np.int64) * n_cts + c_idx.astype(np.int64)

    # Compact per-subchunk Part B accumulators mapped over unique local pair/donor IDs (saves ~25x IPC on Yazar)
    u_pids, inv_pids = np.unique(p_idx, return_inverse=True)
    u_dids, inv_dids = np.unique(d_idx.astype(np.int64), return_inverse=True)
    n_lp = len(u_pids)
    n_ld = len(u_dids)
    n_s = len(entropy_seeds)
    b_res: dict[str, Any] = {
        "u_pids": u_pids,
        "u_dids": u_dids,
        "ct_full": np.zeros((n_lp, 3), dtype=np.float64),
        "pool_full": np.zeros((n_ld, 3), dtype=np.float64),
        "ct_full_gt0": np.zeros((n_lp, 3), dtype=np.float64),
        "pool_full_gt0": np.zeros((n_ld, 3), dtype=np.float64),
        "ct_500": np.zeros((n_lp, 3), dtype=np.float64),
        "pool_500": np.zeros((n_ld, 3), dtype=np.float64),
        "ct_3000": np.zeros((n_lp, 3), dtype=np.float64),
        "pool_3000": np.zeros((n_ld, 3), dtype=np.float64),
        "ct_1500_seeds": np.zeros((n_s, n_lp, 3), dtype=np.float64),
        "pool_1500_seeds": np.zeros((n_s, n_ld, 3), dtype=np.float64),
    }
    if do_500_multiseed:
        b_res["ct_500_seeds"] = np.zeros((n_s, n_lp, 3), dtype=np.float64)
        b_res["pool_500_seeds"] = np.zeros((n_s, n_ld, 3), dtype=np.float64)

    def _accum_metrics(rows: np.ndarray, sh_v: np.ndarray | None, klm_v: np.ndarray, kls_v: np.ndarray, ct_arr_out: np.ndarray, pl_arr_out: np.ndarray) -> None:
        if len(rows) == 0:
            return
        rp = inv_pids[rows]
        rd = inv_dids[rows]
        if sh_v is not None:
            np.add.at(ct_arr_out[:, 0], rp, sh_v)
            np.add.at(pl_arr_out[:, 0], rd, sh_v)
        np.add.at(ct_arr_out[:, 1], rp, klm_v)
        np.add.at(pl_arr_out[:, 1], rd, klm_v)
        np.add.at(ct_arr_out[:, 2], rp, kls_v)
        np.add.at(pl_arr_out[:, 2], rd, kls_v)

    m_gt0 = np.flatnonzero(umi > 0)
    if len(m_gt0):
        mat_gt0 = sub_clean[m_gt0]
        _accum_metrics(m_gt0, shannon(mat_gt0), kl_uniform(mat_gt0, sqrt=False), kl_uniform(mat_gt0, sqrt=True), b_res["ct_full_gt0"], b_res["pool_full_gt0"])
        if partb_full_floor == 1:
            b_res["ct_full"] = b_res["ct_full_gt0"].copy()
            b_res["pool_full"] = b_res["pool_full_gt0"].copy()

    if partb_full_floor > 1:
        m_full = np.flatnonzero(umi >= partb_full_floor)
        if len(m_full):
            mat_fl = sub_clean[m_full]
            _accum_metrics(m_full, shannon(mat_fl), kl_uniform(mat_fl, sqrt=False), kl_uniform(mat_fl, sqrt=True), b_res["ct_full"], b_res["pool_full"])

    # 3000 UMI seed 0 (reused by A1 3000 reference below)
    m_3000 = np.flatnonzero(umi >= 3000)
    ds_3000 = None
    if len(m_3000):
        rng_3k = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 3000, subchunk_id, 0]))
        ds_3000 = multinomial_downsample(sub_clean[m_3000], 3000, rng_3k)
        _accum_metrics(m_3000, shannon(ds_3000), kl_uniform(ds_3000, sqrt=False), kl_uniform(ds_3000, sqrt=True), b_res["ct_3000"], b_res["pool_3000"])

    # 1500 UMI 20 seeds (seed 0 also feeds A1 1500 ref and 500-on-1500 A1 ref)
    m_1500 = np.flatnonzero(umi >= 1500)
    ds_1500_s0 = None
    if len(m_1500):
        mat_15 = sub_clean[m_1500]
        for si in range(n_s):
            rng_15 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 1500, subchunk_id, si]))
            ds = multinomial_downsample(mat_15, 1500, rng_15)
            sh_v = shannon(ds) if si < 5 else None
            klm_v = kl_uniform(ds, sqrt=False)
            kls_v = kl_uniform(ds, sqrt=True)
            _accum_metrics(m_1500, sh_v, klm_v, kls_v, b_res["ct_1500_seeds"][si], b_res["pool_1500_seeds"][si])
            if si == 0:
                ds_1500_s0 = ds
            else:
                del ds

    # 500 UMI (for Part B: all >=500 cells; for A1: >=1500 cells thinned to 500)
    m_500 = np.flatnonzero(umi >= 500)
    if len(m_500):
        mat_500 = sub_clean[m_500]
        seed_count_500 = n_s if do_500_multiseed else 1
        for si in range(seed_count_500):
            rng_5 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 500, subchunk_id, si]))
            ds5 = multinomial_downsample(mat_500, 500, rng_5)
            sh_v = shannon(ds5) if si < 5 else None
            klm_v = kl_uniform(ds5, sqrt=False)
            kls_v = kl_uniform(ds5, sqrt=True)
            if si == 0:
                _accum_metrics(m_500, sh_v, klm_v, kls_v, b_res["ct_500"], b_res["pool_500"])
            if do_500_multiseed:
                _accum_metrics(m_500, sh_v, klm_v, kls_v, b_res["ct_500_seeds"][si], b_res["pool_500_seeds"][si])
            del ds5

    # Compact sparse indicator summation helper returning (u_m, data, indices, indptr)
    def _sparse_pair_csr(
        mat: sp.csr_matrix,
        pair_ids: np.ndarray,
        valid_lookup: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
        mapped = pair_ids if valid_lookup is None else valid_lookup[pair_ids]
        sel = mapped >= 0
        if not np.any(sel):
            return None
        sub_m = mat[sel]
        u_m, inv = np.unique(mapped[sel], return_inverse=True)
        ind = sp.csr_matrix((np.ones(len(inv), dtype=np.float64), (inv, np.arange(len(inv)))), shape=(len(u_m), len(inv)))
        summed = (ind @ sub_m).tocsr()
        summed.sum_duplicates()
        return (u_m.astype(np.int64, copy=False), summed.data, summed.indices, summed.indptr)

    # A1 reference sums (YOUNG=1 & MID=2 for full_1500/1500; YOUNG=1 only for 500_same/3000/full_500)
    a1_ref_chunks: dict[str, Any] = {}
    ym_lk = np.where(np.isin(bin_arr, [1, 2]), np.arange(n_donors), -1)
    y_lk = np.where(bin_arr == 1, np.arange(n_donors), -1)
    if len(m_1500):
        pids_15 = p_idx[m_1500]
        dids_15 = d_idx[m_1500]

        a1_ref_chunks["ct_full_1500"] = _sparse_pair_csr(sub_clean[m_1500], pids_15, ref_pair_ym_lookup)
        a1_ref_chunks["ct_1500"] = _sparse_pair_csr(ds_1500_s0, pids_15, ref_pair_ym_lookup)

        rng_5_same = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 5001500, subchunk_id, 0]))
        ds_500_same = multinomial_downsample(sub_clean[m_1500], 500, rng_5_same)
        a1_ref_chunks["ct_500_same"] = _sparse_pair_csr(ds_500_same, pids_15, ref_pair_y_lookup)

        a1_ref_chunks["pool_full_1500"] = _sparse_pair_csr(sub_clean[m_1500], dids_15, ym_lk)
        a1_ref_chunks["pool_1500"] = _sparse_pair_csr(ds_1500_s0, dids_15, ym_lk)
        a1_ref_chunks["pool_500_same"] = _sparse_pair_csr(ds_500_same, dids_15, y_lk)
        del ds_500_same, ds_1500_s0

    if ds_3000 is not None:
        pids_3k = p_idx[m_3000]
        dids_3k = d_idx[m_3000]
        a1_ref_chunks["ct_3000"] = _sparse_pair_csr(ds_3000, pids_3k, ref_pair_y_lookup)
        a1_ref_chunks["pool_3000"] = _sparse_pair_csr(ds_3000, dids_3k, y_lk)
        del ds_3000

    # Extra Gao A1 >=500 ref if needed (YOUNG only)
    if ds_id_code == 1 and len(m_500):
        pids_5 = p_idx[m_500]
        dids_5 = d_idx[m_500]
        a1_ref_chunks["ct_full_500"] = _sparse_pair_csr(sub_clean[m_500], pids_5, ref_pair_y_lookup)
        a1_ref_chunks["pool_full_500"] = _sparse_pair_csr(sub_clean[m_500], dids_5, y_lk)

    # A2 pseudobulk accumulators
    a2_chunks: dict[str, Any] = {}
    if a2_pb_mode == "gao_all_pairs":
        a2_chunks["ct_all"] = _sparse_pair_csr(sub_clean, p_idx)
        a2_chunks["pool_all"] = _sparse_pair_csr(sub_clean, d_idx)
    elif a2_pb_mode == "yazar_pool_only":
        a2_chunks["pool_all"] = _sparse_pair_csr(sub_clean, d_idx)

    # Gao A3 centroids in Pass 1 (saves an entire pass for Gao!)
    a3_chunks: dict[str, Any] = {}
    if gao_a3_pair_lookup is not None and len(m_1500):
        rng_a3 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 31500, subchunk_id, 0]))
        mat_15 = sub_clean[m_1500]
        r_a3 = _rarefy_csr(mat_15, 1500, rng_a3)
        a3_chunks["ct_1500"] = _sparse_pair_csr(r_a3, p_idx[m_1500], gao_a3_pair_lookup)
        a3_chunks["pool_1500"] = _sparse_pair_csr(r_a3, d_idx[m_1500])
        del r_a3

        inv_tot = 1.0 / np.maximum(umi[m_1500].astype(np.float64), 1.0)
        mat_l1 = mat_15.multiply(inv_tot[:, None]).tocsr()
        a3_chunks["ct_full_l1"] = _sparse_pair_csr(mat_l1, p_idx[m_1500], gao_a3_pair_lookup)
        a3_chunks["pool_full_l1"] = _sparse_pair_csr(mat_l1, d_idx[m_1500])
        del mat_l1

    del sub_clean
    return {"partb": b_res, "a1_refs": a1_ref_chunks, "a2": a2_chunks, "a3": a3_chunks}


def _compute_a3_loo_distances_for_matrix(
    sub_mat: sp.csr_matrix,
    pair_ids: np.ndarray,
    centroid_sums: np.ndarray,
    donor_cell_counts: np.ndarray,
    lookup: np.ndarray,
    is_l1_already: bool = False,
) -> dict[int, tuple[float, int]]:
    """Compute exact LOO cell-to-centroid Hellinger distance sums matching _donor_pseudobulk_and_dispersion."""
    mapped = lookup[pair_ids]
    sel = mapped >= 0
    out: dict[int, tuple[float, int]] = {}
    if not np.any(sel):
        return out
    m = sub_mat[sel].tocsr()
    m_ids = mapped[sel]
    row_sums = np.asarray(m.sum(axis=1)).ravel()
    if not is_l1_already:
        row_sums = np.maximum(row_sums, 1.0)
    # Process grouped by mapped pair so 1 vector lookup per donor-group in chunk
    for g_id in np.unique(m_ids):
        n_r = int(donor_cell_counts[g_id])
        if n_r < 2:
            continue
        rm = m_ids == g_id
        sub_g = m[rm]
        tot_g = centroid_sums[g_id]  # probability sum across all n_r cells (sums to n_r)
        rs_g = 1.0 if is_l1_already else row_sums[rm]
        # Build per-cell p_nz and LOO centroid (tot_g - p_i)/(n_r - 1) identically to _donor_pseudobulk_and_dispersion:124-133
        h_sum = 0.0
        for row_i in range(sub_g.shape[0]):
            s0, s1 = sub_g.indptr[row_i], sub_g.indptr[row_i + 1]
            c_idx = sub_g.indices[s0:s1]
            c_val = sub_g.data[s0:s1].astype(np.float64)
            p_nz = c_val if is_l1_already else (c_val / rs_g[row_i])
            p_loo_nz = np.maximum(tot_g[c_idx] - p_nz, 0.0) / float(n_r - 1)
            bc = float(np.clip(np.sum(np.sqrt(p_nz * p_loo_nz)), 0.0, 1.0))
            h_sum += float(np.sqrt(max(0.0, 1.0 - bc)))
        out[int(g_id)] = (h_sum, int(sub_g.shape[0]))
    return out


def _subchunk_pass2_worker(task: tuple[Any, ...]) -> dict[str, Any]:
    """Compute Pass-2 A1 distances and (when requested) A3 centroids/LOO distances using exact Pass-1 seeds."""
    (
        data,
        indices,
        indptr,
        n_raw_vars,
        d_idx,
        c_idx,
        ds_id_code,
        subchunk_id,
        n_donors,
        n_cts,
        a1_pair_elig_1500,
        a1_pair_elig_3000,
        a1_pair_elig_500,
        a1_pool_elig_1500,
        a1_pool_elig_3000,
        a1_pool_elig_500,
        a3_mode,
        a3_lookup,
        a3_pool_lookup,
        a3_batch_id,
        a3_n_cells_ct,
        a3_n_cells_pool,
    ) = task

    sub = sp.csr_matrix((data.astype(np.float64, copy=False), indices, indptr), shape=(len(d_idx), n_raw_vars))
    sub_clean = (sub @ _W_GPROJ).tocsr()
    del sub
    sub_clean.eliminate_zeros()
    umi = np.rint(np.asarray(sub_clean.sum(axis=1)).ravel()).astype(np.int64)
    p_idx = d_idx * n_cts + c_idx

    out: dict[str, Any] = {"a1": {}, "a3_centroids": {}, "a3_dist": {}}

    # Batch-scoped A3 centroid-only subpass (Yazar Pass 2b part 1)
    if a3_mode == "centroids_only":
        m_15 = np.flatnonzero(umi >= 1500)
        if len(m_15):
            rng_a3 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 31500 + a3_batch_id, subchunk_id, 0]))
            mat_15 = sub_clean[m_15]
            r_a3 = _rarefy_csr(mat_15, 1500, rng_a3)

            def _batch_sum_csr(mat_s: sp.csr_matrix, ids_s: np.ndarray, lkup: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
                mp = lkup[ids_s]
                ok = mp >= 0
                if not np.any(ok):
                    return None
                u_m, inv = np.unique(mp[ok], return_inverse=True)
                ind = sp.csr_matrix((np.ones(len(inv), dtype=np.float64), (inv, np.arange(len(inv)))), shape=(len(u_m), len(inv)))
                summed = (ind @ mat_s[ok]).tocsr()
                summed.sum_duplicates()
                return (u_m.astype(np.int64, copy=False), summed.data, summed.indices, summed.indptr)

            out["a3_centroids"]["ct_1500"] = _batch_sum_csr(r_a3, p_idx[m_15], a3_lookup)
            if a3_pool_lookup is not None:
                out["a3_centroids"]["pool_1500"] = _batch_sum_csr(r_a3, d_idx[m_15], a3_pool_lookup)
            del r_a3
            inv_tot = 1.0 / np.maximum(umi[m_15].astype(np.float64), 1.0)
            mat_l1 = mat_15.multiply(inv_tot[:, None]).tocsr()
            out["a3_centroids"]["ct_full_l1"] = _batch_sum_csr(mat_l1, p_idx[m_15], a3_lookup)
            if a3_pool_lookup is not None:
                out["a3_centroids"]["pool_full_l1"] = _batch_sum_csr(mat_l1, d_idx[m_15], a3_pool_lookup)
            del mat_l1
        del sub_clean
        return out

    # A1 computation helper via imported hellinger_distance & reference_without_donor
    def _eval_a1_arm(mat_arm: sp.csr_matrix, rows_orig: np.ndarray, ref_dict: dict[str, Any], arm_key: str, is_pooled: bool) -> dict[tuple[int, str], tuple[float, int]]:
        res: dict[tuple[int, str], tuple[float, int]] = {}
        if mat_arm.shape[0] == 0 or not ref_dict:
            return res
        elig_arr = a1_pool_elig_1500 if (is_pooled and "3000" not in arm_key and "full_500" not in arm_key) else None
        if is_pooled:
            if "3000" in arm_key:
                elig_arr = a1_pool_elig_3000
            elif "full_500" in arm_key:
                elig_arr = a1_pool_elig_500
            ids = d_idx[rows_orig]
            valid = np.flatnonzero(elig_arr[ids] >= 0)
        else:
            e_arr = a1_pair_elig_3000 if "3000" in arm_key else (a1_pair_elig_500 if "full_500" in arm_key else a1_pair_elig_1500)
            ids = p_idx[rows_orig]
            valid = np.flatnonzero(e_arr[ids] >= 0)
        if len(valid) == 0:
            return res

        sub_v = mat_arm[valid]
        v_rows = rows_orig[valid]
        v_d = d_idx[v_rows]
        v_c = c_idx[v_rows] if not is_pooled else np.zeros(len(v_rows), dtype=np.int32)

        # Group by (donor, ct) using safe stride >= max n_cts
        keys = v_d.astype(np.int64) * 100_000 + v_c.astype(np.int64)
        for kv in np.unique(keys):
            km = keys == kv
            di = int(v_d[km][0])
            ci = int(v_c[km][0])
            g_name = "_POOLED_" if is_pooled else _W_REFS["cell_types"][ci]
            d_str = _W_REFS["donor_ids"][di]
            if g_name not in ref_dict:
                continue
            ref_stats = ref_dict[g_name]
            sub_cells = sub_v[km]
            # LOO young ref distance (primary)
            if "young" in ref_stats:
                r_loo = reference_without_donor(ref_stats["young"], d_str)
                if r_loo["n_cells"] >= MIN_LOO_REF_CELLS and r_loo["counts"].sum() > 0:
                    hd = hellinger_distance(sub_cells, r_loo["counts"])
                    res[(di * (1 if is_pooled else n_cts) + (0 if is_pooled else ci), f"{arm_key}_y_loo")] = (float(hd.sum()), int(len(hd)))
                # Non-LOO self-baseline (YOUNG donors only -> unsubtracted full YOUNG reference, strictly biased low)
                if d_str in ref_stats["young"]["donors"]:
                    r_full = ref_stats["young"]
                    if r_full["n_cells"] >= MIN_LOO_REF_CELLS and r_full["counts"].sum() > 0:
                        hd_s = hellinger_distance(sub_cells, r_full["counts"])
                        res[(di * (1 if is_pooled else n_cts) + (0 if is_pooled else ci), f"{arm_key}_y_self_nonloo")] = (float(hd_s.sum()), int(len(hd_s)))
            # Split-half Half 2 -> Half 1 ref distance (supports intra-cell-type fallback when global H1/H2 empty)
            r_h1 = ref_stats.get("half1")
            if r_h1 is not None and d_str in r_h1.get("h2_donors", set()):
                if r_h1["n_cells"] >= MIN_LOO_REF_CELLS and r_h1["counts"].sum() > 0:
                    hd_sh = hellinger_distance(sub_cells, r_h1["counts"])
                    res[(di * (1 if is_pooled else n_cts) + (0 if is_pooled else ci), f"{arm_key}_h2_to_h1")] = (float(hd_sh.sum()), int(len(hd_sh)))
            # Mature (MID) ref distance (MID + OLD donors)
            if "mid" in ref_stats and _W_REFS["age_bin"][di] in (2, 3):
                r_m = reference_without_donor(ref_stats["mid"], d_str)
                if r_m["n_cells"] >= MIN_LOO_REF_CELLS and r_m["counts"].sum() > 0:
                    hd_m = hellinger_distance(sub_cells, r_m["counts"])
                    res[(di * (1 if is_pooled else n_cts) + (0 if is_pooled else ci), f"{arm_key}_mid_loo")] = (float(hd_m.sum()), int(len(hd_m)))
        return res

    # Run A1 when requested (a3_mode != "distances_only")
    if a3_mode != "distances_only":
        m_15 = np.flatnonzero(umi >= 1500)
        if len(m_15):
            mat_15 = sub_clean[m_15]
            out["a1"].update(_eval_a1_arm(mat_15, m_15, _W_REFS.get("ct_full_1500", {}), "a1_full_1500_ct", False))
            out["a1"].update(_eval_a1_arm(mat_15, m_15, _W_REFS.get("pool_full_1500", {}), "a1_full_1500_pool", True))

            rng_15 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 1500, subchunk_id, 0]))
            ds_15 = multinomial_downsample(mat_15, 1500, rng_15)
            out["a1"].update(_eval_a1_arm(ds_15, m_15, _W_REFS.get("ct_1500", {}), "a1_1500_ct", False))
            out["a1"].update(_eval_a1_arm(ds_15, m_15, _W_REFS.get("pool_1500", {}), "a1_1500_pool", True))
            del ds_15

            rng_5s = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 5001500, subchunk_id, 0]))
            ds_5s = multinomial_downsample(mat_15, 500, rng_5s)
            out["a1"].update(_eval_a1_arm(ds_5s, m_15, _W_REFS.get("ct_500_same", {}), "a1_500_ct", False))
            out["a1"].update(_eval_a1_arm(ds_5s, m_15, _W_REFS.get("pool_500_same", {}), "a1_500_pool", True))
            del ds_5s

        m_3k = np.flatnonzero(umi >= 3000)
        if len(m_3k):
            rng_3k = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 3000, subchunk_id, 0]))
            ds_3k = multinomial_downsample(sub_clean[m_3k], 3000, rng_3k)
            out["a1"].update(_eval_a1_arm(ds_3k, m_3k, _W_REFS.get("ct_3000", {}), "a1_3000_ct", False))
            out["a1"].update(_eval_a1_arm(ds_3k, m_3k, _W_REFS.get("pool_3000", {}), "a1_3000_pool", True))
            del ds_3k

        if ds_id_code == 1 and _W_REFS.get("run_gao_500_a1"):
            m_5 = np.flatnonzero(umi >= 500)
            if len(m_5):
                mat_5 = sub_clean[m_5]
                out["a1"].update(_eval_a1_arm(mat_5, m_5, _W_REFS.get("ct_full_500", {}), "a1_full_500_ct", False))
                out["a1"].update(_eval_a1_arm(mat_5, m_5, _W_REFS.get("pool_full_500", {}), "a1_full_500_pool", True))

    # A3 LOO distances when requested ("with_a1" for Gao in 1 pass, "distances_only" for Yazar batch)
    if a3_mode in ("with_a1", "distances_only"):
        m_15 = np.flatnonzero(umi >= 1500)
        if len(m_15):
            mat_15 = sub_clean[m_15]
            rng_a3 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_id_code, 31500 + a3_batch_id, subchunk_id, 0]))
            r_a3 = _rarefy_csr(mat_15, 1500, rng_a3)
            # Cell-level normalized counts / 1500 so centroid sum is probability sum
            r_prob = r_a3.multiply(1.0 / 1500.0).tocsr()
            out["a3_dist"]["ct_1500"] = _compute_a3_loo_distances_for_matrix(
                r_prob, p_idx[m_15], _W_REFS["a3_ct_1500"], a3_n_cells_ct, a3_lookup, is_l1_already=True
            )
            if a3_pool_lookup is not None:
                out["a3_dist"]["pool_1500"] = _compute_a3_loo_distances_for_matrix(
                    r_prob, d_idx[m_15], _W_REFS["a3_pool_1500"], a3_n_cells_pool, a3_pool_lookup, is_l1_already=True
                )
            del r_a3, r_prob

            inv_tot = 1.0 / np.maximum(umi[m_15].astype(np.float64), 1.0)
            mat_l1 = mat_15.multiply(inv_tot[:, None]).tocsr()
            out["a3_dist"]["ct_full"] = _compute_a3_loo_distances_for_matrix(
                mat_l1, p_idx[m_15], _W_REFS["a3_ct_full"], a3_n_cells_ct, a3_lookup, is_l1_already=True
            )
            if a3_pool_lookup is not None:
                out["a3_dist"]["pool_full"] = _compute_a3_loo_distances_for_matrix(
                    mat_l1, d_idx[m_15], _W_REFS["a3_pool_full"], a3_n_cells_pool, a3_pool_lookup, is_l1_already=True
                )
            del mat_l1

    del sub_clean
    return out


def compute_a2_variance_metrics(
    expr_mat: np.ndarray,
    is_comp: np.ndarray,
    cov_mat_with_intercept: np.ndarray | None,
    batch_col_count: int = 0,
    n_perm: int = A2_N_PERM,
    seed: int = MASTER_SEED,
    eps: float = 1e-8,
    min_leverage_one_minus_h: float = 1e-4,
    custom_perm_orders: np.ndarray | None = None,
    return_perm: bool = False,
) -> dict[str, Any] | tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Compute Item-8 eta2 (Freedman-Lane if covariates given) and LVR (studentized full-model mean(R^2) or depth-only)."""
    n_d, n_g = expr_mat.shape
    is_c = np.asarray(is_comp, dtype=bool)
    n_c = int(is_c.sum())
    n_y = n_d - n_c
    if n_c < 2 or n_y < 2 or n_g == 0:
        raise ValueError(f"Insufficient sample size or genes for A2: n_c={n_c}, n_y={n_y}, n_g={n_g}")

    rng = np.random.default_rng(seed)
    c_float = is_c.astype(np.float64)
    if custom_perm_orders is not None:
        if custom_perm_orders.shape != (n_perm, n_d):
            raise ValueError(f"custom_perm_orders shape {custom_perm_orders.shape} != ({n_perm}, {n_d})")

    if cov_mat_with_intercept is None:
        # Unresidualized / direct Item-8 depth-only input: use exact _stratum_stats observed formulas + vectorized perm
        m_y = expr_mat[~is_c].mean(axis=0)
        m_c = expr_mat[is_c].mean(axis=0)
        m_tot = expr_mat.mean(axis=0)
        ss_between = n_y * (m_y - m_tot) ** 2 + n_c * (m_c - m_tot) ** 2
        var_y = expr_mat[~is_c].var(axis=0, ddof=1)
        var_c = expr_mat[is_c].var(axis=0, ddof=1)
        ss_within = (n_y - 1) * var_y + (n_c - 1) * var_c
        ss_total = np.maximum(ss_between + ss_within, 1e-12)
        eta2_g = ss_between / ss_total
        f_g = (ss_between / 1.0) / np.maximum(ss_within / float(n_d - 2), 1e-12)

        bias_diff = _df_logvar_bias(n_c) - _df_logvar_bias(n_y)
        lvr_g = np.log((var_c + eps) / (var_y + eps)) - bias_diff
        obs_mean_eta2 = float(np.mean(eta2_g))
        obs_med_eta2 = float(np.median(eta2_g))
        obs_mean_f = float(np.mean(f_g))
        obs_mean_lvr = float(np.mean(lvr_g))
        obs_med_lvr = float(np.median(lvr_g))
        obs_tot_lvr = float(np.log(max(float(np.mean(var_c)), 1e-12) / max(float(np.mean(var_y)), 1e-12)))
        obs_frac_gt = float(np.mean(var_c > var_y))
        resid_df_model = n_d - 2

        total_ss_all = np.maximum(np.sum((expr_mat - m_tot[None, :]) ** 2, axis=0), 1e-12)
        sum_all = expr_mat.sum(axis=0)
        sum_sq_all = np.sum(expr_mat ** 2, axis=0)
        perm_eta2 = np.empty(n_perm, dtype=np.float64)
        perm_lvr = np.empty(n_perm, dtype=np.float64)
        perm_tot_lvr = np.empty(n_perm, dtype=np.float64)
        for pi in range(n_perm):
            idx_c = custom_perm_orders[pi][is_c] if custom_perm_orders is not None else rng.choice(n_d, size=n_c, replace=False)
            sub_c = expr_mat[idx_c]
            sc = sub_c.sum(axis=0)
            sy = sum_all - sc
            mc_p = sc / float(n_c)
            my_p = sy / float(n_y)
            p_ssb = n_y * (my_p - m_tot) ** 2 + n_c * (mc_p - m_tot) ** 2
            perm_eta2[pi] = float(np.mean(p_ssb / total_ss_all))

            sqc = np.sum(sub_c ** 2, axis=0)
            sqy = sum_sq_all - sqc
            vc_p = np.maximum((sqc - (sc ** 2) / float(n_c)) / float(n_c - 1), 0.0)
            vy_p = np.maximum((sqy - (sy ** 2) / float(n_y)) / float(n_y - 1), 0.0)
            perm_lvr[pi] = float(np.mean(np.log((vc_p + eps) / (vy_p + eps))) - bias_diff)
            perm_tot_lvr[pi] = float(np.log(max(float(np.mean(vc_p)), 1e-12) / max(float(np.mean(vy_p)), 1e-12)))
    else:
        # Full batch + sex + depth adjusted model: Freedman-Lane partial eta2 + studentized full-model LVR using mean(R^2)
        z_mat = np.asarray(cov_mat_with_intercept, dtype=np.float64)
        q_z, _ = np.linalg.qr(z_mat, mode="reduced")
        # Reduced model residuals Y_Z = (I - H_Z) Y, a_Z = (I - H_Z) is_comp
        e_z = expr_mat - q_z @ (q_z.T @ expr_mat)
        a_z = c_float - q_z @ (q_z.T @ c_float)
        a_z_norm2 = float(np.dot(a_z, a_z))
        if a_z_norm2 <= 1e-10:
            raise ValueError("is_comp is collinear with nuisance covariates Z")

        # Full model [Z, is_comp] residuals & hat diagonal
        x_full = np.column_stack([z_mat, c_float])
        q_full, _ = np.linalg.qr(x_full, mode="reduced")
        h_full = np.sum(q_full ** 2, axis=1)
        one_minus_h = 1.0 - h_full
        if np.any(one_minus_h <= min_leverage_one_minus_h):
            raise ValueError(f"Extreme leverage encountered: min(1-h)={float(one_minus_h.min()):.3e}")

        resid_df_model = n_d - q_full.shape[1]
        e_full = expr_mat - q_full @ (q_full.T @ expr_mat)
        ss_resid_full = np.maximum(np.sum(e_full ** 2, axis=0), 1e-12)
        ss_tot_z = np.maximum(np.sum(e_z ** 2, axis=0), 1e-12)
        ss_age_given_z = np.maximum(ss_tot_z - ss_resid_full, 0.0)
        eta2_g = ss_age_given_z / np.maximum(ss_age_given_z + ss_resid_full, 1e-12)
        f_g = (ss_age_given_z / 1.0) / np.maximum(ss_resid_full / float(resid_df_model), 1e-12)

        obs_mean_eta2 = float(np.mean(eta2_g))
        obs_med_eta2 = float(np.median(eta2_g))
        obs_mean_f = float(np.mean(f_g))

        # Freedman-Lane permutation: permute reduced-model residuals e_z across rows & regress on a_z (= project onto span(is_comp | Z))
        r_mat = e_z  # shape (n_d, n_g)
        e_z_norm2 = np.maximum(np.sum(r_mat ** 2, axis=0), 1e-12)  # ss_tot_z
        a_unit = a_z / np.sqrt(a_z_norm2)
        perm_eta2 = np.empty(n_perm, dtype=np.float64)
        for pi in range(n_perm):
            p_idx = custom_perm_orders[pi] if custom_perm_orders is not None else rng.permutation(n_d)
            e_pi = r_mat[p_idx]
            e_pi_z = e_pi - q_z @ (q_z.T @ e_pi)
            proj_a = a_unit @ e_pi_z  # (n_g,)
            ss_a_pi = proj_a ** 2
            ss_tot_pi = np.maximum(np.sum(e_pi_z ** 2, axis=0), 1e-12)
            perm_eta2[pi] = float(np.mean(ss_a_pi / ss_tot_pi))

        # Studentized full-model residuals R_i = e_full_i / sqrt(1 - h_full_ii), uncentered second moment mean(R[g]^2)
        r_stud = e_full / np.sqrt(one_minus_h)[:, None]
        df_c_eff = float(np.sum(one_minus_h[is_c]))
        df_y_eff = float(np.sum(one_minus_h[~is_c]))
        bias_diff = _df_logvar_bias_from_df(df_c_eff) - _df_logvar_bias_from_df(df_y_eff)

        var_c = np.mean(r_stud[is_c] ** 2, axis=0)
        var_y = np.mean(r_stud[~is_c] ** 2, axis=0)
        lvr_g = np.log((var_c + eps) / (var_y + eps)) - bias_diff
        obs_mean_lvr = float(np.mean(lvr_g))
        obs_med_lvr = float(np.median(lvr_g))
        obs_tot_lvr = float(np.log(max(float(np.mean(var_c)), 1e-12) / max(float(np.mean(var_y)), 1e-12)))
        obs_frac_gt = float(np.mean(var_c > var_y))

        r2_mat = r_stud ** 2
        sum_r2 = np.sum(r2_mat, axis=0)
        perm_lvr = np.empty(n_perm, dtype=np.float64)
        perm_tot_lvr = np.empty(n_perm, dtype=np.float64)
        for pi in range(n_perm):
            idx_c = custom_perm_orders[pi][is_c] if custom_perm_orders is not None else rng.choice(n_d, size=n_c, replace=False)
            sc = np.sum(r2_mat[idx_c], axis=0)
            sy = np.maximum(sum_r2 - sc, 0.0)
            vc_p = sc / float(n_c)
            vy_p = sy / float(n_y)
            perm_lvr[pi] = float(np.mean(np.log((vc_p + eps) / (vy_p + eps))) - bias_diff)
            perm_tot_lvr[pi] = float(np.log(max(float(np.mean(vc_p)), 1e-12) / max(float(np.mean(vy_p)), 1e-12)))

    res = _summarize_a2_perm_stats(
        n_g=int(n_g),
        resid_df_model=int(resid_df_model),
        obs_mean_eta2=obs_mean_eta2,
        obs_med_eta2=obs_med_eta2,
        obs_mean_f=obs_mean_f,
        obs_mean_lvr=obs_mean_lvr,
        obs_med_lvr=obs_med_lvr,
        obs_tot_lvr=obs_tot_lvr,
        obs_frac_gt=obs_frac_gt,
        perm_eta2=perm_eta2,
        perm_lvr=perm_lvr,
        perm_tot_lvr=perm_tot_lvr,
    )
    if return_perm:
        return res, {"perm_eta2": perm_eta2, "perm_lvr": perm_lvr, "perm_tot_lvr": perm_tot_lvr}
    return res


def _summarize_a2_perm_stats(
    n_g: int,
    resid_df_model: int,
    obs_mean_eta2: float,
    obs_med_eta2: float,
    obs_mean_f: float,
    obs_mean_lvr: float,
    obs_med_lvr: float,
    obs_tot_lvr: float,
    obs_frac_gt: float,
    perm_eta2: np.ndarray,
    perm_lvr: np.ndarray,
    perm_tot_lvr: np.ndarray,
) -> dict[str, Any]:
    """Summarize observed A2 eta2/LVR metrics against empirical permutation distributions."""
    n_perm = int(len(perm_eta2))
    p_eta2_perm = float((np.sum(perm_eta2 >= obs_mean_eta2 - 1e-12) + 1) / (n_perm + 1))
    p_lvr_gt = float((np.sum(perm_lvr >= obs_mean_lvr - 1e-12) + 1) / (n_perm + 1))
    p_lvr_lt = float((np.sum(perm_lvr <= obs_mean_lvr + 1e-12) + 1) / (n_perm + 1))
    p_lvr_two = min(1.0, 2.0 * min(p_lvr_gt, p_lvr_lt))
    p_tot_gt = float((np.sum(perm_tot_lvr >= obs_tot_lvr - 1e-12) + 1) / (n_perm + 1))
    p_tot_lt = float((np.sum(perm_tot_lvr <= obs_tot_lvr + 1e-12) + 1) / (n_perm + 1))
    p_tot_two = min(1.0, 2.0 * min(p_tot_gt, p_tot_lt))

    null_eta2_sd = float(np.std(perm_eta2, ddof=1))
    null_lvr_sd = float(np.std(perm_lvr, ddof=1))
    return {
        "n_genes": int(n_g),
        "residual_df_resid_model": int(resid_df_model),
        "mean_eta2": obs_mean_eta2,
        "median_eta2": obs_med_eta2,
        "mean_F": obs_mean_f,
        "mean_log_var_ratio": obs_mean_lvr,
        "median_log_var_ratio": obs_med_lvr,
        "log_total_var_ratio": obs_tot_lvr,
        "frac_genes_comp_gt_young": obs_frac_gt,
        "n_perm": n_perm,
        "perm_null_mean_eta2": float(np.mean(perm_eta2)),
        "perm_null_sd_eta2": null_eta2_sd,
        "z_eta2_perm": (obs_mean_eta2 - float(np.mean(perm_eta2))) / null_eta2_sd if null_eta2_sd > 0 else np.nan,
        "p_eta2_perm": p_eta2_perm,
        "perm_null_mean_lvr": float(np.mean(perm_lvr)),
        "perm_null_sd_lvr": null_lvr_sd,
        "z_lvr_perm": (obs_mean_lvr - float(np.mean(perm_lvr))) / null_lvr_sd if null_lvr_sd > 0 else np.nan,
        "p_lvr_one_sided_comp_gt_young": p_lvr_gt,
        "p_lvr_two_sided": p_lvr_two,
        "p_tot_lvr_one_sided": p_tot_gt,
        "p_tot_lvr_two_sided": p_tot_two,
    }


def _prune_covariate_matrix(
    df: pd.DataFrame,
    nonbatch_covs: list[str],
    batch_col: str | None,
    min_batch_donors: int,
    min_residual_df: int,
    ref_level: str | None = None,
    require_is_comp_collinear_check: str | None = None,
    atomic_batch_block: bool = True,
) -> tuple[pd.DataFrame, list[str], list[str], str]:
    """Select full-rank covariate slice with atomic batch-block df check and singleton exclusion."""
    work = df.copy()
    status = "evaluated"
    n_excluded_singleton = 0

    active_batch_cols: list[str] = []
    omitted: list[str] = []
    if batch_col is not None:
        vc = work[batch_col].value_counts()
        ge_batches = vc[vc >= min_batch_donors].index.tolist()
        if min_batch_donors >= 2:
            keep_m = work[batch_col].isin(ge_batches)
            n_excluded_singleton = int((~keep_m).sum())
            work = work.loc[keep_m].copy()
        vc2 = work[batch_col].value_counts()
        if len(vc2) < 2:
            return df.copy(), [], [*nonbatch_covs, batch_col], "not_evaluable_batch_single_level"
        levels_sorted = [k for k, _ in sorted(vc2.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))]
        ref = ref_level if (ref_level is not None and ref_level in levels_sorted) else levels_sorted[0]
        for lv in levels_sorted:
            if lv != ref:
                cname = f"batch_{str(lv).replace(' ', '_')}"
                work[cname] = (work[batch_col].astype(str) == str(lv)).astype(float)
                active_batch_cols.append(cname)

    # Non-batch covariates first (strictly rank-checked, never silently dropping age)
    cols_final: list[str] = []
    base_cols = [work[require_is_comp_collinear_check].to_numpy(dtype=float)] if require_is_comp_collinear_check else []
    design_cols = [np.ones(len(work), dtype=float), *base_cols]
    for cov in nonbatch_covs:
        vals = work[cov].to_numpy(dtype=float)
        if not np.isfinite(vals).all() or np.unique(vals).size < 2:
            omitted.append(f"{cov}(constant_or_na)")
            continue
        cand = np.column_stack([*design_cols, vals])
        if np.linalg.matrix_rank(cand) == cand.shape[1]:
            design_cols.append(vals)
            cols_final.append(cov)
        else:
            omitted.append(f"{cov}(collinear)")

    if active_batch_cols:
        # Iteratively drop collinear batch dummies (and check residual df / leverage atomically if atomic_batch_block)
        kept_b: list[str] = []
        for bc in active_batch_cols:
            vals = work[bc].to_numpy(dtype=float)
            cand = np.column_stack([*design_cols, *[work[c].to_numpy(dtype=float) for c in kept_b], vals])
            if np.linalg.matrix_rank(cand) == cand.shape[1]:
                kept_b.append(bc)
            else:
                omitted.append(f"{bc}(collinear)")
        n_params_total = len(design_cols) + len(kept_b)
        if (len(work) - n_params_total) < min_residual_df or (atomic_batch_block and len(kept_b) == 0):
            return df.copy(), [c for c in cols_final if c in df.columns], [*omitted, "batch_block(df_limit)"], "not_evaluable_within_batch_df"
        cols_final.extend(kept_b)

    final_mat = np.column_stack([np.ones(len(work), dtype=float), *base_cols, *[work[c].to_numpy(dtype=float) for c in cols_final]])
    if (len(work) - final_mat.shape[1]) < min_residual_df:
        return work, cols_final, omitted, "insufficient_df"
    work.attrs["n_excluded_singleton"] = n_excluded_singleton
    work.attrs["status"] = status
    return work, cols_final, omitted, status


def fit_ols_model_row(
    df_group: pd.DataFrame,
    outcome_col: str,
    predictor_col: str,
    nonbatch_covs: list[str],
    batch_col: str | None,
    min_batch_donors: int,
    min_residual_df: int,
    ref_level: str | None,
    meta_fields: dict[str, Any],
    fallback_to_nobatch: bool = False,
) -> dict[str, Any] | None:
    """Fit donor-level OLS model with standardized slope, MDE80, nested batch R2, and status metadata."""
    work = df_group.dropna(subset=[outcome_col, predictor_col]).copy()
    if len(work) < MIN_DONORS_SLOPE:
        return None
    if work[predictor_col].nunique() < 2:
        return None

    yazar_status = "evaluated"
    sub_df, cols_active, cols_omitted, st = _prune_covariate_matrix(
        work,
        nonbatch_covs=nonbatch_covs,
        batch_col=batch_col,
        min_batch_donors=min_batch_donors,
        min_residual_df=min_residual_df,
        ref_level=ref_level,
        require_is_comp_collinear_check=predictor_col,
        atomic_batch_block=True,
    )
    if st != "evaluated":
        if batch_col is not None and fallback_to_nobatch:
            yazar_status = meta_fields.get("default_fallback_status", "not_evaluable_within_batch_df_no_pool")
            sub_df, cols_active, cols_omitted, st2 = _prune_covariate_matrix(
                work,
                nonbatch_covs=nonbatch_covs,
                batch_col=None,
                min_batch_donors=1,
                min_residual_df=5,
                ref_level=None,
                require_is_comp_collinear_check=predictor_col,
            )
            if st2 != "evaluated":
                return None
            cols_omitted.append("batch_fe(fallback_no_batch)")
        else:
            return None

    y = sub_df[outcome_col].to_numpy(dtype=float)
    outcome_sd = float(np.std(y, ddof=1))
    if not np.isfinite(outcome_sd) or outcome_sd <= 0.0:
        return None

    x_cols = [predictor_col, *cols_active]
    X = sm.add_constant(sub_df[x_cols].astype(float), has_constant="add")
    fit = sm.OLS(y, X).fit()
    df_resid = int(fit.df_resid)
    if df_resid < 1:
        return None

    slope_raw = float(fit.params[predictor_col])
    se_raw = float(fit.bse[predictor_col])
    t975 = float(scipy_stats.t.ppf(0.975, df_resid))
    t80 = float(scipy_stats.t.ppf(0.80, df_resid))
    ci_low_raw = slope_raw - t975 * se_raw
    ci_high_raw = slope_raw + t975 * se_raw
    mde_raw = (t975 + t80) * se_raw
    resid_sd = float(np.sqrt(fit.mse_resid))

    batch_delta_r2 = np.nan
    batch_partial_r2 = np.nan
    has_active_batch = any(c.startswith("batch_") for c in cols_active)
    if has_active_batch:
        nb_cols = [c for c in x_cols if not c.startswith("batch_")]
        X_nb = sm.add_constant(sub_df[nb_cols].astype(float), has_constant="add")
        fit_nb = sm.OLS(y, X_nb).fit()
        batch_delta_r2 = float(fit.rsquared - fit_nb.rsquared)
        sse_nb = float(np.sum(fit_nb.resid ** 2))
        sse_full = float(np.sum(fit.resid ** 2))
        batch_partial_r2 = float((sse_nb - sse_full) / sse_nb) if sse_nb > 1e-12 else 0.0

    meta_clean = {k: v for k, v in meta_fields.items() if k != "default_fallback_status"}
    row = {
        **meta_clean,
        "yazar_pool_fe_status": meta_clean.get("yazar_pool_fe_status", yazar_status),
        "n_donors": int(len(sub_df)),
        "n_cells": int(sub_df["n_cells"].sum()) if "n_cells" in sub_df else np.nan,
        "residual_df": df_resid,
        "covariates_active": ";".join(cols_active) if cols_active else "none",
        "covariates_omitted": ";".join(cols_omitted) if cols_omitted else "none",
        "outcome_sd": outcome_sd,
        "resid_sd": resid_sd,
        "slope_raw_per_decade": slope_raw,
        "se_raw": se_raw,
        "ci95_low_raw": ci_low_raw,
        "ci95_high_raw": ci_high_raw,
        "slope_sd_per_decade": slope_raw / outcome_sd,
        "se_sd": se_raw / outcome_sd,
        "ci95_low_sd": ci_low_raw / outcome_sd,
        "ci95_high_sd": ci_high_raw / outcome_sd,
        "pvalue": float(fit.pvalues[predictor_col]),
        "q_bh_within_slice": np.nan,
        "family_size_m": np.nan,
        "mde80_raw_per_decade": mde_raw,
        "mde80_sd_per_decade": mde_raw / outcome_sd,
        "batch_delta_r2": batch_delta_r2,
        "batch_partial_r2": batch_partial_r2,
    }
    return row


def fit_contrast_adj_diff(
    df_group: pd.DataFrame,
    outcome_col: str,
    bin_a: str,
    bin_b: str,
    nonbatch_covs: list[str],
    batch_col: str | None,
    min_batch_donors: int,
    min_residual_df: int,
    ref_level: str | None,
) -> dict[str, float]:
    """Fit batch-adjusted mean contrast (bin_b minus bin_a, e.g., OLD minus YOUNG or OLD minus MID)."""
    sub = df_group.loc[df_group["bin"].isin([bin_a, bin_b])].dropna(subset=[outcome_col]).copy()
    nan_res = {"diff": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    if int((sub["bin"] == bin_a).sum()) < 2 or int((sub["bin"] == bin_b).sum()) < 2:
        return nan_res
    sub["_is_b"] = (sub["bin"] == bin_b).astype(float)
    row = fit_ols_model_row(
        sub,
        outcome_col=outcome_col,
        predictor_col="_is_b",
        nonbatch_covs=nonbatch_covs,
        batch_col=batch_col,
        min_batch_donors=min_batch_donors,
        min_residual_df=min_residual_df,
        ref_level=ref_level,
        meta_fields={},
        fallback_to_nobatch=True,
    )
    if row is None:
        return nan_res
    return {
        "diff": float(row["slope_raw_per_decade"]),
        "se": float(row["se_raw"]),
        "ci_low": float(row["ci95_low_raw"]),
        "ci_high": float(row["ci95_high_raw"]),
        "p": float(row["pvalue"]),
    }


def evaluate_verdicts(
    models_part_a: pd.DataFrame,
    variance_df: pd.DataFrame,
    spread_df: pd.DataFrame,
    robustness_df: pd.DataFrame,
    part_b_primary_depth: int,
    pooled_name: str,
) -> dict[str, Any]:
    """Programmatically evaluate strict pre-stated Section 9 verdicts without softening."""
    out: dict[str, Any] = {}

    # A1 verdict
    a1_full = models_part_a[
        (models_part_a["part"] == "A1")
        & (models_part_a["depth"] == "full_1500")
        & (models_part_a["reference_type"] == "young_loo")
        & (models_part_a["model_spec"] == "primary")
        & (~models_part_a["is_stratified_global"].astype(bool))
        & (models_part_a["group"] == pooled_name)
    ]
    a1_1500 = models_part_a[
        (models_part_a["part"] == "A1")
        & (models_part_a["depth"] == "1500")
        & (models_part_a["reference_type"] == "young_loo")
        & (models_part_a["model_spec"] == "primary")
        & (~models_part_a["is_stratified_global"].astype(bool))
        & (models_part_a["group"] == pooled_name)
    ]
    if a1_full.empty or a1_1500.empty:
        out["A1"] = {"verdict": "INCONCLUSIVE", "reason": "Missing A1 primary pooled model rows"}
    else:
        rf = a1_full.iloc[0]
        r15 = a1_1500.iloc[0]
        cond_full = float(rf["slope_raw_per_decade"]) > 0.0 and float(rf["pvalue"]) < 0.05
        cond_1500 = float(r15["slope_raw_per_decade"]) > 0.0 and float(r15["pvalue"]) < 0.05
        cond_base = float(rf["old_mean"]) > float(rf["split_half_young_baseline"])
        if cond_full and cond_1500 and cond_base:
            out["A1"] = {"verdict": "PROVES", "reason": "Full and 1500 pooled slopes >0 with p<0.05 and OLD > split-half young baseline"}
        elif float(rf["mde80_sd_per_decade"]) <= 0.5:
            out["A1"] = {
                "verdict": "DOES NOT PROVE",
                "reason": (
                    f"Primary pooled A1 fails PROVES (full slope={rf['slope_raw_per_decade']:.4f} p={rf['pvalue']:.4g}, "
                    f"1500 slope={r15['slope_raw_per_decade']:.4f} p={r15['pvalue']:.4g}, "
                    f"OLD={rf['old_mean']:.4f} vs split_half={rf['split_half_young_baseline']:.4f}) with MDE80={rf['mde80_sd_per_decade']:.3f}<=0.5 SD/dec"
                ),
            }
        else:
            out["A1"] = {"verdict": "INCONCLUSIVE", "reason": f"Primary pooled A1 did not prove drift and MDE80={rf['mde80_sd_per_decade']:.3f}>0.5"}

    # A2 verdict (OLD_vs_YOUNG, eps_used == 1e-8)
    a2_b = variance_df[
        (variance_df["contrast"] == "OLD_vs_YOUNG")
        & (variance_df["arm"] == "batch_residualized")
        & (np.isclose(variance_df["eps_used"].astype(float), 1e-8))
        & (~variance_df["is_stratified_global"].astype(bool))
        & (variance_df["group"] == pooled_name)
    ]
    a2_d = variance_df[
        (variance_df["contrast"] == "OLD_vs_YOUNG")
        & (variance_df["arm"] == "depth_only_residualized")
        & (np.isclose(variance_df["eps_used"].astype(float), 1e-8))
        & (~variance_df["is_stratified_global"].astype(bool))
        & (variance_df["group"] == pooled_name)
    ]
    if a2_b.empty or a2_d.empty:
        out["A2"] = {"verdict": "INCONCLUSIVE", "reason": "Missing A2 OLD_vs_YOUNG pooled rows"}
    else:
        rb = a2_b.iloc[0]
        rd = a2_d.iloc[0]
        sig_b = float(rb["p_eta2_perm"]) < 0.05 or (float(rb["p_lvr_one_sided_comp_gt_young"]) < 0.05 and float(rb["mean_log_var_ratio"]) > 0.0)
        same_dir_d = float(rb["mean_log_var_ratio"]) > 0.0 and float(rd["mean_log_var_ratio"]) > 0.0
        if sig_b and same_dir_d:
            out["A2"] = {"verdict": "PROVES", "reason": "Batch-residualized pooled eta2/LVR perm p<0.05 with OLD>YOUNG and matching positive depth-only direction"}
        elif float(rb["mean_log_var_ratio"]) <= 0.0 and float(rb["p_eta2_perm"]) >= 0.05:
            out["A2"] = {"verdict": "DOES NOT PROVE", "reason": f"Batch-residualized pooled mean_log_var_ratio={rb['mean_log_var_ratio']:.4f}<=0 and p_eta2={rb['p_eta2_perm']:.4g}>=0.05"}
        else:
            out["A2"] = {"verdict": "INCONCLUSIVE", "reason": f"A2 mixed/non-significant (lvr_batch={rb['mean_log_var_ratio']:.4f}, p_lvr={rb['p_lvr_one_sided_comp_gt_young']:.4g}, p_eta2={rb['p_eta2_perm']:.4g}, lvr_depth={rd['mean_log_var_ratio']:.4f})"}

    # A3 verdict
    a3_15 = models_part_a[
        (models_part_a["part"] == "A3")
        & (models_part_a["depth"] == "1500")
        & (models_part_a["model_spec"] == "primary")
        & (~models_part_a["is_stratified_global"].astype(bool))
        & (models_part_a["group"] == pooled_name)
    ]
    a3_fl = models_part_a[
        (models_part_a["part"] == "A3")
        & (models_part_a["depth"] == "full_1500")
        & (models_part_a["model_spec"] == "primary")
        & (~models_part_a["is_stratified_global"].astype(bool))
        & (models_part_a["group"] == pooled_name)
    ]
    if a3_15.empty or a3_fl.empty:
        out["A3"] = {"verdict": "INCONCLUSIVE", "reason": "Missing A3 primary pooled rows"}
    else:
        r15 = a3_15.iloc[0]
        rfl = a3_fl.iloc[0]
        if float(r15["slope_raw_per_decade"]) > 0.0 and float(r15["pvalue"]) < 0.05 and float(rfl["slope_raw_per_decade"]) > 0.0:
            out["A3"] = {"verdict": "PROVES", "reason": f"1500 A3 pooled slope>0 (p={r15['pvalue']:.4g}) and same positive direction at full depth (p={rfl['pvalue']:.4g})"}
        elif float(r15["mde80_sd_per_decade"]) <= 0.5:
            out["A3"] = {"verdict": "DOES NOT PROVE", "reason": f"A3 1500 pooled fails PROVES (slope={r15['slope_raw_per_decade']:.4f}, p={r15['pvalue']:.4g}, full_slope={rfl['slope_raw_per_decade']:.4f}) with MDE80={r15['mde80_sd_per_decade']:.3f}<=0.5"}
        else:
            out["A3"] = {"verdict": "INCONCLUSIVE", "reason": f"A3 fails PROVES and MDE80={r15['mde80_sd_per_decade']:.3f}>0.5"}

    out["Part_A_overall"] = {
        "verdict": out["A1"]["verdict"],
        "subparts": {"A1": out["A1"]["verdict"], "A2": out["A2"]["verdict"], "A3": out["A3"]["verdict"]},
    }

    # Part B verdicts (per metric separately, displayed side-by-side)
    part_b_metrics: dict[str, Any] = {}
    for metric, pos_sign in (("shannon", 1), ("deman_kl_methods", -1), ("deman_kl_sqrt", -1)):
        m_sub = spread_df[
            (spread_df["metric"] == metric)
            & (spread_df["model_spec"] == "primary")
            & (~spread_df["is_stratified_global"].astype(bool))
            & (spread_df["group"] == pooled_name)
        ]
        by_d = {str(r["depth"]): r for _, r in m_sub.iterrows()}
        pri_key = str(part_b_primary_depth)
        rob_sub = robustness_df[
            (robustness_df["metric"] == metric)
            & (robustness_df["depth"] == part_b_primary_depth)
            & (robustness_df["model_spec"] == "primary")
            & (robustness_df["group"] == pooled_name)
        ]
        if pri_key not in by_d or not all(str(d) in by_d for d in DEPTH_LEVELS) or rob_sub.empty:
            part_b_metrics[metric] = {"verdict": "INCONCLUSIVE", "reason": "Missing depth or robustness rows"}
            continue
        rp = by_d[pri_key]
        rob_r = rob_sub.iloc[0]
        depths_ok = all((pos_sign * float(by_d[str(d)]["slope_raw_per_decade"]) > 0.0) and float(by_d[str(d)]["pvalue"]) < 0.05 for d in DEPTH_LEVELS)
        rob_ok = bool(rob_r["passes_robustness"]) and (pos_sign * float(rp["slope_raw_per_decade"]) > 0.0)
        if depths_ok and rob_ok:
            part_b_metrics[metric] = {"verdict": "PROVES", "reason": f"{metric} spread increases at 500/1500/3000 (all p<0.05) and passes seed robustness ({int(rob_r['n_same_sign'])}/{int(rob_r['n_seeds'])})"}
        elif float(rp["mde80_sd_per_decade"]) <= 0.5:
            opp_flag = all((pos_sign * float(by_d[str(d)]["slope_raw_per_decade"]) < 0.0) and float(by_d[str(d)]["pvalue"]) < 0.05 for d in DEPTH_LEVELS)
            extra = " (OPPOSITE: spread decreases significantly across all 3 depths)" if opp_flag else ""
            part_b_metrics[metric] = {
                "verdict": "DOES NOT PROVE",
                "reason": f"{metric} fails higher-spread PROVES at primary depth {pri_key} (slope={rp['slope_raw_per_decade']:.4g}, p={rp['pvalue']:.4g}) with MDE80={rp['mde80_sd_per_decade']:.3f}<=0.5{extra}",
            }
        else:
            part_b_metrics[metric] = {"verdict": "INCONCLUSIVE", "reason": f"{metric} fails PROVES with MDE80={rp['mde80_sd_per_decade']:.3f}>0.5"}

    out["Part_B"] = part_b_metrics
    return out


def apply_slice_bh(df: pd.DataFrame, group_cols: list[str], p_col: str = "pvalue", q_col: str = "q_bh_within_slice", m_col: str = "family_size_m") -> pd.DataFrame:
    """Apply Benjamini-Hochberg separately within each explicit non-pooled/non-stratified family slice."""
    out = df.copy()
    if out.empty:
        return out
    out[q_col] = np.nan
    out[m_col] = np.nan
    eligible_mask = (~out["is_pooled"].astype(bool)) & (~out["is_stratified_global"].astype(bool))
    if "yazar_pool_fe_status" in out.columns and "model_spec" in out.columns:
        pool_fe_arm = out["model_spec"].astype(str).isin(["primary", "pool_fe", "pool_fe_no_depth", "age_ge_20"])
        is_fallback = out["yazar_pool_fe_status"].astype(str).str.contains("not_evaluable")
        eligible_mask = eligible_mask & (~(pool_fe_arm & is_fallback))

    for _, idx in out.loc[eligible_mask].groupby(group_cols, sort=False).groups.items():
        pvals = out.loc[idx, p_col].to_numpy(dtype=float)
        valid = np.isfinite(pvals)
        m = int(valid.sum())
        if m == 0:
            continue
        out.loc[idx, m_col] = m
        valid_idx = np.asarray(list(idx))[valid]
        out.loc[valid_idx, q_col] = multipletests(pvals[valid], method="fdr_bh")[1]

    # Pooled rows carry their own singleton family (m = 1, q = p) cleanly separate from per-cell-type family
    pool_idx = out["is_pooled"].astype(bool) & (~out["is_stratified_global"].astype(bool))
    out.loc[pool_idx, m_col] = 1
    out.loc[pool_idx, q_col] = out.loc[pool_idx, p_col]
    return out


def stream_dataset_pass0(meta: dict[str, Any], g_proj: sp.csc_matrix, executor: ProcessPoolExecutor) -> dict[str, Any]:
    """Pass 0: verify count_gate and build exact post-exclusion depth-eligibility table before modeling."""
    t0 = time.perf_counter()
    h5_path = meta["h5_path"]
    slot = meta["matrix_slot"]
    indptr = meta["indptr"]
    keep_mask = meta["keep_mask"]
    c_donor = meta["cell_donor_idx"]
    c_ct = meta["cell_ct_idx"]
    n_donors = len(meta["donor_ids"])
    n_cts = len(meta["cell_types"])
    n_raw_vars = len(meta["gene_symbols"])

    cnt_mat = np.zeros((n_donors * n_cts, 4), dtype=np.int64)
    pool_mat = np.zeros((n_donors, 4), dtype=np.int64)
    log_sum_ct = np.zeros((n_donors * n_cts, 4), dtype=np.float64)
    log_sum_pool = np.zeros((n_donors, 4), dtype=np.float64)
    gate_tot = {"nonfinite": 0, "negative": 0, "noninteger": 0, "nnz": 0, "cells": 0}

    with h5py.File(h5_path, "r") as f:
        grp = f["layers"]["counts"] if slot == "layers/counts" else (f["raw"]["X"] if slot == "raw/X" else f["X"])
        n_obs = meta["n_obs"]
        futures = []

        def _drain(limit: int) -> None:
            while len(futures) > limit:
                fut = futures.pop(0)
                r = fut.result()
                for k in ("nonfinite", "negative", "noninteger"):
                    gate_tot[k] += r["gate"][k]
                gate_tot["nnz"] += r["nnz"]
                gate_tot["cells"] += r["n_cells"]
                cnt_mat[:] += r["cnt_mat"]
                pool_mat[:] += r["pool_mat"]
                log_sum_ct[:] += r["log_sum_ct"]
                log_sum_pool[:] += r["log_sum_pool"]

        for r0 in range(0, n_obs, CHUNK_ROWS):
            r1 = min(r0 + CHUNK_ROWS, n_obs)
            local_keep = keep_mask[r0:r1]
            if not np.any(local_keep):
                continue
            chunk_csr = read_csr_chunk(grp, r0, r1, n_raw_vars, indptr)
            sub_csr = chunk_csr[local_keep]
            d_sub = c_donor[r0:r1][local_keep]
            c_sub = c_ct[r0:r1][local_keep]
            del chunk_csr
            for s0 in range(0, sub_csr.shape[0], SUBCHUNK_ROWS):
                s1 = min(s0 + SUBCHUNK_ROWS, sub_csr.shape[0])
                sc = sub_csr[s0:s1]
                futures.append(executor.submit(_subchunk_pass0_worker, (sc.data, sc.indices, sc.indptr, n_raw_vars, d_sub[s0:s1], c_sub[s0:s1], n_donors, n_cts)))
                _drain(MAX_WORKERS)
        _drain(0)

    # Tabulate pre-modeling depth eligibility table
    donors_df = meta["donors_df"]
    ct_names = meta["cell_types"]
    cnt_3d = cnt_mat.reshape(n_donors, n_cts, 4)
    depth_rows = []
    thresholds = (1, 500, 1500, 3000)
    for ci, ct in enumerate(ct_names):
        row: dict[str, Any] = {"group": ct, "is_pooled": False}
        for k, thr in enumerate(thresholds):
            col_c = cnt_3d[:, ci, k]
            elig = col_c >= MIN_CELLS_PER_DONOR
            row[f"cells_ge{thr}"] = int(col_c.sum())
            row[f"donors_ge20_at_{thr}"] = int(elig.sum())
            row[f"young_donors_ge20_at_{thr}"] = int((elig & (donors_df["bin"].to_numpy() == "YOUNG")).sum())
            row[f"mid_donors_ge20_at_{thr}"] = int((elig & (donors_df["bin"].to_numpy() == "MID")).sum())
            row[f"old_donors_ge20_at_{thr}"] = int((elig & (donors_df["bin"].to_numpy() == "OLD")).sum())
        depth_rows.append(row)

    prow: dict[str, Any] = {"group": meta["pooled_name"], "is_pooled": True}
    for k, thr in enumerate(thresholds):
        col_c = pool_mat[:, k]
        elig = col_c >= MIN_CELLS_PER_DONOR
        prow[f"cells_ge{thr}"] = int(col_c.sum())
        prow[f"donors_ge20_at_{thr}"] = int(elig.sum())
        prow[f"young_donors_ge20_at_{thr}"] = int((elig & (donors_df["bin"].to_numpy() == "YOUNG")).sum())
        prow[f"mid_donors_ge20_at_{thr}"] = int((elig & (donors_df["bin"].to_numpy() == "MID")).sum())
        prow[f"old_donors_ge20_at_{thr}"] = int((elig & (donors_df["bin"].to_numpy() == "OLD")).sum())
    depth_rows.append(prow)
    depth_df = pd.DataFrame(depth_rows)

    return {
        "gate": gate_tot,
        "cnt_3d": cnt_3d,
        "pool_mat": pool_mat,
        "log_sum_ct": log_sum_ct.reshape(n_donors, n_cts, 4),
        "log_sum_pool": log_sum_pool,
        "depth_table": depth_df,
        "elapsed_s": round(time.perf_counter() - t0, 2),
    }


def _safe_mean(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    out = np.full(num.shape, np.nan, dtype=np.float64)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def _accum_csr_tuple(
    target: np.ndarray,
    csr_tuple: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None,
    row_map: dict[int, int] | None = None,
    scale: float = 1.0,
    as_int32: bool = False,
) -> None:
    """Accumulate compact sparse worker output (u_ids, data, indices, indptr) directly onto dense parent buffer."""
    if csr_tuple is None:
        return
    u_m, data, indices, indptr = csr_tuple
    for k_i, g_id in enumerate(u_m):
        g_int = int(g_id)
        dst = row_map.get(g_int, -1) if row_map is not None else g_int
        if dst < 0:
            continue
        s0, s1 = int(indptr[k_i]), int(indptr[k_i + 1])
        vals = data[s0:s1]
        if scale != 1.0:
            vals = vals * scale
        if as_int32:
            vals = np.rint(vals).astype(np.int32)
        target[dst, indices[s0:s1]] += vals


def stream_dataset_pass1(
    meta: dict[str, Any],
    p0: dict[str, Any],
    g_proj: sp.csc_matrix,
    executor: ProcessPoolExecutor,
    do_500_multiseed: bool,
    partb_full_floor: int,
    run_gao_a3_centroids: bool,
) -> dict[str, Any]:
    """Pass 1: compute all Part B sums, A1 reference vectors, donor mean_log_umi, and A2 pseudobulks."""
    t0 = time.perf_counter()
    h5_path = meta["h5_path"]
    slot = meta["matrix_slot"]
    indptr = meta["indptr"]
    keep_mask = meta["keep_mask"]
    c_donor = meta["cell_donor_idx"]
    c_ct = meta["cell_ct_idx"]
    n_donors = len(meta["donor_ids"])
    n_cts = len(meta["cell_types"])
    n_pairs = n_donors * n_cts
    n_raw_vars = len(meta["gene_symbols"])
    n_genes = g_proj.shape[1]
    ds_code = 1 if meta["dataset"] == "gao" else 2

    donors_df = meta["donors_df"]
    bin_arr = donors_df["bin"].map({"LT20": 0, "YOUNG": 1, "MID": 2, "OLD": 3}).to_numpy(dtype=np.int32)
    sh_arr = donors_df["split_half_group"].to_numpy(dtype=np.int32)

    # Pre-filter observed YOUNG/MID pairs meeting >=20 cell floor cleanly to keep int32 ref tensors lean
    cnt_elig_base = p0["cnt_3d"][:, :, 1 if ds_code == 1 else 2]
    ref_pair_ym_lookup = np.full(n_pairs, -1, dtype=np.int64)
    ref_pair_y_lookup = np.full(n_pairs, -1, dtype=np.int64)
    ref_ym_pairs: list[tuple[int, int]] = []
    ref_y_pairs: list[tuple[int, int]] = []
    for di in range(n_donors):
        b_v = int(bin_arr[di])
        if b_v not in (1, 2):
            continue
        for ci in range(n_cts):
            if cnt_elig_base[di, ci] >= MIN_CELLS_PER_DONOR:
                pid = di * n_cts + ci
                ref_pair_ym_lookup[pid] = len(ref_ym_pairs)
                ref_ym_pairs.append((di, ci))
                if b_v == 1:
                    ref_pair_y_lookup[pid] = len(ref_y_pairs)
                    ref_y_pairs.append((di, ci))

    # Gao A3 pair lookup (only pairs >= 20 cells at >= 1500 UMIs)
    gao_a3_lookup = None
    gao_a3_pairs: list[tuple[int, int]] = []
    if run_gao_a3_centroids:
        cnt_15 = p0["cnt_3d"][:, :, 2]
        gao_a3_lookup = np.full(n_pairs, -1, dtype=np.int64)
        for di in range(n_donors):
            for ci in range(n_cts):
                if cnt_15[di, ci] >= MIN_CELLS_PER_DONOR:
                    gao_a3_lookup[di * n_cts + ci] = len(gao_a3_pairs)
                    gao_a3_pairs.append((di, ci))

    n_s = len(ENTROPY_SEEDS)
    b_sums = {
        "ct_full": np.zeros((n_pairs, 3), dtype=np.float64),
        "pool_full": np.zeros((n_donors, 3), dtype=np.float64),
        "ct_full_gt0": np.zeros((n_pairs, 3), dtype=np.float64),
        "pool_full_gt0": np.zeros((n_donors, 3), dtype=np.float64),
        "ct_500": np.zeros((n_pairs, 3), dtype=np.float64),
        "pool_500": np.zeros((n_donors, 3), dtype=np.float64),
        "ct_3000": np.zeros((n_pairs, 3), dtype=np.float64),
        "pool_3000": np.zeros((n_donors, 3), dtype=np.float64),
        "ct_1500_seeds": np.zeros((n_s, n_pairs, 3), dtype=np.float64),
        "pool_1500_seeds": np.zeros((n_s, n_donors, 3), dtype=np.float64),
    }
    if do_500_multiseed:
        b_sums["ct_500_seeds"] = np.zeros((n_s, n_pairs, 3), dtype=np.float64)
        b_sums["pool_500_seeds"] = np.zeros((n_s, n_donors, 3), dtype=np.float64)

    ref_ct_tensors = {
        "ct_full_1500": np.zeros((len(ref_ym_pairs), n_genes), dtype=np.int32),
        "ct_1500": np.zeros((len(ref_ym_pairs), n_genes), dtype=np.int32),
        "ct_500_same": np.zeros((len(ref_y_pairs), n_genes), dtype=np.int32),
        "ct_3000": np.zeros((len(ref_y_pairs), n_genes), dtype=np.int32),
    }
    if ds_code == 1:
        ref_ct_tensors["ct_full_500"] = np.zeros((len(ref_y_pairs), n_genes), dtype=np.int32)

    pool_ym_donors = sorted(np.flatnonzero(np.isin(bin_arr, [1, 2])).tolist())
    pool_y_donors = sorted(np.flatnonzero(bin_arr == 1).tolist())
    ym_d_to_row = {d: i for i, d in enumerate(pool_ym_donors)}
    y_d_to_row = {d: i for i, d in enumerate(pool_y_donors)}
    ref_pool_tensors = {
        "pool_full_1500": np.zeros((len(pool_ym_donors), n_genes), dtype=np.int32),
        "pool_1500": np.zeros((len(pool_ym_donors), n_genes), dtype=np.int32),
        "pool_500_same": np.zeros((len(pool_y_donors), n_genes), dtype=np.int32),
        "pool_3000": np.zeros((len(pool_y_donors), n_genes), dtype=np.int32),
    }
    if ds_code == 1:
        ref_pool_tensors["pool_full_500"] = np.zeros((len(pool_y_donors), n_genes), dtype=np.int32)

    a2_pb_mode = "gao_all_pairs" if ds_code == 1 else "yazar_pool_only"
    a2_ct_gao = np.zeros((n_pairs, n_genes), dtype=np.float64) if ds_code == 1 else None
    a2_pool = np.zeros((n_donors, n_genes), dtype=np.float64)
    a3_ct_1500 = np.zeros((len(gao_a3_pairs), n_genes), dtype=np.float64) if run_gao_a3_centroids else None
    a3_pool_1500 = np.zeros((n_donors, n_genes), dtype=np.float64) if run_gao_a3_centroids else None
    a3_ct_full = np.zeros((len(gao_a3_pairs), n_genes), dtype=np.float64) if run_gao_a3_centroids else None
    a3_pool_full = np.zeros((n_donors, n_genes), dtype=np.float64) if run_gao_a3_centroids else None

    with h5py.File(h5_path, "r") as f:
        grp = f["layers"]["counts"] if slot == "layers/counts" else (f["raw"]["X"] if slot == "raw/X" else f["X"])
        n_obs = meta["n_obs"]
        futures = []
        subchunk_id = 0

        def _drain(limit: int) -> None:
            while len(futures) > limit:
                r = futures.pop(0).result()
                pb = r["partb"]
                u_p = pb.pop("u_pids")
                u_d = pb.pop("u_dids")
                for bk, bv in pb.items():
                    if bk.startswith("ct_"):
                        if bv.ndim == 3:
                            b_sums[bk][:, u_p, :] += bv
                        else:
                            b_sums[bk][u_p] += bv
                    else:
                        if bv.ndim == 3:
                            b_sums[bk][:, u_d, :] += bv
                        else:
                            b_sums[bk][u_d] += bv
                for rk, csr_t in r["a1_refs"].items():
                    if rk.startswith("ct_"):
                        _accum_csr_tuple(ref_ct_tensors[rk], csr_t, as_int32=True)
                    elif rk.startswith("pool_"):
                        rm = ym_d_to_row if rk in ("pool_full_1500", "pool_1500") else y_d_to_row
                        _accum_csr_tuple(ref_pool_tensors[rk], csr_t, row_map=rm, as_int32=True)
                if "ct_all" in r["a2"] and a2_ct_gao is not None:
                    _accum_csr_tuple(a2_ct_gao, r["a2"]["ct_all"])
                if "pool_all" in r["a2"]:
                    _accum_csr_tuple(a2_pool, r["a2"]["pool_all"])
                if run_gao_a3_centroids and r["a3"]:
                    _accum_csr_tuple(a3_ct_1500, r["a3"].get("ct_1500"), scale=1.0 / 1500.0)
                    _accum_csr_tuple(a3_pool_1500, r["a3"].get("pool_1500"), scale=1.0 / 1500.0)
                    _accum_csr_tuple(a3_ct_full, r["a3"].get("ct_full_l1"))
                    _accum_csr_tuple(a3_pool_full, r["a3"].get("pool_full_l1"))

        for r0 in range(0, n_obs, CHUNK_ROWS):
            r1 = min(r0 + CHUNK_ROWS, n_obs)
            local_keep = keep_mask[r0:r1]
            if not np.any(local_keep):
                continue
            chunk_csr = read_csr_chunk(grp, r0, r1, n_raw_vars, indptr)
            sub_csr = chunk_csr[local_keep]
            d_sub = c_donor[r0:r1][local_keep]
            c_sub = c_ct[r0:r1][local_keep]
            del chunk_csr
            for s0 in range(0, sub_csr.shape[0], SUBCHUNK_ROWS):
                s1 = min(s0 + SUBCHUNK_ROWS, sub_csr.shape[0])
                sc = sub_csr[s0:s1]
                futures.append(
                    executor.submit(
                        _subchunk_pass1_worker,
                        (
                            sc.data,
                            sc.indices,
                            sc.indptr,
                            n_raw_vars,
                            d_sub[s0:s1],
                            c_sub[s0:s1],
                            ds_code,
                            subchunk_id,
                            n_donors,
                            n_cts,
                            ref_pair_ym_lookup,
                            ref_pair_y_lookup,
                            sh_arr,
                            bin_arr,
                            ENTROPY_SEEDS,
                            do_500_multiseed,
                            gao_a3_lookup,
                            a2_pb_mode,
                            partb_full_floor,
                        ),
                    )
                )
                subchunk_id += 1
                _drain(MAX_WORKERS)
        _drain(0)

    # Assemble structured A1 reference dictionaries compatible with reference_without_donor(), popping tensors cleanly (int32 views)
    def _assemble_ref_dict(tensor_key: str, is_pooled: bool, cell_count_col_idx: int) -> dict[str, dict[str, Any]]:
        out_d: dict[str, dict[str, Any]] = {}
        donor_ids = meta["donor_ids"]
        ct_names = meta["cell_types"]

        def _build_subset_entry(items: list[tuple[int, int, np.ndarray]]) -> dict[str, Any]:
            per_d: dict[str, Any] = {}
            tot_cnt = np.zeros(n_genes, dtype=np.int64)
            tot_cells = 0
            for di_loc, n_c_loc, vec_loc in items:
                tot_cnt += vec_loc
                tot_cells += n_c_loc
                per_d[donor_ids[di_loc]] = {"counts": vec_loc, "n_cells": n_c_loc}
            return {"counts": tot_cnt, "n_cells": tot_cells, "donors": per_d}

        if not is_pooled:
            tens = ref_ct_tensors.pop(tensor_key)
            pair_src = ref_ym_pairs if tensor_key in ("ct_full_1500", "ct_1500") else ref_y_pairs
            c_3d = p0["cnt_3d"][:, :, cell_count_col_idx]
            for ci, ct in enumerate(ct_names):
                grp_entry: dict[str, Any] = {}
                ct_items_young: list[tuple[int, int, np.ndarray]] = []
                ct_items_mid: list[tuple[int, int, np.ndarray]] = []
                for r_idx, (di, r_ci) in enumerate(pair_src):
                    if r_ci != ci:
                        continue
                    n_c_elig = int(c_3d[di, ci])
                    if n_c_elig < MIN_CELLS_PER_DONOR:
                        continue
                    vec = tens[r_idx]
                    if bin_arr[di] == 1:
                        ct_items_young.append((di, n_c_elig, vec))
                    elif bin_arr[di] == 2:
                        ct_items_mid.append((di, n_c_elig, vec))

                if ct_items_young:
                    grp_entry["young"] = _build_subset_entry(ct_items_young)
                if ct_items_mid:
                    grp_entry["mid"] = _build_subset_entry(ct_items_mid)

                h1_items = [it for it in ct_items_young if sh_arr[it[0]] == 1]
                h2_items = [it for it in ct_items_young if sh_arr[it[0]] == 2]
                h1_cells = sum(it[1] for it in h1_items)
                fb_used = False
                if (len(h1_items) == 0 or len(h2_items) == 0 or h1_cells < MIN_LOO_REF_CELLS) and len(ct_items_young) >= MIN_YOUNG_A1:
                    ysorted = sorted(ct_items_young, key=lambda it: donor_ids[it[0]])
                    h1_cand = ysorted[::2]
                    h2_cand = ysorted[1::2]
                    if sum(it[1] for it in h1_cand) < MIN_LOO_REF_CELLS and sum(it[1] for it in h2_cand) >= MIN_LOO_REF_CELLS:
                        h1_cand, h2_cand = h2_cand, h1_cand
                    if h1_cand and h2_cand and sum(it[1] for it in h1_cand) >= MIN_LOO_REF_CELLS:
                        h1_items, h2_items = h1_cand, h2_cand
                        fb_used = True
                if h1_items and h2_items:
                    h1_entry = _build_subset_entry(h1_items)
                    h1_entry["h1_n_donors"] = len(h1_items)
                    h1_entry["h1_n_cells"] = int(h1_entry["n_cells"])
                    h1_entry["h2_donors"] = {donor_ids[it[0]] for it in h2_items}
                    h1_entry["fallback_used"] = bool(fb_used)
                    grp_entry["half1"] = h1_entry
                if grp_entry:
                    out_d[ct] = grp_entry
        else:
            tens = ref_pool_tensors.pop(tensor_key)
            if tensor_key in ("pool_full_1500", "pool_1500"):
                d_list, d_map = pool_ym_donors, ym_d_to_row
            else:
                d_list, d_map = pool_y_donors, y_d_to_row
            p_cnt = p0["pool_mat"][:, cell_count_col_idx]
            pool_items_young: list[tuple[int, int, np.ndarray]] = []
            pool_items_mid: list[tuple[int, int, np.ndarray]] = []
            for di in d_list:
                n_c_elig = int(p_cnt[di])
                if n_c_elig < MIN_CELLS_PER_DONOR:
                    continue
                vec = tens[d_map[di]]
                if bin_arr[di] == 1:
                    pool_items_young.append((di, n_c_elig, vec))
                elif bin_arr[di] == 2:
                    pool_items_mid.append((di, n_c_elig, vec))
            grp_entry = {}
            if pool_items_young:
                grp_entry["young"] = _build_subset_entry(pool_items_young)
            if pool_items_mid:
                grp_entry["mid"] = _build_subset_entry(pool_items_mid)
            h1_items = [it for it in pool_items_young if sh_arr[it[0]] == 1]
            h2_items = [it for it in pool_items_young if sh_arr[it[0]] == 2]
            if h1_items and h2_items:
                h1_entry = _build_subset_entry(h1_items)
                h1_entry["h1_n_donors"] = len(h1_items)
                h1_entry["h1_n_cells"] = int(h1_entry["n_cells"])
                h1_entry["h2_donors"] = {donor_ids[it[0]] for it in h2_items}
                h1_entry["fallback_used"] = False
                grp_entry["half1"] = h1_entry
            if grp_entry:
                out_d["_POOLED_"] = grp_entry
        return out_d

    a1_refs = {
        "donor_ids": meta["donor_ids"],
        "cell_types": meta["cell_types"],
        "split_half": sh_arr,
        "age_bin": bin_arr,
        "ct_full_1500": _assemble_ref_dict("ct_full_1500", False, 2),
        "ct_1500": _assemble_ref_dict("ct_1500", False, 2),
        "ct_500_same": _assemble_ref_dict("ct_500_same", False, 2),
        "ct_3000": _assemble_ref_dict("ct_3000", False, 3),
        "pool_full_1500": _assemble_ref_dict("pool_full_1500", True, 2),
        "pool_1500": _assemble_ref_dict("pool_1500", True, 2),
        "pool_500_same": _assemble_ref_dict("pool_500_same", True, 2),
        "pool_3000": _assemble_ref_dict("pool_3000", True, 3),
    }
    if ds_code == 1:
        a1_refs["ct_full_500"] = _assemble_ref_dict("ct_full_500", False, 1)
        a1_refs["pool_full_500"] = _assemble_ref_dict("pool_full_500", True, 1)

    del ref_ct_tensors, ref_pool_tensors
    gc.collect()
    return {
        "partb": b_sums,
        "a1_refs": a1_refs,
        "a2_ct_gao": a2_ct_gao,
        "a2_pool": a2_pool,
        "a3_gao_pairs": gao_a3_pairs,
        "a3_gao_lookup": gao_a3_lookup,
        "a3_ct_1500": a3_ct_1500,
        "a3_pool_1500": a3_pool_1500,
        "a3_ct_full": a3_ct_full,
        "a3_pool_full": a3_pool_full,
        "elapsed_s": round(time.perf_counter() - t0, 2),
    }


def _collect_eligible_a1_masks(p0: dict[str, Any], refs: dict[str, Any], bin_arr: np.ndarray, n_donors: int, n_cts: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Enforce >=20 cells/donor and A1 reference eligibility (>=2 YOUNG & min LOO remaining >=50 cells)."""
    elig_reports: list[dict[str, Any]] = []
    masks: dict[str, Any] = {}
    for arm_key, col_idx, is_pool in (
        ("ct_1500", 2, False),
        ("ct_3000", 3, False),
        ("ct_500", 1, False),
        ("pool_1500", 2, True),
        ("pool_3000", 3, True),
        ("pool_500", 1, True),
    ):
        ref_map = refs.get(f"{'pool' if is_pool else 'ct'}_{'full_1500' if col_idx == 2 else ('3000' if col_idx == 3 else 'full_500')}", {})
        if not is_pool:
            m_arr = np.full(n_donors * n_cts, -1, dtype=np.int64)
            c_3d = p0["cnt_3d"][:, :, col_idx]
            for ci, ct in enumerate(refs["cell_types"]):
                elig_d = np.flatnonzero(c_3d[:, ci] >= MIN_CELLS_PER_DONOR)
                y_donors = [refs["donor_ids"][d] for d in elig_d if bin_arr[d] == 1]
                r_obj = ref_map.get(ct, {}).get("young")
                min_loo = 0
                if r_obj is not None and len(y_donors) >= MIN_YOUNG_A1:
                    min_loo = min(int(r_obj["n_cells"]) - int(r_obj["donors"].get(yd, {}).get("n_cells", 0)) for yd in y_donors)
                ok = len(elig_d) >= MIN_DONORS_SLOPE and len(y_donors) >= MIN_YOUNG_A1 and min_loo >= MIN_LOO_REF_CELLS
                if ok:
                    for di in elig_d:
                        m_arr[di * n_cts + ci] = di * n_cts + ci
                if col_idx in (2, 3):
                    elig_reports.append(
                        {
                            "part": "A1",
                            "depth": "1500" if col_idx == 2 else "3000",
                            "group": ct,
                            "eligible": bool(ok),
                            "n_donors": int(len(elig_d)),
                            "n_young": int(len(y_donors)),
                            "min_loo_ref_cells": int(min_loo),
                            "reason": "eligible" if ok else "failed_min_donors_or_young_or_loo50",
                        }
                    )
            masks[arm_key] = m_arr
        else:
            p_arr = np.full(n_donors, -1, dtype=np.int64)
            p_cnt = p0["pool_mat"][:, col_idx]
            elig_d = np.flatnonzero(p_cnt >= MIN_CELLS_PER_DONOR)
            for di in elig_d:
                p_arr[di] = di
            masks[arm_key] = p_arr
    return masks, elig_reports


def stream_dataset_pass2(
    meta: dict[str, Any],
    p0: dict[str, Any],
    p1: dict[str, Any],
    g_proj: sp.csc_matrix,
    run_gao_500_a1: bool,
) -> dict[str, Any]:
    """Pass 2: stream HDF5 for A1 distances (and A3 LOO centroids/distances via clean donor batches)."""
    t0 = time.perf_counter()
    h5_path = meta["h5_path"]
    slot = meta["matrix_slot"]
    indptr = meta["indptr"]
    keep_mask = meta["keep_mask"]
    c_donor = meta["cell_donor_idx"]
    c_ct = meta["cell_ct_idx"]
    n_donors = len(meta["donor_ids"])
    n_cts = len(meta["cell_types"])
    n_pairs = n_donors * n_cts
    n_raw_vars = len(meta["gene_symbols"])
    n_genes = g_proj.shape[1]
    ds_code = 1 if meta["dataset"] == "gao" else 2
    donors_df = meta["donors_df"]
    bin_arr = donors_df["bin"].map({"LT20": 0, "YOUNG": 1, "MID": 2, "OLD": 3}).to_numpy(dtype=np.int32)

    a1_masks, a1_elig_log = _collect_eligible_a1_masks(p0, p1["a1_refs"], bin_arr, n_donors, n_cts)
    worker_refs = {**p1["a1_refs"], "run_gao_500_a1": run_gao_500_a1}

    a1_sums: dict[str, dict[int, float]] = {}
    a1_cnts: dict[str, dict[int, int]] = {}
    a3_sums: dict[str, dict[int, float]] = {}
    a3_cnts: dict[str, dict[int, int]] = {}

    cnt_15 = p0["cnt_3d"][:, :, 2]
    pool_15 = p0["pool_mat"][:, 2]

    # Gao executes Pass 2 (A1 + A3 LOO distances) in 1 single HDF5 pass using Pass-1 A3 centroids
    if ds_code == 1:
        worker_refs["a3_ct_1500"] = p1["a3_ct_1500"]
        worker_refs["a3_pool_1500"] = p1["a3_pool_1500"]
        worker_refs["a3_ct_full"] = p1["a3_ct_full"]
        worker_refs["a3_pool_full"] = p1["a3_pool_full"]
        a3_n_ct = np.array([int(cnt_15[d, c]) for d, c in p1["a3_gao_pairs"]], dtype=np.int64)
        a3_n_pool = pool_15.astype(np.int64)
        a3_pool_lk = np.where(a3_n_pool >= MIN_CELLS_PER_DONOR, np.arange(n_donors, dtype=np.int64), -1)

        with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, worker_refs)) as ex:
            with h5py.File(h5_path, "r") as f:
                grp = f["layers"]["counts"]
                futures = []
                sub_id = 0

                def _drain(lim: int) -> None:
                    while len(futures) > lim:
                        r = futures.pop(0).result()
                        for (g_idx, ak), (sm_v, ct_v) in r["a1"].items():
                            a1_sums.setdefault(ak, {})[g_idx] = a1_sums.get(ak, {}).get(g_idx, 0.0) + sm_v
                            a1_cnts.setdefault(ak, {})[g_idx] = a1_cnts.get(ak, {}).get(g_idx, 0) + ct_v
                        for a3k, mp in r["a3_dist"].items():
                            for gid, (sm_v, ct_v) in mp.items():
                                a3_sums.setdefault(a3k, {})[gid] = a3_sums.get(a3k, {}).get(gid, 0.0) + sm_v
                                a3_cnts.setdefault(a3k, {})[gid] = a3_cnts.get(a3k, {}).get(gid, 0) + ct_v

                for r0 in range(0, meta["n_obs"], CHUNK_ROWS):
                    r1 = min(r0 + CHUNK_ROWS, meta["n_obs"])
                    lk = keep_mask[r0:r1]
                    if not np.any(lk):
                        continue
                    ch = read_csr_chunk(grp, r0, r1, n_raw_vars, indptr)
                    sc_all = ch[lk]
                    d_s = c_donor[r0:r1][lk]
                    c_s = c_ct[r0:r1][lk]
                    del ch
                    for s0 in range(0, sc_all.shape[0], SUBCHUNK_ROWS):
                        s1 = min(s0 + SUBCHUNK_ROWS, sc_all.shape[0])
                        sc = sc_all[s0:s1]
                        futures.append(
                            ex.submit(
                                _subchunk_pass2_worker,
                                (
                                    sc.data,
                                    sc.indices,
                                    sc.indptr,
                                    n_raw_vars,
                                    d_s[s0:s1],
                                    c_s[s0:s1],
                                    ds_code,
                                    sub_id,
                                    n_donors,
                                    n_cts,
                                    a1_masks["ct_1500"],
                                    a1_masks["ct_3000"],
                                    a1_masks["ct_500"],
                                    a1_masks["pool_1500"],
                                    a1_masks["pool_3000"],
                                    a1_masks["pool_500"],
                                    "with_a1",
                                    p1["a3_gao_lookup"],
                                    a3_pool_lk,
                                    0,
                                    a3_n_ct,
                                    a3_n_pool,
                                ),
                            )
                        )
                        sub_id += 1
                        _drain(MAX_WORKERS)
                _drain(0)

        # Map Gao A3 compressed IDs back to canonical pids for ct
        a3_ct_sums_1500 = {p1["a3_gao_pairs"][gid][0] * n_cts + p1["a3_gao_pairs"][gid][1]: v for gid, v in a3_sums.get("ct_1500", {}).items()}
        a3_ct_cnts_1500 = {p1["a3_gao_pairs"][gid][0] * n_cts + p1["a3_gao_pairs"][gid][1]: v for gid, v in a3_cnts.get("ct_1500", {}).items()}
        a3_ct_sums_full = {p1["a3_gao_pairs"][gid][0] * n_cts + p1["a3_gao_pairs"][gid][1]: v for gid, v in a3_sums.get("ct_full", {}).items()}
        a3_ct_cnts_full = {p1["a3_gao_pairs"][gid][0] * n_cts + p1["a3_gao_pairs"][gid][1]: v for gid, v in a3_cnts.get("ct_full", {}).items()}
        a3_sums["ct_1500"], a3_cnts["ct_1500"] = a3_ct_sums_1500, a3_ct_cnts_1500
        a3_sums["ct_full"], a3_cnts["ct_full"] = a3_ct_sums_full, a3_ct_cnts_full
    else:
        # Yazar Pass 2a: A1 distances cleanly in 1 stream (< 2 GB RAM)
        with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, worker_refs)) as ex:
            with h5py.File(h5_path, "r") as f:
                grp = f["raw"]["X"] if slot == "raw/X" else f["X"]
                futures = []
                sub_id = 0

                def _drain_a1(lim: int) -> None:
                    while len(futures) > lim:
                        r = futures.pop(0).result()
                        for (g_idx, ak), (sm_v, ct_v) in r["a1"].items():
                            a1_sums.setdefault(ak, {})[g_idx] = a1_sums.get(ak, {}).get(g_idx, 0.0) + sm_v
                            a1_cnts.setdefault(ak, {})[g_idx] = a1_cnts.get(ak, {}).get(g_idx, 0) + ct_v

                for r0 in range(0, meta["n_obs"], CHUNK_ROWS):
                    r1 = min(r0 + CHUNK_ROWS, meta["n_obs"])
                    lk = keep_mask[r0:r1]
                    if not np.any(lk):
                        continue
                    ch = read_csr_chunk(grp, r0, r1, n_raw_vars, indptr)
                    sc_all = ch[lk]
                    d_s = c_donor[r0:r1][lk]
                    c_s = c_ct[r0:r1][lk]
                    del ch
                    for s0 in range(0, sc_all.shape[0], SUBCHUNK_ROWS):
                        s1 = min(s0 + SUBCHUNK_ROWS, sc_all.shape[0])
                        sc = sc_all[s0:s1]
                        futures.append(
                            ex.submit(
                                _subchunk_pass2_worker,
                                (
                                    sc.data,
                                    sc.indices,
                                    sc.indptr,
                                    n_raw_vars,
                                    d_s[s0:s1],
                                    c_s[s0:s1],
                                    ds_code,
                                    sub_id,
                                    n_donors,
                                    n_cts,
                                    a1_masks["ct_1500"],
                                    a1_masks["ct_3000"],
                                    a1_masks["ct_500"],
                                    a1_masks["pool_1500"],
                                    a1_masks["pool_3000"],
                                    a1_masks["pool_500"],
                                    "a1_only",
                                    None,
                                    None,
                                    0,
                                    None,
                                    None,
                                ),
                            )
                        )
                        sub_id += 1
                        _drain_a1(MAX_WORKERS)
                _drain_a1(0)

        # Yazar Pass 2b: A3 (1500 + full L1) in 6 ~164-donor batches (< 1 GB centroids/batch, < 3.5 GB total RSS)
        n_batches = 6
        batch_bounds = np.array_split(np.arange(n_donors), n_batches)
        for bid, b_donors in enumerate(batch_bounds):
            b_d_set = set(b_donors.tolist())
            b_d_map = {d: i for i, d in enumerate(b_donors.tolist())}
            b_pool_lk = np.array([b_d_map[d] if (d in b_d_map and pool_15[d] >= MIN_CELLS_PER_DONOR) else -1 for d in range(n_donors)], dtype=np.int64)
            b_pool_nc = np.array([int(pool_15[d]) for d in b_donors.tolist()], dtype=np.int64)

            b_pairs = [(d, c) for d in b_donors.tolist() for c in range(n_cts) if cnt_15[d, c] >= MIN_CELLS_PER_DONOR]
            b_ct_lk = np.full(n_pairs, -1, dtype=np.int64)
            for idx_b, (d, c) in enumerate(b_pairs):
                b_ct_lk[d * n_cts + c] = idx_b
            b_ct_nc = np.array([int(cnt_15[d, c]) for d, c in b_pairs], dtype=np.int64)

            b_cent_ct_15 = np.zeros((len(b_pairs), n_genes), dtype=np.float64)
            b_cent_pool_15 = np.zeros((len(b_donors), n_genes), dtype=np.float64)
            b_cent_ct_fl = np.zeros((len(b_pairs), n_genes), dtype=np.float64)
            b_cent_pool_fl = np.zeros((len(b_donors), n_genes), dtype=np.float64)

            batch_row_mask = keep_mask.copy()
            d_in_batch = np.isin(c_donor, list(b_d_set))
            batch_row_mask = batch_row_mask & d_in_batch

            with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, {})) as ex:
                with h5py.File(h5_path, "r") as f:
                    grp = f["raw"]["X"] if slot == "raw/X" else f["X"]
                    futures = []
                    sub_id = 0

                    def _drain_cent(lim: int) -> None:
                        while len(futures) > lim:
                            r = futures.pop(0).result()
                            _accum_csr_tuple(b_cent_ct_15, r["a3_centroids"].get("ct_1500"), scale=1.0 / 1500.0)
                            _accum_csr_tuple(b_cent_pool_15, r["a3_centroids"].get("pool_1500"), scale=1.0 / 1500.0)
                            _accum_csr_tuple(b_cent_ct_fl, r["a3_centroids"].get("ct_full_l1"))
                            _accum_csr_tuple(b_cent_pool_fl, r["a3_centroids"].get("pool_full_l1"))

                    for r0 in range(0, meta["n_obs"], CHUNK_ROWS):
                        r1 = min(r0 + CHUNK_ROWS, meta["n_obs"])
                        lk = batch_row_mask[r0:r1]
                        if not np.any(lk):
                            continue
                        ch = read_csr_chunk(grp, r0, r1, n_raw_vars, indptr)
                        sc_all = ch[lk]
                        d_s = c_donor[r0:r1][lk]
                        c_s = c_ct[r0:r1][lk]
                        del ch
                        for s0 in range(0, sc_all.shape[0], SUBCHUNK_ROWS):
                            s1 = min(s0 + SUBCHUNK_ROWS, sc_all.shape[0])
                            sc = sc_all[s0:s1]
                            futures.append(
                                ex.submit(
                                    _subchunk_pass2_worker,
                                    (
                                        sc.data,
                                        sc.indices,
                                        sc.indptr,
                                        n_raw_vars,
                                        d_s[s0:s1],
                                        c_s[s0:s1],
                                        ds_code,
                                        sub_id,
                                        n_donors,
                                        n_cts,
                                        None,
                                        None,
                                        None,
                                        None,
                                        None,
                                        None,
                                        "centroids_only",
                                        b_ct_lk,
                                        b_pool_lk,
                                        bid + 1,
                                        None,
                                        None,
                                    ),
                                )
                            )
                            sub_id += 1
                            _drain_cent(MAX_WORKERS)
                    _drain_cent(0)

            cent_refs = {
                "a3_ct_1500": b_cent_ct_15,
                "a3_pool_1500": b_cent_pool_15,
                "a3_ct_full": b_cent_ct_fl,
                "a3_pool_full": b_cent_pool_fl,
            }
            with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, cent_refs)) as ex:
                with h5py.File(h5_path, "r") as f:
                    grp = f["raw"]["X"] if slot == "raw/X" else f["X"]
                    futures = []
                    sub_id = 0

                    def _drain_dist(lim: int) -> None:
                        while len(futures) > lim:
                            r = futures.pop(0).result()
                            for gid, (sm_v, ct_v) in r["a3_dist"].get("ct_1500", {}).items():
                                pid = b_pairs[gid][0] * n_cts + b_pairs[gid][1]
                                a3_sums.setdefault("ct_1500", {})[pid] = a3_sums.get("ct_1500", {}).get(pid, 0.0) + sm_v
                                a3_cnts.setdefault("ct_1500", {})[pid] = a3_cnts.get("ct_1500", {}).get(pid, 0) + ct_v
                            for did_l, (sm_v, ct_v) in r["a3_dist"].get("pool_1500", {}).items():
                                did_g = int(b_donors[did_l])
                                a3_sums.setdefault("pool_1500", {})[did_g] = a3_sums.get("pool_1500", {}).get(did_g, 0.0) + sm_v
                                a3_cnts.setdefault("pool_1500", {})[did_g] = a3_cnts.get("pool_1500", {}).get(did_g, 0) + ct_v
                            for gid, (sm_v, ct_v) in r["a3_dist"].get("ct_full", {}).items():
                                pid = b_pairs[gid][0] * n_cts + b_pairs[gid][1]
                                a3_sums.setdefault("ct_full", {})[pid] = a3_sums.get("ct_full", {}).get(pid, 0.0) + sm_v
                                a3_cnts.setdefault("ct_full", {})[pid] = a3_cnts.get("ct_full", {}).get(pid, 0) + ct_v
                            for did_l, (sm_v, ct_v) in r["a3_dist"].get("pool_full", {}).items():
                                did_g = int(b_donors[did_l])
                                a3_sums.setdefault("pool_full", {})[did_g] = a3_sums.get("pool_full", {}).get(did_g, 0.0) + sm_v
                                a3_cnts.setdefault("pool_full", {})[did_g] = a3_cnts.get("pool_full", {}).get(did_g, 0) + ct_v

                    for r0 in range(0, meta["n_obs"], CHUNK_ROWS):
                        r1 = min(r0 + CHUNK_ROWS, meta["n_obs"])
                        lk = batch_row_mask[r0:r1]
                        if not np.any(lk):
                            continue
                        ch = read_csr_chunk(grp, r0, r1, n_raw_vars, indptr)
                        sc_all = ch[lk]
                        d_s = c_donor[r0:r1][lk]
                        c_s = c_ct[r0:r1][lk]
                        del ch
                        for s0 in range(0, sc_all.shape[0], SUBCHUNK_ROWS):
                            s1 = min(s0 + SUBCHUNK_ROWS, sc_all.shape[0])
                            sc = sc_all[s0:s1]
                            futures.append(
                                ex.submit(
                                    _subchunk_pass2_worker,
                                    (
                                        sc.data,
                                        sc.indices,
                                        sc.indptr,
                                        n_raw_vars,
                                        d_s[s0:s1],
                                        c_s[s0:s1],
                                        ds_code,
                                        sub_id,
                                        n_donors,
                                        n_cts,
                                        None,
                                        None,
                                        None,
                                        None,
                                        None,
                                        None,
                                        "distances_only",
                                        b_ct_lk,
                                        b_pool_lk,
                                        bid + 1,
                                        b_ct_nc,
                                        b_pool_nc,
                                    ),
                                )
                            )
                            sub_id += 1
                            _drain_dist(MAX_WORKERS)
                    _drain_dist(0)

            del b_cent_ct_15, b_cent_pool_15, b_cent_ct_fl, b_cent_pool_fl
            gc.collect()

    return {
        "a1_sums": a1_sums,
        "a1_cnts": a1_cnts,
        "a1_masks": a1_masks,
        "a1_elig_log": a1_elig_log,
        "a3_sums": a3_sums,
        "a3_cnts": a3_cnts,
        "elapsed_s": round(time.perf_counter() - t0, 2),
    }


def _build_donors_table_and_partb_cache(
    meta: dict[str, Any],
    p0: dict[str, Any],
    p1: dict[str, Any],
    p2: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build exact Part A long donor-by-group DataFrame and compact Part B donor metric aggregates."""
    donors_df = meta["donors_df"]
    donor_ids = meta["donor_ids"]
    ct_names = meta["cell_types"]
    n_donors = len(donor_ids)
    n_cts = len(ct_names)
    ds_name = meta["dataset"]
    pooled_name = meta["pooled_name"]
    cnt_3d = p0["cnt_3d"]
    pool_mat = p0["pool_mat"]
    mlu_ct = _safe_mean(p0["log_sum_ct"], cnt_3d)
    mlu_pool = _safe_mean(p0["log_sum_pool"], pool_mat)

    def _mean_map(sdict: dict[int, float], cdict: dict[int, int], idx: int) -> float:
        c = cdict.get(idx, 0)
        return float(sdict[idx] / c) if c > 0 else np.nan

    rows: list[dict[str, Any]] = []
    for gi, grp in enumerate([*ct_names, pooled_name]):
        is_pool = grp == pooled_name
        for di, did in enumerate(donor_ids):
            dr = donors_df.iloc[di]
            idx = di if is_pool else (di * n_cts + gi)
            n_all = int(pool_mat[di, 0] if is_pool else cnt_3d[di, gi, 0])
            n_500 = int(pool_mat[di, 1] if is_pool else cnt_3d[di, gi, 1])
            n_1500 = int(pool_mat[di, 2] if is_pool else cnt_3d[di, gi, 2])
            n_3000 = int(pool_mat[di, 3] if is_pool else cnt_3d[di, gi, 3])
            if n_all == 0:
                continue

            suffix = "pool" if is_pool else "ct"
            dist_fl = _mean_map(p2["a1_sums"].get(f"a1_full_1500_{suffix}_y_loo", {}), p2["a1_cnts"].get(f"a1_full_1500_{suffix}_y_loo", {}), idx)
            dist_15 = _mean_map(p2["a1_sums"].get(f"a1_1500_{suffix}_y_loo", {}), p2["a1_cnts"].get(f"a1_1500_{suffix}_y_loo", {}), idx)
            dist_5 = _mean_map(p2["a1_sums"].get(f"a1_500_{suffix}_y_loo", {}), p2["a1_cnts"].get(f"a1_500_{suffix}_y_loo", {}), idx)
            dist_3k = _mean_map(p2["a1_sums"].get(f"a1_3000_{suffix}_y_loo", {}), p2["a1_cnts"].get(f"a1_3000_{suffix}_y_loo", {}), idx)
            dist_mid_fl = _mean_map(p2["a1_sums"].get(f"a1_full_1500_{suffix}_mid_loo", {}), p2["a1_cnts"].get(f"a1_full_1500_{suffix}_mid_loo", {}), idx)
            dist_mid_15 = _mean_map(p2["a1_sums"].get(f"a1_1500_{suffix}_mid_loo", {}), p2["a1_cnts"].get(f"a1_1500_{suffix}_mid_loo", {}), idx)
            dist_self = _mean_map(p2["a1_sums"].get(f"a1_full_1500_{suffix}_y_self_nonloo", {}), p2["a1_cnts"].get(f"a1_full_1500_{suffix}_y_self_nonloo", {}), idx)
            dist_sh_h2 = _mean_map(p2["a1_sums"].get(f"a1_full_1500_{suffix}_h2_to_h1", {}), p2["a1_cnts"].get(f"a1_full_1500_{suffix}_h2_to_h1", {}), idx)
            dist_5_fl = _mean_map(p2["a1_sums"].get(f"a1_full_500_{suffix}_y_loo", {}), p2["a1_cnts"].get(f"a1_full_500_{suffix}_y_loo", {}), idx)

            a3_15 = _mean_map(p2["a3_sums"].get("pool_1500" if is_pool else "ct_1500", {}), p2["a3_cnts"].get("pool_1500" if is_pool else "ct_1500", {}), idx)
            a3_fl = _mean_map(p2["a3_sums"].get("pool_full" if is_pool else "ct_full", {}), p2["a3_cnts"].get("pool_full" if is_pool else "ct_full", {}), idx)

            elig_15 = n_1500 >= MIN_CELLS_PER_DONOR
            elig_3k = n_3000 >= MIN_CELLS_PER_DONOR
            elig_all = n_all >= MIN_CELLS_PER_DONOR
            reasons = []
            if not elig_all:
                reasons.append("lt20_all_cells")
            elif not elig_15:
                reasons.append("lt20_1500_cells")
            elif not elig_3k:
                reasons.append("lt20_3000_cells")

            rec = {
                "dataset": ds_name,
                "donor_id": str(did),
                "age_years": float(dr["age_years"]),
                "age_per_decade": float(dr["age_per_decade"]),
                "bin": str(dr["bin"]),
                "sex": str(dr["sex"]),
                "sex_male": int(dr["sex_male"]),
                "batch_label": "study" if ds_name == "gao" else "pool_number",
                "study_or_pool": str(dr["study"] if ds_name == "gao" else int(dr["pool_number"])),
                "technology": str(dr["technology"]) if ds_name == "gao" else "10x 3' v2",
                "split_half_group": int(dr["split_half_group"]),
                "group": grp,
                "is_pooled": bool(is_pool),
                "n_cells_all": n_all,
                "n_cells_ge500": n_500,
                "n_cells_ge1500": n_1500,
                "n_cells_ge3000": n_3000,
                "mean_log_umi_all": float(mlu_pool[di, 0] if is_pool else mlu_ct[di, gi, 0]),
                "mean_log_umi_ge500": float(mlu_pool[di, 1] if is_pool else mlu_ct[di, gi, 1]),
                "mean_log_umi_ge1500": float(mlu_pool[di, 2] if is_pool else mlu_ct[di, gi, 2]),
                "mean_log_umi_ge3000": float(mlu_pool[di, 3] if is_pool else mlu_ct[di, gi, 3]),
                "a1_dist_full_1500": dist_fl,
                "a1_dist_1500": dist_15,
                "a1_dist_500": dist_5,
                "a1_dist_3000": dist_3k,
                "a1_dist_full_500_sens": dist_5_fl,
                "a1_dist_mature_ref_mid_loo": dist_mid_fl,
                "a1_dist_mature_ref_1500": dist_mid_15,
                "a1_self_non_loo": dist_self,
                "a1_split_half_h2_to_h1": dist_sh_h2,
                "a3_disp_1500": a3_15,
                "a3_disp_full_1500": a3_fl,
                "eligible_a1_1500": bool(elig_15 and np.isfinite(dist_fl)),
                "eligible_a1_3000": bool(elig_3k and np.isfinite(dist_3k)),
                "eligible_a1_young_ref": bool(elig_15 and str(dr["bin"]) == "YOUNG"),
                "eligible_a1_mid_ref": bool(elig_15 and str(dr["bin"]) == "MID"),
                "eligible_a2_oy": bool(elig_all and str(dr["bin"]) in ("YOUNG", "OLD")),
                "eligible_a2_my": bool(elig_all and str(dr["bin"]) in ("YOUNG", "MID")),
                "eligible_a3_1500": bool(elig_15 and np.isfinite(a3_15)),
                "yazar_pool_fe_all_evaluable": bool(grp not in meta.get("not_evaluable_all_umi_csv", [])) if ds_name == "yazar" else True,
                "exclusion_reason": ";".join(reasons) if reasons else "none",
            }
            if ds_name == "gao":
                for col in STUDY_DUMMY_MAP.values():
                    rec[col] = float(dr[col])
            rows.append(rec)

    donors_long = pd.DataFrame(rows)
    # Depth-specific A1 split-half & self-non-LOO donor means so every depth evaluates against its matched depth baseline
    a1_depth_baselines: dict[tuple[str, str], dict[str, float]] = {}
    depth_key_map = (
        ("full_1500", "a1_full_1500", 2),
        ("1500", "a1_1500", 2),
        ("500", "a1_500", 2),
        ("3000", "a1_3000", 3),
        ("full_500_sens", "a1_full_500", 1),
    )
    for gi, grp in enumerate([*ct_names, pooled_name]):
        is_pool = grp == pooled_name
        suffix = "pool" if is_pool else "ct"
        for d_lbl, arm_prefix, cnt_k in depth_key_map:
            sh_vals = []
            self_vals = []
            for di in range(n_donors):
                n_c_req = int(pool_mat[di, cnt_k] if is_pool else cnt_3d[di, gi, cnt_k])
                if n_c_req < MIN_CELLS_PER_DONOR:
                    continue
                idx = di if is_pool else (di * n_cts + gi)
                v_sh = _mean_map(p2["a1_sums"].get(f"{arm_prefix}_{suffix}_h2_to_h1", {}), p2["a1_cnts"].get(f"{arm_prefix}_{suffix}_h2_to_h1", {}), idx)
                if np.isfinite(v_sh):
                    sh_vals.append(v_sh)
                if str(donors_df.iloc[di]["bin"]) == "YOUNG":
                    v_self = _mean_map(p2["a1_sums"].get(f"{arm_prefix}_{suffix}_y_self_nonloo", {}), p2["a1_cnts"].get(f"{arm_prefix}_{suffix}_y_self_nonloo", {}), idx)
                    if np.isfinite(v_self):
                        self_vals.append(v_self)
            a1_depth_baselines[(d_lbl, grp)] = {
                "split_half_young_baseline": float(np.mean(sh_vals)) if sh_vals else np.nan,
                "self_baseline_biased_low": float(np.mean(self_vals)) if self_vals else np.nan,
            }

    partb_ckpt = {
        "b_sums": p1["partb"],
        "cnt_3d": cnt_3d,
        "pool_mat": pool_mat,
        "mlu_ct": mlu_ct,
        "mlu_pool": mlu_pool,
        "cell_types": ct_names,
        "pooled_name": pooled_name,
        "donor_ids": donor_ids,
        "a1_depth_baselines": a1_depth_baselines,
    }
    return donors_long, partb_ckpt


def _fit_stratified_global_row(
    strat_df: pd.DataFrame,
    elig_cts: list[str],
    out_col: str,
    cell_col: str,
    depth_col: str,
    ds_name: str,
    part: str,
    metric: str,
    depth: str,
    ref_type: str,
    include_part_a_cols: bool,
) -> dict[str, Any] | None:
    """Fit donor-clustered stratified OLS across z-scored eligible cell types."""
    if len(elig_cts) < 2 or strat_df.empty:
        return None
    parts = []
    for _, cpart in strat_df[strat_df["group"].isin(elig_cts)].dropna(subset=[out_col]).groupby("group"):
        if len(cpart) < MIN_DONORS_SLOPE or cpart["age_per_decade"].nunique() < 2:
            continue
        sd_y = float(cpart[out_col].std(ddof=1))
        if not np.isfinite(sd_y) or sd_y <= 0.0:
            continue
        cp = cpart.copy()
        cp["y_std"] = (cp[out_col] - cp[out_col].mean()) / sd_y
        cp["sex_c"] = cp["sex_male"] - cp["sex_male"].mean()
        cp["depth_c"] = cp[depth_col] - cp[depth_col].mean()
        parts.append(cp)
    if len(parts) < 2:
        return None
    df_st = pd.concat(parts, ignore_index=True)
    if df_st["donor_id"].nunique() < MIN_DONORS_SLOPE:
        return None
    ct_dum = pd.get_dummies(df_st["group"], prefix="ct", drop_first=True, dtype=float)
    cov_blocks = [df_st[["age_per_decade"]]]
    cov_names = []
    # Rank-check centered nuisance covariates & CT dummies jointly so clustered OLS is always full-rank
    design_cols = [np.ones(len(df_st), dtype=float), df_st["age_per_decade"].to_numpy(dtype=float)]
    for cand_col in ("sex_c", "depth_c"):
        vals = df_st[cand_col].to_numpy(dtype=float)
        if np.unique(vals).size >= 2 and np.linalg.matrix_rank(np.column_stack([*design_cols, vals])) == len(design_cols) + 1:
            design_cols.append(vals)
            cov_blocks.append(df_st[[cand_col]])
            cov_names.append(cand_col)
    cov_blocks.append(ct_dum)
    cov_names.append("cell_type_dummies")
    X_st = sm.add_constant(pd.concat(cov_blocks, axis=1), has_constant="add")
    fit_st = sm.OLS(df_st["y_std"], X_st).fit(cov_type="cluster", cov_kwds={"groups": df_st["donor_id"]}, use_t=True)
    df_c = max(1, int(df_st["donor_id"].nunique() - 1))
    t975_c = float(scipy_stats.t.ppf(0.975, df_c))
    t80_c = float(scipy_stats.t.ppf(0.80, df_c))
    sl = float(fit_st.params["age_per_decade"])
    se = float(fit_st.bse["age_per_decade"])
    row: dict[str, Any] = {
        "dataset": ds_name,
        "part": part,
        "metric": metric,
        "depth": depth,
        "reference_type": ref_type,
        "model_spec": "stratified_global_primary",
        "group": "Stratified_CellType_Global",
        "is_pooled": True,
        "is_stratified_global": True,
        "eligibility_status": "sensitivity_stratified_global",
        "yazar_pool_fe_status": "not_applicable",
        "n_donors": int(df_st["donor_id"].nunique()),
        "n_cells": int(df_st[cell_col].sum()),
        "residual_df": df_c,
        "covariates_active": ";".join(cov_names),
        "covariates_omitted": "none",
        "outcome_sd": 1.0,
        "resid_sd": float(np.sqrt(fit_st.mse_resid)),
        "slope_raw_per_decade": sl,
        "se_raw": se,
        "ci95_low_raw": sl - t975_c * se,
        "ci95_high_raw": sl + t975_c * se,
        "slope_sd_per_decade": sl,
        "se_sd": se,
        "ci95_low_sd": sl - t975_c * se,
        "ci95_high_sd": sl + t975_c * se,
        "pvalue": float(fit_st.pvalues["age_per_decade"]),
        "q_bh_within_slice": np.nan,
        "family_size_m": np.nan,
        "mde80_raw_per_decade": (t975_c + t80_c) * se,
        "mde80_sd_per_decade": (t975_c + t80_c) * se,
        "batch_delta_r2": np.nan,
        "batch_partial_r2": np.nan,
    }
    if include_part_a_cols:
        row.update(
            {
                "split_half_young_baseline": np.nan,
                "split_half_n_donors_h1": 0,
                "split_half_n_cells_h1": 0,
                "split_half_n_donors_h2": 0,
                "split_half_fallback_flag": False,
                "self_baseline_biased_low": np.nan,
                "young_loo_mean": np.nan,
                "mid_mean": np.nan,
                "old_mean": np.nan,
                "old_minus_young_diff_adj": np.nan,
                "old_minus_young_se_adj": np.nan,
                "old_minus_young_ci95_low_adj": np.nan,
                "old_minus_young_ci95_high_adj": np.nan,
                "old_minus_young_p_adj": np.nan,
                "old_minus_mid_diff_adj": np.nan,
                "old_minus_mid_se_adj": np.nan,
                "old_minus_mid_ci95_low_adj": np.nan,
                "old_minus_mid_ci95_high_adj": np.nan,
                "old_minus_mid_p_adj": np.nan,
            }
        )
    return row


def _run_a1_a3_models_from_donors(
    ds_name: str,
    donors_long: pd.DataFrame,
    refs_summary: dict[str, Any],
    part_a_models: list[dict[str, Any]],
    extra_baselines: dict[tuple[str, str], dict[str, float]] | None = None,
) -> None:
    """Fit all Part A1 and Part A3 OLS + stratified global models from donor-by-group table."""
    is_gao = ds_name == "gao"
    batch_col = "study_or_pool"
    ref_lv = "ENCODE v4 (Snyder)" if is_gao else None

    def _model_specs(base_df: pd.DataFrame) -> list[tuple[str, pd.DataFrame, list[str], str | None, int, int, str | None]]:
        if is_gao:
            d_10x = base_df[base_df["technology"] != "Dropseq"].copy()
            d_ge8 = base_df[~base_df["study_or_pool"].isin(["Kuppe 2022", "Sim 2021"])].copy()
            d_enc = base_df[base_df["study_or_pool"] == "ENCODE v4 (Snyder)"].copy()
            return [
                ("primary", base_df, ["sex_male", "_depth"], batch_col, 2, 5, ref_lv),
                ("batch_adj_no_depth", base_df, ["sex_male"], batch_col, 2, 5, ref_lv),
                ("age_only", base_df, [], None, 1, 5, None),
                ("age_plus_sex", base_df, ["sex_male"], None, 1, 5, None),
                ("10x_only_no_penn", d_10x, ["sex_male", "_depth"], batch_col, 2, 5, ref_lv),
                ("studies_ge8", d_ge8, ["sex_male", "_depth"], batch_col, 2, 5, ref_lv),
                ("encode_v4_only", d_enc, ["sex_male", "_depth"], None, 1, 5, None),
            ]
        return [
            ("primary", base_df, ["sex_male", "_depth"], batch_col, 2, 10, None),
            ("pool_fe_no_depth", base_df, ["sex_male"], batch_col, 2, 10, None),
            ("no_pool", base_df, ["sex_male", "_depth"], None, 1, 5, None),
            ("age_ge_20", base_df[base_df["age_years"] >= 20.0].copy(), ["sex_male", "_depth"], batch_col, 2, 10, None),
        ]

    # Define A1 and A3 arms
    arms = [
        ("A1", "hellinger_young_ref", "full_1500", "young_loo", "a1_dist_full_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_young_ref", "1500", "young_loo", "a1_dist_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_young_ref", "500", "young_loo", "a1_dist_500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_young_ref", "3000", "young_loo", "a1_dist_3000", "n_cells_ge3000", "mean_log_umi_ge3000"),
        ("A1", "hellinger_mature_ref", "full_1500", "mid_loo", "a1_dist_mature_ref_mid_loo", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A1", "hellinger_mature_ref", "1500", "mid_loo", "a1_dist_mature_ref_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A3", "within_donor_loo_disp", "1500", "none", "a3_disp_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
        ("A3", "within_donor_loo_disp", "full_1500", "none", "a3_disp_full_1500", "n_cells_ge1500", "mean_log_umi_ge1500"),
    ]
    if is_gao and "a1_dist_full_500_sens" in donors_long.columns and donors_long["a1_dist_full_500_sens"].notna().any():
        arms.append(("A1", "hellinger_young_ref", "full_500_sens", "young_loo", "a1_dist_full_500_sens", "n_cells_ge500", "mean_log_umi_ge500"))

    groups = sorted(donors_long["group"].unique().tolist(), key=lambda g: (g in ("All_LV_Nuclei", "All_Cells"), g))
    for part, metric, depth, ref_type, out_col, cell_col, depth_col in arms:
        for grp in groups:
            gdf = donors_long[donors_long["group"] == grp].copy()
            gdf["n_cells"] = gdf[cell_col]
            gdf["_depth"] = gdf[depth_col]
            elig_col = "eligible_a1_3000" if depth == "3000" else ("eligible_a3_1500" if part == "A3" else "eligible_a1_1500")
            if depth == "full_500_sens":
                gdf = gdf[gdf["n_cells_ge500"] >= MIN_CELLS_PER_DONOR].copy()
            else:
                gdf = gdf[gdf[elig_col]].copy()
            if ref_type == "mid_loo":
                gdf = gdf[gdf["age_years"] >= 40.0].copy()
            if gdf.empty:
                continue

            is_pool = bool(gdf["is_pooled"].iloc[0])
            fb_status = "not_evaluable_within_pool_df" if (not is_gao and not bool(gdf["yazar_pool_fe_all_evaluable"].iloc[0])) else "not_evaluable_within_pool_df_slice"

            # Precompute unadjusted depth-matched baselines & split-half counts for this (group, depth, ref_type)
            ref_arm_key = f"{'pool' if is_pool else 'ct'}_{'500_same' if depth == '500' else ('full_500' if depth == 'full_500_sens' else depth)}"
            ref_grp_key = "_POOLED_" if is_pool else grp
            h1_info = refs_summary.get(ref_arm_key, {}).get(ref_grp_key, {}).get("half1", {}) if (part == "A1" and ref_type == "young_loo") else {}

            sh_h2 = gdf.loc[gdf["a1_split_half_h2_to_h1"].notna(), "a1_split_half_h2_to_h1"] if (part == "A1" and ref_type == "young_loo") else pd.Series(dtype=float)
            sh_base_fallback = float(sh_h2.mean()) if not sh_h2.empty else np.nan
            self_fallback = float(gdf.loc[gdf["bin"] == "YOUNG", "a1_self_non_loo"].mean()) if (part == "A1" and ref_type == "young_loo") else np.nan
            ext = extra_baselines.get((depth, grp), {}) if (extra_baselines is not None and part == "A1" and ref_type == "young_loo") else {}
            sh_base = float(ext.get("split_half_young_baseline", sh_base_fallback)) if (part == "A1" and ref_type == "young_loo") else np.nan
            self_nonloo = float(ext.get("self_baseline_biased_low", self_fallback)) if (part == "A1" and ref_type == "young_loo") else np.nan

            y_loo_mean = float(gdf.loc[gdf["bin"] == "YOUNG", out_col].mean()) if (part == "A1" and ref_type == "young_loo") else np.nan
            m_mean = float(gdf.loc[gdf["bin"] == "MID", out_col].mean()) if (part == "A1") else np.nan
            o_mean = float(gdf.loc[gdf["bin"] == "OLD", out_col].mean()) if (part == "A1") else np.nan
            if h1_info:
                sh_h1_n = int(h1_info.get("h1_n_donors", 0))
                sh_h1_cells = int(h1_info.get("h1_n_cells", 0))
                sh_h2_n = int(len(h1_info.get("h2_donors", ())))
                sh_fb_flag = bool(h1_info.get("fallback_used", False))
            else:
                sh_h1_n = int(((gdf["bin"] == "YOUNG") & (gdf["split_half_group"] == 1)).sum()) if (part == "A1" and ref_type == "young_loo") else 0
                sh_h2_n = int(len(sh_h2))
                sh_h1_cells = int(gdf.loc[(gdf["bin"] == "YOUNG") & (gdf["split_half_group"] == 1), cell_col].sum()) if (part == "A1" and ref_type == "young_loo") else 0
                sh_fb_flag = False

            for spec_name, sub_s, nbc, bcol, min_bd, min_df, r_lv in _model_specs(gdf):
                row = fit_ols_model_row(
                    sub_s,
                    outcome_col=out_col,
                    predictor_col="age_per_decade",
                    nonbatch_covs=nbc,
                    batch_col=bcol,
                    min_batch_donors=min_bd,
                    min_residual_df=min_df,
                    ref_level=r_lv,
                    meta_fields={
                        "dataset": ds_name,
                        "part": part,
                        "metric": metric,
                        "depth": depth,
                        "reference_type": ref_type,
                        "model_spec": spec_name,
                        "group": grp,
                        "is_pooled": is_pool,
                        "is_stratified_global": False,
                        "eligibility_status": "eligible",
                        "default_fallback_status": fb_status,
                    },
                    fallback_to_nobatch=(not is_gao),
                )
                if row is None:
                    continue
                if part == "A1" and ref_type == "young_loo":
                    oy = fit_contrast_adj_diff(sub_s, out_col, "YOUNG", "OLD", nbc, bcol, min_bd, min_df, r_lv)
                    row.update(
                        {
                            "split_half_young_baseline": sh_base,
                            "split_half_n_donors_h1": sh_h1_n,
                            "split_half_n_cells_h1": sh_h1_cells,
                            "split_half_n_donors_h2": sh_h2_n,
                            "split_half_fallback_flag": sh_fb_flag,
                            "self_baseline_biased_low": self_nonloo,
                            "young_loo_mean": y_loo_mean,
                            "mid_mean": m_mean,
                            "old_mean": o_mean,
                            "old_minus_young_diff_adj": oy["diff"],
                            "old_minus_young_se_adj": oy["se"],
                            "old_minus_young_ci95_low_adj": oy["ci_low"],
                            "old_minus_young_ci95_high_adj": oy["ci_high"],
                            "old_minus_young_p_adj": oy["p"],
                            "old_minus_mid_diff_adj": np.nan,
                            "old_minus_mid_se_adj": np.nan,
                            "old_minus_mid_ci95_low_adj": np.nan,
                            "old_minus_mid_ci95_high_adj": np.nan,
                            "old_minus_mid_p_adj": np.nan,
                        }
                    )
                elif part == "A1" and ref_type == "mid_loo":
                    om = fit_contrast_adj_diff(sub_s, out_col, "MID", "OLD", nbc, bcol, min_bd, min_df, r_lv)
                    row.update(
                        {
                            "split_half_young_baseline": np.nan,
                            "split_half_n_donors_h1": 0,
                            "split_half_n_cells_h1": 0,
                            "split_half_n_donors_h2": 0,
                            "split_half_fallback_flag": False,
                            "self_baseline_biased_low": np.nan,
                            "young_loo_mean": np.nan,
                            "mid_mean": m_mean,
                            "old_mean": o_mean,
                            "old_minus_young_diff_adj": np.nan,
                            "old_minus_young_se_adj": np.nan,
                            "old_minus_young_ci95_low_adj": np.nan,
                            "old_minus_young_ci95_high_adj": np.nan,
                            "old_minus_young_p_adj": np.nan,
                            "old_minus_mid_diff_adj": om["diff"],
                            "old_minus_mid_se_adj": om["se"],
                            "old_minus_mid_ci95_low_adj": om["ci_low"],
                            "old_minus_mid_ci95_high_adj": om["ci_high"],
                            "old_minus_mid_p_adj": om["p"],
                        }
                    )
                else:
                    row.update(
                        {
                            "split_half_young_baseline": np.nan,
                            "split_half_n_donors_h1": 0,
                            "split_half_n_cells_h1": 0,
                            "split_half_n_donors_h2": 0,
                            "split_half_fallback_flag": False,
                            "self_baseline_biased_low": np.nan,
                            "young_loo_mean": np.nan,
                            "mid_mean": np.nan,
                            "old_mean": np.nan,
                            "old_minus_young_diff_adj": np.nan,
                            "old_minus_young_se_adj": np.nan,
                            "old_minus_young_ci95_low_adj": np.nan,
                            "old_minus_young_ci95_high_adj": np.nan,
                            "old_minus_young_p_adj": np.nan,
                            "old_minus_mid_diff_adj": np.nan,
                            "old_minus_mid_se_adj": np.nan,
                            "old_minus_mid_ci95_low_adj": np.nan,
                            "old_minus_mid_ci95_high_adj": np.nan,
                            "old_minus_mid_p_adj": np.nan,
                        }
                    )
                part_a_models.append(row)

        # Stratified donor-clustered global sensitivity across non-pooled eligible cell types
        ct_rows = [
            r
            for r in part_a_models
            if r["dataset"] == ds_name and r["part"] == part and r["depth"] == depth and r["reference_type"] == ref_type and r["model_spec"] == "primary" and not r["is_pooled"]
        ]
        elig_cts = [r["group"] for r in ct_rows]
        elig_col = "eligible_a1_3000" if depth == "3000" else ("eligible_a3_1500" if part == "A3" else "eligible_a1_1500")
        strat_src = donors_long[donors_long["n_cells_ge500"] >= MIN_CELLS_PER_DONOR] if depth == "full_500_sens" else donors_long[donors_long[elig_col]]
        if ref_type == "mid_loo":
            strat_src = strat_src[strat_src["age_years"] >= 40.0]
        strat_row = _fit_stratified_global_row(
            strat_src,
            elig_cts=elig_cts,
            out_col=out_col,
            cell_col=cell_col,
            depth_col=depth_col,
            ds_name=ds_name,
            part=part,
            metric=metric,
            depth=depth,
            ref_type=ref_type,
            include_part_a_cols=True,
        )
        if strat_row is not None:
            part_a_models.append(strat_row)


def run_a2_pipeline(meta: dict[str, Any], donors_long: pd.DataFrame, a2_inputs: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run Item-8 variance/eta2 comparisons (OLD_vs_YOUNG & MID_vs_YOUNG) across both residualization arms."""
    ds_name = meta["dataset"]
    is_gao = ds_name == "gao"
    donors_df = meta["donors_df"]
    donor_ids = meta["donor_ids"]
    ct_names = meta["cell_types"]
    pooled_name = meta["pooled_name"]
    n_cts = len(ct_names)
    d_to_i = {d: i for i, d in enumerate(donor_ids)}

    rows: list[dict[str, Any]] = []
    inelig: list[dict[str, Any]] = []

    if not is_gao:
        manifest = json.loads((YAZAR_CACHE / "pb_manifest.json").read_text(encoding="utf-8"))
        # Gene order alignment check
        sample_pkl = next(iter(manifest.values()))
        with open(YAZAR_CACHE / sample_pkl, "rb") as pf:
            pb_s = pickle.load(pf)
        proj_info = build_human_gene_projection(pb_s.index.astype(str).tolist())
        cache_genes = proj_info["audit"]["retained_symbols"]
        assert cache_genes == a2_inputs["retained_symbols"], "Yazar cache genes differ from H5AD retained genes!"
        g_proj_cache = proj_info["proj"]
    else:
        manifest = {}
        g_proj_cache = None

    def _load_group_pb_counts(grp: str, sub_donors: list[str]) -> np.ndarray:
        if is_gao:
            if grp == pooled_name:
                idxs = [d_to_i[d] for d in sub_donors]
                return a2_inputs["a2_pool"][idxs]
            ci = ct_names.index(grp)
            idxs = [d_to_i[d] * n_cts + ci for d in sub_donors]
            return a2_inputs["a2_ct_gao"][idxs]
        else:
            if grp == pooled_name:
                # Doublet-free All_Cells accumulated in Pass 1 (all 1,248,849 singlets!)
                idxs = [d_to_i[d] for d in sub_donors]
                return a2_inputs["a2_pool"][idxs]
            with open(YAZAR_CACHE / manifest[grp], "rb") as pf:
                raw_df = pickle.load(pf)
            raw_df.columns = raw_df.columns.astype(str)
            sub_raw = raw_df[sub_donors].to_numpy(dtype=np.float64).T  # n_donors x n_raw_genes
            return (sp.csr_matrix(sub_raw) @ g_proj_cache).toarray()

    all_groups = [*ct_names, pooled_name]
    for contrast_idx, (contrast_name, comp_bin) in enumerate((("OLD_vs_YOUNG", "OLD"), ("MID_vs_YOUNG", "MID"))):
        elig_col = "eligible_a2_oy" if contrast_name == "OLD_vs_YOUNG" else "eligible_a2_my"
        u_donors = sorted(donors_long.loc[donors_long["bin"].isin(["YOUNG", comp_bin]), "donor_id"].astype(str).unique().tolist())
        u_to_pos = {d: i for i, d in enumerate(u_donors)}
        rng_sync = np.random.default_rng(MASTER_SEED + 500 + contrast_idx)
        global_ranks = np.array([rng_sync.permutation(len(u_donors)) for _ in range(A2_N_PERM)], dtype=np.int32)
        strat_accum: dict[tuple[str, float], list[dict[str, Any]]] = {}

        for grp in all_groups:
            gdf = donors_long[(donors_long["group"] == grp) & (donors_long["bin"].isin(["YOUNG", comp_bin]))].copy()
            gdf = gdf[gdf["n_cells_all"] >= MIN_CELLS_PER_DONOR].copy()
            if (not is_gao) and (grp != pooled_name) and (grp not in manifest):
                inelig.append({"dataset": ds_name, "contrast": contrast_name, "group": grp, "reason": "not_in_yazar_cache"})
                continue
            n_y = int((gdf["bin"] == "YOUNG").sum())
            n_c = int((gdf["bin"] == comp_bin).sum())
            if n_y < MIN_DONORS_A2_BIN or n_c < MIN_DONORS_A2_BIN:
                inelig.append({"dataset": ds_name, "contrast": contrast_name, "group": grp, "n_young": n_y, "n_comp": n_c, "reason": "lt6_in_young_or_comp"})
                continue

            d_order = gdf["donor_id"].astype(str).tolist()
            pb_raw = _load_group_pb_counts(grp, d_order)
            det_mask = np.sum(pb_raw > 0, axis=0) >= (0.5 * len(d_order))
            pb_det = pb_raw[:, det_mask]
            tot = pb_det.sum(axis=1, keepdims=True)
            expr_log = np.log1p(pb_det / np.maximum(tot, 1.0) * 1e4)
            is_comp_arr = (gdf["bin"] == comp_bin).to_numpy(dtype=bool)
            mlu_arr = gdf["mean_log_umi_all"].to_numpy(dtype=float)

            # Arm 1: Depth-only (Item 8 direct port via imported _residualize_against_depth)
            expr_depth = _residualize_against_depth(expr_log, mlu_arr)

            # Arm 2: Batch-residualized with dynamic >=2 donor batch check, leverage guard, and fallback
            gdf["_is_comp"] = is_comp_arr.astype(float)
            sub_b, cols_b, _, st_b = _prune_covariate_matrix(
                gdf,
                nonbatch_covs=["sex_male", "mean_log_umi_all"],
                batch_col="study_or_pool",
                min_batch_donors=2,
                min_residual_df=5 if is_gao else 10,
                ref_level="ENCODE v4 (Snyder)" if is_gao else None,
                require_is_comp_collinear_check="_is_comp",
                atomic_batch_block=True,
            )
            pool_status = "evaluated"
            n_excl_sing = 0
            if st_b == "evaluated":
                keep_pos = [d_order.index(d) for d in sub_b["donor_id"].astype(str).tolist()]
                is_comp_b = is_comp_arr[keep_pos]
                if int((~is_comp_b).sum()) < MIN_DONORS_A2_BIN or int(is_comp_b.sum()) < MIN_DONORS_A2_BIN:
                    st_b = "insufficient_after_singleton_drop"
                else:
                    n_excl_sing = len(d_order) - len(keep_pos)
                    pb_sub = pb_raw[keep_pos]
                    det_b = np.sum(pb_sub > 0, axis=0) >= (0.5 * len(keep_pos))
                    pb_bd = pb_sub[:, det_b]
                    tot_b = pb_bd.sum(axis=1, keepdims=True)
                    expr_b = np.log1p(pb_bd / np.maximum(tot_b, 1.0) * 1e4)
                    z_b = np.column_stack([np.ones(len(keep_pos), dtype=float), sub_b[cols_b].to_numpy(dtype=float)])
                    q_chk, _ = np.linalg.qr(np.column_stack([z_b, is_comp_b.astype(float)]), mode="reduced")
                    if np.any((1.0 - np.sum(q_chk ** 2, axis=1)) <= 1e-4):
                        st_b = "extreme_leverage_guard"
            if st_b != "evaluated":
                pool_status = f"not_evaluable_within_{'study' if is_gao else 'pool'}_df_slice"
                sub_b2, cols_b2, _, _ = _prune_covariate_matrix(
                    gdf,
                    nonbatch_covs=["sex_male", "mean_log_umi_all"],
                    batch_col=None,
                    min_batch_donors=1,
                    min_residual_df=5,
                    ref_level=None,
                    require_is_comp_collinear_check="_is_comp",
                )
                keep_pos = [d_order.index(d) for d in sub_b2["donor_id"].astype(str).tolist()]
                is_comp_b = is_comp_arr
                expr_b = expr_log
                z_b = np.column_stack([np.ones(len(d_order), dtype=float), sub_b2[cols_b2].to_numpy(dtype=float)])

            d_order_b = [d_order[k] for k in keep_pos]
            orders_d = np.argsort(global_ranks[:, np.array([u_to_pos[d] for d in d_order], dtype=np.int32)], axis=1)
            orders_b = orders_d if len(d_order_b) == len(d_order) and d_order_b == d_order else np.argsort(global_ranks[:, np.array([u_to_pos[d] for d in d_order_b], dtype=np.int32)], axis=1)
            is_pool = grp == pooled_name
            n_cells_tot = int(gdf["n_cells_all"].sum())
            n_cells_b = int(gdf.iloc[keep_pos]["n_cells_all"].sum())

            for eps_val in (1e-8, 1e-3):
                st_batch, perm_b = compute_a2_variance_metrics(
                    expr_b, is_comp_b, cov_mat_with_intercept=z_b, n_perm=A2_N_PERM,
                    seed=MASTER_SEED + (0 if eps_val < 1e-5 else 1), eps=eps_val,
                    custom_perm_orders=None if is_pool else orders_b, return_perm=True,
                )
                st_depth, perm_d = compute_a2_variance_metrics(
                    expr_depth, is_comp_arr, cov_mat_with_intercept=None, n_perm=A2_N_PERM,
                    seed=MASTER_SEED + (2 if eps_val < 1e-5 else 3), eps=eps_val,
                    custom_perm_orders=None if is_pool else orders_d, return_perm=True,
                )
                disagree = (st_batch["mean_log_var_ratio"] * st_depth["mean_log_var_ratio"]) < 0
                row_b = {
                    "dataset": ds_name,
                    "contrast": contrast_name,
                    "arm": "batch_residualized",
                    "eps_used": eps_val,
                    "group": grp,
                    "is_pooled": is_pool,
                    "is_stratified_global": False,
                    "pool_fe_status": pool_status,
                    "n_donors_ref_young": int((~is_comp_b).sum()),
                    "n_donors_comp": int(is_comp_b.sum()),
                    "n_donors_excluded_singleton_batch": int(n_excl_sing),
                    "n_cells_total": n_cells_b,
                    **st_batch,
                    "gene_count": st_batch.pop("n_genes"),
                    "q_bh_eta2": np.nan,
                    "q_bh_lvr_one_sided": np.nan,
                    "q_bh_lvr_two_sided": np.nan,
                    "q_bh_tot_lvr_one_sided": np.nan,
                    "family_size_m": np.nan,
                    "arm_direction_disagreement_flag": bool(disagree),
                }
                row_d = {
                    "dataset": ds_name,
                    "contrast": contrast_name,
                    "arm": "depth_only_residualized",
                    "eps_used": eps_val,
                    "group": grp,
                    "is_pooled": is_pool,
                    "is_stratified_global": False,
                    "pool_fe_status": "depth_only_all_donors",
                    "n_donors_ref_young": n_y,
                    "n_donors_comp": n_c,
                    "n_donors_excluded_singleton_batch": 0,
                    "n_cells_total": n_cells_tot,
                    **st_depth,
                    "gene_count": st_depth.pop("n_genes"),
                    "q_bh_eta2": np.nan,
                    "q_bh_lvr_one_sided": np.nan,
                    "q_bh_lvr_two_sided": np.nan,
                    "q_bh_tot_lvr_one_sided": np.nan,
                    "family_size_m": np.nan,
                    "arm_direction_disagreement_flag": bool(disagree),
                }
                rows.extend([row_b, row_d])
                if not is_pool:
                    if pool_status == "evaluated":
                        strat_accum.setdefault(("batch_residualized", eps_val), []).append(
                            {"row": row_b, "perm": perm_b, "y_set": {d_order_b[i] for i, v in enumerate(is_comp_b) if not v}, "c_set": {d_order_b[i] for i, v in enumerate(is_comp_b) if v}}
                        )
                    strat_accum.setdefault(("depth_only_residualized", eps_val), []).append(
                        {"row": row_d, "perm": perm_d, "y_set": {d_order[i] for i, v in enumerate(is_comp_arr) if not v}, "c_set": {d_order[i] for i, v in enumerate(is_comp_arr) if v}}
                    )

        # Stratified across-cell-type synchronized permutation global rows per contrast and eps
        for eps_val in (1e-8, 1e-3):
            b_items = strat_accum.get(("batch_residualized", eps_val), [])
            d_items = strat_accum.get(("depth_only_residualized", eps_val), [])
            b_stat = None
            d_stat = None
            for arm_lbl, items in (("batch_residualized", b_items), ("depth_only_residualized", d_items)):
                if len(items) < 2:
                    continue
                st_joint = _summarize_a2_perm_stats(
                    n_g=int(round(np.mean([it["row"]["gene_count"] for it in items]))),
                    resid_df_model=int(round(np.mean([it["row"]["residual_df_resid_model"] for it in items]))),
                    obs_mean_eta2=float(np.mean([it["row"]["mean_eta2"] for it in items])),
                    obs_med_eta2=float(np.median([it["row"]["mean_eta2"] for it in items])),
                    obs_mean_f=float(np.mean([it["row"]["mean_F"] for it in items])),
                    obs_mean_lvr=float(np.mean([it["row"]["mean_log_var_ratio"] for it in items])),
                    obs_med_lvr=float(np.median([it["row"]["mean_log_var_ratio"] for it in items])),
                    obs_tot_lvr=float(np.mean([it["row"]["log_total_var_ratio"] for it in items])),
                    obs_frac_gt=float(np.mean([it["row"]["frac_genes_comp_gt_young"] for it in items])),
                    perm_eta2=np.mean([it["perm"]["perm_eta2"] for it in items], axis=0),
                    perm_lvr=np.mean([it["perm"]["perm_lvr"] for it in items], axis=0),
                    perm_tot_lvr=np.mean([it["perm"]["perm_tot_lvr"] for it in items], axis=0),
                )
                n_g_val = st_joint.pop("n_genes")
                u_y = set().union(*(it["y_set"] for it in items))
                u_c = set().union(*(it["c_set"] for it in items))
                r_obj = {
                    "dataset": ds_name,
                    "contrast": contrast_name,
                    "arm": arm_lbl,
                    "eps_used": eps_val,
                    "group": "Stratified_CellType_Global",
                    "is_pooled": True,
                    "is_stratified_global": True,
                    "pool_fe_status": f"stratified_{len(items)}_cell_types",
                    "n_donors_ref_young": int(len(u_y)),
                    "n_donors_comp": int(len(u_c)),
                    "n_donors_excluded_singleton_batch": 0,
                    "n_cells_total": int(sum(it["row"]["n_cells_total"] for it in items)),
                    **st_joint,
                    "gene_count": n_g_val,
                    "q_bh_eta2": np.nan,
                    "q_bh_lvr_one_sided": np.nan,
                    "q_bh_lvr_two_sided": np.nan,
                    "q_bh_tot_lvr_one_sided": np.nan,
                    "family_size_m": np.nan,
                    "arm_direction_disagreement_flag": False,
                }
                if arm_lbl == "batch_residualized":
                    b_stat = r_obj
                else:
                    d_stat = r_obj
            if b_stat is not None and d_stat is not None:
                dis_s = (b_stat["mean_log_var_ratio"] * d_stat["mean_log_var_ratio"]) < 0
                b_stat["arm_direction_disagreement_flag"] = bool(dis_s)
                d_stat["arm_direction_disagreement_flag"] = bool(dis_s)
            if b_stat is not None:
                rows.append(b_stat)
            if d_stat is not None:
                rows.append(d_stat)

    var_df = pd.DataFrame(rows)
    if not var_df.empty:
        non_strat = (~var_df["is_pooled"]) & (~var_df["is_stratified_global"].astype(bool)) & (~var_df["pool_fe_status"].str.startswith("not_evaluable"))
        for (contrast, arm, eps_v), sub_idx in var_df[non_strat].groupby(["contrast", "arm", "eps_used"], sort=False).groups.items():
            idx_l = list(sub_idx)
            var_df.loc[idx_l, "family_size_m"] = len(idx_l)
            for p_col, q_col in (
                ("p_eta2_perm", "q_bh_eta2"),
                ("p_lvr_one_sided_comp_gt_young", "q_bh_lvr_one_sided"),
                ("p_lvr_two_sided", "q_bh_lvr_two_sided"),
                ("p_tot_lvr_one_sided", "q_bh_tot_lvr_one_sided"),
            ):
                var_df.loc[idx_l, q_col] = multipletests(var_df.loc[idx_l, p_col].to_numpy(float), method="fdr_bh")[1]
        p_idx = var_df["is_pooled"].astype(bool) & (~var_df["is_stratified_global"].astype(bool))
        var_df.loc[p_idx, "family_size_m"] = 1
        for p_col, q_col in (
            ("p_eta2_perm", "q_bh_eta2"),
            ("p_lvr_one_sided_comp_gt_young", "q_bh_lvr_one_sided"),
            ("p_lvr_two_sided", "q_bh_lvr_two_sided"),
            ("p_tot_lvr_one_sided", "q_bh_tot_lvr_one_sided"),
        ):
            var_df.loc[p_idx, q_col] = var_df.loc[p_idx, p_col]

    return var_df, {"ineligible": inelig}


def run_partb_pipeline(
    meta: dict[str, Any],
    donors_long: pd.DataFrame,
    partb_ckpt: dict[str, Any],
    primary_depth: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit all Part B Spread models (500, 1500, 3000, full_depth) + 20-seed KL / 5-seed Shannon robustness."""
    ds_name = meta["dataset"]
    is_gao = ds_name == "gao"
    donors_df = meta["donors_df"]
    ct_names = meta["cell_types"]
    pooled_name = meta["pooled_name"]
    n_donors = len(meta["donor_ids"])
    n_cts = len(ct_names)
    b_sums = partb_ckpt["b_sums"]
    cnt_3d = partb_ckpt["cnt_3d"]
    pool_mat = partb_ckpt["pool_mat"]
    mlu_ct = partb_ckpt["mlu_ct"]
    mlu_pool = partb_ckpt["mlu_pool"]

    # Map depth strings -> (cell count idx, mean_log_umi idx, numerator matrix key)
    depth_configs = [
        ("500", 1, 1, "ct_500", "pool_500"),
        ("1500", 2, 2, "ct_1500_seeds", "pool_1500_seeds"),
        ("3000", 3, 3, "ct_3000", "pool_3000"),
        (f"full_{primary_depth}", 1 if primary_depth == 500 else 2, 1 if primary_depth == 500 else 2, "ct_full", "pool_full"),
        ("full_gt0", 0, 0, "ct_full_gt0", "pool_full_gt0"),
    ]

    model_rows: list[dict[str, Any]] = []
    robust_rows: list[dict[str, Any]] = []

    groups = [*ct_names, pooled_name]
    for depth_lbl, cnt_k, mlu_k, ct_key, pl_key in depth_configs:
        strat_parts_b: dict[str, list[pd.DataFrame]] = {m_name: [] for m_name in ("shannon", "deman_kl_methods", "deman_kl_sqrt")}
        strat_elig_b: dict[str, list[str]] = {m_name: [] for m_name in ("shannon", "deman_kl_methods", "deman_kl_sqrt")}
        for gi, grp in enumerate(groups):
            is_pool = grp == pooled_name
            c_vec = pool_mat[:, cnt_k] if is_pool else cnt_3d[:, gi, cnt_k]
            elig_idx = np.flatnonzero(c_vec >= MIN_CELLS_PER_DONOR)
            if len(elig_idx) < MIN_DONORS_SLOPE:
                continue

            sub_df = donors_df.iloc[elig_idx].copy()
            sub_df["n_cells"] = c_vec[elig_idx]
            sub_df["_depth"] = (mlu_pool[elig_idx, mlu_k] if is_pool else mlu_ct[elig_idx, gi, mlu_k])
            sub_df["group"] = grp
            if is_gao:
                sub_df["study_or_pool"] = sub_df["study"]
            else:
                sub_df["study_or_pool"] = sub_df["pool_number"].astype(str)

            fb_status = "not_evaluable_within_pool_df" if (not is_gao and grp in meta.get("not_evaluable_all_umi_csv", [])) else "not_evaluable_within_pool_df_slice"

            specs = (
                [
                    ("primary", sub_df, ["sex_male", "_depth"], "study_or_pool", 2, 5, "ENCODE v4 (Snyder)"),
                    ("batch_adj_no_depth", sub_df, ["sex_male"], "study_or_pool", 2, 5, "ENCODE v4 (Snyder)"),
                    ("age_only", sub_df, [], None, 1, 5, None),
                    ("age_plus_sex", sub_df, ["sex_male"], None, 1, 5, None),
                    ("10x_only_no_penn", sub_df[sub_df["technology"] != "Dropseq"].copy(), ["sex_male", "_depth"], "study_or_pool", 2, 5, "ENCODE v4 (Snyder)"),
                    ("studies_ge8", sub_df[~sub_df["study"].isin(["Kuppe 2022", "Sim 2021"])].copy(), ["sex_male", "_depth"], "study_or_pool", 2, 5, "ENCODE v4 (Snyder)"),
                    ("encode_v4_only", sub_df[sub_df["study"] == "ENCODE v4 (Snyder)"].copy(), ["sex_male", "_depth"], None, 1, 5, None),
                ]
                if is_gao
                else [
                    ("primary", sub_df, ["sex_male", "_depth"], "study_or_pool", 2, 10, None),
                    ("pool_fe_no_depth", sub_df, ["sex_male"], "study_or_pool", 2, 10, None),
                    ("no_pool", sub_df, ["sex_male", "_depth"], None, 1, 5, None),
                    ("age_ge_20", sub_df[sub_df["age_years"] >= 20.0].copy(), ["sex_male", "_depth"], "study_or_pool", 2, 10, None),
                ]
            )

            # Base numerator array (seed 0 for 1500)
            if depth_lbl == "1500":
                arr_ct = b_sums[ct_key][0].reshape(n_donors, n_cts, 3)[:, gi, :] if not is_pool else b_sums[pl_key][0]
            else:
                arr_ct = b_sums[ct_key].reshape(n_donors, n_cts, 3)[:, gi, :] if not is_pool else b_sums[pl_key]
            for m_i, m_name in enumerate(("shannon", "deman_kl_methods", "deman_kl_sqrt")):
                sub_df[f"_{m_name}"] = arr_ct[elig_idx, m_i] / np.maximum(c_vec[elig_idx].astype(np.float64), 1.0)
                for spec_name, sp_df, nbc, bcol, min_bd, min_df, r_lv in specs:
                    sp_copy = sp_df.copy()
                    sp_copy[f"_{m_name}"] = sub_df.loc[sp_copy.index, f"_{m_name}"]
                    row = fit_ols_model_row(
                        sp_copy,
                        outcome_col=f"_{m_name}",
                        predictor_col="age_per_decade",
                        nonbatch_covs=nbc,
                        batch_col=bcol,
                        min_batch_donors=min_bd,
                        min_residual_df=min_df,
                        ref_level=r_lv,
                        meta_fields={
                            "dataset": ds_name,
                            "part": "B",
                            "metric": m_name,
                            "depth": depth_lbl,
                            "reference_type": "none",
                            "model_spec": spec_name,
                            "group": grp,
                            "is_pooled": is_pool,
                            "is_stratified_global": False,
                            "eligibility_status": "eligible",
                            "default_fallback_status": fb_status,
                        },
                        fallback_to_nobatch=(not is_gao),
                    )
                    if row is not None:
                        model_rows.append(row)
                        if spec_name == "primary" and not is_pool:
                            strat_elig_b[m_name].append(grp)
                            strat_parts_b[m_name].append(sub_df[["donor_id", "group", "age_per_decade", "sex_male", "n_cells", "_depth", f"_{m_name}"]].copy())

        for m_name in ("shannon", "deman_kl_methods", "deman_kl_sqrt"):
            if len(strat_elig_b[m_name]) >= 2 and strat_parts_b[m_name]:
                st_df_b = pd.concat(strat_parts_b[m_name], ignore_index=True)
                s_row = _fit_stratified_global_row(
                    st_df_b,
                    elig_cts=strat_elig_b[m_name],
                    out_col=f"_{m_name}",
                    cell_col="n_cells",
                    depth_col="_depth",
                    ds_name=ds_name,
                    part="B",
                    metric=m_name,
                    depth=depth_lbl,
                    ref_type="none",
                    include_part_a_cols=False,
                )
                if s_row is not None:
                    model_rows.append(s_row)

    # Seed robustness sweeps (at 1500 and at 500 if Gao fallback triggered)
    seed_sweeps = [(1500, 2, "ct_1500_seeds", "pool_1500_seeds")]
    if is_gao and "ct_500_seeds" in b_sums:
        seed_sweeps.append((500, 1, "ct_500_seeds", "pool_500_seeds"))

    for d_int, cnt_k, ct_sk, pl_sk in seed_sweeps:
        for gi, grp in enumerate(groups):
            is_pool = grp == pooled_name
            c_vec = pool_mat[:, cnt_k] if is_pool else cnt_3d[:, gi, cnt_k]
            elig_idx = np.flatnonzero(c_vec >= MIN_CELLS_PER_DONOR)
            if len(elig_idx) < MIN_DONORS_SLOPE:
                continue
            sub_df = donors_df.iloc[elig_idx].copy()
            sub_df["n_cells"] = c_vec[elig_idx]
            sub_df["_depth"] = (mlu_pool[elig_idx, cnt_k] if is_pool else mlu_ct[elig_idx, gi, cnt_k])
            sub_df["study_or_pool"] = sub_df["study"] if is_gao else sub_df["pool_number"].astype(str)

            for m_i, m_name, n_seeds_req in ((0, "shannon", 5), (1, "deman_kl_methods", 20), (2, "deman_kl_sqrt", 20)):
                slopes = []
                slopes_sd = []
                pvals = []
                for si in range(n_seeds_req):
                    arr_s = b_sums[pl_sk][si] if is_pool else b_sums[ct_sk][si].reshape(n_donors, n_cts, 3)[:, gi, :]
                    sub_df["_y"] = arr_s[elig_idx, m_i] / np.maximum(c_vec[elig_idx].astype(np.float64), 1.0)
                    r = fit_ols_model_row(
                        sub_df,
                        outcome_col="_y",
                        predictor_col="age_per_decade",
                        nonbatch_covs=["sex_male", "_depth"],
                        batch_col="study_or_pool",
                        min_batch_donors=2,
                        min_residual_df=5 if is_gao else 10,
                        ref_level="ENCODE v4 (Snyder)" if is_gao else None,
                        meta_fields={},
                        fallback_to_nobatch=(not is_gao),
                    )
                    if r is None:
                        break
                    slopes.append(float(r["slope_raw_per_decade"]))
                    slopes_sd.append(float(r["slope_sd_per_decade"]))
                    pvals.append(float(r["pvalue"]))
                if len(slopes) != n_seeds_req:
                    continue
                pri_s = slopes[0]
                pri_sign = 1 if pri_s >= 0 else -1
                n_same = int(sum((1 if s >= 0 else -1) == pri_sign for s in slopes))
                n_p05 = int(sum(p < 0.05 for p in pvals))
                n_both = int(sum(((1 if s >= 0 else -1) == pri_sign) and (p < 0.05) for s, p in zip(slopes, pvals)))
                min_req = 4 if n_seeds_req == 5 else 18
                robust_rows.append(
                    {
                        "dataset": ds_name,
                        "depth": d_int,
                        "metric": m_name,
                        "model_spec": "primary",
                        "group": grp,
                        "is_pooled": is_pool,
                        "primary_seed": int(ENTROPY_SEEDS[0]),
                        "primary_slope_raw": pri_s,
                        "primary_slope_sd": float(slopes_sd[0]),
                        "primary_pvalue": float(pvals[0]),
                        "n_seeds": n_seeds_req,
                        "n_same_sign": n_same,
                        "frac_same_sign": float(n_same / n_seeds_req),
                        "n_p_lt_05": n_p05,
                        "n_same_sign_and_p_lt_05": n_both,
                        "min_slope_raw": float(min(slopes)),
                        "median_slope_raw": float(np.median(slopes)),
                        "max_slope_raw": float(max(slopes)),
                        "min_pvalue": float(min(pvals)),
                        "median_pvalue": float(np.median(pvals)),
                        "max_pvalue": float(max(pvals)),
                        "passes_robustness": bool(n_same >= min_req),
                    }
                )

    mod_df = pd.DataFrame(model_rows)
    if not mod_df.empty:
        mod_df = apply_slice_bh(mod_df, ["dataset", "part", "metric", "depth", "model_spec"])
    rob_df = pd.DataFrame(robust_rows)
    return mod_df, rob_df


def _log_stage(ds_name: str, msg: str) -> None:
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(f"[{ds_name}] {msg} | self_max_rss={rss_mb:.1f} MB", flush=True)


def run_dataset_extract_and_part_a(ds_name: str, skip_part_b_stream: bool = False) -> dict[str, Any]:
    """Execute Pass 0, Pass 1, Pass 2, Part A1/A2/A3 models & export Part A + cache for Part B."""
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    meta = reconcile_gao_metadata(GAO_H5AD) if ds_name == "gao" else reconcile_yazar_metadata(YAZAR_H5AD, YAZAR_CACHE)
    proj_info = build_human_gene_projection(meta["gene_symbols"])
    g_proj = proj_info["proj"]
    _log_stage(ds_name, f"reconciled metadata & gene proj ({proj_info['audit']['n_genes_retained']} retained)")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, {})) as ex:
        p0 = stream_dataset_pass0(meta, g_proj, ex)
    _log_stage(ds_name, f"Pass 0 complete in {p0['elapsed_s']}s")

    depth_df = p0["depth_table"]
    if ds_name == "gao":
        ct_depth = depth_df[~depth_df["is_pooled"]]
        few_1500 = ct_depth[ct_depth["donors_ge20_at_1500"] < MIN_DONORS_SLOPE]["group"].tolist()
        partb_primary_depth = 500 if len(few_1500) > 0 else 1500
        do_500_ms = partb_primary_depth == 500
        run_gao_500_a1 = int((ct_depth["donors_ge20_at_1500"] < MIN_DONORS_SLOPE).sum()) > 2
    else:
        few_1500 = []
        partb_primary_depth = 1500
        do_500_ms = False
        run_gao_500_a1 = False

    with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, {})) as ex:
        p1 = stream_dataset_pass1(
            meta,
            p0,
            g_proj,
            executor=ex,
            do_500_multiseed=do_500_ms,
            partb_full_floor=partb_primary_depth,
            run_gao_a3_centroids=(ds_name == "gao"),
        )
    _log_stage(ds_name, f"Pass 1 complete in {p1['elapsed_s']}s")

    p2 = stream_dataset_pass2(meta, p0, p1, g_proj, run_gao_500_a1=run_gao_500_a1)
    _log_stage(ds_name, f"Pass 2 complete in {p2['elapsed_s']}s")
    donors_long, partb_ckpt = _build_donors_table_and_partb_cache(meta, p0, p1, p2)

    part_a_rows: list[dict[str, Any]] = []
    _run_a1_a3_models_from_donors(ds_name, donors_long, p1["a1_refs"], part_a_rows, extra_baselines=partb_ckpt.get("a1_depth_baselines"))
    part_a_df = apply_slice_bh(pd.DataFrame(part_a_rows), ["dataset", "part", "metric", "depth", "reference_type", "model_spec"])
    _log_stage(ds_name, f"Part A1/A3 models complete ({len(part_a_df)} rows)")

    # Free A1/A3 intermediate refs before running A2
    a2_inputs = {
        "retained_symbols": proj_info["audit"]["retained_symbols"],
        "a2_ct_gao": p1.pop("a2_ct_gao"),
        "a2_pool": p1.pop("a2_pool"),
    }
    a1_elig_log = p2["a1_elig_log"]
    pass_timings = {"pass0": p0["elapsed_s"], "pass1": p1["elapsed_s"], "pass2": p2["elapsed_s"]}
    p1_gate = p0["gate"]
    p1.clear()
    p2.clear()
    gc.collect()
    var_df, var_audit = run_a2_pipeline(meta, donors_long, a2_inputs)
    del a2_inputs
    gc.collect()
    _log_stage(ds_name, f"Part A2 complete ({len(var_df)} rows)")

    # Export Commit 1 / Commit 2 CSVs immediately
    write_clean_csv(donors_long, RESULTS_DIR / f"followup-between-subject-{ds_name}-donors.csv")
    write_clean_csv(part_a_df, RESULTS_DIR / f"followup-between-subject-{ds_name}-models.csv")
    write_clean_csv(var_df, RESULTS_DIR / f"followup-between-subject-variance-{ds_name}.csv")

    # Save compact (<50 MB) Part B & QC checkpoint outside git
    ckpt = {
        "meta_summary": {
            "dataset": ds_name,
            "reconcile": meta["reconcile"],
            "gene_audit": {k: v for k, v in proj_info["audit"].items() if k != "retained_symbols"},
            "count_gate_pass0": p1_gate,
            "depth_table": depth_df.to_dict(orient="records"),
            "gao_partb_fallback_triggered": partb_primary_depth == 500,
            "gao_partb_primary_depth": partb_primary_depth,
            "gao_partb_1500_dropped_groups": few_1500,
            "a1_elig_log": a1_elig_log,
            "a2_ineligible": var_audit["ineligible"],
            "timings_s": pass_timings,
        },
        "meta_light": {
            "dataset": ds_name,
            "donors_df": meta["donors_df"],
            "donor_ids": meta["donor_ids"],
            "cell_types": meta["cell_types"],
            "pooled_name": meta["pooled_name"],
            "not_evaluable_all_umi_csv": meta.get("not_evaluable_all_umi_csv", []),
        },
        "donors_long": donors_long,
        "partb_ckpt": partb_ckpt,
        "partb_primary_depth": partb_primary_depth,
    }
    with open(CHECKPOINT_DIR / f"{ds_name}_ckpt.pkl", "wb") as pf:
        pickle.dump(ckpt, pf, protocol=pickle.HIGHEST_PROTOCOL)
    return ckpt["meta_summary"]


def run_part_b_and_qc() -> dict[str, Any]:
    """Load Part A/B external checkpoints, fit Part B spread & robustness, compute verdicts & QC JSON."""
    qc_out: dict[str, Any] = {"datasets": {}, "verdicts": {}, "seeds": {"master_seed": MASTER_SEED, "entropy_seeds": ENTROPY_SEEDS, "a2_n_perm": A2_N_PERM}}
    rob_parts = []
    for ds_name in ("gao", "yazar"):
        ckpt_path = CHECKPOINT_DIR / f"{ds_name}_ckpt.pkl"
        with open(ckpt_path, "rb") as pf:
            ckpt = pickle.load(pf)
        mod_b, rob_b = run_partb_pipeline(ckpt["meta_light"], ckpt["donors_long"], ckpt["partb_ckpt"], ckpt["partb_primary_depth"])
        write_clean_csv(mod_b, RESULTS_DIR / f"followup-spread-{ds_name}-models.csv")
        rob_parts.append(rob_b)

        part_a_df = pd.read_csv(RESULTS_DIR / f"followup-between-subject-{ds_name}-models.csv")
        var_df = pd.read_csv(RESULTS_DIR / f"followup-between-subject-variance-{ds_name}.csv")
        verdicts = evaluate_verdicts(part_a_df, var_df, mod_b, rob_b, ckpt["partb_primary_depth"], ckpt["meta_light"]["pooled_name"])
        qc_out["datasets"][ds_name] = ckpt["meta_summary"]
        qc_out["verdicts"][ds_name] = verdicts

    rob_all = pd.concat(rob_parts, ignore_index=True) if rob_parts else pd.DataFrame()
    write_clean_csv(rob_all, RESULTS_DIR / "followup-between-subject-robustness.csv")

    qc_text = json.dumps(qc_out, indent=2) + "\n"
    assert "/home/reza" not in qc_text, "Forbidden home path detected in QC JSON"
    (RESULTS_DIR / "followup-between-subject-gao-yazar-qc.json").write_text(qc_text, encoding="utf-8")
    return qc_out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["gao", "yazar", "all"], default="all")
    parser.add_argument("--stage", choices=["part-a", "part-b", "all"], default="all")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if args.stage in ("part-a", "all"):
        targets = ["gao", "yazar"] if args.dataset == "all" else [args.dataset]
        for ds in targets:
            print(f"=== Running {ds.upper()} Part A (+ Part B stream extraction) ===", flush=True)
            summary = run_dataset_extract_and_part_a(ds)
            print(json.dumps(summary["timings_s"], indent=2), flush=True)
    if args.stage in ("part-b", "all"):
        print("=== Running Part B Spread + Robustness + Verdicts + QC JSON ===", flush=True)
        qc = run_part_b_and_qc()
        print(json.dumps(qc["verdicts"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
