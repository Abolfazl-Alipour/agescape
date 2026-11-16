#!/usr/bin/env python3
"""Full-depth cell-level drift (Part 1), spread (Part 2), and ectopic identity (Part 3) in Gao Heart and Yazar Blood."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gc
import json
import os
from pathlib import Path
import pickle
import re
import resource
import shutil
import sys
import time
from typing import Any
import zlib

for _env_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_env_key, "1")

import h5py
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats as scipy_stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_cross_species_expression_vector import reference_without_donor, technical_gene_mask
from scripts.run_followup_between_subject_gao_yazar import (
    STUDY_DUMMY_MAP,
    _prune_covariate_matrix,
    apply_slice_bh,
    assign_age_bin,
    build_human_gene_projection,
    fit_contrast_adj_diff,
    fit_ols_model_row,
    read_csr_chunk,
    reconcile_gao_metadata,
    reconcile_yazar_metadata,
)
from scripts.run_followup_de_man_entropy import SEEDS as ENTROPY_SEEDS
from scripts.run_followup_de_man_entropy import kl_uniform, multinomial_downsample, shannon
from scripts.run_h3_expression_vector import hellinger_distance
from src.followup_muscle import count_gate
from src.followup_organs import write_clean_csv

DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
GAO_H5AD = DATA_ROOT / "organ17_gao_heart" / "GSE290367_integrated_RNA_adata.h5ad"
YAZAR_H5AD = DATA_ROOT / "hcata" / "yazar_1e44db10.h5ad"
YAZAR_CACHE = DATA_ROOT / "hcata" / "yazar_cache"
CACHE_DIR = DATA_ROOT / "fulldepth_cache"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"

MASTER_SEED = 42
SH_SEEDS = [20261116 + i for i in range(5)]
SUB20_SEEDS = [20262000 + i for i in range(20)]
SPLINE_KNOTS = (np.log(1500.0), np.log(2500.0), np.log(4500.0), np.log(9000.0))
DEPTH_STRATA = (
    ("s1_1500_2500", 1500, 2500),
    ("s2_2500_4500", 2500, 4500),
    ("s3_4500_9000", 4500, 9000),
    ("s4_9000_18000", 9000, 18000),
)
N_SHARDS = 16
CHUNK_ROWS = 16000
SUBCHUNK_ROWS = 1000
MAX_WORKERS = 6
SHARD_WORKERS = 3
MIN_CELLS_PER_DONOR = 20
MIN_DONORS_SLOPE = 6
MIN_YOUNG_REF = 2
MIN_LOO_REF_CELLS = 50
MIN_LOO_YOUNG_E2_STRATUM = 200
MAX_E2_GLOBAL_POOL_FRAC = 0.01

# Exact biological ambient/RBC/Ig/platelet/myeloid/IEG masks (sparing IGF*, HBEGF, HBS1L, etc.)
IG_RE = re.compile(r"^(IGH[ACDEGMVJ]|IGK[CVJ]|IGL[CVJ])")
HG_EXACT = frozenset({"HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"})
BLOOD_AMBIENT_EXACT = frozenset({"JCHAIN", "MZB1", "PPBP", "PF4", "S100A8", "S100A9", "LYZ"})
IEG_EXACT = frozenset({
    "FOS", "FOSB", "JUN", "JUNB", "JUND",
    "EGR1", "EGR2", "EGR3", "IER2", "IER3", "ZFP36", "ATF3", "DUSP1",
})
HEART_CM_CURATED = frozenset({"TTN", "MYH6", "MYH7", "RYR2", "TNNT2", "ACTN2", "MB", "CKM", "NPPA"})

# Worker globals for lightweight Pass 1/2/3 shard execution
_W_GPROJ: sp.csc_matrix | None = None
_W_CTX: dict[str, Any] = {}


def _log(tag: str, msg: str) -> None:
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(f"[{tag}] {msg} | self_max_rss={rss_mb:.1f} MB", flush=True)


def _worker_init(g_proj: sp.csc_matrix | None = None, ctx: dict[str, Any] | None = None) -> None:
    global _W_GPROJ, _W_CTX
    if g_proj is not None:
        _W_GPROJ = g_proj
    if ctx is not None:
        _W_CTX = ctx
    gc.freeze()


def build_stringent_ambient_mask(
    retained_symbols: list[str],
    ds_name: str,
    young_cm_counts: np.ndarray | None = None,
) -> np.ndarray:
    """Exact ambient/RBC/Ig/platelet/IEG exclusion (+ YOUNG CM top-50 ambient for Gao non-CM)."""
    mask = np.zeros(len(retained_symbols), dtype=bool)
    for j, sym in enumerate(retained_symbols):
        up = str(sym).upper().strip()
        if IG_RE.match(up) or up in HG_EXACT or up in BLOOD_AMBIENT_EXACT or up in IEG_EXACT:
            mask[j] = True
        if ds_name == "gao" and up in HEART_CM_CURATED:
            mask[j] = True
    if ds_name == "gao" and young_cm_counts is not None and young_cm_counts.sum() > 0:
        top50 = np.argsort(young_cm_counts)[::-1][:50]
        mask[top50] = True
    return mask


def assign_depth_stratum(umi: np.ndarray) -> np.ndarray:
    """Assign 4 pre-specified bounded depth strata (-1 if >= 18000 or < 1500)."""
    out = np.full(len(umi), -1, dtype=np.int32)
    for si, (_, lo, hi) in enumerate(DEPTH_STRATA):
        m = (umi >= lo) & (umi < hi)
        out[m] = si
    return out


def spline_design_matrix(umi: np.ndarray) -> np.ndarray:
    """Fixed-knot natural cubic spline on log(umi) (4 knots -> 3 basis cols) + 1/umi + intercept (5 cols)."""
    t = np.log(np.maximum(umi.astype(np.float64), 1.0))
    k1, k2, k3, k4 = SPLINE_KNOTS

    def _d(k: float) -> np.ndarray:
        return (np.maximum(t - k, 0.0) ** 3 - np.maximum(t - k4, 0.0) ** 3) / (k4 - k)

    b1 = _d(k1) - _d(k3)
    b2 = _d(k2) - _d(k3)
    inv_u = 1.0 / np.maximum(umi.astype(np.float64), 1.0)
    return np.column_stack([np.ones(len(t), dtype=np.float64), t, b1, b2, inv_u])


def fast_entropy_from_csr(mat: sp.csr_matrix, n_genes: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact bit-equivalent Shannon, KL_methods (= log(G)-Shannon), and KL_sqrt from CSR without triple copies."""
    n_rows = mat.shape[0]
    sh = np.full(n_rows, np.nan, dtype=np.float64)
    kl_m = np.full(n_rows, np.nan, dtype=np.float64)
    kl_s = np.full(n_rows, np.nan, dtype=np.float64)
    if n_rows == 0:
        return sh, kl_m, kl_s
    row_ids = np.repeat(np.arange(n_rows, dtype=np.int32), np.diff(mat.indptr))
    d = mat.data.astype(np.float64, copy=False)
    pos = d > 0.0
    if not np.any(pos):
        return sh, kl_m, kl_s
    r_p = row_ids[pos]
    v_p = d[pos]
    tots = np.bincount(r_p, weights=v_p, minlength=n_rows)
    valid = tots > 0.0
    p = v_p / tots[r_p]
    sh_vals = -np.bincount(r_p, weights=p * np.log(p), minlength=n_rows)
    sh[valid] = sh_vals[valid]
    kl_m[valid] = np.log(float(n_genes)) - sh_vals[valid]

    sq_v = np.sqrt(v_p)
    tots_sq = np.bincount(r_p, weights=sq_v, minlength=n_rows)
    q = sq_v / tots_sq[r_p]
    kl_s_vals = np.bincount(r_p, weights=q * np.log(q * float(n_genes)), minlength=n_rows)
    kl_s[valid] = kl_s_vals[valid]
    return sh, kl_m, kl_s


def fast_hellinger_csr_vec(mat: sp.csr_matrix, ref_sqrt: np.ndarray, already_l1: bool = False) -> np.ndarray:
    """Vectorized per-row Hellinger distance against one normalized sqrt-reference vector (atol < 1e-14 vs hellinger_distance)."""
    n_rows = mat.shape[0]
    if n_rows == 0:
        return np.zeros(0, dtype=np.float64)
    row_ids = np.repeat(np.arange(n_rows, dtype=np.int32), np.diff(mat.indptr))
    d = mat.data.astype(np.float64, copy=False)
    if not already_l1:
        tots = np.maximum(np.bincount(row_ids, weights=np.maximum(d, 0.0), minlength=n_rows), 1.0)
        p = np.maximum(d, 0.0) / tots[row_ids]
    else:
        p = np.maximum(d, 0.0)
    bc = np.bincount(row_ids, weights=np.sqrt(p) * ref_sqrt[mat.indices], minlength=n_rows)
    return np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc, 0.0, 1.0)))


def _shard_dir(ds_name: str) -> Path:
    return CACHE_DIR / f"{ds_name}_shards"


def _pass1_worker(task: tuple[Any, ...]) -> dict[str, Any]:
    """Pass 1 worker: compute Part 2 Full/1500/3000/5-SH/4-strata/spline + YOUNG counts & det + compact >=1500 CSR."""
    (
        data,
        indices,
        indptr,
        n_raw_vars,
        d_idx,
        c_idx,
        bin_arr,
        n_donors,
        n_cts,
        subchunk_id,
        ds_code,
        shard_cell_mask_bits,
    ) = task
    sub = sp.csr_matrix((data.astype(np.float64, copy=False), indices, indptr), shape=(len(d_idx), n_raw_vars))
    sub_clean = (sub @ _W_GPROJ).tocsr()
    del sub
    sub_clean.eliminate_zeros()

    n_genes = sub_clean.shape[1]
    umi = np.rint(np.asarray(sub_clean.sum(axis=1)).ravel()).astype(np.int64)
    m_15 = np.flatnonzero(umi >= 1500)
    out: dict[str, Any] = {"gate_n_all": len(d_idx), "n_15": len(m_15)}
    if len(m_15) == 0:
        return out

    mat_15 = sub_clean[m_15].tocsr()
    del sub_clean
    d_15 = d_idx[m_15].astype(np.int64)
    c_15 = c_idx[m_15].astype(np.int64)
    u_15 = umi[m_15]
    p_15 = d_15 * n_cts + c_15
    lu_15 = np.log(u_15.astype(np.float64))
    z_sp = spline_design_matrix(u_15)
    st_15 = assign_depth_stratum(u_15)

    u_pids, inv_p = np.unique(p_15, return_inverse=True)
    u_dids, inv_d = np.unique(d_15, return_inverse=True)
    n_lp, n_ld = len(u_pids), len(u_dids)

    # 1. Full-depth Part 2 (no seeds)
    sh_f, klm_f, kls_f = fast_entropy_from_csr(mat_15, n_genes)

    def _agg_vec(arr: np.ndarray, inv: np.ndarray, length: int) -> np.ndarray:
        return np.bincount(inv, weights=arr, minlength=length)

    pb: dict[str, Any] = {
        "u_pids": u_pids,
        "u_dids": u_dids,
        "cnt_ct": np.bincount(inv_p, minlength=n_lp),
        "cnt_pl": np.bincount(inv_d, minlength=n_ld),
        "lu_ct": _agg_vec(lu_15, inv_p, n_lp),
        "lu_pl": _agg_vec(lu_15, inv_d, n_ld),
        "invu_ct": _agg_vec(1.0 / u_15.astype(np.float64), inv_p, n_lp),
        "invu_pl": _agg_vec(1.0 / u_15.astype(np.float64), inv_d, n_ld),
        "full_ct": np.column_stack([_agg_vec(sh_f, inv_p, n_lp), _agg_vec(klm_f, inv_p, n_lp), _agg_vec(kls_f, inv_p, n_lp)]),
        "full_pl": np.column_stack([_agg_vec(sh_f, inv_d, n_ld), _agg_vec(klm_f, inv_d, n_ld), _agg_vec(kls_f, inv_d, n_ld)]),
    }

    # Spline sufficient stats (CT level Z^T Z & Z^T y + donor mean z vector across the 5 cols)
    u_cts, inv_c = np.unique(c_15, return_inverse=True)
    n_lc = len(u_cts)
    ztz_c = np.zeros((n_lc, 5, 5), dtype=np.float64)
    zty_c = np.zeros((n_lc, 3, 5), dtype=np.float64)
    for k_c in range(n_lc):
        mk = inv_c == k_c
        zk = z_sp[mk]
        ztz_c[k_c] = zk.T @ zk
        zty_c[k_c, 0] = zk.T @ sh_f[mk]
        zty_c[k_c, 1] = zk.T @ klm_f[mk]
        zty_c[k_c, 2] = zk.T @ kls_f[mk]
    zsum_ct = np.column_stack([_agg_vec(z_sp[:, j], inv_p, n_lp) for j in range(5)])
    zsum_pl = np.column_stack([_agg_vec(z_sp[:, j], inv_d, n_ld) for j in range(5)])
    pb["u_cts"] = u_cts
    pb["ztz_c"] = ztz_c
    pb["zty_c"] = zty_c
    pb["zsum_ct"] = zsum_ct
    pb["zsum_pl"] = zsum_pl

    # 4 depth strata Part 2 accumulators
    st_cnt_ct = np.zeros((4, n_lp), dtype=np.int64)
    st_lu_ct = np.zeros((4, n_lp), dtype=np.float64)
    st_sum_ct = np.zeros((4, n_lp, 3), dtype=np.float64)
    for si in range(4):
        sm = st_15 == si
        if np.any(sm):
            ip_s = inv_p[sm]
            st_cnt_ct[si] = np.bincount(ip_s, minlength=n_lp)
            st_lu_ct[si] = np.bincount(ip_s, weights=lu_15[sm], minlength=n_lp)
            st_sum_ct[si, :, 0] = np.bincount(ip_s, weights=sh_f[sm], minlength=n_lp)
            st_sum_ct[si, :, 1] = np.bincount(ip_s, weights=klm_f[sm], minlength=n_lp)
            st_sum_ct[si, :, 2] = np.bincount(ip_s, weights=kls_f[sm], minlength=n_lp)
    pb["st_cnt_ct"] = st_cnt_ct
    pb["st_lu_ct"] = st_lu_ct
    pb["st_sum_ct"] = st_sum_ct

    # 2. Vectorized 5-seed Binomial Split-Half for Part 2 (without copying CSR index structure!)
    row_ids_15 = np.repeat(np.arange(len(m_15), dtype=np.int32), np.diff(mat_15.indptr))
    d_int = np.rint(mat_15.data).astype(np.int64)
    sh_h1_ct = np.zeros((len(SH_SEEDS), n_lp, 3), dtype=np.float64)
    sh_h2_ct = np.zeros((len(SH_SEEDS), n_lp, 3), dtype=np.float64)
    sh_h1_pl = np.zeros((len(SH_SEEDS), n_ld, 3), dtype=np.float64)
    sh_h2_pl = np.zeros((len(SH_SEEDS), n_ld, 3), dtype=np.float64)

    def _entropy_half(vals: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        pos = vals > 0
        rp = row_ids_15[pos]
        vp = vals[pos].astype(np.float64)
        tot = np.maximum(np.bincount(rp, weights=vp, minlength=len(m_15)), 1.0)
        p = vp / tot[rp]
        s_v = -np.bincount(rp, weights=p * np.log(p), minlength=len(m_15))
        km_v = np.log(float(n_genes)) - s_v
        sq = np.sqrt(vp)
        tot_sq = np.maximum(np.bincount(rp, weights=sq, minlength=len(m_15)), 1e-12)
        q = sq / tot_sq[rp]
        ks_v = np.bincount(rp, weights=q * np.log(q * float(n_genes)), minlength=len(m_15))
        return s_v, km_v, ks_v

    for si_sh, seed_sh in enumerate(SH_SEEDS):
        rng_sh = np.random.default_rng(np.random.SeedSequence([int(seed_sh), ds_code, subchunk_id]))
        h1_v = rng_sh.binomial(d_int, 0.5)
        h2_v = d_int - h1_v
        s1, km1, ks1 = _entropy_half(h1_v)
        s2, km2, ks2 = _entropy_half(h2_v)
        for arr_c, arr_p, v_tup in ((sh_h1_ct[si_sh], sh_h1_pl[si_sh], (s1, km1, ks1)), (sh_h2_ct[si_sh], sh_h2_pl[si_sh], (s2, km2, ks2))):
            for m_i in range(3):
                arr_c[:, m_i] = _agg_vec(v_tup[m_i], inv_p, n_lp)
                arr_p[:, m_i] = _agg_vec(v_tup[m_i], inv_d, n_ld)
    pb["sh_h1_ct"] = sh_h1_ct
    pb["sh_h2_ct"] = sh_h2_ct
    pb["sh_h1_pl"] = sh_h1_pl
    pb["sh_h2_pl"] = sh_h2_pl

    # 3. 1500-UMI rarefaction sensitivity (5 Sh / 20 KL seeds) + 3000-UMI selection-bias comparison
    m_3k = np.flatnonzero(u_15 >= 3000)
    n_s = len(ENTROPY_SEEDS)
    r15_ct = np.zeros((n_s, n_lp, 3), dtype=np.float64)
    r15_pl = np.zeros((n_s, n_ld, 3), dtype=np.float64)
    r15_on_3k_ct = np.zeros((n_lp, 3), dtype=np.float64)
    r15_on_3k_pl = np.zeros((n_ld, 3), dtype=np.float64)
    mat_15_i64 = mat_15.astype(np.int64)
    for si in range(n_s):
        rng_15 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, 1500, subchunk_id, si]))
        ds15 = multinomial_downsample(mat_15_i64, 1500, rng_15)
        ds15.eliminate_zeros()  # downsample keeps dropped genes as explicit stored zeros; drop them before any stored-entry counting
        s_v, km_v, ks_v = fast_entropy_from_csr(ds15, n_genes)
        if si < 5:
            r15_ct[si, :, 0] = _agg_vec(s_v, inv_p, n_lp)
            r15_pl[si, :, 0] = _agg_vec(s_v, inv_d, n_ld)
        r15_ct[si, :, 1] = _agg_vec(km_v, inv_p, n_lp)
        r15_pl[si, :, 1] = _agg_vec(km_v, inv_d, n_ld)
        r15_ct[si, :, 2] = _agg_vec(ks_v, inv_p, n_lp)
        r15_pl[si, :, 2] = _agg_vec(ks_v, inv_d, n_ld)
        if si == 0 and len(m_3k):
            r15_on_3k_ct[:, 0] = np.bincount(inv_p[m_3k], weights=s_v[m_3k], minlength=n_lp)
            r15_on_3k_pl[:, 0] = np.bincount(inv_d[m_3k], weights=s_v[m_3k], minlength=n_ld)
            r15_on_3k_ct[:, 1] = np.bincount(inv_p[m_3k], weights=km_v[m_3k], minlength=n_lp)
            r15_on_3k_pl[:, 1] = np.bincount(inv_d[m_3k], weights=km_v[m_3k], minlength=n_ld)
            r15_on_3k_ct[:, 2] = np.bincount(inv_p[m_3k], weights=ks_v[m_3k], minlength=n_lp)
            r15_on_3k_pl[:, 2] = np.bincount(inv_d[m_3k], weights=ks_v[m_3k], minlength=n_ld)
        del ds15
    pb["r15_ct"] = r15_ct
    pb["r15_pl"] = r15_pl
    pb["r15_on_3k_ct"] = r15_on_3k_ct
    pb["r15_on_3k_pl"] = r15_on_3k_pl

    # 3000-UMI rarefaction sensitivity
    r3k_ct = np.zeros((n_s, n_lp, 3), dtype=np.float64)
    r3k_pl = np.zeros((n_s, n_ld, 3), dtype=np.float64)
    pb["cnt3k_ct"] = np.bincount(inv_p[m_3k], minlength=n_lp) if len(m_3k) else np.zeros(n_lp, dtype=np.int64)
    pb["cnt3k_pl"] = np.bincount(inv_d[m_3k], minlength=n_ld) if len(m_3k) else np.zeros(n_ld, dtype=np.int64)
    pb["lu3k_ct"] = np.bincount(inv_p[m_3k], weights=lu_15[m_3k], minlength=n_lp) if len(m_3k) else np.zeros(n_lp, dtype=np.float64)
    pb["lu3k_pl"] = np.bincount(inv_d[m_3k], weights=lu_15[m_3k], minlength=n_ld) if len(m_3k) else np.zeros(n_ld, dtype=np.float64)
    if len(m_3k):
        mat_3k_i64 = mat_15_i64[m_3k]
        ip_3k = inv_p[m_3k]
        id_3k = inv_d[m_3k]
        for si in range(n_s):
            rng_3k = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, 3000, subchunk_id, si]))
            ds3k = multinomial_downsample(mat_3k_i64, 3000, rng_3k)
            ds3k.eliminate_zeros()  # downsample keeps dropped genes as explicit stored zeros; drop them before any stored-entry counting
            s_v, km_v, ks_v = fast_entropy_from_csr(ds3k, n_genes)
            if si < 5:
                r3k_ct[si, :, 0] = _agg_vec(s_v, ip_3k, n_lp)
                r3k_pl[si, :, 0] = _agg_vec(s_v, id_3k, n_ld)
            r3k_ct[si, :, 1] = _agg_vec(km_v, ip_3k, n_lp)
            r3k_pl[si, :, 1] = _agg_vec(km_v, id_3k, n_ld)
            r3k_ct[si, :, 2] = _agg_vec(ks_v, ip_3k, n_lp)
            r3k_pl[si, :, 2] = _agg_vec(ks_v, id_3k, n_ld)
            del ds3k
    pb["r3k_ct"] = r3k_ct
    pb["r3k_pl"] = r3k_pl

    # 4. Accumulate YOUNG cell counts & gene detections for Part 1 & Part 3 references (CSR tuples)
    y_mask = bin_arr[d_15] == 1
    ref_chunks: dict[str, Any] = {}
    if np.any(y_mask):
        mat_y = mat_15[y_mask]
        p_y = p_15[y_mask]
        u_py, inv_py = np.unique(p_y, return_inverse=True)
        ind_y = sp.csr_matrix((np.ones(len(inv_py), dtype=np.float64), (inv_py, np.arange(len(inv_py)))), shape=(len(u_py), len(inv_py)))
        sum_cnt_y = (ind_y @ mat_y).tocsr()
        sum_cnt_y.sum_duplicates()
        mat_y_bin = mat_y.copy()
        mat_y_bin.data = np.ones_like(mat_y_bin.data)
        sum_det_y = (ind_y @ mat_y_bin).tocsr()
        sum_det_y.sum_duplicates()
        ref_chunks["y_full_cnt"] = (u_py, np.rint(sum_cnt_y.data).astype(np.int64), sum_cnt_y.indices, sum_cnt_y.indptr)
        ref_chunks["y_full_det"] = (u_py, np.rint(sum_det_y.data).astype(np.int64), sum_det_y.indices, sum_det_y.indptr)
        ref_chunks["y_n_cells"] = (u_py, np.bincount(inv_py))

        # 4 bounded depth strata YOUNG dets for E2_stratum
        for si in range(4):
            ym_s = y_mask & (st_15 == si)
            if np.any(ym_s):
                my_s = mat_15[ym_s].copy()
                my_s.data = np.ones_like(my_s.data)
                u_ps, inv_ps = np.unique(p_15[ym_s], return_inverse=True)
                ind_s = sp.csr_matrix((np.ones(len(inv_ps), dtype=np.float64), (inv_ps, np.arange(len(inv_ps)))), shape=(len(u_ps), len(inv_ps)))
                sd = (ind_s @ my_s).tocsr()
                sd.sum_duplicates()
                ref_chunks[f"y_st{si}_det"] = (u_ps, np.rint(sd.data).astype(np.int64), sd.indices, sd.indptr)
                ref_chunks[f"y_st{si}_n"] = (u_ps, np.bincount(inv_ps))

        # Seed-0 1500 YOUNG counts & 5-seed 1500 YOUNG dets for Part 1 & Part 3 rarefied references
        my_i64 = mat_15_i64[y_mask]
        for si_e in range(5):
            rng_ye = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, 1500, subchunk_id, si_e]))
            ds_ye = multinomial_downsample(my_i64, 1500, rng_ye)
            ds_ye.eliminate_zeros()  # downsample keeps dropped genes as explicit stored zeros; drop them before detection binarization below
            if si_e == 0:
                sc15 = (ind_y @ ds_ye).tocsr()
                sc15.sum_duplicates()
                ref_chunks["y_1500_cnt0"] = (u_py, np.rint(sc15.data).astype(np.int64), sc15.indices, sc15.indptr)
            ds_ye.data = np.ones_like(ds_ye.data)
            sd15 = (ind_y @ ds_ye).tocsr()
            sd15.sum_duplicates()
            ref_chunks[f"y_1500_det_{si_e}"] = (u_py, np.rint(sd15.data).astype(np.int64), sd15.indices, sd15.indptr)

    # 5. Prepare serialized K=16 shard payloads directly (int32 data, uint16 indices)
    sh_payloads: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    shards = d_15 % N_SHARDS
    for s_id in np.unique(shards):
        sm = shards == s_id
        sub_s = mat_15[sm].tocsr()
        nnz_per_row = np.diff(sub_s.indptr).astype(np.int64)
        cells_hdr = np.column_stack([d_15[sm], c_15[sm], u_15[sm], nnz_per_row]).astype(np.int64)
        sh_payloads[int(s_id)] = (cells_hdr, sub_s.indices.astype(np.uint16), np.rint(sub_s.data).astype(np.int32))

    out["partb"] = pb
    out["refs"] = ref_chunks
    out["shards"] = sh_payloads
    out["subchunk_id"] = subchunk_id
    return out


def _scatter_csr_rows(target: np.ndarray, csr_tuple: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None, row_map: dict[int, int]) -> None:
    if csr_tuple is None:
        return
    u_m, data, indices, indptr = csr_tuple
    for k_i, g_id in enumerate(u_m):
        dst = row_map.get(int(g_id), -1)
        if dst < 0:
            continue
        s0, s1 = int(indptr[k_i]), int(indptr[k_i + 1])
        target[dst, indices[s0:s1]] += data[s0:s1].astype(target.dtype, copy=False)


def run_pass1_and_export_part2(meta: dict[str, Any], ds_name: str) -> dict[str, Any]:
    """Single HDF5 scan per dataset: writes K=16 CSR shards, Part 1 & 3 YOUNG refs, and exports Part 2 CSV."""
    t0 = time.perf_counter()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    sdir = _shard_dir(ds_name)
    if sdir.exists():
        shutil.rmtree(sdir)
    sdir.mkdir(parents=True, exist_ok=True)

    proj_info = build_human_gene_projection(meta["gene_symbols"])
    g_proj = proj_info["proj"]
    n_genes = g_proj.shape[1]
    assert n_genes < 65535, f"n_genes={n_genes} exceeds uint16 index bound!"
    ds_code = 1 if ds_name == "gao" else 2

    donors_df = meta["donors_df"]
    n_donors = len(meta["donor_ids"])
    n_cts = len(meta["cell_types"])
    n_pairs = n_donors * n_cts
    bin_arr = donors_df["bin"].map({"LT20": 0, "YOUNG": 1, "MID": 2, "OLD": 3}).to_numpy(dtype=np.int32)

    # Identify all (YOUNG d, T) pairs upfront cleanly (at most 20*13=260 Gao or 103*28=2,884 Yazar)
    y_donors = np.flatnonzero(bin_arr == 1).tolist()
    y_pair_map: dict[int, int] = {}
    y_pairs_list: list[tuple[int, int]] = []
    for d in y_donors:
        for c in range(n_cts):
            pid = d * n_cts + c
            y_pair_map[pid] = len(y_pairs_list)
            y_pairs_list.append((d, c))
    n_yp = len(y_pairs_list)

    n_s = len(ENTROPY_SEEDS)
    n_sh = len(SH_SEEDS)
    acc = {
        "cnt_ct": np.zeros(n_pairs, dtype=np.int64),
        "cnt_pl": np.zeros(n_donors, dtype=np.int64),
        "cnt3k_ct": np.zeros(n_pairs, dtype=np.int64),
        "cnt3k_pl": np.zeros(n_donors, dtype=np.int64),
        "lu_ct": np.zeros(n_pairs, dtype=np.float64),
        "lu_pl": np.zeros(n_donors, dtype=np.float64),
        "lu3k_ct": np.zeros(n_pairs, dtype=np.float64),
        "lu3k_pl": np.zeros(n_donors, dtype=np.float64),
        "invu_ct": np.zeros(n_pairs, dtype=np.float64),
        "invu_pl": np.zeros(n_donors, dtype=np.float64),
        "full_ct": np.zeros((n_pairs, 3), dtype=np.float64),
        "full_pl": np.zeros((n_donors, 3), dtype=np.float64),
        "ztz_c": np.zeros((n_cts, 5, 5), dtype=np.float64),
        "zty_c": np.zeros((n_cts, 3, 5), dtype=np.float64),
        "zsum_ct": np.zeros((n_pairs, 5), dtype=np.float64),
        "zsum_pl": np.zeros((n_donors, 5), dtype=np.float64),
        "st_cnt_ct": np.zeros((4, n_pairs), dtype=np.int64),
        "st_lu_ct": np.zeros((4, n_pairs), dtype=np.float64),
        "st_sum_ct": np.zeros((4, n_pairs, 3), dtype=np.float64),
        "sh_h1_ct": np.zeros((n_sh, n_pairs, 3), dtype=np.float64),
        "sh_h2_ct": np.zeros((n_sh, n_pairs, 3), dtype=np.float64),
        "sh_h1_pl": np.zeros((n_sh, n_donors, 3), dtype=np.float64),
        "sh_h2_pl": np.zeros((n_sh, n_donors, 3), dtype=np.float64),
        "r15_ct": np.zeros((n_s, n_pairs, 3), dtype=np.float64),
        "r15_pl": np.zeros((n_s, n_donors, 3), dtype=np.float64),
        "r15_on_3k_ct": np.zeros((n_pairs, 3), dtype=np.float64),
        "r15_on_3k_pl": np.zeros((n_donors, 3), dtype=np.float64),
        "r3k_ct": np.zeros((n_s, n_pairs, 3), dtype=np.float64),
        "r3k_pl": np.zeros((n_s, n_donors, 3), dtype=np.float64),
    }

    y_full_cnt = np.zeros((n_yp, n_genes), dtype=np.int64)
    y_full_det = np.zeros((n_yp, n_genes), dtype=np.int64)
    y_n_cells = np.zeros(n_yp, dtype=np.int64)
    y_st_det = np.zeros((4, n_yp, n_genes), dtype=np.int32)
    y_st_n = np.zeros((4, n_yp), dtype=np.int64)
    y_1500_cnt0 = np.zeros((n_yp, n_genes), dtype=np.int32)
    y_1500_det = np.zeros((5, n_yp, n_genes), dtype=np.int32)

    f_cells = [open(sdir / f"shard_{s}_cells.bin", "wb") for s in range(N_SHARDS)]
    f_idx = [open(sdir / f"shard_{s}_idx.bin", "wb") for s in range(N_SHARDS)]
    f_dat = [open(sdir / f"shard_{s}_dat.bin", "wb") for s in range(N_SHARDS)]
    shard_cell_counts = [0] * N_SHARDS
    total_ge1500 = 0

    try:
        h5_path = meta["h5_path"]
        slot = meta["matrix_slot"]
        indptr = meta["indptr"]
        keep_mask = meta["keep_mask"]
        c_donor = meta["cell_donor_idx"]
        c_ct = meta["cell_ct_idx"]
        n_raw_vars = len(meta["gene_symbols"])

        with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init, initargs=(g_proj, {})) as ex:
            with h5py.File(h5_path, "r") as f:
                grp = f["layers"]["counts"] if slot == "layers/counts" else (f["raw"]["X"] if slot == "raw/X" else f["X"])
                n_obs = meta["n_obs"]
                futures = []
                sub_id = 0

                def _drain(lim: int) -> None:
                    nonlocal total_ge1500
                    while len(futures) > lim:
                        r = futures.pop(0).result()
                        if r["n_15"] == 0:
                            continue
                        total_ge1500 += r["n_15"]
                        pb = r["partb"]
                        up = pb["u_pids"]
                        ud = pb["u_dids"]
                        uc = pb["u_cts"]
                        for k in ("cnt_ct", "cnt3k_ct", "lu_ct", "lu3k_ct", "invu_ct", "full_ct", "zsum_ct", "r15_on_3k_ct"):
                            acc[k][up] += pb[k]
                        for k in ("cnt_pl", "cnt3k_pl", "lu_pl", "lu3k_pl", "invu_pl", "full_pl", "zsum_pl", "r15_on_3k_pl"):
                            acc[k][ud] += pb[k]
                        acc["ztz_c"][uc] += pb["ztz_c"]
                        acc["zty_c"][uc] += pb["zty_c"]
                        for k in ("st_cnt_ct", "st_lu_ct", "st_sum_ct", "sh_h1_ct", "sh_h2_ct", "r15_ct", "r3k_ct"):
                            acc[k][:, up] += pb[k]
                        for k in ("sh_h1_pl", "sh_h2_pl", "r15_pl", "r3k_pl"):
                            acc[k][:, ud] += pb[k]

                        rf = r["refs"]
                        if "y_full_cnt" in rf:
                            _scatter_csr_rows(y_full_cnt, rf["y_full_cnt"], y_pair_map)
                            _scatter_csr_rows(y_full_det, rf["y_full_det"], y_pair_map)
                            u_py, nc_py = rf["y_n_cells"]
                            for g_id, cnt_v in zip(u_py, nc_py):
                                y_n_cells[y_pair_map[int(g_id)]] += int(cnt_v)
                            for si in range(4):
                                if f"y_st{si}_det" in rf:
                                    _scatter_csr_rows(y_st_det[si], rf[f"y_st{si}_det"], y_pair_map)
                                    u_ps, nc_ps = rf[f"y_st{si}_n"]
                                    for g_id, cnt_v in zip(u_ps, nc_ps):
                                        y_st_n[si, y_pair_map[int(g_id)]] += int(cnt_v)
                            _scatter_csr_rows(y_1500_cnt0, rf["y_1500_cnt0"], y_pair_map)
                            for si_e in range(5):
                                _scatter_csr_rows(y_1500_det[si_e], rf[f"y_1500_det_{si_e}"], y_pair_map)

                        for s_id, (hdr, idx_a, dat_a) in r["shards"].items():
                            hdr.tofile(f_cells[s_id])
                            idx_a.tofile(f_idx[s_id])
                            dat_a.tofile(f_dat[s_id])
                            shard_cell_counts[s_id] += len(hdr)

                for r0 in range(0, n_obs, CHUNK_ROWS):
                    r1 = min(r0 + CHUNK_ROWS, n_obs)
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
                                _pass1_worker,
                                (
                                    sc.data,
                                    sc.indices,
                                    sc.indptr,
                                    n_raw_vars,
                                    d_s[s0:s1],
                                    c_s[s0:s1],
                                    bin_arr,
                                    n_donors,
                                    n_cts,
                                    sub_id,
                                    ds_code,
                                    0,
                                ),
                            )
                        )
                        sub_id += 1
                        _drain(MAX_WORKERS)
                _drain(0)
    finally:
        for flist in (f_cells, f_idx, f_dat):
            for fh in flist:
                fh.close()

    elapsed = round(time.perf_counter() - t0, 2)
    _log(ds_name, f"Pass 1 single HDF5 scan done in {elapsed}s (kept >=1500 cells: {total_ge1500})")

    # Save compact Part 1 & Part 3 references externally (never touching HDF5 again!)
    part1_refs = {
        "y_pairs_list": y_pairs_list,
        "y_pair_map": y_pair_map,
        "y_full_cnt": y_full_cnt,
        "y_1500_cnt0": y_1500_cnt0,
        "y_n_cells": y_n_cells,
        "cnt_ct_15": acc["cnt_ct"].reshape(n_donors, n_cts),
        "cnt_pl_15": acc["cnt_pl"],
        "lu_ct_15": acc["lu_ct"].reshape(n_donors, n_cts),
        "lu_pl_15": acc["lu_pl"],
    }
    part3_refs = {
        "gene_audit": proj_info["audit"],
        "y_pairs_list": y_pairs_list,
        "y_pair_map": y_pair_map,
        "y_full_cnt": y_full_cnt,
        "y_full_det": y_full_det,
        "y_n_cells": y_n_cells,
        "y_st_det": y_st_det,
        "y_st_n": y_st_n,
        "y_1500_cnt0": y_1500_cnt0,
        "y_1500_det": y_1500_det,
        "cnt_ct_15": acc["cnt_ct"].reshape(n_donors, n_cts),
        "lu_ct_15": acc["lu_ct"].reshape(n_donors, n_cts),
    }
    part2_ckpt = {"acc": acc, "shard_cell_counts": shard_cell_counts, "total_ge1500": total_ge1500, "elapsed_pass1_s": elapsed}

    with open(CACHE_DIR / f"{ds_name}_pass1_part1_refs.pkl", "wb") as pf:
        pickle.dump(part1_refs, pf, protocol=pickle.HIGHEST_PROTOCOL)
    with open(CACHE_DIR / f"{ds_name}_pass1_part3_refs.pkl", "wb") as pf:
        pickle.dump(part3_refs, pf, protocol=pickle.HIGHEST_PROTOCOL)
    with open(CACHE_DIR / f"{ds_name}_pass1_part2.pkl", "wb") as pf:
        pickle.dump(part2_ckpt, pf, protocol=pickle.HIGHEST_PROTOCOL)

    del y_full_cnt, y_full_det, y_st_det, y_1500_cnt0, y_1500_det, part1_refs, part3_refs
    gc.collect()
    return part2_ckpt


def _model_specs_for_dataset(base_df: pd.DataFrame, ds_name: str, extra_log_ncells: bool = False) -> list[tuple[str, pd.DataFrame, list[str], str | None, int, int, str | None]]:
    """Canonical sensitivity specs reused uniformly across Parts 1, 2, 3."""
    if ds_name == "gao":
        bcol, ref_lv, min_bd, min_df = "study_or_pool", "ENCODE v4 (Snyder)", 2, 5
        specs = [
            ("primary", base_df, ["sex_male", "_depth"], bcol, min_bd, min_df, ref_lv),
            ("batch_adj_no_depth", base_df, ["sex_male"], bcol, min_bd, min_df, ref_lv),
            ("nonlinear_depth", base_df, ["sex_male", "_depth", "_depth2", "_inv_umi"], bcol, min_bd, min_df, ref_lv),
            ("age_only", base_df, [], None, 1, 5, None),
            ("age_plus_sex", base_df, ["sex_male"], None, 1, 5, None),
            ("10x_only_no_penn", base_df[base_df["technology"] != "Dropseq"].copy(), ["sex_male", "_depth"], bcol, min_bd, min_df, ref_lv),
            ("studies_ge8", base_df[~base_df["study_or_pool"].isin(["Kuppe 2022", "Sim 2021"])].copy(), ["sex_male", "_depth"], bcol, min_bd, min_df, ref_lv),
            ("encode_v4_only", base_df[base_df["study_or_pool"] == "ENCODE v4 (Snyder)"].copy(), ["sex_male", "_depth"], None, 1, 5, None),
        ]
    else:
        bcol, ref_lv, min_bd, min_df = "study_or_pool", None, 2, 10
        specs = [
            ("primary", base_df, ["sex_male", "_depth"], bcol, min_bd, min_df, None),
            ("pool_fe_no_depth", base_df, ["sex_male"], bcol, min_bd, min_df, None),
            ("nonlinear_depth", base_df, ["sex_male", "_depth", "_depth2", "_inv_umi"], bcol, min_bd, min_df, None),
            ("no_pool", base_df, ["sex_male", "_depth"], None, 1, 5, None),
            ("age_ge_20", base_df[base_df["age_years"] >= 20.0].copy(), ["sex_male", "_depth"], bcol, min_bd, min_df, None),
        ]
    if extra_log_ncells:
        specs.append(("plus_log_ncells", base_df, ["sex_male", "_depth", "_log_n_cells"], bcol, min_bd, min_df, ref_lv))
    return specs


def fit_stratified_global_with_batch(
    strat_df: pd.DataFrame,
    elig_cts: list[str],
    out_col: str,
    cell_col: str,
    depth_col: str,
    ds_name: str,
    part: str,
    metric: str,
    depth_mode: str,
    analysis_role: str,
    model_spec: str = "stratified_global_primary",
    gene_set: str = "none",
    strata_col: str | None = None,
    batch_col: str | None = "study_or_pool",
    ref_level: str | None = "ENCODE v4 (Snyder)",
    extra_covs: list[str] | None = None,
    depth_x_ct: bool = False,
    custom_scale_map: dict[Any, float] | None = None,
    fallback_exclude_groups: set[str] | None = None,
) -> tuple[dict[str, Any] | None, pd.DataFrame | None]:
    """Repaired donor-clustered Stratified Global OLS INCLUDING study/pool FE and OLD_vs_YOUNG contrast."""
    if len(elig_cts) < 2 or strat_df.empty:
        return None, None
    is_gao = ds_name == "gao"
    use_cts = [c for c in elig_cts if fallback_exclude_groups is None or c not in fallback_exclude_groups]
    if len(use_cts) < 2:
        return None, None

    group_key = ["group", strata_col] if strata_col else ["group"]
    parts = []
    scale_used: dict[Any, float] = {}
    for gk, cpart in strat_df[strat_df["group"].isin(use_cts)].dropna(subset=[out_col]).groupby(group_key):
        if len(cpart) < MIN_DONORS_SLOPE or cpart["age_per_decade"].nunique() < 2:
            continue
        sd_y = custom_scale_map[gk] if (custom_scale_map is not None and gk in custom_scale_map) else float(cpart[out_col].std(ddof=1))
        if not np.isfinite(sd_y) or sd_y <= 0.0:
            continue
        scale_used[gk] = sd_y
        cp = cpart.copy()
        cp["y_std"] = (cp[out_col] - cp[out_col].mean()) / sd_y
        cp["sex_c"] = cp["sex_male"] - cp["sex_male"].mean()
        cp["depth_c"] = cp[depth_col] - cp[depth_col].mean()
        cp["_gkey"] = "|".join(gk) if isinstance(gk, tuple) else str(gk)
        if extra_covs:
            for ec in extra_covs:
                cp[f"{ec}_c"] = cp[ec] - cp[ec].mean()
        parts.append(cp)
    if len(parts) < 2:
        return None, None
    df_st = pd.concat(parts, ignore_index=True)
    if df_st["donor_id"].nunique() < MIN_DONORS_SLOPE:
        return None, None

    # Enforce >=2 donors/batch cleanly over unique participating donors
    active_covs: list[str] = []
    batch_cols: list[str] = []
    if batch_col is not None and batch_col in df_st.columns:
        d_batches = df_st[["donor_id", batch_col]].drop_duplicates()
        vc = d_batches[batch_col].value_counts()
        ge_b = vc[vc >= 2].index.tolist()
        df_st = df_st[df_st[batch_col].isin(ge_b)].copy()
        vc2 = df_st[["donor_id", batch_col]].drop_duplicates()[batch_col].value_counts()
        if len(vc2) >= 2:
            lv_sorted = [k for k, _ in sorted(vc2.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))]
            r_lv = ref_level if (is_gao and ref_level in lv_sorted) else lv_sorted[0]
            for lv in lv_sorted:
                if lv != r_lv:
                    cname = f"batch_{str(lv).replace(' ', '_')}"
                    df_st[cname] = (df_st[batch_col].astype(str) == str(lv)).astype(float)
                    batch_cols.append(cname)

    g_dum = pd.get_dummies(df_st["_gkey"], prefix="g", drop_first=True, dtype=float)
    design_list = [np.ones(len(df_st), dtype=float), df_st["age_per_decade"].to_numpy(dtype=float)]
    x_blocks = [df_st[["age_per_decade"]]]

    cand_covs = ["sex_c", "depth_c"] + ([f"{ec}_c" for ec in extra_covs] if extra_covs else [])
    for cc in cand_covs:
        vals = df_st[cc].to_numpy(dtype=float)
        if np.unique(vals).size >= 2 and np.linalg.matrix_rank(np.column_stack([*design_list, vals])) == len(design_list) + 1:
            design_list.append(vals)
            x_blocks.append(df_st[[cc]])
            active_covs.append(cc)

    if depth_x_ct and "depth_c" in active_covs and len(g_dum.columns) > 0:
        d_vals = df_st["depth_c"].to_numpy(dtype=float)
        for col in g_dum.columns:
            inter = d_vals * g_dum[col].to_numpy(dtype=float)
            if np.linalg.matrix_rank(np.column_stack([*design_list, inter])) == len(design_list) + 1:
                design_list.append(inter)
                x_blocks.append(pd.DataFrame({f"dx_{col}": inter}, index=df_st.index))
        active_covs.append("depth_x_ct")

    for gc_col in g_dum.columns:
        vals = g_dum[gc_col].to_numpy(dtype=float)
        if np.linalg.matrix_rank(np.column_stack([*design_list, vals])) == len(design_list) + 1:
            design_list.append(vals)
            x_blocks.append(g_dum[[gc_col]])
    active_covs.append("ct_stratum_dummies")

    kept_batch = []
    for bc in batch_cols:
        vals = df_st[bc].to_numpy(dtype=float)
        if np.linalg.matrix_rank(np.column_stack([*design_list, vals])) == len(design_list) + 1:
            design_list.append(vals)
            x_blocks.append(df_st[[bc]])
            kept_batch.append(bc)
    if kept_batch:
        active_covs.append(f"batch_fe({len(kept_batch)})")

    X_st = sm.add_constant(pd.concat(x_blocks, axis=1), has_constant="add")
    fit_st = sm.OLS(df_st["y_std"], X_st).fit(cov_type="cluster", cov_kwds={"groups": df_st["donor_id"]}, use_t=True)
    df_c = max(1, int(df_st["donor_id"].nunique() - 1))
    t975 = float(scipy_stats.t.ppf(0.975, df_c))
    t80 = float(scipy_stats.t.ppf(0.80, df_c))
    sl = float(fit_st.params["age_per_decade"])
    se = float(fit_st.bse["age_per_decade"])

    # Fit exact matched OLD_vs_YOUNG contrast on this same design
    oy_diff, oy_se, oy_p = np.nan, np.nan, np.nan
    oy_mask = df_st["bin"].isin(["YOUNG", "OLD"])
    sub_oy = df_st[oy_mask].copy()
    if sub_oy[sub_oy["bin"] == "YOUNG"]["donor_id"].nunique() >= 2 and sub_oy[sub_oy["bin"] == "OLD"]["donor_id"].nunique() >= 2:
        is_old_s = pd.Series((sub_oy["bin"] == "OLD").astype(float), index=sub_oy.index, name="_is_old")
        x_oy_pool = X_st.loc[sub_oy.index]
        oy_cols = [c for c in x_oy_pool.columns if c not in ("const", "age_per_decade")]
        keep_oy_c: list[str] = []
        d_oy = [np.ones(len(sub_oy), dtype=float), is_old_s.to_numpy(dtype=float)]
        for c in oy_cols:
            v = x_oy_pool[c].to_numpy(dtype=float)
            if np.unique(v).size >= 2 and np.linalg.matrix_rank(np.column_stack([*d_oy, v])) == len(d_oy) + 1:
                d_oy.append(v)
                keep_oy_c.append(c)
        X_oy = sm.add_constant(pd.concat([is_old_s, x_oy_pool[keep_oy_c].astype(float)], axis=1), has_constant="add")
        fit_oy = sm.OLS(sub_oy["y_std"], X_oy).fit(cov_type="cluster", cov_kwds={"groups": sub_oy["donor_id"]}, use_t=True)
        oy_diff = float(fit_oy.params["_is_old"])
        oy_se = float(fit_oy.bse["_is_old"])
        oy_p = float(fit_oy.pvalues["_is_old"])

    row = {
        "dataset": ds_name,
        "part": part,
        "metric": metric,
        "gene_set": gene_set,
        "depth_mode": depth_mode,
        "analysis_role": analysis_role,
        "model_spec": model_spec,
        "group": "Stratified_CellType_Global",
        "is_pooled": True,
        "is_stratified_global": True,
        "yazar_pool_fe_status": "evaluated" if kept_batch else "no_batch",
        "n_donors": int(df_st["donor_id"].nunique()),
        "n_cells": int(df_st[cell_col].sum()),
        "residual_df": df_c,
        "covariates_active": ";".join(active_covs),
        "covariates_omitted": "none",
        "outcome_sd": 1.0,
        "resid_sd": float(np.sqrt(fit_st.mse_resid)),
        "slope_raw_per_decade": sl,
        "se_raw": se,
        "ci95_low_raw": sl - t975 * se,
        "ci95_high_raw": sl + t975 * se,
        "slope_sd_per_decade": sl,
        "se_sd": se,
        "ci95_low_sd": sl - t975 * se,
        "ci95_high_sd": sl + t975 * se,
        "pvalue": float(fit_st.pvalues["age_per_decade"]),
        "q_bh_within_slice": np.nan,
        "family_size_m": np.nan,
        "mde80_raw_per_decade": (t975 + t80) * se,
        "mde80_sd_per_decade": (t975 + t80) * se,
        "old_minus_young_diff_adj": oy_diff,
        "old_minus_young_se_adj": oy_se,
        "old_minus_young_p_adj": oy_p,
    }
    return row, scale_used


def compute_sh_and_quenouille_stats(
    df_group: pd.DataFrame,
    obs_col: str,
    h1_cols: list[str],
    h2_cols: list[str],
    nonbatch_covs: list[str],
    batch_col: str | None,
    min_batch_donors: int,
    min_residual_df: int,
    ref_level: str | None,
    is_part3_count_metric: bool,
    apply_quenouille: bool,
    fallback_to_nobatch: bool = False,
) -> dict[str, float]:
    """Compute exact Sampling SE (factor 0.5 for scale-invariant, 1.0 for raw count) and Quenouille Bias (Parts 1-2 only)."""
    nan_d = {
        "sh_rms_delta_donor": np.nan,
        "sh_se_samp_slope": np.nan,
        "sh_se_samp_oy": np.nan,
        "sh_rho_tech2": np.nan,
        "sh_snr_slope": np.nan,
        "sh_snr_oy": np.nan,
        "quenouille_bias_slope": np.nan,
        "quenouille_bias_frac": np.nan,
        "quenouille_corr_slope": np.nan,
        "quenouille_corr_p": np.nan,
    }
    sub_base = df_group.dropna(subset=[obs_col]).copy()
    if len(sub_base) < MIN_DONORS_SLOPE:
        return nan_d
    # Fit primary once cleanly to obtain residual mse_resid and obs_slope
    pri_row = fit_ols_model_row(
        sub_base, obs_col, "age_per_decade", nonbatch_covs, batch_col, min_batch_donors, min_residual_df, ref_level, {}, fallback_to_nobatch
    )
    if pri_row is None:
        return nan_d
    mse_res = float(pri_row["resid_sd"]) ** 2
    obs_sl = float(pri_row["slope_raw_per_decade"])
    oy_obs = fit_contrast_adj_diff(sub_base, obs_col, "YOUNG", "OLD", nonbatch_covs, batch_col, min_batch_donors, min_residual_df, ref_level)

    factor = 1.0 if is_part3_count_metric else 0.5
    var_scale = 1.0 if is_part3_count_metric else 0.25

    se_sl_sq, se_oy_sq, rms_sq, bias_mats = [], [], [], []
    for h1c, h2c in zip(h1_cols, h2_cols):
        if h1c not in sub_base.columns or h2c not in sub_base.columns:
            continue
        w = sub_base.dropna(subset=[h1c, h2c]).copy()
        if len(w) < MIN_DONORS_SLOPE:
            continue
        delta = w[h1c].to_numpy(float) - w[h2c].to_numpy(float)
        rms_sq.append(float(np.mean(delta ** 2)))
        w["_delta"] = delta
        if np.std(delta, ddof=1) > 1e-15:
            # OLS directly on design for delta_sh
            sub_p, cols_a, _, st_p = _prune_covariate_matrix(w, nonbatch_covs, batch_col, min_batch_donors, min_residual_df, ref_level, "age_per_decade", True)
            if st_p != "evaluated" and batch_col is not None and fallback_to_nobatch:
                sub_p, cols_a, _, st_p = _prune_covariate_matrix(w, nonbatch_covs, None, 1, 5, None, "age_per_decade", False)
            if st_p == "evaluated":
                X = sm.add_constant(sub_p[["age_per_decade", *cols_a]].astype(float), has_constant="add")
                fit_d = sm.OLS(sub_p["_delta"].astype(float), X).fit()
                se_sl_sq.append(float(fit_d.bse["age_per_decade"]) ** 2)
                # Contrast SE
                w_oy = sub_p[sub_p["bin"].isin(["YOUNG", "OLD"])].copy()
                if w_oy["bin"].eq("YOUNG").sum() >= 2 and w_oy["bin"].eq("OLD").sum() >= 2:
                    w_oy["_is_old"] = (w_oy["bin"] == "OLD").astype(float)
                    sp_o, ca_o, _, st_o = _prune_covariate_matrix(w_oy, nonbatch_covs, None, 1, 5, None, "_is_old", False)
                    if st_o == "evaluated":
                        X_o = sm.add_constant(sp_o[["_is_old", *ca_o]].astype(float), has_constant="add")
                        fit_do = sm.OLS(sp_o["_delta"].astype(float), X_o).fit()
                        se_oy_sq.append(float(fit_do.bse["_is_old"]) ** 2)
        else:
            se_sl_sq.append(0.0)
            se_oy_sq.append(0.0)

        if apply_quenouille:
            b_s = 0.5 * (w[h1c].to_numpy(float) + w[h2c].to_numpy(float)) - w[obs_col].to_numpy(float)
            bias_mats.append(pd.Series(b_s, index=w.index))

    rms_delta = float(np.sqrt(np.mean(rms_sq))) if rms_sq else np.nan
    se_samp_sl = float(factor * np.sqrt(np.mean(se_sl_sq))) if se_sl_sq else np.nan
    se_samp_oy = float(factor * np.sqrt(np.mean(se_oy_sq))) if se_oy_sq else np.nan
    rho_t2 = float((var_scale * (rms_delta ** 2)) / mse_res) if (np.isfinite(rms_delta) and mse_res > 1e-18) else np.nan
    snr_sl = float(abs(obs_sl) / se_samp_sl) if (np.isfinite(se_samp_sl) and se_samp_sl > 1e-18) else np.nan
    snr_oy = float(abs(oy_obs["diff"]) / se_samp_oy) if (np.isfinite(oy_obs["diff"]) and np.isfinite(se_samp_oy) and se_samp_oy > 1e-18) else np.nan

    q_sl, q_frac, q_csl, q_cp = np.nan, np.nan, np.nan, np.nan
    if apply_quenouille and bias_mats:
        b_df = pd.concat(bias_mats, axis=1).mean(axis=1)
        wb = sub_base.loc[b_df.index].copy()
        wb["_bias"] = b_df.to_numpy(float)
        wb["_corr"] = wb[obs_col].to_numpy(float) - wb["_bias"].to_numpy(float)
        r_b = fit_ols_model_row(wb, "_bias", "age_per_decade", nonbatch_covs, batch_col, min_batch_donors, min_residual_df, ref_level, {}, fallback_to_nobatch)
        r_c = fit_ols_model_row(wb, "_corr", "age_per_decade", nonbatch_covs, batch_col, min_batch_donors, min_residual_df, ref_level, {}, fallback_to_nobatch)
        if r_b is not None:
            q_sl = float(r_b["slope_raw_per_decade"])
            q_frac = float(abs(q_sl) / abs(obs_sl)) if abs(obs_sl) > 1e-18 else np.nan
        if r_c is not None:
            q_csl = float(r_c["slope_raw_per_decade"])
            q_cp = float(r_c["pvalue"])

    return {
        "sh_rms_delta_donor": rms_delta,
        "sh_se_samp_slope": se_samp_sl,
        "sh_se_samp_oy": se_samp_oy,
        "sh_rho_tech2": rho_t2,
        "sh_snr_slope": snr_sl,
        "sh_snr_oy": snr_oy,
        "quenouille_bias_slope": q_sl,
        "quenouille_bias_frac": q_frac,
        "quenouille_corr_slope": q_csl,
        "quenouille_corr_p": q_cp,
    }


def _build_part2_table_and_fit(meta: dict[str, Any], part2_ckpt: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Build Part 2 donor aggregates, fit Full-Depth Primary + 1500/3000/Spline/Strata models + Split-Half/Quenouille."""
    ds_name = meta["dataset"]
    is_gao = ds_name == "gao"
    donors_df = meta["donors_df"]
    ct_names = meta["cell_types"]
    pooled_name = "All_LV_Nuclei" if is_gao else "All_Cells"
    n_donors = len(meta["donor_ids"])
    n_cts = len(ct_names)
    acc = part2_ckpt["acc"]
    not_eval_set = set(meta.get("not_evaluable_all_umi_csv", []))

    # Solve per-CT cell-level fixed-knot spline(df=4) + 1/umi beta via exact normal equations
    cnt_ct = acc["cnt_ct"].reshape(n_donors, n_cts)
    cnt_pl = acc["cnt_pl"]
    lu_ct = acc["lu_ct"].reshape(n_donors, n_cts)
    lu_pl = acc["lu_pl"]
    invu_ct = acc["invu_ct"].reshape(n_donors, n_cts)
    invu_pl = acc["invu_pl"]
    full_ct = acc["full_ct"].reshape(n_donors, n_cts, 3)
    full_pl = acc["full_pl"]
    zsum_ct = acc["zsum_ct"].reshape(n_donors, n_cts, 5)
    zsum_pl = acc["zsum_pl"]

    spline_beta_ct = np.zeros((n_cts, 3, 5), dtype=np.float64)
    for ci in range(n_cts):
        ztz = acc["ztz_c"][ci]
        if ztz[0, 0] >= 20:
            for m_i in range(3):
                spline_beta_ct[ci, m_i] = np.linalg.pinv(ztz) @ acc["zty_c"][ci, m_i]

    ztz_all = acc["ztz_c"].sum(axis=0)
    spline_beta_pool = np.zeros((3, 5), dtype=np.float64)
    if ztz_all[0, 0] >= 20:
        for m_i in range(3):
            spline_beta_pool[m_i] = np.linalg.pinv(ztz_all) @ acc["zty_c"][:, m_i, :].sum(axis=0)

    # Assemble donor-by-group long table for Part 2
    d_rows: list[dict[str, Any]] = []
    for gi, grp in enumerate([*ct_names, pooled_name]):
        is_pool = grp == pooled_name
        c_vec = cnt_pl if is_pool else cnt_ct[:, gi]
        c3k_vec = acc["cnt3k_pl"] if is_pool else acc["cnt3k_ct"].reshape(n_donors, n_cts)[:, gi]
        lu_vec = lu_pl if is_pool else lu_ct[:, gi]
        lu3k_vec = acc["lu3k_pl"] if is_pool else acc["lu3k_ct"].reshape(n_donors, n_cts)[:, gi]
        inv_vec = invu_pl if is_pool else invu_ct[:, gi]
        f_mat = full_pl if is_pool else full_ct[:, gi, :]
        z_mat = zsum_pl if is_pool else zsum_ct[:, gi, :]
        r15_mat = acc["r15_pl"] if is_pool else acc["r15_ct"].reshape(len(ENTROPY_SEEDS), n_donors, n_cts, 3)[:, :, gi, :]
        r3k_mat = acc["r3k_pl"] if is_pool else acc["r3k_ct"].reshape(len(ENTROPY_SEEDS), n_donors, n_cts, 3)[:, :, gi, :]
        r15_on3k = acc["r15_on_3k_pl"] if is_pool else acc["r15_on_3k_ct"].reshape(n_donors, n_cts, 3)[:, gi, :]
        h1_mat = acc["sh_h1_pl"] if is_pool else acc["sh_h1_ct"].reshape(len(SH_SEEDS), n_donors, n_cts, 3)[:, :, gi, :]
        h2_mat = acc["sh_h2_pl"] if is_pool else acc["sh_h2_ct"].reshape(len(SH_SEEDS), n_donors, n_cts, 3)[:, :, gi, :]

        for di, did in enumerate(meta["donor_ids"]):
            n_c = int(c_vec[di])
            if n_c == 0:
                continue
            dr = donors_df.iloc[di]
            mlu = float(lu_vec[di] / n_c)
            rec: dict[str, Any] = {
                "dataset": ds_name,
                "donor_id": str(did),
                "age_years": float(dr["age_years"]),
                "age_per_decade": float(dr["age_per_decade"]),
                "bin": str(dr["bin"]),
                "sex": str(dr["sex"]),
                "sex_male": int(dr["sex_male"]),
                "study_or_pool": str(dr["study"] if is_gao else int(dr["pool_number"])),
                "technology": str(dr["technology"]) if is_gao else "10x 3' v2",
                "group": grp,
                "is_pooled": is_pool,
                "yazar_pool_fe_all_evaluable": True if is_gao else (grp not in not_eval_set),
                "n_cells": n_c,
                "n_cells_3000": int(c3k_vec[di]),
                "_depth": mlu,
                "_depth2": mlu ** 2,
                "_inv_umi": float(inv_vec[di] / n_c),
                "_depth_3k": float(lu3k_vec[di] / c3k_vec[di]) if c3k_vec[di] > 0 else np.nan,
            }
            for m_i, m_nm in enumerate(("shannon", "deman_kl_methods", "deman_kl_sqrt")):
                obs_v = float(f_mat[di, m_i] / n_c)
                rec[f"full_{m_nm}"] = obs_v
                beta_sp = spline_beta_pool[m_i] if is_pool else spline_beta_ct[gi, m_i]
                z_mean = z_mat[di] / float(n_c)
                rec[f"spline_resid_{m_nm}"] = obs_v - float(z_mean @ beta_sp)
                rec[f"r15_s0_{m_nm}"] = float(r15_mat[0, di, m_i] / n_c)
                rec[f"r15_on3k_s0_{m_nm}"] = float(r15_on3k[di, m_i] / c3k_vec[di]) if c3k_vec[di] > 0 else np.nan
                rec[f"r3k_s0_{m_nm}"] = float(r3k_mat[0, di, m_i] / c3k_vec[di]) if c3k_vec[di] > 0 else np.nan
                for si_sh in range(len(SH_SEEDS)):
                    rec[f"h1_s{si_sh}_{m_nm}"] = float(h1_mat[si_sh, di, m_i] / n_c)
                    rec[f"h2_s{si_sh}_{m_nm}"] = float(h2_mat[si_sh, di, m_i] / n_c)
                n_req = 5 if m_i == 0 else 20
                for si_k in range(n_req):
                    rec[f"r15_s{si_k}_{m_nm}"] = float(r15_mat[si_k, di, m_i] / n_c)
                    rec[f"r3k_s{si_k}_{m_nm}"] = float(r3k_mat[si_k, di, m_i] / c3k_vec[di]) if c3k_vec[di] > 0 else np.nan
            d_rows.append(rec)

    df_d = pd.DataFrame(d_rows)
    model_rows: list[dict[str, Any]] = []
    rob_rows: list[dict[str, Any]] = []

    depth_modes = [
        ("full_depth", "primary_full_depth", "n_cells", "_depth", "full_"),
        ("spline_resid_full", "sensitivity_spline_resid", "n_cells", "_depth", "spline_resid_"),
        ("rarefied_1500", "sensitivity_rarefied_1500", "n_cells", "_depth", "r15_s0_"),
        ("rarefied_1500_on_ge3000", "sensitivity_selection_check", "n_cells_3000", "_depth_3k", "r15_on3k_s0_"),
        ("rarefied_3000", "sensitivity_rarefied_3000", "n_cells_3000", "_depth_3k", "r3k_s0_"),
    ]
    bcol = "study_or_pool"
    ref_lv = "ENCODE v4 (Snyder)" if is_gao else None
    min_df = 5 if is_gao else 10

    for d_mode, a_role, c_col, dep_col, pref in depth_modes:
        for m_nm in ("shannon", "deman_kl_methods", "deman_kl_sqrt"):
            out_c = f"{pref}{m_nm}"
            elig_cts = []
            scale_map_full: dict[str, float] = {}
            for grp in [*ct_names, pooled_name]:
                sub = df_d[(df_d["group"] == grp) & (df_d[c_col] >= MIN_CELLS_PER_DONOR) & df_d[out_c].notna()].copy()
                if len(sub) < MIN_DONORS_SLOPE:
                    continue
                is_pool = grp == pooled_name
                sub["n_cells"] = sub[c_col]
                sub["_depth"] = sub[dep_col]
                fb_status = "not_evaluable_within_pool_df" if (not is_gao and not bool(sub["yazar_pool_fe_all_evaluable"].iloc[0])) else "not_evaluable_within_pool_df_slice"

                sh_stats = {}
                if d_mode == "full_depth":
                    h1_cs = [f"h1_s{i}_{m_nm}" for i in range(len(SH_SEEDS))]
                    h2_cs = [f"h2_s{i}_{m_nm}" for i in range(len(SH_SEEDS))]
                    sh_stats = compute_sh_and_quenouille_stats(sub, out_c, h1_cs, h2_cs, ["sex_male", "_depth"], bcol, 2, min_df, ref_lv, False, True, (not is_gao))

                for spec_name, sp_df, nbc, sp_bc, mbd, mdf, rlv in _model_specs_for_dataset(sub, ds_name):
                    row = fit_ols_model_row(
                        sp_df, out_c, "age_per_decade", nbc, sp_bc, mbd, mdf, rlv,
                        meta_fields={
                            "dataset": ds_name, "part": "2", "metric": m_nm, "gene_set": "none",
                            "depth_mode": d_mode, "analysis_role": a_role, "model_spec": spec_name,
                            "group": grp, "is_pooled": is_pool, "is_stratified_global": False,
                            "default_fallback_status": fb_status,
                        },
                        fallback_to_nobatch=(not is_gao),
                    )
                    if row is None:
                        continue
                    oy = fit_contrast_adj_diff(sp_df, out_c, "YOUNG", "OLD", nbc, sp_bc, mbd, mdf, rlv)
                    row.update({"old_minus_young_diff_adj": oy["diff"], "old_minus_young_se_adj": oy["se"], "old_minus_young_p_adj": oy["p"]})
                    if d_mode == "full_depth" and spec_name == "primary":
                        row.update(sh_stats)
                    model_rows.append(row)
                    if spec_name == "primary" and not is_pool and row["yazar_pool_fe_status"] == "evaluated":
                        elig_cts.append(grp)
                        scale_map_full[grp] = float(sub[out_c].std(ddof=1))

            # Stratified_CellType_Global WITH batch FE + depth_x_ct sensitivity
            if len(elig_cts) >= 2:
                st_sub = df_d[(df_d["group"].isin(elig_cts)) & (df_d[c_col] >= MIN_CELLS_PER_DONOR) & df_d[out_c].notna()].copy()
                st_sub["n_cells"] = st_sub[c_col]
                st_sub["_depth"] = st_sub[dep_col]
                for sp_nm, dx in (("stratified_global_primary", False), ("stratified_global_depth_x_ct", True)):
                    srow, _ = fit_stratified_global_with_batch(
                        st_sub, elig_cts, out_c, c_col, "_depth", ds_name, "2", m_nm, d_mode, a_role,
                        model_spec=sp_nm, depth_x_ct=dx,
                    )
                    if srow is not None:
                        if d_mode == "full_depth" and sp_nm == "stratified_global_primary":
                            # Stratified split-half Sampling SE + Quenouille Bias scaled by EXACT SD(M_full)!
                            se_sh_sq, se_oy_sq, rms_v, b_list = [], [], [], []
                            for si_sh in range(len(SH_SEEDS)):
                                h1c, h2c = f"h1_s{si_sh}_{m_nm}", f"h2_s{si_sh}_{m_nm}"
                                st_sub["_d_sh"] = st_sub[h1c] - st_sub[h2c]
                                st_sub["_b_sh"] = 0.5 * (st_sub[h1c] + st_sub[h2c]) - st_sub[out_c]
                                b_list.append(st_sub["_b_sh"])
                                rms_v.append(float(np.mean((st_sub["_d_sh"] / st_sub["group"].map(scale_map_full)) ** 2)))
                                sr_sh, _ = fit_stratified_global_with_batch(
                                    st_sub, elig_cts, "_d_sh", c_col, "_depth", ds_name, "2", m_nm, d_mode, "split_half_se",
                                    custom_scale_map=scale_map_full,
                                )
                                if sr_sh is not None:
                                    se_sh_sq.append(float(sr_sh["se_raw"]) ** 2)
                                    if np.isfinite(sr_sh["old_minus_young_se_adj"]):
                                        se_oy_sq.append(float(sr_sh["old_minus_young_se_adj"]) ** 2)
                            se_samp = float(0.5 * np.sqrt(np.mean(se_sh_sq))) if se_sh_sq else np.nan
                            se_samp_oy = float(0.5 * np.sqrt(np.mean(se_oy_sq))) if se_oy_sq else np.nan
                            st_sub["_b_mean"] = pd.concat(b_list, axis=1).mean(axis=1)
                            st_sub["_q_corr"] = st_sub[out_c] - st_sub["_b_mean"]
                            sr_b, _ = fit_stratified_global_with_batch(st_sub, elig_cts, "_b_mean", c_col, "_depth", ds_name, "2", m_nm, d_mode, "quenouille_bias", custom_scale_map=scale_map_full)
                            sr_c, _ = fit_stratified_global_with_batch(st_sub, elig_cts, "_q_corr", c_col, "_depth", ds_name, "2", m_nm, d_mode, "quenouille_corr", custom_scale_map=scale_map_full)
                            rms_sd = float(np.sqrt(np.mean(rms_v))) if rms_v else np.nan
                            srow.update({
                                "sh_rms_delta_donor": rms_sd,
                                "sh_se_samp_slope": se_samp,
                                "sh_se_samp_oy": se_samp_oy,
                                "sh_rho_tech2": float(0.25 * (rms_sd ** 2) / (srow["resid_sd"] ** 2)) if np.isfinite(rms_sd) else np.nan,
                                "sh_snr_slope": float(abs(srow["slope_raw_per_decade"]) / se_samp) if np.isfinite(se_samp) and se_samp > 0 else np.nan,
                                "sh_snr_oy": float(abs(srow["old_minus_young_diff_adj"]) / se_samp_oy) if np.isfinite(se_samp_oy) and se_samp_oy > 0 and np.isfinite(srow["old_minus_young_diff_adj"]) else np.nan,
                                "quenouille_bias_slope": float(sr_b["slope_raw_per_decade"]) if sr_b else np.nan,
                                "quenouille_bias_frac": float(abs(sr_b["slope_raw_per_decade"]) / abs(srow["slope_raw_per_decade"])) if (sr_b and abs(srow["slope_raw_per_decade"]) > 1e-18) else np.nan,
                                "quenouille_corr_slope": float(sr_c["slope_raw_per_decade"]) if sr_c else np.nan,
                                "quenouille_corr_p": float(sr_c["pvalue"]) if sr_c else np.nan,
                            })
                        model_rows.append(srow)

    # Multi-seed robustness sweeps for 1500 & 3000 rarefaction sensitivities (5 Sh, 20 KL seeds)
    for d_mode, c_col, dep_col, pref in (("rarefied_1500", "n_cells", "_depth", "r15"), ("rarefied_3000", "n_cells_3000", "_depth_3k", "r3k")):
        for m_nm, n_req, min_ok in (("shannon", 5, 4), ("deman_kl_methods", 20, 18), ("deman_kl_sqrt", 20, 18)):
            for grp in [*ct_names, pooled_name]:
                sub = df_d[(df_d["group"] == grp) & (df_d[c_col] >= MIN_CELLS_PER_DONOR)].copy()
                if len(sub) < MIN_DONORS_SLOPE:
                    continue
                is_pool = grp == pooled_name
                sub["n_cells"] = sub[c_col]
                sub["_depth"] = sub[dep_col]
                slopes, pvals = [], []
                for si_k in range(n_req):
                    col_k = f"{pref}_s{si_k}_{m_nm}"
                    r_k = fit_ols_model_row(sub, col_k, "age_per_decade", ["sex_male", "_depth"], bcol, 2, min_df, ref_lv, {}, fallback_to_nobatch=(not is_gao))
                    if r_k is None:
                        break
                    slopes.append(float(r_k["slope_raw_per_decade"]))
                    pvals.append(float(r_k["pvalue"]))
                if len(slopes) != n_req:
                    continue
                pri_s = slopes[0]
                sig = 1 if pri_s >= 0 else -1
                n_same = int(sum((1 if s >= 0 else -1) == sig for s in slopes))
                rob_rows.append({
                    "dataset": ds_name, "part": "2", "depth_mode": d_mode, "metric": m_nm, "gene_set": "none",
                    "group": grp, "is_pooled": is_pool, "n_seeds": n_req, "n_same_sign": n_same,
                    "frac_same_sign": n_same / n_req, "n_p_lt_05": int(sum(p < 0.05 for p in pvals)),
                    "min_slope_raw": float(min(slopes)), "median_slope_raw": float(np.median(slopes)), "max_slope_raw": float(max(slopes)),
                    "median_pvalue": float(np.median(pvals)), "passes_robustness": bool(n_same >= min_ok),
                })

    # Add within-stratum Part 2 models ([1500,2500), [2500,4500), [4500,9000), [9000,18000))
    for si, (s_lbl, _, _) in enumerate(DEPTH_STRATA):
        c_st = acc["st_cnt_ct"][si].reshape(n_donors, n_cts)
        lu_st = acc["st_lu_ct"][si].reshape(n_donors, n_cts)
        s_mat = acc["st_sum_ct"][si].reshape(n_donors, n_cts, 3)
        for m_i, m_nm in enumerate(("shannon", "deman_kl_methods", "deman_kl_sqrt")):
            st_rows_ct = []
            for ci, ct in enumerate(ct_names):
                elig_d = np.flatnonzero(c_st[:, ci] >= MIN_CELLS_PER_DONOR)
                if len(elig_d) < MIN_DONORS_SLOPE:
                    continue
                sdf = donors_df.iloc[elig_d].copy()
                sdf["n_cells"] = c_st[elig_d, ci]
                sdf["_depth"] = lu_st[elig_d, ci] / c_st[elig_d, ci]
                sdf["_y"] = s_mat[elig_d, ci, m_i] / c_st[elig_d, ci]
                sdf["study_or_pool"] = sdf["study"] if is_gao else sdf["pool_number"].astype(str)
                sdf["group"] = ct
                row = fit_ols_model_row(
                    sdf, "_y", "age_per_decade", ["sex_male", "_depth"], bcol, 2, min_df, ref_lv,
                    meta_fields={
                        "dataset": ds_name, "part": "2", "metric": m_nm, "gene_set": "none",
                        "depth_mode": s_lbl, "analysis_role": "sensitivity_within_stratum",
                        "model_spec": "primary", "group": ct, "is_pooled": False,
                        "is_stratified_global": False, "default_fallback_status": "not_evaluable_within_pool_df_slice",
                    },
                    fallback_to_nobatch=(not is_gao),
                )
                if row is not None:
                    model_rows.append(row)
                    if row["yazar_pool_fe_status"] == "evaluated":
                        st_rows_ct.append(sdf)
            if len(st_rows_ct) >= 2:
                st_cat = pd.concat(st_rows_ct, ignore_index=True)
                srow, _ = fit_stratified_global_with_batch(
                    st_cat, sorted(st_cat["group"].unique()), "_y", "n_cells", "_depth", ds_name, "2", m_nm, s_lbl, "sensitivity_within_stratum",
                )
                if srow is not None:
                    model_rows.append(srow)

    mod_df = apply_slice_bh(pd.DataFrame(model_rows), ["dataset", "part", "metric", "gene_set", "depth_mode", "model_spec"])
    rob_df = pd.DataFrame(rob_rows)
    return mod_df, rob_df, {"part2_rows": len(mod_df)}


def _part1_shard_worker(shard_id: int) -> dict[str, Any]:
    """Process 1 contiguous donor shard from disk for Part 1 (center_shift + scatter + 5-SH + 1500 + EqualRef + 20cell)."""
    ds_name = _W_CTX["ds_name"]
    ds_code = 1 if ds_name == "gao" else 2
    n_genes = _W_CTX["n_genes"]
    n_cts = _W_CTX["n_cts"]
    n_donors = _W_CTX["n_donors"]
    bin_arr = _W_CTX["bin_arr"]
    sqrt_ref_full = _W_CTX["sqrt_ref_full"]          # {ct_i: float64[n_genes]}
    sqrt_ref_loo = _W_CTX["sqrt_ref_loo"]            # {(d, ct_i): float64[n_genes]}
    sqrt_ref_15_full = _W_CTX["sqrt_ref_15_full"]
    sqrt_ref_15_loo = _W_CTX["sqrt_ref_15_loo"]
    q_loo_mats = _W_CTX["q_loo_mats"]                # {ct_i: (n_genes, n_Y_T), young_d_list}

    sdir = _shard_dir(ds_name)
    c_path = sdir / f"shard_{shard_id}_cells.bin"
    if not c_path.exists() or c_path.stat().st_size == 0:
        return {"rows": [], "pooled_rows": [], "ct_ztz": np.zeros((n_cts, 5, 5)), "ct_zty": np.zeros((n_cts, 2, 5))}

    cells = np.fromfile(c_path, dtype=np.int64).reshape(-1, 4)
    indices = np.fromfile(sdir / f"shard_{shard_id}_idx.bin", dtype=np.uint16).astype(np.int32)
    data_i32 = np.fromfile(sdir / f"shard_{shard_id}_dat.bin", dtype=np.int32)
    indptr = np.concatenate([[0], np.cumsum(cells[:, 3], dtype=np.int64)])
    csr_i32 = sp.csr_matrix((data_i32, indices, indptr), shape=(len(cells), n_genes))

    d_col, c_col, u_col = cells[:, 0], cells[:, 1], cells[:, 2]
    row_rows: list[dict[str, Any]] = []
    ct_ztz = np.zeros((n_cts, 5, 5), dtype=np.float64)
    ct_zty = np.zeros((n_cts, 2, 5), dtype=np.float64)

    def _process_subset(sub_rows: np.ndarray, di: int, ci: int, is_pool: bool) -> dict[str, Any] | None:
        nc = len(sub_rows)
        if nc < MIN_CELLS_PER_DONOR:
            return None
        sub_i32 = csr_i32[sub_rows].tocsr()
        sub_m = sub_i32.astype(np.float64)
        sub_i = sub_i32.astype(np.int64)
        sub_u = u_col[sub_rows]
        z_sp = spline_design_matrix(sub_u)

        # Normalized L1 probabilities per cell
        inv_u = 1.0 / np.maximum(sub_u.astype(np.float64), 1.0)
        l1_m = sub_m.multiply(inv_u[:, None]).tocsr()
        cent_sum = np.asarray(l1_m.sum(axis=0)).ravel()

        # Vectorized internal LOO scatter Hellinger: p_loo = (cent_sum - p_i) / (nc - 1)
        r_ids = np.repeat(np.arange(nc, dtype=np.int32), np.diff(l1_m.indptr))
        p_val = l1_m.data
        p_loo = np.maximum(cent_sum[l1_m.indices] - p_val, 0.0) / float(nc - 1)
        bc_sc = np.bincount(r_ids, weights=np.sqrt(p_val * p_loo), minlength=nc)
        sc_full = np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc_sc, 0.0, 1.0)))

        # Center shift to YOUNG LOO ref
        cs_full = np.full(nc, np.nan, dtype=np.float64)
        ref_sq = sqrt_ref_loo.get((di, ci)) if bin_arr[di] == 1 else sqrt_ref_full.get(ci)
        pyth_h2_cent, equal_ref_mean = np.nan, np.nan
        if ref_sq is not None:
            cs_full = fast_hellinger_csr_vec(l1_m, ref_sq, already_l1=True)
            cent_mean = cent_sum / float(nc)
            pyth_h2_cent = float(max(0.0, 1.0 - float(np.sum(np.sqrt(cent_mean) * ref_sq))))
        if (not is_pool) and ci in q_loo_mats:
            q_mat, y_d_list = q_loo_mats[ci]
            p_sq_m = l1_m.copy()
            p_sq_m.data = np.sqrt(p_sq_m.data)
            bc_eq = p_sq_m @ q_mat
            d_eq = np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc_eq, 0.0, 1.0)))
            if di in y_d_list:
                j_self = y_d_list.index(di)
                d_eq = np.delete(d_eq, j_self, axis=1)
            equal_ref_mean = float(d_eq.mean()) if d_eq.shape[1] > 0 else np.nan

        if not is_pool:
            ct_ztz[ci] += z_sp.T @ z_sp
            ct_zty[ci, 1] += z_sp.T @ sc_full
            if ref_sq is not None:
                ct_zty[ci, 0] += z_sp.T @ cs_full

        # 1500-UMI rarefaction sensitivity (seed 0)
        rng_15 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, 15001, di, ci if not is_pool else 999]))
        ds15 = multinomial_downsample(sub_i, 1500, rng_15)
        ds15.eliminate_zeros()  # downsample keeps dropped genes as explicit stored zeros; drop them before any stored-entry counting
        l1_15 = ds15.multiply(1.0 / 1500.0).tocsr()
        cent15 = np.asarray(l1_15.sum(axis=0)).ravel()
        r15_ids = np.repeat(np.arange(nc, dtype=np.int32), np.diff(l1_15.indptr))
        pl15 = np.maximum(cent15[l1_15.indices] - l1_15.data, 0.0) / float(nc - 1)
        bc15_sc = np.bincount(r15_ids, weights=np.sqrt(l1_15.data * pl15), minlength=nc)
        sc_15 = np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc15_sc, 0.0, 1.0)))
        ref15_sq = sqrt_ref_15_loo.get((di, ci)) if bin_arr[di] == 1 else sqrt_ref_15_full.get(ci)
        cs_15_mean = float(fast_hellinger_csr_vec(l1_15, ref15_sq, already_l1=True).mean()) if ref15_sq is not None else np.nan

        # 5-seed Binomial Split-Half (against exact full YOUNG ref for center_shift, internal H1/H2 centroids for scatter!)
        sh_cs_h1, sh_cs_h2, sh_sc_h1, sh_sc_h2 = [], [], [], []
        d_raw = np.rint(sub_m.data).astype(np.int64)
        r_m = np.repeat(np.arange(nc, dtype=np.int32), np.diff(sub_m.indptr))
        c_m = sub_m.indices
        for s_sh in SH_SEEDS:
            rng_sh = np.random.default_rng(np.random.SeedSequence([int(s_sh), ds_code, di, ci if not is_pool else 999]))
            h1_d = rng_sh.binomial(d_raw, 0.5)
            h2_d = d_raw - h1_d
            res_halves = []
            for h_arr in (h1_d, h2_d):
                pos = h_arr > 0
                rp, cp, vp = r_m[pos], c_m[pos], h_arr[pos].astype(np.float64)
                tot_h = np.maximum(np.bincount(rp, weights=vp, minlength=nc), 1.0)
                pv = vp / tot_h[rp]
                ch_m = sp.csr_matrix((pv, (rp, cp)), shape=(nc, n_genes))
                cs_v = float(fast_hellinger_csr_vec(ch_m, ref_sq, already_l1=True).mean()) if ref_sq is not None else np.nan
                csum_h = np.asarray(ch_m.sum(axis=0)).ravel()
                pl_h = np.maximum(csum_h[cp] - pv, 0.0) / float(nc - 1)
                bc_h = np.bincount(rp, weights=np.sqrt(pv * pl_h), minlength=nc)
                sc_v = float(np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc_h, 0.0, 1.0))).mean())
                res_halves.append((cs_v, sc_v))
            sh_cs_h1.append(res_halves[0][0])
            sh_cs_h2.append(res_halves[1][0])
            sh_sc_h1.append(res_halves[0][1])
            sh_sc_h2.append(res_halves[1][1])

        # 20-seed subsampled_20cells scatter (exact n_cells=20 control using contiguous l1_m)
        if nc == 20:
            sc_20_mean = float(sc_full.mean())
        else:
            s20_vals = []
            for s20 in SUB20_SEEDS:
                rng20 = np.random.default_rng(np.random.SeedSequence([int(s20), di, ci if not is_pool else 999]))
                pick = rng20.choice(nc, size=20, replace=False)
                sub20 = l1_m[pick]
                c20 = np.asarray(sub20.sum(axis=0)).ravel()
                r20 = np.repeat(np.arange(20, dtype=np.int32), np.diff(sub20.indptr))
                pl20 = np.maximum(c20[sub20.indices] - sub20.data, 0.0) / 19.0
                bc20 = np.bincount(r20, weights=np.sqrt(sub20.data * pl20), minlength=20)
                s20_vals.append(float(np.sqrt(np.maximum(0.0, 1.0 - np.clip(bc20, 0.0, 1.0))).mean()))
            sc_20_mean = float(np.mean(s20_vals))

        rec = {
            "donor_idx": di, "ct_idx": ci, "is_pooled": is_pool, "n_cells": nc,
            "mean_log_umi": float(np.log(sub_u.astype(np.float64)).mean()),
            "mean_inv_umi": float(inv_u.mean()),
            "z_mean": z_sp.mean(axis=0),
            "center_shift": float(cs_full.mean()) if np.all(np.isfinite(cs_full)) else np.nan,
            "center_shift_1500": cs_15_mean,
            "center_shift_eqref": equal_ref_mean,
            "pyth_h2_centroid": pyth_h2_cent,
            "scatter": float(sc_full.mean()),
            "pyth_h2_scatter": float((sc_full ** 2).mean()),
            "scatter_1500": float(sc_15.mean()),
            "scatter_20cells": sc_20_mean,
        }
        for si_sh in range(len(SH_SEEDS)):
            rec[f"h1_s{si_sh}_center_shift"] = sh_cs_h1[si_sh]
            rec[f"h2_s{si_sh}_center_shift"] = sh_cs_h2[si_sh]
            rec[f"h1_s{si_sh}_scatter"] = sh_sc_h1[si_sh]
            rec[f"h2_s{si_sh}_scatter"] = sh_sc_h2[si_sh]
        return rec

    for di in np.unique(d_col):
        dm = d_col == di
        d_rows_idx = np.flatnonzero(dm)
        pr = _process_subset(d_rows_idx, int(di), -1, True)
        if pr is not None:
            row_rows.append(pr)
        for ci in np.unique(c_col[dm]):
            c_rows_idx = d_rows_idx[c_col[d_rows_idx] == ci]
            cr = _process_subset(c_rows_idx, int(di), int(ci), False)
            if cr is not None:
                row_rows.append(cr)

    return {"rows": row_rows, "ct_ztz": ct_ztz, "ct_zty": ct_zty}


def run_pass2_part1_and_export(meta: dict[str, Any], ds_name: str, part2_mod_df: pd.DataFrame) -> dict[str, Any]:
    """Pass 2 (ZERO HDF5 reads): Part 1 center_shift & scatter + EqualRef + Pythagorean + Part A cross-check -> CSVs."""
    t0 = time.perf_counter()
    with open(CACHE_DIR / f"{ds_name}_pass1_part1_refs.pkl", "rb") as pf:
        p1r = pickle.load(pf)

    donors_df = meta["donors_df"]
    ct_names = meta["cell_types"]
    pooled_name = "All_LV_Nuclei" if ds_name == "gao" else "All_Cells"
    is_gao = ds_name == "gao"
    n_donors, n_cts = len(meta["donor_ids"]), len(ct_names)
    n_genes = p1r["y_full_cnt"].shape[1]
    bin_arr = donors_df["bin"].map({"LT20": 0, "YOUNG": 1, "MID": 2, "OLD": 3}).to_numpy(dtype=np.int32)
    not_eval_set = set(meta.get("not_evaluable_all_umi_csv", []))

    cnt_ct_15 = p1r["cnt_ct_15"]
    cnt_pl_15 = p1r["cnt_pl_15"]
    y_pair_map = p1r["y_pair_map"]
    y_full_cnt = p1r["y_full_cnt"]
    y_1500_cnt0 = p1r["y_1500_cnt0"]
    y_n_cells = p1r["y_n_cells"]

    # Build normalized sqrt(YOUNG ref) vectors cleanly (ge20 YOUNG donors, LOO min >= 50 cells, matching Part A!)
    def _build_ref_lookups(cnt_tensor: np.ndarray) -> tuple[dict[int, np.ndarray], dict[tuple[int, int], np.ndarray], dict[int, tuple[np.ndarray, list[int]]]]:
        full_sq: dict[int, np.ndarray] = {}
        loo_sq: dict[tuple[int, int], np.ndarray] = {}
        q_mats: dict[int, tuple[np.ndarray, list[int]]] = {}
        for ci in range(n_cts):
            y_elig = [di for di in range(n_donors) if bin_arr[di] == 1 and cnt_ct_15[di, ci] >= MIN_CELLS_PER_DONOR]
            if len(y_elig) < MIN_YOUNG_REF:
                continue
            tot_cells = sum(int(y_n_cells[y_pair_map[di * n_cts + ci]]) for di in y_elig)
            min_loo = min(tot_cells - int(y_n_cells[y_pair_map[di * n_cts + ci]]) for di in y_elig)
            if min_loo < MIN_LOO_REF_CELLS:
                continue
            d_vecs = {di: cnt_tensor[y_pair_map[di * n_cts + ci]].astype(np.float64) for di in y_elig}
            tot_v = sum(d_vecs.values())
            if tot_v.sum() <= 0:
                continue
            full_sq[ci] = np.sqrt(tot_v / tot_v.sum())
            cols_q = []
            for di in y_elig:
                lv = tot_v - d_vecs[di]
                sq_l = np.sqrt(lv / lv.sum())
                loo_sq[(di, ci)] = sq_l
                cols_q.append(sq_l)
            q_mats[ci] = (np.column_stack(cols_q), y_elig)
        # Pooled YOUNG ref (ci = -1)
        yp_elig = [di for di in range(n_donors) if bin_arr[di] == 1 and cnt_pl_15[di] >= MIN_CELLS_PER_DONOR]
        if len(yp_elig) >= MIN_YOUNG_REF:
            p_vecs = {}
            p_nc = {}
            for di in yp_elig:
                pv = np.zeros(n_genes, dtype=np.float64)
                nc_p = 0
                for ci in range(n_cts):
                    idx_r = y_pair_map[di * n_cts + ci]
                    pv += cnt_tensor[idx_r]
                    nc_p += int(y_n_cells[idx_r])
                p_vecs[di] = pv
                p_nc[di] = nc_p
            if min(sum(p_nc.values()) - p_nc[d] for d in yp_elig) >= MIN_LOO_REF_CELLS:
                tot_p = sum(p_vecs.values())
                full_sq[-1] = np.sqrt(tot_p / tot_p.sum())
                for di in yp_elig:
                    lv = tot_p - p_vecs[di]
                    loo_sq[(di, -1)] = np.sqrt(lv / lv.sum())
        return full_sq, loo_sq, q_mats

    sq_full, sq_loo, q_mats = _build_ref_lookups(y_full_cnt)
    sq15_full, sq15_loo, _ = _build_ref_lookups(y_1500_cnt0)
    del p1r, y_full_cnt, y_1500_cnt0, cnt_ct_15, cnt_pl_15, y_pair_map, y_n_cells
    gc.collect()

    global _W_CTX
    _W_CTX = {
        "ds_name": ds_name, "n_genes": n_genes, "n_cts": n_cts, "n_donors": n_donors,
        "bin_arr": bin_arr, "sqrt_ref_full": sq_full, "sqrt_ref_loo": sq_loo,
        "sqrt_ref_15_full": sq15_full, "sqrt_ref_15_loo": sq15_loo, "q_loo_mats": q_mats,
    }
    all_rows: list[dict[str, Any]] = []
    ct_ztz = np.zeros((n_cts, 5, 5), dtype=np.float64)
    ct_zty = np.zeros((n_cts, 2, 5), dtype=np.float64)
    try:
        with ProcessPoolExecutor(max_workers=SHARD_WORKERS, initializer=_worker_init, initargs=()) as ex:
            for res in ex.map(_part1_shard_worker, list(range(N_SHARDS))):
                all_rows.extend(res["rows"])
                ct_ztz += res["ct_ztz"]
                ct_zty += res["ct_zty"]
    finally:
        _W_CTX = {}
        gc.collect()

    beta_sp = np.zeros((n_cts, 2, 5), dtype=np.float64)
    for ci in range(n_cts):
        if ct_ztz[ci, 0, 0] >= 20:
            for k_m in range(2):
                beta_sp[ci, k_m] = np.linalg.pinv(ct_ztz[ci]) @ ct_zty[ci, k_m]

    # Build donor-level Part 1 dataframe & compute spline residuals
    df1_rows = []
    for r in all_rows:
        di, ci, is_pool = r["donor_idx"], r["ct_idx"], r["is_pooled"]
        grp = pooled_name if is_pool else ct_names[ci]
        dr = donors_df.iloc[di]
        z_m = r["z_mean"]
        rec = {
            "dataset": ds_name,
            "donor_id": str(dr["donor_id"]),
            "age_years": float(dr["age_years"]),
            "age_per_decade": float(dr["age_per_decade"]),
            "bin": str(dr["bin"]),
            "sex": str(dr["sex"]),
            "sex_male": int(dr["sex_male"]),
            "study_or_pool": str(dr["study"] if is_gao else int(dr["pool_number"])),
            "technology": str(dr["technology"]) if is_gao else "10x 3' v2",
            "group": grp,
            "is_pooled": is_pool,
            "yazar_pool_fe_all_evaluable": True if is_gao else (grp not in not_eval_set),
            "n_cells": int(r["n_cells"]),
            "_depth": float(r["mean_log_umi"]),
            "_depth2": float(r["mean_log_umi"]) ** 2,
            "_inv_umi": float(r["mean_inv_umi"]),
            "_log_n_cells": float(np.log(max(r["n_cells"], 1))),
            "center_shift": r["center_shift"],
            "center_shift_1500": r["center_shift_1500"],
            "center_shift_eqref": r["center_shift_eqref"],
            "spline_resid_center_shift": (r["center_shift"] - float(z_m @ beta_sp[ci, 0])) if (not is_pool and np.isfinite(r["center_shift"])) else np.nan,
            "pyth_h2_centroid": r["pyth_h2_centroid"],
            "scatter": r["scatter"],
            "scatter_1500": r["scatter_1500"],
            "scatter_20cells": r["scatter_20cells"],
            "spline_resid_scatter": (r["scatter"] - float(z_m @ beta_sp[ci, 1])) if (not is_pool) else np.nan,
            "pyth_h2_scatter": r["pyth_h2_scatter"],
        }
        for si_sh in range(len(SH_SEEDS)):
            for m_nm in ("center_shift", "scatter"):
                rec[f"h1_s{si_sh}_{m_nm}"] = r[f"h1_s{si_sh}_{m_nm}"]
                rec[f"h2_s{si_sh}_{m_nm}"] = r[f"h2_s{si_sh}_{m_nm}"]
        df1_rows.append(rec)

    df1 = pd.DataFrame(df1_rows)

    # Bit-exact Part A cross-check (against results/followup-between-subject-{ds_name}-donors.csv when present!)
    pa_path = RESULTS_DIR / f"followup-between-subject-{ds_name}-donors.csv"
    cross_chk: dict[str, Any] = {"checked": False}
    if pa_path.exists():
        pa = pd.read_csv(pa_path)
        mrg = df1.merge(pa[["donor_id", "group", "n_cells_ge1500", "a1_dist_full_1500", "a3_disp_full_1500", "eligible_a1_1500", "eligible_a3_1500"]], on=["donor_id", "group"], how="inner")
        assert (mrg["n_cells"] == mrg["n_cells_ge1500"]).all(), "Part 1 n_cells_ge1500 mismatch vs Part A CSV!"
        m_cs = mrg["eligible_a1_1500"].astype(bool) & mrg["center_shift"].notna() & mrg["a1_dist_full_1500"].notna()
        max_cs_diff = float(np.abs(mrg.loc[m_cs, "center_shift"] - mrg.loc[m_cs, "a1_dist_full_1500"]).max()) if m_cs.any() else 0.0
        m_sc = mrg["eligible_a3_1500"].astype(bool) & mrg["scatter"].notna() & mrg["a3_disp_full_1500"].notna()
        max_sc_diff = float(np.abs(mrg.loc[m_sc, "scatter"] - mrg.loc[m_sc, "a3_disp_full_1500"]).max()) if m_sc.any() else 0.0
        assert max_cs_diff <= 1e-10 and max_sc_diff <= 1e-10, f"Part A mismatch: cs={max_cs_diff}, sc={max_sc_diff}"
        cross_chk = {"checked": True, "n_matched_rows": int(len(mrg)), "max_center_shift_diff": max_cs_diff, "max_scatter_diff": max_sc_diff}

    # Fit all Part 1 OLS + Stratified_CellType_Global + BH + Split-Half/Quenouille
    model_rows: list[dict[str, Any]] = []
    bcol, ref_lv = "study_or_pool", ("ENCODE v4 (Snyder)" if is_gao else None)
    min_df = 5 if is_gao else 10

    arms_p1 = [
        ("center_shift", "full_depth", "primary_full_depth", "center_shift", False),
        ("center_shift", "equal_ref_size", "sensitivity_equal_ref", "center_shift_eqref", False),
        ("center_shift", "spline_resid_full", "sensitivity_spline_resid", "spline_resid_center_shift", False),
        ("center_shift", "pythagorean_centroid_h2", "sensitivity_pythagorean", "pyth_h2_centroid", False),
        ("center_shift", "rarefied_1500", "sensitivity_rarefied_1500", "center_shift_1500", False),
        ("scatter", "full_depth", "primary_full_depth", "scatter", True),
        ("scatter", "subsampled_20cells", "sensitivity_20cells", "scatter_20cells", True),
        ("scatter", "spline_resid_full", "sensitivity_spline_resid", "spline_resid_scatter", True),
        ("scatter", "pythagorean_scatter_h2", "sensitivity_pythagorean", "pyth_h2_scatter", True),
        ("scatter", "rarefied_1500", "sensitivity_rarefied_1500", "scatter_1500", True),
    ]
    for m_nm, d_mode, a_role, out_c, add_log_nc in arms_p1:
        elig_cts = []
        scale_full: dict[str, float] = {}
        for grp in [*ct_names, pooled_name]:
            sub = df1[(df1["group"] == grp) & df1[out_c].notna()].copy()
            if len(sub) < MIN_DONORS_SLOPE:
                continue
            is_pool = grp == pooled_name
            fb_status = "not_evaluable_within_pool_df" if (not is_gao and not bool(sub["yazar_pool_fe_all_evaluable"].iloc[0])) else "not_evaluable_within_pool_df_slice"

            sh_stats = {}
            if d_mode == "full_depth":
                h1_cs = [f"h1_s{i}_{m_nm}" for i in range(len(SH_SEEDS))]
                h2_cs = [f"h2_s{i}_{m_nm}" for i in range(len(SH_SEEDS))]
                sh_stats = compute_sh_and_quenouille_stats(sub, out_c, h1_cs, h2_cs, ["sex_male", "_depth"], bcol, 2, min_df, ref_lv, False, True, (not is_gao))

            y_mean = float(sub.loc[sub["bin"] == "YOUNG", out_c].mean())
            o_mean = float(sub.loc[sub["bin"] == "OLD", out_c].mean())
            for spec_name, sp_df, nbc, sp_bc, mbd, mdf, rlv in _model_specs_for_dataset(sub, ds_name, extra_log_ncells=add_log_nc):
                row = fit_ols_model_row(
                    sp_df, out_c, "age_per_decade", nbc, sp_bc, mbd, mdf, rlv,
                    meta_fields={
                        "dataset": ds_name, "part": "1", "metric": m_nm, "gene_set": "none",
                        "depth_mode": d_mode, "analysis_role": a_role, "model_spec": spec_name,
                        "group": grp, "is_pooled": is_pool, "is_stratified_global": False,
                        "default_fallback_status": fb_status,
                    },
                    fallback_to_nobatch=(not is_gao),
                )
                if row is None:
                    continue
                oy = fit_contrast_adj_diff(sp_df, out_c, "YOUNG", "OLD", nbc, sp_bc, mbd, mdf, rlv)
                row.update({"young_loo_mean": y_mean, "old_mean": o_mean, "old_minus_young_diff_adj": oy["diff"], "old_minus_young_se_adj": oy["se"], "old_minus_young_p_adj": oy["p"]})
                if d_mode == "full_depth" and spec_name == "primary":
                    row.update(sh_stats)
                model_rows.append(row)
                if spec_name == "primary" and not is_pool and row["yazar_pool_fe_status"] == "evaluated":
                    elig_cts.append(grp)
                    scale_full[grp] = float(sub[out_c].std(ddof=1))

        if len(elig_cts) >= 2:
            st_sub = df1[(df1["group"].isin(elig_cts)) & df1[out_c].notna()].copy()
            for sp_nm, dx in (("stratified_global_primary", False), ("stratified_global_depth_x_ct", True)):
                srow, _ = fit_stratified_global_with_batch(
                    st_sub, elig_cts, out_c, "n_cells", "_depth", ds_name, "1", m_nm, d_mode, a_role, model_spec=sp_nm, depth_x_ct=dx,
                )
                if srow is not None:
                    y_st = float(st_sub.loc[st_sub["bin"] == "YOUNG", out_c].mean())
                    o_st = float(st_sub.loc[st_sub["bin"] == "OLD", out_c].mean())
                    srow.update({"young_loo_mean": y_st, "old_mean": o_st})
                    if d_mode == "full_depth" and sp_nm == "stratified_global_primary":
                        se_sh_sq, se_oy_sq, rms_v, b_list = [], [], [], []
                        for si_sh in range(len(SH_SEEDS)):
                            h1c, h2c = f"h1_s{si_sh}_{m_nm}", f"h2_s{si_sh}_{m_nm}"
                            st_sub["_d_sh"] = st_sub[h1c] - st_sub[h2c]
                            st_sub["_b_sh"] = 0.5 * (st_sub[h1c] + st_sub[h2c]) - st_sub[out_c]
                            b_list.append(st_sub["_b_sh"])
                            rms_v.append(float(np.mean((st_sub["_d_sh"] / st_sub["group"].map(scale_full)) ** 2)))
                            sr_sh, _ = fit_stratified_global_with_batch(st_sub, elig_cts, "_d_sh", "n_cells", "_depth", ds_name, "1", m_nm, d_mode, "split_half_se", custom_scale_map=scale_full)
                            if sr_sh is not None:
                                se_sh_sq.append(float(sr_sh["se_raw"]) ** 2)
                                if np.isfinite(sr_sh["old_minus_young_se_adj"]):
                                    se_oy_sq.append(float(sr_sh["old_minus_young_se_adj"]) ** 2)
                        se_samp = float(0.5 * np.sqrt(np.mean(se_sh_sq))) if se_sh_sq else np.nan
                        se_samp_oy = float(0.5 * np.sqrt(np.mean(se_oy_sq))) if se_oy_sq else np.nan
                        st_sub["_b_mean"] = pd.concat(b_list, axis=1).mean(axis=1)
                        st_sub["_q_corr"] = st_sub[out_c] - st_sub["_b_mean"]
                        sr_b, _ = fit_stratified_global_with_batch(st_sub, elig_cts, "_b_mean", "n_cells", "_depth", ds_name, "1", m_nm, d_mode, "quenouille_bias", custom_scale_map=scale_full)
                        sr_c, _ = fit_stratified_global_with_batch(st_sub, elig_cts, "_q_corr", "n_cells", "_depth", ds_name, "1", m_nm, d_mode, "quenouille_corr", custom_scale_map=scale_full)
                        rms_sd = float(np.sqrt(np.mean(rms_v))) if rms_v else np.nan
                        srow.update({
                            "sh_rms_delta_donor": rms_sd,
                            "sh_se_samp_slope": se_samp,
                            "sh_se_samp_oy": se_samp_oy,
                            "sh_rho_tech2": float(0.25 * (rms_sd ** 2) / (srow["resid_sd"] ** 2)) if np.isfinite(rms_sd) else np.nan,
                            "sh_snr_slope": float(abs(srow["slope_raw_per_decade"]) / se_samp) if np.isfinite(se_samp) and se_samp > 0 else np.nan,
                            "sh_snr_oy": float(abs(srow["old_minus_young_diff_adj"]) / se_samp_oy) if np.isfinite(se_samp_oy) and se_samp_oy > 0 and np.isfinite(srow["old_minus_young_diff_adj"]) else np.nan,
                            "quenouille_bias_slope": float(sr_b["slope_raw_per_decade"]) if sr_b else np.nan,
                            "quenouille_bias_frac": float(abs(sr_b["slope_raw_per_decade"]) / abs(srow["slope_raw_per_decade"])) if (sr_b and abs(srow["slope_raw_per_decade"]) > 1e-18) else np.nan,
                            "quenouille_corr_slope": float(sr_c["slope_raw_per_decade"]) if sr_c else np.nan,
                            "quenouille_corr_p": float(sr_c["pvalue"]) if sr_c else np.nan,
                        })
                    model_rows.append(srow)

    mod1_df = apply_slice_bh(pd.DataFrame(model_rows), ["dataset", "part", "metric", "gene_set", "depth_mode", "model_spec"])
    write_clean_csv(mod1_df, RESULTS_DIR / f"followup-fulldepth-{ds_name}-cell-drift.csv")
    write_clean_csv(part2_mod_df, RESULTS_DIR / f"followup-fulldepth-{ds_name}-spread.csv")

    elapsed = round(time.perf_counter() - t0, 2)
    _log(ds_name, f"Pass 2 (Part 1 + CSV exports) complete in {elapsed}s")
    return {"cross_check": cross_chk, "elapsed_pass2_s": elapsed, "n_part1_rows": len(mod1_df)}


def build_e1_offprogram_lookup(
    y_full_cnt: np.ndarray,
    y_full_det: np.ndarray,
    y_n_cells: np.ndarray,
    y_pair_map: dict[int, int],
    cnt_ct_15: np.ndarray,
    bin_arr: np.ndarray,
    n_donors: int,
    n_cts: int,
    n_genes: int,
    ge20_only_ref: bool = True,
    exclude_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    """Fixed H3 off-program marker builder checking EACH other CT jointly with invariant LOO other sources."""
    y_donors = [d for d in range(n_donors) if bin_arr[d] == 1]

    def _donor_contrib(d: int, ci: int) -> bool:
        if ge20_only_ref:
            return cnt_ct_15[d, ci] >= MIN_CELLS_PER_DONOR
        return int(y_n_cells[y_pair_map[d * n_cts + ci]]) > 0

    # Tot stats per cell type + invariant source CT eligibility (min_{d in YOUNG} other_loo >= 50)
    tot_cnt = np.zeros((n_cts, n_genes), dtype=np.float64)
    tot_det = np.zeros((n_cts, n_genes), dtype=np.float64)
    tot_lib = np.zeros(n_cts, dtype=np.float64)
    tot_n = np.zeros(n_cts, dtype=np.int64)
    invariant_source_cts: list[int] = []
    for ci in range(n_cts):
        contrib_d = [d for d in y_donors if _donor_contrib(d, ci)]
        nc_list = [int(y_n_cells[y_pair_map[d * n_cts + ci]]) for d in contrib_d]
        for d in contrib_d:
            idx = y_pair_map[d * n_cts + ci]
            tot_cnt[ci] += y_full_cnt[idx]
            tot_det[ci] += y_full_det[idx]
            tot_lib[ci] += float(y_full_cnt[idx].sum())
            tot_n[ci] += int(y_n_cells[idx])
        if len(contrib_d) >= MIN_YOUNG_REF and (tot_n[ci] - max(nc_list)) >= MIN_LOO_REF_CELLS and tot_lib[ci] > 0:
            invariant_source_cts.append(ci)

    def _compute_markers(target_ci: int, excl_d: int | None) -> tuple[np.ndarray, list[str]]:
        t_n = int(tot_n[target_ci])
        t_lib = float(tot_lib[target_ci])
        t_cnt = tot_cnt[target_ci].copy()
        t_det = tot_det[target_ci].copy()
        if excl_d is not None and bin_arr[excl_d] == 1 and _donor_contrib(excl_d, target_ci):
            idx = y_pair_map[excl_d * n_cts + target_ci]
            t_n -= int(y_n_cells[idx])
            t_cnt -= y_full_cnt[idx]
            t_det -= y_full_det[idx]
            t_lib -= float(y_full_cnt[idx].sum())
        if t_n < MIN_LOO_REF_CELLS or t_lib <= 0:
            return np.zeros(0, dtype=np.int32), []

        t_det_frac = t_det / float(t_n)
        t_cp10k = t_cnt / t_lib * 1e4
        t_low = t_det_frac <= 0.05

        qual_any = np.zeros(n_genes, dtype=bool)
        best_cp10k = np.zeros(n_genes, dtype=np.float64)
        best_src = np.full(n_genes, -1, dtype=np.int32)
        for oci in invariant_source_cts:
            if oci == target_ci:
                continue
            o_n = int(tot_n[oci])
            o_lib = float(tot_lib[oci])
            o_cnt = tot_cnt[oci]
            o_det = tot_det[oci]
            if excl_d is not None and bin_arr[excl_d] == 1 and _donor_contrib(excl_d, oci):
                idx_o = y_pair_map[excl_d * n_cts + oci]
                o_n -= int(y_n_cells[idx_o])
                o_cnt = o_cnt - y_full_cnt[idx_o]
                o_det = o_det - y_full_det[idx_o]
                o_lib -= float(y_full_cnt[idx_o].sum())
            if o_n < MIN_LOO_REF_CELLS or o_lib <= 0:
                continue
            o_det_frac = o_det / float(o_n)
            o_cp10k = o_cnt / o_lib * 1e4
            pass_o = t_low & (o_det_frac >= 0.20) & (o_cp10k >= 20.0) & (o_cp10k >= 5.0 * (t_cp10k + 1.0))
            if exclude_mask is not None:
                pass_o &= ~exclude_mask
            qual_any |= pass_o
            upd = pass_o & (o_cp10k > best_cp10k)
            best_cp10k[upd] = o_cp10k[upd]
            best_src[upd] = oci
        g_idx = np.flatnonzero(qual_any).astype(np.int32)
        return g_idx, [str(best_src[g]) for g in g_idx]

    e1_full: dict[int, np.ndarray] = {}
    e1_loo: dict[tuple[int, int], np.ndarray] = {}
    summary: dict[int, int] = {}
    for ci in range(n_cts):
        if ci not in invariant_source_cts:
            continue
        g_f, _ = _compute_markers(ci, None)
        e1_full[ci] = g_f
        summary[ci] = len(g_f)
        for d in y_donors:
            g_l, _ = _compute_markers(ci, d)
            e1_loo[(d, ci)] = g_l
    return {"full": e1_full, "loo": e1_loo, "summary": summary, "invariant_cts": invariant_source_cts}


def build_e2_silent_lookup(
    det_tensor: np.ndarray,
    n_cells_vec: np.ndarray,
    y_pair_map: dict[int, int],
    cnt_ct_elig: np.ndarray,
    bin_arr: np.ndarray,
    n_donors: int,
    n_cts: int,
    n_genes: int,
    min_loo_young_cells: int,
    ge20_only_ref: bool = True,
    exclude_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    """Build < 1% YOUNG-SILENT LOO gene sets (stratum-matched >=200 or global/1500 >=50)."""
    y_donors = [d for d in range(n_donors) if bin_arr[d] == 1]
    e2_full: dict[int, np.ndarray] = {}
    e2_loo: dict[tuple[int, int], np.ndarray] = {}
    summary: dict[int, int] = {}
    for ci in range(n_cts):
        contrib = [d for d in y_donors if (cnt_ct_elig[d, ci] >= MIN_CELLS_PER_DONOR if ge20_only_ref else int(n_cells_vec[y_pair_map[d * n_cts + ci]]) > 0)]
        if len(contrib) < MIN_YOUNG_REF:
            continue
        tot_n = sum(int(n_cells_vec[y_pair_map[d * n_cts + ci]]) for d in contrib)
        min_loo = min(tot_n - int(n_cells_vec[y_pair_map[d * n_cts + ci]]) for d in contrib)
        if min_loo < min_loo_young_cells:
            continue
        tot_d = sum(det_tensor[y_pair_map[d * n_cts + ci]].astype(np.int64) for d in contrib)
        m_f = (tot_d.astype(np.float64) / float(tot_n)) < MAX_E2_GLOBAL_POOL_FRAC
        if exclude_mask is not None:
            m_f &= ~exclude_mask
        e2_full[ci] = np.flatnonzero(m_f).astype(np.int32)
        summary[ci] = len(e2_full[ci])
        for d in y_donors:
            idx = y_pair_map[d * n_cts + ci]
            if d in contrib:
                dn = tot_n - int(n_cells_vec[idx])
                dd = tot_d - det_tensor[idx].astype(np.int64)
            else:
                dn, dd = tot_n, tot_d
            m_l = (dd.astype(np.float64) / float(dn)) < MAX_E2_GLOBAL_POOL_FRAC
            if exclude_mask is not None:
                m_l &= ~exclude_mask
            e2_loo[(d, ci)] = np.flatnonzero(m_l).astype(np.int32)
    return {"full": e2_full, "loo": e2_loo, "summary": summary}


def compute_binary_ectopic_metrics(
    bin_mat: sp.csr_matrix,
    rng_perm: np.random.Generator,
    n_genes_set: int,
    run_excess: bool = True,
) -> dict[str, float]:
    """Row-blocked (block=512) Jaccard (all, union>0, both>0), singleton_fraction (0/0->NaN), evenness/log(|E|), & excess_jaccard."""
    nc = bin_mat.shape[0]
    if nc < 2 or n_genes_set == 0:
        return {"jaccard_all_pairs": np.nan, "jaccard_union_pos": np.nan, "jaccard_both_pos": np.nan, "singleton_fraction": np.nan, "top10_ectopic_share": np.nan, "ectopic_pool_evenness": np.nan, "excess_jaccard": np.nan, "n_pos_cells": 0, "has_any_ectopic": 0}
    b_csr = bin_mat.tocsr(copy=True)
    b_csr.data = np.ones_like(b_csr.data, dtype=np.float64)
    col_cnts = np.asarray(b_csr.sum(axis=0)).ravel().astype(np.int64)
    n_det_genes = int((col_cnts >= 1).sum())
    n_sing = int((col_cnts == 1).sum())
    tot_hits = int(col_cnts.sum())
    sing_frac = float(n_sing / n_det_genes) if n_det_genes > 0 else np.nan
    top10_sh = float(np.sort(col_cnts)[::-1][:10].sum() / tot_hits) if tot_hits > 0 else np.nan
    if tot_hits > 0 and n_genes_set > 1:
        nz = col_cnts[col_cnts > 0].astype(np.float64) / float(tot_hits)
        ev = float(-np.sum(nz * np.log(nz)) / np.log(float(n_genes_set)))
    else:
        ev = np.nan

    deg = np.diff(b_csr.indptr).astype(np.int64)
    pos_rows = np.flatnonzero(deg > 0)
    n_pos = len(pos_rows)
    p_all = nc * (nc - 1) // 2
    n_zero = nc - n_pos
    p_zero_zero = n_zero * (n_zero - 1) // 2
    p_union_pos = p_all - p_zero_zero
    p_both_pos = n_pos * (n_pos - 1) // 2
    if n_pos < 2:
        return {"jaccard_all_pairs": 0.0, "jaccard_union_pos": 0.0 if p_union_pos > 0 else np.nan, "jaccard_both_pos": np.nan, "singleton_fraction": sing_frac, "top10_ectopic_share": top10_sh, "ectopic_pool_evenness": ev, "excess_jaccard": np.nan, "n_pos_cells": n_pos, "has_any_ectopic": int(n_det_genes > 0)}

    b_pos = b_csr[pos_rows].tocsr()
    deg_pos = deg[pos_rows].astype(np.float64)
    j_sum = 0.0
    blk = 512
    for r0 in range(0, n_pos, blk):
        r1 = min(r0 + blk, n_pos)
        inter = (b_pos[r0:r1] @ b_pos.T).toarray()
        # Keep only strictly upper triangle cols > r0 + local_r
        for lr in range(r1 - r0):
            gr = r0 + lr
            if gr + 1 >= n_pos:
                break
            iv = inter[lr, gr + 1 :]
            un = deg_pos[gr] + deg_pos[gr + 1 :] - iv
            j_sum += float(np.sum(iv / np.maximum(un, 1.0)))

    j_both = float(j_sum / p_both_pos)
    j_union = float(j_sum / p_union_pos)
    j_all = float(j_sum / p_all)

    exc = np.nan
    if run_excess:
        # Preserve gene marginal counts k_g under P=20 random permutations
        nz_cols = np.flatnonzero(col_cnts > 0)
        b_sub = b_csr[:, nz_cols].tocsc()
        k_nz = col_cnts[nz_cols]
        null_js = []
        n_perms = 20
        if n_pos <= 500:
            row_pool = np.arange(nc, dtype=np.int32)
            for _ in range(n_perms):
                rows_l, cols_l = [], []
                for cj, kc in enumerate(k_nz):
                    if kc >= nc:
                        chosen = row_pool
                    elif kc == 1:
                        chosen = np.array([int(rng_perm.integers(0, nc))], dtype=np.int32)
                    else:
                        chosen = rng_perm.choice(nc, size=int(kc), replace=False)
                    rows_l.append(chosen)
                    cols_l.append(np.full(len(chosen), cj, dtype=np.int32))
                bn = sp.csr_matrix((np.ones(sum(len(x) for x in rows_l), dtype=np.float64), (np.concatenate(rows_l), np.concatenate(cols_l))), shape=(nc, len(nz_cols)))
                dn = np.diff(bn.indptr)
                pn = np.flatnonzero(dn > 0)
                if len(pn) >= 2:
                    bp = bn[pn]
                    dp = dn[pn].astype(np.float64)
                    im = (bp @ bp.T).toarray()
                    iu = np.triu_indices(len(pn), k=1)
                    iv = im[iu]
                    un = dp[iu[0]] + dp[iu[1]] - iv
                    null_js.append(float(np.mean(iv / np.maximum(un, 1.0))))
        else:
            # Monte Carlo 20,000 active-pair approximation when n_pos > 500
            # Expected both-pos Jaccard under independent Bernoulli placement at per-gene rate k_g / nc:
            # Simulate 20 draws of 20k random cell pairs preserving per-gene marginals concisely
            p_sample = min(20000, p_both_pos)
            sets = [set() for _ in range(nc)]
            # Sample marginal vectors via exact Binomial/choice into donor sparsity, then evaluate 20k random active pairs
            row_pool = np.arange(nc, dtype=np.int32)
            for _ in range(n_perms):
                rows_l, cols_l = [], []
                for cj, kc in enumerate(k_nz):
                    chosen = row_pool if kc >= nc else rng_perm.choice(nc, size=int(kc), replace=False)
                    rows_l.append(chosen)
                    cols_l.append(np.full(len(chosen), cj, dtype=np.int32))
                bn = sp.csr_matrix((np.ones(sum(len(x) for x in rows_l), dtype=np.float64), (np.concatenate(rows_l), np.concatenate(cols_l))), shape=(nc, len(nz_cols)))
                dn = np.diff(bn.indptr)
                pn = np.flatnonzero(dn > 0)
                if len(pn) >= 2:
                    i1 = pn[rng_perm.integers(0, len(pn), size=p_sample)]
                    i2 = pn[rng_perm.integers(0, len(pn), size=p_sample)]
                    neq = i1 != i2
                    if np.any(neq):
                        i1, i2 = i1[neq], i2[neq]
                        inter_mc = np.asarray(bn[i1].multiply(bn[i2]).sum(axis=1)).ravel()
                        un_mc = dn[i1] + dn[i2] - inter_mc
                        null_js.append(float(np.mean(inter_mc / np.maximum(un_mc, 1.0))))
        if null_js:
            exc = j_both - float(np.mean(null_js))

    return {
        "jaccard_all_pairs": j_all,
        "jaccard_union_pos": j_union,
        "jaccard_both_pos": j_both,
        "singleton_fraction": sing_frac,
        "top10_ectopic_share": top10_sh,
        "ectopic_pool_evenness": ev,
        "excess_jaccard": exc,
        "n_pos_cells": n_pos,
        "has_any_ectopic": int(n_det_genes > 0),
    }


def _part3_shard_worker(shard_id: int) -> dict[str, Any]:
    """Read shard_{id} CSR cleanly and compute all Part 3 E1, E2, stringent, 20-cell, strata & 5-SH donor metrics."""
    ds_name = _W_CTX["ds_name"]
    ds_code = 1 if ds_name == "gao" else 2
    n_genes = _W_CTX["n_genes"]
    n_cts = _W_CTX["n_cts"]
    bin_arr = _W_CTX["bin_arr"]
    sets = _W_CTX["sets"]  # E1_full, E1_str, E1_15, E2_g, E2_str, E2_15, E2_strata_0..3

    sdir = _shard_dir(ds_name)
    c_path = sdir / f"shard_{shard_id}_cells.bin"
    if not c_path.exists() or c_path.stat().st_size == 0:
        return {"rows": [], "strata_rows": [], "ztz": np.zeros((n_cts, 2, 5, 5)), "zty": np.zeros((n_cts, 2, 5))}

    cells = np.fromfile(c_path, dtype=np.int64).reshape(-1, 4)
    indices = np.fromfile(sdir / f"shard_{shard_id}_idx.bin", dtype=np.uint16).astype(np.int32)
    data_i32 = np.fromfile(sdir / f"shard_{shard_id}_dat.bin", dtype=np.int32)
    indptr = np.concatenate([[0], np.cumsum(cells[:, 3], dtype=np.int64)])
    csr_i32 = sp.csr_matrix((data_i32, indices, indptr), shape=(len(cells), n_genes))

    d_col, c_col, u_col = cells[:, 0], cells[:, 1], cells[:, 2]
    rows: list[dict[str, Any]] = []
    st_rows: list[dict[str, Any]] = []
    ztz = np.zeros((n_cts, 2, 5, 5), dtype=np.float64)
    zty = np.zeros((n_cts, 2, 5), dtype=np.float64)

    def _get_gidx(ref_obj: dict[str, Any], di: int, ci: int) -> np.ndarray | None:
        if bin_arr[di] == 1:
            return ref_obj["loo"].get((di, ci))
        return ref_obj["full"].get(ci)

    def _score_gene_set(sub_i: sp.csr_matrix, sub_u: np.ndarray, di: int, ci: int, g_arr: np.ndarray, key_tag: str, run_sh: bool = False, run_20c: bool = False) -> dict[str, Any]:
        nc = sub_i.shape[0]
        if len(g_arr) == 0:
            return {}
        tag_crc = zlib.crc32(key_tag.encode("ascii")) & 0xFFFF
        b_full = sub_i[:, g_arr].tocsr()
        cnt_per_cell = np.diff(b_full.indptr).astype(np.float64)
        amb_umi_frac = float(np.asarray(b_full.sum(axis=1)).ravel().sum() / np.maximum(sub_u.sum(), 1))
        rng_p = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, di, ci, tag_crc]))
        bm = compute_binary_ectopic_metrics(b_full, rng_p, len(g_arr), run_excess=True)

        # 20-seed subsampled_20cells for singleton_fraction, evenness, top10_share
        s20_sings, s20_ev, s20_top10, s20_jac = [], [], [], []
        if run_20c:
            for si20 in SUB20_SEEDS:
                r20 = np.random.default_rng(np.random.SeedSequence([int(si20), di, ci, tag_crc]))
                pick = np.arange(nc) if nc == 20 else r20.choice(nc, size=20, replace=False)
                bm20 = compute_binary_ectopic_metrics(b_full[pick], r20, len(g_arr), run_excess=False)
                if np.isfinite(bm20["singleton_fraction"]):
                    s20_sings.append(bm20["singleton_fraction"])
                if np.isfinite(bm20["ectopic_pool_evenness"]):
                    s20_ev.append(bm20["ectopic_pool_evenness"])
                if np.isfinite(bm20["top10_ectopic_share"]):
                    s20_top10.append(bm20["top10_ectopic_share"])
                if np.isfinite(bm20["jaccard_both_pos"]):
                    s20_jac.append(bm20["jaccard_both_pos"])

        res = {
            f"{key_tag}_ectopic_count": float(cnt_per_cell.mean()),
            f"{key_tag}_ectopic_umi_frac": amb_umi_frac,
            f"{key_tag}_jaccard_both_pos": bm["jaccard_both_pos"],
            f"{key_tag}_jaccard_all_pairs": bm["jaccard_all_pairs"],
            f"{key_tag}_jaccard_union_pos": bm["jaccard_union_pos"],
            f"{key_tag}_excess_jaccard": bm["excess_jaccard"],
            f"{key_tag}_singleton_fraction": bm["singleton_fraction"],
            f"{key_tag}_top10_share": bm["top10_ectopic_share"],
            f"{key_tag}_evenness": bm["ectopic_pool_evenness"],
            f"{key_tag}_has_any": float(bm["has_any_ectopic"]),
            f"{key_tag}_singleton_20cells": float(np.mean(s20_sings)) if s20_sings else np.nan,
            f"{key_tag}_evenness_20cells": float(np.mean(s20_ev)) if s20_ev else np.nan,
            f"{key_tag}_top10_20cells": float(np.mean(s20_top10)) if s20_top10 else np.nan,
            f"{key_tag}_jaccard_both_pos_20cells": float(np.mean(s20_jac)) if len(s20_jac) >= 10 else np.nan,
            f"{key_tag}_n_valid_20seeds_sing": len(s20_sings),
            f"{key_tag}_n_genes_set": len(g_arr),
            "_cell_cnt_vec": cnt_per_cell,
        }

        if run_sh:
            r_b = np.repeat(np.arange(nc, dtype=np.int32), np.diff(b_full.indptr))
            c_b = b_full.indices
            d_b = np.rint(b_full.data).astype(np.int64)
            n_g_b = b_full.shape[1]
            # Within-cell x_cg >= 2 repeatability + donor pseudobulk det(H1) vs det(H2) Jaccard + 5 SH deltas
            ge2_m = d_b >= 2
            wc_jac_ge2_list, pb_sh_jac_list = [], []
            for si_sh, seed_sh in enumerate(SH_SEEDS):
                rng_s = np.random.default_rng(np.random.SeedSequence([int(seed_sh), ds_code, di, ci, tag_crc]))
                h1_d = rng_s.binomial(d_b, 0.5)
                h2_d = d_b - h1_d
                h1_b = sp.csr_matrix((np.ones(int((h1_d > 0).sum()), dtype=np.float64), (r_b[h1_d > 0], c_b[h1_d > 0])), shape=(nc, n_g_b))
                h2_b = sp.csr_matrix((np.ones(int((h2_d > 0).sum()), dtype=np.float64), (r_b[h2_d > 0], c_b[h2_d > 0])), shape=(nc, n_g_b))
                cnt1 = float(np.diff(h1_b.indptr).mean())
                cnt2 = float(np.diff(h2_b.indptr).mean())
                bm1 = compute_binary_ectopic_metrics(h1_b, rng_s, n_g_b, run_excess=False)
                bm2 = compute_binary_ectopic_metrics(h2_b, rng_s, n_g_b, run_excess=False)
                res[f"h1_s{si_sh}_{key_tag}_ectopic_count"] = cnt1
                res[f"h2_s{si_sh}_{key_tag}_ectopic_count"] = cnt2
                res[f"h1_s{si_sh}_{key_tag}_jaccard_both_pos"] = bm1["jaccard_both_pos"]
                res[f"h2_s{si_sh}_{key_tag}_jaccard_both_pos"] = bm2["jaccard_both_pos"]
                res[f"h1_s{si_sh}_{key_tag}_singleton_fraction"] = bm1["singleton_fraction"]
                res[f"h2_s{si_sh}_{key_tag}_singleton_fraction"] = bm2["singleton_fraction"]
                # Pseudobulk gene detection Jaccard across the donor's 2 halves
                det1 = np.asarray(h1_b.sum(axis=0)).ravel() > 0
                det2 = np.asarray(h2_b.sum(axis=0)).ravel() > 0
                u12 = int((det1 | det2).sum())
                if u12 > 0:
                    pb_sh_jac_list.append(float((det1 & det2).sum() / u12))
                if si_sh == 0 and np.any(ge2_m):
                    # Restricted strictly to entries x_{cg} >= 2 where split into BOTH halves is possible!
                    inter_2 = ((h1_d[ge2_m] > 0) & (h2_d[ge2_m] > 0)).sum()
                    wc_jac_ge2_list.append(float(inter_2 / int(ge2_m.sum())))
            res[f"{key_tag}_sh_pb_jaccard"] = float(np.mean(pb_sh_jac_list)) if pb_sh_jac_list else np.nan
            res[f"{key_tag}_sh_ge2_repeatability"] = float(np.mean(wc_jac_ge2_list)) if wc_jac_ge2_list else np.nan
        return res

    for di in np.unique(d_col):
        dm = d_col == di
        d_idx = np.flatnonzero(dm)
        for ci in np.unique(c_col[dm]):
            sub_idx = d_idx[c_col[d_idx] == ci]
            nc = len(sub_idx)
            if nc < MIN_CELLS_PER_DONOR:
                continue
            sub_i = csr_i32[sub_idx].tocsr().astype(np.int64)
            sub_u = u_col[sub_idx]
            z_sp = spline_design_matrix(sub_u)
            rec: dict[str, Any] = {
                "donor_idx": int(di), "ct_idx": int(ci), "n_cells": nc,
                "mean_log_umi": float(np.log(sub_u.astype(np.float64)).mean()),
                "mean_inv_umi": float((1.0 / np.maximum(sub_u.astype(np.float64), 1.0)).mean()),
                "log_n_cells": float(np.log(float(nc))),
                "z_mean": z_sp.mean(axis=0),
            }
            for tag, s_key, r_sh, r_20, sp_idx in (
                ("e1", "E1_full", True, True, 0),
                ("e1_str", "E1_str", False, True, -1),
                ("e2_g", "E2_global", True, True, 1),
                ("e2_str", "E2_str", False, True, -1),
            ):
                g_arr = _get_gidx(sets[s_key], int(di), int(ci))
                if g_arr is not None:
                    sc_r = _score_gene_set(sub_i, sub_u, int(di), int(ci), g_arr, tag, run_sh=r_sh, run_20c=r_20)
                    c_vec = sc_r.pop("_cell_cnt_vec", None)
                    rec.update(sc_r)
                    if sp_idx >= 0 and c_vec is not None:
                        ztz[ci, sp_idx] += z_sp.T @ z_sp
                        zty[ci, sp_idx] += z_sp.T @ c_vec

            # 1500-UMI rarefaction sensitivity (E1 & E2) across 5 seeds
            rng_15 = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, ds_code, 15003, int(di), int(ci)]))
            ds15 = multinomial_downsample(sub_i, 1500, rng_15)
            ds15.eliminate_zeros()  # downsample keeps dropped genes as explicit stored zeros; _score_gene_set counts stored entries
            ds15 = ds15.astype(np.int64)
            u15_arr = np.full(nc, 1500, dtype=np.int64)
            for tag, s_key in (("e1_r15", "E1_1500"), ("e2_r15", "E2_1500")):
                g_arr = _get_gidx(sets[s_key], int(di), int(ci))
                if g_arr is not None:
                    sc15 = _score_gene_set(ds15, u15_arr, int(di), int(ci), g_arr, tag, run_sh=False, run_20c=True)
                    sc15.pop("_cell_cnt_vec", None)
                    rec.update(sc15)

            # Bounded depth strata E2_stratum (>=20 cells within stratum s AND min_loo_young_s >= 200)
            st_ids = assign_depth_stratum(sub_u)
            for si, (s_lbl, _, _) in enumerate(DEPTH_STRATA):
                g_arr = _get_gidx(sets[f"E2_st_{si}"], int(di), int(ci))
                if g_arr is None:
                    continue
                sm = st_ids == si
                if int(sm.sum()) < MIN_CELLS_PER_DONOR:
                    continue
                sub_st = sub_i[sm].tocsr()
                u_st = sub_u[sm]
                sc_s = _score_gene_set(sub_st, u_st, int(di), int(ci), g_arr, "e2_st", run_sh=True, run_20c=True)
                sc_s.pop("_cell_cnt_vec", None)
                st_rows.append({
                    "donor_idx": int(di), "ct_idx": int(ci), "stratum": s_lbl,
                    "n_cells": int(sm.sum()), "mean_log_umi": float(np.log(u_st.astype(np.float64)).mean()),
                    "log_n_cells": float(np.log(float(sm.sum()))),
                    **sc_s,
                })

            rows.append(rec)
    return {"rows": rows, "strata_rows": st_rows, "ztz": ztz, "zty": zty}


def run_part3_ectopic_and_qc(all_datasets: list[str] | None = None) -> dict[str, Any]:
    """Run Part 3 Ectopic Identity on Gao & Yazar (reloading Part 1/3 caches, zero HDF5 reads) -> CSVs + QC JSON."""
    t0 = time.perf_counter()
    targets = all_datasets or ["gao", "yazar"]
    qc_summary: dict[str, Any] = {"datasets": {}, "verdicts": {}, "seeds": {"master": MASTER_SEED, "split_half": SH_SEEDS, "sub20": SUB20_SEEDS}}
    rob_all_parts: list[pd.DataFrame] = []

    for ds_name in targets:
        meta = reconcile_gao_metadata(GAO_H5AD) if ds_name == "gao" else reconcile_yazar_metadata(YAZAR_H5AD, YAZAR_CACHE)
        with open(CACHE_DIR / f"{ds_name}_pass1_part3_refs.pkl", "rb") as pf:
            p3r = pickle.load(pf)

        donors_df = meta["donors_df"]
        ct_names = meta["cell_types"]
        n_donors, n_cts = len(meta["donor_ids"]), len(ct_names)
        gene_audit = p3r["gene_audit"]
        n_genes = int(gene_audit["n_genes_retained"])
        bin_arr = donors_df["bin"].map({"LT20": 0, "YOUNG": 1, "MID": 2, "OLD": 3}).to_numpy(dtype=np.int32)
        is_gao = ds_name == "gao"
        not_eval_set = set(meta.get("not_evaluable_all_umi_csv", []))

        # Compute YOUNG CM top-50 ambient mask strictly from YOUNG Cardiomyocytes on Gao
        young_cm_counts = None
        if is_gao and "Cardiomyocyte" in ct_names:
            cm_ci = ct_names.index("Cardiomyocyte")
            young_cm_counts = np.zeros(n_genes, dtype=np.float64)
            for di in range(n_donors):
                if bin_arr[di] == 1:
                    young_cm_counts += p3r["y_full_cnt"][p3r["y_pair_map"][di * n_cts + cm_ci]]
        str_mask = build_stringent_ambient_mask(gene_audit["retained_symbols"], ds_name, young_cm_counts)

        # Build E1 (Full, Stringent, 1500) & E2 (Global, Stringent, 1500-5seed pooled, 4 Strata)
        y_1500_det_mean = np.rint(p3r["y_1500_det"].mean(axis=0)).astype(np.int64)
        sets: dict[str, Any] = {
            "E1_full": build_e1_offprogram_lookup(p3r["y_full_cnt"], p3r["y_full_det"], p3r["y_n_cells"], p3r["y_pair_map"], p3r["cnt_ct_15"], bin_arr, n_donors, n_cts, n_genes, True, None),
            "E1_str": build_e1_offprogram_lookup(p3r["y_full_cnt"], p3r["y_full_det"], p3r["y_n_cells"], p3r["y_pair_map"], p3r["cnt_ct_15"], bin_arr, n_donors, n_cts, n_genes, True, str_mask),
            "E1_1500": build_e1_offprogram_lookup(p3r["y_1500_cnt0"], y_1500_det_mean, p3r["y_n_cells"], p3r["y_pair_map"], p3r["cnt_ct_15"], bin_arr, n_donors, n_cts, n_genes, True, None),
            "E2_global": build_e2_silent_lookup(p3r["y_full_det"], p3r["y_n_cells"], p3r["y_pair_map"], p3r["cnt_ct_15"], bin_arr, n_donors, n_cts, n_genes, MIN_LOO_REF_CELLS, True, None),
            "E2_str": build_e2_silent_lookup(p3r["y_full_det"], p3r["y_n_cells"], p3r["y_pair_map"], p3r["cnt_ct_15"], bin_arr, n_donors, n_cts, n_genes, MIN_LOO_REF_CELLS, True, str_mask),
            "E2_1500": build_e2_silent_lookup(y_1500_det_mean, p3r["y_n_cells"], p3r["y_pair_map"], p3r["cnt_ct_15"], bin_arr, n_donors, n_cts, n_genes, MIN_LOO_REF_CELLS, True, None),
        }
        dummy_cnt_all = np.full((n_donors, n_cts), MIN_CELLS_PER_DONOR, dtype=np.int64)
        for si in range(4):
            sets[f"E2_st_{si}"] = build_e2_silent_lookup(
                p3r["y_st_det"][si], p3r["y_st_n"][si], p3r["y_pair_map"],
                dummy_cnt_all,
                bin_arr, n_donors, n_cts, n_genes, MIN_LOO_YOUNG_E2_STRATUM, False, None,
            )

        del p3r, y_1500_det_mean, young_cm_counts
        gc.collect()

        global _W_CTX
        _W_CTX = {"ds_name": ds_name, "n_genes": n_genes, "n_cts": n_cts, "bin_arr": bin_arr, "sets": sets}
        raw_rows, raw_st_rows = [], []
        ztz = np.zeros((n_cts, 2, 5, 5), dtype=np.float64)
        zty = np.zeros((n_cts, 2, 5), dtype=np.float64)
        try:
            with ProcessPoolExecutor(max_workers=SHARD_WORKERS, initializer=_worker_init, initargs=()) as ex:
                for res in ex.map(_part3_shard_worker, list(range(N_SHARDS))):
                    raw_rows.extend(res["rows"])
                    raw_st_rows.extend(res["strata_rows"])
                    ztz += res["ztz"]
                    zty += res["zty"]
        finally:
            _W_CTX = {}
            gc.collect()

        beta_p3 = np.zeros((n_cts, 2, 5), dtype=np.float64)
        for ci in range(n_cts):
            for k_m in range(2):
                if ztz[ci, k_m, 0, 0] >= 20:
                    beta_p3[ci, k_m] = np.linalg.pinv(ztz[ci, k_m]) @ zty[ci, k_m]

        def _enrich_meta(rlist: list[dict[str, Any]], add_spline: bool = False) -> pd.DataFrame:
            out_r = []
            for r in rlist:
                di, ci = r["donor_idx"], r["ct_idx"]
                dr = donors_df.iloc[di]
                grp = ct_names[ci]
                rec = {
                    "dataset": ds_name, "donor_id": str(dr["donor_id"]),
                    "age_years": float(dr["age_years"]), "age_per_decade": float(dr["age_per_decade"]),
                    "bin": str(dr["bin"]), "sex_male": int(dr["sex_male"]),
                    "study_or_pool": str(dr["study"] if is_gao else int(dr["pool_number"])),
                    "technology": str(dr["technology"]) if is_gao else "10x 3' v2",
                    "group": grp, "is_pooled": False,
                    "yazar_pool_fe_all_evaluable": True if is_gao else (grp not in not_eval_set),
                    "n_cells": r["n_cells"],
                    "_depth": r["mean_log_umi"], "_depth2": r["mean_log_umi"] ** 2,
                    "_inv_umi": r.get("mean_inv_umi", np.nan), "_log_n_cells": r["log_n_cells"],
                    **r,
                }
                if add_spline and "z_mean" in r:
                    zm = r["z_mean"]
                    if "e1_ectopic_count" in r:
                        rec["e1_spline_count"] = r["e1_ectopic_count"] - float(zm @ beta_p3[ci, 0])
                    if "e2_g_ectopic_count" in r:
                        rec["e2_g_spline_count"] = r["e2_g_ectopic_count"] - float(zm @ beta_p3[ci, 1])
                out_r.append(rec)
            return pd.DataFrame(out_r)

        df_p3 = _enrich_meta(raw_rows, add_spline=True)
        df_st3 = _enrich_meta(raw_st_rows, add_spline=False)

        bcol, ref_lv, min_df = "study_or_pool", ("ENCODE v4 (Snyder)" if is_gao else None), (5 if is_gao else 10)
        p3_models: list[dict[str, Any]] = []

        # Fit models across E1 and E2 outcomes
        metric_specs = [
            ("E1", "mean_ectopic_count", "full_depth", "primary_full_depth", "e1_ectopic_count", True, True),
            ("E1", "mean_ectopic_count", "spline_resid_full", "sensitivity_spline_resid", "e1_spline_count", False, False),
            ("E1", "mean_ectopic_count", "rarefied_1500", "sensitivity_rarefied_1500", "e1_r15_ectopic_count", False, False),
            ("E1_stringent", "mean_ectopic_count", "full_depth", "sensitivity_stringent", "e1_str_ectopic_count", False, False),
            ("E1", "singleton_fraction_20cells", "subsampled_20cells", "primary_20cells", "e1_singleton_20cells", False, False),
            ("E1", "singleton_fraction", "full_depth", "sensitivity_full_cells", "e1_singleton_fraction", True, False),
            ("E1", "jaccard_both_pos", "full_depth", "primary_full_depth", "e1_jaccard_both_pos", True, False),
            ("E1", "excess_jaccard", "full_depth", "sensitivity_excess_jaccard", "e1_excess_jaccard", False, False),
            ("E1", "jaccard_all_pairs", "full_depth", "sensitivity_all_pairs", "e1_jaccard_all_pairs", False, False),
            ("E1", "ectopic_pool_evenness_20cells", "subsampled_20cells", "sensitivity_evenness", "e1_evenness_20cells", False, False),
            ("E1", "top10_ectopic_share_20cells", "subsampled_20cells", "sensitivity_top10", "e1_top10_20cells", False, False),
            ("E1", "has_any_ectopic", "full_depth", "diagnostic_missingness", "e1_has_any", False, False),
            ("E2_global", "mean_ectopic_count", "full_depth", "sensitivity_global_depth_vuln", "e2_g_ectopic_count", True, True),
            ("E2_global", "mean_ectopic_count", "rarefied_1500", "sensitivity_rarefied_1500", "e2_r15_ectopic_count", False, False),
            ("E2_global_stringent", "mean_ectopic_count", "full_depth", "sensitivity_stringent", "e2_str_ectopic_count", False, False),
            ("E2_global", "singleton_fraction_20cells", "subsampled_20cells", "sensitivity_global_depth_vuln", "e2_g_singleton_20cells", False, False),
            ("E2_global", "jaccard_both_pos", "full_depth", "sensitivity_global_depth_vuln", "e2_g_jaccard_both_pos", True, False),
            ("E2_global", "excess_jaccard", "full_depth", "sensitivity_excess_jaccard", "e2_g_excess_jaccard", False, False),
            ("E2_global", "ectopic_pool_evenness_20cells", "subsampled_20cells", "sensitivity_evenness", "e2_g_evenness_20cells", False, False),
            ("E2_global", "top10_ectopic_share_20cells", "subsampled_20cells", "sensitivity_top10", "e2_g_top10_20cells", False, False),
            ("E2_rarefied_1500", "singleton_fraction_20cells", "rarefied_1500", "primary_rarefied_1500", "e2_r15_singleton_20cells", False, False),
            ("E2_rarefied_1500", "jaccard_both_pos", "rarefied_1500", "primary_rarefied_1500", "e2_r15_jaccard_both_pos", False, False),
            ("E2_rarefied_1500", "excess_jaccard", "rarefied_1500", "sensitivity_excess_jaccard", "e2_r15_excess_jaccard", False, False),
        ]

        for g_set, m_nm, d_mode, a_role, col_nm, do_sh, is_cnt in metric_specs:
            if col_nm not in df_p3.columns:
                continue
            elig_cts = []
            scale_map: dict[str, float] = {}
            for grp in ct_names:
                sub = df_p3[(df_p3["group"] == grp) & df_p3[col_nm].notna()].copy()
                if len(sub) < MIN_DONORS_SLOPE:
                    continue
                fb_status = "not_evaluable_within_pool_df" if (not is_gao and not bool(sub["yazar_pool_fe_all_evaluable"].iloc[0])) else "not_evaluable_within_pool_df_slice"
                sh_s = {}
                if do_sh:
                    tag_p = "e1" if col_nm.startswith("e1") else "e2_g"
                    suf_p = "ectopic_count" if is_cnt else ("jaccard_both_pos" if "jaccard" in m_nm else "singleton_fraction")
                    h1cs = [f"h1_s{i}_{tag_p}_{suf_p}" for i in range(len(SH_SEEDS))]
                    h2cs = [f"h2_s{i}_{tag_p}_{suf_p}" for i in range(len(SH_SEEDS))]
                    sh_s = compute_sh_and_quenouille_stats(sub, col_nm, h1cs, h2cs, ["sex_male", "_depth"], bcol, 2, min_df, ref_lv, is_cnt, False, (not is_gao))

                for spec_name, sp_df, nbc, sp_bc, mbd, mdf, rlv in _model_specs_for_dataset(sub, ds_name, extra_log_ncells=True):
                    row = fit_ols_model_row(
                        sp_df, col_nm, "age_per_decade", nbc, sp_bc, mbd, mdf, rlv,
                        meta_fields={
                            "dataset": ds_name, "part": "3", "metric": m_nm, "gene_set": g_set,
                            "depth_mode": d_mode, "analysis_role": a_role, "model_spec": spec_name,
                            "group": grp, "is_pooled": False, "is_stratified_global": False,
                            "default_fallback_status": fb_status,
                        },
                        fallback_to_nobatch=(not is_gao),
                    )
                    if row is None:
                        continue
                    oy = fit_contrast_adj_diff(sp_df, col_nm, "YOUNG", "OLD", nbc, sp_bc, mbd, mdf, rlv)
                    row.update({"old_minus_young_diff_adj": oy["diff"], "old_minus_young_se_adj": oy["se"], "old_minus_young_p_adj": oy["p"]})
                    if do_sh and spec_name == "primary":
                        row.update(sh_s)
                    p3_models.append(row)
                    if spec_name == "primary" and row["yazar_pool_fe_status"] == "evaluated":
                        elig_cts.append(grp)
                        scale_map[grp] = float(sub[col_nm].std(ddof=1))

            if len(elig_cts) >= 2:
                st_sub = df_p3[(df_p3["group"].isin(elig_cts)) & df_p3[col_nm].notna()].copy()
                srow, _ = fit_stratified_global_with_batch(
                    st_sub, elig_cts, col_nm, "n_cells", "_depth", ds_name, "3", m_nm, d_mode, a_role, gene_set=g_set,
                )
                if srow is not None:
                    if do_sh:
                        tag_p = "e1" if col_nm.startswith("e1") else "e2_g"
                        suf_p = "ectopic_count" if is_cnt else ("jaccard_both_pos" if "jaccard" in m_nm else "singleton_fraction")
                        fac = 1.0 if is_cnt else 0.5
                        v_sc = 1.0 if is_cnt else 0.25
                        se_sq, se_oysq, rms_v = [], [], []
                        for si_sh in range(len(SH_SEEDS)):
                            h1c, h2c = f"h1_s{si_sh}_{tag_p}_{suf_p}", f"h2_s{si_sh}_{tag_p}_{suf_p}"
                            wsh = st_sub.dropna(subset=[h1c, h2c]).copy()
                            if wsh.empty:
                                continue
                            wsh["_d_sh"] = wsh[h1c] - wsh[h2c]
                            rms_v.append(float(np.mean((wsh["_d_sh"] / wsh["group"].map(scale_map)) ** 2)))
                            sr_sh, _ = fit_stratified_global_with_batch(wsh, elig_cts, "_d_sh", "n_cells", "_depth", ds_name, "3", m_nm, d_mode, "split_half_se", custom_scale_map=scale_map)
                            if sr_sh is not None:
                                se_sq.append(float(sr_sh["se_raw"]) ** 2)
                                if np.isfinite(sr_sh["old_minus_young_se_adj"]):
                                    se_oysq.append(float(sr_sh["old_minus_young_se_adj"]) ** 2)
                        se_samp = float(fac * np.sqrt(np.mean(se_sq))) if se_sq else np.nan
                        se_samp_oy = float(fac * np.sqrt(np.mean(se_oysq))) if se_oysq else np.nan
                        rms_sd = float(np.sqrt(np.mean(rms_v))) if rms_v else np.nan
                        srow.update({
                            "sh_rms_delta_donor": rms_sd, "sh_se_samp_slope": se_samp, "sh_se_samp_oy": se_samp_oy,
                            "sh_rho_tech2": float(v_sc * (rms_sd ** 2) / (srow["resid_sd"] ** 2)) if np.isfinite(rms_sd) else np.nan,
                            "sh_snr_slope": float(abs(srow["slope_raw_per_decade"]) / se_samp) if np.isfinite(se_samp) and se_samp > 0 else np.nan,
                        })
                    p3_models.append(srow)

        # Primary depth-matched E2_stratum models ((CT x stratum) intercepts + z-scoring + within-stratum log_umi + batch FE!)
        if not df_st3.empty:
            for m_nm, col_s, d_m in (
                ("mean_ectopic_count", "e2_st_ectopic_count", "bounded_depth_strata"),
                ("singleton_fraction_20cells", "e2_st_singleton_20cells", "bounded_depth_strata_20cells"),
                ("jaccard_both_pos", "e2_st_jaccard_both_pos", "bounded_depth_strata"),
                ("excess_jaccard", "e2_st_excess_jaccard", "bounded_depth_strata"),
                ("ectopic_pool_evenness_20cells", "e2_st_evenness_20cells", "bounded_depth_strata_20cells"),
                ("top10_ectopic_share_20cells", "e2_st_top10_20cells", "bounded_depth_strata_20cells"),
            ):
                if col_s not in df_st3.columns:
                    continue
                sub_st = df_st3.dropna(subset=[col_s]).copy()
                elig_cts_s = []
                st_scale_map: dict[Any, float] = {}
                for (grp, st_nm), gpart in sub_st.groupby(["group", "stratum"]):
                    if len(gpart) >= MIN_DONORS_SLOPE and (grp not in not_eval_set if not is_gao else True):
                        sd_g = float(gpart[col_s].std(ddof=1))
                        if np.isfinite(sd_g) and sd_g > 0.0:
                            elig_cts_s.append(grp)
                            st_scale_map[(grp, st_nm)] = sd_g
                elig_cts_u = sorted(set(elig_cts_s))
                if len(elig_cts_u) >= 2:
                    srow, _ = fit_stratified_global_with_batch(
                        sub_st, elig_cts_u, col_s, "n_cells", "_depth", ds_name, "3", m_nm, d_m, "primary_depth_matched_e2",
                        model_spec="stratified_global_primary", gene_set="E2_stratum", strata_col="stratum",
                    )
                    if srow is not None:
                        if m_nm in ("mean_ectopic_count", "jaccard_both_pos"):
                            is_cnt_s = m_nm == "mean_ectopic_count"
                            suf_s = "ectopic_count" if is_cnt_s else "jaccard_both_pos"
                            fac_s, v_scs = (1.0, 1.0) if is_cnt_s else (0.5, 0.25)
                            se_sq, se_oysq, rms_v = [], [], []
                            for si_sh in range(len(SH_SEEDS)):
                                h1c, h2c = f"h1_s{si_sh}_e2_st_{suf_s}", f"h2_s{si_sh}_e2_st_{suf_s}"
                                if h1c not in sub_st.columns or h2c not in sub_st.columns:
                                    continue
                                wsh = sub_st.dropna(subset=[h1c, h2c]).copy()
                                if wsh.empty:
                                    continue
                                wsh["_d_sh"] = wsh[h1c] - wsh[h2c]
                                sc_keys = list(zip(wsh["group"], wsh["stratum"]))
                                sc_vals = np.array([st_scale_map.get(k, np.nan) for k in sc_keys], dtype=np.float64)
                                valid_sc = np.isfinite(sc_vals) & (sc_vals > 0)
                                if np.any(valid_sc):
                                    rms_v.append(float(np.mean((wsh["_d_sh"].to_numpy()[valid_sc] / sc_vals[valid_sc]) ** 2)))
                                sr_sh, _ = fit_stratified_global_with_batch(
                                    wsh, elig_cts_u, "_d_sh", "n_cells", "_depth", ds_name, "3", m_nm, d_m, "split_half_se",
                                    gene_set="E2_stratum", strata_col="stratum", custom_scale_map=st_scale_map,
                                )
                                if sr_sh is not None:
                                    se_sq.append(float(sr_sh["se_raw"]) ** 2)
                                    if np.isfinite(sr_sh["old_minus_young_se_adj"]):
                                        se_oysq.append(float(sr_sh["old_minus_young_se_adj"]) ** 2)
                            se_samp = float(fac_s * np.sqrt(np.mean(se_sq))) if se_sq else np.nan
                            se_samp_oy = float(fac_s * np.sqrt(np.mean(se_oysq))) if se_oysq else np.nan
                            rms_sd = float(np.sqrt(np.mean(rms_v))) if rms_v else np.nan
                            srow.update({
                                "sh_rms_delta_donor": rms_sd, "sh_se_samp_slope": se_samp, "sh_se_samp_oy": se_samp_oy,
                                "sh_rho_tech2": float(v_scs * (rms_sd ** 2) / (srow["resid_sd"] ** 2)) if np.isfinite(rms_sd) else np.nan,
                                "sh_snr_slope": float(abs(srow["slope_raw_per_decade"]) / se_samp) if np.isfinite(se_samp) and se_samp > 0 else np.nan,
                            })
                        p3_models.append(srow)

        mod3_df = apply_slice_bh(pd.DataFrame(p3_models), ["dataset", "part", "metric", "gene_set", "depth_mode", "model_spec"])
        write_clean_csv(mod3_df, RESULTS_DIR / f"followup-fulldepth-{ds_name}-ectopic.csv")

        rob2_path = CACHE_DIR / f"{ds_name}_rob2.pkl"
        if rob2_path.exists():
            with open(rob2_path, "rb") as pf:
                rob2_df = pickle.load(pf)
            if isinstance(rob2_df, pd.DataFrame) and not rob2_df.empty:
                rob_all_parts.append(rob2_df)

        # Load Part 1 & Part 2 CSVs for verdict evaluation + QC JSON
        p1_df = pd.read_csv(RESULTS_DIR / f"followup-fulldepth-{ds_name}-cell-drift.csv")
        p2_df = pd.read_csv(RESULTS_DIR / f"followup-fulldepth-{ds_name}-spread.csv")
        qc_summary["datasets"][ds_name] = {
            "gene_audit": {k: v for k, v in gene_audit.items() if k != "retained_symbols"},
            "n_e1_cts": len(sets["E1_full"]["summary"]),
            "e1_median_markers": float(np.median(list(sets["E1_full"]["summary"].values()))) if sets["E1_full"]["summary"] else 0.0,
            "e2_global_median_markers": float(np.median(list(sets["E2_global"]["summary"].values()))) if sets["E2_global"]["summary"] else 0.0,
            "e1_sh_pb_jaccard_median": float(df_p3["e1_sh_pb_jaccard"].median()) if "e1_sh_pb_jaccard" in df_p3.columns and df_p3["e1_sh_pb_jaccard"].notna().any() else np.nan,
            "e2_g_sh_pb_jaccard_median": float(df_p3["e2_g_sh_pb_jaccard"].median()) if "e2_g_sh_pb_jaccard" in df_p3.columns and df_p3["e2_g_sh_pb_jaccard"].notna().any() else np.nan,
            "e1_sh_ge2_repeat_median": float(df_p3["e1_sh_ge2_repeatability"].median()) if "e1_sh_ge2_repeatability" in df_p3.columns and df_p3["e1_sh_ge2_repeatability"].notna().any() else np.nan,
            "e2_st_ectopic_umi_frac_median": float(df_st3["e2_st_ectopic_umi_frac"].median()) if (not df_st3.empty and "e2_st_ectopic_umi_frac" in df_st3.columns and df_st3["e2_st_ectopic_umi_frac"].notna().any()) else np.nan,
        }
        qc_summary["verdicts"][ds_name] = evaluate_fulldepth_verdicts(p1_df, p2_df, mod3_df)

    rob_all = pd.concat(rob_all_parts, ignore_index=True) if rob_all_parts else pd.DataFrame()
    write_clean_csv(rob_all, RESULTS_DIR / "followup-fulldepth-robustness.csv")
    qc_txt = json.dumps(qc_summary, indent=2) + "\n"
    assert ("/home/" + "reza") not in qc_txt and "\u2013" not in qc_txt and "\u2014" not in qc_txt
    (RESULTS_DIR / "followup-fulldepth-qc.json").write_text(qc_txt, encoding="utf-8")
    _log("all", f"Part 3 + Verdicts + QC JSON complete in {round(time.perf_counter() - t0, 2)}s")
    return qc_summary


def evaluate_fulldepth_verdicts(p1_df: pd.DataFrame, p2_df: pd.DataFrame, p3_df: pd.DataFrame) -> dict[str, Any]:
    """Strict, uncollapsed Verdict engine per dataset & per part (and E1/E2 separately)."""
    out: dict[str, Any] = {}

    def _g_row(df: pd.DataFrame, part: str, metric: str, depth_mode: str, spec: str = "stratified_global_primary", g_set: str = "none") -> pd.Series | None:
        m = (df["part"].astype(str) == str(part)) & (df["metric"] == metric) & (df["depth_mode"] == depth_mode) & (df["model_spec"] == spec) & (df["gene_set"] == g_set)
        sub = df[m]
        return None if sub.empty else sub.iloc[0]

    # Part 1: center_shift and scatter strictly separate
    for metric in ("center_shift", "scatter"):
        rp = _g_row(p1_df, "1", metric, "full_depth")
        rs = _g_row(p1_df, "1", metric, "spline_resid_full")
        if rp is None:
            out[f"Part1_{metric}"] = {"verdict": "INCONCLUSIVE", "reason": "Missing global stratified row"}
            continue
        p_ok = float(rp["slope_raw_per_decade"]) > 0.0 and float(rp["pvalue"]) < 0.05
        oy_ok = float(rp["old_minus_young_diff_adj"]) > 0.0 and float(rp["old_minus_young_p_adj"]) < 0.05
        snr_ok = np.isfinite(rp.get("sh_snr_slope", np.nan)) and float(rp["sh_snr_slope"]) > 2.0
        q_ok = (not np.isfinite(rp.get("quenouille_bias_frac", np.nan))) or (float(rp["quenouille_bias_frac"]) < 0.25 and float(rp.get("quenouille_corr_slope", 0.0)) > 0.0)
        sp_ok = (rs is None) or (float(rs["slope_raw_per_decade"]) > 0.0 and float(rs["pvalue"]) < 0.05)
        extra_ok = True if metric == "scatter" else (float(rp.get("old_mean", np.nan)) > float(rp.get("young_loo_mean", np.nan)))
        if p_ok and oy_ok and snr_ok and q_ok and sp_ok and extra_ok:
            out[f"Part1_{metric}"] = {"verdict": "PROVES", "reason": f"Global {metric} >0 (p={rp['pvalue']:.3e}, OY p={rp['old_minus_young_p_adj']:.3e}, SNR={rp['sh_snr_slope']:.2f})"}
        elif float(rp["mde80_sd_per_decade"]) <= 0.5:
            tier = "high_power<=0.25" if float(rp["mde80_sd_per_decade"]) <= 0.25 else "moderate_power<=0.5"
            out[f"Part1_{metric}"] = {"verdict": "DOES NOT PROVE", "reason": f"Global {metric} fails PROVES (slope={rp['slope_raw_per_decade']:.4f}, p={rp['pvalue']:.3g}, MDE80={rp['mde80_sd_per_decade']:.3f} [{tier}])"}
        else:
            out[f"Part1_{metric}"] = {"verdict": "INCONCLUSIVE", "reason": f"MDE80={rp['mde80_sd_per_decade']:.3f}>0.5"}

    # Part 2: spread (per metric)
    p2_out: dict[str, Any] = {}
    for metric, sign in (("shannon", 1), ("deman_kl_methods", -1), ("deman_kl_sqrt", -1)):
        rp = _g_row(p2_df, "2", metric, "full_depth")
        r3k = _g_row(p2_df, "2", metric, "rarefied_3000")
        rs = _g_row(p2_df, "2", metric, "spline_resid_full")
        if rp is None:
            p2_out[metric] = {"verdict": "INCONCLUSIVE", "reason": "Missing full-depth global row"}
            continue
        dir_f = (sign * float(rp["slope_raw_per_decade"]) > 0.0) and float(rp["pvalue"]) < 0.05
        dir_3k = (r3k is not None) and (sign * float(r3k["slope_raw_per_decade"]) > 0.0)
        snr_ok = np.isfinite(rp.get("sh_snr_slope", np.nan)) and float(rp["sh_snr_slope"]) > 2.0
        q_ok = (not np.isfinite(rp.get("quenouille_bias_frac", np.nan))) or (float(rp["quenouille_bias_frac"]) < 0.25 and sign * float(rp.get("quenouille_corr_slope", 0.0)) > 0.0)
        sp_ok = (rs is None) or (sign * float(rs["slope_raw_per_decade"]) > 0.0 and float(rs["pvalue"]) < 0.05)
        if dir_f and dir_3k and snr_ok and q_ok and sp_ok:
            p2_out[metric] = {"verdict": "PROVES", "reason": f"Full-depth global {metric} more spread (p={rp['pvalue']:.3e}, SNR={rp['sh_snr_slope']:.2f}) & same direction at 3000"}
        elif float(rp["mde80_sd_per_decade"]) <= 0.5:
            tier = "high_power<=0.25" if float(rp["mde80_sd_per_decade"]) <= 0.25 else "moderate_power<=0.5"
            p2_out[metric] = {"verdict": "DOES NOT PROVE", "reason": f"{metric} fails PROVES (slope={rp['slope_raw_per_decade']:.4f}, p={rp['pvalue']:.3g}, MDE80={rp['mde80_sd_per_decade']:.3f} [{tier}])"}
        else:
            p2_out[metric] = {"verdict": "INCONCLUSIVE", "reason": f"MDE80={rp['mde80_sd_per_decade']:.3f}>0.5"}
    out["Part2"] = p2_out

    # Part 3: Ectopic Identity (E1 and E2 separately with Holm check across {E1, E2})
    p3_out: dict[str, Any] = {}
    for g_tag, cnt_d, sin_g, sin_d, jac_g, jac_d, r15_g, str_g in (
        ("E1", "full_depth", "E1", "subsampled_20cells", "E1", "full_depth", "E1", "E1_stringent"),
        ("E2", "bounded_depth_strata", "E2_stratum", "bounded_depth_strata_20cells", "E2_stratum", "bounded_depth_strata", "E2_global", "E2_global_stringent"),
    ):
        g_cnt_src = "E1" if g_tag == "E1" else "E2_stratum"
        rc = _g_row(p3_df, "3", "mean_ectopic_count", cnt_d, g_set=g_cnt_src)
        rsn = _g_row(p3_df, "3", "singleton_fraction_20cells", sin_d, g_set=sin_g)
        rjc = _g_row(p3_df, "3", "jaccard_both_pos", jac_d, g_set=jac_g)
        rexc = _g_row(p3_df, "3", "excess_jaccard", jac_d, g_set=jac_g)
        r15 = _g_row(p3_df, "3", "mean_ectopic_count", "rarefied_1500", g_set=r15_g)
        rstr = _g_row(p3_df, "3", "mean_ectopic_count", "full_depth", g_set=str_g)
        if rc is None:
            p3_out[g_tag] = {"verdict": "INCONCLUSIVE", "reason": "Missing global ectopic count row"}
            continue
        cnt_ok = float(rc["slope_raw_per_decade"]) > 0.0 and float(rc["pvalue"]) < 0.05
        snr_cnt_ok = np.isfinite(rc.get("sh_snr_slope", np.nan)) and float(rc["sh_snr_slope"]) > 2.0
        r15_ok = (r15 is not None) and float(r15["slope_raw_per_decade"]) > 0.0 and float(r15["pvalue"]) < 0.05
        str_ok = (rstr is not None) and float(rstr["slope_raw_per_decade"]) > 0.0 and float(rstr["pvalue"]) < 0.05
        sin_up = (rsn is not None) and float(rsn["slope_raw_per_decade"]) > 0.0 and float(rsn["pvalue"]) < 0.05
        jac_dn = (rjc is not None) and float(rjc["slope_raw_per_decade"]) < 0.0 and float(rjc["pvalue"]) < 0.05
        exc_up = (rexc is not None) and float(rexc["slope_raw_per_decade"]) > 0.0 and float(rexc["pvalue"]) < 0.05
        if cnt_ok and snr_cnt_ok and r15_ok and str_ok and (sin_up or jac_dn) and (not exc_up):
            p3_out[g_tag] = {"verdict": "PROVES (RANDOMNESS)", "reason": f"{g_tag} count UP (SNR={rc['sh_snr_slope']:.2f}) + 1500/stringent + randomness dispersion confirmed"}
        elif cnt_ok and snr_cnt_ok and exc_up:
            p3_out[g_tag] = {"verdict": "DOES NOT PROVE RANDOMNESS (SHARED_INDUCTION_OR_AMBIENT)", "reason": f"{g_tag} count UP (SNR={rc['sh_snr_slope']:.2f}) but excess_jaccard RISES significantly (p={rexc['pvalue']:.3e})"}
        elif float(rc["mde80_sd_per_decade"]) <= 0.5:
            p3_out[g_tag] = {"verdict": "DOES NOT PROVE", "reason": f"{g_tag} fails randomness PROVES (cnt p={rc['pvalue']:.3g}, SNR={rc.get('sh_snr_slope', np.nan):.2f}, MDE80={rc['mde80_sd_per_decade']:.3f}<=0.5)"}
        else:
            p3_out[g_tag] = {"verdict": "INCONCLUSIVE", "reason": f"{g_tag} MDE80={rc['mde80_sd_per_decade']:.3f}>0.5"}
    out["Part3"] = p3_out
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["gao", "yazar", "all"], default="all")
    parser.add_argument("--stage", choices=["parts12", "part3", "all"], default="all")
    parser.add_argument("--reuse-pass1", action="store_true", help="Reuse cached Pass 1 shards and Part 2 summary if complete")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    targets = ["gao", "yazar"] if args.dataset == "all" else [args.dataset]
    if args.stage in ("parts12", "all"):
        for ds in targets:
            _log(ds, "Starting Parts 1 + 2 full-depth run")
            meta = reconcile_gao_metadata(GAO_H5AD) if ds == "gao" else reconcile_yazar_metadata(YAZAR_H5AD, YAZAR_CACHE)
            p2_cache = CACHE_DIR / f"{ds}_pass1_part2.pkl"
            p1_cache = CACHE_DIR / f"{ds}_pass1_part1_refs.pkl"
            if args.reuse_pass1 and p2_cache.exists() and p1_cache.exists() and all((_shard_dir(ds) / f"shard_{s}_cells.bin").exists() for s in range(N_SHARDS)):
                _log(ds, "Reusing cached Pass 1 shards and Part 2 checkpoint")
                with open(p2_cache, "rb") as pf:
                    p2_ckpt = pickle.load(pf)
            else:
                p2_ckpt = run_pass1_and_export_part2(meta, ds)
            mod2_df, rob2_df, _ = _build_part2_table_and_fit(meta, p2_ckpt)
            write_clean_csv(mod2_df, RESULTS_DIR / f"followup-fulldepth-{ds}-spread.csv")
            # Save Part 2 dataset robustness into cache for Commit C robustness.csv
            with open(CACHE_DIR / f"{ds}_rob2.pkl", "wb") as pf:
                pickle.dump(rob2_df, pf, protocol=pickle.HIGHEST_PROTOCOL)
            run_pass2_part1_and_export(meta, ds, mod2_df)
    if args.stage in ("part3", "all"):
        run_part3_ectopic_and_qc(targets)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
