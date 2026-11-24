#!/usr/bin/env python3
"""Phase 2 Audit Item 5 Batch A: mature (30-plus) follow-up organ transports.

Cohorts in order:
  1. Gao heart (GSE290367)
  2. Yazar blood (HCATA OneK1K PBMC, GSE196830)
  3. Perez blood controls (HCATA SLE/Control PBMC, GSE174188)
  4. Guilliams liver sorted cells (HCATA Liver Cell Atlas, GSE192742)

Reuses the exact original follow-up statistical rules (pseudobulk log2(CPM+1),
technical-gene removal, CPM>=1 in >=50% mature donors, original covariate models,
2000-draw expression-decile-matched null with seed 20261105, 1000-perm concordance
with seed 20261105), restricted strictly to donors aged >= 30.0 years BEFORE
pseudobulk aggregation and detection filtering, and scored with mature P1m and P0m.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import re
import resource
import shutil
import subprocess
import sys
import time
from typing import Any

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse as sp
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_organ_17_heart_gao import (
    H5AD_PATH as GAO_H5AD,
    STUDY_DUMMY_MAP as GAO_STUDY_DUMMY_MAP,
    build_donor_audit as gao_build_donor_audit,
    build_donor_meta_table as gao_build_donor_meta_table,
    stream_pseudobulk_counts as gao_stream_pseudobulk_counts,
)
from src.followup_muscle import (
    count_gate,
    expression_deciles,
    fit_ols,
    matched_null_effect_pvalue,
    random_signature_scores as orig_random_signature_scores,
)
from src.followup_organs import (
    assert_no_forbidden_dashes,
    finalize_pseudobulk_group,
    run_family_a_clr,
    sanitize_dashes,
    select_design_columns,
    write_clean_csv,
)
from src.reviewer_stats import pairwise_pearson

RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
HCATA_DIR = DATA_ROOT / "tmp_organs" / "hcata"
YAZAR_H5AD = HCATA_DIR / "yazar_1e44db10.h5ad"
YAZAR_CACHE = HCATA_DIR / "yazar_cache"
TMP_WORK_DIR = DATA_ROOT / "tmp_organs" / "hcata_mature_batch_a_tmp"

PEREZ_CZ_URL = "https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad"
PEREZ_EXPECTED_BYTES = 12218105530
GUILLIAMS_CZ_URL = "https://datasets.cellxgene.cziscience.com/a2ebec68-7e59-4c44-a88f-c9ef87139d2e.h5ad"
GUILLIAMS_EXPECTED_BYTES = 1026110286

MIN_AGE = 30.0
MIN_CELLS_PER_DONOR = 20
ORIG_NULL_SEED = 20261105
FROZEN_DEFAULT_SEED = 20261108
P1_ITERATIONS = 2000
CONCORDANCE_PERMS = 1000
MDE80_NULL_FACTOR = 2.4865
ALPHA = 0.05
MIN_FREE_DISK_GIB = 30.0
MIN_AVAIL_RAM_GIB = 12.0
RAM_WAIT_TIMEOUT_SEC = 3600
MODEL_USED = os.environ.get("AGESPACE_MODEL_USED", "direct execution (no model)")
YAZAR_DOUBLETS_LABEL = "peripheral blood mononuclear cell"
CHUNK_ROWS = 25000
MAX_GROUP_BATCH = 6

PREFIX = "followup-mature-organ-transports-a"
RUNTIME_CSV = RESULTS_DIR / f"{PREFIX}-runtime.csv"
SUMMARY_CSV = RESULTS_DIR / f"{PREFIX}-summary.csv"
BH_CHECK_CSV = RESULTS_DIR / f"{PREFIX}-bh-check.csv"

REQUIRED_COHORT_COLS = [
    "cohort",
    "cell_type",
    "endpoint_family",
    "n_donors",
    "n_donors_30plus_excluded",
    "original_p",
    "original_q",
    "mature_p",
    "mature_q",
    "mature_empirical_p",
    "mature_empirical_q",
    "verdict",
    "notes",
]

ORIGINAL_RESULT_FILES = [
    "results/followup-mature-p1m-genes.csv",
    "results/followup-mature-p1m-p0m-genes.csv",
    "results/followup-organ-17-heart-gao-donors.csv",
    "results/followup-organ-17-heart-gao-cellcounts.csv",
    "results/followup-organ-17-heart-gao-fractions.csv",
    "results/followup-organ-17-heart-gao-genes-summary.csv",
    "results/followup-organ-17-heart-gao-top-genes.csv",
    "results/followup-organ-17-heart-gao-concordance.csv",
    "results/followup-organ-17-heart-gao-p1-signature.csv",
    "results/followup-organ-17-heart-gao-qc.csv",
    "results/followup-organ-17-heart-gao-metadata-summary.json",
    "results/followup-hcata-yazar-donors.csv",
    "results/followup-hcata-yazar-fractions.csv",
    "results/followup-hcata-yazar-genes-summary.csv",
    "results/followup-hcata-yazar-concordance.csv",
    "results/followup-hcata-yazar-p1-signature.csv",
    "results/followup-hcata-yazar-pool-fe-eligibility.csv",
    "results/followup-hcata-yazar-qc.csv",
    "results/followup-hcata-pool-sensitivity-fractions.csv",
    "results/followup-hcata-pool-sensitivity-genes-summary.csv",
    "results/followup-hcata-pool-sensitivity-concordance.csv",
    "results/followup-hcata-pool-sensitivity-p1-signature.csv",
    "results/followup-hcata-metadata-summary.json",
    "results/followup-hcata-perez-donors.csv",
    "results/followup-hcata-perez-fractions.csv",
    "results/followup-hcata-perez-genes-summary.csv",
    "results/followup-hcata-perez-concordance.csv",
    "results/followup-hcata-perez-p1-signature.csv",
    "results/followup-hcata-perez-qc.csv",
    "results/followup-hcata-guilliams-donors.csv",
    "results/followup-hcata-guilliams-genes-summary.csv",
    "results/followup-hcata-guilliams-p1-signature.csv",
    "results/followup-hcata-guilliams-qc.csv",
    "results/followup-hcata-guilliams-sensitivities-genes.csv",
    "results/followup-hcata-guilliams-sensitivities-p1.csv",
    "results/followup-11pt-block5-p0-transport.csv",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def peak_rss_mb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0


def available_ram_gib() -> float:
    with open("/proc/meminfo", encoding="ascii") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024.0 ** 2)
    raise RuntimeError("Could not read MemAvailable from /proc/meminfo")


def wait_for_memory(stage: str, min_gib: float = MIN_AVAIL_RAM_GIB, timeout_sec: int = RAM_WAIT_TIMEOUT_SEC) -> float:
    deadline = time.monotonic() + timeout_sec
    while True:
        avail = available_ram_gib()
        if avail >= min_gib:
            print(f"[RAM OK] {stage}: {avail:.2f} GiB available (>= {min_gib:.1f} GiB)", flush=True)
            return avail
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise MemoryError(f"{stage}: available RAM {avail:.2f} GiB stayed below {min_gib:.1f} GiB for {timeout_sec}s")
        sleep_s = min(30.0, remaining)
        print(f"[RAM WAIT] {stage}: {avail:.2f} GiB available < {min_gib:.1f} GiB; sleeping {sleep_s:.0f}s", flush=True)
        time.sleep(sleep_s)


def append_runtime(cohort: str, step: str, start_utc: str, t0: float, status: str, notes: str = "") -> None:
    row = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "step": step,
                "start_utc": start_utc,
                "end_utc": utc_now(),
                "duration_seconds": round(time.monotonic() - t0, 3),
                "peak_rss_mb": round(peak_rss_mb(), 2),
                "model_used": MODEL_USED,
                "status": status,
                "notes": sanitize_dashes(notes),
            }
        ]
    )
    if RUNTIME_CSV.exists():
        prev = pd.read_csv(RUNTIME_CSV)
        prev = prev[~((prev["cohort"] == cohort) & (prev["step"] == step))]
        out = pd.concat([prev, row], ignore_index=True)
    else:
        out = row
    write_clean_csv(out, RUNTIME_CSV)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def record_sha256_manifest() -> dict[str, str]:
    manifest: dict[str, str] = {}
    for rel in ORIGINAL_RESULT_FILES:
        p = REPO_ROOT / rel
        if not p.exists():
            raise FileNotFoundError(f"Expected original input missing: {p}")
        manifest[rel] = sha256_file(p)
    TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
    snap_path = TMP_WORK_DIR / "original_sha256_before.json"
    if not snap_path.exists():
        snap_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def bh_explicit(p_values: np.ndarray | list[float] | pd.Series) -> np.ndarray:
    p = np.asarray(p_values, dtype=np.float64)
    n = p.size
    if n == 0:
        return np.array([], dtype=np.float64)
    out = np.full(n, np.nan, dtype=np.float64)
    valid = np.isfinite(p)
    if not np.any(valid):
        return out
    pv = p[valid]
    m = pv.size
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    adj = ranked * float(m) / np.arange(1, m + 1, dtype=np.float64)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    qv = np.empty(m, dtype=np.float64)
    qv[order] = adj
    out[valid] = qv
    return out


def bh_checked(p_values: np.ndarray | pd.Series) -> tuple[np.ndarray, np.ndarray, float]:
    p = np.asarray(p_values, dtype=np.float64)
    q_exp = bh_explicit(p)
    q_sm = np.full_like(q_exp, np.nan)
    valid = np.isfinite(p)
    if np.any(valid):
        q_sm[valid] = multipletests(p[valid], method="fdr_bh")[1]
    diff = np.nanmax(np.abs(q_exp - q_sm)) if np.any(valid) else 0.0
    if diff > 1e-12:
        raise AssertionError(f"BH explicit vs statsmodels discrepancy {diff:.3e} > 1e-12")
    return q_exp, q_sm, float(diff)


def mde80_from_se(se: float, residual_df: int) -> float:
    if not np.isfinite(se) or se <= 0 or residual_df < 1:
        return float("nan")
    factor = float(scipy_stats.t.ppf(0.975, residual_df) + scipy_stats.t.ppf(0.80, residual_df))
    return factor * float(se)


def mde80_from_null_sd(null_sd: float) -> float:
    if not np.isfinite(null_sd) or null_sd <= 0:
        return float("nan")
    return MDE80_NULL_FACTOR * float(null_sd)


def load_verified_mature_programs() -> dict[str, dict[str, list[str]]]:
    """Load and assert exact counts for P1m and P0m against derivation report."""
    p1m_df = pd.read_csv(RESULTS_DIR / "followup-mature-p1m-genes.csv")
    p0m_df = pd.read_csv(RESULTS_DIR / "followup-mature-p1m-p0m-genes.csv")

    assert len(p1m_df) == 1914, f"P1m mouse count {len(p1m_df)} != 1914"
    assert int((p1m_df["direction"] == "UP").sum()) == 1119
    assert int((p1m_df["direction"] == "DOWN").sum()) == 795

    p1m_hm = p1m_df[p1m_df["human_ortholog"].notna() & (p1m_df["human_ortholog"].astype(str).str.strip() != "")].copy()
    p1m_hm["human_ortholog"] = p1m_hm["human_ortholog"].astype(str).str.strip()
    assert len(p1m_hm) == 1724, f"P1m human mapped {len(p1m_hm)} != 1724"
    assert p1m_hm["human_ortholog"].nunique() == 1724
    p1m_up = sorted(p1m_hm.loc[p1m_hm["direction"] == "UP", "human_ortholog"].tolist())
    p1m_dn = sorted(p1m_hm.loc[p1m_hm["direction"] == "DOWN", "human_ortholog"].tolist())
    assert len(p1m_up) == 1013 and len(p1m_dn) == 711

    assert len(p0m_df) == 2598, f"P0m mouse count {len(p0m_df)} != 2598"
    assert int((p0m_df["direction"] == "UP").sum()) == 1474
    assert int((p0m_df["direction"] == "DOWN").sum()) == 1124

    p0m_hm = p0m_df[p0m_df["human_ortholog"].notna() & (p0m_df["human_ortholog"].astype(str).str.strip() != "")].copy()
    p0m_hm["human_ortholog"] = p0m_hm["human_ortholog"].astype(str).str.strip()
    assert len(p0m_hm) == 2049, f"P0m human mapped {len(p0m_hm)} != 2049"
    assert p0m_hm["human_ortholog"].nunique() == 2049
    p0m_up = sorted(p0m_hm.loc[p0m_hm["direction"] == "UP", "human_ortholog"].tolist())
    p0m_dn = sorted(p0m_hm.loc[p0m_hm["direction"] == "DOWN", "human_ortholog"].tolist())
    assert len(p0m_up) == 1157 and len(p0m_dn) == 892

    p1m_set = set(p1m_up) | set(p1m_dn)
    p0m_set = set(p0m_up) | set(p0m_dn)
    assert not p1m_set.issubset(p0m_set), "Expected P1m not to be a strict subset of P0m"
    print(
        f"[VERIFIED] P1m={len(p1m_set)} ({len(p1m_up)} UP, {len(p1m_dn)} DOWN); "
        f"P0m={len(p0m_set)} ({len(p0m_up)} UP, {len(p0m_dn)} DOWN); "
        f"union={len(p1m_set | p0m_set)}",
        flush=True,
    )
    return {
        "P1m": {"up": p1m_up, "down": p1m_dn},
        "P0m": {"up": p0m_up, "down": p0m_dn},
    }


def fast_random_signature_scores(
    matrix: np.ndarray,
    up: np.ndarray,
    down: np.ndarray,
    iterations: int = P1_ITERATIONS,
    seed: int = ORIG_NULL_SEED,
    verify_parity: bool = False,
) -> tuple[np.ndarray, int]:
    """Bit-for-bit parity with src.followup_muscle.random_signature_scores, precomputing pools."""
    rng = np.random.default_rng(seed)
    deciles = expression_deciles(matrix)
    available = np.ones(matrix.shape[0], dtype=bool)
    available[np.concatenate([up, down])] = False

    plan: list[list[tuple[np.ndarray, int]]] = []
    min_margin = 10 ** 9
    for program_genes in (up, down):
        vc_items = list(pd.Series(deciles[program_genes]).value_counts().items())
        dir_plan: list[tuple[np.ndarray, int]] = []
        for decile, count in vc_items:
            cnt = int(count)
            pool = np.flatnonzero((deciles == decile) & available)
            margin = len(pool) - cnt
            min_margin = min(min_margin, margin)
            if len(pool) < cnt:
                raise ValueError(
                    f"Insufficient expression-matched random-signature pool in decile {decile}: "
                    f"need {cnt}, have {len(pool)}"
                )
            dir_plan.append((pool, cnt))
        plan.append(dir_plan)

    n_up = int(len(up))
    n_dn = int(len(down))
    n_donors = matrix.shape[1]
    scores = np.empty((iterations, n_donors), dtype=np.float64)
    up_buf = np.empty(n_up, dtype=np.int64)
    dn_buf = np.empty(n_dn, dtype=np.int64)

    for it in range(iterations):
        pos = 0
        for pool, cnt in plan[0]:
            up_buf[pos : pos + cnt] = rng.choice(pool, size=cnt, replace=False)
            pos += cnt
        pos = 0
        for pool, cnt in plan[1]:
            dn_buf[pos : pos + cnt] = rng.choice(pool, size=cnt, replace=False)
            pos += cnt
        scores[it, :] = matrix[up_buf, :].mean(axis=0) - matrix[dn_buf, :].mean(axis=0)

    if verify_parity:
        ref = orig_random_signature_scores(matrix, up, down, iterations=min(5, iterations), seed=seed)
        np.testing.assert_allclose(scores[: ref.shape[0]], ref, rtol=0.0, atol=1e-15)

    return scores, int(min_margin)


def run_family_c_mature(
    expr_by_group: dict[str, pd.DataFrame],
    meta_by_group: dict[str, pd.DataFrame],
    program: dict[str, list[str]],
    age_col: str = "age_per_decade",
    candidate_covariates: tuple[str, ...] = ("sex_male", "log_lib"),
    min_donors: int = 4,
    min_cpm: float = 1.0,
    min_expr_frac: float = 0.5,
    min_residual_df: int = 2,
    iterations: int = P1_ITERATIONS,
    seed: int = ORIG_NULL_SEED,
    verify_first_parity: bool = True,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Run Family C scoring + expression-decile-matched null with explicit untestable reasons."""
    up_set = set(program["up"])
    down_set = set(program["down"])
    log2_thresh = float(np.log2(min_cpm + 1.0))

    rows: list[dict[str, Any]] = []
    skipped: dict[str, str] = {}
    verified_once = False

    for group, log2cpm_df in expr_by_group.items():
        meta = meta_by_group[group]
        common_donors = [d for d in meta.index if d in log2cpm_df.columns]
        if len(common_donors) < min_donors:
            skipped[group] = f"n_donors={len(common_donors)} < min_donors={min_donors}"
            continue
        mat_df = log2cpm_df[common_donors]
        m_meta = meta.loc[common_donors]
        assert float(m_meta["age_years"].min()) >= MIN_AGE, f"Under-30 donor leaked into {group}"
        min_req = max(3, int(np.ceil(min_expr_frac * len(common_donors))))
        keep_mask = (mat_df.to_numpy(dtype=np.float64) >= log2_thresh).sum(axis=1) >= min_req
        mat_sub = mat_df.loc[keep_mask]
        if mat_sub.shape[0] < 200:
            skipped[group] = f"n_expressed_pool={mat_sub.shape[0]} < 200"
            continue

        genes = mat_sub.index.astype(str).tolist()
        up_idx = np.array([i for i, g in enumerate(genes) if g in up_set], dtype=int)
        down_idx = np.array([i for i, g in enumerate(genes) if g in down_set], dtype=int)
        if len(up_idx) == 0 or len(down_idx) == 0:
            skipped[group] = f"detected up={len(up_idx)}, down={len(down_idx)}"
            continue

        try:
            design, names, omitted = select_design_columns(
                m_meta,
                age_col=age_col,
                candidate_covariates=candidate_covariates,
                min_residual_df=min_residual_df,
            )
        except ValueError as exc:
            skipped[group] = f"design_error: {exc}"
            continue

        mat = mat_sub.to_numpy(dtype=np.float64)
        up_score = mat[up_idx, :].mean(axis=0)
        down_score = mat[down_idx, :].mean(axis=0)
        comp_score = up_score - down_score
        ols = fit_ols(np.column_stack([up_score, down_score, comp_score]), design, effect_column=1)

        try:
            null_scores, min_decile_margin = fast_random_signature_scores(
                mat,
                up_idx,
                down_idx,
                iterations=iterations,
                seed=seed,
                verify_parity=(verify_first_parity and not verified_once),
            )
            verified_once = True
        except ValueError as exc:
            skipped[group] = f"matched_null_pool_error: {exc}"
            continue

        null_ols = fit_ols(null_scores.T, design, effect_column=1)
        null_effects = null_ols["effect"]
        p_emp = matched_null_effect_pvalue(float(ols["effect"][2]), null_effects)

        rows.append(
            {
                "group": str(group),
                "n_donors": int(len(common_donors)),
                "age_min": float(m_meta["age_years"].min()),
                "age_max": float(m_meta["age_years"].max()),
                "n_expressed_pool": int(len(genes)),
                "n_up_detected": int(len(up_idx)),
                "n_down_detected": int(len(down_idx)),
                "n_total_detected": int(len(up_idx) + len(down_idx)),
                "min_decile_margin": int(min_decile_margin),
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "up_effect": float(ols["effect"][0]),
                "up_pvalue": float(ols["pvalue"][0]),
                "down_effect": float(ols["effect"][1]),
                "down_pvalue": float(ols["pvalue"][1]),
                "composite_effect": float(ols["effect"][2]),
                "composite_se": float(ols["se"][2]),
                "composite_ci95_low": float(ols["ci95_low"][2]),
                "composite_ci95_high": float(ols["ci95_high"][2]),
                "composite_pvalue": float(ols["pvalue"][2]),
                "null_mean_effect": float(null_effects.mean()),
                "null_sd_effect": float(null_effects.std(ddof=1)),
                "null_iterations": int(iterations),
                "empirical_null_pvalue": float(p_emp),
                "residual_df": int(ols["residual_df"][2]),
            }
        )

    if not rows:
        return pd.DataFrame(), skipped
    res = pd.DataFrame(rows)
    res["family_c_size"] = int(len(res))
    q_comp, _, _ = bh_checked(res["composite_pvalue"])
    q_emp, _, _ = bh_checked(res["empirical_null_pvalue"])
    res["q_bh_composite"] = q_comp
    res["q_bh_empirical"] = q_emp
    return res.sort_values("empirical_null_pvalue", ascending=True).reset_index(drop=True), skipped


def run_family_b_mature(
    expr_by_group: dict[str, pd.DataFrame],
    meta_by_group: dict[str, pd.DataFrame],
    age_col: str = "age_per_decade",
    candidate_covariates: tuple[str, ...] = ("sex_male", "log_lib"),
    min_donors: int = 4,
    min_cpm: float = 1.0,
    min_expr_frac: float = 0.5,
    min_residual_df: int = 2,
    run_concordance: bool = True,
    concordance_perms: int = CONCORDANCE_PERMS,
    seed: int = ORIG_NULL_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run Family B donor-level gene OLS + FWL-verified exact permutation concordance with null stats."""
    log2_thresh = float(np.log2(min_cpm + 1.0))
    gene_frames: list[pd.DataFrame] = []
    group_summaries: list[dict[str, Any]] = []
    t_series_by_group: dict[str, pd.Series] = {}
    designs_by_group: dict[str, tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]] = {}

    for group, log2cpm_df in expr_by_group.items():
        meta = meta_by_group[group]
        common_donors = [d for d in meta.index if d in log2cpm_df.columns]
        if len(common_donors) < min_donors:
            group_summaries.append(
                {
                    "group": group,
                    "status": f"untestable_low_donors({len(common_donors)}<{min_donors})",
                    "n_donors": int(len(common_donors)),
                    "n_genes_tested": 0,
                }
            )
            continue
        mat_df = log2cpm_df[common_donors]
        m_meta = meta.loc[common_donors]
        assert float(m_meta["age_years"].min()) >= MIN_AGE, f"Under-30 donor leaked into {group}"
        min_req = max(3, int(np.ceil(min_expr_frac * len(common_donors))))
        keep_mask = (mat_df.to_numpy(dtype=np.float64) >= log2_thresh).sum(axis=1) >= min_req
        if int(keep_mask.sum()) == 0:
            group_summaries.append(
                {
                    "group": group,
                    "status": "untestable_no_expressed_genes",
                    "n_donors": int(len(common_donors)),
                    "n_genes_tested": 0,
                }
            )
            continue
        mat_sub = mat_df.loc[keep_mask]
        try:
            design, names, omitted = select_design_columns(
                m_meta,
                age_col=age_col,
                candidate_covariates=candidate_covariates,
                min_residual_df=min_residual_df,
            )
        except ValueError as exc:
            group_summaries.append(
                {
                    "group": group,
                    "status": f"untestable_design({exc})",
                    "n_donors": int(len(common_donors)),
                    "n_genes_tested": 0,
                }
            )
            continue

        outcomes = mat_sub.to_numpy(dtype=np.float64).T
        ols = fit_ols(outcomes, design, effect_column=1)
        t_stats = np.where(ols["se"] > 0, ols["effect"] / ols["se"], np.nan)
        t_series_by_group[group] = pd.Series(t_stats, index=mat_sub.index)
        designs_by_group[group] = (mat_sub, m_meta, design, names)

        df_val = int(ols["residual_df"][0])
        mde_factor = float(scipy_stats.t.ppf(0.975, df_val) + scipy_stats.t.ppf(0.80, df_val))
        mde_vec = mde_factor * ols["se"]

        q_within, _, _ = bh_checked(ols["pvalue"])
        df_g = pd.DataFrame(
            {
                "group": str(group),
                "gene": mat_sub.index.astype(str),
                "n_donors": int(len(common_donors)),
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "mean_log2cpm": outcomes.mean(axis=0),
                "effect": ols["effect"],
                "se": ols["se"],
                "mde80": mde_vec,
                "t_stat": t_stats,
                "pvalue": ols["pvalue"],
                "residual_df": ols["residual_df"].astype(int),
                "q_bh_within_group": q_within,
            }
        )
        gene_frames.append(df_g)
        min_idx = int(np.nanargmin(ols["pvalue"]))
        group_summaries.append(
            {
                "group": str(group),
                "status": "tested",
                "n_donors": int(len(common_donors)),
                "age_min": float(m_meta["age_years"].min()),
                "age_max": float(m_meta["age_years"].max()),
                "n_genes_tested": int(len(df_g)),
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "residual_df": df_val,
                "median_se": float(np.nanmedian(ols["se"])),
                "median_mde80": float(np.nanmedian(mde_vec)),
                "top_gene": str(mat_sub.index[min_idx]),
                "top_gene_effect": float(ols["effect"][min_idx]),
                "top_gene_se": float(ols["se"][min_idx]),
            }
        )

    if not gene_frames:
        return pd.DataFrame(), pd.DataFrame(group_summaries), pd.DataFrame()

    all_genes = pd.concat(gene_frames, ignore_index=True)
    all_genes["family_b_global_size"] = int(len(all_genes))
    q_glob, _, bh_diff = bh_checked(all_genes["pvalue"])
    all_genes["q_bh_global"] = q_glob

    summary_df = pd.DataFrame(group_summaries)
    within_hits = all_genes.groupby("group")["q_bh_within_group"].apply(lambda s: int((s < ALPHA).sum())).to_dict()
    global_hits = all_genes.groupby("group")["q_bh_global"].apply(lambda s: int((s < ALPHA).sum())).to_dict()
    min_p_map = all_genes.groupby("group")["pvalue"].min().to_dict()
    min_q_map = all_genes.groupby("group")["q_bh_global"].min().to_dict()
    summary_df["n_sig_within_q05"] = summary_df["group"].map(within_hits).fillna(0).astype(int)
    summary_df["n_sig_global_q05"] = summary_df["group"].map(global_hits).fillna(0).astype(int)
    summary_df["min_pvalue"] = summary_df["group"].map(min_p_map)
    summary_df["min_q_bh_global"] = summary_df["group"].map(min_q_map)
    summary_df["family_b_global_size"] = int(len(all_genes))
    summary_df["bh_max_diff"] = float(bh_diff)

    conc_df = pd.DataFrame()
    if run_concordance and len(t_series_by_group) >= 2 and concordance_perms > 0:
        conc_df = compute_concordance_fwl(
            t_series_by_group=t_series_by_group,
            designs_by_group=designs_by_group,
            age_col=age_col,
            concordance_perms=concordance_perms,
            seed=seed,
        )

    all_genes = all_genes.sort_values("pvalue", ascending=True).reset_index(drop=True)
    return all_genes, summary_df, conc_df


def compute_concordance_fwl(
    t_series_by_group: dict[str, pd.Series],
    designs_by_group: dict[str, tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]],
    age_col: str,
    concordance_perms: int,
    seed: int,
) -> pd.DataFrame:
    """Exact permutation concordance via FWL orthogonal projection (verified against fit_ols)."""
    groups = sorted(t_series_by_group.keys())
    common_genes = sorted(set.intersection(*(set(t_series_by_group[g].index) for g in groups)))
    if len(common_genes) < 50:
        return pd.DataFrame(
            [
                {
                    "n_cell_types": int(len(groups)),
                    "n_common_genes": int(len(common_genes)),
                    "status": "untestable_few_common_genes",
                }
            ]
        )

    t_mat = np.vstack([t_series_by_group[g].loc[common_genes].to_numpy(dtype=np.float64) for g in groups])
    t_mat = np.nan_to_num(t_mat, nan=0.0)
    pair_r = pairwise_pearson(t_mat)
    obs_mean_r = float(np.mean(pair_r)) if len(pair_r) else np.nan
    obs_median_r = float(np.median(pair_r)) if len(pair_r) else np.nan

    pair_spear = []
    pair_sign = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            rho, _ = scipy_stats.spearmanr(t_mat[i], t_mat[j])
            pair_spear.append(float(rho))
            pair_sign.append(float(np.mean(np.sign(t_mat[i]) == np.sign(t_mat[j]))))

    all_donors = sorted(set.union(*(set(designs_by_group[g][1].index) for g in groups)))
    donor_age: dict[str, float] = {}
    for g in groups:
        _, meta_g, _, _ = designs_by_group[g]
        for d, val in meta_g[age_col].items():
            donor_age[d] = float(val)
    age_vec = np.array([donor_age[d] for d in all_donors], dtype=np.float64)
    donor_to_pos = {d: idx for idx, d in enumerate(all_donors)}

    # Precompute FWL orthogonal bases Q_Z and residualized Y_perp for each group
    fwl_cache: dict[str, dict[str, Any]] = {}
    for g in groups:
        mat_sub, meta_g, design_g, _ = designs_by_group[g]
        Y = mat_sub.loc[common_genes].to_numpy(dtype=np.float64).T  # (n, G)
        n_g, p_g = design_g.shape
        df_g = n_g - p_g
        z_cols = [0] + list(range(2, p_g))
        Z = design_g[:, z_cols]
        Q, _ = np.linalg.qr(Z, mode="reduced")
        Y_perp = Y - Q @ (Q.T @ Y)
        ss_y = np.sum(Y_perp * Y_perp, axis=0)
        d_idx = np.array([donor_to_pos[d] for d in meta_g.index], dtype=np.int64)

        # Verify FWL reproduces observed t-stats exactly
        a_obs = design_g[:, 1]
        a_perp = a_obs - Q @ (Q.T @ a_obs)
        s_a = float(a_perp @ a_perp)
        dot_obs = a_perp @ Y_perp
        rss_obs = np.maximum(0.0, ss_y - (dot_obs * dot_obs) / s_a)
        se_obs = np.sqrt((rss_obs / df_g) / s_a)
        t_fwl = np.where(se_obs > 0, (dot_obs / s_a) / np.where(se_obs > 0, se_obs, 1.0), 0.0)
        t_ref = np.nan_to_num(t_series_by_group[g].loc[common_genes].to_numpy(dtype=np.float64), nan=0.0)
        np.testing.assert_allclose(t_fwl, t_ref, rtol=1e-7, atol=1e-7)

        fwl_cache[g] = {
            "Q": Q,
            "Y_perp": Y_perp,
            "ss_y": ss_y,
            "df": float(df_g),
            "d_idx": d_idx,
            "Y_raw": Y,
            "design_raw": design_g,
        }

    rng = np.random.default_rng(seed)
    perm_means: list[float] = []
    verified_perm0 = False

    for _ in range(concordance_perms):
        perm_ages = rng.permutation(age_vec)
        perm_t_rows = []
        ok = True
        for g in groups:
            c = fwl_cache[g]
            a = perm_ages[c["d_idx"]]
            Q = c["Q"]
            a_perp = a - Q @ (Q.T @ a)
            s_a = float(a_perp @ a_perp)
            a_norm_sq = float(a @ a)
            if s_a <= 1e-20 * max(a_norm_sq, 1.0):
                ok = False
                break
            dot = a_perp @ c["Y_perp"]
            rss = np.maximum(0.0, c["ss_y"] - (dot * dot) / s_a)
            se = np.sqrt((rss / c["df"]) / s_a)
            t_perm = np.where(se > 0, (dot / s_a) / np.where(se > 0, se, 1.0), 0.0)
            t_perm = np.nan_to_num(t_perm, nan=0.0, posinf=0.0, neginf=0.0)

            if not verified_perm0 and g == groups[0]:
                d_check = c["design_raw"].copy()
                d_check[:, 1] = a
                ols_ref = fit_ols(c["Y_raw"], d_check, effect_column=1)
                t_ref = np.where(
                    ols_ref["se"] > 0,
                    ols_ref["effect"] / np.where(ols_ref["se"] > 0, ols_ref["se"], 1.0),
                    0.0,
                )
                np.testing.assert_allclose(t_perm, np.nan_to_num(t_ref, nan=0.0), rtol=1e-7, atol=1e-7)
                verified_perm0 = True

            perm_t_rows.append(t_perm)
        if ok:
            perm_r = pairwise_pearson(np.vstack(perm_t_rows))
            perm_means.append(float(np.mean(perm_r)))

    perm_arr = np.asarray(perm_means, dtype=np.float64)
    perm_null_mean = float(np.mean(perm_arr)) if len(perm_arr) > 0 else np.nan
    perm_null_sd = float(np.std(perm_arr, ddof=1)) if len(perm_arr) > 1 else np.nan
    perm_p_one = float((1 + np.sum(perm_arr >= obs_mean_r)) / (len(perm_arr) + 1)) if len(perm_arr) > 0 else np.nan
    perm_p_two = float((1 + np.sum(np.abs(perm_arr) >= abs(obs_mean_r))) / (len(perm_arr) + 1)) if len(perm_arr) > 0 else np.nan

    return pd.DataFrame(
        [
            {
                "status": "tested",
                "n_cell_types": int(len(groups)),
                "n_pairs": int(len(pair_r)),
                "n_common_genes": int(len(common_genes)),
                "n_donors_union": int(len(all_donors)),
                "mean_pairwise_pearson_r": obs_mean_r,
                "median_pairwise_pearson_r": obs_median_r,
                "min_pairwise_pearson_r": float(np.min(pair_r)),
                "max_pairwise_pearson_r": float(np.max(pair_r)),
                "mean_pairwise_spearman_rho": float(np.mean(pair_spear)),
                "mean_sign_agreement": float(np.mean(pair_sign)),
                "n_valid_perms": int(len(perm_arr)),
                "perm_null_mean": perm_null_mean,
                "perm_null_sd": perm_null_sd,
                "mde80_null": mde80_from_null_sd(perm_null_sd),
                "perm_p_one_sided": perm_p_one,
                "perm_p_two_sided": perm_p_two,
                "cell_types": ";".join(groups),
            }
        ]
    )


def transport_class(comp_q: float | None, emp_q: float | None, testable: bool = True) -> str:
    if not testable:
        return "UNTESTABLE"
    c_pass = bool(pd.notna(comp_q) and float(comp_q) < ALPHA)
    e_pass = bool(pd.notna(emp_q) and float(emp_q) < ALPHA)
    if c_pass and e_pass:
        return "BOTH_PASS"
    if e_pass:
        return "EMPIRICAL_ONLY"
    if c_pass:
        return "COMPOSITE_ONLY"
    return "FAIL_BOTH"


def verdict_standard(
    orig_sig: bool | None,
    mature_sig: bool,
    orig_eff: float | None = None,
    mature_eff: float | None = None,
    testable: bool = True,
) -> tuple[str, str]:
    if not testable:
        return "UNTESTABLE", ""
    if orig_sig is None:
        return ("NEW" if mature_sig else "NULL"), "no_original_baseline"
    if orig_sig and mature_sig:
        if orig_eff is not None and mature_eff is not None and np.isfinite(orig_eff) and np.isfinite(mature_eff):
            if np.sign(orig_eff) != np.sign(mature_eff):
                return "DIED", "sign_flipped_significant"
        return "SURVIVED", ""
    if orig_sig and not mature_sig:
        return "DIED", ""
    if (not orig_sig) and mature_sig:
        return "NEW", ""
    return "NULL", ""


def anova_batch_age(donor_meta: pd.DataFrame, col: str, single_only: bool = True) -> dict[str, Any]:
    if col not in donor_meta.columns:
        return {"column": col, "present": False}
    vals = donor_meta[col].astype(str)
    valid = ~vals.str.contains(";") if single_only else pd.Series(True, index=donor_meta.index)
    sub = donor_meta.loc[valid]
    groups = [g["age_years"].to_numpy(dtype=np.float64) for _, g in sub.groupby(vals.loc[valid]) if len(g) >= 2]
    if len(groups) < 2:
        return {
            "column": col,
            "present": True,
            "n_single_donors": int(valid.sum()),
            "n_groups_ge2": int(len(groups)),
            "anova_r2": np.nan,
            "anova_pvalue": np.nan,
        }
    f_stat, p_val = scipy_stats.f_oneway(*groups)
    y = np.concatenate(groups)
    gm = float(y.mean())
    ss_tot = float(np.sum((y - gm) ** 2))
    ss_bet = float(sum(len(g) * (float(g.mean()) - gm) ** 2 for g in groups))
    return {
        "column": col,
        "present": True,
        "n_single_donors": int(valid.sum()),
        "n_donors_in_ge2": int(len(y)),
        "n_groups_ge2": int(len(groups)),
        "anova_r2": ss_bet / ss_tot if ss_tot > 0 else 0.0,
        "anova_f": float(f_stat),
        "anova_pvalue": float(p_val),
    }


def make_row(
    cohort: str,
    cell_type: str,
    endpoint_family: str,
    n_donors: int | float,
    n_excluded: int | float,
    original_p: float | None,
    original_q: float | None,
    mature_p: float | None,
    mature_q: float | None,
    mature_empirical_p: float | None,
    mature_empirical_q: float | None,
    verdict: str,
    notes: str,
    **extra: Any,
) -> dict[str, Any]:
    base: dict[str, Any] = {
        "cohort": cohort,
        "cell_type": sanitize_dashes(str(cell_type)),
        "endpoint_family": endpoint_family,
        "n_donors": n_donors,
        "n_donors_30plus_excluded": n_excluded,
        "original_p": original_p if original_p is not None else np.nan,
        "original_q": original_q if original_q is not None else np.nan,
        "mature_p": mature_p if mature_p is not None else np.nan,
        "mature_q": mature_q if mature_q is not None else np.nan,
        "mature_empirical_p": mature_empirical_p if mature_empirical_p is not None else np.nan,
        "mature_empirical_q": mature_empirical_q if mature_empirical_q is not None else np.nan,
        "verdict": verdict,
        "notes": sanitize_dashes(str(notes)),
    }
    for k, v in extra.items():
        base[k] = sanitize_dashes(v) if isinstance(v, str) else v
    return base


def order_cohort_df(rows: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    extra_cols = [c for c in df.columns if c not in REQUIRED_COHORT_COLS]
    df = df[REQUIRED_COHORT_COLS + extra_cols]
    dup_mask = df.duplicated(subset=["cohort", "cell_type", "endpoint_family"], keep=False)
    if bool(dup_mask.any()):
        dups = df.loc[dup_mask, ["cohort", "cell_type", "endpoint_family"]].to_dict("records")
        raise AssertionError(f"Duplicate (cohort, cell_type, endpoint_family) rows: {dups[:5]}")
    return df


def save_cohort_payload(cohort_key: str, payload: dict[str, Any]) -> None:
    TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
    with (TMP_WORK_DIR / f"payload_{cohort_key}.pkl").open("wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_cohort_payload(cohort_key: str) -> dict[str, Any]:
    with (TMP_WORK_DIR / f"payload_{cohort_key}.pkl").open("rb") as fh:
        return pickle.load(fh)


def check_disk_for_download(path: Path, expected_bytes: int) -> None:
    free_bytes = shutil.disk_usage(path).free
    reserve_bytes = int(MIN_FREE_DISK_GIB * (1024 ** 3))
    if free_bytes - expected_bytes < reserve_bytes:
        raise RuntimeError(
            f"Insufficient disk space at {path}: free={free_bytes / 1024**3:.2f} GiB, "
            f"download={expected_bytes / 1024**3:.2f} GiB, reserve={MIN_FREE_DISK_GIB:.1f} GiB"
        )


def download_h5ad_verified(url: str, dest: Path, expected_bytes: int) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size == expected_bytes:
        print(f"[DOWNLOAD] Reusing verified local file {dest} ({expected_bytes} bytes)", flush=True)
        return dest
    check_disk_for_download(dest.parent, expected_bytes)
    cmd = [
        "curl",
        "--fail",
        "--location",
        "--continue-at",
        "-",
        "--retry",
        "5",
        "--retry-delay",
        "5",
        "--retry-all-errors",
        "--connect-timeout",
        "30",
        "--speed-limit",
        "10240",
        "--speed-time",
        "120",
        "--output",
        str(dest),
        url,
    ]
    print(f"[DOWNLOAD] Fetching {url} -> {dest} ({expected_bytes / 1024**3:.2f} GiB)", flush=True)
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"curl failed with exit code {proc.returncode} for {url}")
    actual_bytes = dest.stat().st_size
    if actual_bytes != expected_bytes:
        raise RuntimeError(f"Downloaded size mismatch for {dest}: got {actual_bytes}, expected {expected_bytes}")
    print(f"[DOWNLOAD OK] Verified exact size {actual_bytes} bytes", flush=True)
    return dest


def stream_hcata_mature_from_h5ad(
    study: str,
    h5_path: Path,
    expected_n_obs: int,
    expected_total_nnz: int,
    expected_raw_source: str,
    expected_orig_donors_csv: Path,
    min_donors_extract: int = 4,
) -> dict[str, Any]:
    """Audit full H5AD against original metadata, then stream strictly 30-plus raw integer pseudobulks."""
    orig_donors_df = pd.read_csv(expected_orig_donors_csv)
    orig_donors_df["donor_id"] = orig_donors_df["donor_id"].astype(str)

    with h5py.File(h5_path, "r") as f:
        obs = read_elem(f["obs"])
        if "raw" in f and "var" in f["raw"]:
            raw_grp = f["raw/X"]
            var = read_elem(f["raw/var"])
            raw_source = "raw/X"
        else:
            raw_grp = f["X"]
            var = read_elem(f["var"])
            raw_source = "X"

        shape = tuple(int(x) for x in raw_grp.attrs["shape"])
        n_obs, n_vars = shape
        assert n_obs == expected_n_obs, f"{study}: n_obs {n_obs} != expected {expected_n_obs}"
        assert raw_source == expected_raw_source, f"{study}: raw_source {raw_source} != {expected_raw_source}"

        data_ds = raw_grp["data"]
        total_nnz = int(data_ds.shape[0])
        assert total_nnz == expected_total_nnz, f"{study}: total_nnz {total_nnz} != expected {expected_total_nnz}"

        # Audit full-object normal exact-age donors across ALL assays (to reconcile audit notes e.g. Guilliams)
        full_normal_exact = obs[
            (obs["disease"].astype(str) == "normal")
            & (obs["development_stage"].astype(str).str.endswith("-year-old stage"))
        ].copy()
        full_normal_exact["_age"] = full_normal_exact["development_stage"].astype(str).str.split("-").str[0].astype(float)
        full_obj_donors = (
            full_normal_exact.groupby("donor_id", observed=True)
            .agg(
                age_years=("_age", "first"),
                assays=("assay", lambda s: ";".join(sorted(set(s.astype(str))))),
                suspensions=("suspension_type", lambda s: ";".join(sorted(set(s.astype(str))))),
            )
            .reset_index()
        )
        full_obj_under30 = full_obj_donors[full_obj_donors["age_years"] < MIN_AGE].to_dict("records")

        # Apply exact original study cell mask
        if study == "guilliams":
            assay_target = "10x 3' v3"
        else:
            assay_target = "10x 3' v2"

        mask_all = (
            (obs["disease"].astype(str) == "normal")
            & (obs["assay"].astype(str) == assay_target)
            & (obs["suspension_type"].astype(str) == "cell")
            & (obs["development_stage"].astype(str).str.endswith("-year-old stage"))
        )
        sub_all = obs.loc[mask_all].copy()
        sub_all["_donor"] = sub_all["donor_id"].astype(str)
        sub_all["_age"] = sub_all["development_stage"].astype(str).str.split("-").str[0].astype(float)
        sub_all["_sex"] = sub_all["sex"].astype(str)

        all_donor_rows = []
        for did, dsub in sub_all.groupby("_donor", observed=True):
            age_vals = sorted(dsub["_age"].unique().tolist())
            sex_vals = sorted(dsub["_sex"].unique().tolist())
            assert len(age_vals) == 1 and len(sex_vals) == 1
            rec: dict[str, Any] = {
                "donor_id": did,
                "age_years": float(age_vals[0]),
                "age_per_decade": float(age_vals[0]) / 10.0,
                "sex": sex_vals[0],
                "sex_male": int(sex_vals[0] == "male"),
                "n_cells": int(len(dsub)),
            }
            for extra in (
                "reported_diseases",
                "diet",
                "suspension_enrichment_factors",
                "typeSample",
                "self_reported_ethnicity",
                "pool_number",
                "Processing_Cohort",
            ):
                if extra in dsub.columns:
                    rec[extra] = ";".join(sorted(set(dsub[extra].astype(str))))
            all_donor_rows.append(rec)

        donor_meta_all = pd.DataFrame(all_donor_rows).sort_values(["age_years", "donor_id"]).set_index("donor_id")
        assert set(donor_meta_all.index) == set(orig_donors_df["donor_id"]), (
            f"{study}: donor ID set mismatch vs {expected_orig_donors_csv}"
        )
        for row in orig_donors_df.itertuples(index=False):
            did = str(row.donor_id)
            assert float(donor_meta_all.loc[did, "age_years"]) == float(row.age_years)
            assert int(donor_meta_all.loc[did, "n_cells"]) == int(row.n_cells)

        # Restrict strictly to age >= 30.0 BEFORE cell counting, pseudobulk aggregation, and detection filtering
        donor_meta_30 = donor_meta_all[donor_meta_all["age_years"] >= MIN_AGE].copy()
        mature_donor_set = set(donor_meta_30.index.tolist())
        obs_donors = obs["donor_id"].astype(str).to_numpy()
        keep_30_bool = mask_all.to_numpy() & np.isin(obs_donors, list(mature_donor_set))
        assert int(keep_30_bool.sum()) == int(donor_meta_30["n_cells"].sum())

        sub30 = obs.loc[keep_30_bool]
        ct_counts_30 = sub30.groupby(["donor_id", "cell_type"], observed=True).size().unstack(fill_value=0)
        ct_counts_30.index = ct_counts_30.index.astype(str)
        ct_counts_30.columns = ct_counts_30.columns.astype(str)
        ct_counts_30 = ct_counts_30.reindex(donor_meta_30.index, fill_value=0)

        cd45_bool_30 = None
        cd45_counts_30 = None
        if study == "guilliams":
            cd45_bool_30 = keep_30_bool & (obs["suspension_enrichment_factors"].astype(str).to_numpy() == "CD45+")
            cd45_counts_30 = obs.loc[cd45_bool_30].groupby("donor_id", observed=True).size()
            cd45_counts_30.index = cd45_counts_30.index.astype(str)
            cd45_counts_30 = cd45_counts_30.reindex(donor_meta_30.index, fill_value=0).astype(int)

        elig_cts = [
            ct
            for ct in sorted(ct_counts_30.columns)
            if int((ct_counts_30[ct] >= MIN_CELLS_PER_DONOR).sum()) >= min_donors_extract
        ]
        groups_to_extract = list(elig_cts) + ["All_Cells"]
        if study == "guilliams":
            groups_to_extract.append("All_CD45pos_Cells")

        gene_names = (
            var["feature_name"].astype(str).to_numpy()
            if "feature_name" in var.columns
            else var.index.astype(str).to_numpy()
        )
        ct_arr = obs["cell_type"].astype(str).to_numpy()
        unique_donors = list(donor_meta_30.index)
        donor_idx_map = {d: i for i, d in enumerate(unique_donors)}
        donor_code_arr = np.array([donor_idx_map.get(d, -1) for d in obs_donors], dtype=np.int32)

        indptr_ds = raw_grp["indptr"]
        indices_ds = raw_grp["indices"]

        gate_30 = {"nonfinite": 0, "negative": 0, "noninteger": 0, "total_nnz_30": 0, "cells_30": 0}
        raw_pb_by_group: dict[str, pd.DataFrame] = {}

        for b0 in range(0, len(groups_to_extract), MAX_GROUP_BATCH):
            batch_groups = groups_to_extract[b0 : b0 + MAX_GROUP_BATCH]
            print(f"[{study}] Streaming 30-plus batch {b0 // MAX_GROUP_BATCH + 1}: {batch_groups}", flush=True)
            acc = {g: np.zeros((len(unique_donors), n_vars), dtype=np.float64) for g in batch_groups}

            for r0 in range(0, n_obs, CHUNK_ROWS):
                r1 = min(r0 + CHUNK_ROWS, n_obs)
                cmask = keep_30_bool[r0:r1]
                if not np.any(cmask):
                    continue
                ptrs = indptr_ds[r0 : r1 + 1]
                p0 = int(ptrs[0])
                p1 = int(ptrs[-1])
                mat = sp.csr_matrix(
                    (data_ds[p0:p1].astype(np.float64, copy=False), indices_ds[p0:p1], ptrs - p0),
                    shape=(r1 - r0, n_vars),
                )
                sub_rows = np.flatnonzero(cmask)
                mat_sub = mat[sub_rows]
                if b0 == 0:
                    g_chk = count_gate(mat_sub.data)
                    gate_30["nonfinite"] += g_chk["nonfinite"]
                    gate_30["negative"] += g_chk["negative"]
                    gate_30["noninteger"] += g_chk["noninteger"]
                    gate_30["total_nnz_30"] += int(len(mat_sub.data))
                    gate_30["cells_30"] += int(mat_sub.shape[0])
                    if g_chk["nonfinite"] or g_chk["negative"] or g_chk["noninteger"]:
                        raise ValueError(f"{study} 30-plus count_gate failed in {r0}:{r1}: {g_chk}")

                d_sub = donor_code_arr[r0:r1][sub_rows]
                c_sub = ct_arr[r0:r1][sub_rows]

                for grp in batch_groups:
                    if grp == "All_Cells":
                        rows_g = np.arange(len(sub_rows))
                    elif grp == "All_CD45pos_Cells":
                        rows_g = np.flatnonzero(cd45_bool_30[r0:r1][sub_rows])
                    else:
                        rows_g = np.flatnonzero(c_sub == grp)
                    if len(rows_g) == 0:
                        continue
                    dg = d_sub[rows_g]
                    mg = mat_sub[rows_g]
                    uniq_d, inv_d = np.unique(dg, return_inverse=True)
                    ind = sp.csr_matrix(
                        (np.ones(len(rows_g), dtype=np.float64), (inv_d, np.arange(len(rows_g)))),
                        shape=(len(uniq_d), len(rows_g)),
                    )
                    acc[grp][uniq_d, :] += (ind @ mg).toarray()

            for grp in batch_groups:
                if grp == "All_Cells":
                    elig_donors = unique_donors
                elif grp == "All_CD45pos_Cells":
                    elig_donors = [d for d in unique_donors if int(cd45_counts_30.loc[d]) >= MIN_CELLS_PER_DONOR]
                else:
                    elig_donors = [d for d in unique_donors if int(ct_counts_30.loc[d, grp]) >= MIN_CELLS_PER_DONOR]
                d_idx = [donor_idx_map[d] for d in elig_donors]
                raw_pb_by_group[grp] = pd.DataFrame(acc[grp][d_idx, :].T, index=gene_names, columns=elig_donors)
            del acc
            gc.collect()

    return {
        "study": study,
        "n_obs_total": int(n_obs),
        "total_nnz_file": int(total_nnz),
        "raw_source": raw_source,
        "gate_30": gate_30,
        "n_donors_orig": int(len(donor_meta_all)),
        "n_donors_30": int(len(donor_meta_30)),
        "n_donors_excluded": int(len(donor_meta_all) - len(donor_meta_30)),
        "full_object_donors_total": int(len(full_obj_donors)),
        "full_object_under30_donors": full_obj_under30,
        "donor_meta_30": donor_meta_30,
        "cell_counts_30": ct_counts_30,
        "cd45_counts_30": cd45_counts_30,
        "raw_pb_by_group": raw_pb_by_group,
    }


# =============================================================================
# COHORT 1: GAO HEART (GSE290367)
# =============================================================================

def run_cohort_gao(programs: dict[str, dict[str, list[str]]]) -> None:
    cohort_name = "Gao_Heart"
    t0_all = time.monotonic()
    s0_all = utc_now()
    wait_for_memory("Gao_Heart stream")

    print("\n================ COHORT 1: GAO HEART (GSE290367) ================", flush=True)
    gao_cache_pkl = TMP_WORK_DIR / "gao_30plus_raw_pb.pkl"

    if gao_cache_pkl.exists():
        print(f"[Gao] Loading verified 30-plus raw pseudobulk cache {gao_cache_pkl}", flush=True)
        with gao_cache_pkl.open("rb") as fh:
            cached = pickle.load(fh)
        raw_pb_by_group = cached["raw_pb_by_group"]
        cell_counts_df = cached["cell_counts_df"]
        donor_meta_30 = cached["donor_meta_30"]
        n_donors_orig = cached["n_donors_orig"]
        gate_totals = cached["gate_totals"]
        assert float(donor_meta_30["age_years"].min()) >= MIN_AGE
    else:
        t0_stream = time.monotonic()
        s0_stream = utc_now()
        with h5py.File(GAO_H5AD, "r") as handle:
            obs = read_elem(handle["obs"])
            var = read_elem(handle["var"])
        gene_names = var.index.astype(str).to_numpy()
        audit_df = gao_build_donor_audit(obs)
        donor_meta_all = gao_build_donor_meta_table(audit_df)
        n_donors_orig = int(len(donor_meta_all))
        assert n_donors_orig == 117, f"Gao original primary donors {n_donors_orig} != 117"

        donor_meta_30 = donor_meta_all[donor_meta_all["age_years"] >= MIN_AGE].copy()
        assert len(donor_meta_30) == 109, f"Gao 30-plus donors {len(donor_meta_30)} != 109"
        inc30 = set(donor_meta_30.index.tolist())
        obs_donors = obs["donor_id"].astype(str).to_numpy()
        obs_cts = obs["final_cell_type"].astype(str).to_numpy()
        keep_mask = np.isin(obs_donors, list(inc30))
        print(f"[Gao] Streaming {len(inc30)} mature donors, {int(keep_mask.sum())} nuclei...", flush=True)

        raw_pb_by_group, cell_counts_df, gate_totals = gao_stream_pseudobulk_counts(
            GAO_H5AD, keep_mask, obs_donors, obs_cts, gene_names
        )
        del obs, var, obs_donors, obs_cts, keep_mask
        gc.collect()
        TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
        tmp_p = gao_cache_pkl.with_suffix(".tmp")
        with tmp_p.open("wb") as fh:
            pickle.dump(
                {
                    "raw_pb_by_group": raw_pb_by_group,
                    "cell_counts_df": cell_counts_df,
                    "donor_meta_30": donor_meta_30,
                    "n_donors_orig": n_donors_orig,
                    "gate_totals": gate_totals,
                },
                fh,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        tmp_p.replace(gao_cache_pkl)
        append_runtime(
            cohort_name,
            "stream_h5ad_30plus",
            s0_stream,
            t0_stream,
            "ok",
            f"109 mature donors, {gate_totals['cells_checked']} nuclei, noninteger={gate_totals['noninteger']}",
        )

    n_donors_30 = int(len(donor_meta_30))
    n_excluded_cohort = n_donors_orig - n_donors_30
    celltype_names = sorted([c for c in cell_counts_df.columns if c != "Total_Cells"])
    assert len(celltype_names) == 13
    active_studies = set(donor_meta_30["study"].tolist())
    study_cols = [col for st, col in GAO_STUDY_DUMMY_MAP.items() if st in active_studies]
    cov_a = ("sex_male", *study_cols)
    cov_bc = ("sex_male", "log_lib", *study_cols)

    study_anova = anova_batch_age(donor_meta_30, "study", single_only=False)

    # Load original results for comparison
    orig_frac = pd.read_csv(RESULTS_DIR / "followup-organ-17-heart-gao-fractions.csv").set_index("cell_type")
    orig_gsum = pd.read_csv(RESULTS_DIR / "followup-organ-17-heart-gao-genes-summary.csv").set_index("group")
    orig_top = pd.read_csv(RESULTS_DIR / "followup-organ-17-heart-gao-top-genes.csv")
    orig_conc = pd.read_csv(RESULTS_DIR / "followup-organ-17-heart-gao-concordance.csv").iloc[0]
    orig_p1 = pd.read_csv(RESULTS_DIR / "followup-organ-17-heart-gao-p1-signature.csv").set_index("group")
    orig_b5 = pd.read_csv(RESULTS_DIR / "followup-11pt-block5-p0-transport.csv")
    orig_p0 = orig_b5[
        (orig_b5["cohort"] == "Gao_Heart") & (orig_b5["arm"] == "primary") & (orig_b5["gene_set"] == "P0_Raw257")
    ].set_index("unit")

    # Family A: CLR fractions
    t0_a = time.monotonic()
    s0_a = utc_now()
    fam_a = run_family_a_clr(
        cell_counts_df.loc[donor_meta_30.index, celltype_names],
        donor_meta_30,
        age_col="age_per_decade",
        candidate_covariates=cov_a,
        min_donors_present=3,
        min_total_cells=MIN_CELLS_PER_DONOR,
        min_residual_df=2,
    )
    q_a_exp, q_a_sm, bh_diff_a = bh_checked(fam_a["pvalue"])
    fam_a["q_bh_family_a"] = q_a_exp
    append_runtime(cohort_name, "family_a_fractions", s0_a, t0_a, "ok", f"n_tested={len(fam_a)}")

    # Finalize pseudobulks strictly on 30-plus donors
    expr_ct: dict[str, pd.DataFrame] = {}
    meta_ct: dict[str, pd.DataFrame] = {}
    for grp_name, raw_df in raw_pb_by_group.items():
        cols = [d for d in donor_meta_30.index if d in raw_df.columns]
        if len(cols) < 4:
            continue
        c_ser = (
            cell_counts_df.loc[cols, "Total_Cells"]
            if grp_name == "All_LV_Nuclei"
            else cell_counts_df.loc[cols, grp_name]
        )
        l2, gm, _ = finalize_pseudobulk_group(raw_df[cols], donor_meta_30.loc[cols], c_ser)
        expr_ct[grp_name] = l2
        meta_ct[grp_name] = gm

    # Family B: Genes + Concordance
    wait_for_memory("Gao_Heart Family B")
    t0_b = time.monotonic()
    s0_b = utc_now()
    expr_only_ct = {k: v for k, v in expr_ct.items() if k != "All_LV_Nuclei"}
    meta_only_ct = {k: v for k, v in meta_ct.items() if k != "All_LV_Nuclei"}
    genes_ct, sum_ct, conc_df = run_family_b_mature(
        expr_only_ct,
        meta_only_ct,
        age_col="age_per_decade",
        candidate_covariates=cov_bc,
        min_donors=4,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        run_concordance=True,
        concordance_perms=CONCORDANCE_PERMS,
        seed=ORIG_NULL_SEED,
    )
    genes_all_lv, sum_all_lv, _ = run_family_b_mature(
        {"All_LV_Nuclei": expr_ct["All_LV_Nuclei"]},
        {"All_LV_Nuclei": meta_ct["All_LV_Nuclei"]},
        age_col="age_per_decade",
        candidate_covariates=cov_bc,
        min_donors=4,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        run_concordance=False,
        seed=ORIG_NULL_SEED,
    )
    genes_combined = pd.concat([genes_ct, genes_all_lv], ignore_index=True)
    q_all14_exp, q_all14_sm, bh_diff_all14 = bh_checked(genes_combined["pvalue"])
    genes_combined["q_bh_all14"] = q_all14_exp
    genes_combined["q_bh_celltypes13"] = np.nan
    ct_mask = genes_combined["group"] != "All_LV_Nuclei"
    q_ct13_exp, q_ct13_sm, bh_diff_ct13 = bh_checked(genes_combined.loc[ct_mask, "pvalue"])
    genes_combined.loc[ct_mask, "q_bh_celltypes13"] = q_ct13_exp
    sum_b_all = pd.concat([sum_ct, sum_all_lv], ignore_index=True)
    all14_hits_map = genes_combined.groupby("group")["q_bh_all14"].apply(lambda s: int((s < ALPHA).sum())).to_dict()
    sum_b_all["n_sig_q05_all14"] = sum_b_all["group"].map(all14_hits_map).fillna(0).astype(int)
    append_runtime(
        cohort_name,
        "family_b_genes_concordance",
        s0_b,
        t0_b,
        "ok",
        f"ct13_tests={int(ct_mask.sum())}, ct13_sig={int((q_ct13_exp < ALPHA).sum())}, all14_sig={int((q_all14_exp < ALPHA).sum())}",
    )

    # Family C: P1m and P0m transport (14 groups, matching original family_c_size=14)
    wait_for_memory("Gao_Heart Family C")
    t0_c = time.monotonic()
    s0_c = utc_now()
    fam_c_p1m, skip_p1m = run_family_c_mature(
        expr_ct,
        meta_ct,
        program=programs["P1m"],
        candidate_covariates=cov_bc,
        min_donors=4,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=ORIG_NULL_SEED,
    )
    fam_c_p0m, skip_p0m = run_family_c_mature(
        expr_ct,
        meta_ct,
        program=programs["P0m"],
        candidate_covariates=cov_bc,
        min_donors=4,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=ORIG_NULL_SEED,
    )
    append_runtime(
        cohort_name,
        "family_c_transport",
        s0_c,
        t0_c,
        "ok",
        f"P1m_tested={len(fam_c_p1m)}, P0m_tested={len(fam_c_p0m)}, seed={ORIG_NULL_SEED}, draws={P1_ITERATIONS}",
    )

    # Assemble cohort CSV rows & BH check rows
    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []

    # 1. Fractions rows
    for idx_r, r in fam_a.iterrows():
        ct = str(r["cell_type"])
        orow = orig_frac.loc[ct]
        orig_sig = bool(float(orow["q_bh_family_a"]) < ALPHA)
        mat_sig = bool(float(r["q_bh_family_a"]) < ALPHA)
        verd, vnote = verdict_standard(orig_sig, mat_sig, float(orow["effect_clr"]), float(r["effect_clr"]))
        mde = mde80_from_se(float(r["se"]), int(r["residual_df"]))
        out_rows.append(
            make_row(
                cohort_name,
                ct,
                "fractions_clr",
                int(r["n_donors"]),
                int(orow["n_donors"]) - int(r["n_donors"]),
                float(orow["pvalue"]),
                float(orow["q_bh_family_a"]),
                float(r["pvalue"]),
                float(r["q_bh_family_a"]),
                None,
                None,
                verd,
                f"study_adjusted;n_present={int(r['n_donors_present'])};{vnote}".strip(";"),
                arm="pooled_7study_primary",
                residual_df=int(r["residual_df"]),
                original_effect=float(orow["effect_clr"]),
                mature_effect=float(r["effect_clr"]),
                mature_se=float(r["se"]),
                null_sd=np.nan,
                mde80=mde,
                mde80_method="t_quantile_se",
                covariates=str(r["covariates"]),
                omitted_covariates=str(r["omitted_covariates"]),
            )
        )
        bh_rows.append(
            {
                "cohort": cohort_name,
                "family": "gao_fractions_clr_13",
                "unit": ct,
                "family_size": int(len(fam_a)),
                "raw_p": float(r["pvalue"]),
                "q_bh_explicit": float(q_a_exp[idx_r]),
                "q_bh_statsmodels": float(q_a_sm[idx_r]),
                "abs_diff": float(abs(q_a_exp[idx_r] - q_a_sm[idx_r])),
                "family_max_abs_diff": bh_diff_a,
            }
        )

    # 2. Per-cell-type gene summaries + individual original/mature hit rows
    for r in sum_b_all.itertuples(index=False):
        grp = str(r.group)
        is_pool = grp == "All_LV_Nuclei"
        orow = orig_gsum.loc[grp]
        n_ex = int(orow["n_donors"]) - int(r.n_donors)
        if not is_pool:
            orig_nsig = int(orow["n_sig_global_q05"])
            mat_nsig = int(r.n_sig_global_q05)
            orig_minq = float(orow["min_q_bh_global"])
            mat_minq = float(r.min_q_bh_global)
            fam_label = "gao_genes_ct13"
            q_col = "q_bh_celltypes13"
        else:
            orig_nsig = int(orow["n_sig_q05_all14"])
            mat_nsig = int(r.n_sig_q05_all14)
            orig_minq = float(orig_top[orig_top["group"] == "All_LV_Nuclei"]["q_bh_all14"].min())
            mat_minq = float(genes_combined[genes_combined["group"] == "All_LV_Nuclei"]["q_bh_all14"].min())
            fam_label = "gao_genes_all14"
            q_col = "q_bh_all14"
        verd, _ = verdict_standard(orig_nsig > 0, mat_nsig > 0)
        out_rows.append(
            make_row(
                cohort_name,
                grp,
                "genes_celltype",
                int(r.n_donors),
                n_ex,
                float(orow["min_pvalue"]),
                orig_minq,
                float(r.min_pvalue),
                mat_minq,
                None,
                None,
                verd,
                f"n_genes={int(r.n_genes_tested)};orig_sig={orig_nsig};mature_sig={mat_nsig};top={r.top_gene};all14_sig={int(r.n_sig_q05_all14)}",
                arm="pooled_7study_primary",
                residual_df=int(r.residual_df),
                original_effect=np.nan,
                mature_effect=float(r.top_gene_effect),
                mature_se=float(r.median_se),
                null_sd=np.nan,
                mde80=float(r.median_mde80),
                mde80_method="median_t_quantile_se",
                n_genes_tested=int(r.n_genes_tested),
                n_sig_q05=mat_nsig,
                covariates=str(r.covariates),
                omitted_covariates=str(r.omitted_covariates),
            )
        )
        # Record top (min-p) gene in bh-check
        g_sub = genes_combined[genes_combined["group"] == grp].sort_values("pvalue").iloc[0]
        p_val = float(g_sub["pvalue"])
        if not is_pool:
            idx_in_ct = int(np.flatnonzero((genes_combined.loc[ct_mask, "group"] == grp) & (genes_combined.loc[ct_mask, "gene"] == g_sub["gene"]))[0])
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": "gao_genes_ct13",
                    "unit": f"{grp}::{g_sub['gene']}",
                    "family_size": int(ct_mask.sum()),
                    "raw_p": p_val,
                    "q_bh_explicit": float(q_ct13_exp[idx_in_ct]),
                    "q_bh_statsmodels": float(q_ct13_sm[idx_in_ct]),
                    "abs_diff": float(abs(q_ct13_exp[idx_in_ct] - q_ct13_sm[idx_in_ct])),
                    "family_max_abs_diff": bh_diff_ct13,
                }
            )
        idx_all = int(np.flatnonzero((genes_combined["group"] == grp) & (genes_combined["gene"] == g_sub["gene"]))[0])
        bh_rows.append(
            {
                "cohort": cohort_name,
                "family": "gao_genes_all14",
                "unit": f"{grp}::{g_sub['gene']}",
                "family_size": int(len(genes_combined)),
                "raw_p": p_val,
                "q_bh_explicit": float(q_all14_exp[idx_all]),
                "q_bh_statsmodels": float(q_all14_sm[idx_all]),
                "abs_diff": float(abs(q_all14_exp[idx_all] - q_all14_sm[idx_all])),
                "family_max_abs_diff": bh_diff_all14,
            }
        )

    # Individual gene hit tracking: all 12 original all14 hits (which includes all 8 ct13 hits) + all mature hits
    orig_sig_all14 = orig_top[orig_top["q_bh_all14"] < ALPHA].copy()
    orig_sig_map = {(str(r.group), str(r.gene)): r for r in orig_sig_all14.itertuples(index=False)}
    mat_sig_mask = (genes_combined["q_bh_celltypes13"] < ALPHA) | (genes_combined["q_bh_all14"] < ALPHA)
    mat_sig_pairs = set(zip(genes_combined.loc[mat_sig_mask, "group"], genes_combined.loc[mat_sig_mask, "gene"]))
    track_pairs = sorted(set(orig_sig_map.keys()) | mat_sig_pairs)

    genes_idx = genes_combined.set_index(["group", "gene"])
    ct_genes_pos = {
        (str(g), str(gn)): i
        for i, (g, gn) in enumerate(zip(genes_combined.loc[ct_mask, "group"], genes_combined.loc[ct_mask, "gene"]))
    }
    all_genes_pos = {
        (str(g), str(gn)): i for i, (g, gn) in enumerate(zip(genes_combined["group"], genes_combined["gene"]))
    }

    for grp, gene in track_pairs:
        is_pool = grp == "All_LV_Nuclei"
        orig_rec = orig_sig_map.get((grp, gene))
        if (grp, gene) in genes_idx.index:
            mrec = genes_idx.loc[(grp, gene)]
            mat_p = float(mrec["pvalue"])
            mat_q = float(mrec["q_bh_all14"]) if is_pool else float(mrec["q_bh_celltypes13"])
            mat_q14 = float(mrec["q_bh_all14"])
            mat_eff = float(mrec["effect"])
            mat_se = float(mrec["se"])
            mat_df = int(mrec["residual_df"])
            mat_n = int(mrec["n_donors"])
            mat_mde = float(mrec["mde80"])
            testable = True
        else:
            mat_p = np.nan
            mat_q = np.nan
            mat_q14 = np.nan
            mat_eff = np.nan
            mat_se = np.nan
            mat_df = np.nan
            mat_n = int(sum_b_all.set_index("group").loc[grp, "n_donors"])
            mat_mde = np.nan
            testable = False

        if orig_rec is not None:
            orig_p = float(orig_rec.pvalue)
            orig_q = float(orig_rec.q_bh_all14) if is_pool else float(orig_rec.q_bh_celltypes13)
            orig_sig_ct = bool(pd.notna(orig_rec.q_bh_celltypes13) and float(orig_rec.q_bh_celltypes13) < ALPHA)
            orig_sig_14 = bool(float(orig_rec.q_bh_all14) < ALPHA)
            orig_eff = float(orig_rec.effect)
            n_ex = int(orig_rec.n_donors) - mat_n
        else:
            orig_p = np.nan
            orig_q = np.nan
            orig_sig_ct = False
            orig_sig_14 = False
            orig_eff = np.nan
            n_ex = int(orig_gsum.loc[grp, "n_donors"]) - mat_n

        orig_sig_use = orig_sig_14 if is_pool else orig_sig_ct
        mat_sig_use = bool(testable and mat_q < ALPHA)
        if not testable:
            verd, vnote = "UNTESTABLE", "gene_failed_mature_50pct_detection_filter"
        else:
            verd, vnote = verdict_standard(orig_sig_use, mat_sig_use, orig_eff, mat_eff)
        note_str = (
            f"orig_ct13_sig={orig_sig_ct};orig_all14_sig={orig_sig_14};"
            f"mature_q_all14={mat_q14:.6g}" if np.isfinite(mat_q14) else f"orig_ct13_sig={orig_sig_ct};orig_all14_sig={orig_sig_14};unfiltered_out"
        )
        if vnote:
            note_str += f";{vnote}"
        out_rows.append(
            make_row(
                cohort_name,
                f"{grp}::{gene}",
                "genes_hit",
                mat_n,
                n_ex,
                orig_p,
                orig_q,
                mat_p,
                mat_q,
                None,
                None,
                verd,
                note_str,
                arm="pooled_7study_primary",
                residual_df=mat_df,
                original_effect=orig_eff,
                mature_effect=mat_eff,
                mature_se=mat_se,
                null_sd=np.nan,
                mde80=mat_mde,
                mde80_method="t_quantile_se",
            )
        )
        if testable:
            if not is_pool:
                pos_ct = ct_genes_pos[(grp, gene)]
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": "gao_genes_ct13",
                        "unit": f"{grp}::{gene}",
                        "family_size": int(ct_mask.sum()),
                        "raw_p": mat_p,
                        "q_bh_explicit": float(q_ct13_exp[pos_ct]),
                        "q_bh_statsmodels": float(q_ct13_sm[pos_ct]),
                        "abs_diff": float(abs(q_ct13_exp[pos_ct] - q_ct13_sm[pos_ct])),
                        "family_max_abs_diff": bh_diff_ct13,
                    }
                )
            pos_all = all_genes_pos[(grp, gene)]
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": "gao_genes_all14",
                    "unit": f"{grp}::{gene}",
                    "family_size": int(len(genes_combined)),
                    "raw_p": mat_p,
                    "q_bh_explicit": float(q_all14_exp[pos_all]),
                    "q_bh_statsmodels": float(q_all14_sm[pos_all]),
                    "abs_diff": float(abs(q_all14_exp[pos_all] - q_all14_sm[pos_all])),
                    "family_max_abs_diff": bh_diff_all14,
                }
            )

    # 3. Concordance row
    crow = conc_df.iloc[0]
    orig_c_p = float(orig_conc["perm_p_two_sided"])
    mat_c_p = float(crow["perm_p_two_sided"])
    verd_c, _ = verdict_standard(orig_c_p < ALPHA, mat_c_p < ALPHA, float(orig_conc["mean_pairwise_pearson_r"]), float(crow["mean_pairwise_pearson_r"]))
    out_rows.append(
        make_row(
            cohort_name,
            "All_13_CellTypes",
            "concordance",
            int(crow["n_donors_union"]),
            n_donors_orig - int(crow["n_donors_union"]),
            orig_c_p,
            orig_c_p,
            mat_c_p,
            mat_c_p,
            mat_c_p,
            mat_c_p,
            verd_c,
            f"n_pairs={int(crow['n_pairs'])};n_common_genes={int(crow['n_common_genes'])};median_r={float(crow['median_pairwise_pearson_r']):.4f};sign_agree={float(crow['mean_sign_agreement']):.4f};perms={int(crow['n_valid_perms'])};seed={ORIG_NULL_SEED}",
            arm="pooled_7study_primary",
            residual_df=np.nan,
            original_effect=float(orig_conc["mean_pairwise_pearson_r"]),
            mature_effect=float(crow["mean_pairwise_pearson_r"]),
            mature_se=np.nan,
            null_sd=float(crow["perm_null_sd"]),
            mde80=float(crow["mde80_null"]),
            mde80_method="2.4865_times_null_sd",
        )
    )

    # 4. Family C: transport_P1m and transport_P0m
    for tier, df_c, skips, orig_ref, fam_key in [
        ("P1m", fam_c_p1m, skip_p1m, orig_p1, "transport_P1m"),
        ("P0m", fam_c_p0m, skip_p0m, orig_p0, "transport_P0m"),
    ]:
        q_c_exp, q_c_sm, diff_c = bh_checked(df_c["composite_pvalue"])
        q_e_exp, q_e_sm, diff_e = bh_checked(df_c["empirical_null_pvalue"])
        df_c_idx = df_c.set_index("group")
        pos_map = {str(g): i for i, g in enumerate(df_c["group"])}

        for grp in sorted(orig_p1.index.tolist()):
            orig_n = int(orig_p1.loc[grp, "n_donors"])
            if tier == "P1m":
                op_c = float(orig_ref.loc[grp, "composite_pvalue"])
                oq_c = float(orig_ref.loc[grp, "q_bh_composite"])
                op_e = float(orig_ref.loc[grp, "empirical_null_pvalue"])
                oq_e = float(orig_ref.loc[grp, "q_bh_empirical"])
                oeff = float(orig_ref.loc[grp, "composite_effect"])
            else:
                op_c = float(orig_ref.loc[grp, "composite_pvalue"])
                oq_c = float(orig_ref.loc[grp, "composite_q"])
                op_e = float(orig_ref.loc[grp, "empirical_pvalue"])
                oq_e = float(orig_ref.loc[grp, "empirical_q"])
                oeff = float(orig_ref.loc[grp, "composite_slope"])
            orig_cls = transport_class(oq_c, oq_e, True)
            orig_both = orig_cls == "BOTH_PASS"

            if grp not in df_c_idx.index:
                reason = skips.get(grp, "ineligible_at_30plus")
                out_rows.append(
                    make_row(
                        cohort_name,
                        grp,
                        fam_key,
                        0,
                        orig_n,
                        op_c,
                        oq_c,
                        None,
                        None,
                        None,
                        None,
                        "UNTESTABLE",
                        f"reason={reason};orig_class={orig_cls}",
                        arm="pooled_7study_primary",
                        original_effect=oeff,
                        original_empirical_p=op_e,
                        original_empirical_q=oq_e,
                        orig_transport_class=orig_cls,
                        mature_transport_class="UNTESTABLE",
                    )
                )
                continue

            r = df_c_idx.loc[grp]
            pos = pos_map[grp]
            mat_cls = transport_class(float(r["q_bh_composite"]), float(r["q_bh_empirical"]), True)
            mat_both = mat_cls == "BOTH_PASS"
            verd, vnote = verdict_standard(orig_both, mat_both, oeff, float(r["composite_effect"]))
            mde_null = mde80_from_null_sd(float(r["null_sd_effect"]))
            mde_se = mde80_from_se(float(r["composite_se"]), int(r["residual_df"]))
            note_parts = [
                f"orig_class={orig_cls}",
                f"mature_class={mat_cls}",
                f"up_eff={float(r['up_effect']):+.4f}(p={float(r['up_pvalue']):.3g})",
                f"dn_eff={float(r['down_effect']):+.4f}(p={float(r['down_pvalue']):.3g})",
                f"mde80_se={mde_se:.4f}",
                f"draws={P1_ITERATIONS}",
                f"seed={ORIG_NULL_SEED}",
                f"min_decile_margin={int(r['min_decile_margin'])}",
            ]
            if vnote:
                note_parts.append(vnote)
            out_rows.append(
                make_row(
                    cohort_name,
                    grp,
                    fam_key,
                    int(r["n_donors"]),
                    orig_n - int(r["n_donors"]),
                    op_c,
                    oq_c,
                    float(r["composite_pvalue"]),
                    float(r["q_bh_composite"]),
                    float(r["empirical_null_pvalue"]),
                    float(r["q_bh_empirical"]),
                    verd,
                    ";".join(note_parts),
                    arm="pooled_7study_primary",
                    residual_df=int(r["residual_df"]),
                    original_effect=oeff,
                    mature_effect=float(r["composite_effect"]),
                    mature_se=float(r["composite_se"]),
                    null_sd=float(r["null_sd_effect"]),
                    mde80=mde_null,
                    mde80_method="2.4865_times_null_sd",
                    original_empirical_p=op_e,
                    original_empirical_q=oq_e,
                    orig_transport_class=orig_cls,
                    mature_transport_class=mat_cls,
                    n_up_detected=int(r["n_up_detected"]),
                    n_down_detected=int(r["n_down_detected"]),
                    covariates=str(r["covariates"]),
                    omitted_covariates=str(r["omitted_covariates"]),
                )
            )
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": f"gao_{fam_key}_composite_14",
                    "unit": grp,
                    "family_size": int(len(df_c)),
                    "raw_p": float(r["composite_pvalue"]),
                    "q_bh_explicit": float(q_c_exp[pos]),
                    "q_bh_statsmodels": float(q_c_sm[pos]),
                    "abs_diff": float(abs(q_c_exp[pos] - q_c_sm[pos])),
                    "family_max_abs_diff": diff_c,
                }
            )
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": f"gao_{fam_key}_empirical_14",
                    "unit": grp,
                    "family_size": int(len(df_c)),
                    "raw_p": float(r["empirical_null_pvalue"]),
                    "q_bh_explicit": float(q_e_exp[pos]),
                    "q_bh_statsmodels": float(q_e_sm[pos]),
                    "abs_diff": float(abs(q_e_exp[pos] - q_e_sm[pos])),
                    "family_max_abs_diff": diff_e,
                }
            )

    gao_df = order_cohort_df(out_rows)
    write_clean_csv(gao_df, RESULTS_DIR / f"{PREFIX}-gao.csv")

    summary_rows = [
        {
            "cohort": cohort_name,
            "endpoint_family": "fractions_clr",
            "arm": "pooled_7study_primary",
            "donors_original": n_donors_orig,
            "donors_mature": n_donors_30,
            "donors_excluded": n_excluded_cohort,
            "original_headline": "0/13 CLR hits (min q=0.3128, Endothelial p=0.0362)",
            "mature_headline": f"{int((fam_a['q_bh_family_a'] < ALPHA).sum())}/{len(fam_a)} CLR hits (min q={float(fam_a['q_bh_family_a'].min()):.4f}, top={fam_a.iloc[0]['cell_type']} p={float(fam_a.iloc[0]['pvalue']):.4g})",
            "mde80_summary": f"median={float(gao_df[gao_df['endpoint_family']=='fractions_clr']['mde80'].median()):.4f} (range {float(gao_df[gao_df['endpoint_family']=='fractions_clr']['mde80'].min()):.4f} to {float(gao_df[gao_df['endpoint_family']=='fractions_clr']['mde80'].max()):.4f}) CLR/dec (t_quantile_se)",
        },
        {
            "cohort": cohort_name,
            "endpoint_family": "genes",
            "arm": "pooled_7study_primary",
            "donors_original": n_donors_orig,
            "donors_mature": n_donors_30,
            "donors_excluded": n_excluded_cohort,
            "original_headline": "8/145873 ct13 hits; 12/159052 all14 hits",
            "mature_headline": f"{int((q_ct13_exp < ALPHA).sum())}/{int(ct_mask.sum())} ct13 hits; {int((q_all14_exp < ALPHA).sum())}/{len(genes_combined)} all14 hits",
            "mde80_summary": f"median cell-type median MDE80={float(sum_ct['median_mde80'].median()):.4f} (range {float(sum_ct['median_mde80'].min()):.4f} to {float(sum_ct['median_mde80'].max()):.4f}) log2CPM/dec (t_quantile_se)",
        },
        {
            "cohort": cohort_name,
            "endpoint_family": "concordance",
            "arm": "pooled_7study_primary",
            "donors_original": n_donors_orig,
            "donors_mature": int(crow["n_donors_union"]),
            "donors_excluded": n_donors_orig - int(crow["n_donors_union"]),
            "original_headline": f"null (mean r={float(orig_conc['mean_pairwise_pearson_r']):.4f}, perm_p={orig_c_p:.4f})",
            "mature_headline": f"{'sig' if mat_c_p < ALPHA else 'null'} (mean r={float(crow['mean_pairwise_pearson_r']):.4f}, perm_p={mat_c_p:.4f}, common_genes={int(crow['n_common_genes'])})",
            "mde80_summary": f"MDE80={float(crow['mde80_null']):.4f} (2.4865 * null_sd={float(crow['perm_null_sd']):.4f})",
        },
    ]
    for tier, df_c, fam_key in [("P1m", fam_c_p1m, "transport_P1m"), ("P0m", fam_c_p0m, "transport_P0m")]:
        bp = int(((df_c["q_bh_composite"] < ALPHA) & (df_c["q_bh_empirical"] < ALPHA)).sum())
        eo = int(((df_c["q_bh_composite"] >= ALPHA) & (df_c["q_bh_empirical"] < ALPHA)).sum())
        co = int(((df_c["q_bh_composite"] < ALPHA) & (df_c["q_bh_empirical"] >= ALPHA)).sum())
        endo = df_c[df_c["group"] == "Endocardial"].iloc[0]
        mde_s = gao_df[gao_df["endpoint_family"] == fam_key]["mde80"]
        summary_rows.append(
            {
                "cohort": cohort_name,
                "endpoint_family": fam_key,
                "arm": "pooled_7study_primary",
                "donors_original": n_donors_orig,
                "donors_mature": n_donors_30,
                "donors_excluded": n_excluded_cohort,
                "original_headline": (
                    "1/14 both-pass (Endocardial q_comp=0.0150, q_emp=0.0093); 6 empirical-only"
                    if tier == "P1m"
                    else f"{int(((orig_p0['composite_q'] < ALPHA) & (orig_p0['empirical_q'] < ALPHA)).sum())}/14 both-pass (Block 5 P0_Raw257)"
                ),
                "mature_headline": f"{bp}/{len(df_c)} both-pass, {eo} empirical-only, {co} composite-only; Endocardial eff={float(endo['composite_effect']):+.4f}, q_comp={float(endo['q_bh_composite']):.4f}, q_emp={float(endo['q_bh_empirical']):.4f}",
                "mde80_summary": f"median={float(mde_s.median()):.4f} (range {float(mde_s.min()):.4f} to {float(mde_s.max()):.4f}) per dec (2.4865 * null_sd)",
            }
        )

    save_cohort_payload(
        "gao",
        {
            "summary_rows": summary_rows,
            "bh_rows": bh_rows,
            "diagnostics": {
                "n_donors_orig": n_donors_orig,
                "n_donors_30": n_donors_30,
                "n_excluded": n_excluded_cohort,
                "age_min": float(donor_meta_30["age_years"].min()),
                "age_max": float(donor_meta_30["age_years"].max()),
                "n_male": int(donor_meta_30["sex_male"].sum()),
                "n_female": int((1 - donor_meta_30["sex_male"]).sum()),
                "study_counts": donor_meta_30["study"].value_counts().to_dict(),
                "study_anova": study_anova,
                "gate_totals": gate_totals,
                "bh_diff_a": bh_diff_a,
                "bh_diff_ct13": bh_diff_ct13,
                "bh_diff_all14": bh_diff_all14,
            },
        },
    )
    append_runtime(cohort_name, "cohort_total", s0_all, t0_all, "ok", f"rows={len(gao_df)}")
    print(f"[Gao DONE] rows={len(gao_df)}, elapsed={time.monotonic() - t0_all:.1f}s", flush=True)


# =============================================================================
# COHORT 2: YAZAR BLOOD (HCATA OneK1K PBMC, GSE196830)
# =============================================================================

def add_yazar_pool_dummies_mature(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    out = df.copy()
    vc = out["pool_number"].astype(int).value_counts()
    pools_sorted = [int(p) for p, _ in sorted(vc.items(), key=lambda kv: (-int(kv[1]), int(kv[0])))]
    dummy_cols: list[str] = []
    for p in pools_sorted[1:]:
        col = f"pool_{p}"
        out[col] = (out["pool_number"].astype(int) == p).astype(float)
        dummy_cols.append(col)
    return out, dummy_cols


def run_cohort_yazar(programs: dict[str, dict[str, list[str]]]) -> None:
    cohort_name = "Yazar_Blood"
    t0_all = time.monotonic()
    s0_all = utc_now()
    wait_for_memory("Yazar_Blood stream")

    print("\n================ COHORT 2: YAZAR BLOOD (HCATA OneK1K) ================", flush=True)
    yazar_cache_pkl = TMP_WORK_DIR / "yazar_30plus_raw_pb.pkl"

    if yazar_cache_pkl.exists():
        print(f"[Yazar] Loading verified 30-plus raw pseudobulk cache {yazar_cache_pkl}", flush=True)
        with yazar_cache_pkl.open("rb") as fh:
            ydata = pickle.load(fh)
        assert float(ydata["donor_meta_30"]["age_years"].min()) >= MIN_AGE
    else:
        t0_stream = time.monotonic()
        s0_stream = utc_now()
        ydata = stream_hcata_mature_from_h5ad(
            study="yazar",
            h5_path=YAZAR_H5AD,
            expected_n_obs=1248980,
            expected_total_nnz=1183532685,
            expected_raw_source="X",
            expected_orig_donors_csv=RESULTS_DIR / "followup-hcata-yazar-donors.csv",
            min_donors_extract=6,
        )
        # Cross-verify bit-for-bit parity against existing raw yazar_cache for 30-plus donors
        manifest = json.loads((YAZAR_CACHE / "pb_manifest.json").read_text(encoding="utf-8"))
        for grp in ("All_Cells", "natural killer cell", "mucosal-associated invariant T cell"):
            with open(YAZAR_CACHE / manifest[grp], "rb") as pf:
                old_raw = pickle.load(pf)
            old_raw.columns = old_raw.columns.astype(str)
            new_raw = ydata["raw_pb_by_group"][grp]
            cols = list(new_raw.columns)
            np.testing.assert_array_equal(new_raw.to_numpy(), old_raw[cols].to_numpy())
            del old_raw
        print("[Yazar] Verified bit-for-bit raw count parity with yazar_cache on 30-plus donors.", flush=True)
        gc.collect()
        TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
        tmp_p = yazar_cache_pkl.with_suffix(".tmp")
        with tmp_p.open("wb") as fh:
            pickle.dump(ydata, fh, protocol=pickle.HIGHEST_PROTOCOL)
        tmp_p.replace(yazar_cache_pkl)
        append_runtime(
            cohort_name,
            "stream_h5ad_30plus",
            s0_stream,
            t0_stream,
            "ok",
            f"932 mature donors, cells_30={ydata['gate_30']['cells_30']}, bit_parity_yazar_cache=PASS",
        )

    donor_meta_30: pd.DataFrame = ydata["donor_meta_30"].copy()
    donor_meta_30["pool_number"] = donor_meta_30["pool_number"].astype(int)
    cell_counts_30: pd.DataFrame = ydata["cell_counts_30"].copy()
    raw_pb_by_group: dict[str, pd.DataFrame] = ydata["raw_pb_by_group"]
    n_donors_orig = int(ydata["n_donors_orig"])
    n_donors_30 = int(ydata["n_donors_30"])
    n_excluded_cohort = int(ydata["n_donors_excluded"])
    assert n_donors_orig == 981 and n_donors_30 == 932

    # Pool-age confounding stats at 30-plus
    pool_anova = anova_batch_age(donor_meta_30, "pool_number", single_only=False)
    pool_vc = donor_meta_30["pool_number"].value_counts()
    pool_stats = {
        "n_pools": int(len(pool_vc)),
        "min_donors_per_pool": int(pool_vc.min()),
        "median_donors_per_pool": float(pool_vc.median()),
        "max_donors_per_pool": int(pool_vc.max()),
        "n_pools_lt2": int((pool_vc < 2).sum()),
        "anova": pool_anova,
    }

    # Drop Azimuth doublets for Family A CLR (exact 28 cell types)
    fam_a_counts = cell_counts_30.drop(columns=[YAZAR_DOUBLETS_LABEL], errors="ignore")
    assert fam_a_counts.shape[1] == 28, f"Expected 28 non-doublet cell types, got {fam_a_counts.shape[1]}"

    # Load original pool-FE & unadjusted references
    orig_pool_frac = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-fractions.csv")
    orig_pool_frac = orig_pool_frac[orig_pool_frac["study"] == "yazar"].set_index("cell_type")
    orig_pool_gsum = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-genes-summary.csv")
    orig_pool_gsum = orig_pool_gsum[orig_pool_gsum["study"] == "yazar"].set_index("group")
    orig_pool_conc = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-concordance.csv")
    orig_pool_conc = orig_pool_conc[orig_pool_conc["study"] == "yazar"].iloc[0]
    orig_pool_p1 = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-p1-signature.csv")
    orig_pool_p1 = orig_pool_p1[orig_pool_p1["study"] == "yazar"].set_index("group")
    orig_pool_elig = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-pool-fe-eligibility.csv").set_index("group")

    orig_unadj_frac = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-fractions.csv")
    orig_unadj_frac = orig_unadj_frac[orig_unadj_frac["arm"] == "primary"].set_index("cell_type")
    orig_unadj_gsum = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-genes-summary.csv").set_index("group")
    orig_unadj_conc = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-concordance.csv")
    orig_unadj_conc = orig_unadj_conc[orig_unadj_conc["arm"] == "primary"].iloc[0]
    orig_unadj_p1 = pd.read_csv(RESULTS_DIR / "followup-hcata-yazar-p1-signature.csv")
    orig_unadj_p1 = orig_unadj_p1[orig_unadj_p1["arm"] == "primary"].set_index("group")

    orig_b5 = pd.read_csv(RESULTS_DIR / "followup-11pt-block5-p0-transport.csv")
    orig_p0_pool = orig_b5[
        (orig_b5["cohort"] == "Yazar_Blood") & (orig_b5["arm"] == "primary") & (orig_b5["gene_set"] == "P0_Raw257")
    ].set_index("unit")

    # Family A: pool_fe (primary) and unadjusted (sensitivity)
    t0_a = time.monotonic()
    s0_a = utc_now()
    donor_clr, clr_pool_cols = add_yazar_pool_dummies_mature(donor_meta_30)
    fam_a_fe = run_family_a_clr(
        fam_a_counts,
        donor_clr,
        age_col="age_per_decade",
        candidate_covariates=("sex_male", *clr_pool_cols),
        min_donors_present=3,
        min_total_cells=MIN_CELLS_PER_DONOR,
        min_residual_df=10,
    )
    q_afe_exp, q_afe_sm, bh_diff_afe = bh_checked(fam_a_fe["pvalue"])
    fam_a_fe["q_bh_family_a"] = q_afe_exp

    fam_a_unadj = run_family_a_clr(
        fam_a_counts,
        donor_meta_30,
        age_col="age_per_decade",
        candidate_covariates=("sex_male",),
        min_donors_present=6,
        min_total_cells=100,
        min_residual_df=2,
    )
    q_aun_exp, q_aun_sm, bh_diff_aun = bh_checked(fam_a_unadj["pvalue"])
    fam_a_unadj["q_bh_family_a"] = q_aun_exp
    append_runtime(cohort_name, "family_a_fractions", s0_a, t0_a, "ok", f"pool_fe={len(fam_a_fe)}, unadj={len(fam_a_unadj)}")

    # Finalize 30-plus pseudobulks (trim at >=35% for large N>=50 just like run_followup_hcata_studies to bound RAM,
    # strictly looser than >=50% detection filter inside Family B and C)
    wait_for_memory("Yazar_Blood finalize pseudobulks")
    expr_30: dict[str, pd.DataFrame] = {}
    meta_30: dict[str, pd.DataFrame] = {}
    thresh = float(np.log2(2.0))
    for grp, raw_df in raw_pb_by_group.items():
        cols = list(raw_df.columns)
        if grp == "All_Cells":
            g_counts = fam_a_counts.loc[cols].sum(axis=1)
        else:
            g_counts = cell_counts_30.loc[cols, grp]
        l2, gm, _ = finalize_pseudobulk_group(raw_df, donor_meta_30.loc[cols], g_counts)
        if len(cols) >= 50:
            min_keep = max(3, int(np.floor(0.35 * len(cols))))
            keep_m = (l2.to_numpy(dtype=np.float64) >= thresh).sum(axis=1) >= min_keep
            l2 = l2.loc[keep_m].copy()
        expr_30[grp] = l2
        meta_30[grp] = gm
    del raw_pb_by_group, ydata
    gc.collect()

    # Build pool-FE subsets (pools with >=2 mature donors in group, est_df >= 10)
    all_pools = sorted(donor_meta_30["pool_number"].unique().tolist())
    all_pool_cols = [f"pool_{p}" for p in all_pools]
    fe_expr: dict[str, pd.DataFrame] = {}
    fe_meta: dict[str, pd.DataFrame] = {}
    pool_elig_rows: list[dict[str, Any]] = []

    for grp in sorted(meta_30.keys()):
        m_df = meta_30[grp]
        vc = m_df["pool_number"].value_counts()
        ge2 = vc[vc >= 2]
        n_lt2 = int((vc < 2).sum())
        keep_donors = m_df.index[m_df["pool_number"].isin(ge2.index)].tolist()
        n_ge2 = len(keep_donors)
        n_pools_ge2 = int(len(ge2))
        est_df = n_ge2 - (4 + max(0, n_pools_ge2 - 1))
        if n_pools_ge2 < 2 or est_df < 10:
            pool_elig_rows.append(
                {
                    "group": grp,
                    "n_donors_30plus": int(len(m_df)),
                    "n_donors_ge2_pools": n_ge2,
                    "n_pools_total_present": int(len(vc)),
                    "n_pools_ge2": n_pools_ge2,
                    "n_pools_singleton_dropped": n_lt2,
                    "estimated_residual_df": int(est_df),
                    "status": "UNTESTABLE",
                }
            )
            continue
        sub_m = m_df.loc[keep_donors].copy()
        ref_p = int(ge2.index[0])
        for p in all_pools:
            sub_m[f"pool_{p}"] = (sub_m["pool_number"] == p).astype(float) if p != ref_p else 0.0
        fe_meta[grp] = sub_m
        fe_expr[grp] = expr_30[grp][keep_donors]
        pool_elig_rows.append(
            {
                "group": grp,
                "n_donors_30plus": int(len(m_df)),
                "n_donors_ge2_pools": n_ge2,
                "n_pools_total_present": int(len(vc)),
                "n_pools_ge2": n_pools_ge2,
                "n_pools_singleton_dropped": n_lt2,
                "estimated_residual_df": int(est_df),
                "status": "TESTED",
            }
        )
    pool_elig_df = pd.DataFrame(pool_elig_rows).set_index("group")

    # Family B: pool_fe (primary) and unadjusted (sensitivity summary)
    wait_for_memory("Yazar_Blood Family B pool_fe")
    t0_b = time.monotonic()
    s0_b = utc_now()
    cov_fe = ("sex_male", "log_lib", *all_pool_cols)
    fe_ct_expr = {g: fe_expr[g] for g in fe_expr if g != "All_Cells"}
    fe_ct_meta = {g: fe_meta[g] for g in fe_meta if g != "All_Cells"}

    genes_fe_ct, sum_fe_ct, conc_fe = run_family_b_mature(
        fe_ct_expr,
        fe_ct_meta,
        age_col="age_per_decade",
        candidate_covariates=cov_fe,
        min_donors=4,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=10,
        run_concordance=True,
        concordance_perms=CONCORDANCE_PERMS,
        seed=ORIG_NULL_SEED,
    )
    genes_fe_pool, sum_fe_pool, _ = run_family_b_mature(
        {"All_Cells": fe_expr["All_Cells"]},
        {"All_Cells": fe_meta["All_Cells"]},
        age_col="age_per_decade",
        candidate_covariates=cov_fe,
        min_donors=4,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=10,
        run_concordance=False,
        seed=ORIG_NULL_SEED,
    )
    q_b_ct_exp, q_b_ct_sm, bh_diff_b_ct = bh_checked(genes_fe_ct["pvalue"])
    q_b_pl_exp, q_b_pl_sm, bh_diff_b_pl = bh_checked(genes_fe_pool["pvalue"])
    n_fe_ct_tests = int(len(genes_fe_ct))
    n_fe_ct_sig = int((q_b_ct_exp < ALPHA).sum())
    n_fe_pool_tests = int(len(genes_fe_pool))
    n_fe_pool_sig = int((q_b_pl_exp < ALPHA).sum())

    # Capture top min-p gene per group for BH check BEFORE freeing giant genes_fe_ct DataFrame
    yazar_gene_bh_entries: list[dict[str, Any]] = []
    for grp in sorted(fe_ct_expr.keys()):
        g_rows = genes_fe_ct[genes_fe_ct["group"] == grp]
        idx_min = g_rows["pvalue"].idxmin()
        row_min = genes_fe_ct.loc[idx_min]
        pos_min = int(genes_fe_ct.index.get_loc(idx_min))
        yazar_gene_bh_entries.append(
            {
                "cohort": cohort_name,
                "family": "yazar_pool_fe_genes_ct",
                "unit": f"{grp}::{row_min['gene']}",
                "family_size": n_fe_ct_tests,
                "raw_p": float(row_min["pvalue"]),
                "q_bh_explicit": float(q_b_ct_exp[pos_min]),
                "q_bh_statsmodels": float(q_b_ct_sm[pos_min]),
                "abs_diff": float(abs(q_b_ct_exp[pos_min] - q_b_ct_sm[pos_min])),
                "family_max_abs_diff": bh_diff_b_ct,
            }
        )
    idx_min_p = genes_fe_pool["pvalue"].idxmin()
    row_min_p = genes_fe_pool.loc[idx_min_p]
    pos_min_p = int(genes_fe_pool.index.get_loc(idx_min_p))
    yazar_gene_bh_entries.append(
        {
            "cohort": cohort_name,
            "family": "yazar_pool_fe_genes_all_cells",
            "unit": f"All_Cells::{row_min_p['gene']}",
            "family_size": n_fe_pool_tests,
            "raw_p": float(row_min_p["pvalue"]),
            "q_bh_explicit": float(q_b_pl_exp[pos_min_p]),
            "q_bh_statsmodels": float(q_b_pl_sm[pos_min_p]),
            "abs_diff": float(abs(q_b_pl_exp[pos_min_p] - q_b_pl_sm[pos_min_p])),
            "family_max_abs_diff": bh_diff_b_pl,
        }
    )
    del genes_fe_ct, genes_fe_pool
    gc.collect()
    append_runtime(
        cohort_name,
        "family_b_pool_fe",
        s0_b,
        t0_b,
        "ok",
        f"ct_tests={n_fe_ct_tests}, ct_sig={n_fe_ct_sig}, pool_sig={n_fe_pool_sig}, conc_r={float(conc_fe.iloc[0]['mean_pairwise_pearson_r']):.4f}",
    )

    # Family B unadjusted diagnostic (for pool-confounding attenuation reporting)
    wait_for_memory("Yazar_Blood Family B unadjusted")
    t0_bu = time.monotonic()
    s0_bu = utc_now()
    unadj_ct_expr = {g: expr_30[g] for g in expr_30 if g != "All_Cells"}
    unadj_ct_meta = {g: meta_30[g] for g in meta_30 if g != "All_Cells"}
    genes_un_ct, sum_un_ct, conc_un = run_family_b_mature(
        unadj_ct_expr,
        unadj_ct_meta,
        age_col="age_per_decade",
        candidate_covariates=("sex_male", "log_lib"),
        min_donors=6,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        run_concordance=True,
        concordance_perms=CONCORDANCE_PERMS,
        seed=ORIG_NULL_SEED,
    )
    n_un_ct_tests = int(len(genes_un_ct))
    n_un_ct_sig = int((genes_un_ct["q_bh_global"] < ALPHA).sum())
    del genes_un_ct
    gc.collect()
    append_runtime(
        cohort_name,
        "family_b_unadjusted",
        s0_bu,
        t0_bu,
        "ok",
        f"unadj_ct_tests={n_un_ct_tests}, unadj_ct_sig={n_un_ct_sig}, conc_r={float(conc_un.iloc[0]['mean_pairwise_pearson_r']):.4f}",
    )

    # Family C: pool_fe P1m and P0m (cell types separate from All_Cells, matching original + fixed all-groups BH)
    wait_for_memory("Yazar_Blood Family C pool_fe")
    t0_c = time.monotonic()
    s0_c = utc_now()
    fam_c_results: dict[str, tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]] = {}
    for tier in ("P1m", "P0m"):
        df_ct, sk_ct = run_family_c_mature(
            fe_ct_expr,
            fe_ct_meta,
            program=programs[tier],
            candidate_covariates=cov_fe,
            min_donors=4,
            min_residual_df=10,
            iterations=P1_ITERATIONS,
            seed=ORIG_NULL_SEED,
        )
        df_pl, sk_pl = run_family_c_mature(
            {"All_Cells": fe_expr["All_Cells"]},
            {"All_Cells": fe_meta["All_Cells"]},
            program=programs[tier],
            candidate_covariates=cov_fe,
            min_donors=4,
            min_residual_df=10,
            iterations=P1_ITERATIONS,
            seed=ORIG_NULL_SEED,
        )
        fam_c_results[tier] = (df_ct, df_pl, {**sk_ct, **sk_pl})
    append_runtime(
        cohort_name,
        "family_c_pool_fe",
        s0_c,
        t0_c,
        "ok",
        f"P1m_ct={len(fam_c_results['P1m'][0])}, P0m_ct={len(fam_c_results['P0m'][0])}, seed={ORIG_NULL_SEED}",
    )

    # Assemble Yazar cohort CSV rows and BH check rows
    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = list(yazar_gene_bh_entries)

    # 1. Fractions (pool_fe primary)
    for idx_r, r in fam_a_fe.iterrows():
        ct = str(r["cell_type"])
        orow = orig_pool_frac.loc[ct]
        orig_sig = bool(float(orow["q_bh_family_a"]) < ALPHA)
        mat_sig = bool(float(r["q_bh_family_a"]) < ALPHA)
        verd, vnote = verdict_standard(orig_sig, mat_sig, float(orow["effect_clr"]), float(r["effect_clr"]))
        mde = mde80_from_se(float(r["se"]), int(r["residual_df"]))
        n_pool_act = str(r["covariates"]).count("pool_")
        out_rows.append(
            make_row(
                cohort_name,
                ct,
                "fractions_clr",
                int(r["n_donors"]),
                int(orow["n_donors"]) - int(r["n_donors"]),
                float(orow["pvalue"]),
                float(orow["q_bh_family_a"]),
                float(r["pvalue"]),
                float(r["q_bh_family_a"]),
                None,
                None,
                verd,
                f"pool_fe({n_pool_act}_dummies);n_present={int(r['n_donors_present'])};{vnote}".strip(";"),
                arm="pool_fe",
                residual_df=int(r["residual_df"]),
                original_effect=float(orow["effect_clr"]),
                mature_effect=float(r["effect_clr"]),
                mature_se=float(r["se"]),
                null_sd=np.nan,
                mde80=mde,
                mde80_method="t_quantile_se",
                covariates=f"sex_male;pool_fe({n_pool_act}_dummies)",
                omitted_covariates="ref_or_unused_pools_only",
            )
        )
        bh_rows.append(
            {
                "cohort": cohort_name,
                "family": "yazar_pool_fe_fractions_clr_28",
                "unit": ct,
                "family_size": int(len(fam_a_fe)),
                "raw_p": float(r["pvalue"]),
                "q_bh_explicit": float(q_afe_exp[idx_r]),
                "q_bh_statsmodels": float(q_afe_sm[idx_r]),
                "abs_diff": float(abs(q_afe_exp[idx_r] - q_afe_sm[idx_r])),
                "family_max_abs_diff": bh_diff_afe,
            }
        )

    # 2. Genes per cell type (pool_fe: 16 tested + All_Cells + 3 untestable)
    sum_fe_all = pd.concat([sum_fe_ct, sum_fe_pool], ignore_index=True).set_index("group")
    for grp in sorted(orig_pool_elig.index.tolist()):
        erow = pool_elig_df.loc[grp]
        orig_n_ge2 = int(orig_pool_elig.loc[grp, "n_donors_ge2_pools"])
        if erow["status"] != "TESTED" or grp not in sum_fe_all.index:
            out_rows.append(
                make_row(
                    cohort_name,
                    grp,
                    "genes_celltype",
                    int(erow["n_donors_ge2_pools"]),
                    int(orig_pool_elig.loc[grp, "n_donors_original"]) - int(erow["n_donors_30plus"]),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    "UNTESTABLE",
                    f"not_evaluable_within_pool_df;n30={int(erow['n_donors_30plus'])};n_ge2={int(erow['n_donors_ge2_pools'])};pools_ge2={int(erow['n_pools_ge2'])};est_df={int(erow['estimated_residual_df'])}<10",
                    arm="pool_fe",
                    residual_df=int(erow["estimated_residual_df"]),
                )
            )
            continue
        r = sum_fe_all.loc[grp]
        orow = orig_pool_gsum.loc[grp]
        orig_nsig = int(orow["n_sig_global_q05"])
        mat_nsig = int(r["n_sig_global_q05"])
        verd, _ = verdict_standard(orig_nsig > 0, mat_nsig > 0)
        n_pool_act = str(r["covariates"]).count("pool_")
        out_rows.append(
            make_row(
                cohort_name,
                grp,
                "genes_celltype",
                int(r["n_donors"]),
                orig_n_ge2 - int(r["n_donors"]),
                float(orow["min_pvalue"]),
                float(orow["min_q_bh_global"]),
                float(r["min_pvalue"]),
                float(r["min_q_bh_global"]),
                None,
                None,
                verd,
                f"n_genes={int(r['n_genes_tested'])};orig_sig={orig_nsig};mature_sig={mat_nsig};top={r['top_gene']};pools_ge2={int(erow['n_pools_ge2'])};singletons_dropped={int(erow['n_pools_singleton_dropped'])}",
                arm="pool_fe",
                residual_df=int(r["residual_df"]),
                original_effect=np.nan,
                mature_effect=float(r["top_gene_effect"]),
                mature_se=float(r["median_se"]),
                null_sd=np.nan,
                mde80=float(r["median_mde80"]),
                mde80_method="median_t_quantile_se",
                n_genes_tested=int(r["n_genes_tested"]),
                n_sig_q05=mat_nsig,
                covariates=f"sex_male;log_lib;pool_fe({n_pool_act}_dummies)",
                omitted_covariates="ref_or_unused_pools_only",
            )
        )

    # 3. Concordance (pool_fe primary + unadjusted sensitivity)
    for label, cdf, o_cdf, fam_k, arm_k in [
        ("All_16_PoolFE_CellTypes", conc_fe.iloc[0], orig_pool_conc, "concordance", "pool_fe"),
        ("All_19_Unadjusted_CellTypes", conc_un.iloc[0], orig_unadj_conc, "concordance_unadjusted", "unadjusted"),
    ]:
        op = float(o_cdf["perm_p_two_sided"])
        mp = float(cdf["perm_p_two_sided"])
        verd, _ = verdict_standard(op < ALPHA, mp < ALPHA, float(o_cdf["mean_pairwise_pearson_r"]), float(cdf["mean_pairwise_pearson_r"]))
        out_rows.append(
            make_row(
                cohort_name,
                label,
                fam_k,
                int(cdf["n_donors_union"]),
                n_donors_orig - int(cdf["n_donors_union"]),
                op,
                op,
                mp,
                mp,
                mp,
                mp,
                verd,
                f"n_ct={int(cdf['n_cell_types'])};n_pairs={int(cdf['n_pairs'])};n_common_genes={int(cdf['n_common_genes'])};median_r={float(cdf['median_pairwise_pearson_r']):.4f};sign_agree={float(cdf['mean_sign_agreement']):.4f};perms={int(cdf['n_valid_perms'])};seed={ORIG_NULL_SEED}",
                arm=arm_k,
                original_effect=float(o_cdf["mean_pairwise_pearson_r"]),
                mature_effect=float(cdf["mean_pairwise_pearson_r"]),
                null_sd=float(cdf["perm_null_sd"]),
                mde80=float(cdf["mde80_null"]),
                mde80_method="2.4865_times_null_sd",
            )
        )

    # 4. Family C: transport_P1m and transport_P0m (pool_fe)
    for tier, fam_key in [("P1m", "transport_P1m"), ("P0m", "transport_P0m")]:
        df_ct, df_pl, skips = fam_c_results[tier]
        q_c_exp, q_c_sm, diff_c = bh_checked(df_ct["composite_pvalue"])
        q_e_exp, q_e_sm, diff_e = bh_checked(df_ct["empirical_null_pvalue"])
        pos_ct = {str(g): i for i, g in enumerate(df_ct["group"])}

        # Combined 17-group BH (fixing original add_bh bug!)
        df_both = pd.concat([df_ct.assign(is_pooled=False), df_pl.assign(is_pooled=True)], ignore_index=True)
        q_all_c_exp, q_all_c_sm, diff_all_c = bh_checked(df_both["composite_pvalue"])
        q_all_e_exp, q_all_e_sm, diff_all_e = bh_checked(df_both["empirical_null_pvalue"])
        df_both["q_bh_all_composite"] = q_all_c_exp
        df_both["q_bh_all_empirical"] = q_all_e_exp
        df_both_idx = df_both.set_index("group")
        pos_all = {str(g): i for i, g in enumerate(df_both["group"])}

        for grp in sorted(orig_pool_elig.index.tolist()):
            erow = pool_elig_df.loc[grp]
            if erow["status"] != "TESTED" or grp not in df_both_idx.index:
                out_rows.append(
                    make_row(
                        cohort_name,
                        grp,
                        fam_key,
                        int(erow["n_donors_ge2_pools"]),
                        int(orig_pool_elig.loc[grp, "n_donors_original"]) - int(erow["n_donors_30plus"]),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        "UNTESTABLE",
                        f"not_evaluable_within_pool_df;n30={int(erow['n_donors_30plus'])};n_ge2={int(erow['n_donors_ge2_pools'])};pools_ge2={int(erow['n_pools_ge2'])};est_df={int(erow['estimated_residual_df'])}<10",
                        arm="pool_fe",
                        orig_transport_class="UNTESTABLE",
                        mature_transport_class="UNTESTABLE",
                    )
                )
                continue

            r = df_both_idx.loc[grp]
            is_pool = bool(r["is_pooled"])
            orig_n_ge2 = int(orig_pool_elig.loc[grp, "n_donors_ge2_pools"])
            if tier == "P1m":
                oref = orig_pool_p1.loc[grp]
                op_c = float(oref["composite_pvalue"])
                oq_c = float(oref["q_bh_composite"]) if not is_pool else op_c
                op_e = float(oref["empirical_null_pvalue"])
                oq_e = float(oref["q_bh_empirical"]) if not is_pool else op_e
                oeff = float(oref["composite_effect"])
            else:
                oref = orig_p0_pool.loc[grp]
                op_c = float(oref["composite_pvalue"])
                oq_c = float(oref["composite_q"])
                op_e = float(oref["empirical_pvalue"])
                oq_e = float(oref["empirical_q"])
                oeff = float(oref["composite_slope"])

            mat_q_c = float(r["q_bh_all_composite"])
            mat_q_e = float(r["q_bh_all_empirical"])
            orig_cls = transport_class(oq_c, oq_e, True)
            mat_cls = transport_class(mat_q_c, mat_q_e, True)
            verd, vnote = verdict_standard(orig_cls == "BOTH_PASS", mat_cls == "BOTH_PASS", oeff, float(r["composite_effect"]))
            mde_null = mde80_from_null_sd(float(r["null_sd_effect"]))
            mde_se = mde80_from_se(float(r["composite_se"]), int(r["residual_df"]))
            n_pool_act = str(r["covariates"]).count("pool_")
            notes = (
                f"orig_class={orig_cls};mature_class={mat_cls};"
                f"q_all17_comp={float(r['q_bh_all_composite']):.4g};q_all17_emp={float(r['q_bh_all_empirical']):.4g};"
                f"up_eff={float(r['up_effect']):+.4f}(p={float(r['up_pvalue']):.3g});"
                f"dn_eff={float(r['down_effect']):+.4f}(p={float(r['down_pvalue']):.3g});"
                f"mde80_se={mde_se:.4f};draws={P1_ITERATIONS};seed={ORIG_NULL_SEED};"
                f"pools_ge2={int(erow['n_pools_ge2'])};min_decile_margin={int(r['min_decile_margin'])}"
            )
            if vnote:
                notes += f";{vnote}"
            out_rows.append(
                make_row(
                    cohort_name,
                    grp,
                    fam_key,
                    int(r["n_donors"]),
                    orig_n_ge2 - int(r["n_donors"]),
                    op_c,
                    oq_c,
                    float(r["composite_pvalue"]),
                    mat_q_c,
                    float(r["empirical_null_pvalue"]),
                    mat_q_e,
                    verd,
                    notes,
                    arm="pool_fe",
                    residual_df=int(r["residual_df"]),
                    original_effect=oeff,
                    mature_effect=float(r["composite_effect"]),
                    mature_se=float(r["composite_se"]),
                    null_sd=float(r["null_sd_effect"]),
                    mde80=mde_null,
                    mde80_method="2.4865_times_null_sd",
                    original_empirical_p=op_e,
                    original_empirical_q=oq_e,
                    orig_transport_class=orig_cls,
                    mature_transport_class=mat_cls,
                    n_up_detected=int(r["n_up_detected"]),
                    n_down_detected=int(r["n_down_detected"]),
                    covariates=f"sex_male;log_lib;pool_fe({n_pool_act}_dummies)",
                    omitted_covariates="ref_or_unused_pools_only",
                )
            )
            if not is_pool:
                p_i = pos_ct[grp]
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"yazar_pool_fe_{fam_key}_ct16_composite",
                        "unit": grp,
                        "family_size": int(len(df_ct)),
                        "raw_p": float(r["composite_pvalue"]),
                        "q_bh_explicit": float(q_c_exp[p_i]),
                        "q_bh_statsmodels": float(q_c_sm[p_i]),
                        "abs_diff": float(abs(q_c_exp[p_i] - q_c_sm[p_i])),
                        "family_max_abs_diff": diff_c,
                    }
                )
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"yazar_pool_fe_{fam_key}_ct16_empirical",
                        "unit": grp,
                        "family_size": int(len(df_ct)),
                        "raw_p": float(r["empirical_null_pvalue"]),
                        "q_bh_explicit": float(q_e_exp[p_i]),
                        "q_bh_statsmodels": float(q_e_sm[p_i]),
                        "abs_diff": float(abs(q_e_exp[p_i] - q_e_sm[p_i])),
                        "family_max_abs_diff": diff_e,
                    }
                )
            pa = pos_all[grp]
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": f"yazar_pool_fe_{fam_key}_all17_composite",
                    "unit": grp,
                    "family_size": int(len(df_both)),
                    "raw_p": float(r["composite_pvalue"]),
                    "q_bh_explicit": float(q_all_c_exp[pa]),
                    "q_bh_statsmodels": float(q_all_c_sm[pa]),
                    "abs_diff": float(abs(q_all_c_exp[pa] - q_all_c_sm[pa])),
                    "family_max_abs_diff": diff_all_c,
                }
            )
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": f"yazar_pool_fe_{fam_key}_all17_empirical",
                    "unit": grp,
                    "family_size": int(len(df_both)),
                    "raw_p": float(r["empirical_null_pvalue"]),
                    "q_bh_explicit": float(q_all_e_exp[pa]),
                    "q_bh_statsmodels": float(q_all_e_sm[pa]),
                    "abs_diff": float(abs(q_all_e_exp[pa] - q_all_e_sm[pa])),
                    "family_max_abs_diff": diff_all_e,
                }
            )

    yazar_df = order_cohort_df(out_rows)
    write_clean_csv(yazar_df, RESULTS_DIR / f"{PREFIX}-yazar.csv")

    summary_rows = [
        {
            "cohort": cohort_name,
            "endpoint_family": "fractions_clr",
            "arm": "pool_fe",
            "donors_original": n_donors_orig,
            "donors_mature": n_donors_30,
            "donors_excluded": n_excluded_cohort,
            "original_headline": "21/28 pool-FE CLR hits (naive CD8 q=3.11e-92)",
            "mature_headline": f"{int((fam_a_fe['q_bh_family_a'] < ALPHA).sum())}/{len(fam_a_fe)} pool-FE CLR hits (unadj: {int((fam_a_unadj['q_bh_family_a'] < ALPHA).sum())}/{len(fam_a_unadj)}); top={fam_a_fe.iloc[0]['cell_type']} q={float(fam_a_fe.iloc[0]['q_bh_family_a']):.3e}",
            "mde80_summary": f"median={float(yazar_df[yazar_df['endpoint_family']=='fractions_clr']['mde80'].median()):.4f} (range {float(yazar_df[yazar_df['endpoint_family']=='fractions_clr']['mde80'].min()):.4f} to {float(yazar_df[yazar_df['endpoint_family']=='fractions_clr']['mde80'].max()):.4f}) CLR/dec (t_quantile_se)",
        },
        {
            "cohort": cohort_name,
            "endpoint_family": "genes",
            "arm": "pool_fe",
            "donors_original": n_donors_orig,
            "donors_mature": n_donors_30,
            "donors_excluded": n_excluded_cohort,
            "original_headline": "1741/134036 pool-FE ct hits (16 types); 6241/13103 All_Cells",
            "mature_headline": f"{n_fe_ct_sig}/{n_fe_ct_tests} pool-FE ct hits (16 types, 3 untestable; unadj={n_un_ct_sig}/{n_un_ct_tests}); {n_fe_pool_sig}/{n_fe_pool_tests} All_Cells",
            "mde80_summary": f"median cell-type median MDE80={float(sum_fe_ct['median_mde80'].median()):.4f} (range {float(sum_fe_ct['median_mde80'].min()):.4f} to {float(sum_fe_ct['median_mde80'].max()):.4f}) log2CPM/dec (t_quantile_se)",
        },
        {
            "cohort": cohort_name,
            "endpoint_family": "concordance",
            "arm": "pool_fe",
            "donors_original": n_donors_orig,
            "donors_mature": int(conc_fe.iloc[0]["n_donors_union"]),
            "donors_excluded": n_donors_orig - int(conc_fe.iloc[0]["n_donors_union"]),
            "original_headline": f"pool-FE mean r={float(orig_pool_conc['mean_pairwise_pearson_r']):.4f}, perm_p={float(orig_pool_conc['perm_p_two_sided']):.4f} (16 types)",
            "mature_headline": f"pool-FE mean r={float(conc_fe.iloc[0]['mean_pairwise_pearson_r']):.4f}, perm_p={float(conc_fe.iloc[0]['perm_p_two_sided']):.4f} (16 types, {int(conc_fe.iloc[0]['n_common_genes'])} genes; unadj r={float(conc_un.iloc[0]['mean_pairwise_pearson_r']):.4f})",
            "mde80_summary": f"MDE80={float(conc_fe.iloc[0]['mde80_null']):.4f} (2.4865 * null_sd={float(conc_fe.iloc[0]['perm_null_sd']):.4f})",
        },
    ]
    for tier, fam_key in [("P1m", "transport_P1m"), ("P0m", "transport_P0m")]:
        df_ct, df_pl, _ = fam_c_results[tier]
        df_both_sum = pd.concat([df_ct.assign(is_pooled=False), df_pl.assign(is_pooled=True)], ignore_index=True)
        q_sum_c, _, _ = bh_checked(df_both_sum["composite_pvalue"])
        q_sum_e, _, _ = bh_checked(df_both_sum["empirical_null_pvalue"])
        both_sum = (q_sum_c < ALPHA) & (q_sum_e < ALPHA)
        bp = int(both_sum.sum())
        eo = int(((q_sum_c >= ALPHA) & (q_sum_e < ALPHA)).sum())
        co = int(((q_sum_c < ALPHA) & (q_sum_e >= ALPHA)).sum())
        bp_pos = int((both_sum & (df_both_sum["composite_effect"] > 0)).sum())
        pl_r = df_pl.iloc[0]
        mde_s = yazar_df[(yazar_df["endpoint_family"] == fam_key) & (yazar_df["verdict"] != "UNTESTABLE")]["mde80"]
        summary_rows.append(
            {
                "cohort": cohort_name,
                "endpoint_family": fam_key,
                "arm": "pool_fe",
                "donors_original": n_donors_orig,
                "donors_mature": n_donors_30,
                "donors_excluded": n_excluded_cohort,
                "original_headline": (
                    "5/16 positive pool-FE both-pass (All_Cells null emp_p=0.5287)"
                    if tier == "P1m"
                    else f"{int(((orig_p0_pool.drop(index=['All_Cells'], errors='ignore')['composite_q'] < ALPHA) & (orig_p0_pool.drop(index=['All_Cells'], errors='ignore')['empirical_q'] < ALPHA)).sum())}/16 both-pass (Block 5 P0_Raw257)"
                ),
                "mature_headline": f"{bp}/{len(df_both_sum)} both-pass ({bp_pos} positive), {eo} empirical-only, {co} composite-only (3 untestable); All_Cells eff={float(pl_r['composite_effect']):+.4f}, p_comp={float(pl_r['composite_pvalue']):.4g}, p_emp={float(pl_r['empirical_null_pvalue']):.4g}",
                "mde80_summary": f"median={float(mde_s.median()):.4f} (range {float(mde_s.min()):.4f} to {float(mde_s.max()):.4f}) per dec (2.4865 * null_sd)",
            }
        )

    save_cohort_payload(
        "yazar",
        {
            "summary_rows": summary_rows,
            "bh_rows": bh_rows,
            "diagnostics": {
                "n_donors_orig": n_donors_orig,
                "n_donors_30": n_donors_30,
                "n_excluded": n_excluded_cohort,
                "age_min": float(donor_meta_30["age_years"].min()),
                "age_max": float(donor_meta_30["age_years"].max()),
                "n_male": int(donor_meta_30["sex_male"].sum()),
                "n_female": int((1 - donor_meta_30["sex_male"]).sum()),
                "pool_stats": pool_stats,
                "pool_eligibility": pool_elig_df.reset_index().to_dict("records"),
                "unadj_clr_sig": int((fam_a_unadj["q_bh_family_a"] < ALPHA).sum()),
                "unadj_genes_sig": n_un_ct_sig,
                "unadj_genes_tested": n_un_ct_tests,
                "bh_diff_afe": bh_diff_afe,
                "bh_diff_b_ct": bh_diff_b_ct,
                "bh_diff_b_pl": bh_diff_b_pl,
            },
        },
    )
    append_runtime(cohort_name, "cohort_total", s0_all, t0_all, "ok", f"rows={len(yazar_df)}")
    print(f"[Yazar DONE] rows={len(yazar_df)}, elapsed={time.monotonic() - t0_all:.1f}s", flush=True)


# =============================================================================
# COHORT 3: PEREZ BLOOD CONTROLS (HCATA SLE/Control PBMC, GSE174188)
# =============================================================================

def prepare_perez_design_for_arm(
    meta: pd.DataFrame,
    arm: str,
    include_log_lib: bool,
    min_residual_df: int,
) -> tuple[pd.DataFrame, tuple[str, ...], str | None]:
    """Build Perez cohort parameterization and strictly test identifiability instead of silent dropping."""
    sub = meta.copy()
    sub["Processing_Cohort"] = sub["Processing_Cohort"].fillna("").astype(str)
    base_covs = ["sex_male", "log_lib"] if include_log_lib else ["sex_male"]

    if arm == "unadjusted":
        return sub, tuple(base_covs), None

    if arm == "female_only":
        sub = sub[sub["sex_male"] == 0].copy()
        covs = ("log_lib",) if include_log_lib else ()
        if len(sub) < 4:
            return sub, covs, f"female_only has {len(sub)} donors < 4"
        return sub, covs, None

    if arm == "cohort_multihot":
        tokens = sub["Processing_Cohort"].apply(lambda s: [t.strip() for t in s.split(";") if t.strip()])
        mh_cols = []
        for k in ("1.0", "2.0", "3.0", "4.0"):
            cname = "cohort_" + k.replace(".", "_")
            sub[cname] = tokens.apply(lambda ts, kk=k: float(kk in ts))
            mh_cols.append(cname)
        covs = tuple(base_covs + mh_cols)
        # Check strict identifiability: >=1 cohort token must enter non-collinearly, residual_df >= min_residual_df,
        # and age_per_decade must not be collinear with the cohort + base covariate span.
        try:
            design, names, omitted = select_design_columns(
                sub, age_col="age_per_decade", candidate_covariates=covs, min_residual_df=min_residual_df
            )
        except ValueError as exc:
            return sub, covs, f"cohort_multihot_unidentifiable: {exc}"
        active_cohort_cols = [c for c in names if c.startswith("cohort_")]
        if len(active_cohort_cols) == 0:
            return sub, covs, f"cohort_multihot_unidentifiable: zero active cohort indicators (omitted={omitted})"
        df_res = len(sub) - design.shape[1]
        if df_res < min_residual_df:
            return sub, covs, f"cohort_multihot_unidentifiable: residual_df={df_res} < {min_residual_df}"
        # Check if any cohort dummy was dropped due to df_limit (meaning cohort adjustment was silently truncated!)
        df_dropped = [o for o in omitted if o.startswith("cohort_") and "(df_limit)" in o]
        if df_dropped:
            return sub, covs, f"cohort_multihot_unidentifiable: cohort indicators dropped by df_limit ({df_dropped})"
        return sub, covs, None

    if arm == "cohort_level_ge2":
        vc = sub["Processing_Cohort"].value_counts()
        ge2_levels = sorted(vc[vc >= 2].index.tolist())
        dropped_singletons = {k: int(v) for k, v in vc[vc < 2].items()}
        sub = sub[sub["Processing_Cohort"].isin(ge2_levels)].copy()
        if len(ge2_levels) < 2:
            return (
                sub,
                (),
                f"cohort_level_ge2_unidentifiable: only {len(ge2_levels)} exact level(s) with >=2 donors (dropped_singletons={dropped_singletons})",
            )
        # Reference level is '4.0' if present (matching original), else largest remaining level
        ref_lvl = "4.0" if "4.0" in ge2_levels else str(vc.loc[ge2_levels].idxmax())
        lvl_cols = []
        for lv in ge2_levels:
            if lv == ref_lvl:
                continue
            cname = "lvl_" + lv.replace(".", "").replace(";", "_")
            sub[cname] = (sub["Processing_Cohort"] == lv).astype(float)
            lvl_cols.append(cname)
        covs = tuple(base_covs + lvl_cols)
        try:
            design, names, omitted = select_design_columns(
                sub, age_col="age_per_decade", candidate_covariates=covs, min_residual_df=min_residual_df
            )
        except ValueError as exc:
            return sub, covs, f"cohort_level_ge2_unidentifiable: {exc} (n={len(sub)}, levels={len(ge2_levels)}, dropped={dropped_singletons})"
        dropped_lvl = [o for o in omitted if o.startswith("lvl_")]
        if dropped_lvl:
            return (
                sub,
                covs,
                f"cohort_level_ge2_unidentifiable: level dummy dropped ({dropped_lvl}); n={len(sub)}, levels={ge2_levels}",
            )
        df_res = len(sub) - design.shape[1]
        if df_res < min_residual_df:
            return sub, covs, f"cohort_level_ge2_unidentifiable: residual_df={df_res} < {min_residual_df}"
        sub.attrs["ref_lvl"] = ref_lvl
        sub.attrs["dropped_singletons"] = dropped_singletons
        sub.attrs["n_levels_ge2"] = len(ge2_levels)
        return sub, covs, None

    raise ValueError(f"Unknown Perez arm {arm!r}")


def run_cohort_perez(programs: dict[str, dict[str, list[str]]]) -> None:
    cohort_name = "Perez_Blood_Controls"
    t0_all = time.monotonic()
    s0_all = utc_now()
    wait_for_memory("Perez_Blood_Controls start")

    print("\n================ COHORT 3: PEREZ BLOOD CONTROLS (GSE174188) ================", flush=True)
    perez_cache_pkl = TMP_WORK_DIR / "perez_30plus_raw_pb.pkl"

    if perez_cache_pkl.exists():
        print(f"[Perez] Loading verified 30-plus raw pseudobulk cache {perez_cache_pkl}", flush=True)
        with perez_cache_pkl.open("rb") as fh:
            pdata = pickle.load(fh)
        assert float(pdata["donor_meta_30"]["age_years"].min()) >= MIN_AGE
    else:
        t0_dl = time.monotonic()
        s0_dl = utc_now()
        tmp_h5 = TMP_WORK_DIR / "perez_218acb0f_tmp.h5ad"
        try:
            download_h5ad_verified(PEREZ_CZ_URL, tmp_h5, PEREZ_EXPECTED_BYTES)
            append_runtime(cohort_name, "download_h5ad", s0_dl, t0_dl, "ok", f"bytes={PEREZ_EXPECTED_BYTES}")
            wait_for_memory("Perez_Blood_Controls stream")
            t0_st = time.monotonic()
            s0_st = utc_now()
            pdata = stream_hcata_mature_from_h5ad(
                study="perez",
                h5_path=tmp_h5,
                expected_n_obs=1263676,
                expected_total_nnz=1000691133,
                expected_raw_source="raw/X",
                expected_orig_donors_csv=RESULTS_DIR / "followup-hcata-perez-donors.csv",
                min_donors_extract=4,
            )
            TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
            tmp_p = perez_cache_pkl.with_suffix(".tmp")
            with tmp_p.open("wb") as fh:
                pickle.dump(pdata, fh, protocol=pickle.HIGHEST_PROTOCOL)
            tmp_p.replace(perez_cache_pkl)
            append_runtime(
                cohort_name,
                "stream_h5ad_30plus",
                s0_st,
                t0_st,
                "ok",
                f"61 mature donors, cells_30={pdata['gate_30']['cells_30']}",
            )
        finally:
            if tmp_h5.exists():
                tmp_h5.unlink()
                print(f"[Perez] Deleted temporary H5AD {tmp_h5} to preserve disk space.", flush=True)

    donor_meta_30: pd.DataFrame = pdata["donor_meta_30"].copy()
    cell_counts_30: pd.DataFrame = pdata["cell_counts_30"].copy()
    raw_pb_by_group: dict[str, pd.DataFrame] = pdata["raw_pb_by_group"]
    n_donors_orig = int(pdata["n_donors_orig"])
    n_donors_30 = int(pdata["n_donors_30"])
    n_excluded_cohort = int(pdata["n_donors_excluded"])
    assert n_donors_orig == 99 and n_donors_30 == 61

    cohort_anova = anova_batch_age(donor_meta_30, "Processing_Cohort", single_only=True)

    # Finalize 30-plus pseudobulks
    expr_30: dict[str, pd.DataFrame] = {}
    meta_30: dict[str, pd.DataFrame] = {}
    for grp, raw_df in raw_pb_by_group.items():
        cols = list(raw_df.columns)
        g_counts = cell_counts_30.loc[cols].sum(axis=1) if grp == "All_Cells" else cell_counts_30.loc[cols, grp]
        l2, gm, _ = finalize_pseudobulk_group(raw_df, donor_meta_30.loc[cols], g_counts)
        expr_30[grp] = l2
        meta_30[grp] = gm
    del raw_pb_by_group, pdata
    gc.collect()

    # Load original Perez references
    orig_pool_frac = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-fractions.csv")
    orig_pool_frac = orig_pool_frac[orig_pool_frac["study"] == "perez"]
    orig_pool_gsum = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-genes-summary.csv")
    orig_pool_gsum = orig_pool_gsum[orig_pool_gsum["study"] == "perez"]
    orig_pool_conc = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-concordance.csv")
    orig_pool_conc = orig_pool_conc[orig_pool_conc["study"] == "perez"]
    orig_pool_p1 = pd.read_csv(RESULTS_DIR / "followup-hcata-pool-sensitivity-p1-signature.csv")
    orig_pool_p1 = orig_pool_p1[orig_pool_p1["study"] == "perez"]

    orig_unadj_frac = pd.read_csv(RESULTS_DIR / "followup-hcata-perez-fractions.csv")
    orig_unadj_frac = orig_unadj_frac[orig_unadj_frac["arm"] == "primary"].set_index("cell_type")
    orig_unadj_gsum = pd.read_csv(RESULTS_DIR / "followup-hcata-perez-genes-summary.csv").set_index("group")
    orig_unadj_conc = pd.read_csv(RESULTS_DIR / "followup-hcata-perez-concordance.csv")
    orig_unadj_conc = orig_unadj_conc[orig_unadj_conc["arm"] == "primary"].iloc[0]
    orig_unadj_p1 = pd.read_csv(RESULTS_DIR / "followup-hcata-perez-p1-signature.csv")
    orig_unadj_p1 = orig_unadj_p1[orig_unadj_p1["arm"] == "primary"].set_index("group")

    orig_clr_roster = sorted(orig_unadj_frac.index.tolist())
    orig_group_roster = sorted(orig_unadj_gsum.index.tolist())
    orig_ct_roster = [g for g in orig_group_roster if g != "All_Cells"]

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    arm_diag: dict[str, Any] = {}

    arms_spec = [
        ("cohort_multihot", 10, 3, 20, 4),
        ("cohort_level_ge2", 10, 3, 20, 4),
        ("unadjusted", 2, 6, 100, 6),
    ]

    for arm_name, min_df, min_pres_a, min_cells_a, min_don_bc in arms_spec:
        wait_for_memory(f"Perez {arm_name}")
        t0_arm = time.monotonic()
        s0_arm = utc_now()
        print(f"[Perez] Running arm={arm_name} (min_df={min_df})...", flush=True)

        if arm_name == "unadjusted":
            o_frac_df = orig_unadj_frac
            o_gsum_df = orig_unadj_gsum
            o_conc_row = orig_unadj_conc
            o_p1_df = orig_unadj_p1
        else:
            o_frac_df = orig_pool_frac[orig_pool_frac["arm"] == arm_name].set_index("cell_type")
            o_gsum_df = orig_pool_gsum[orig_pool_gsum["arm"] == arm_name].set_index("group")
            o_conc_row = orig_pool_conc[orig_pool_conc["arm"] == arm_name].iloc[0]
            o_p1_df = orig_pool_p1[orig_pool_p1["arm"] == arm_name].set_index("group")

        # Family A
        meta_a, covs_a, err_a = prepare_perez_design_for_arm(donor_meta_30, arm_name, include_log_lib=False, min_residual_df=min_df)
        if err_a is not None:
            for ct in orig_clr_roster:
                orow = o_frac_df.loc[ct]
                out_rows.append(
                    make_row(
                        cohort_name,
                        ct,
                        f"fractions_clr_{arm_name}",
                        int(len(meta_a)),
                        int(orow["n_donors"]) - int(len(meta_a)),
                        float(orow["pvalue"]),
                        float(orow["q_bh_family_a"]),
                        None,
                        None,
                        None,
                        None,
                        "UNTESTABLE",
                        err_a,
                        arm=arm_name,
                    )
                )
            fam_a = pd.DataFrame()
        else:
            fam_a = run_family_a_clr(
                cell_counts_30.loc[meta_a.index, orig_clr_roster],
                meta_a,
                age_col="age_per_decade",
                candidate_covariates=covs_a,
                min_donors_present=min_pres_a,
                min_total_cells=min_cells_a,
                min_residual_df=min_df,
            )
            q_a_exp, q_a_sm, diff_a = bh_checked(fam_a["pvalue"])
            fam_a["q_bh_family_a"] = q_a_exp
            fam_a_idx = fam_a.set_index("cell_type")
            pos_a = {str(c): i for i, c in enumerate(fam_a["cell_type"])}
            for ct in orig_clr_roster:
                orow = o_frac_df.loc[ct]
                if ct not in fam_a_idx.index:
                    out_rows.append(
                        make_row(
                            cohort_name,
                            ct,
                            f"fractions_clr_{arm_name}",
                            int(len(meta_a)),
                            int(orow["n_donors"]) - int(len(meta_a)),
                            float(orow["pvalue"]),
                            float(orow["q_bh_family_a"]),
                            None,
                            None,
                            None,
                            None,
                            "UNTESTABLE",
                            f"failed_clr_cell_or_donor_floor(min_present={min_pres_a},min_cells={min_cells_a})",
                            arm=arm_name,
                        )
                    )
                    continue
                r = fam_a_idx.loc[ct]
                pi = pos_a[ct]
                verd, vnote = verdict_standard(
                    float(orow["q_bh_family_a"]) < ALPHA,
                    float(r["q_bh_family_a"]) < ALPHA,
                    float(orow["effect_clr"]),
                    float(r["effect_clr"]),
                )
                mde = mde80_from_se(float(r["se"]), int(r["residual_df"]))
                out_rows.append(
                    make_row(
                        cohort_name,
                        ct,
                        f"fractions_clr_{arm_name}",
                        int(r["n_donors"]),
                        int(orow["n_donors"]) - int(r["n_donors"]),
                        float(orow["pvalue"]),
                        float(orow["q_bh_family_a"]),
                        float(r["pvalue"]),
                        float(r["q_bh_family_a"]),
                        None,
                        None,
                        verd,
                        f"n_present={int(r['n_donors_present'])};covs={r['covariates']};omitted={r['omitted_covariates']};{vnote}".strip(";"),
                        arm=arm_name,
                        residual_df=int(r["residual_df"]),
                        original_effect=float(orow["effect_clr"]),
                        mature_effect=float(r["effect_clr"]),
                        mature_se=float(r["se"]),
                        mde80=mde,
                        mde80_method="t_quantile_se",
                        covariates=str(r["covariates"]),
                        omitted_covariates=str(r["omitted_covariates"]),
                    )
                )
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"perez_{arm_name}_fractions_clr",
                        "unit": ct,
                        "family_size": int(len(fam_a)),
                        "raw_p": float(r["pvalue"]),
                        "q_bh_explicit": float(q_a_exp[pi]),
                        "q_bh_statsmodels": float(q_a_sm[pi]),
                        "abs_diff": float(abs(q_a_exp[pi] - q_a_sm[pi])),
                        "family_max_abs_diff": diff_a,
                    }
                )

        # Prepare per-group designs for Family B & C with explicit identifiability checks
        arm_expr: dict[str, pd.DataFrame] = {}
        arm_meta: dict[str, pd.DataFrame] = {}
        arm_covs_by_group: dict[str, tuple[str, ...]] = {}
        untestable_bc: dict[str, str] = {}

        for grp in orig_group_roster:
            if grp not in meta_30:
                untestable_bc[grp] = f"fewer than 4 mature donors with >= {MIN_CELLS_PER_DONOR} cells"
                continue
            m_g, c_g, err_g = prepare_perez_design_for_arm(meta_30[grp], arm_name, include_log_lib=True, min_residual_df=min_df)
            if err_g is not None:
                untestable_bc[grp] = err_g
                continue
            if len(m_g) < min_don_bc:
                untestable_bc[grp] = f"n_donors={len(m_g)} < min_donors={min_don_bc}"
                continue
            arm_meta[grp] = m_g
            arm_expr[grp] = expr_30[grp][list(m_g.index)]
            arm_covs_by_group[grp] = c_g

        # Family B: Because cohort_level_ge2 can have group-specific active level columns if a level has <2 donors in a group,
        # run OLS per group cleanly with its verified covariate tuple, then combine global cell-type BH and concordance!
        ct_genes_frames: list[pd.DataFrame] = []
        ct_sum_rows: list[dict[str, Any]] = []
        t_series_arm: dict[str, pd.Series] = {}
        designs_arm: dict[str, tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]] = {}
        log2_thresh = float(np.log2(2.0))

        for grp in orig_ct_roster:
            if grp in untestable_bc:
                continue
            mat_df = arm_expr[grp]
            m_df = arm_meta[grp]
            covs_g = arm_covs_by_group[grp]
            min_req = max(3, int(np.ceil(0.5 * len(m_df))))
            keep_m = (mat_df.to_numpy(dtype=np.float64) >= log2_thresh).sum(axis=1) >= min_req
            mat_sub = mat_df.loc[keep_m]
            design, names, omitted = select_design_columns(
                m_df, age_col="age_per_decade", candidate_covariates=covs_g, min_residual_df=min_df
            )
            outcomes = mat_sub.to_numpy(dtype=np.float64).T
            ols = fit_ols(outcomes, design, effect_column=1)
            t_stats = np.where(ols["se"] > 0, ols["effect"] / ols["se"], np.nan)
            t_series_arm[grp] = pd.Series(t_stats, index=mat_sub.index)
            designs_arm[grp] = (mat_sub, m_df, design, names)
            df_val = int(ols["residual_df"][0])
            mde_factor = float(scipy_stats.t.ppf(0.975, df_val) + scipy_stats.t.ppf(0.80, df_val))
            mde_vec = mde_factor * ols["se"]
            q_w, _, _ = bh_checked(ols["pvalue"])
            ct_genes_frames.append(
                pd.DataFrame(
                    {
                        "group": grp,
                        "gene": mat_sub.index.astype(str),
                        "pvalue": ols["pvalue"],
                        "effect": ols["effect"],
                        "se": ols["se"],
                        "q_bh_within_group": q_w,
                    }
                )
            )
            min_i = int(np.nanargmin(ols["pvalue"]))
            ct_sum_rows.append(
                {
                    "group": grp,
                    "n_donors": int(len(m_df)),
                    "n_genes_tested": int(len(mat_sub)),
                    "residual_df": df_val,
                    "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                    "omitted_covariates": ";".join(omitted) if omitted else "none",
                    "median_se": float(np.nanmedian(ols["se"])),
                    "median_mde80": float(np.nanmedian(mde_vec)),
                    "top_gene": str(mat_sub.index[min_i]),
                    "top_gene_effect": float(ols["effect"][min_i]),
                    "min_pvalue": float(ols["pvalue"][min_i]),
                }
            )

        if ct_genes_frames:
            all_ct_g = pd.concat(ct_genes_frames, ignore_index=True)
            q_g_exp, q_g_sm, diff_g_ct = bh_checked(all_ct_g["pvalue"])
            all_ct_g["q_bh_global"] = q_g_exp
            n_ct_tests_arm = int(len(all_ct_g))
            n_ct_sig_arm = int((q_g_exp < ALPHA).sum())
            sig_by_g = all_ct_g.groupby("group")["q_bh_global"].apply(lambda s: int((s < ALPHA).sum())).to_dict()
            minq_by_g = all_ct_g.groupby("group")["q_bh_global"].min().to_dict()
            ct_sum_df = pd.DataFrame(ct_sum_rows).set_index("group")
            for grp in ct_sum_df.index:
                sub_g = all_ct_g[all_ct_g["group"] == grp]
                idx_m = sub_g["pvalue"].idxmin()
                rm = all_ct_g.loc[idx_m]
                pm = int(all_ct_g.index.get_loc(idx_m))
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"perez_{arm_name}_genes_ct",
                        "unit": f"{grp}::{rm['gene']}",
                        "family_size": n_ct_tests_arm,
                        "raw_p": float(rm["pvalue"]),
                        "q_bh_explicit": float(q_g_exp[pm]),
                        "q_bh_statsmodels": float(q_g_sm[pm]),
                        "abs_diff": float(abs(q_g_exp[pm] - q_g_sm[pm])),
                        "family_max_abs_diff": diff_g_ct,
                    }
                )
        else:
            all_ct_g = pd.DataFrame()
            ct_sum_df = pd.DataFrame()
            n_ct_tests_arm = 0
            n_ct_sig_arm = 0
            sig_by_g = {}
            minq_by_g = {}
            diff_g_ct = 0.0

        # Pooled All_Cells genes
        pool_sig_arm = 0
        pool_tests_arm = 0
        pool_sum_rec = None
        if "All_Cells" not in untestable_bc:
            mat_df = arm_expr["All_Cells"]
            m_df = arm_meta["All_Cells"]
            covs_g = arm_covs_by_group["All_Cells"]
            min_req = max(3, int(np.ceil(0.5 * len(m_df))))
            keep_m = (mat_df.to_numpy(dtype=np.float64) >= log2_thresh).sum(axis=1) >= min_req
            mat_sub = mat_df.loc[keep_m]
            design, names, omitted = select_design_columns(
                m_df, age_col="age_per_decade", candidate_covariates=covs_g, min_residual_df=min_df
            )
            outcomes = mat_sub.to_numpy(dtype=np.float64).T
            ols = fit_ols(outcomes, design, effect_column=1)
            q_p_exp, q_p_sm, diff_p_g = bh_checked(ols["pvalue"])
            df_val = int(ols["residual_df"][0])
            mde_factor = float(scipy_stats.t.ppf(0.975, df_val) + scipy_stats.t.ppf(0.80, df_val))
            mde_vec = mde_factor * ols["se"]
            min_i = int(np.nanargmin(ols["pvalue"]))
            pool_tests_arm = int(len(mat_sub))
            pool_sig_arm = int((q_p_exp < ALPHA).sum())
            pool_sum_rec = {
                "n_donors": int(len(m_df)),
                "n_genes_tested": pool_tests_arm,
                "residual_df": df_val,
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "median_se": float(np.nanmedian(ols["se"])),
                "median_mde80": float(np.nanmedian(mde_vec)),
                "top_gene": str(mat_sub.index[min_i]),
                "top_gene_effect": float(ols["effect"][min_i]),
                "min_pvalue": float(ols["pvalue"][min_i]),
                "min_q": float(q_p_exp[min_i]),
            }
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": f"perez_{arm_name}_genes_all_cells",
                    "unit": f"All_Cells::{mat_sub.index[min_i]}",
                    "family_size": pool_tests_arm,
                    "raw_p": float(ols["pvalue"][min_i]),
                    "q_bh_explicit": float(q_p_exp[min_i]),
                    "q_bh_statsmodels": float(q_p_sm[min_i]),
                    "abs_diff": float(abs(q_p_exp[min_i] - q_p_sm[min_i])),
                    "family_max_abs_diff": diff_p_g,
                }
            )

        # Emit genes_celltype rows
        for grp in orig_group_roster:
            orow = o_gsum_df.loc[grp]
            orig_n = int(orow["n_donors"])
            orig_nsig = int(orow["n_sig_global_q05"])
            if grp in untestable_bc:
                out_rows.append(
                    make_row(
                        cohort_name,
                        grp,
                        f"genes_celltype_{arm_name}",
                        int(len(meta_30[grp])) if grp in meta_30 else 0,
                        orig_n - (int(len(meta_30[grp])) if grp in meta_30 else 0),
                        float(orow["min_pvalue"]),
                        float(orow["min_q_bh_global"]),
                        None,
                        None,
                        None,
                        None,
                        "UNTESTABLE",
                        untestable_bc[grp],
                        arm=arm_name,
                    )
                )
                continue
            if grp == "All_Cells":
                rec = pool_sum_rec
                mat_nsig = pool_sig_arm
                mat_minq = float(rec["min_q"])
            else:
                rec = ct_sum_df.loc[grp]
                mat_nsig = int(sig_by_g.get(grp, 0))
                mat_minq = float(minq_by_g[grp])
            verd, _ = verdict_standard(orig_nsig > 0, mat_nsig > 0)
            out_rows.append(
                make_row(
                    cohort_name,
                    grp,
                    f"genes_celltype_{arm_name}",
                    int(rec["n_donors"]),
                    orig_n - int(rec["n_donors"]),
                    float(orow["min_pvalue"]),
                    float(orow["min_q_bh_global"]),
                    float(rec["min_pvalue"]),
                    mat_minq,
                    None,
                    None,
                    verd,
                    f"n_genes={int(rec['n_genes_tested'])};orig_sig={orig_nsig};mature_sig={mat_nsig};top={rec['top_gene']}",
                    arm=arm_name,
                    residual_df=int(rec["residual_df"]),
                    mature_effect=float(rec["top_gene_effect"]),
                    mature_se=float(rec["median_se"]),
                    mde80=float(rec["median_mde80"]),
                    mde80_method="median_t_quantile_se",
                    n_genes_tested=int(rec["n_genes_tested"]),
                    n_sig_q05=mat_nsig,
                    covariates=str(rec["covariates"]),
                    omitted_covariates=str(rec["omitted_covariates"]),
                )
            )

        # Concordance
        if len(t_series_arm) >= 2:
            conc_arm = compute_concordance_fwl(
                t_series_arm, designs_arm, "age_per_decade", CONCORDANCE_PERMS, ORIG_NULL_SEED
            ).iloc[0]
            op_c = float(o_conc_row["perm_p_two_sided"])
            mp_c = float(conc_arm["perm_p_two_sided"])
            verd_c, _ = verdict_standard(
                op_c < ALPHA,
                mp_c < ALPHA,
                float(o_conc_row["mean_pairwise_pearson_r"]),
                float(conc_arm["mean_pairwise_pearson_r"]),
            )
            out_rows.append(
                make_row(
                    cohort_name,
                    f"All_{int(conc_arm['n_cell_types'])}_CellTypes",
                    f"concordance_{arm_name}",
                    int(conc_arm["n_donors_union"]),
                    int(o_gsum_df.loc["All_Cells", "n_donors"]) - int(conc_arm["n_donors_union"]),
                    op_c,
                    op_c,
                    mp_c,
                    mp_c,
                    mp_c,
                    mp_c,
                    verd_c,
                    f"n_ct={int(conc_arm['n_cell_types'])};n_pairs={int(conc_arm['n_pairs'])};n_common_genes={int(conc_arm['n_common_genes'])};median_r={float(conc_arm['median_pairwise_pearson_r']):.4f};untestable_cts={sorted([g for g in orig_ct_roster if g in untestable_bc])}",
                    arm=arm_name,
                    original_effect=float(o_conc_row["mean_pairwise_pearson_r"]),
                    mature_effect=float(conc_arm["mean_pairwise_pearson_r"]),
                    null_sd=float(conc_arm["perm_null_sd"]),
                    mde80=float(conc_arm["mde80_null"]),
                    mde80_method="2.4865_times_null_sd",
                )
            )
        else:
            conc_arm = None
            out_rows.append(
                make_row(
                    cohort_name,
                    "All_CellTypes",
                    f"concordance_{arm_name}",
                    0,
                    n_donors_orig,
                    float(o_conc_row["perm_p_two_sided"]),
                    float(o_conc_row["perm_p_two_sided"]),
                    None,
                    None,
                    None,
                    None,
                    "UNTESTABLE",
                    f"fewer than 2 testable cell types under {arm_name}",
                    arm=arm_name,
                )
            )

        # Family C: P1m and P0m per group (with group-specific verified covariates!)
        tier_stats: dict[str, dict[str, Any]] = {}
        for tier in ("P1m", "P0m"):
            fam_key = f"transport_{tier}_{arm_name}"
            ct_p_rows = []
            pl_p_row = None
            tier_skips = dict(untestable_bc)

            for grp in orig_group_roster:
                if grp in untestable_bc:
                    continue
                res_one, sk_one = run_family_c_mature(
                    {grp: arm_expr[grp]},
                    {grp: arm_meta[grp]},
                    program=programs[tier],
                    candidate_covariates=arm_covs_by_group[grp],
                    min_donors=min_don_bc,
                    min_residual_df=min_df,
                    iterations=P1_ITERATIONS,
                    seed=ORIG_NULL_SEED,
                    verify_first_parity=(grp == "B cell" and tier == "P1m"),
                )
                if res_one.empty:
                    tier_skips[grp] = sk_one.get(grp, "failed_family_c_gate")
                    continue
                if grp == "All_Cells":
                    pl_p_row = res_one.iloc[0].to_dict()
                else:
                    ct_p_rows.append(res_one.iloc[0].to_dict())

            df_ct_p = pd.DataFrame(ct_p_rows) if ct_p_rows else pd.DataFrame()
            if not df_ct_p.empty:
                q_cc_exp, q_cc_sm, diff_cc = bh_checked(df_ct_p["composite_pvalue"])
                q_ce_exp, q_ce_sm, diff_ce = bh_checked(df_ct_p["empirical_null_pvalue"])
                df_ct_p["q_bh_composite"] = q_cc_exp
                df_ct_p["q_bh_empirical"] = q_ce_exp
                df_ct_idx = df_ct_p.set_index("group")
                pos_ct_p = {str(g): i for i, g in enumerate(df_ct_p["group"])}
            else:
                df_ct_idx = pd.DataFrame()
                pos_ct_p = {}
                q_cc_exp = q_cc_sm = q_ce_exp = q_ce_sm = np.array([])
                diff_cc = diff_ce = 0.0

            # Also compute all-groups BH (cell types + All_Cells) cleanly
            all_p_list = list(ct_p_rows) + ([pl_p_row] if pl_p_row is not None else [])
            df_all_p = pd.DataFrame(all_p_list) if all_p_list else pd.DataFrame()
            if not df_all_p.empty:
                q_ac_exp, q_ac_sm, diff_ac = bh_checked(df_all_p["composite_pvalue"])
                q_ae_exp, q_ae_sm, diff_ae = bh_checked(df_all_p["empirical_null_pvalue"])
                pos_all_p = {str(g): i for i, g in enumerate(df_all_p["group"])}
            else:
                q_ac_exp = q_ac_sm = q_ae_exp = q_ae_sm = np.array([])
                diff_ac = diff_ae = 0.0
                pos_all_p = {}

            for grp in orig_group_roster:
                oref = o_p1_df.loc[grp]
                orig_n = int(oref["n_donors"])
                is_pool = grp == "All_Cells"
                if tier == "P1m":
                    op_c = float(oref["composite_pvalue"])
                    oq_c = float(oref["q_bh_composite"]) if (not is_pool and pd.notna(oref["q_bh_composite"])) else op_c
                    op_e = float(oref["empirical_null_pvalue"])
                    oq_e = float(oref["q_bh_empirical"]) if (not is_pool and pd.notna(oref["q_bh_empirical"])) else op_e
                    oeff = float(oref["composite_effect"])
                    orig_cls = transport_class(oq_c, oq_e, True)
                    orig_sig_flag: bool | None = orig_cls == "BOTH_PASS"
                else:
                    op_c = oq_c = op_e = oq_e = oeff = np.nan
                    orig_cls = "NO_ORIG_P0"
                    orig_sig_flag = None

                if grp in tier_skips or (not is_pool and grp not in df_ct_idx.index) or (is_pool and pl_p_row is None):
                    reason = tier_skips.get(grp, "untestable_at_30plus")
                    out_rows.append(
                        make_row(
                            cohort_name,
                            grp,
                            fam_key,
                            int(len(meta_30[grp])) if grp in meta_30 else 0,
                            orig_n - (int(len(meta_30[grp])) if grp in meta_30 else 0),
                            op_c if np.isfinite(op_c) else None,
                            oq_c if np.isfinite(oq_c) else None,
                            None,
                            None,
                            None,
                            None,
                            "UNTESTABLE",
                            f"{reason};orig_class={orig_cls}",
                            arm=arm_name,
                            orig_transport_class=orig_cls,
                            mature_transport_class="UNTESTABLE",
                        )
                    )
                    continue

                r = pl_p_row if is_pool else df_ct_idx.loc[grp]
                mat_qc = float(q_ac_exp[pos_all_p[grp]])
                mat_qe = float(q_ae_exp[pos_all_p[grp]])
                mat_cls = transport_class(mat_qc, mat_qe, True)
                verd, vnote = verdict_standard(
                    orig_sig_flag,
                    mat_cls == "BOTH_PASS",
                    oeff if np.isfinite(oeff) else None,
                    float(r["composite_effect"]),
                )
                mde_null = mde80_from_null_sd(float(r["null_sd_effect"]))
                mde_se = mde80_from_se(float(r["composite_se"]), int(r["residual_df"]))
                pa = pos_all_p[grp]
                notes = (
                    f"orig_class={orig_cls};mature_class={mat_cls};"
                    f"q_all_comp={float(q_ac_exp[pa]):.4g};q_all_emp={float(q_ae_exp[pa]):.4g};"
                    f"mde80_se={mde_se:.4f};draws={P1_ITERATIONS};seed={ORIG_NULL_SEED};"
                    f"min_decile_margin={int(r['min_decile_margin'])}"
                )
                if vnote:
                    notes += f";{vnote}"
                out_rows.append(
                    make_row(
                        cohort_name,
                        grp,
                        fam_key,
                        int(r["n_donors"]),
                        orig_n - int(r["n_donors"]),
                        op_c if np.isfinite(op_c) else None,
                        oq_c if np.isfinite(oq_c) else None,
                        float(r["composite_pvalue"]),
                        mat_qc,
                        float(r["empirical_null_pvalue"]),
                        mat_qe,
                        verd,
                        notes,
                        arm=arm_name,
                        residual_df=int(r["residual_df"]),
                        original_effect=oeff,
                        mature_effect=float(r["composite_effect"]),
                        mature_se=float(r["composite_se"]),
                        null_sd=float(r["null_sd_effect"]),
                        mde80=mde_null,
                        mde80_method="2.4865_times_null_sd",
                        original_empirical_p=op_e,
                        original_empirical_q=oq_e,
                        orig_transport_class=orig_cls,
                        mature_transport_class=mat_cls,
                        n_up_detected=int(r["n_up_detected"]),
                        n_down_detected=int(r["n_down_detected"]),
                        covariates=str(r["covariates"]),
                        omitted_covariates=str(r["omitted_covariates"]),
                    )
                )
                if not is_pool:
                    pi = pos_ct_p[grp]
                    bh_rows.append(
                        {
                            "cohort": cohort_name,
                            "family": f"perez_{arm_name}_{tier}_ct_composite",
                            "unit": grp,
                            "family_size": int(len(df_ct_p)),
                            "raw_p": float(r["composite_pvalue"]),
                            "q_bh_explicit": float(q_cc_exp[pi]),
                            "q_bh_statsmodels": float(q_cc_sm[pi]),
                            "abs_diff": float(abs(q_cc_exp[pi] - q_cc_sm[pi])),
                            "family_max_abs_diff": diff_cc,
                        }
                    )
                    bh_rows.append(
                        {
                            "cohort": cohort_name,
                            "family": f"perez_{arm_name}_{tier}_ct_empirical",
                            "unit": grp,
                            "family_size": int(len(df_ct_p)),
                            "raw_p": float(r["empirical_null_pvalue"]),
                            "q_bh_explicit": float(q_ce_exp[pi]),
                            "q_bh_statsmodels": float(q_ce_sm[pi]),
                            "abs_diff": float(abs(q_ce_exp[pi] - q_ce_sm[pi])),
                            "family_max_abs_diff": diff_ce,
                        }
                    )
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"perez_{arm_name}_{tier}_all_composite",
                        "unit": grp,
                        "family_size": int(len(df_all_p)),
                        "raw_p": float(r["composite_pvalue"]),
                        "q_bh_explicit": float(q_ac_exp[pa]),
                        "q_bh_statsmodels": float(q_ac_sm[pa]),
                        "abs_diff": float(abs(q_ac_exp[pa] - q_ac_sm[pa])),
                        "family_max_abs_diff": diff_ac,
                    }
                )
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"perez_{arm_name}_{tier}_all_empirical",
                        "unit": grp,
                        "family_size": int(len(df_all_p)),
                        "raw_p": float(r["empirical_null_pvalue"]),
                        "q_bh_explicit": float(q_ae_exp[pa]),
                        "q_bh_statsmodels": float(q_ae_sm[pa]),
                        "abs_diff": float(abs(q_ae_exp[pa] - q_ae_sm[pa])),
                        "family_max_abs_diff": diff_ae,
                    }
                )

            bp = int(((q_ac_exp < ALPHA) & (q_ae_exp < ALPHA)).sum()) if not df_all_p.empty else 0
            eo = int(((q_ac_exp >= ALPHA) & (q_ae_exp < ALPHA)).sum()) if not df_all_p.empty else 0
            co = int(((q_ac_exp < ALPHA) & (q_ae_exp >= ALPHA)).sum()) if not df_all_p.empty else 0
            tier_stats[tier] = {
                "n_ct_tested": int(len(df_ct_p)),
                "n_ct_untestable": int(len(orig_ct_roster) - len(df_ct_p)),
                "both_pass": bp,
                "emp_only": eo,
                "comp_only": co,
                "all_cells": pl_p_row,
            }

        # Summary rows for this Perez arm
        n_a_sig = int((fam_a["q_bh_family_a"] < ALPHA).sum()) if not fam_a.empty else 0
        n_a_orig_sig = int((o_frac_df["q_bh_family_a"] < ALPHA).sum())
        summary_rows.append(
            {
                "cohort": cohort_name,
                "endpoint_family": f"fractions_clr_{arm_name}",
                "arm": arm_name,
                "donors_original": int(o_frac_df["n_donors"].iloc[0]),
                "donors_mature": int(len(meta_a)),
                "donors_excluded": int(o_frac_df["n_donors"].iloc[0]) - int(len(meta_a)),
                "original_headline": f"{n_a_orig_sig}/{len(o_frac_df)} CLR hits (min q={float(o_frac_df['q_bh_family_a'].min()):.4f})",
                "mature_headline": (
                    f"UNTESTABLE ({err_a})"
                    if err_a is not None
                    else f"{n_a_sig}/{len(fam_a)} CLR hits (min q={float(fam_a['q_bh_family_a'].min()):.4f}, top={fam_a.iloc[0]['cell_type']} p={float(fam_a.iloc[0]['pvalue']):.4g})"
                ),
                "mde80_summary": (
                    "UNTESTABLE"
                    if fam_a.empty
                    else f"median={float(pd.Series([mde80_from_se(float(r.se), int(r.residual_df)) for r in fam_a.itertuples()]).median()):.4f} CLR/dec (t_quantile_se)"
                ),
            }
        )
        orig_ct_nsig_tot = int(o_gsum_df.drop(index=["All_Cells"], errors="ignore")["n_sig_global_q05"].sum())
        orig_ct_ntest_tot = int(o_gsum_df.drop(index=["All_Cells"], errors="ignore")["n_genes_tested"].sum())
        summary_rows.append(
            {
                "cohort": cohort_name,
                "endpoint_family": f"genes_{arm_name}",
                "arm": arm_name,
                "donors_original": int(o_gsum_df.loc["All_Cells", "n_donors"]),
                "donors_mature": int(pool_sum_rec["n_donors"]) if pool_sum_rec else 0,
                "donors_excluded": int(o_gsum_df.loc["All_Cells", "n_donors"]) - (int(pool_sum_rec["n_donors"]) if pool_sum_rec else 0),
                "original_headline": f"{orig_ct_nsig_tot}/{orig_ct_ntest_tot} ct hits; {int(o_gsum_df.loc['All_Cells', 'n_sig_global_q05'])} All_Cells",
                "mature_headline": f"{n_ct_sig_arm}/{n_ct_tests_arm} ct hits ({len(ct_sum_df)}/9 types tested, {9 - len(ct_sum_df)} untestable); {pool_sig_arm}/{pool_tests_arm} All_Cells",
                "mde80_summary": (
                    f"median cell-type median MDE80={float(ct_sum_df['median_mde80'].median()):.4f} log2CPM/dec (t_quantile_se)"
                    if not ct_sum_df.empty
                    else "UNTESTABLE"
                ),
            }
        )
        summary_rows.append(
            {
                "cohort": cohort_name,
                "endpoint_family": f"concordance_{arm_name}",
                "arm": arm_name,
                "donors_original": int(o_gsum_df.loc["All_Cells", "n_donors"]),
                "donors_mature": int(conc_arm["n_donors_union"]) if conc_arm is not None else 0,
                "donors_excluded": int(o_gsum_df.loc["All_Cells", "n_donors"]) - (int(conc_arm["n_donors_union"]) if conc_arm is not None else 0),
                "original_headline": f"mean r={float(o_conc_row['mean_pairwise_pearson_r']):.4f}, perm_p={float(o_conc_row['perm_p_two_sided']):.4f}",
                "mature_headline": (
                    "UNTESTABLE"
                    if conc_arm is None
                    else f"mean r={float(conc_arm['mean_pairwise_pearson_r']):.4f}, perm_p={float(conc_arm['perm_p_two_sided']):.4f} ({int(conc_arm['n_cell_types'])} types, {int(conc_arm['n_common_genes'])} genes)"
                ),
                "mde80_summary": (
                    "UNTESTABLE"
                    if conc_arm is None
                    else f"MDE80={float(conc_arm['mde80_null']):.4f} (2.4865 * null_sd={float(conc_arm['perm_null_sd']):.4f})"
                ),
            }
        )
        for tier in ("P1m", "P0m"):
            ts = tier_stats[tier]
            ac = ts["all_cells"]
            ac_str = (
                f"All_Cells eff={float(ac['composite_effect']):+.4f}, p_comp={float(ac['composite_pvalue']):.4g}, p_emp={float(ac['empirical_null_pvalue']):.4g}"
                if ac is not None
                else "All_Cells UNTESTABLE"
            )
            sub_mde = [
                r["mde80"]
                for r in out_rows
                if r["endpoint_family"] == f"transport_{tier}_{arm_name}" and r["verdict"] != "UNTESTABLE"
            ]
            summary_rows.append(
                {
                    "cohort": cohort_name,
                    "endpoint_family": f"transport_{tier}_{arm_name}",
                    "arm": arm_name,
                    "donors_original": int(o_p1_df.loc["All_Cells", "n_donors"]),
                    "donors_mature": int(ac["n_donors"]) if ac is not None else 0,
                    "donors_excluded": int(o_p1_df.loc["All_Cells", "n_donors"]) - (int(ac["n_donors"]) if ac is not None else 0),
                    "original_headline": (
                        f"0/9 per-type empirical passes; All_Cells eff={float(o_p1_df.loc['All_Cells', 'composite_effect']):+.4f}, emp_p={float(o_p1_df.loc['All_Cells', 'empirical_null_pvalue']):.4g}"
                        if tier == "P1m"
                        else "No original P0 run"
                    ),
                    "mature_headline": f"{ts['both_pass']}/{ts['n_ct_tested'] + (1 if ts['all_cells'] is not None else 0)} both-pass, {ts['emp_only']} emp-only, {ts['comp_only']} comp-only ({ts['n_ct_untestable']} untestable); {ac_str}",
                    "mde80_summary": (
                        f"median={float(np.nanmedian(sub_mde)):.4f} (range {float(np.nanmin(sub_mde)):.4f} to {float(np.nanmax(sub_mde)):.4f}) per dec (2.4865 * null_sd)"
                        if sub_mde
                        else "UNTESTABLE"
                    ),
                }
            )
        arm_diag[arm_name] = {
            "err_a": err_a,
            "untestable_bc": untestable_bc,
            "n_ct_tests": n_ct_tests_arm,
            "n_ct_sig": n_ct_sig_arm,
        }
        append_runtime(cohort_name, f"arm_{arm_name}", s0_arm, t0_arm, "ok", f"ct_sig={n_ct_sig_arm}, untestable_bc={len(untestable_bc)}")

    perez_df = order_cohort_df(out_rows)
    write_clean_csv(perez_df, RESULTS_DIR / f"{PREFIX}-perez.csv")

    save_cohort_payload(
        "perez",
        {
            "summary_rows": summary_rows,
            "bh_rows": bh_rows,
            "diagnostics": {
                "n_donors_orig": n_donors_orig,
                "n_donors_30": n_donors_30,
                "n_excluded": n_excluded_cohort,
                "age_min": float(donor_meta_30["age_years"].min()),
                "age_max": float(donor_meta_30["age_years"].max()),
                "n_male": int(donor_meta_30["sex_male"].sum()),
                "n_female": int((1 - donor_meta_30["sex_male"]).sum()),
                "cohort_counts": donor_meta_30["Processing_Cohort"].value_counts().to_dict(),
                "cohort_anova": cohort_anova,
                "arms": arm_diag,
            },
        },
    )
    append_runtime(cohort_name, "cohort_total", s0_all, t0_all, "ok", f"rows={len(perez_df)}")
    print(f"[Perez DONE] rows={len(perez_df)}, elapsed={time.monotonic() - t0_all:.1f}s", flush=True)


# =============================================================================
# COHORT 4: GUILLIAMS LIVER SORTED CELLS (HCATA GSE192742)
# =============================================================================

def run_cohort_guilliams(programs: dict[str, dict[str, list[str]]]) -> None:
    cohort_name = "Guilliams_Liver_Sorted"
    t0_all = time.monotonic()
    s0_all = utc_now()
    wait_for_memory("Guilliams_Liver_Sorted start")

    print("\n================ COHORT 4: GUILLIAMS LIVER SORTED CELLS (GSE192742) ================", flush=True)
    gui_cache_pkl = TMP_WORK_DIR / "guilliams_30plus_raw_pb.pkl"

    if gui_cache_pkl.exists():
        print(f"[Guilliams] Loading verified 30-plus raw pseudobulk cache {gui_cache_pkl}", flush=True)
        with gui_cache_pkl.open("rb") as fh:
            gdata = pickle.load(fh)
        assert float(gdata["donor_meta_30"]["age_years"].min()) >= MIN_AGE
    else:
        t0_dl = time.monotonic()
        s0_dl = utc_now()
        tmp_h5 = TMP_WORK_DIR / "guilliams_e84f2780_tmp.h5ad"
        try:
            download_h5ad_verified(GUILLIAMS_CZ_URL, tmp_h5, GUILLIAMS_EXPECTED_BYTES)
            append_runtime(cohort_name, "download_h5ad", s0_dl, t0_dl, "ok", f"bytes={GUILLIAMS_EXPECTED_BYTES}")
            wait_for_memory("Guilliams_Liver_Sorted stream")
            t0_st = time.monotonic()
            s0_st = utc_now()
            gdata = stream_hcata_mature_from_h5ad(
                study="guilliams",
                h5_path=tmp_h5,
                expected_n_obs=167598,
                expected_total_nnz=215769434,
                expected_raw_source="raw/X",
                expected_orig_donors_csv=RESULTS_DIR / "followup-hcata-guilliams-donors.csv",
                min_donors_extract=4,
            )
            TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
            tmp_p = gui_cache_pkl.with_suffix(".tmp")
            with tmp_p.open("wb") as fh:
                pickle.dump(gdata, fh, protocol=pickle.HIGHEST_PROTOCOL)
            tmp_p.replace(gui_cache_pkl)
            append_runtime(
                cohort_name,
                "stream_h5ad_30plus",
                s0_st,
                t0_st,
                "ok",
                f"mature_donors={gdata['n_donors_30']}, cells_30={gdata['gate_30']['cells_30']}",
            )
        finally:
            if tmp_h5.exists():
                tmp_h5.unlink()
                print(f"[Guilliams] Deleted temporary H5AD {tmp_h5} to preserve disk space.", flush=True)

    donor_meta_30: pd.DataFrame = gdata["donor_meta_30"].copy()
    donor_meta_30["is_citeseq"] = donor_meta_30["typeSample"].astype(str).str.contains("citeSeq").astype(float)
    cell_counts_30: pd.DataFrame = gdata["cell_counts_30"].copy()
    cd45_counts_30: pd.Series = gdata["cd45_counts_30"].copy()
    raw_pb_by_group: dict[str, pd.DataFrame] = gdata["raw_pb_by_group"]
    n_donors_orig = int(gdata["n_donors_orig"])
    n_donors_30 = int(gdata["n_donors_30"])
    n_excluded_cohort = int(gdata["n_donors_excluded"])
    full_obj_donors_total = int(gdata["full_object_donors_total"])
    full_obj_under30 = gdata["full_object_under30_donors"]
    assert n_donors_orig == 10 and n_donors_30 == 9 and n_excluded_cohort == 1

    # Finalize 30-plus pseudobulks (keeping all >=4-donor groups so we can report exact donor counts for <6 untestable!)
    expr_30: dict[str, pd.DataFrame] = {}
    meta_30: dict[str, pd.DataFrame] = {}
    for grp, raw_df in raw_pb_by_group.items():
        cols = list(raw_df.columns)
        if grp == "All_Cells":
            g_counts = donor_meta_30.loc[cols, "n_cells"].astype(int)
        elif grp == "All_CD45pos_Cells":
            g_counts = cd45_counts_30.loc[cols].astype(int)
        else:
            g_counts = cell_counts_30.loc[cols, grp].astype(int)
        l2, gm, _ = finalize_pseudobulk_group(raw_df, donor_meta_30.loc[cols], g_counts)
        expr_30[grp] = l2
        meta_30[grp] = gm
    del raw_pb_by_group
    gc.collect()

    orig_gsum = pd.read_csv(RESULTS_DIR / "followup-hcata-guilliams-genes-summary.csv").set_index("group")
    orig_p1 = pd.read_csv(RESULTS_DIR / "followup-hcata-guilliams-p1-signature.csv").set_index("group")
    orig_sens_g = pd.read_csv(RESULTS_DIR / "followup-hcata-guilliams-sensitivities-genes.csv")
    orig_sens_p = pd.read_csv(RESULTS_DIR / "followup-hcata-guilliams-sensitivities-p1.csv")

    orig_ct_roster = sorted([g for g in orig_gsum.index if not str(g).startswith("All_")])
    orig_pool_roster = sorted([g for g in orig_gsum.index if str(g).startswith("All_")])
    orig_all_roster = orig_ct_roster + orig_pool_roster

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    # Explicitly document sorted-cell exclusion for Family A and concordance
    out_rows.append(
        make_row(
            cohort_name,
            "All_Sorted_Cells",
            "fractions_clr",
            n_donors_30,
            n_excluded_cohort,
            None,
            None,
            None,
            None,
            None,
            None,
            "UNTESTABLE",
            "FACS_sorted_CD45pos_CD45neg_enrichment_distorts_whole_liver_proportions;not_run_in_original",
            arm="primary",
        )
    )
    out_rows.append(
        make_row(
            cohort_name,
            "All_Sorted_Cells",
            "concordance",
            n_donors_30,
            n_excluded_cohort,
            None,
            None,
            None,
            None,
            None,
            None,
            "UNTESTABLE",
            "FACS_sorted_exclusion;no_cross_celltype_concordance_in_original",
            arm="primary",
        )
    )

    # Primary + 4 original sensitivity arms (using exact original MIN_DONORS = 6 and min_residual_df = 2!)
    all_ids = set(donor_meta_30.index)
    arm_specs = [
        ("primary", all_ids, ("sex_male", "log_lib")),
        ("no_malignant_resection", all_ids - {"H16", "H21"}, ("sex_male", "log_lib")),
        ("no_any_tumor", all_ids - {"H13", "H16", "H21"}, ("log_lib",)),
        ("female_only", set(donor_meta_30.index[donor_meta_30["sex_male"] == 0]), ("log_lib",)),
        ("adjust_citeseq", all_ids, ("sex_male", "log_lib", "is_citeseq")),
    ]
    orig_sens_name_map = {
        "no_malignant_resection": "no_malignant_resection_n8",
        "no_any_tumor": "no_any_tumor_n7",
        "female_only": "female_only_n8",
        "adjust_citeseq": "adjust_citeseq_n10",
    }

    arm_diagnostics: dict[str, Any] = {}
    for arm_name, dset, covs in arm_specs:
        t0_arm = time.monotonic()
        s0_arm = utc_now()
        suffix = "" if arm_name == "primary" else f"_{arm_name}"

        if arm_name == "primary":
            o_g = orig_gsum
            o_p = orig_p1
        else:
            osn = orig_sens_name_map[arm_name]
            o_g = orig_sens_g[orig_sens_g["sensitivity"] == osn].set_index("group")
            o_p = orig_sens_p[orig_sens_p["sensitivity"] == osn].set_index("group")

        ct_e: dict[str, pd.DataFrame] = {}
        ct_m: dict[str, pd.DataFrame] = {}
        pl_e: dict[str, pd.DataFrame] = {}
        pl_m: dict[str, pd.DataFrame] = {}
        elig_skip: dict[str, str] = {}

        for grp in orig_all_roster:
            if grp not in expr_30:
                n_here = int((cell_counts_30.loc[list(dset), grp] >= MIN_CELLS_PER_DONOR).sum()) if grp in cell_counts_30.columns else 0
                elig_skip[grp] = f"n_donors={n_here} < MIN_DONORS=6"
                continue
            cols = [d for d in expr_30[grp].columns if d in dset]
            if len(cols) < 6:
                elig_skip[grp] = f"n_donors={len(cols)} < MIN_DONORS=6"
                continue
            if grp.startswith("All_"):
                pl_e[grp] = expr_30[grp][cols]
                pl_m[grp] = meta_30[grp].loc[cols]
            else:
                ct_e[grp] = expr_30[grp][cols]
                ct_m[grp] = meta_30[grp].loc[cols]

        # Family B
        g_ct, s_ct, _ = run_family_b_mature(
            ct_e, ct_m, candidate_covariates=covs, min_donors=6, min_residual_df=2, run_concordance=False, seed=ORIG_NULL_SEED
        )
        g_pl, s_pl, _ = run_family_b_mature(
            pl_e, pl_m, candidate_covariates=covs, min_donors=6, min_residual_df=2, run_concordance=False, seed=ORIG_NULL_SEED
        )

        q_ct_exp, q_ct_sm, diff_ct = bh_checked(g_ct["pvalue"]) if not g_ct.empty else (np.array([]), np.array([]), 0.0)
        q_pl_exp, q_pl_sm, diff_pl = bh_checked(g_pl["pvalue"]) if not g_pl.empty else (np.array([]), np.array([]), 0.0)

        s_ct_idx = s_ct.set_index("group") if not s_ct.empty else pd.DataFrame()
        s_pl_idx = s_pl.set_index("group") if not s_pl.empty else pd.DataFrame()

        for grp in orig_all_roster:
            is_pool = grp.startswith("All_")
            orig_present = grp in o_g.index
            orig_n = int(o_g.loc[grp, "n_donors"]) if orig_present else int(orig_gsum.loc[grp, "n_donors"])
            orig_nsig = int(o_g.loc[grp, "n_sig_global_q05"]) if orig_present else 0
            orig_minp = float(o_g.loc[grp, "min_pvalue"]) if orig_present else np.nan
            orig_minq = float(o_g.loc[grp, "min_q_bh_global"]) if orig_present else np.nan

            s_idx = s_pl_idx if is_pool else s_ct_idx
            if grp in elig_skip or grp not in s_idx.index or str(s_idx.loc[grp, "status"]) != "tested":
                reason = elig_skip.get(grp, str(s_idx.loc[grp, "status"]) if grp in s_idx.index else "untestable")
                n_mat = int((cell_counts_30.loc[list(dset), grp] >= MIN_CELLS_PER_DONOR).sum()) if grp in cell_counts_30.columns else 0
                out_rows.append(
                    make_row(
                        cohort_name,
                        grp,
                        f"genes_celltype{suffix}",
                        n_mat,
                        orig_n - n_mat,
                        orig_minp if np.isfinite(orig_minp) else None,
                        orig_minq if np.isfinite(orig_minq) else None,
                        None,
                        None,
                        None,
                        None,
                        "UNTESTABLE",
                        reason,
                        arm=arm_name,
                    )
                )
                continue

            r = s_idx.loc[grp]
            mat_nsig = int(r["n_sig_global_q05"])
            verd, _ = verdict_standard(orig_nsig > 0 if orig_present else None, mat_nsig > 0)
            out_rows.append(
                make_row(
                    cohort_name,
                    grp,
                    f"genes_celltype{suffix}",
                    int(r["n_donors"]),
                    orig_n - int(r["n_donors"]),
                    orig_minp if np.isfinite(orig_minp) else None,
                    orig_minq if np.isfinite(orig_minq) else None,
                    float(r["min_pvalue"]),
                    float(r["min_q_bh_global"]),
                    None,
                    None,
                    verd,
                    f"n_genes={int(r['n_genes_tested'])};orig_sig={orig_nsig};mature_sig={mat_nsig};top={r['top_gene']}",
                    arm=arm_name,
                    residual_df=int(r["residual_df"]),
                    mature_effect=float(r["top_gene_effect"]),
                    mature_se=float(r["median_se"]),
                    mde80=float(r["median_mde80"]),
                    mde80_method="median_t_quantile_se",
                    n_genes_tested=int(r["n_genes_tested"]),
                    n_sig_q05=mat_nsig,
                    covariates=str(r["covariates"]),
                    omitted_covariates=str(r["omitted_covariates"]),
                )
            )
            g_src = g_pl if is_pool else g_ct
            q_e_src = q_pl_exp if is_pool else q_ct_exp
            q_s_src = q_pl_sm if is_pool else q_ct_sm
            diff_src = diff_pl if is_pool else diff_ct
            sub_g = g_src[g_src["group"] == grp]
            idx_m = sub_g["pvalue"].idxmin()
            rm = g_src.loc[idx_m]
            pm = int(g_src.index.get_loc(idx_m))
            bh_rows.append(
                {
                    "cohort": cohort_name,
                    "family": f"guilliams_{arm_name}_{'pool' if is_pool else 'ct'}_genes",
                    "unit": f"{grp}::{rm['gene']}",
                    "family_size": int(len(g_src)),
                    "raw_p": float(rm["pvalue"]),
                    "q_bh_explicit": float(q_e_src[pm]),
                    "q_bh_statsmodels": float(q_s_src[pm]),
                    "abs_diff": float(abs(q_e_src[pm] - q_s_src[pm])),
                    "family_max_abs_diff": diff_src,
                }
            )
            # Also log any mature gene hit with q < 0.05 in bh-check
            hits_g = sub_g[sub_g["q_bh_global"] < ALPHA]
            for h_idx, h_row in hits_g.iterrows():
                if h_idx == idx_m:
                    continue
                hp = int(g_src.index.get_loc(h_idx))
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"guilliams_{arm_name}_{'pool' if is_pool else 'ct'}_genes",
                        "unit": f"{grp}::{h_row['gene']}",
                        "family_size": int(len(g_src)),
                        "raw_p": float(h_row["pvalue"]),
                        "q_bh_explicit": float(q_e_src[hp]),
                        "q_bh_statsmodels": float(q_s_src[hp]),
                        "abs_diff": float(abs(q_e_src[hp] - q_s_src[hp])),
                        "family_max_abs_diff": diff_src,
                    }
                )

        # Family C: P1m and P0m
        tier_res: dict[str, Any] = {}
        for tier in ("P1m", "P0m"):
            fam_key = f"transport_{tier}{suffix}"
            df_c_ct, sk_c_ct = run_family_c_mature(
                ct_e,
                ct_m,
                program=programs[tier],
                candidate_covariates=covs,
                min_donors=6,
                min_residual_df=2,
                iterations=P1_ITERATIONS,
                seed=ORIG_NULL_SEED,
                verify_first_parity=(arm_name == "primary" and tier == "P1m"),
            )
            df_c_pl, sk_c_pl = run_family_c_mature(
                pl_e,
                pl_m,
                program=programs[tier],
                candidate_covariates=covs,
                min_donors=6,
                min_residual_df=2,
                iterations=P1_ITERATIONS,
                seed=ORIG_NULL_SEED,
                verify_first_parity=False,
            )
            sk_all = {**elig_skip, **sk_c_ct, **sk_c_pl}
            q_cc_e, q_cc_s, diff_cc = bh_checked(df_c_ct["composite_pvalue"]) if not df_c_ct.empty else (np.array([]), np.array([]), 0.0)
            q_ce_e, q_ce_s, diff_ce = bh_checked(df_c_ct["empirical_null_pvalue"]) if not df_c_ct.empty else (np.array([]), np.array([]), 0.0)
            df_all_c = pd.concat([df_c_ct.assign(is_pooled=False), df_c_pl.assign(is_pooled=True)], ignore_index=True)
            q_ac_e, q_ac_s, diff_ac = bh_checked(df_all_c["composite_pvalue"]) if not df_all_c.empty else (np.array([]), np.array([]), 0.0)
            q_ae_e, q_ae_s, diff_ae = bh_checked(df_all_c["empirical_null_pvalue"]) if not df_all_c.empty else (np.array([]), np.array([]), 0.0)

            df_ct_idx = df_c_ct.set_index("group") if not df_c_ct.empty else pd.DataFrame()
            df_pl_idx = df_c_pl.set_index("group") if not df_c_pl.empty else pd.DataFrame()
            pos_ct_c = {str(g): i for i, g in enumerate(df_c_ct["group"])} if not df_c_ct.empty else {}
            pos_all_c = {str(g): i for i, g in enumerate(df_all_c["group"])} if not df_all_c.empty else {}

            for grp in orig_all_roster:
                is_pool = grp.startswith("All_")
                orig_present = grp in o_p.index
                orig_n = int(o_p.loc[grp, "n_donors"]) if orig_present else int(orig_p1.loc[grp, "n_donors"])
                if tier == "P1m" and orig_present:
                    oref = o_p.loc[grp]
                    op_c = float(oref["composite_pvalue"])
                    oq_c = float(oref["q_bh_composite"]) if (not is_pool and pd.notna(oref["q_bh_composite"])) else op_c
                    op_e = float(oref["empirical_null_pvalue"])
                    oq_e = float(oref["q_bh_empirical"]) if (not is_pool and pd.notna(oref["q_bh_empirical"])) else op_e
                    oeff = float(oref["composite_effect"])
                    orig_cls = transport_class(oq_c, oq_e, True)
                    orig_sig_flag: bool | None = orig_cls == "BOTH_PASS"
                else:
                    op_c = oq_c = op_e = oq_e = oeff = np.nan
                    orig_cls = "NO_ORIG_P0" if tier == "P0m" else "NOT_IN_ORIG_SENS"
                    orig_sig_flag = None

                src_idx = df_pl_idx if is_pool else df_ct_idx
                if grp in sk_all or grp not in src_idx.index:
                    reason = sk_all.get(grp, "untestable_at_30plus")
                    n_mat = int((cell_counts_30.loc[list(dset), grp] >= MIN_CELLS_PER_DONOR).sum()) if grp in cell_counts_30.columns else 0
                    out_rows.append(
                        make_row(
                            cohort_name,
                            grp,
                            fam_key,
                            n_mat,
                            orig_n - n_mat,
                            op_c if np.isfinite(op_c) else None,
                            oq_c if np.isfinite(oq_c) else None,
                            None,
                            None,
                            None,
                            None,
                            "UNTESTABLE",
                            f"{reason};orig_class={orig_cls}",
                            arm=arm_name,
                            orig_transport_class=orig_cls,
                            mature_transport_class="UNTESTABLE",
                        )
                    )
                    continue

                r = src_idx.loc[grp]
                mat_qc = float(q_ac_e[pos_all_c[grp]])
                mat_qe = float(q_ae_e[pos_all_c[grp]])
                mat_cls = transport_class(mat_qc, mat_qe, True)
                verd, vnote = verdict_standard(
                    orig_sig_flag,
                    mat_cls == "BOTH_PASS",
                    oeff if np.isfinite(oeff) else None,
                    float(r["composite_effect"]),
                )
                mde_null = mde80_from_null_sd(float(r["null_sd_effect"]))
                mde_se = mde80_from_se(float(r["composite_se"]), int(r["residual_df"]))
                pa = pos_all_c[grp]
                notes = (
                    f"orig_class={orig_cls};mature_class={mat_cls};"
                    f"q_all_comp={float(q_ac_e[pa]):.4g};q_all_emp={float(q_ae_e[pa]):.4g};"
                    f"mde80_se={mde_se:.4f};draws={P1_ITERATIONS};seed={ORIG_NULL_SEED};"
                    f"min_decile_margin={int(r['min_decile_margin'])}"
                )
                if vnote:
                    notes += f";{vnote}"
                out_rows.append(
                    make_row(
                        cohort_name,
                        grp,
                        fam_key,
                        int(r["n_donors"]),
                        orig_n - int(r["n_donors"]),
                        op_c if np.isfinite(op_c) else None,
                        oq_c if np.isfinite(oq_c) else None,
                        float(r["composite_pvalue"]),
                        mat_qc,
                        float(r["empirical_null_pvalue"]),
                        mat_qe,
                        verd,
                        notes,
                        arm=arm_name,
                        residual_df=int(r["residual_df"]),
                        original_effect=oeff,
                        mature_effect=float(r["composite_effect"]),
                        mature_se=float(r["composite_se"]),
                        null_sd=float(r["null_sd_effect"]),
                        mde80=mde_null,
                        mde80_method="2.4865_times_null_sd",
                        original_empirical_p=op_e,
                        original_empirical_q=oq_e,
                        orig_transport_class=orig_cls,
                        mature_transport_class=mat_cls,
                        n_up_detected=int(r["n_up_detected"]),
                        n_down_detected=int(r["n_down_detected"]),
                        covariates=str(r["covariates"]),
                        omitted_covariates=str(r["omitted_covariates"]),
                    )
                )
                if not is_pool:
                    pi = pos_ct_c[grp]
                    bh_rows.append(
                        {
                            "cohort": cohort_name,
                            "family": f"guilliams_{arm_name}_{tier}_ct_composite",
                            "unit": grp,
                            "family_size": int(len(df_c_ct)),
                            "raw_p": float(r["composite_pvalue"]),
                            "q_bh_explicit": float(q_cc_e[pi]),
                            "q_bh_statsmodels": float(q_cc_s[pi]),
                            "abs_diff": float(abs(q_cc_e[pi] - q_cc_s[pi])),
                            "family_max_abs_diff": diff_cc,
                        }
                    )
                    bh_rows.append(
                        {
                            "cohort": cohort_name,
                            "family": f"guilliams_{arm_name}_{tier}_ct_empirical",
                            "unit": grp,
                            "family_size": int(len(df_c_ct)),
                            "raw_p": float(r["empirical_null_pvalue"]),
                            "q_bh_explicit": float(q_ce_e[pi]),
                            "q_bh_statsmodels": float(q_ce_s[pi]),
                            "abs_diff": float(abs(q_ce_e[pi] - q_ce_s[pi])),
                            "family_max_abs_diff": diff_ce,
                        }
                    )
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"guilliams_{arm_name}_{tier}_all_composite",
                        "unit": grp,
                        "family_size": int(len(df_all_c)),
                        "raw_p": float(r["composite_pvalue"]),
                        "q_bh_explicit": float(q_ac_e[pa]),
                        "q_bh_statsmodels": float(q_ac_s[pa]),
                        "abs_diff": float(abs(q_ac_e[pa] - q_ac_s[pa])),
                        "family_max_abs_diff": diff_ac,
                    }
                )
                bh_rows.append(
                    {
                        "cohort": cohort_name,
                        "family": f"guilliams_{arm_name}_{tier}_all_empirical",
                        "unit": grp,
                        "family_size": int(len(df_all_c)),
                        "raw_p": float(r["empirical_null_pvalue"]),
                        "q_bh_explicit": float(q_ae_e[pa]),
                        "q_bh_statsmodels": float(q_ae_s[pa]),
                        "abs_diff": float(abs(q_ae_e[pa] - q_ae_s[pa])),
                        "family_max_abs_diff": diff_ae,
                    }
                )

            bp = int(((q_ac_e < ALPHA) & (q_ae_e < ALPHA)).sum()) if not df_all_c.empty else 0
            eo = int(((q_ac_e >= ALPHA) & (q_ae_e < ALPHA)).sum()) if not df_all_c.empty else 0
            co = int(((q_ac_e < ALPHA) & (q_ae_e >= ALPHA)).sum()) if not df_all_c.empty else 0
            tier_res[tier] = {
                "n_ct_tested": int(len(df_c_ct)),
                "both_pass": bp,
                "emp_only": eo,
                "comp_only": co,
                "df_pl": df_c_pl.to_dict("records") if not df_c_pl.empty else [],
            }

        n_ct_tests_g = int(len(g_ct))
        n_ct_sig_g = int((q_ct_exp < ALPHA).sum()) if len(q_ct_exp) else 0
        n_pl_tests_g = int(len(g_pl))
        n_pl_sig_g = int((q_pl_exp < ALPHA).sum()) if len(q_pl_exp) else 0
        arm_diagnostics[arm_name] = {
            "n_donors": int(len(dset)),
            "elig_skip": elig_skip,
            "ct_tests": n_ct_tests_g,
            "ct_sig": n_ct_sig_g,
            "pl_tests": n_pl_tests_g,
            "pl_sig": n_pl_sig_g,
            "tier_res": tier_res,
        }

        if arm_name == "primary":
            summary_rows.append(
                {
                    "cohort": cohort_name,
                    "endpoint_family": "fractions_clr",
                    "arm": "primary",
                    "donors_original": n_donors_orig,
                    "donors_mature": n_donors_30,
                    "donors_excluded": n_excluded_cohort,
                    "original_headline": "Not run (FACS-sorted CD45+/CD45- exclusion)",
                    "mature_headline": "UNTESTABLE (FACS-sorted CD45+/CD45- exclusion retained verbatim)",
                    "mde80_summary": "NA (sorted cells)",
                }
            )
            summary_rows.append(
                {
                    "cohort": cohort_name,
                    "endpoint_family": "concordance",
                    "arm": "primary",
                    "donors_original": n_donors_orig,
                    "donors_mature": n_donors_30,
                    "donors_excluded": n_excluded_cohort,
                    "original_headline": "Not run (FACS-sorted exclusion)",
                    "mature_headline": "UNTESTABLE (FACS-sorted exclusion retained verbatim)",
                    "mde80_summary": "NA (sorted cells)",
                }
            )
            summary_rows.append(
                {
                    "cohort": cohort_name,
                    "endpoint_family": "genes",
                    "arm": "primary",
                    "donors_original": n_donors_orig,
                    "donors_mature": n_donors_30,
                    "donors_excluded": n_excluded_cohort,
                    "original_headline": "0/81881 ct hits (7 types, min q=0.6100); 0/26220 pooled",
                    "mature_headline": f"{n_ct_sig_g}/{n_ct_tests_g} ct hits ({len(ct_e)}/7 types tested, {7 - len(ct_e)} untestable; min q={float(q_ct_exp.min()) if len(q_ct_exp) else np.nan:.4f}); {n_pl_sig_g}/{n_pl_tests_g} pooled",
                    "mde80_summary": f"median cell-type median MDE80={float(s_ct[s_ct['status']=='tested']['median_mde80'].median()):.4f} log2CPM/dec (t_quantile_se)",
                }
            )
            for tier in ("P1m", "P0m"):
                tr = tier_res[tier]
                pl_desc = "; ".join(
                    f"{p['group']} eff={float(p['composite_effect']):+.4f}(p_comp={float(p['composite_pvalue']):.4g},p_emp={float(p['empirical_null_pvalue']):.4g})"
                    for p in tr["df_pl"]
                )
                sub_mde = [
                    r["mde80"]
                    for r in out_rows
                    if r["endpoint_family"] == f"transport_{tier}" and r["verdict"] != "UNTESTABLE"
                ]
                summary_rows.append(
                    {
                        "cohort": cohort_name,
                        "endpoint_family": f"transport_{tier}",
                        "arm": "primary",
                        "donors_original": n_donors_orig,
                        "donors_mature": n_donors_30,
                        "donors_excluded": n_excluded_cohort,
                        "original_headline": (
                            "0/7 ct comp or emp passes (min q_emp=0.1259); pooled All_CD45pos emp_p=0.0085, All_Cells emp_p=0.0165 negative"
                            if tier == "P1m"
                            else "No original P0 run"
                        ),
                        "mature_headline": f"{tr['both_pass']}/{tr['n_ct_tested'] + len(tr['df_pl'])} both-pass, {tr['emp_only']} emp-only, {tr['comp_only']} comp-only ({7 - tr['n_ct_tested']} untestable); {pl_desc}",
                        "mde80_summary": f"median={float(np.nanmedian(sub_mde)):.4f} (range {float(np.nanmin(sub_mde)):.4f} to {float(np.nanmax(sub_mde)):.4f}) per dec (2.4865 * null_sd)",
                    }
                )
        append_runtime(cohort_name, f"arm_{arm_name}", s0_arm, t0_arm, "ok", f"n={len(dset)}, ct_sig={n_ct_sig_g}")

    gui_df = order_cohort_df(out_rows)
    write_clean_csv(gui_df, RESULTS_DIR / f"{PREFIX}-guilliams.csv")

    save_cohort_payload(
        "guilliams",
        {
            "summary_rows": summary_rows,
            "bh_rows": bh_rows,
            "diagnostics": {
                "n_donors_orig": n_donors_orig,
                "n_donors_30": n_donors_30,
                "n_excluded": n_excluded_cohort,
                "age_min": float(donor_meta_30["age_years"].min()),
                "age_max": float(donor_meta_30["age_years"].max()),
                "n_male": int(donor_meta_30["sex_male"].sum()),
                "n_female": int((1 - donor_meta_30["sex_male"]).sum()),
                "n_citeseq": int(donor_meta_30["is_citeseq"].sum()),
                "full_object_donors_total": full_obj_donors_total,
                "full_object_under30_donors": full_obj_under30,
                "arms": arm_diagnostics,
            },
        },
    )
    append_runtime(cohort_name, "cohort_total", s0_all, t0_all, "ok", f"rows={len(gui_df)}")
    print(f"[Guilliams DONE] rows={len(gui_df)}, elapsed={time.monotonic() - t0_all:.1f}s", flush=True)


# =============================================================================
# FINALIZE: AGGREGATE SUMMARY, BH CHECK, SHA-256 INTEGRITY & ASSERTIONS
# =============================================================================

def run_finalize() -> None:
    t0 = time.monotonic()
    s0 = utc_now()
    print("\n================ FINALIZING BATCH A OUTPUTS ================", flush=True)

    summary_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    diag_all: dict[str, Any] = {}

    for key in ("gao", "yazar", "perez", "guilliams"):
        payload = load_cohort_payload(key)
        summary_rows.extend(payload["summary_rows"])
        bh_rows.extend(payload["bh_rows"])
        diag_all[key] = payload["diagnostics"]

    summary_df = pd.DataFrame(summary_rows)
    write_clean_csv(summary_df, SUMMARY_CSV)

    bh_df = pd.DataFrame(bh_rows).drop_duplicates(subset=["cohort", "family", "unit"], keep="first").reset_index(drop=True)
    max_diff = float(bh_df["abs_diff"].max()) if not bh_df.empty else 0.0
    max_fam_diff = float(bh_df["family_max_abs_diff"].max()) if not bh_df.empty else 0.0
    assert max_diff <= 1e-12 and max_fam_diff <= 1e-12, f"BH check failed: max_diff={max_diff}, max_fam_diff={max_fam_diff}"
    write_clean_csv(bh_df, BH_CHECK_CSV)

    # Verify SHA-256 of EVERY original result file read has NOT changed
    before_path = TMP_WORK_DIR / "original_sha256_before.json"
    before_map = json.loads(before_path.read_text(encoding="utf-8"))
    after_map: dict[str, str] = {}
    for rel, h_before in before_map.items():
        h_after = sha256_file(REPO_ROOT / rel)
        after_map[rel] = h_after
        if h_after != h_before:
            raise AssertionError(f"Original file modified! {rel}: before={h_before}, after={h_after}")
    (TMP_WORK_DIR / "original_sha256_verified.json").write_text(
        json.dumps({"verified_utc": utc_now(), "all_unchanged": True, "files": after_map}, indent=2) + "\n",
        encoding="utf-8",
    )

    # Verify all 4 cohort CSVs schema + zero forbidden dashes + cross-check q values against bh-check
    for key in ("gao", "yazar", "perez", "guilliams"):
        p = RESULTS_DIR / f"{PREFIX}-{key}.csv"
        df = pd.read_csv(p)
        assert list(df.columns[: len(REQUIRED_COHORT_COLS)]) == REQUIRED_COHORT_COLS, f"Schema mismatch in {p}"
        assert_no_forbidden_dashes(p.read_text(encoding="utf-8"), label=str(p))

    (TMP_WORK_DIR / "diagnostics_all.json").write_text(
        json.dumps(diag_all, indent=2, default=str) + "\n", encoding="utf-8"
    )

    append_runtime(
        "Batch_A_All",
        "finalize_and_verify",
        s0,
        t0,
        "ok",
        f"summary_rows={len(summary_df)}, bh_rows={len(bh_df)}, max_bh_diff={max_diff:.3e}, sha256_checked={len(after_map)}",
    )
    print(
        f"[FINALIZE OK] summary={len(summary_df)} rows, bh_check={len(bh_df)} rows "
        f"(max_diff={max_diff:.3e}, max_fam_diff={max_fam_diff:.3e}), "
        f"sha256_unchanged={len(after_map)}/{len(after_map)}",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cohort",
        required=True,
        choices=["gao", "yazar", "perez", "guilliams", "finalize", "all"],
    )
    args = parser.parse_args()

    record_sha256_manifest()
    if args.cohort == "finalize":
        run_finalize()
        return 0

    programs = load_verified_mature_programs()
    if args.cohort == "gao":
        run_cohort_gao(programs)
    elif args.cohort == "yazar":
        run_cohort_yazar(programs)
    elif args.cohort == "perez":
        run_cohort_perez(programs)
    elif args.cohort == "guilliams":
        run_cohort_guilliams(programs)
    elif args.cohort == "all":
        run_cohort_gao(programs)
        gc.collect()
        run_cohort_yazar(programs)
        gc.collect()
        run_cohort_perez(programs)
        gc.collect()
        run_cohort_guilliams(programs)
        gc.collect()
        run_finalize()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
