"""Shared statistical, resource-guard, and RDS-extraction helpers for Datasets 01 to 20."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd
import rdata
from rdata.parser import RObjectType
from scipy import sparse as sp
from scipy import stats as scipy_stats

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import resource_guard
from src.aging_programs import load_program_tiers
from src.followup_muscle import (
    NULL_ITERATIONS,
    NULL_SEED,
    add_bh,
    clr_transform,
    count_gate,
    fit_ols,
    matched_null_effect_pvalue,
    random_signature_scores,
    technical_gene_mask,
)
from src.reviewer_stats import pairwise_pearson

DATA_TMP_ROOT = Path("/home/reza/projects/agescape-revision-data/tmp_organs")
MIN_FREE_GB = 30.0
FORBIDDEN_DASH_RE = re.compile(r"[\u2013\u2014]")
ERCC_RE = re.compile(r"^ERCC[-_]", re.IGNORECASE)


def assert_no_forbidden_dashes(text: str, label: str = "output") -> None:
    """Reject em dashes and en dashes in generated text and tables."""
    match = FORBIDDEN_DASH_RE.search(text)
    if match:
        raise ValueError(f"Forbidden Unicode dash {match.group(0)!r} found in {label}")


def write_clean_markdown(text: str, path: Path) -> Path:
    """Validate text for dash rules and write markdown atomically."""
    assert_no_forbidden_dashes(text, label=str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def write_clean_csv(frame: pd.DataFrame, path: Path) -> Path:
    """Validate CSV serialization for dash rules and write to disk."""
    csv_text = frame.to_csv(index=False)
    assert_no_forbidden_dashes(csv_text, label=str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(csv_text, encoding="utf-8")
    return path


def sanitize_dashes(text: str) -> str:
    """Replace external Unicode en/em dashes from metadata with ASCII hyphen-minus or 'to'."""
    return text.replace("\u2014", " - ").replace("\u2013", "-")


def download_with_guard(url: str, dest: Path, max_file_gb: float = 20.0) -> Path:
    """Download a file using scripts/resource_guard.py with fallback when HEAD lacks size."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    resource_guard.check_free_space(dest.parent, MIN_FREE_GB)
    args = argparse.Namespace(
        url=url,
        output=dest,
        min_free_gb=MIN_FREE_GB,
        max_file_gb=float(max_file_gb),
    )
    try:
        rc = resource_guard.download_guarded(args)
        if rc != 0:
            raise RuntimeError(f"curl exited with status {rc} downloading {url}")
        return dest
    except resource_guard.ResourceGuardError as exc:
        msg = str(exc)
        header_missing = "cannot retrieve download headers" in msg or "missing or invalid Content-Length" in msg
        if not header_missing:
            raise
    with resource_guard.heavy_job_lock():
        resource_guard.check_free_space(dest.parent, MIN_FREE_GB)
        command = resource_guard.build_curl_command(
            url, dest, resource_guard.gb_to_bytes(float(max_file_gb))
        )
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"curl fallback exited with status {completed.returncode} downloading {url}")
    resource_guard.check_free_space(dest.parent, MIN_FREE_GB)
    return dest


def load_p1_clean_human(repo_root: Path = REPO_ROOT) -> dict[str, list[str]]:
    """Return canonical mapped human P1_Clean UP (76) and DOWN (42) gene symbols."""
    return load_program_tiers(repo_root, "human")["P1_Clean"]


def extended_technical_mask(gene_names: np.ndarray | pd.Index) -> np.ndarray:
    """Mask MT-/RP[LS]/MRP[LS]/sex-chromosome/ERCC genes before library-size normalization."""
    arr = np.asarray(gene_names, dtype=str)
    base_mask = technical_gene_mask(arr)
    ercc_mask = np.array([bool(ERCC_RE.match(g.strip())) for g in arr], dtype=bool)
    return base_mask | ercc_mask


def prepare_counts_matrix(
    counts_df: pd.DataFrame,
    allow_fractional_expected: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate raw counts, collapse duplicate symbols by sum, and strip technical genes."""
    clean_idx = counts_df.index.astype(str).str.strip()
    valid_name = (clean_idx != "") & (clean_idx.str.lower() != "nan") & (clean_idx.str.lower() != "na")
    work = counts_df.loc[valid_name].copy()
    work.index = clean_idx[valid_name]
    if work.index.has_duplicates:
        work = work.groupby(work.index, sort=False).sum()

    vals = work.to_numpy(dtype=np.float64)
    gate = count_gate(vals)
    if gate["nonfinite"] > 0 or gate["negative"] > 0:
        raise ValueError(f"Non-finite or negative values in raw count matrix: {gate}")
    if not allow_fractional_expected and gate["noninteger"] > 0:
        raise ValueError(f"Non-integer values in integer count matrix: {gate}")

    tech_mask = extended_technical_mask(work.index)
    n_tech = int(tech_mask.sum())
    filtered = work.loc[~tech_mask].copy()
    lib_sizes = filtered.sum(axis=0)
    stats = {
        "n_genes_input": int(counts_df.shape[0]),
        "n_genes_unique": int(work.shape[0]),
        "n_technical_removed": n_tech,
        "n_genes_nontech": int(filtered.shape[0]),
        "count_gate": gate,
        "min_lib_size": float(lib_sizes.min()),
        "median_lib_size": float(lib_sizes.median()),
        "max_lib_size": float(lib_sizes.max()),
    }
    return filtered, stats


def counts_to_log2cpm(counts_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Compute log2(CPM + 1) and library sizes from non-technical raw counts."""
    lib_sizes = counts_df.sum(axis=0).astype(float)
    safe_lib = np.maximum(lib_sizes.to_numpy(dtype=float), 1.0)
    cpm = counts_df.to_numpy(dtype=np.float64) / safe_lib[None, :] * 1e6
    log2cpm = np.log2(cpm + 1.0)
    return pd.DataFrame(log2cpm, index=counts_df.index, columns=counts_df.columns), lib_sizes


def select_design_columns(
    donor_meta: pd.DataFrame,
    age_col: str = "age_per_decade",
    candidate_covariates: tuple[str, ...] = ("sex_male", "log_lib"),
    min_residual_df: int = 2,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Build full-rank OLS design [intercept, age_col, non-collinear covariates...]."""
    if age_col not in donor_meta:
        raise ValueError(f"Missing age column {age_col!r}")
    age_vals = donor_meta[age_col].to_numpy(dtype=float)
    if not np.isfinite(age_vals).all() or np.unique(age_vals).size < 2:
        raise ValueError(f"Age predictor {age_col!r} must be finite with >= 2 distinct values")

    n = len(donor_meta)
    columns = [np.ones(n, dtype=float), age_vals]
    names = ["intercept", age_col]
    omitted: list[str] = []

    for cov in candidate_covariates:
        if (n - (len(columns) + 1)) < min_residual_df:
            omitted.append(f"{cov}(df_limit)")
            continue
        if cov not in donor_meta:
            omitted.append(f"{cov}(missing)")
            continue
        cov_vals = donor_meta[cov].to_numpy(dtype=float)
        if not np.isfinite(cov_vals).all() or np.unique(cov_vals).size < 2:
            omitted.append(f"{cov}(constant_or_na)")
            continue
        candidate = np.column_stack([*columns, cov_vals])
        if np.linalg.matrix_rank(candidate) == candidate.shape[1]:
            columns.append(cov_vals)
            names.append(cov)
        else:
            omitted.append(f"{cov}(collinear)")

    design = np.column_stack(columns)
    if (n - design.shape[1]) < min_residual_df:
        raise ValueError(
            f"Insufficient donors (n={n}, params={design.shape[1]}) for min_residual_df={min_residual_df}"
        )
    return design, names, omitted


def run_family_a_clr(
    cell_counts_df: pd.DataFrame,
    donor_meta: pd.DataFrame,
    age_col: str = "age_per_decade",
    candidate_covariates: tuple[str, ...] = ("sex_male",),
    min_donors_present: int = 3,
    min_total_cells: int = 20,
    min_residual_df: int = 2,
) -> pd.DataFrame:
    """Family A: donor-level CLR fraction OLS with BH correction across K tested cell types."""
    common_donors = [d for d in donor_meta.index if d in cell_counts_df.index]
    counts = cell_counts_df.loc[common_donors].copy()
    meta = donor_meta.loc[common_donors].copy()
    totals = counts.sum(axis=1).astype(float)
    props = counts.div(np.maximum(totals, 1.0), axis=0)
    clr = clr_transform(counts, pseudocount=0.5)

    present_donors = (counts > 0).sum(axis=0)
    total_ct_cells = counts.sum(axis=0)
    tested_cts = [
        ct
        for ct in counts.columns
        if int(present_donors[ct]) >= min_donors_present and int(total_ct_cells[ct]) >= min_total_cells
    ]
    if not tested_cts:
        return pd.DataFrame()

    design, names, omitted = select_design_columns(
        meta,
        age_col=age_col,
        candidate_covariates=candidate_covariates,
        min_residual_df=min_residual_df,
    )
    outcomes = clr[tested_cts].to_numpy(dtype=float)
    ols = fit_ols(outcomes, design, effect_column=1)

    rows = []
    for idx, ct in enumerate(tested_cts):
        ct_props = props[ct].to_numpy(dtype=float)
        rows.append(
            {
                "cell_type": str(ct),
                "n_donors": int(len(common_donors)),
                "n_donors_present": int(present_donors[ct]),
                "total_cells": int(total_ct_cells[ct]),
                "median_fraction": float(np.median(ct_props)),
                "min_fraction": float(np.min(ct_props)),
                "max_fraction": float(np.max(ct_props)),
                "predictor": age_col,
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "family_a_size": int(len(tested_cts)),
                "effect_clr": float(ols["effect"][idx]),
                "se": float(ols["se"][idx]),
                "ci95_low": float(ols["ci95_low"][idx]),
                "ci95_high": float(ols["ci95_high"][idx]),
                "t_stat": float(ols["effect"][idx] / ols["se"][idx]) if ols["se"][idx] > 0 else np.nan,
                "pvalue": float(ols["pvalue"][idx]),
                "residual_df": int(ols["residual_df"][idx]),
                "r2": float(ols["r2"][idx]),
            }
        )
    res = pd.DataFrame(rows)
    res = add_bh(res, ["family_a_size"], p_column="pvalue", q_column="q_bh_family_a")
    res = res.sort_values("pvalue", ascending=True).reset_index(drop=True)
    return res


def run_family_b_genes_and_concordance(
    expr_by_group: dict[str, pd.DataFrame],
    donor_meta_by_group: dict[str, pd.DataFrame],
    age_col: str = "age_per_decade",
    candidate_covariates: tuple[str, ...] = ("sex_male", "log_lib"),
    min_donors: int = 4,
    min_cpm: float = 1.0,
    min_expr_frac: float = 0.5,
    min_residual_df: int = 2,
    is_unbiased_single_cell: bool = True,
    concordance_perms: int = 1000,
    seed: int = NULL_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Family B: donor-level gene age OLS (with global BH) plus unbiased cross-cell-type concordance."""
    log2_thresh = float(np.log2(min_cpm + 1.0))
    gene_frames: list[pd.DataFrame] = []
    group_summaries: list[dict[str, Any]] = []
    t_series_by_group: dict[str, pd.Series] = {}
    designs_by_group: dict[str, tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]] = {}

    for group, log2cpm_df in expr_by_group.items():
        meta = donor_meta_by_group[group]
        common_donors = [d for d in meta.index if d in log2cpm_df.columns]
        if len(common_donors) < min_donors:
            group_summaries.append(
                {
                    "group": group,
                    "status": "skipped_low_donors",
                    "n_donors": int(len(common_donors)),
                    "n_genes_tested": 0,
                }
            )
            continue
        mat_df = log2cpm_df[common_donors]
        m_meta = meta.loc[common_donors]
        min_req = max(3, int(np.ceil(min_expr_frac * len(common_donors))))
        keep_mask = (mat_df.to_numpy(dtype=float) >= log2_thresh).sum(axis=1) >= min_req
        if int(keep_mask.sum()) == 0:
            group_summaries.append(
                {
                    "group": group,
                    "status": "skipped_no_expressed_genes",
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
                    "status": f"skipped_{exc}",
                    "n_donors": int(len(common_donors)),
                    "n_genes_tested": 0,
                }
            )
            continue

        outcomes = mat_sub.to_numpy(dtype=float).T
        ols = fit_ols(outcomes, design, effect_column=1)
        t_stats = np.where(ols["se"] > 0, ols["effect"] / ols["se"], np.nan)
        t_series_by_group[group] = pd.Series(t_stats, index=mat_sub.index)
        designs_by_group[group] = (mat_sub, m_meta, design, names)

        df_g = pd.DataFrame(
            {
                "group": str(group),
                "gene": mat_sub.index.astype(str),
                "n_donors": int(len(common_donors)),
                "predictor": age_col,
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "mean_log2cpm": outcomes.mean(axis=0),
                "effect": ols["effect"],
                "se": ols["se"],
                "ci95_low": ols["ci95_low"],
                "ci95_high": ols["ci95_high"],
                "t_stat": t_stats,
                "pvalue": ols["pvalue"],
                "residual_df": ols["residual_df"].astype(int),
                "r2": ols["r2"],
            }
        )
        df_g = add_bh(df_g, ["group"], p_column="pvalue", q_column="q_bh_within_group")
        gene_frames.append(df_g)
        group_summaries.append(
            {
                "group": str(group),
                "status": "tested",
                "n_donors": int(len(common_donors)),
                "n_genes_tested": int(len(df_g)),
                "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
                "omitted_covariates": ";".join(omitted) if omitted else "none",
                "residual_df": int(ols["residual_df"][0]),
            }
        )

    if not gene_frames:
        return pd.DataFrame(), pd.DataFrame(group_summaries), pd.DataFrame()

    all_genes = pd.concat(gene_frames, ignore_index=True)
    all_genes["family_b_global_size"] = int(len(all_genes))
    all_genes = add_bh(
        all_genes, ["family_b_global_size"], p_column="pvalue", q_column="q_bh_global"
    )

    summary_df = pd.DataFrame(group_summaries)
    within_hits = (
        all_genes.groupby("group")["q_bh_within_group"].apply(lambda s: int((s < 0.05).sum())).to_dict()
    )
    global_hits = (
        all_genes.groupby("group")["q_bh_global"].apply(lambda s: int((s < 0.05).sum())).to_dict()
    )
    min_p_map = all_genes.groupby("group")["pvalue"].min().to_dict()
    min_q_map = all_genes.groupby("group")["q_bh_global"].min().to_dict()
    summary_df["n_sig_within_q05"] = summary_df["group"].map(within_hits).fillna(0).astype(int)
    summary_df["n_sig_global_q05"] = summary_df["group"].map(global_hits).fillna(0).astype(int)
    summary_df["min_pvalue"] = summary_df["group"].map(min_p_map)
    summary_df["min_q_bh_global"] = summary_df["group"].map(min_q_map)
    summary_df["family_b_global_size"] = int(len(all_genes))

    concordance_df = pd.DataFrame()
    if is_unbiased_single_cell and len(t_series_by_group) >= 2:
        concordance_df = _compute_concordance_table(
            t_series_by_group=t_series_by_group,
            designs_by_group=designs_by_group,
            age_col=age_col,
            concordance_perms=concordance_perms,
            seed=seed,
        )

    all_genes = all_genes.sort_values("pvalue", ascending=True).reset_index(drop=True)
    return all_genes, summary_df, concordance_df


def _compute_concordance_table(
    t_series_by_group: dict[str, pd.Series],
    designs_by_group: dict[str, tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]],
    age_col: str,
    concordance_perms: int,
    seed: int,
) -> pd.DataFrame:
    """Compute pairwise Pearson concordance across cell types plus donor-age permutation p."""
    groups = sorted(t_series_by_group.keys())
    common_genes = sorted(set.intersection(*(set(t_series_by_group[g].index) for g in groups)))
    if len(common_genes) < 50:
        return pd.DataFrame(
            [
                {
                    "n_cell_types": len(groups),
                    "n_common_genes": len(common_genes),
                    "status": "skipped_few_common_genes",
                }
            ]
        )

    t_mat = np.vstack([t_series_by_group[g].loc[common_genes].to_numpy(dtype=float) for g in groups])
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

    rng = np.random.default_rng(seed)
    all_donors = sorted(set.union(*(set(designs_by_group[g][1].index) for g in groups)))
    donor_age = {}
    for g in groups:
        _, meta_g, _, _ = designs_by_group[g]
        for d, val in meta_g[age_col].items():
            donor_age[d] = float(val)
    age_vec = np.array([donor_age[d] for d in all_donors], dtype=float)

    common_arr_by_group = {}
    donor_indices_by_group = {}
    donor_to_pos = {d: idx for idx, d in enumerate(all_donors)}
    for g in groups:
        mat_sub, meta_g, design_g, _ = designs_by_group[g]
        common_arr_by_group[g] = (mat_sub.loc[common_genes].to_numpy(dtype=float).T, design_g)
        donor_indices_by_group[g] = np.array([donor_to_pos[d] for d in meta_g.index], dtype=int)

    perm_means = []
    for _ in range(concordance_perms):
        perm_ages = rng.permutation(age_vec)
        perm_t_rows = []
        ok = True
        for g in groups:
            expr_mat, design_g = common_arr_by_group[g]
            d_g = design_g.copy()
            d_g[:, 1] = perm_ages[donor_indices_by_group[g]]
            if np.linalg.matrix_rank(d_g) < d_g.shape[1]:
                ok = False
                break
            ols_p = fit_ols(expr_mat, d_g, effect_column=1)
            perm_t_rows.append(np.nan_to_num(ols_p["effect"] / ols_p["se"], nan=0.0))
        if ok:
            perm_r = pairwise_pearson(np.vstack(perm_t_rows))
            perm_means.append(float(np.mean(perm_r)))

    perm_arr = np.asarray(perm_means, dtype=float)
    perm_p_one_sided = (
        float((1 + np.sum(perm_arr >= obs_mean_r)) / (len(perm_arr) + 1))
        if len(perm_arr) > 0
        else np.nan
    )
    perm_p_two_sided = (
        float((1 + np.sum(np.abs(perm_arr) >= abs(obs_mean_r))) / (len(perm_arr) + 1))
        if len(perm_arr) > 0
        else np.nan
    )

    return pd.DataFrame(
        [
            {
                "n_cell_types": int(len(groups)),
                "n_pairs": int(len(pair_r)),
                "n_common_genes": int(len(common_genes)),
                "mean_pairwise_pearson_r": obs_mean_r,
                "median_pairwise_pearson_r": obs_median_r,
                "min_pairwise_pearson_r": float(np.min(pair_r)),
                "max_pairwise_pearson_r": float(np.max(pair_r)),
                "mean_pairwise_spearman_rho": float(np.mean(pair_spear)),
                "mean_sign_agreement": float(np.mean(pair_sign)),
                "n_valid_perms": int(len(perm_arr)),
                "perm_p_one_sided": perm_p_one_sided,
                "perm_p_two_sided": perm_p_two_sided,
                "cell_types": ";".join(groups),
            }
        ]
    )


def run_family_c_p1_signature(
    expr_by_group: dict[str, pd.DataFrame],
    donor_meta_by_group: dict[str, pd.DataFrame],
    p1_program: dict[str, list[str]] | None = None,
    age_col: str = "age_per_decade",
    candidate_covariates: tuple[str, ...] = ("sex_male", "log_lib"),
    min_donors: int = 4,
    min_cpm: float = 1.0,
    min_expr_frac: float = 0.5,
    min_residual_df: int = 2,
    iterations: int = NULL_ITERATIONS,
    seed: int = NULL_SEED,
) -> pd.DataFrame:
    """Family C: P1_Clean mouse aging signature OLS + 200-draw matched random gene-set null + BH."""
    if p1_program is None:
        p1_program = load_p1_clean_human(REPO_ROOT)
    up_set = set(p1_program["up"])
    down_set = set(p1_program["down"])
    log2_thresh = float(np.log2(min_cpm + 1.0))

    rows = []
    for group, log2cpm_df in expr_by_group.items():
        meta = donor_meta_by_group[group]
        common_donors = [d for d in meta.index if d in log2cpm_df.columns]
        if len(common_donors) < min_donors:
            continue
        mat_df = log2cpm_df[common_donors]
        m_meta = meta.loc[common_donors]
        min_req = max(3, int(np.ceil(min_expr_frac * len(common_donors))))
        keep_mask = (mat_df.to_numpy(dtype=float) >= log2_thresh).sum(axis=1) >= min_req
        mat_sub = mat_df.loc[keep_mask]
        if mat_sub.shape[0] < 200:
            continue

        genes = mat_sub.index.astype(str).tolist()
        up_idx = np.array([i for i, g in enumerate(genes) if g in up_set], dtype=int)
        down_idx = np.array([i for i, g in enumerate(genes) if g in down_set], dtype=int)
        if len(up_idx) == 0 or len(down_idx) == 0:
            continue

        try:
            design, names, omitted = select_design_columns(
                m_meta,
                age_col=age_col,
                candidate_covariates=candidate_covariates,
                min_residual_df=min_residual_df,
            )
        except ValueError:
            continue

        mat = mat_sub.to_numpy(dtype=float)
        up_score = mat[up_idx, :].mean(axis=0)
        down_score = mat[down_idx, :].mean(axis=0)
        comp_score = up_score - down_score
        ols = fit_ols(np.column_stack([up_score, down_score, comp_score]), design, effect_column=1)

        null_scores = random_signature_scores(
            mat, up_idx, down_idx, iterations=iterations, seed=seed
        )
        null_ols = fit_ols(null_scores.T, design, effect_column=1)
        null_effects = null_ols["effect"]
        p_emp = matched_null_effect_pvalue(float(ols["effect"][2]), null_effects)

        rows.append(
            {
                "group": str(group),
                "n_donors": int(len(common_donors)),
                "n_expressed_pool": int(len(genes)),
                "n_p1_up_detected": int(len(up_idx)),
                "n_p1_down_detected": int(len(down_idx)),
                "n_p1_total_detected": int(len(up_idx) + len(down_idx)),
                "predictor": age_col,
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
        return pd.DataFrame()
    res = pd.DataFrame(rows)
    res["family_c_size"] = int(len(res))
    res = add_bh(res, ["family_c_size"], p_column="composite_pvalue", q_column="q_bh_composite")
    res = add_bh(
        res, ["family_c_size"], p_column="empirical_null_pvalue", q_column="q_bh_empirical"
    )
    return res.sort_values("empirical_null_pvalue", ascending=True).reset_index(drop=True)


def _robject_to_python(obj: Any, converter: rdata.conversion.SimpleConverter) -> Any:
    """Convert an AST RObject using rdata SimpleConverter with dereferenced refs."""
    if getattr(obj.info, "type", None) == RObjectType.REF and obj.referenced_object is not None:
        return _robject_to_python(obj.referenced_object, converter)
    return converter.convert(obj)


def _robject_attrs_dict(obj: Any, converter: rdata.conversion.SimpleConverter) -> dict[str, Any]:
    """Return raw AST attribute mapping without forcing deep conversion of every slot."""
    if getattr(obj.info, "type", None) == RObjectType.REF and obj.referenced_object is not None:
        obj = obj.referenced_object
    raw_attrs = getattr(obj, "attributes", None)
    out: dict[str, Any] = {}
    while raw_attrs is not None and getattr(raw_attrs.info, "type", None) == RObjectType.LIST:
        tag = getattr(raw_attrs, "tag", None)
        if tag is not None:
            key = _robject_to_python(tag, converter)
            if isinstance(key, str):
                out[key] = raw_attrs.value[0]
        raw_attrs = raw_attrs.value[1]
    return out


def _dgc_from_robject(obj: Any, converter: rdata.conversion.SimpleConverter) -> sp.csc_matrix:
    """Convert an R dgCMatrix AST node into a SciPy csc_matrix without loading sibling assays."""
    attrs = _robject_attrs_dict(obj, converter)
    i_arr = np.asarray(_robject_to_python(attrs["i"], converter), dtype=np.int32)
    p_arr = np.asarray(_robject_to_python(attrs["p"], converter), dtype=np.int32)
    x_arr = np.asarray(_robject_to_python(attrs["x"], converter), dtype=np.float64)
    dim = tuple(int(v) for v in _robject_to_python(attrs["Dim"], converter))
    return sp.csc_matrix((x_arr, i_arr, p_arr), shape=dim)


def _dimnames_from_robject(
    obj: Any, converter: rdata.conversion.SimpleConverter
) -> tuple[np.ndarray, np.ndarray]:
    """Extract row and column Dimnames from an R dgCMatrix AST node."""
    attrs = _robject_attrs_dict(obj, converter)
    dimnames = _robject_to_python(attrs["Dimnames"], converter)
    rows = np.asarray(dimnames[0], dtype=str)
    cols = np.asarray(dimnames[1], dtype=str)
    return rows, cols


def extract_seurat_rds_selective(rds_path: Path) -> dict[str, Any]:
    """Selectively parse Seurat v2/v3/v4/v5 RDS meta.data, ident, and raw counts using rdata."""
    parsed = rdata.parser.parse_file(str(rds_path))
    conv = rdata.conversion.SimpleConverter(default_encoding="utf-8")
    root = parsed.object
    root_attrs = _robject_attrs_dict(root, conv)

    result: dict[str, Any] = {"classes": _robject_to_python(root_attrs.get("class"), conv) if "class" in root_attrs else None}
    if "meta.data" in root_attrs:
        result["meta.data"] = _robject_to_python(root_attrs["meta.data"], conv)

    for ident_key in ("active.ident", "ident"):
        if ident_key in root_attrs:
            result[ident_key] = _robject_to_python(root_attrs[ident_key], conv)

    # Seurat v2 stores @raw.data directly; v3/v4 stores @assays$RNA@counts
    if "raw.data" in root_attrs:
        raw_node = root_attrs["raw.data"]
        result["raw_counts"] = _dgc_from_robject(raw_node, conv)
        result["genes"], result["barcodes"] = _dimnames_from_robject(raw_node, conv)
    elif "assays" in root_attrs:
        assays_node = root_attrs["assays"]
        assays_names = _robject_to_python(_robject_attrs_dict(assays_node, conv)["names"], conv)
        assay_idx = list(assays_names).index("RNA") if "RNA" in list(assays_names) else 0
        rna_node = assays_node.value[assay_idx]
        rna_attrs = _robject_attrs_dict(rna_node, conv)
        if "counts" in rna_attrs:
            counts_node = rna_attrs["counts"]
            result["raw_counts"] = _dgc_from_robject(counts_node, conv)
            result["genes"], result["barcodes"] = _dimnames_from_robject(counts_node, conv)
        elif "layers" in rna_attrs:
            layers_node = rna_attrs["layers"]
            layer_names = list(_robject_to_python(_robject_attrs_dict(layers_node, conv)["names"], conv))
            count_layers = [i for i, ln in enumerate(layer_names) if ln == "counts" or ln.startswith("counts.")]
            if count_layers:
                counts_node = layers_node.value[count_layers[0]]
                result["raw_counts"] = _dgc_from_robject(counts_node, conv)
                result["genes"], result["barcodes"] = _dimnames_from_robject(counts_node, conv)
    return result


def clean_temp_dir(path: Path) -> None:
    """Remove temporary downloaded raw files after donor-level pseudobulk extraction."""
    if path.exists() and DATA_TMP_ROOT in path.resolve().parents:
        shutil.rmtree(path)
