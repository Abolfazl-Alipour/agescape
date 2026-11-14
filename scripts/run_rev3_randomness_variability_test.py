#!/usr/bin/env python3
"""Item 8: Cross-sectional 18m vs 30m mouse TMS randomness and variability test.

Fundamental limitation: Tabula Muris Senis is cross-sectional (four 18m mice:
2 females, 2 males; four 30m mice: 4 males only). Within-group variance at 18m
combines inter-individual variance with residual sexual dimorphism, whereas 30m
is all-male and reflects survivorship to 30 months.
"""

from __future__ import annotations

import gc
import itertools
from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.special import digamma
from scipy.stats import f as f_dist
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
H5AD_CANDIDATES = [
    REPO / "data" / "interim" / "m1-qc-cleaned.h5ad",
    Path("/home/reza/projects/agescape/data/interim/m1-qc-cleaned.h5ad"),
    Path("/home/reza/projects/agescape-revision-data/m1-qc-cleaned.h5ad"),
]
SEX_GENES = {"Xist", "Tsix", "Ddx3y", "Eif2s3y", "Uty", "Kdm5d"}
RAREFY_UMI = 1500
SEED = sum(map(ord, "agescape-rev3-randomness-variability-v1"))


def _resolve_h5ad() -> Path:
    for cand in H5AD_CANDIDATES:
        if cand.exists():
            return cand
    raise FileNotFoundError(f"Could not find m1-qc-cleaned.h5ad in {H5AD_CANDIDATES}")


def _excluded_gene_mask(gene_names: pd.Index) -> np.ndarray:
    genes = gene_names.astype(str)
    return np.array(
        [bool(re.match(r"^(mt-|Rp[ls]|Mrp[ls])", g)) or g in SEX_GENES for g in genes],
        dtype=bool,
    )


def _rarefy_csr(matrix: sp.csr_matrix, target_umi: int, rng: np.random.Generator) -> sp.csr_matrix:
    """Exact without-replacement subsampling to target_umi UMIs per qualifying cell."""
    n_rows, n_cols = matrix.shape
    indptr = np.zeros(n_rows + 1, dtype=np.int64)
    idx_parts: list[np.ndarray] = []
    val_parts: list[np.ndarray] = []
    cur = 0
    for i in range(n_rows):
        start, end = matrix.indptr[i], matrix.indptr[i + 1]
        cols = matrix.indices[start:end]
        vals = np.rint(matrix.data[start:end]).astype(np.int64)
        tot = int(vals.sum())
        if tot == target_umi:
            keep = vals > 0
            c_keep, v_keep = cols[keep], vals[keep].astype(np.float32)
        elif tot > target_umi:
            sampled = rng.multivariate_hypergeometric(vals, target_umi)
            keep = sampled > 0
            c_keep, v_keep = cols[keep], sampled[keep].astype(np.float32)
        else:
            c_keep = np.empty(0, dtype=np.int32)
            v_keep = np.empty(0, dtype=np.float32)
        idx_parts.append(c_keep)
        val_parts.append(v_keep)
        cur += len(c_keep)
        indptr[i + 1] = cur
    if cur == 0:
        return sp.csr_matrix((n_rows, n_cols), dtype=np.float32)
    return sp.csr_matrix(
        (np.concatenate(val_parts), np.concatenate(idx_parts), indptr),
        shape=(n_rows, n_cols),
        dtype=np.float32,
    )


def _donor_pseudobulk_and_dispersion(
    ct_mat: sp.csr_matrix,
    ct_rare: sp.csr_matrix | None,
    donor_vec: np.ndarray,
    age_vec: np.ndarray,
    sex_vec: np.ndarray,
    umi_vec: np.ndarray,
    donors_keep: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]:
    """Compute pseudobulk counts, detected mask counts, depth, age, sex, rarefied pseudobulk, LOO cell Hellinger dispersion."""
    n_d = len(donors_keep)
    n_g = ct_mat.shape[1]
    pb_full = np.zeros((n_d, n_g), dtype=np.float64)
    pb_rare = np.zeros((n_d, n_g), dtype=np.float64) if ct_rare is not None else None
    is_30 = np.zeros(n_d, dtype=bool)
    is_male = np.zeros(n_d, dtype=bool)
    mean_log_umi = np.zeros(n_d, dtype=np.float64)
    disp_loo = np.full(n_d, np.nan, dtype=np.float64)

    for di, d in enumerate(donors_keep):
        m = donor_vec == d
        sub = ct_mat[m]
        pb_full[di] = np.asarray(sub.sum(axis=0), dtype=np.float64).ravel()
        is_30[di] = age_vec[m][0] == "30m"
        is_male[di] = sex_vec[m][0].lower() == "male"
        mean_log_umi[di] = float(np.mean(np.log1p(umi_vec[m])))
        if ct_rare is not None:
            rsub = ct_rare[m]
            row_sums = np.asarray(rsub.sum(axis=1)).ravel()
            valid = np.isclose(row_sums, RAREFY_UMI)
            n_r = int(valid.sum())
            if n_r >= 2:
                rsub_v = rsub[valid]
                pb_r = np.asarray(rsub_v.sum(axis=0), dtype=np.float64).ravel()
                pb_rare[di] = pb_r
                if n_r >= 5:
                    # Leave-one-out cell centroid Hellinger distance via sparse nonzero indices
                    p_tot = np.asarray(rsub_v.sum(axis=0), dtype=np.float64).ravel() / float(RAREFY_UMI)
                    p_nz = rsub_v.data.astype(np.float64) / float(RAREFY_UMI)
                    p_loo_nz = np.maximum(p_tot[rsub_v.indices] - p_nz, 0.0) / float(n_r - 1)
                    bc_csr = sp.csr_matrix(
                        (np.sqrt(p_nz * p_loo_nz), rsub_v.indices, rsub_v.indptr),
                        shape=rsub_v.shape,
                    )
                    bc = np.clip(np.asarray(bc_csr.sum(axis=1)).ravel(), 0.0, 1.0)
                    hell = np.sqrt(np.maximum(0.0, 1.0 - bc))
                    disp_loo[di] = float(np.mean(hell))

    return pb_full, is_30, is_male, mean_log_umi, disp_loo, pb_rare, np.array(donors_keep, dtype=object)


def _residualize_against_depth(expr: np.ndarray, mean_log_umi: np.ndarray) -> np.ndarray:
    """Linearly residualize each gene against mean_log_umi across donors in the stratum."""
    x = mean_log_umi - np.mean(mean_log_umi)
    denom = float(np.dot(x, x))
    if denom <= 1e-12:
        return expr - np.mean(expr, axis=0, keepdims=True)
    y_c = expr - np.mean(expr, axis=0, keepdims=True)
    beta = (x[:, None] * y_c).sum(axis=0) / denom
    resid = y_c - x[:, None] * beta[None, :]
    return resid


def _df_logvar_bias(n: int) -> float:
    """Exact Gaussian E[log(s^2)] - log(sigma^2) = psi(df/2) - log(df/2) for sample size n."""
    df = float(n - 1)
    return float(digamma(0.5 * df) - np.log(0.5 * df))


def _stratum_stats(expr: np.ndarray, mask30: np.ndarray, eps: float = 1e-3) -> dict:
    """Compute per-gene ANOVA F, eta^2, digamma-debiased log(s30^2 / s18^2), and exact animal-label permutations."""
    n_d, n_g = expr.shape
    n30 = int(mask30.sum())
    n18 = n_d - n30
    x18 = expr[~mask30]
    x30 = expr[mask30]
    bias_diff = _df_logvar_bias(n30) - _df_logvar_bias(n18)

    m18 = x18.mean(axis=0)
    m30 = x30.mean(axis=0)
    mtot = expr.mean(axis=0)
    ss_between = n18 * (m18 - mtot) ** 2 + n30 * (m30 - mtot) ** 2
    var18 = x18.var(axis=0, ddof=1)
    var30 = x30.var(axis=0, ddof=1)
    ss_within = (n18 - 1) * var18 + (n30 - 1) * var30
    ss_total = np.maximum(ss_between + ss_within, 1e-8)

    eta2_g = ss_between / ss_total
    f_g = (ss_between / 1.0) / np.maximum(ss_within / float(n_d - 2), 1e-8)
    lvr_g = np.log((var30 + eps) / (var18 + eps)) - bias_diff

    p_f_g = f_dist.sf(f_g, 1, n_d - 2)
    q_f_g = multipletests(p_f_g, method="fdr_bh")[1]

    obs_mean_eta2 = float(np.mean(eta2_g))
    obs_med_eta2 = float(np.median(eta2_g))
    obs_mean_f = float(np.mean(f_g))
    obs_mean_lvr = float(np.mean(lvr_g))
    obs_med_lvr = float(np.median(lvr_g))
    obs_tot_lvr = float(np.log(max(float(np.mean(var30)), 1e-8) / max(float(np.mean(var18)), 1e-8)))
    obs_frac_30_gt_18 = float(np.mean(var30 > var18))

    # Exact animal-label permutations across all binom(n_d, n30) combinations
    total_ss_all = np.maximum(np.sum((expr - mtot[None, :]) ** 2, axis=0), 1e-8)
    sum_all = expr.sum(axis=0)
    sum_sq_all = np.sum(expr ** 2, axis=0)
    perm_eta2_vals = []
    perm_lvr_vals = []
    perm_tot_lvr_vals = []
    for combo in itertools.combinations(range(n_d), n30):
        idx30 = list(combo)
        sub30 = expr[idx30]
        s30 = sub30.sum(axis=0)
        s18 = sum_all - s30
        p_m30 = s30 / float(n30)
        p_m18 = s18 / float(n18)
        p_ssb = n18 * (p_m18 - mtot) ** 2 + n30 * (p_m30 - mtot) ** 2
        perm_eta2_vals.append(float(np.mean(p_ssb / total_ss_all)))

        sq30 = np.sum(sub30 ** 2, axis=0)
        sq18 = sum_sq_all - sq30
        p_v30 = np.maximum((sq30 - s30 ** 2 / float(n30)) / float(n30 - 1), 0.0)
        p_v18 = np.maximum((sq18 - s18 ** 2 / float(n18)) / float(n18 - 1), 0.0)
        perm_lvr_vals.append(float(np.mean(np.log((p_v30 + eps) / (p_v18 + eps))) - bias_diff))
        perm_tot_lvr_vals.append(float(np.log(max(float(np.mean(p_v30)), 1e-8) / max(float(np.mean(p_v18)), 1e-8))))

    perm_eta2 = np.asarray(perm_eta2_vals, dtype=float)
    perm_lvr = np.asarray(perm_lvr_vals, dtype=float)
    perm_tot_lvr = np.asarray(perm_tot_lvr_vals, dtype=float)
    n_perm = len(perm_eta2)

    p_eta2_exact = float(np.mean(perm_eta2 >= obs_mean_eta2 - 1e-12))
    p_lvr_gt0_exact = float(np.mean(perm_lvr >= obs_mean_lvr - 1e-12))
    p_lvr_lt0_exact = float(np.mean(perm_lvr <= obs_mean_lvr + 1e-12))
    p_lvr_two_exact = min(1.0, 2.0 * min(p_lvr_gt0_exact, p_lvr_lt0_exact))
    p_tot_lvr_gt0 = float(np.mean(perm_tot_lvr >= obs_tot_lvr - 1e-12))
    p_tot_lvr_lt0 = float(np.mean(perm_tot_lvr <= obs_tot_lvr + 1e-12))
    p_tot_lvr_two = min(1.0, 2.0 * min(p_tot_lvr_gt0, p_tot_lvr_lt0))
    null_eta2_sd = float(np.std(perm_eta2, ddof=1)) if n_perm > 1 else np.nan
    null_lvr_sd = float(np.std(perm_lvr, ddof=1)) if n_perm > 1 else np.nan

    return {
        "n_genes": int(n_g),
        "mean_eta2": obs_mean_eta2,
        "median_eta2": obs_med_eta2,
        "mean_F": obs_mean_f,
        "frac_genes_F_param_q05": float(np.mean(q_f_g < 0.05)),
        "mean_log_var_ratio_30_vs_18": obs_mean_lvr,
        "median_log_var_ratio_30_vs_18": obs_med_lvr,
        "log_total_var_ratio_30_vs_18": obs_tot_lvr,
        "frac_genes_s30_gt_s18": obs_frac_30_gt_18,
        "n_exact_perm": int(n_perm),
        "perm_null_mean_eta2": float(np.mean(perm_eta2)),
        "perm_null_sd_eta2": null_eta2_sd,
        "z_eta2_perm": (obs_mean_eta2 - float(np.mean(perm_eta2))) / null_eta2_sd if null_eta2_sd > 0 else np.nan,
        "p_eta2_exact": p_eta2_exact,
        "perm_null_mean_lvr": float(np.mean(perm_lvr)),
        "perm_null_sd_lvr": null_lvr_sd,
        "z_lvr_perm": (obs_mean_lvr - float(np.mean(perm_lvr))) / null_lvr_sd if null_lvr_sd > 0 else np.nan,
        "p_lvr_gt0_exact": p_lvr_gt0_exact,
        "p_lvr_lt0_exact": p_lvr_lt0_exact,
        "p_lvr_two_sided_exact": p_lvr_two_exact,
        "p_tot_lvr_gt0_exact": p_tot_lvr_gt0,
        "p_tot_lvr_two_sided_exact": p_tot_lvr_two,
    }


def _disp_exact_perm(disp_vals: np.ndarray, mask30: np.ndarray) -> dict:
    valid = np.isfinite(disp_vals)
    vals = disp_vals[valid]
    m30 = mask30[valid]
    n30 = int(m30.sum())
    n18 = len(vals) - n30
    if n30 < 3 or n18 < 3:
        return {}
    diff = float(np.mean(vals[m30]) - np.mean(vals[~m30]))
    perms = []
    tot = float(np.sum(vals))
    for combo in itertools.combinations(range(len(vals)), n30):
        s30 = float(np.sum(vals[list(combo)]))
        s18 = tot - s30
        perms.append((s30 / float(n30)) - (s18 / float(n18)))
    perms_arr = np.asarray(perms, dtype=float)
    p_gt = float(np.mean(perms_arr >= diff - 1e-12))
    p_lt = float(np.mean(perms_arr <= diff + 1e-12))
    return {
        "n18_disp": int(n18),
        "n30_disp": int(n30),
        "mean_disp_18m": float(np.mean(vals[~m30])),
        "mean_disp_30m": float(np.mean(vals[m30])),
        "diff_30m_minus_18m_disp": diff,
        "p_disp_gt0_exact": p_gt,
        "p_disp_two_sided_exact": min(1.0, 2.0 * min(p_gt, p_lt)),
    }


def _apply_bh_within_arm(df: pd.DataFrame) -> pd.DataFrame:
    out_parts = []
    for arm, part in df.groupby("arm", sort=False):
        part = part.copy()
        for pcol, qcol in [
            ("p_eta2_exact", "q_bh_eta2_within_arm"),
            ("p_lvr_gt0_exact", "q_bh_lvr_gt0_within_arm"),
            ("p_lvr_lt0_exact", "q_bh_lvr_lt0_within_arm"),
            ("p_lvr_two_sided_exact", "q_bh_lvr_two_within_arm"),
            ("p_tot_lvr_two_sided_exact", "q_bh_tot_lvr_two_within_arm"),
            ("p_disp_gt0_exact", "q_bh_disp_gt0_within_arm"),
            ("p_disp_two_sided_exact", "q_bh_disp_two_within_arm"),
        ]:
            part[qcol] = np.nan
            if pcol in part.columns:
                valid = part[pcol].notna()
                if valid.sum() > 0:
                    part.loc[valid, qcol] = multipletests(part.loc[valid, pcol], method="fdr_bh")[1]
        out_parts.append(part)
    return pd.concat(out_parts, ignore_index=True)


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    h5_path = _resolve_h5ad()
    print(f"[item8] Reading metadata from {h5_path}", flush=True)
    adata = ad.read_h5ad(h5_path, backed="r")
    obs = adata.obs.copy()
    keep_gene_mask = ~_excluded_gene_mask(adata.var_names)

    obs["age_str"] = obs["age"].astype(str)
    obs["mouse_str"] = obs["mouse.id"].astype(str)
    obs["sex_str"] = obs["sex"].astype(str)
    obs["tissue_str"] = obs["tissue"].astype(str)
    obs["ct_str"] = obs["cell_ontology_class"].astype(str)
    obs["pred_dbl"] = obs["predicted_doublet"].astype(str) == "True"

    base_mask = (
        obs["age_str"].isin(["18m", "30m"])
        & ~obs["pred_dbl"]
        & (obs["n_counts"] >= 500)
        & (obs["ct_str"] != "unknown")
    )
    sub_obs = obs.loc[base_mask].copy()

    # Determine all eligible strata (>=3 mice 18m and >=3 mice 30m at >=10 cells/mouse)
    all_18_mice = sorted(sub_obs.loc[sub_obs["age_str"] == "18m", "mouse_str"].unique())
    all_30_mice = sorted(sub_obs.loc[sub_obs["age_str"] == "30m", "mouse_str"].unique())
    all_8_mice = all_18_mice + all_30_mice

    cnt_tbl = (
        sub_obs.groupby(["tissue_str", "ct_str", "mouse_str"], observed=True)
        .size()
        .unstack("mouse_str", fill_value=0)
        .reindex(columns=all_8_mice, fill_value=0)
    )
    elig10 = cnt_tbl[((cnt_tbl[all_18_mice] >= 10).sum(1) >= 3) & ((cnt_tbl[all_30_mice] >= 10).sum(1) >= 3)]
    tissues_elig = sorted(elig10.index.get_level_values(0).unique())
    print(f"[item8] Eligible >=10c strata: {len(elig10)} across {tissues_elig}", flush=True)

    rows: list[dict] = []
    # Store prepared matrices per primary stratum for exact joint 8-mouse global permutations
    joint_cache: list[tuple[str, str, list[str], np.ndarray]] = []
    joint_cache_rare: list[tuple[str, str, list[str], np.ndarray]] = []

    rng_global = np.random.default_rng(SEED)

    for tissue in tissues_elig:
        t_pos = np.flatnonzero((base_mask & (obs["tissue_str"] == tissue)).to_numpy())
        print(f"[item8] Loading {tissue}: {len(t_pos)} cells into sparse CSR", flush=True)
        raw_sub = adata[t_pos, :].to_memory()
        mat_full = sp.csr_matrix(raw_sub.X, dtype=np.float32)[:, keep_gene_mask]
        t_obs = obs.iloc[t_pos].copy()
        del raw_sub
        gc.collect()

        # Cell filtered UMI sums & rarefied 1500-UMI sparse matrix for cells with filtered UMI >= 1500
        cell_umi_filt = np.asarray(mat_full.sum(axis=1)).ravel()
        mat_rare = _rarefy_csr(mat_full, RAREFY_UMI, rng_global)
        cell_umi_raw = t_obs["n_counts"].to_numpy(dtype=float)

        strata_t = [ct for (ti, ct) in elig10.index if ti == tissue]
        for ct in strata_t:
            ct_mask = (t_obs["ct_str"] == ct).to_numpy()
            ct_obs = t_obs.iloc[ct_mask]
            ct_mat = mat_full[ct_mask]
            ct_rare = mat_rare[ct_mask]
            ct_d = ct_obs["mouse_str"].to_numpy()
            ct_a = ct_obs["age_str"].to_numpy()
            ct_s = ct_obs["sex_str"].to_numpy()
            ct_u = cell_umi_raw[ct_mask]
            ct_u_filt = cell_umi_filt[ct_mask]

            # Check donor eligibility for >=10 and >=20 cells
            d_counts = pd.Series(ct_d).value_counts()
            for min_c, arm_prefix in [(10, "primary_10cells"), (20, "sensitivity_20cells")]:
                keep_d = [d for d in all_8_mice if d_counts.get(d, 0) >= min_c]
                n18_d = sum(d in all_18_mice for d in keep_d)
                n30_d = sum(d in all_30_mice for d in keep_d)
                if n18_d < 3 or n30_d < 3:
                    continue

                pb_full, is_30, is_male, mlu, disp_loo, pb_rare, d_arr = _donor_pseudobulk_and_dispersion(
                    ct_mat, ct_rare if min_c == 10 else None, ct_d, ct_a, ct_s, ct_u, keep_d
                )
                det_mask = np.sum(pb_full > 0, axis=0) >= (0.5 * len(keep_d))
                pb_det = pb_full[:, det_mask]
                tot = pb_det.sum(axis=1, keepdims=True)
                expr_log = np.log1p(pb_det / np.maximum(tot, 1.0) * 1e4)
                expr_resid = _residualize_against_depth(expr_log, mlu)

                base_meta = {
                    "tissue": tissue,
                    "cell_type": ct,
                    "min_cells_per_mouse": int(min_c),
                    "n_donors_18m": int(n18_d),
                    "n_donors_30m": int(n30_d),
                    "n_cells_total": int(sum(d_counts[d] for d in keep_d)),
                }

                # Residualized full-depth arm (primary when min_c==10)
                st_res = _stratum_stats(expr_resid, is_30)
                rows.append({"arm": f"{arm_prefix}_full_residualized", **base_meta, **st_res})

                if min_c == 10:
                    joint_cache.append((tissue, ct, keep_d, expr_resid))
                    # Unresidualized full-depth sensitivity
                    st_un = _stratum_stats(expr_log, is_30)
                    rows.append({"arm": "sensitivity_10cells_full_unresidualized", **base_meta, **st_un})

                    # Male-only sensitivity (18-M-52, 18-M-53 vs >=3 30m males)
                    male_idx = np.flatnonzero(is_male)
                    if int((~is_30[male_idx]).sum()) >= 2 and int((is_30[male_idx]).sum()) >= 3:
                        pb_m = pb_full[male_idx]
                        det_m = np.sum(pb_m > 0, axis=0) >= (0.5 * len(male_idx))
                        pb_md = pb_m[:, det_m]
                        tot_m = pb_md.sum(axis=1, keepdims=True)
                        expr_m = np.log1p(pb_md / np.maximum(tot_m, 1.0) * 1e4)
                        expr_m_res = _residualize_against_depth(expr_m, mlu[male_idx])
                        st_m = _stratum_stats(expr_m_res, is_30[male_idx])
                        rows.append({
                            "arm": "sensitivity_10cells_male_only_residualized",
                            "tissue": tissue,
                            "cell_type": ct,
                            "min_cells_per_mouse": 10,
                            "n_donors_18m": int((~is_30[male_idx]).sum()),
                            "n_donors_30m": int((is_30[male_idx]).sum()),
                            "n_cells_total": int(sum(d_counts[keep_d[i]] for i in male_idx)),
                            **st_m,
                        })

                    # Rarefied 1500-UMI pseudobulk & within-mouse cell-to-centroid Hellinger dispersion
                    r_valid_counts = pd.Series(ct_d[ct_u_filt >= RAREFY_UMI]).value_counts()
                    r_keep_idx = [i for i, d in enumerate(keep_d) if r_valid_counts.get(d, 0) >= 10]
                    if (
                        len(r_keep_idx) >= 6
                        and int((~is_30[r_keep_idx]).sum()) >= 3
                        and int((is_30[r_keep_idx]).sum()) >= 3
                    ):
                        r_keep_d = [keep_d[i] for i in r_keep_idx]
                        pbr = pb_rare[r_keep_idx]
                        det_r = np.sum(pbr > 0, axis=0) >= (0.5 * len(r_keep_idx))
                        pbr_d = pbr[:, det_r]
                        tot_r = pbr_d.sum(axis=1, keepdims=True)
                        expr_r = np.log1p(pbr_d / np.maximum(tot_r, 1.0) * 1e4)
                        st_r = _stratum_stats(expr_r, is_30[r_keep_idx])
                        disp_meta = _disp_exact_perm(disp_loo[r_keep_idx], is_30[r_keep_idx])
                        rows.append({
                            "arm": "sensitivity_10cells_rarefied_1500_pseudobulk",
                            "tissue": tissue,
                            "cell_type": ct,
                            "min_cells_per_mouse": 10,
                            "n_donors_18m": int((~is_30[r_keep_idx]).sum()),
                            "n_donors_30m": int((is_30[r_keep_idx]).sum()),
                            "n_cells_total": int(sum(r_valid_counts[d] for d in r_keep_d)),
                            **st_r,
                            **disp_meta,
                        })
                        joint_cache_rare.append((tissue, ct, r_keep_d, expr_r))

        del mat_full, mat_rare, t_obs
        gc.collect()

    adata.file.close()
    res_df = _apply_bh_within_arm(pd.DataFrame(rows))
    res_df.to_csv(RESULTS / "rev3-randomness-variability-test.csv", index=False)

    # Joint 8-mouse global exact permutation over all binom(8, 4) = 70 global age partitions
    joint_rows = []
    for label, cache in [("primary_full_residualized_34", joint_cache), ("rarefied_1500_pseudobulk", joint_cache_rare)]:
        obs_m_eta2 = []
        obs_m_lvr = []
        for _, _, keep_d, expr in cache:
            m30 = np.array([d in all_30_mice for d in keep_d])
            st = _stratum_stats(expr, m30)
            obs_m_eta2.append(st["mean_eta2"])
            obs_m_lvr.append(st["mean_log_var_ratio_30_vs_18"])
        perm_global_eta2 = []
        perm_global_lvr = []
        eps_lvr = 1e-3
        for combo in itertools.combinations(range(len(all_8_mice)), 4):
            fake30 = {all_8_mice[i] for i in combo}
            eta_list = []
            lvr_list = []
            for _, _, keep_d, expr in cache:
                m30 = np.array([d in fake30 for d in keep_d])
                n30 = int(m30.sum())
                n18 = len(m30) - n30
                if n30 < 2 or n18 < 2:
                    continue
                x18, x30 = expr[~m30], expr[m30]
                mtot = expr.mean(axis=0)
                ss_b = n18 * (x18.mean(0) - mtot) ** 2 + n30 * (x30.mean(0) - mtot) ** 2
                v18, v30 = x18.var(0, ddof=1), x30.var(0, ddof=1)
                ss_t = np.maximum(ss_b + (n18 - 1) * v18 + (n30 - 1) * v30, 1e-8)
                eta_list.append(float(np.mean(ss_b / ss_t)))
                b_diff = _df_logvar_bias(n30) - _df_logvar_bias(n18)
                lvr_list.append(float(np.mean(np.log((v30 + eps_lvr) / (v18 + eps_lvr))) - b_diff))
            perm_global_eta2.append(float(np.mean(eta_list)))
            perm_global_lvr.append(float(np.mean(lvr_list)))

        p_g_eta2 = np.asarray(perm_global_eta2, dtype=float)
        p_g_lvr = np.asarray(perm_global_lvr, dtype=float)
        obs_e = float(np.mean(obs_m_eta2))
        obs_l = float(np.mean(obs_m_lvr))
        p_l_gt = float(np.mean(p_g_lvr >= obs_l - 1e-12))
        p_l_lt = float(np.mean(p_g_lvr <= obs_l + 1e-12))
        joint_rows.append({
            "arm": label,
            "n_strata": len(cache),
            "n_global_mouse_perms": len(p_g_eta2),
            "observed_mean_eta2": obs_e,
            "joint_null_mean_eta2": float(np.mean(p_g_eta2)),
            "joint_null_sd_eta2": float(np.std(p_g_eta2, ddof=1)),
            "joint_p_eta2_exact": float(np.mean(p_g_eta2 >= obs_e - 1e-12)),
            "observed_mean_log_var_ratio": obs_l,
            "joint_null_mean_lvr": float(np.mean(p_g_lvr)),
            "joint_null_sd_lvr": float(np.std(p_g_lvr, ddof=1)),
            "joint_p_lvr_gt0_exact": p_l_gt,
            "joint_p_lvr_lt0_exact": p_l_lt,
            "joint_p_lvr_two_sided_exact": min(1.0, 2.0 * min(p_l_gt, p_l_lt)),
        })

    joint_df = pd.DataFrame(joint_rows)
    joint_df.to_csv(RESULTS / "rev3-randomness-variability-joint-global.csv", index=False)
    print("\n=== ITEM 8 ARM SUMMARY ===", flush=True)
    for arm, g in res_df.groupby("arm", sort=False):
        print(
            f"{arm}: n_strata={len(g)}, med_eta2={g['mean_eta2'].median():.4f}, "
            f"min_q_eta2={g['q_bh_eta2_within_arm'].min():.4f}, "
            f"med_lvr={g['mean_log_var_ratio_30_vs_18'].median():+.4f}, "
            f"frac_s30_gt_s18={g['frac_genes_s30_gt_s18'].median():.3f}, "
            f"min_q_lvr_gt0={g['q_bh_lvr_gt0_within_arm'].min():.4f}, "
            f"min_q_lvr_lt0={g['q_bh_lvr_lt0_within_arm'].min():.4f}",
            flush=True,
        )
    print("\n=== JOINT 8-MOUSE PERMUTATION ===", flush=True)
    print(joint_df.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
