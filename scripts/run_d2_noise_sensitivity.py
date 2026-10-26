#!/usr/bin/env python3
"""
scripts/run_d2_noise_sensitivity.py

D2 Sensitivity Analysis: Single-cell transcriptional noise under
WITHOUT-REPLACEMENT subsampling to exactly 1,500 UMIs per cell.

Compares Distance-to-Median (DM) and multivariate dispersion between:
1. Multinomial resampling (with replacement, adds sampling variance)
2. Hypergeometric subsampling (without replacement, exact count subsetting)

Outputs:
  - results/d2-noise-sensitivity.csv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import anndata as ad
import scanpy as sc

REPO = Path(__file__).resolve().parent.parent

def compute_dm(mat, cell_type_name):
    # mat is cell x gene counts (dense or sparse)
    if sp.issparse(mat):
        mat = mat.toarray()
    n_cells, n_genes = mat.shape
    if n_cells < 20:
        return np.nan, np.nan

    means = np.mean(mat, axis=0)
    vars = np.var(mat, axis=0, ddof=1)
    
    # Filter expressed genes (detection >= 5%, mean >= 0.01)
    det_rate = np.mean(mat > 0, axis=0)
    valid_genes = (det_rate >= 0.05) & (means >= 0.01) & (vars > 0)
    if valid_genes.sum() < 20:
        return np.nan, np.nan

    m_val = means[valid_genes]
    v_val = vars[valid_genes]
    cv2 = v_val / (m_val ** 2)

    log_m = np.log10(m_val)
    log_cv2 = np.log10(cv2)

    # 25-bin running median trend
    bins = np.linspace(log_m.min(), log_m.max(), 26)
    bin_idx = np.digitize(log_m, bins) - 1
    trend = np.zeros_like(log_cv2)
    for b in range(25):
        in_b = bin_idx == b
        if in_b.sum() > 0:
            trend[in_b] = np.median(log_cv2[in_b])
        else:
            trend[in_b] = np.nan

    # Fill any empty bins by interpolation
    s_trend = pd.Series(trend).interpolate(method="linear").bfill().ffill().values
    dm_genes = log_cv2 - s_trend
    mean_dm = np.mean(dm_genes)

    # Multivariate dispersion: mean distance to donor centroid in log2(counts + 1)
    log_mat = np.log2(mat + 1.0)
    centroid = np.mean(log_mat, axis=0)
    dists = np.linalg.norm(log_mat - centroid, axis=1)
    mean_disp = np.mean(dists)

    return mean_dm, mean_disp

def main():
    print("[D2 Noise Sensitivity] Loading QC-cleaned data from data/interim/m1-qc-cleaned.h5ad ...")
    raw_path = REPO / "data" / "interim" / "m1-qc-cleaned.h5ad"
    adata = ad.read_h5ad(raw_path, backed="r")

    obs = adata.obs
    lung_mask = (obs["tissue"] == "Lung") & (~obs["predicted_doublet"]) & (obs["n_counts"] >= 1500)
    lung_adata = adata[lung_mask].to_memory()
    print(f"Qualified cells: {lung_adata.shape[0]}")

    # 1. Multinomial resampling (matching run_h3_entropy.py)
    print("Performing multinomial downsampling to 1,500 UMIs...")
    rng = np.random.default_rng(42)
    mat = lung_adata.X.tocsr()
    multi_data, multi_indices, multi_indptr = [], [], [0]
    for i in range(mat.shape[0]):
        row = mat.getrow(i)
        idx = row.indices
        d = row.data.astype(int)
        tot = d.sum()
        sampled = rng.multinomial(1500, d / tot)
        nz = sampled > 0
        multi_indices.extend(idx[nz])
        multi_data.extend(sampled[nz])
        multi_indptr.append(len(multi_indices))
    multi_csr = sp.csr_matrix((multi_data, multi_indices, multi_indptr), shape=mat.shape)

    # 2. Without-replacement hypergeometric subsampling
    print("Performing without-replacement hypergeometric subsampling to 1,500 UMIs...")
    wo_adata = lung_adata.copy()
    sc.pp.downsample_counts(wo_adata, counts_per_cell=1500, replace=False, random_state=42)
    wo_csr = wo_adata.X.tocsr()

    # Age mapping
    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    lung_adata.obs["age_months"] = lung_adata.obs["age"].astype(str).map(age_map)
    lung_adata.obs["sex_bin"] = (lung_adata.obs["sex"].astype(str).str.lower() == "male").astype(float)

    results_rows = []
    cell_types = [
        "fibroblast of lung",
        "classical monocyte",
        "bronchial smooth muscle cell",
        "alveolar macrophage",
        "B cell",
        "NK cell",
        "CD8-positive, alpha-beta T cell",
        "vein endothelial cell"
    ]

    for ct in cell_types:
        ct_mask = lung_adata.obs["cell_ontology_class"] == ct
        sub_obs = lung_adata.obs[ct_mask]
        donors = sub_obs["mouse.id"].unique()

        donor_records_multi = []
        donor_records_wo = []

        for d in donors:
            d_mask = (sub_obs["mouse.id"] == d).values
            if d_mask.sum() < 20:
                continue

            sub_idx = np.where(ct_mask)[0][d_mask]
            
            # Multi
            m_mat = multi_csr[sub_idx, :]
            dm_m, disp_m = compute_dm(m_mat, ct)
            if not np.isnan(dm_m):
                donor_records_multi.append({
                    "donor": d,
                    "age": sub_obs.loc[sub_obs["mouse.id"] == d, "age_months"].iloc[0],
                    "sex": sub_obs.loc[sub_obs["mouse.id"] == d, "sex_bin"].iloc[0],
                    "dm": dm_m,
                    "disp": disp_m
                })

            # Without replacement
            w_mat = wo_csr[sub_idx, :]
            dm_w, disp_w = compute_dm(w_mat, ct)
            if not np.isnan(dm_w):
                donor_records_wo.append({
                    "donor": d,
                    "age": sub_obs.loc[sub_obs["mouse.id"] == d, "age_months"].iloc[0],
                    "sex": sub_obs.loc[sub_obs["mouse.id"] == d, "sex_bin"].iloc[0],
                    "dm": dm_w,
                    "disp": disp_w
                })

        # Regress Multi
        if len(donor_records_multi) >= 6:
            df_m = pd.DataFrame(donor_records_multi)
            X_m = np.column_stack([np.ones(len(df_m)), df_m["age"], df_m["sex"]])
            b_dm_m, *_ = np.linalg.lstsq(X_m, df_m["dm"], rcond=None)
            res_dm_m = df_m["dm"] - X_m @ b_dm_m
            se_dm_m = np.sqrt(np.diag(((res_dm_m @ res_dm_m)/(len(df_m)-3))*np.linalg.inv(X_m.T @ X_m)))[1]
            p_dm_m = 2 * stats.t.sf(abs(b_dm_m[1]/se_dm_m), len(df_m)-3)

            b_dp_m, *_ = np.linalg.lstsq(X_m, df_m["disp"], rcond=None)
            res_dp_m = df_m["disp"] - X_m @ b_dp_m
            se_dp_m = np.sqrt(np.diag(((res_dp_m @ res_dp_m)/(len(df_m)-3))*np.linalg.inv(X_m.T @ X_m)))[1]
            p_dp_m = 2 * stats.t.sf(abs(b_dp_m[1]/se_dp_m), len(df_m)-3)
        else:
            b_dm_m = [0, np.nan]; se_dm_m = np.nan; p_dm_m = np.nan
            b_dp_m = [0, np.nan]; se_dp_m = np.nan; p_dp_m = np.nan
            df_m = pd.DataFrame()

        # Regress Without Replacement
        if len(donor_records_wo) >= 6:
            df_w = pd.DataFrame(donor_records_wo)
            X_w = np.column_stack([np.ones(len(df_w)), df_w["age"], df_w["sex"]])
            b_dm_w, *_ = np.linalg.lstsq(X_w, df_w["dm"], rcond=None)
            res_dm_w = df_w["dm"] - X_w @ b_dm_w
            se_dm_w = np.sqrt(np.diag(((res_dm_w @ res_dm_w)/(len(df_w)-3))*np.linalg.inv(X_w.T @ X_w)))[1]
            p_dm_w = 2 * stats.t.sf(abs(b_dm_w[1]/se_dm_w), len(df_w)-3)

            b_dp_w, *_ = np.linalg.lstsq(X_w, df_w["disp"], rcond=None)
            res_dp_w = df_w["disp"] - X_w @ b_dp_w
            se_dp_w = np.sqrt(np.diag(((res_dp_w @ res_dp_w)/(len(df_w)-3))*np.linalg.inv(X_w.T @ X_w)))[1]
            p_dp_w = 2 * stats.t.sf(abs(b_dp_w[1]/se_dp_w), len(df_w)-3)
        else:
            b_dm_w = [0, np.nan]; se_dm_w = np.nan; p_dm_w = np.nan
            b_dp_w = [0, np.nan]; se_dp_w = np.nan; p_dp_w = np.nan
            df_w = pd.DataFrame()

        results_rows.append({
            "cell_type": ct,
            "n_donors": len(df_w),
            "dm_slope_multinomial": b_dm_m[1],
            "dm_p_multinomial": p_dm_m,
            "dm_slope_without_replacement": b_dm_w[1],
            "dm_se_without_replacement": se_dm_w,
            "dm_p_without_replacement": p_dm_w,
            "disp_slope_multinomial": b_dp_m[1],
            "disp_p_multinomial": p_dp_m,
            "disp_slope_without_replacement": b_dp_w[1],
            "disp_se_without_replacement": se_dp_w,
            "disp_p_without_replacement": p_dp_w
        })

    out_df = pd.DataFrame(results_rows)
    # Compute FDR q
    out_df["dm_q_without_replacement_bh"] = multipletests(out_df["dm_p_without_replacement"].fillna(1.0), method="fdr_bh")[1]

    out_path = REPO / "results" / "d2-noise-sensitivity.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved D2 noise sensitivity to {out_path}")
    print(out_df[["cell_type", "n_donors", "dm_slope_multinomial", "dm_p_multinomial", "dm_slope_without_replacement", "dm_p_without_replacement", "dm_q_without_replacement_bh"]].to_string())

if __name__ == "__main__":
    main()
