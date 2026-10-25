#!/usr/bin/env python3
"""
scripts/run_h3_entropy.py

Tests Hypothesis 3: Cellular Entropy and Transcriptional Noise in Aging.
Quantifies within-cell-type cell-to-cell variability across the lifespan
using Distance-to-Median (DM; Martinez-Jimenez et al. Science 2017)
and multivariate dispersion (Anderson betadisper).

Outputs:
  - results/h3-entropy-mouse-lung.csv
  - figures/h3-transcriptional-noise-vs-age.png
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import binned_statistic
import statsmodels.api as sm
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


def compute_donor_dispersion(mat, baseline_func, n_top_genes=2000):
    """
    Computes:
    1. Distance-to-Median (DM): deviation of log(CV^2) from the common reference baseline.
    2. Multivariate distance-to-centroid in normalized expression space.
    3. Shannon entropy per cell (mean across cells).
    """
    n_cells, n_genes = mat.shape
    if n_cells < 20:
        return None

    # Normalization: counts per 10k
    totals = mat.sum(axis=1, keepdims=True)
    norm_mat = (mat / np.maximum(totals, 1.0)) * 1e4
    log_mat = np.log1p(norm_mat)

    # 1. Gene-level statistics
    means = norm_mat.mean(axis=0)
    vars_ = norm_mat.var(axis=0, ddof=1)
    cv2 = vars_ / np.maximum(means**2, 1e-6)

    # Filter expressed genes
    det_rate = (mat > 0).mean(axis=0)
    expr_mask = (det_rate >= 0.05) & (means > 0.01)
    if expr_mask.sum() < 50:
        return None

    sub_means = means[expr_mask]
    sub_cv2 = cv2[expr_mask]
    log_means = np.log10(sub_means)
    log_cv2 = np.log10(np.maximum(sub_cv2, 1e-4))

    # Evaluate against common baseline
    expected_log_cv2 = baseline_func(log_means)
    dm = log_cv2 - expected_log_cv2
    mean_dm = float(np.mean(dm))
    median_dm = float(np.median(dm))

    # 2. Multivariate distance to centroid (in log space)
    centroid = log_mat[:, expr_mask].mean(axis=0)
    dists = np.linalg.norm(log_mat[:, expr_mask] - centroid, axis=1)
    mean_dist = float(dists.mean())

    # 3. Cell-level Shannon entropy of gene expression
    p_mat = norm_mat[:, expr_mask] / np.maximum(norm_mat[:, expr_mask].sum(axis=1, keepdims=True), 1e-6)
    log_p = np.log2(np.maximum(p_mat, 1e-12))
    shannon_per_cell = -np.sum(p_mat * log_p, axis=1)
    mean_shannon = float(shannon_per_cell.mean())

    return {
        "mean_dm": mean_dm,
        "median_dm": median_dm,
        "mean_euclidean_dist": mean_dist,
        "shannon_entropy": mean_shannon,
        "n_cells": n_cells,
        "mean_umi_per_cell": float(totals.mean())
    }


def main():
    print("[H3 Entropy] Starting Phase 2: Cellular Entropy & Transcriptional Noise...")
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    import anndata as ad
    print("[H3 Entropy] Loading mouse lung single-cell data from data/interim/m1-processed.h5ad ...")
    adata = ad.read_h5ad(REPO_ROOT / "data" / "interim" / "m1-processed.h5ad", backed="r")
    obs = adata.obs
    lung_obs = obs[obs["tissue"] == "Lung"].copy()

    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    lung_obs["age_months"] = lung_obs["age"].astype(str).map(age_map).astype(float)
    lung_obs["sex_bin"] = (lung_obs["sex"].astype(str).str.lower() == "male").astype(float)

    # Candidate cell types with good cell counts across donors
    ct_counts = lung_obs["cell_ontology_class"].value_counts()
    target_cts = [
        "classical monocyte",
        "intermediate monocyte",
        "bronchial smooth muscle cell",
        "alveolar macrophage",
        "lung macrophage",
        "B cell",
        "NK cell",
        "CD8-positive, alpha-beta T cell",
        "fibroblast of lung",
        "vein endothelial cell"
    ]

    all_celltype_results = []
    donor_records = []

    for ct in target_cts:
        print(f"\n[H3 Entropy] Processing cell type: {ct} ...")
        ct_obs = lung_obs[lung_obs["cell_ontology_class"] == ct]
        mice_with_cells = ct_obs["mouse.id"].value_counts()[lambda x: x >= 25].index.tolist()
        if len(mice_with_cells) < 8:
            print(f"  Skipping {ct}: only {len(mice_with_cells)} mice with >= 25 cells (need >= 8).")
            continue

        # Step 1: Fit common baseline curve across all young mice (1m, 3m) or all mice if young sparse
        young_mice = ct_obs[ct_obs["age_months"] <= 3.0]["mouse.id"].unique()
        ref_mice = young_mice if len(young_mice) >= 3 else mice_with_cells

        ref_cells = ct_obs[ct_obs["mouse.id"].isin(ref_mice)].index
        ref_row_idx = np.where(obs.index.isin(ref_cells))[0]
        
        # Subsample up to 2000 cells for reference baseline
        if len(ref_row_idx) > 2000:
            np.random.seed(42)
            ref_row_idx = np.random.choice(ref_row_idx, size=2000, replace=False)

        sub_X = adata.X[ref_row_idx]
        ref_mat = sub_X.toarray() if sp.issparse(sub_X) else np.array(sub_X)
        ref_totals = ref_mat.sum(axis=1, keepdims=True)
        ref_norm = (ref_mat / np.maximum(ref_totals, 1.0)) * 1e4
        ref_means = ref_norm.mean(axis=0)
        ref_vars = ref_norm.var(axis=0, ddof=1)
        ref_cv2 = ref_vars / np.maximum(ref_means**2, 1e-6)

        expr_mask = (ref_means > 0.01) & ((ref_mat > 0).mean(axis=0) >= 0.05)
        log_m = np.log10(ref_means[expr_mask])
        log_c = np.log10(np.maximum(ref_cv2[expr_mask], 1e-4))

        # 15-bin running median
        bin_meds, bin_edges, bin_num = binned_statistic(log_m, log_c, statistic="median", bins=15)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        def baseline_func(query_log_m):
            return np.interp(query_log_m, bin_centers, bin_meds, left=bin_meds[0], right=bin_meds[-1])

        # Step 2: Calculate dispersion per donor
        ct_donor_rows = []
        for m in mice_with_cells:
            donor_cells = ct_obs[ct_obs["mouse.id"] == m].index
            row_idx = np.where(obs.index.isin(donor_cells))[0]
            # Subsample up to 500 cells to standardize donor size influence
            if len(row_idx) > 500:
                np.random.seed(42)
                row_idx = np.random.choice(row_idx, size=500, replace=False)

            donor_X = adata.X[row_idx]
            donor_mat = donor_X.toarray() if sp.issparse(donor_X) else np.array(donor_X)

            disp = compute_donor_dispersion(donor_mat, baseline_func)
            if disp is None:
                continue

            age_val = ct_obs[ct_obs["mouse.id"] == m]["age_months"].iloc[0]
            sex_val = ct_obs[ct_obs["mouse.id"] == m]["sex_bin"].iloc[0]

            disp.update({
                "cell_type": ct,
                "mouse_id": m,
                "age_months": age_val,
                "sex_bin": sex_val
            })
            ct_donor_rows.append(disp)
            donor_records.append(disp)

        donor_df = pd.DataFrame(ct_donor_rows)
        if len(donor_df) < 8:
            continue

        # Step 3: Regress dispersion metrics on Age + Sex
        # Distance-to-Median (DM)
        X = sm.add_constant(donor_df[["age_months", "sex_bin"]])
        res_dm = sm.OLS(donor_df["mean_dm"], X).fit()
        res_dist = sm.OLS(donor_df["mean_euclidean_dist"], X).fit()
        res_shannon = sm.OLS(donor_df["shannon_entropy"], X).fit()

        all_celltype_results.append({
            "cell_type": ct,
            "n_donors": len(donor_df),
            "dm_slope_per_month": res_dm.params["age_months"],
            "dm_se": res_dm.bse["age_months"],
            "dm_pval": res_dm.pvalues["age_months"],
            "dist_slope_per_month": res_dist.params["age_months"],
            "dist_pval": res_dist.pvalues["age_months"],
            "shannon_slope_per_month": res_shannon.params["age_months"],
            "shannon_pval": res_shannon.pvalues["age_months"]
        })

    res_df = pd.DataFrame(all_celltype_results)
    from statsmodels.stats.multitest import multipletests
    _, qvals, _, _ = multipletests(res_df["dm_pval"], method="fdr_bh")
    res_df["dm_q_bh"] = qvals
    res_df = res_df.sort_values("dm_pval").reset_index(drop=True)

    out_csv = results_dir / "h3-entropy-mouse-lung.csv"
    res_df.to_csv(out_csv, index=False)
    print(f"\n[H3 Entropy] Saved mouse lung entropy results to {out_csv}")
    print(res_df[["cell_type", "n_donors", "dm_slope_per_month", "dm_pval", "dm_q_bh", "dist_slope_per_month", "dist_pval"]])

    # -------------------------------------------------------------
    # Generate Publication Figure
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Forest plot of DM slopes
    y_pos = np.arange(len(res_df))
    colors = ["crimson" if q < 0.05 and s > 0 else "gray" for q, s in zip(res_df["dm_q_bh"], res_df["dm_slope_per_month"])]
    axes[0].barh(y_pos, res_df["dm_slope_per_month"], color=colors, alpha=0.8)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(res_df["cell_type"], fontsize=9)
    axes[0].axvline(0, color="black", linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("Distance-to-Median Slope (per month)")
    axes[0].set_title("A. Within-Cell-Type Transcriptional Noise Slope\n(Red: FDR q < 0.05)")

    # Panel B: Scatter plot of top cell type (monocyte or AT2) across age
    donor_all_df = pd.DataFrame(donor_records)
    top_ct = res_df.iloc[0]["cell_type"]
    sub_top = donor_all_df[donor_all_df["cell_type"] == top_ct]
    axes[1].scatter(sub_top["age_months"], sub_top["mean_dm"], color="royalblue", s=60, edgecolors="black", alpha=0.9)
    
    # Fit regression line
    x_grid = np.linspace(sub_top["age_months"].min(), sub_top["age_months"].max(), 50)
    m_slope = res_df.iloc[0]["dm_slope_per_month"]
    intercept = sub_top["mean_dm"].mean() - m_slope * sub_top["age_months"].mean()
    axes[1].plot(x_grid, intercept + m_slope * x_grid, color="firebrick", linewidth=2)
    axes[1].set_xlabel("Donor Age (months)")
    axes[1].set_ylabel("Distance-to-Median (Excess Noise)")
    axes[1].set_title(f"B. {top_ct} (p = {res_df.iloc[0]['dm_pval']:.3f}, q = {res_df.iloc[0]['dm_q_bh']:.3f})")

    plt.tight_layout()
    fig_path = figures_dir / "h3-transcriptional-noise-vs-age.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[H3 Entropy] Saved publication figure to {fig_path}")


if __name__ == "__main__":
    main()
