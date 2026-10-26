#!/usr/bin/env python3
"""
scripts/run_h3_entropy.py

Tests Hypothesis 3: Cellular Entropy and Transcriptional Noise in Aging.
Quantifies within-cell-type cell-to-cell variability across the mouse lifespan
using Distance-to-Median (DM; Martinez-Jimenez et al. Science 2017)
and multivariate dispersion from genuine 1,500-UMI multinomial rarefied counts.

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
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


def main():
    print("[H3 Entropy] Starting Phase 2: Cellular Entropy & Transcriptional Noise...")
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    import anndata as ad
    raw_path = REPO_ROOT / "data" / "interim" / "m1-qc-cleaned.h5ad"
    print(f"[H3 Entropy] Loading QC-cleaned AnnData (raw counts) from {raw_path}...")
    adata = ad.read_h5ad(raw_path, backed="r")

    # Filter to singlet Lung cells with at least 1,500 UMIs
    obs = adata.obs
    lung_mask = (obs["tissue"] == "Lung") & (~obs["predicted_doublet"]) & (obs["n_counts"] >= 1500)
    lung_adata = adata[lung_mask].to_memory()
    print(f"[H3 Entropy] Extracted {lung_adata.shape[0]} qualified lung cells.")

    # -------------------------------------------------------------
    # 1. Exact Multinomial Rarefaction to 1,500 UMIs per cell
    # -------------------------------------------------------------
    print("[H3 Entropy] Performing multinomial downsampling to exactly 1,500 UMIs/cell...")
    rng = np.random.default_rng(42)
    mat = lung_adata.X.tocsr()
    rarefied_data = []
    rarefied_indices = []
    rarefied_indptr = [0]

    for i in range(mat.shape[0]):
        row = mat.getrow(i)
        indices = row.indices
        data_vals = row.data.astype(int)
        tot = data_vals.sum()
        p = data_vals / tot
        sampled = rng.multinomial(1500, p)
        nz = sampled > 0
        rarefied_indices.extend(indices[nz])
        rarefied_data.extend(sampled[nz])
        rarefied_indptr.append(len(rarefied_indices))

    rarefied_csr = sp.csr_matrix((rarefied_data, rarefied_indices, rarefied_indptr), shape=mat.shape)
    lung_adata.layers["rarefied_1500"] = rarefied_csr

    # Verify rarefaction invariant
    umi_sums = np.array(rarefied_csr.sum(axis=1)).flatten()
    assert np.all(umi_sums == 1500), f"Rarefaction invariant failed: {np.unique(umi_sums)}"
    print(f"[H3 Entropy] Rarefaction confirmed: all {len(umi_sums)} cells have exactly 1,500 UMIs.")

    # Demographics mapping
    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    lung_adata.obs["age_months"] = lung_adata.obs["age"].astype(str).map(age_map).astype(float)
    lung_adata.obs["sex_bin"] = (lung_adata.obs["sex"].astype(str).str.lower() == "male").astype(float)

    # Candidate cell types
    target_cts = [
        "fibroblast of lung",
        "bronchial smooth muscle cell",
        "alveolar macrophage",
        "classical monocyte",
        "intermediate monocyte",
        "B cell",
        "NK cell",
        "CD8-positive, alpha-beta T cell",
        "vein endothelial cell"
    ]

    all_celltype_results = []
    donor_records = []

    for ct in target_cts:
        ct_adata = lung_adata[lung_adata.obs["cell_ontology_class"] == ct]
        if ct_adata.n_obs < 100:
            continue

        donors = ct_adata.obs["mouse.id"].unique()
        pooled_mat = ct_adata.layers["rarefied_1500"].toarray()
        p_means = pooled_mat.mean(axis=0)
        p_vars = pooled_mat.var(axis=0, ddof=1)
        p_cv2 = np.where(p_means > 0, p_vars / (p_means**2), np.nan)
        p_det = (pooled_mat > 0).mean(axis=0)

        # Expressed genes: detected in >= 5% of cells and mean >= 0.01 UMI
        g_mask = (p_det >= 0.05) & (p_means >= 0.01) & np.isfinite(p_cv2)
        if g_mask.sum() < 100:
            continue

        sub_log_means = np.log10(p_means[g_mask])
        sub_log_cv2 = np.log10(p_cv2[g_mask])

        # Running median baseline across 25 bins
        bin_medians, bin_edges, _ = binned_statistic(sub_log_means, sub_log_cv2, statistic="median", bins=25)
        bin_mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        valid_bins = ~np.isnan(bin_medians)
        baseline_func = lambda x: np.interp(x, bin_mids[valid_bins], bin_medians[valid_bins])

        d_rows = []
        for d in donors:
            d_cells = ct_adata[ct_adata.obs["mouse.id"] == d]
            if len(d_cells) < 20:
                continue

            d_mat = d_cells.layers["rarefied_1500"].toarray()
            d_means = d_mat.mean(axis=0)
            d_vars = d_mat.var(axis=0, ddof=1)
            d_cv2 = np.where(d_means > 0, d_vars / (d_means**2), np.nan)

            dg_mask = g_mask & (d_means > 0) & np.isfinite(d_cv2) & (d_cv2 > 0)
            dm = np.log10(d_cv2[dg_mask]) - baseline_func(np.log10(d_means[dg_mask]))

            # Multivariate distance to centroid in log2(count + 1) space
            log_mat = np.log2(d_mat[:, g_mask] + 1)
            centroid = log_mat.mean(axis=0)
            dists = np.linalg.norm(log_mat - centroid, axis=1)

            age_m = d_cells.obs["age_months"].iloc[0]
            sex_b = d_cells.obs["sex_bin"].iloc[0]
            d_rows.append({
                "cell_type": ct,
                "donor": d,
                "n_cells": len(d_cells),
                "age_months": age_m,
                "sex_bin": sex_b,
                "mean_dm": float(np.mean(dm)),
                "median_dm": float(np.median(dm)),
                "mean_euclidean_dist": float(np.mean(dists))
            })

        df_d = pd.DataFrame(d_rows)
        if len(df_d) >= 5:
            donor_records.extend(d_rows)
            X = sm.add_constant(df_d[["age_months", "sex_bin"]])

            ols_dm = sm.OLS(df_d["mean_dm"], X).fit()
            ols_med = sm.OLS(df_d["median_dm"], X).fit()
            ols_disp = sm.OLS(df_d["mean_euclidean_dist"], X).fit()

            all_celltype_results.append({
                "cell_type": ct,
                "n_donors": len(df_d),
                "total_cells": int(df_d["n_cells"].sum()),
                "dm_slope_per_month": ols_dm.params["age_months"],
                "dm_se": ols_dm.bse["age_months"],
                "dm_pval": ols_dm.pvalues["age_months"],
                "dm_r2": ols_dm.rsquared,
                "med_dm_slope_per_month": ols_med.params["age_months"],
                "med_dm_pval": ols_med.pvalues["age_months"],
                "dispersion_slope_per_month": ols_disp.params["age_months"],
                "dispersion_se": ols_disp.bse["age_months"],
                "dispersion_pval": ols_disp.pvalues["age_months"],
                "dispersion_r2": ols_disp.rsquared
            })

    res_df = pd.DataFrame(all_celltype_results)

    # Benjamini-Hochberg FDR correction
    res_df["dm_q_bh"] = multipletests(res_df["dm_pval"], method="fdr_bh")[1]
    res_df["dispersion_q_bh"] = multipletests(res_df["dispersion_pval"], method="fdr_bh")[1]

    # Save results
    out_csv = results_dir / "h3-entropy-mouse-lung.csv"
    res_df.to_csv(out_csv, index=False)
    print(f"\n[H3 Entropy] Saved full cell-type entropy results to {out_csv}")
    print(res_df[["cell_type", "n_donors", "dm_slope_per_month", "dm_pval", "dm_q_bh", "dispersion_slope_per_month", "dispersion_pval", "dispersion_q_bh"]].to_string(index=False))

    # -------------------------------------------------------------
    # 2. Publication Figure
    # -------------------------------------------------------------
    df_donors = pd.DataFrame(donor_records)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    plot_cts = ["fibroblast of lung", "classical monocyte", "bronchial smooth muscle cell"]
    colors = {"fibroblast of lung": "crimson", "classical monocyte": "navy", "bronchial smooth muscle cell": "darkorange"}

    for idx, ct in enumerate(plot_cts):
        sub = df_donors[df_donors["cell_type"] == ct]
        ax = axes[idx]
        if len(sub) > 0:
            ax.scatter(sub["age_months"], sub["mean_dm"], color=colors[ct], s=60, alpha=0.8, edgecolors="black")
            fit = np.poly1d(np.polyfit(sub["age_months"], sub["mean_dm"], 1))
            xs = np.linspace(sub["age_months"].min() - 1, sub["age_months"].max() + 1, 100)
            row = res_df[res_df["cell_type"] == ct].iloc[0]
            ax.plot(xs, fit(xs), color=colors[ct], linewidth=2.0,
                    label=f"Slope: {row['dm_slope_per_month']:+.4f}/mo\n(p = {row['dm_pval']:.4f}, q = {row['dm_q_bh']:.4f})")
            ax.set_xlabel("Donor Age (Months)")
            ax.set_ylabel("Transcriptional Noise (Mean DM)")
            ax.set_title(f"{ct.title()} (N = {row['n_donors']} mice)")
            ax.legend(frameon=True, fontsize=9)
            ax.grid(alpha=0.3)

    plt.tight_layout()
    fig_path = figures_dir / "h3-transcriptional-noise-vs-age.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[H3 Entropy] Saved publication figure to {fig_path}")


if __name__ == "__main__":
    main()
