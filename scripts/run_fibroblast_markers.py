#!/usr/bin/env python3
"""
scripts/run_fibroblast_markers.py

Formalized, committed audit of candidate DNA repair, chromatin maintenance,
and cellular senescence markers in mouse lung fibroblasts across lifespan donors.

Regresses donor pseudobulk expression (CP10K) against Age + Sex for 11 targeted markers:
  Parp1, Atm, Xrcc5, Sirt6, Cdkn1a, Cdkn2a, Stat1, Nfkb1, Ccl2, Cxcl12, B2m.
Applies Benjamini-Hochberg FDR correction across all 11 markers.

Outputs:
  - results/fibroblast-markers-audit.csv
  - figures/fibroblast-candidate-markers.png
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


def main():
    print("[Markers Audit] Starting Fibroblast Candidate Markers Evaluation...")
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    import anndata as ad
    raw_path = REPO_ROOT / "data" / "interim" / "m1-qc-cleaned.h5ad"
    adata = ad.read_h5ad(raw_path, backed="r")

    obs = adata.obs
    fibro_mask = (obs["tissue"] == "Lung") & (obs["cell_ontology_class"] == "fibroblast of lung") & (~obs["predicted_doublet"])
    fibro_adata = adata[fibro_mask].to_memory()
    print(f"[Markers Audit] Loaded {fibro_adata.shape[0]} mouse lung fibroblasts.")

    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    fibro_adata.obs["age_months"] = fibro_adata.obs["age"].astype(str).map(age_map).astype(float)
    fibro_adata.obs["sex_bin"] = (fibro_adata.obs["sex"].astype(str).str.lower() == "male").astype(float)

    # Normalize to CP10K for pseudobulk mean
    totals = np.array(fibro_adata.X.sum(axis=1)).flatten()
    norm_x = (fibro_adata.X.toarray() / np.maximum(totals[:, None], 1.0)) * 1e4
    fibro_adata.layers["cp10k"] = norm_x

    markers = ["Cdkn1a", "Sirt6", "Ccl2", "Stat1", "Cxcl12", "Parp1", "B2m", "Nfkb1", "Atm", "Xrcc5", "Cdkn2a"]
    donors = fibro_adata.obs["mouse.id"].unique()

    donor_means = []
    donor_meta = []
    for d in donors:
        d_mask = fibro_adata.obs["mouse.id"] == d
        d_cells = fibro_adata[d_mask]
        if len(d_cells) >= 20:
            d_mean = d_cells.layers["cp10k"].mean(axis=0)
            donor_means.append(d_mean)
            donor_meta.append({
                "donor": d,
                "n_cells": len(d_cells),
                "age_months": d_cells.obs["age_months"].iloc[0],
                "sex_bin": d_cells.obs["sex_bin"].iloc[0]
            })

    df_meta = pd.DataFrame(donor_meta)
    mat_means = np.array(donor_means)
    var_idx = {g: list(fibro_adata.var_names).index(g) for g in markers}

    rows = []
    for g in markers:
        idx = var_idx[g]
        y = mat_means[:, idx]
        X = sm.add_constant(df_meta[["age_months", "sex_bin"]])
        model = sm.OLS(y, X).fit()
        rows.append({
            "gene": g,
            "n_donors": len(df_meta),
            "mean_expression_cp10k": float(np.mean(y)),
            "slope_per_month": float(model.params["age_months"]),
            "se": float(model.bse["age_months"]),
            "t_stat": float(model.tvalues["age_months"]),
            "p_value": float(model.pvalues["age_months"]),
            "r2": float(model.rsquared)
        })

    df_res = pd.DataFrame(rows)
    df_res["q_value_bh"] = multipletests(df_res["p_value"], method="fdr_bh")[1]

    out_csv = results_dir / "fibroblast-markers-audit.csv"
    df_res.to_csv(out_csv, index=False)
    print(f"\n[Markers Audit] Saved audit results to {out_csv}")
    print(df_res[["gene", "mean_expression_cp10k", "slope_per_month", "p_value", "q_value_bh", "r2"]].to_string(index=False))

    # Generate Figure for Cdkn1a and Sirt6
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Cdkn1a
    y_cdkn1a = mat_means[:, var_idx["Cdkn1a"]]
    axes[0].scatter(df_meta["age_months"], y_cdkn1a, color="crimson", s=60, alpha=0.8, edgecolors="black")
    fit_c = np.poly1d(np.polyfit(df_meta["age_months"], y_cdkn1a, 1))
    xs = np.linspace(df_meta["age_months"].min() - 1, df_meta["age_months"].max() + 1, 100)
    row_c = df_res[df_res["gene"] == "Cdkn1a"].iloc[0]
    axes[0].plot(xs, fit_c(xs), color="darkred", linewidth=2.0,
                 label=f"Slope: {row_c['slope_per_month']:+.4f}/mo\n(p = {row_c['p_value']:.4f}, q = {row_c['q_value_bh']:.4f}, R² = {row_c['r2']:.2f})")
    axes[0].set_xlabel("Donor Age (Months)")
    axes[0].set_ylabel("Cdkn1a Expression (CP10K)")
    axes[0].set_title("Mouse Lung Fibroblasts: Cdkn1a (p21)")
    axes[0].legend(frameon=True)
    axes[0].grid(alpha=0.3)

    # Sirt6
    y_sirt6 = mat_means[:, var_idx["Sirt6"]]
    axes[1].scatter(df_meta["age_months"], y_sirt6, color="royalblue", s=60, alpha=0.8, edgecolors="black")
    fit_s = np.poly1d(np.polyfit(df_meta["age_months"], y_sirt6, 1))
    row_s = df_res[df_res["gene"] == "Sirt6"].iloc[0]
    axes[1].plot(xs, fit_s(xs), color="navy", linewidth=2.0,
                 label=f"Slope: {row_s['slope_per_month']:+.4f}/mo\n(p = {row_s['p_value']:.4f}, q = {row_s['q_value_bh']:.4f}, R² = {row_s['r2']:.2f})")
    axes[1].set_xlabel("Donor Age (Months)")
    axes[1].set_ylabel("Sirt6 Expression (CP10K)")
    axes[1].set_title("Mouse Lung Fibroblasts: Sirt6")
    axes[1].legend(frameon=True)
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    fig_path = figures_dir / "fibroblast-candidate-markers.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[Markers Audit] Saved publication figure to {fig_path}")


if __name__ == "__main__":
    main()
