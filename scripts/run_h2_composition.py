#!/usr/bin/env python3
"""
scripts/run_h2_composition.py

Tests Hypothesis 2: Cell-Type Compositional Remodeling in Mouse and Human Lung.
Uses Aitchison Centred Log-Ratio (CLR) transformation on donor-level cell counts
to avoid closed-simplex negative correlative bias.

Outputs:
  - results/h2-composition-mouse-lung.csv
  - results/h2-composition-human-lung.csv
  - results/h2-composition-cross-species.csv
  - figures/h2-composition-shifts-mouse-vs-human.png
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr, binomtest
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Homology mapping: mouse cell_ontology_class -> human cell_type
HOMOLOGY_MAP = {
    "classical monocyte": "classical monocyte",
    "non-classical monocyte": "non-classical monocyte",
    "alveolar macrophage": "alveolar macrophage",
    "lung macrophage": "lung macrophage",
    "B cell": "B cell",
    "CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
    "CD4-positive, alpha-beta T cell": "CD4-positive, alpha-beta T cell",
    "NK cell": "natural killer cell",
    "type II pneumocyte": "pulmonary alveolar type 2 cell",
    "vein endothelial cell": "vein endothelial cell",
    "fibroblast of lung": "alveolar type 1 fibroblast cell",
}


def compute_clr_table(count_matrix: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    """
    Computes Centred Log-Ratio (CLR) for each row (donor):
    CLR(p_ik) = ln(p_ik) - mean(ln(p_ij))
    Strict mathematical invariant: row sums of CLR must equal 0.
    """
    counts = count_matrix.values.astype(float) + pseudocount
    totals = counts.sum(axis=1, keepdims=True)
    props = counts / totals
    log_props = np.log(props)
    geom_means = log_props.mean(axis=1, keepdims=True)
    clr = log_props - geom_means
    clr_df = pd.DataFrame(clr, index=count_matrix.index, columns=count_matrix.columns)
    # Verify simplex closure invariant
    np.testing.assert_allclose(clr_df.sum(axis=1).values, 0.0, atol=1e-10)
    return clr_df


def run_composition_regression(clr_df: pd.DataFrame, meta_df: pd.DataFrame, age_col: str = "age", sex_col: str = "sex") -> pd.DataFrame:
    """
    Fits OLS for each cell type: CLR ~ Age + Sex.
    """
    rows = []
    donors = clr_df.index.intersection(meta_df.index)
    clr_sub = clr_df.loc[donors]
    meta_sub = meta_df.loc[donors]

    for ct in clr_sub.columns:
        y = clr_sub[ct].values
        X_df = pd.DataFrame({
            "const": 1.0,
            "age": meta_sub[age_col].values,
            "sex": meta_sub[sex_col].values
        })
        model = sm.OLS(y, X_df).fit()
        slope = model.params["age"]
        se = model.bse["age"]
        pval = model.pvalues["age"]
        ci_lo, ci_hi = model.conf_int().loc["age"]

        rows.append({
            "cell_type": ct,
            "n_donors": len(donors),
            "slope_per_year_or_month": slope,
            "se": se,
            "pvalue": pval,
            "ci_lo": ci_lo,
            "ci_hi": ci_hi,
            "r2": model.rsquared
        })

    res = pd.DataFrame(rows)
    # Benjamini-Hochberg FDR
    from statsmodels.stats.multitest import multipletests
    _, qvals, _, _ = multipletests(res["pvalue"], method="fdr_bh")
    res["q_bh"] = qvals
    res = res.sort_values("pvalue").reset_index(drop=True)
    return res


def main():
    print("[H2 Composition] Starting Phase 1: Cell-Type Compositional Remodeling...")
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    # -------------------------------------------------------------
    # 1. Mouse Arm: Tabula Muris Senis Lung Droplet
    # -------------------------------------------------------------
    print("[H2 Composition] Processing Mouse Lung from data/raw/tms_droplet_raw.h5ad ...")
    import anndata as ad
    adata = ad.read_h5ad(REPO_ROOT / "data" / "raw" / "tms_droplet_raw.h5ad", backed="r")
    obs = adata.obs
    mouse_lung = obs[obs["tissue"] == "Lung"].copy()
    
    # Map age strings to numeric months
    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    mouse_lung["age_months"] = mouse_lung["age"].map(age_map)
    mouse_lung["sex_bin"] = (mouse_lung["sex"].str.lower() == "male").astype(float)

    # Filter cell types with at least 50 total cells
    ct_counts = mouse_lung["cell_ontology_class"].value_counts()
    valid_cts = ct_counts[ct_counts >= 50].index.tolist()

    # Pivot count matrix: mouse.id x cell_type
    mouse_counts = (
        mouse_lung[mouse_lung["cell_ontology_class"].isin(valid_cts)]
        .groupby(["mouse.id", "cell_ontology_class"], observed=True)
        .size()
        .unstack(fill_value=0)
    )

    # Mouse metadata
    mouse_meta = (
        mouse_lung.groupby("mouse.id", observed=True)
        .agg({"age_months": "first", "sex_bin": "first"})
        .rename(columns={"age_months": "age", "sex_bin": "sex"})
    )

    # Compute CLR
    mouse_clr = compute_clr_table(mouse_counts)
    mouse_res = run_composition_regression(mouse_clr, mouse_meta)
    mouse_res.to_csv(results_dir / "h2-composition-mouse-lung.csv", index=False)
    print(f"[H2 Composition] Saved mouse results ({len(mouse_res)} cell types) to results/h2-composition-mouse-lung.csv")
    print(mouse_res.head(6)[["cell_type", "slope_per_year_or_month", "pvalue", "q_bh"]])

    # -------------------------------------------------------------
    # 2. Human Arm: CELLxGENE Lung_Integrated (137 donors)
    # -------------------------------------------------------------
    print("\n[H2 Composition] Processing Human Lung from results/m9-donor-scores.csv ...")
    m9_df = pd.read_csv(results_dir / "m9-donor-scores.csv")
    human_lung = m9_df[m9_df["cohort"] == "Lung_Integrated"].copy()

    # Filter cell types present across at least 25 donors
    h_ct_counts = human_lung["cell_type"].value_counts()
    valid_h_cts = h_ct_counts[h_ct_counts >= 25].index.tolist()

    # Pivot count matrix: donor_id x cell_type
    human_counts = (
        human_lung[human_lung["cell_type"].isin(valid_h_cts)]
        .set_index(["donor_id", "cell_type"])["n_cells"]
        .unstack(fill_value=0)
    )

    human_meta = (
        human_lung.groupby("donor_id")
        .agg({"age": "first", "sex": "first"})
    )
    human_meta["sex"] = (human_meta["sex"].str.lower() == "male").astype(float)

    # Compute CLR
    human_clr = compute_clr_table(human_counts)
    human_res = run_composition_regression(human_clr, human_meta)
    human_res.to_csv(results_dir / "h2-composition-human-lung.csv", index=False)
    print(f"[H2 Composition] Saved human results ({len(human_res)} cell types) to results/h2-composition-human-lung.csv")
    print(human_res.head(6)[["cell_type", "slope_per_year_or_month", "pvalue", "q_bh"]])

    # -------------------------------------------------------------
    # 3. Cross-Species Homology & Concordance
    # -------------------------------------------------------------
    print("\n[H2 Composition] Evaluating Cross-Species Concordance across Homologous Cell Types...")
    comp_rows = []
    mouse_dict = mouse_res.set_index("cell_type").to_dict("index")
    human_dict = human_res.set_index("cell_type").to_dict("index")

    for m_ct, h_ct in HOMOLOGY_MAP.items():
        if m_ct in mouse_dict and h_ct in human_dict:
            m_stat = mouse_dict[m_ct]
            h_stat = human_dict[h_ct]
            same_sign = (m_stat["slope_per_year_or_month"] * h_stat["slope_per_year_or_month"]) > 0

            comp_rows.append({
                "homologous_pair": f"{m_ct} <-> {h_ct}",
                "mouse_cell_type": m_ct,
                "human_cell_type": h_ct,
                "mouse_slope_per_month": m_stat["slope_per_year_or_month"],
                "mouse_pval": m_stat["pvalue"],
                "mouse_q_bh": m_stat["q_bh"],
                "human_slope_per_year": h_stat["slope_per_year_or_month"],
                "human_pval": h_stat["pvalue"],
                "human_q_bh": h_stat["q_bh"],
                "sign_concordance": same_sign
            })

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(results_dir / "h2-composition-cross-species.csv", index=False)
    print(f"[H2 Composition] Saved cross-species concordance ({len(comp_df)} pairs) to results/h2-composition-cross-species.csv")

    # Statistical concordance tests
    rho, rho_p = spearmanr(comp_df["mouse_slope_per_month"], comp_df["human_slope_per_year"])
    n_concordant = comp_df["sign_concordance"].sum()
    n_total = len(comp_df)
    binom_res = binomtest(n_concordant, n_total, p=0.5, alternative="greater")

    print("\n--- Cross-Species Compositional Remodeling Verdict ---")
    print(f"Homologous pairs tested: {n_total}")
    print(f"Sign concordant pairs: {n_concordant}/{n_total} ({n_concordant/n_total:.1%})")
    print(f"Binomial sign test p-value: {binom_res.pvalue:.4f}")
    print(f"Spearman rank correlation rho: {rho:.3f} (p = {rho_p:.4f})")

    # -------------------------------------------------------------
    # 4. Generate Publication-Quality Figure
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Forest plot of mouse CLR slopes
    m_sorted = mouse_res.sort_values("slope_per_year_or_month")
    y_pos = np.arange(len(m_sorted))
    axes[0].barh(y_pos, m_sorted["slope_per_year_or_month"], color=["crimson" if q < 0.05 else "gray" for q in m_sorted["q_bh"]], alpha=0.8)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(m_sorted["cell_type"], fontsize=9)
    axes[0].axvline(0, color="black", linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("Mouse Lung CLR Slope (per month)")
    axes[0].set_title(f"A. Mouse Lung Compositional Shifts (N={len(mouse_meta)} mice)\nRed: FDR q < 0.05")

    # Panel B: Cross-species concordance scatter
    for _, row in comp_df.iterrows():
        color = "crimson" if row["sign_concordance"] else "navy"
        axes[1].scatter(row["mouse_slope_per_month"], row["human_slope_per_year"], color=color, s=70, zorder=3)
        # Label points
        label = row["mouse_cell_type"][:18]
        axes[1].annotate(label, (row["mouse_slope_per_month"], row["human_slope_per_year"]), fontsize=8, xytext=(4, 2), textcoords="offset points")

    axes[1].axhline(0, color="black", linestyle="--", linewidth=0.8)
    axes[1].axvline(0, color="black", linestyle="--", linewidth=0.8)
    axes[1].set_xlabel("Mouse CLR Slope (per month)")
    axes[1].set_ylabel("Human CLR Slope (per year)")
    axes[1].set_title(f"B. Cross-Species Concordance (Spearman rho = {rho:.2f}, p = {rho_p:.3f})\nSign Concordance: {n_concordant}/{n_total} ({n_concordant/n_total:.0%})")

    plt.tight_layout()
    fig_path = figures_dir / "h2-composition-shifts-mouse-vs-human.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[H2 Composition] Saved publication figure to {fig_path}")


if __name__ == "__main__":
    main()
