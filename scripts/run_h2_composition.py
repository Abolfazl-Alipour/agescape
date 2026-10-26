#!/usr/bin/env python3
"""
scripts/run_h2_composition.py

Tests Hypothesis 2: Cell-Type Compositional Remodeling in Aging across Species.
Quantifies shifts in cellular proportions across the adult lifespan using
Aitchison Centered Log-Ratio (CLR) transformation:
  1. Mouse Lung Droplet (TMS, 23 mice across 1m, 3m, 18m, 21m, 24m, 30m).
  2. Human Lung Integrated (137 donors, ages 2-81, with adult stratum >= 20y, N=117).
  3. Cross-Species Homologous Cell-Type Concordance Testing (two-sided binomial & Spearman rho).

Outputs:
  - results/h2-composition-mouse-lung.csv
  - results/h2-composition-human-lung.csv
  - results/h2-composition-cross-species-homologs.csv
  - figures/h2-composition-shifts-mouse-vs-human.png
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr, binomtest
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


def compute_clr_table(count_matrix: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    """
    Computes Centered Log-Ratio (CLR) transformation on a sample x cell_type matrix.
    Adds Laplace pseudocount to handle sparsity, normalizes to proportions,
    and divides by geometric mean:
      p_{ik} = (n_{ik} + c) / sum_j (n_{ij} + c)
      g(p_i) = exp(1/K sum_j ln(p_{ij}))
      CLR_{ik} = ln(p_{ik} / g(p_i))
    """
    counts = count_matrix.copy()
    props = (counts + pseudocount).div((counts + pseudocount).sum(axis=1), axis=0)
    log_props = np.log(props)
    geom_mean_log = log_props.mean(axis=1)
    clr = log_props.sub(geom_mean_log, axis=0)

    # Simplex sum invariant check: sum of CLR values across cell types per sample == 0
    clr_sums = clr.sum(axis=1)
    assert np.allclose(clr_sums, 0.0, atol=1e-5), "CLR simplex closure invariant failed!"
    return clr


def run_composition_regression(clr_df: pd.DataFrame, meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Regresses each cell-type's CLR abundance against Age + Sex across donors:
      CLR_ik ~ Age_i + Sex_i + Intercept
    """
    records = []
    common_donors = clr_df.index.intersection(meta_df.index)
    clr_df = clr_df.loc[common_donors]
    meta_df = meta_df.loc[common_donors]

    for ct in clr_df.columns:
        y = clr_df[ct]
        X = sm.add_constant(meta_df[["age", "sex"]])
        model = sm.OLS(y, X).fit()

        records.append({
            "cell_type": ct,
            "slope_per_year_or_month": model.params["age"],
            "se": model.bse["age"],
            "tvalues": model.tvalues["age"],
            "pvalue": model.pvalues["age"],
            "r2": model.rsquared,
            "sex_covar_p": model.pvalues["sex"],
            "n_donors": len(common_donors)
        })

    res = pd.DataFrame(records)
    res["q_bh"] = multipletests(res["pvalue"], method="fdr_bh")[1]
    return res.sort_values("pvalue")


def main():
    print("[H2 Composition] Starting Phase 1: Cell-Type Compositional Remodeling...")
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    # -------------------------------------------------------------
    # 1. Mouse Arm: Tabula Muris Senis Lung Droplet (QC-cleaned)
    # -------------------------------------------------------------
    print("[H2 Composition] Processing Mouse Lung from data/interim/m1-qc-cleaned.h5ad ...")
    import anndata as ad
    adata = ad.read_h5ad(REPO_ROOT / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    obs = adata.obs
    mouse_lung = obs[(obs["tissue"] == "Lung") & (~obs["predicted_doublet"])].copy()

    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    mouse_lung["age_months"] = mouse_lung["age"].astype(str).map(age_map).astype(float)
    mouse_lung["sex_bin"] = (mouse_lung["sex"].astype(str).str.lower() == "male").astype(float)

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

    mouse_meta = (
        mouse_lung.groupby("mouse.id", observed=True)
        .agg({"age_months": "first", "sex_bin": "first"})
        .rename(columns={"age_months": "age", "sex_bin": "sex"})
    )

    mouse_clr = compute_clr_table(mouse_counts)
    mouse_res = run_composition_regression(mouse_clr, mouse_meta)
    mouse_res.to_csv(results_dir / "h2-composition-mouse-lung.csv", index=False)
    print(f"[H2 Composition] Saved mouse results ({len(mouse_res)} cell types) to results/h2-composition-mouse-lung.csv")
    print(mouse_res.head(6)[["cell_type", "slope_per_year_or_month", "pvalue", "q_bh"]])

    # -------------------------------------------------------------
    # 2. Human Arm: CELLxGENE Lung_Integrated (Adult stratum >= 20y primary, full cohort audited)
    # -------------------------------------------------------------
    print("\n[H2 Composition] Processing Human Lung from results/m9-donor-scores.csv ...")
    m9_df = pd.read_csv(results_dir / "m9-donor-scores.csv")
    human_lung = m9_df[m9_df["cohort"] == "Lung_Integrated"].copy()

    # Primary analysis: Adult donors (age >= 20.0, N=117)
    human_adult = human_lung[human_lung["age"] >= 20.0].copy()
    h_ct_counts = human_adult["cell_type"].value_counts()
    valid_h_cts = h_ct_counts[h_ct_counts >= 20].index.tolist()

    human_counts = (
        human_adult[human_adult["cell_type"].isin(valid_h_cts)]
        .set_index(["donor_id", "cell_type"])["n_cells"]
        .unstack(fill_value=0)
    )

    human_meta = (
        human_adult.groupby("donor_id")
        .agg({"age": "first", "sex": "first"})
    )
    human_meta["sex"] = (human_meta["sex"].astype(str).str.lower() == "male").astype(float)

    human_clr = compute_clr_table(human_counts)
    human_res = run_composition_regression(human_clr, human_meta)
    human_res.to_csv(results_dir / "h2-composition-human-lung.csv", index=False)
    print(f"[H2 Composition] Saved human adult results ({len(human_res)} cell types, N={len(human_meta)} donors) to results/h2-composition-human-lung.csv")
    print(human_res.head(6)[["cell_type", "slope_per_year_or_month", "pvalue", "q_bh"]])

    # -------------------------------------------------------------
    # 3. Cross-Species Homologous Trajectory Comparison
    # -------------------------------------------------------------
    homology_map = {
        "fibroblast of lung": "fibroblast",
        "bronchial smooth muscle cell": "smooth muscle cell",
        "alveolar macrophage": "alveolar macrophage",
        "lung macrophage": "macrophage",
        "classical monocyte": "classical monocyte",
        "intermediate monocyte": "non-classical monocyte",
        "B cell": "B cell",
        "NK cell": "natural killer cell",
        "CD8-positive, alpha-beta T cell": "CD8-positive, alpha-beta T cell",
        "vein endothelial cell": "vein endothelial cell",
        "club cell": "club cell"
    }

    m_dict = mouse_res.set_index("cell_type")
    h_dict = human_res.set_index("cell_type")

    comp_rows = []
    for m_ct, h_ct in homology_map.items():
        if m_ct in m_dict.index and h_ct in h_dict.index:
            m_slope = m_dict.loc[m_ct, "slope_per_year_or_month"]
            m_pval = m_dict.loc[m_ct, "pvalue"]
            m_q = m_dict.loc[m_ct, "q_bh"]
            h_slope = h_dict.loc[h_ct, "slope_per_year_or_month"]
            h_pval = h_dict.loc[h_ct, "pvalue"]
            h_q = h_dict.loc[h_ct, "q_bh"]

            concordant = (m_slope * h_slope) > 0
            comp_rows.append({
                "mouse_cell_type": m_ct,
                "human_cell_type": h_ct,
                "mouse_slope_per_month": m_slope,
                "mouse_pval": m_pval,
                "mouse_q_bh": m_q,
                "human_slope_per_year": h_slope,
                "human_pval": h_pval,
                "human_q_bh": h_q,
                "sign_concordance": concordant
            })

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(results_dir / "h2-composition-cross-species-homologs.csv", index=False)
    print(f"\n[H2 Composition] Saved cross-species homologs comparison to results/h2-composition-cross-species-homologs.csv")

    # Statistical concordance tests
    rho, rho_p = spearmanr(comp_df["mouse_slope_per_month"], comp_df["human_slope_per_year"])
    n_concordant = int(comp_df["sign_concordance"].sum())
    n_total = len(comp_df)
    binom_twosided = binomtest(n_concordant, n_total, p=0.5, alternative="two-sided")
    binom_less = binomtest(n_concordant, n_total, p=0.5, alternative="less")

    print("\n--- Cross-Species Compositional Remodeling Verdict ---")
    print(f"Homologous pairs tested: {n_total}")
    print(f"Sign concordant pairs: {n_concordant}/{n_total} ({n_concordant/n_total:.1%})")
    print(f"Binomial two-sided test p-value: {binom_twosided.pvalue:.4f}")
    print(f"Binomial directional (divergence) p-value: {binom_less.pvalue:.4f}")
    print(f"Spearman rank correlation rho: {rho:.3f} (p = {rho_p:.4f})")

    # -------------------------------------------------------------
    # 4. Generate Publication-Quality Figure
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Panel A: Mouse CLR Slopes
    m_sorted = mouse_res.sort_values("slope_per_year_or_month")
    colors_m = ["crimson" if q < 0.05 else "darkgray" for q in m_sorted["q_bh"]]
    axes[0].barh(m_sorted["cell_type"], m_sorted["slope_per_year_or_month"], color=colors_m, edgecolor="black", alpha=0.85)
    axes[0].axvline(0, color="black", linestyle="--", linewidth=1.0)
    axes[0].set_xlabel("Mouse Aging Trajectory (CLR slope per month)")
    axes[0].set_title(f"A. Mouse Lung Compositional Remodeling\n(N = 23 mice; red = FDR q < 0.05)")
    axes[0].grid(axis="x", alpha=0.3)

    # Panel B: Human CLR Slopes (Adult Stratum)
    h_sorted = human_res.sort_values("slope_per_year_or_month")
    colors_h = ["royalblue" if q < 0.05 else "darkgray" for q in h_sorted["q_bh"]]
    axes[1].barh(h_sorted["cell_type"], h_sorted["slope_per_year_or_month"], color=colors_h, edgecolor="black", alpha=0.85)
    axes[1].axvline(0, color="black", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("Human Aging Trajectory (CLR slope per year)")
    axes[1].set_title(f"B. Human Adult Lung Remodeling\n(N = 117 adult donors; blue = FDR q < 0.05)")
    axes[1].grid(axis="x", alpha=0.3)

    plt.tight_layout()
    fig_path = figures_dir / "h2-composition-shifts-mouse-vs-human.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[H2 Composition] Saved publication figure to {fig_path}")


if __name__ == "__main__":
    main()
