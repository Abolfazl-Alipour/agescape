#!/usr/bin/env python3
"""Estimate mouse and human lung cell-type composition changes with donor-level counts."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import binomtest, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


def compute_clr_table(count_matrix: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    adjusted = count_matrix.astype(float) + pseudocount
    props = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_props = np.log(props)
    clr = log_props.sub(log_props.mean(axis=1), axis=0)
    assert np.allclose(clr.sum(axis=1), 0.0, atol=1e-5)
    return clr


def run_composition_regression(clr_df: pd.DataFrame, meta_df: pd.DataFrame, observed_counts=None) -> pd.DataFrame:
    records = []
    common_donors = clr_df.index.intersection(meta_df.index)
    clr_df = clr_df.loc[common_donors]
    meta_df = meta_df.loc[common_donors]
    for cell_type in clr_df.columns:
        if observed_counts is None:
            donors = common_donors
        else:
            donors = common_donors[observed_counts.loc[common_donors, cell_type] > 0]
        if len(donors) < 4:
            continue
        y = clr_df.loc[donors, cell_type]
        X = sm.add_constant(meta_df.loc[donors, ["age", "sex"]])
        model = sm.OLS(y, X).fit()
        records.append({
            "cell_type": cell_type,
            "slope_per_year_or_month": model.params["age"],
            "se": model.bse["age"],
            "tvalues": model.tvalues["age"],
            "pvalue": model.pvalues["age"],
            "r2": model.rsquared,
            "sex_covar_p": model.pvalues["sex"],
            "n_donors": len(donors),
        })
    res = pd.DataFrame(records)
    res["q_bh"] = multipletests(res["pvalue"], method="fdr_bh")[1]
    return res.sort_values("pvalue").reset_index(drop=True)


def presence_analysis(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell_type in counts.columns:
        present = (counts[cell_type] > 0).astype(int)
        if present.nunique() < 2:
            r, p = np.nan, np.nan
        else:
            r, p = pearsonr(meta.loc[counts.index, "age"], present)
        rows.append({
            "cell_type": cell_type,
            "n_present": int(present.sum()),
            "n_absent": int((1 - present).sum()),
            "presence_fraction": present.mean(),
            "presence_age_r": r,
            "presence_age_p": p,
            "n_donors": len(present),
        })
    result = pd.DataFrame(rows)
    valid = result["presence_age_p"].notna()
    result["q_bh"] = np.nan
    result.loc[valid, "q_bh"] = multipletests(result.loc[valid, "presence_age_p"], method="fdr_bh")[1]
    return result.sort_values("presence_age_p", na_position="last").reset_index(drop=True)


def load_human_counts(results_dir: Path):
    source = results_dir / "m9-donor-cell-counts-prethreshold.csv"
    if not source.exists():
        raise FileNotFoundError(f"Missing complete pre-threshold counts: {source}")
    data = pd.read_csv(source)
    human = data[(data["cohort"] == "Lung_Integrated") & (data["age"] >= 20)].copy()
    eligible = human["cell_type"].value_counts()
    valid_cell_types = eligible[eligible >= 20].index.tolist()
    human = human[human["cell_type"].isin(valid_cell_types)]
    counts = (human.set_index(["donor_id", "cell_type"])["n_cells"]
              .unstack(fill_value=0).sort_index())
    meta = human.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"))
    meta["sex"] = (meta["sex"].astype(str).str.lower() == "male").astype(float)
    return counts, meta


def main():
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    import anndata as ad
    adata = ad.read_h5ad(REPO_ROOT / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    obs = adata.obs
    mouse_lung = obs[(obs["tissue"] == "Lung") & (~obs["predicted_doublet"])].copy()
    age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
    mouse_lung["age_months"] = mouse_lung["age"].astype(str).map(age_map).astype(float)
    mouse_lung["sex_bin"] = (mouse_lung["sex"].astype(str).str.lower() == "male").astype(float)
    ct_counts = mouse_lung["cell_ontology_class"].value_counts()
    valid_cts = ct_counts[ct_counts >= 50].index.tolist()
    mouse_counts = (mouse_lung[mouse_lung["cell_ontology_class"].isin(valid_cts)]
                    .groupby(["mouse.id", "cell_ontology_class"], observed=True)
                    .size().unstack(fill_value=0))
    mouse_meta = (mouse_lung.groupby("mouse.id", observed=True)
                  .agg(age_months=("age_months", "first"), sex_bin=("sex_bin", "first"))
                  .rename(columns={"age_months": "age", "sex_bin": "sex"}))
    mouse_clr = compute_clr_table(mouse_counts)
    mouse_res = run_composition_regression(mouse_clr, mouse_meta)
    mouse_res.to_csv(results_dir / "h2-composition-mouse-lung.csv", index=False)

    human_counts, human_meta = load_human_counts(results_dir)
    human_clr = compute_clr_table(human_counts)
    human_res = run_composition_regression(human_clr, human_meta)
    human_res.to_csv(results_dir / "h2-composition-human-lung.csv", index=False)
    observed_res = run_composition_regression(human_clr, human_meta, observed_counts=human_counts)
    observed_res.to_csv(results_dir / "h2-composition-human-lung-observed-only.csv", index=False)
    presence_res = presence_analysis(human_counts, human_meta)
    presence_res.to_csv(results_dir / "h2-composition-human-lung-presence.csv", index=False)

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
        "club cell": "club cell",
    }
    mouse_dict = mouse_res.set_index("cell_type")
    human_dict = human_res.set_index("cell_type")
    rows = []
    for mouse_ct, human_ct in homology_map.items():
        if mouse_ct not in mouse_dict.index or human_ct not in human_dict.index:
            continue
        m = mouse_dict.loc[mouse_ct]
        h = human_dict.loc[human_ct]
        rows.append({
            "mouse_cell_type": mouse_ct, "human_cell_type": human_ct,
            "mouse_slope_per_month": m["slope_per_year_or_month"], "mouse_pval": m["pvalue"],
            "mouse_q_bh": m["q_bh"], "human_slope_per_year": h["slope_per_year_or_month"],
            "human_pval": h["pvalue"], "human_q_bh": h["q_bh"],
            "sign_concordance": bool(m["slope_per_year_or_month"] * h["slope_per_year_or_month"] > 0),
        })
    comparison = pd.DataFrame(rows)
    comparison.to_csv(results_dir / "h2-composition-cross-species.csv", index=False)
    comparison.to_csv(results_dir / "h2-composition-cross-species-homologs.csv", index=False)
    rho, rho_p = spearmanr(comparison["mouse_slope_per_month"], comparison["human_slope_per_year"])
    concordant = int(comparison["sign_concordance"].sum())
    total = len(comparison)
    test_row = pd.DataFrame([{
        "n_pairs": total, "n_sign_concordant": concordant,
        "sign_concordance_fraction": concordant / total,
        "binomial_two_sided_p": binomtest(concordant, total, p=0.5, alternative="two-sided").pvalue,
        "binomial_less_p": binomtest(concordant, total, p=0.5, alternative="less").pvalue,
        "spearman_rho": rho, "spearman_p": rho_p,
    }])
    test_row.to_csv(results_dir / "h2-composition-cross-species-test.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, result, title, color in (
        (axes[0], mouse_res, f"A. Mouse lung (N = {len(mouse_meta)} mice)", "crimson"),
        (axes[1], human_res, f"B. Human adult lung (N = {len(human_meta)} donors)", "royalblue"),
    ):
        ordered = result.sort_values("slope_per_year_or_month")
        colors = [color if q < 0.05 else "darkgray" for q in ordered["q_bh"]]
        ax.barh(ordered["cell_type"], ordered["slope_per_year_or_month"], color=colors, edgecolor="black", alpha=0.85)
        ax.axvline(0, color="black", linestyle="--", linewidth=1.0)
        ax.set_xlabel("CLR slope per month" if ax is axes[0] else "CLR slope per year")
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures_dir / "h2-composition-shifts-mouse-vs-human.png", dpi=200)
    plt.close(fig)

    print(f"Mouse: {len(mouse_res)} cell types, {len(mouse_meta)} donors")
    print(f"Human primary complete counts: {len(human_res)} cell types, {len(human_meta)} donors")
    print(human_res.head(10).to_string(index=False))
    print("Presence versus age:")
    print(presence_res.head(10).to_string(index=False))
    print("Observed-record-only sensitivity:")
    print(observed_res.head(10).to_string(index=False))
    print(test_row.to_string(index=False))


if __name__ == "__main__":
    main()
