#!/usr/bin/env python3
"""
scripts/generate_m10_stratified_summary.py

Generates the committed supplementary table for all 27 PsychAD DLPFC cell types,
reporting pooled covariate-adjusted and adult-stratum (>=20y) regressions alongside
within-family Benjamini-Hochberg FDR q-values.

Source: results/m10-confounder-audit.csv
Output: results/m10-celltype-stratified-summary.csv
"""
from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent

def main():
    audit_file = REPO / "results" / "m10-confounder-audit.csv"
    if not audit_file.exists():
        raise FileNotFoundError(f"Missing {audit_file}")

    df = pd.read_csv(audit_file)

    # Recompute / verify within-family Benjamini-Hochberg FDR
    q_covar = multipletests(df["p_covar_adj"], method="fdr_bh")[1]
    q_adult = multipletests(df["p_adult_ge20"], method="fdr_bh")[1]

    out_df = pd.DataFrame({
        "cell_type": df["cell_type"],
        "n_pooled": df["n_covar_adj"],
        "slope_pooled_adjusted": df["slope_covar_adj"],
        "se_pooled_adjusted": df["se_covar_adj"],
        "p_pooled_adjusted": df["p_covar_adj"],
        "q_pooled_adjusted_bh": q_covar,
        "n_adult": df["n_adult_ge20"],
        "slope_adult_stratum": df["slope_adult_ge20"],
        "se_adult_stratum": df["se_adult_ge20"],
        "p_adult_stratum": df["p_adult_ge20"],
        "q_adult_stratum_bh": q_adult
    })

    # Sort by p_pooled_adjusted
    out_df = out_df.sort_values("p_pooled_adjusted")

    out_path = REPO / "results" / "m10-celltype-stratified-summary.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved 27-cell-type stratified summary to {out_path}")
    print(out_df[["cell_type", "slope_pooled_adjusted", "p_pooled_adjusted", "q_pooled_adjusted_bh", "slope_adult_stratum", "p_adult_stratum", "q_adult_stratum_bh"]].head(10).to_string())

if __name__ == "__main__":
    main()
