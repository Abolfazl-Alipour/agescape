#!/usr/bin/env python3
"""Regenerate the November 17 pooled donor-level depth models."""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
OUTPUT = RESULTS / "csv-tests-2026-11-17" / "depth-vs-age" / "depth-vs-age-results.csv"
COHORTS = {
    "Gao_heart_GSE290367": "gao",
    "Yazar_blood_GSE196830": "yazar",
}
OUTCOMES = ["mean_log_umi_ge1500", "mean_log_umi_all", "log_n_cells_all", "frac_cells_ge1500_EXPLORATORY"]


def fit_model(frame, outcome, adjusted):
    work = frame.dropna(subset=[outcome, "age_per_decade", "sex_male", "study_or_pool"]).copy()
    design = pd.DataFrame({
        "intercept": np.ones(len(work)),
        "age_per_decade": work["age_per_decade"].to_numpy(dtype=float),
    }, index=work.index)
    if adjusted:
        design["sex_male"] = work["sex_male"].to_numpy(dtype=float)
        batches = pd.get_dummies(work["study_or_pool"].astype(str), drop_first=True, dtype=float)
        design = pd.concat([design, batches], axis=1)
    fit = sm.OLS(work[outcome].to_numpy(dtype=float), design.astype(float)).fit()
    ci = fit.conf_int().loc["age_per_decade"]
    return {
        "slope_per_decade": fit.params["age_per_decade"],
        "ci_low": ci.iloc[0],
        "ci_high": ci.iloc[1],
        "p_value": fit.pvalues["age_per_decade"],
        "n_donors": len(work),
        "r_squared": fit.rsquared,
    }


def main():
    rows = []
    for dataset, short_name in COHORTS.items():
        donors = pd.read_csv(RESULTS / f"followup-between-subject-{short_name}-donors.csv")
        pooled = donors["is_pooled"].astype(str).str.lower().eq("true")
        frame = donors.loc[pooled].copy()
        if frame["donor_id"].nunique() != len(frame):
            raise ValueError(f"Expected one pooled row per donor for {dataset}")
        frame["log_n_cells_all"] = np.log(frame["n_cells_all"])
        frame["frac_cells_ge1500_EXPLORATORY"] = frame["n_cells_ge1500"] / frame["n_cells_all"]
        for outcome in OUTCOMES:
            for model, adjusted in [("A_unadjusted", False), ("B_adjusted", True)]:
                rows.append({
                    "dataset": dataset,
                    "outcome": outcome,
                    "model": model,
                    **fit_model(frame, outcome, adjusted),
                    "method": "statsmodels OLS, classical SE",
                })
    columns = ["dataset", "outcome", "model", "slope_per_decade", "ci_low", "ci_high", "p_value", "n_donors", "r_squared", "method"]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows)[columns].to_csv(OUTPUT, index=False)
    print(f"Wrote {len(rows)} donor-level fits to {OUTPUT}")


if __name__ == "__main__":
    main()
