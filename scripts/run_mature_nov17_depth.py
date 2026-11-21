#!/usr/bin/env python3
"""Refit the November 17 pooled depth models using donors aged 30-plus."""

from datetime import datetime, timezone
from pathlib import Path
import json
import time

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
OUTPUT = RESULTS / "followup-mature-nov17-depth.csv"
RUNTIME = RESULTS / "followup-mature-nov17-runtime.csv"
COHORTS = {
    "Gao_heart_GSE290367": "gao",
    "Yazar_blood_GSE196830": "yazar",
}
OUTCOMES = [
    "mean_log_umi_ge1500",
    "mean_log_umi_all",
    "log_n_cells_all",
    "frac_cells_ge1500_EXPLORATORY",
]


def utc_now():
    return datetime.now(timezone.utc).isoformat()


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
    residual_df = float(fit.df_resid)
    se = float(fit.bse["age_per_decade"])
    mde80 = float((stats.t.ppf(0.975, residual_df) + stats.t.ppf(0.80, residual_df)) * se)
    return {
        "slope_per_decade": float(fit.params["age_per_decade"]),
        "standard_error": se,
        "ci_low": float(ci.iloc[0]),
        "ci_high": float(ci.iloc[1]),
        "p_value": float(fit.pvalues["age_per_decade"]),
        "n_donors": int(len(work)),
        "r_squared": float(fit.rsquared),
        "residual_df": residual_df,
        "design_rank": int(fit.model.rank),
        "mde80_per_decade": mde80,
        "mde80_outcome_sd_per_decade": float(mde80 / work[outcome].std(ddof=1)),
        "min_age": float(work["age_years"].min()),
        "max_age": float(work["age_years"].max()),
        "n_batches": int(work["study_or_pool"].nunique()),
    }


def main():
    started = time.monotonic()
    start_utc = utc_now()
    rows = []
    summary = {}
    for dataset, short_name in COHORTS.items():
        donors = pd.read_csv(RESULTS / f"followup-between-subject-{short_name}-donors.csv")
        pooled = donors["is_pooled"].astype(str).str.lower().eq("true")
        full = donors.loc[pooled].copy()
        if full["donor_id"].nunique() != len(full):
            raise ValueError(f"Expected one pooled row per donor for {dataset}")
        frame = full.loc[full["age_years"] >= 30].copy()
        if frame["age_years"].min() < 30:
            raise ValueError(f"Age restriction failed for {dataset}")
        frame["log_n_cells_all"] = np.log(frame["n_cells_all"])
        frame["frac_cells_ge1500_EXPLORATORY"] = frame["n_cells_ge1500"] / frame["n_cells_all"]
        summary[dataset] = {
            "n_original": int(len(full)),
            "n_mature": int(len(frame)),
            "n_lost": int(len(full) - len(frame)),
            "age_min": float(frame["age_years"].min()),
            "age_max": float(frame["age_years"].max()),
        }
        for outcome in OUTCOMES:
            for model, adjusted in [("A_unadjusted", False), ("B_adjusted", True)]:
                rows.append({
                    "dataset": dataset,
                    "outcome": outcome,
                    "model": model,
                    "n_original_donors": int(len(full)),
                    "n_donors_lost_to_age_restriction": int(len(full) - len(frame)),
                    **fit_model(frame, outcome, adjusted),
                    "method": "statsmodels OLS, classical SE; MDE80 two-sided alpha 0.05",
                })
    columns = [
        "dataset", "outcome", "model", "slope_per_decade", "standard_error",
        "ci_low", "ci_high", "p_value", "n_donors", "n_original_donors",
        "n_donors_lost_to_age_restriction", "r_squared", "residual_df",
        "design_rank", "mde80_per_decade", "mde80_outcome_sd_per_decade",
        "min_age", "max_age", "n_batches", "method",
    ]
    pd.DataFrame(rows)[columns].to_csv(OUTPUT, index=False)
    runtime_row = pd.DataFrame([{
        "pipeline": "mature_nov17_depth",
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "duration_seconds": time.monotonic() - started,
        "model_used": "kimi-code/k3.1-pro",
        "status": "completed",
    }])
    if RUNTIME.exists():
        runtime_row = pd.concat([pd.read_csv(RUNTIME), runtime_row], ignore_index=True)
    runtime_row.to_csv(RUNTIME, index=False)
    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(rows)[columns].to_string(index=False))
    print(f"Wrote {len(rows)} mature donor-level fits to {OUTPUT}")


if __name__ == "__main__":
    main()
