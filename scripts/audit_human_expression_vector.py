#!/usr/bin/env python3
"""Audit design, replication, and UMI-depth dependence of human vector associations."""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def fit(frame, outcome, predictors):
    X = sm.add_constant(frame[predictors].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        predictors = [name for name in predictors if name != "sex_male"]
        X = sm.add_constant(frame[predictors].astype(float), has_constant="add")
    return sm.OLS(frame[outcome].astype(float), X).fit()


def main():
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    donors = donors[donors["species"] == "Human"].copy()
    tests = pd.read_csv(RESULTS / "h3-cross-species-vector-tests.csv")
    tests = tests[
        (tests["species"] == "Human")
        & (tests["metric"] == "profile_hellinger")
    ]
    rows = []
    for (cohort, cell_type), frame in donors.groupby(["cohort", "cell_type"]):
        if len(frame) < 6:
            continue
        depth_age = fit(frame, "mean_log_umi", ["age", "sex_male"])
        distance_depth = fit(frame, "profile_hellinger", ["age", "sex_male", "mean_log_umi"])
        unit_tests = tests[(tests["cohort"] == cohort) & (tests["cell_type"] == cell_type)]
        primary = unit_tests[unit_tests["model"] == "unadjusted"].iloc[0]
        adjusted = unit_tests[unit_tests["model"] == "depth_adjusted"].iloc[0]
        rows.append({
            "cohort": cohort,
            "cell_type": cell_type,
            "n_donors": len(frame),
            "age_sd": frame["age"].std(),
            "age_min": frame["age"].min(),
            "age_max": frame["age"].max(),
            "median_cells": frame["n_cells"].median(),
            "mean_log_umi": frame["mean_log_umi"].mean(),
            "depth_age_slope_per_year": depth_age.params["age"],
            "depth_age_pvalue": depth_age.pvalues["age"],
            "distance_depth_slope": distance_depth.params["mean_log_umi"],
            "distance_depth_pvalue": distance_depth.pvalues["mean_log_umi"],
            "primary_slope": primary["slope_per_year"],
            "primary_pvalue": primary["pvalue"],
            "primary_q": primary["q_bh"],
            "depth_adjusted_slope": adjusted["slope_per_year"],
            "depth_adjusted_pvalue": adjusted["pvalue"],
            "depth_adjusted_q": adjusted["q_bh"],
            "primary_fdr_hit": primary["q_bh"] < 0.05,
            "depth_adjusted_fdr_hit": adjusted["q_bh"] < 0.05,
        })
    audit = pd.DataFrame(rows)
    audit.to_csv(RESULTS / "h3-human-vector-design-audit.csv", index=False)

    replication_rows = []
    labels = [
        "GABAergic neuron", "VIP GABAergic cortical interneuron",
        "sst GABAergic cortical interneuron", "pvalb GABAergic cortical interneuron",
    ]
    brain = donors[donors["cell_type"].isin(labels) & donors["cohort"].isin(["MSSM_DLPFC", "HBCC_DLPFC"])].copy()
    for cell_type, frame in brain.groupby("cell_type"):
        frame["mssm"] = (frame["cohort"] == "MSSM_DLPFC").astype(float)
        frame["age_mssm"] = frame["age"] * frame["mssm"]
        X = sm.add_constant(frame[["age", "sex_male", "mssm", "age_mssm"]], has_constant="add")
        fit_interaction = sm.OLS(frame["profile_hellinger"], X).fit(
            cov_type="cluster", cov_kwds={"groups": frame["cohort"] + "|" + frame["donor"]}, use_t=True
        )
        for cohort, part in frame.groupby("cohort"):
            age_fit = fit(part, "profile_hellinger", ["age", "sex_male"])
            replication_rows.append({
                "cell_type": cell_type, "cohort": cohort, "n_donors": len(part),
                "slope_per_year": age_fit.params["age"], "pvalue": age_fit.pvalues["age"],
                "cohort_by_age_interaction": fit_interaction.params["age_mssm"],
                "interaction_pvalue": fit_interaction.pvalues["age_mssm"],
            })
    replication = pd.DataFrame(replication_rows)
    replication.to_csv(RESULTS / "h3-human-brain-cohort-replication.csv", index=False)
    print(audit[audit["primary_fdr_hit"]].to_string(index=False), flush=True)
    print(replication.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
