#!/usr/bin/env python3

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RESULTS = Path(__file__).resolve().parent.parent / "results"


def fit_age(part, extra=None):
    predictors = ["age", "mean_log_umi"] + (extra or [])
    X = sm.add_constant(part[predictors].astype(float), has_constant="add")
    fit = sm.OLS(part["profile_hellinger"].astype(float), X).fit()
    return fit


def permutation_rows(donors, analysis):
    rows = []
    for cell_type, part in donors.groupby("cell_type"):
        part = part.copy().reset_index(drop=True)
        if len(part) < 6 or len(part) > 7:
            continue
        observed = fit_age(part)
        observed_t = observed.tvalues["age"]
        permuted_t = []
        for indices in itertools.permutations(part.index):
            shuffled = part.copy()
            shuffled["age"] = part["age"].to_numpy()[list(indices)]
            permuted_t.append(fit_age(shuffled).tvalues["age"])
        permutation_p = float(np.mean(np.abs(permuted_t) >= abs(observed_t) - 1e-12))
        leave_one_out = []
        for donor in part["donor"]:
            leave_one_out.append(fit_age(part[part["donor"] != donor]).params["age"])
        rows.append({
            "analysis": analysis, "cell_type": cell_type, "n_fish": len(part),
            "slope_per_month": observed.params["age"], "parametric_p": observed.pvalues["age"],
            "exact_permutation_p": permutation_p,
            "leave_one_out_min": min(leave_one_out), "leave_one_out_max": max(leave_one_out),
        })
    result = pd.DataFrame(rows)
    result["permutation_q_bh"] = multipletests(result["exact_permutation_p"], method="fdr_bh")[1]
    return result


def main():
    fine = pd.read_csv(RESULTS / "extended-zebrafish-vector-donors-fine.csv")
    permutation_rows(fine, "fine_5_to_22_months").to_csv(
        RESULTS / "extended-zebrafish-vector-permutation.csv", index=False
    )
    broad = pd.read_csv(RESULTS / "extended-zebrafish-vector-donors-broad.csv")
    broad["incomplete_spermatogenesis"] = (
        broad["donor"].isin(["mo27_fish1", "mo27_fish2", "mo27_fish4"])
    ).astype(float)
    rows = []
    for cell_type, part in broad.groupby("cell_type"):
        fit = fit_age(part, ["incomplete_spermatogenesis"])
        rows.append({
            "cell_type": cell_type, "n_fish": len(part),
            "age_slope_per_month": fit.params["age"], "age_p": fit.pvalues["age"],
            "incomplete_effect": fit.params["incomplete_spermatogenesis"],
            "incomplete_p": fit.pvalues["incomplete_spermatogenesis"],
        })
    result = pd.DataFrame(rows)
    result["age_q_bh"] = multipletests(result["age_p"], method="fdr_bh")[1]
    result.to_csv(RESULTS / "extended-zebrafish-vector-spermatogenesis-sensitivity.csv", index=False)
    print(pd.read_csv(RESULTS / "extended-zebrafish-vector-permutation.csv").to_string(index=False))
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
