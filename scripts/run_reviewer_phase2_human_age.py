#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_proliferation_moderator import core_composition

RESULTS = REPO / "results"


def fit_comparison(part, outcome, depth_adjusted):
    part = part.copy()
    spline = np.asarray(dmatrix("bs(age, df=3, degree=3, include_intercept=False)", part))
    linear_columns = [part["age"].to_numpy()]
    spline_columns = [spline[:, index] for index in range(spline.shape[1])]
    for columns in (linear_columns, spline_columns):
        if part["sex_male"].nunique() > 1:
            columns.append(part["sex_male"].to_numpy())
        if depth_adjusted:
            columns.append(part["mean_log_umi"].to_numpy())
    linear = sm.OLS(part[outcome], sm.add_constant(np.column_stack(linear_columns))).fit()
    nonlinear = sm.OLS(part[outcome], sm.add_constant(np.column_stack(spline_columns))).fit()
    fvalue, pvalue, _ = nonlinear.compare_f_test(linear)
    return linear, nonlinear, fvalue, pvalue


def analyze(data, outcome, analysis, depth_adjusted):
    rows = []
    for (cohort, tissue, cell_type), part in data.groupby(["cohort", "tissue", "cell_type"]):
        if len(part) < 8:
            continue
        linear, nonlinear, fvalue, pvalue = fit_comparison(part, outcome, depth_adjusted)
        rows.append({
            "analysis": analysis, "cohort": cohort, "tissue": tissue, "cell_type": cell_type,
            "n_donors": len(part), "age_min": part["age"].min(), "age_max": part["age"].max(),
            "linear_slope": linear.params[1], "linear_p": linear.pvalues[1],
            "nonlinear_f": fvalue, "nonlinear_p": pvalue,
        })
    result = pd.DataFrame(rows)
    result["q_linear_bh"] = np.nan
    result["q_nonlinear_bh"] = np.nan
    for _, indices in result.groupby("cohort").groups.items():
        result.loc[indices, "q_linear_bh"] = multipletests(result.loc[indices, "linear_p"], method="fdr_bh")[1]
        result.loc[indices, "q_nonlinear_bh"] = multipletests(result.loc[indices, "nonlinear_p"], method="fdr_bh")[1]
    return result


def main():
    vector = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    vector = vector[vector["species"] == "Human"].copy()
    vector_results = analyze(vector, "profile_hellinger", "expression_vector_depth_adjusted", True)
    vector_results.to_csv(RESULTS / "reviewer-phase2-human-vector-nonlinear.csv", index=False)
    composition = core_composition("all_donors")
    composition = composition[composition["species"] == "Human"].copy()
    composition_results = analyze(composition, "outcome", "cell_fraction", False)
    composition_results.to_csv(RESULTS / "reviewer-phase2-human-composition-nonlinear.csv", index=False)
    distributions = []
    for analysis, data in (("expression_vector", vector), ("cell_fraction", composition)):
        for cohort, part in data.groupby("cohort"):
            metadata = part.groupby("donor")["age"].first()
            distributions.append({
                "analysis": analysis, "cohort": cohort, "n_donors": len(metadata),
                "age_min": metadata.min(), "age_q1": metadata.quantile(0.25),
                "age_median": metadata.median(), "age_q3": metadata.quantile(0.75),
                "age_max": metadata.max(), "all_adult": bool((metadata >= 20).all()),
            })
    pd.DataFrame(distributions).to_csv(RESULTS / "reviewer-phase2-human-age-distributions.csv", index=False)
    for result in (vector_results, composition_results):
        print(result.groupby("cohort").agg(n=("cell_type", "size"),
              linear_fdr=("q_linear_bh", lambda x: int((x < 0.05).sum())),
              nonlinear_fdr=("q_nonlinear_bh", lambda x: int((x < 0.05).sum()))).to_string(), flush=True)


if __name__ == "__main__":
    main()
