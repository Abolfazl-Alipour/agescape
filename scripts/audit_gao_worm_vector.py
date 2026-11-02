#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RESULTS = Path(__file__).resolve().parent.parent / "results"
REPRODUCTIVE = {"Germline", "Embroynic_cells", "Sperms"}


def fit_unit(part):
    X = sm.add_constant(part[["age", "mean_log_umi", "platform_nextseq"]], has_constant="add")
    return sm.OLS(part["profile_hellinger"], X).fit()


def permuted_age(part, rng):
    result = part.copy()
    result["age"] = result["age"].to_numpy()
    for _, indices in result.groupby("platform_nextseq").groups.items():
        values = result.loc[indices, "age"].to_numpy().copy()
        result.loc[indices, "age"] = rng.permutation(values)
    return result


def unit_tests(donors, suffix):
    rng = np.random.default_rng(20261102)
    rows = []
    for cell_type, part in donors.groupby("cell_type"):
        if len(part) < 6:
            continue
        fit = fit_unit(part)
        permuted = [fit_unit(permuted_age(part, rng)).tvalues["age"] for _ in range(5000)]
        row = {
            "cell_type": cell_type, "n_pools": len(part),
            "age_slope_per_day": fit.params["age"], "age_p": fit.pvalues["age"],
            "platform_effect": fit.params["platform_nextseq"],
            "platform_p": fit.pvalues["platform_nextseq"],
            "stratified_permutation_p": float(
                (1 + np.sum(np.abs(permuted) >= abs(fit.tvalues["age"]) - 1e-12)) / 5001
            ),
        }
        for platform, selected in part.groupby("platform_nextseq"):
            X = sm.add_constant(selected[["age", "mean_log_umi"]], has_constant="add")
            platform_fit = sm.OLS(selected["profile_hellinger"], X).fit()
            row[f"age_slope_platform_{int(platform)}"] = platform_fit.params["age"]
        rows.append(row)
    result = pd.DataFrame(rows)
    result["q_bh"] = multipletests(result["age_p"], method="fdr_bh")[1]
    result["permutation_q_bh"] = multipletests(
        result["stratified_permutation_p"], method="fdr_bh"
    )[1]
    result.to_csv(RESULTS / f"extended-worm-vector-platform-audit{suffix}.csv", index=False)


def prepare_global(donors):
    data = donors.copy()
    data["unit"] = data["cell_type"]
    standardized = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 6 or part["profile_hellinger"].std(ddof=1) <= 0:
            continue
        part["outcome"] = ((part["profile_hellinger"] - part["profile_hellinger"].mean())
                           / part["profile_hellinger"].std(ddof=1))
        part["age_scaled"] = (part["age"] - 1) / 14
        for column in ("mean_log_umi", "platform_nextseq"):
            part[f"{column}_centered"] = part[column] - part[column].mean()
        standardized.append(part)
    return pd.concat(standardized, ignore_index=True)


def global_fit(data, cluster=True):
    dummies = pd.get_dummies(data["unit"], drop_first=True, dtype=float)
    predictors = ["age_scaled", "mean_log_umi_centered", "platform_nextseq_centered"]
    X = sm.add_constant(pd.concat([data[predictors], dummies], axis=1), has_constant="add")
    if cluster:
        return sm.OLS(data["outcome"], X).fit(
            cov_type="cluster", cov_kwds={"groups": data["donor"]}, use_t=True
        )
    return sm.OLS(data["outcome"], X).fit()


def permuted_global_age(data, rng):
    result = data.copy()
    donor_meta = result[["donor", "platform_nextseq", "age"]].drop_duplicates("donor")
    for _, part in donor_meta.groupby("platform_nextseq"):
        shuffled = rng.permutation(part["age"].to_numpy())
        mapping = dict(zip(part["donor"], shuffled))
        selected = result["donor"].isin(part["donor"])
        result.loc[selected, "age_scaled"] = (
            result.loc[selected, "donor"].map(mapping) - 1
        ) / 14
    return result


def global_tests(donors, suffix):
    rng = np.random.default_rng(20261103)
    rows = []
    for analysis, selected in (
        ("all_annotated_cells", donors),
        ("somatic_only", donors[~donors["cell_type"].isin(REPRODUCTIVE)]),
    ):
        data = prepare_global(selected)
        fit = global_fit(data)
        ordinary = global_fit(data, cluster=False)
        permuted_t = []
        for _ in range(2000):
            shuffled = permuted_global_age(data, rng)
            permuted_t.append(global_fit(shuffled, cluster=False).tvalues["age_scaled"])
        rows.append({
            "analysis": analysis, "n_pools": data["donor"].nunique(),
            "n_units": data["unit"].nunique(), "effect_sd_per_14_days": fit.params["age_scaled"],
            "se": fit.bse["age_scaled"], "pvalue": fit.pvalues["age_scaled"],
            "platform_stratified_permutation_p": float(
                (1 + np.sum(np.abs(permuted_t) >= abs(ordinary.tvalues["age_scaled"]) - 1e-12)) / 2001
            ),
        })
    pd.DataFrame(rows).to_csv(RESULTS / f"extended-worm-vector-global-platform-audit{suffix}.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--suffix", default="")
    parser.add_argument("--global-only", action="store_true")
    args = parser.parse_args()
    donors = pd.read_csv(RESULTS / f"extended-worm-vector-donors{args.suffix}.csv")
    if not args.global_only:
        unit_tests(donors, args.suffix)
    global_tests(donors, args.suffix)


if __name__ == "__main__":
    main()
