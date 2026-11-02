#!/usr/bin/env python3

import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def clr_table(counts, pseudocount=0.5):
    adjusted = counts.astype(float) + pseudocount
    proportions = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_proportions = np.log(proportions)
    return log_proportions.sub(log_proportions.mean(axis=1), axis=0)


def composition_counts(obs, sample_column, cell_column):
    counts = (obs.groupby([sample_column, cell_column], observed=True).size()
              .unstack(fill_value=0).sort_index())
    return counts.loc[:, (counts > 0).sum(axis=0) >= 6]


def composition_tests(counts, metadata, species, organ, cohort, age_unit,
                      sensitivity="all_samples", observed_only=False, min_samples=6,
                      adjust_total=False, extra_covariates=None):
    clr = clr_table(counts)
    metadata = metadata.copy()
    metadata["log_total_cells"] = np.log(counts.sum(axis=1))
    rows = []
    for cell_type in counts.columns:
        donors = counts.index
        if observed_only:
            donors = donors[counts.loc[donors, cell_type] > 0]
        if len(donors) < min_samples:
            continue
        predictors = ["age"]
        if metadata.loc[donors, "sex_male"].nunique() > 1:
            predictors.append("sex_male")
        if adjust_total:
            metadata.loc[donors, "log_total_cells"] -= metadata.loc[donors, "log_total_cells"].mean()
            predictors.append("log_total_cells")
        if extra_covariates:
            predictors.extend(extra_covariates)
        X = sm.add_constant(metadata.loc[donors, predictors].astype(float), has_constant="add")
        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            continue
        fit = sm.OLS(clr.loc[donors, cell_type], X).fit()
        rows.append({
            "species": species, "organ": organ, "cohort": cohort, "cell_type": cell_type,
            "sensitivity": sensitivity, "effect_type": f"CLR_slope_per_{age_unit}",
            "effect": fit.params["age"], "se": fit.bse["age"],
            "ci95_low": fit.conf_int().loc["age", 0],
            "ci95_high": fit.conf_int().loc["age", 1],
            "pvalue": fit.pvalues["age"], "n_samples": len(donors),
            "n_present": int((counts.loc[donors, cell_type] > 0).sum()),
            "total_cells": int(counts.loc[donors, cell_type].sum()),
        })
    if not rows:
        return []
    qvalues = multipletests([row["pvalue"] for row in rows], method="fdr_bh")[1]
    for row, qvalue in zip(rows, qvalues):
        row["q_bh"] = qvalue
    return rows


def afca_sample_id(cell_names):
    return pd.Index(cell_names.astype(str)).str.split("_", n=1).str[1]


def zebrafish_sample_metadata(filename):
    match = re.search(r"mo(\d+)_fish(\d+)_rep(\d+)", Path(filename).name)
    if not match:
        raise ValueError(f"Unrecognized zebrafish sample name: {filename}")
    age, fish, replicate = map(int, match.groups())
    return {
        "sample": f"mo{age}_fish{fish}_rep{replicate}",
        "donor": f"mo{age}_fish{fish}",
        "age": age,
        "replicate": replicate,
    }


def write_vector_results(donors, age_scales, prefix, age_units=None):
    from scripts.run_cross_species_expression_vector import fit_global, fit_unit

    result_rows = []
    global_rows = []
    for (species, tissue, cohort, cell_type), part in donors.groupby(
        ["species", "tissue", "cohort", "cell_type"]
    ):
        for model in ("unadjusted", "depth_adjusted"):
            row = fit_unit(part, "profile_hellinger", species, tissue, cohort, cell_type, model)
            if row:
                result_rows.append(row)
    for (species, tissue, cohort), part in donors.groupby(["species", "tissue", "cohort"]):
        for model in ("unadjusted", "depth_adjusted"):
            row = fit_global(part, "profile_hellinger", species, tissue, cohort, model,
                             age_scales[species], "all_samples")
            if row:
                global_rows.append(row)
    for species, part in donors.groupby("species"):
        for model in ("unadjusted", "depth_adjusted"):
            row = fit_global(part, "profile_hellinger", species, "All tested tissues", "All cohorts",
                             model, age_scales[species], "all_samples")
            if row:
                global_rows.append(row)
    results = pd.DataFrame(result_rows)
    if not results.empty:
        results["q_bh"] = np.nan
        for _, indices in results.groupby(["cohort", "metric", "model", "species"]).groups.items():
            results.loc[indices, "q_bh"] = multipletests(
                results.loc[indices, "pvalue"], method="fdr_bh"
            )[1]
    globals_out = pd.DataFrame(global_rows)
    results = results.rename(columns={"slope_per_year": "slope_per_age_unit"})
    globals_out = globals_out.rename(columns={"age_scale_years": "age_scale"})
    if age_units:
        results["age_unit"] = results["species"].map(age_units)
        globals_out["age_unit"] = globals_out["species"].map(age_units)
    results.to_csv(RESULTS / f"{prefix}-tests.csv", index=False)
    globals_out.to_csv(RESULTS / f"{prefix}-global.csv", index=False)
    return results, globals_out
