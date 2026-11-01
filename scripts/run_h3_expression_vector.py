#!/usr/bin/env python3
"""Compare complete young and old expression vectors without UMI downsampling."""

from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_h3_transcriptional_disorder import (
    AGE_MAP,
    TARGET_CELL_TYPES,
    excluded_gene_mask,
    reference_stats,
    without_donor,
)

RESULTS = REPO / "results"


def hellinger_distance(matrix, reference_proportions):
    """Hellinger distance from each count row to a reference probability vector."""
    matrix = sp.csr_matrix(matrix, dtype=float)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError("Hellinger distance requires positive row sums")
    reference_proportions = np.asarray(reference_proportions, dtype=float)
    if reference_proportions.sum() <= 0:
        raise ValueError("Reference profile has no counts")
    reference_proportions = reference_proportions / reference_proportions.sum()
    root_matrix = matrix.copy()
    root_matrix.data = np.sqrt(root_matrix.data)
    affinity = np.asarray(root_matrix @ np.sqrt(reference_proportions)).ravel() / np.sqrt(totals)
    return np.sqrt(np.maximum(0.0, 1.0 - affinity))


def raw_euclidean_distance(matrix, reference_mean_counts):
    """Euclidean distance from each raw count row to the young mean count vector."""
    matrix = sp.csr_matrix(matrix, dtype=float)
    reference_mean_counts = np.asarray(reference_mean_counts, dtype=float)
    row_squares = np.asarray(matrix.multiply(matrix).sum(axis=1)).ravel()
    cross_products = np.asarray(matrix @ reference_mean_counts).ravel()
    squared = row_squares + np.dot(reference_mean_counts, reference_mean_counts) - 2 * cross_products
    return np.sqrt(np.maximum(0.0, squared))


def donor_regression(records, metric, cell_type, donor_set, adjust_depth):
    predictors = ["age_months", "sex_male"]
    model = "age_plus_sex"
    if adjust_depth:
        predictors.append("mean_log_umi")
        model += "_plus_log_umi"
    X = sm.add_constant(records[predictors].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        predictors.remove("sex_male")
        model = model.replace("age_plus_sex", "age_only_rank_deficient_sex")
        X = sm.add_constant(records[predictors].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    fit = sm.OLS(records[metric].astype(float), X).fit()
    ci = fit.conf_int().loc["age_months"]
    return {
        "cell_type": cell_type,
        "donor_set": donor_set,
        "metric": metric,
        "model": model,
        "n_donors": len(records),
        "n_cells": int(records["n_cells"].sum()),
        "slope_per_month": fit.params["age_months"],
        "se": fit.bse["age_months"],
        "ci95_low": ci.iloc[0],
        "ci95_high": ci.iloc[1],
        "pvalue": fit.pvalues["age_months"],
        "r2": fit.rsquared,
        "log_umi_slope": fit.params.get("mean_log_umi", np.nan),
        "log_umi_pvalue": fit.pvalues.get("mean_log_umi", np.nan),
    }


def global_regression(donors, metric, donor_set, adjust_depth):
    rows = []
    for cell_type, part in donors.groupby("cell_type"):
        part = part.dropna(subset=[metric]).copy()
        if len(part) < 6:
            continue
        scale = part[metric].std(ddof=1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        part["standardized_outcome"] = (part[metric] - part[metric].mean()) / scale
        part["age_contrast"] = (part["age_months"] - 3.0) / 15.0
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        part["depth_centered"] = part["mean_log_umi"] - part["mean_log_umi"].mean()
        part["cell_type"] = cell_type
        rows.append(part)
    data = pd.concat(rows, ignore_index=True)
    if donor_set == "exclude_pooled_ids":
        pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                     .query("is_pooled_id")["mouse.id"].astype(str))
        data = data[~data["donor"].isin(pooled)].copy()
    cell_dummies = pd.get_dummies(data["cell_type"], prefix="cell_type", drop_first=True, dtype=float)
    predictors = ["age_contrast", "sex_centered"]
    model = "cell_type_age_sex_clustered_by_animal"
    if adjust_depth:
        predictors.append("depth_centered")
        model += "_plus_log_umi"
    X = pd.concat([data[predictors].astype(float), cell_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    fit = sm.OLS(data["standardized_outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": data["donor"]}
    )
    ci = fit.conf_int().loc["age_contrast"]
    return {
        "donor_set": donor_set,
        "metric": metric,
        "model": model,
        "n_donors": data["donor"].nunique(),
        "n_cell_types": data["cell_type"].nunique(),
        "n_animal_cell_types": len(data),
        "effect_15_months_sd": fit.params["age_contrast"],
        "se": fit.bse["age_contrast"],
        "ci95_low": ci.iloc[0],
        "ci95_high": ci.iloc[1],
        "pvalue": fit.pvalues["age_contrast"],
    }


def main():
    adata = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    lung_mask = (
        (adata.obs["tissue"] == "Lung")
        & ~adata.obs["predicted_doublet"]
        & (adata.obs["n_counts"] >= 1500)
        & (adata.obs["cell_ontology_class"].astype(str) != "unknown")
    )
    lung = adata[lung_mask].to_memory()
    lung.obs["age_months"] = lung.obs["age"].astype(str).map(AGE_MAP)
    lung.obs["sex_male"] = (lung.obs["sex"].astype(str).str.lower() == "male").astype(float)
    lung.obs["mouse.id"] = lung.obs["mouse.id"].astype(str)
    gene_names = lung.var_names.astype(str).to_numpy()
    technical, _ = excluded_gene_mask(gene_names)
    lung = lung[:, ~technical].copy()
    lung.X = lung.X.tocsr()
    print(f"Loaded {lung.n_obs} cells and {lung.n_vars} eligible genes without downsampling", flush=True)

    young = lung[lung.obs["age_months"] <= 3].copy()
    stats = reference_stats(young.X, young.obs)
    donor_rows = []
    for cell_type in TARGET_CELL_TYPES:
        ct_mask = (lung.obs["cell_ontology_class"] == cell_type).to_numpy()
        ct_obs = lung.obs.loc[ct_mask]
        ct_matrix = lung.X[ct_mask]
        for donor in ct_obs["mouse.id"].unique():
            donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
            if donor_mask.sum() < 20:
                continue
            donor_obs = ct_obs.loc[donor_mask]
            donor_matrix = ct_matrix[donor_mask]
            reference = without_donor(stats[cell_type], str(donor))
            if reference["n_cells"] < 50 or reference["library"] <= 0:
                continue
            profile_distance = hellinger_distance(donor_matrix, reference["counts"])
            raw_distance = raw_euclidean_distance(
                donor_matrix, reference["counts"] / reference["n_cells"]
            )
            totals = np.asarray(donor_matrix.sum(axis=1)).ravel()
            donor_rows.append({
                "cell_type": cell_type,
                "donor": donor,
                "age_months": float(donor_obs["age_months"].iloc[0]),
                "sex_male": float(donor_obs["sex_male"].iloc[0]),
                "n_cells": int(donor_matrix.shape[0]),
                "mean_umi": float(totals.mean()),
                "median_umi": float(np.median(totals)),
                "mean_log_umi": float(np.log(totals).mean()),
                "profile_hellinger": float(profile_distance.mean()),
                "raw_count_euclidean": float(raw_distance.mean()),
                "reference_young_animals": int(young.obs.loc[
                    young.obs["cell_ontology_class"] == cell_type, "mouse.id"
                ].nunique() - (donor in set(young.obs["mouse.id"]))),
            })
        print(f"Completed {cell_type}", flush=True)

    donors = pd.DataFrame(donor_rows)
    donors.to_csv(RESULTS / "h3-expression-vector-donors.csv", index=False)
    pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                 .query("is_pooled_id")["mouse.id"].astype(str))
    result_rows = []
    global_rows = []
    for cell_type in TARGET_CELL_TYPES:
        ct_donors = donors[donors["cell_type"] == cell_type]
        for donor_set, selected in (
            ("all_donor_ids", ct_donors),
            ("exclude_pooled_ids", ct_donors[~ct_donors["donor"].isin(pooled)]),
        ):
            for metric in ("profile_hellinger", "raw_count_euclidean"):
                for adjust_depth in (False, True):
                    row = donor_regression(selected, metric, cell_type, donor_set, adjust_depth)
                    if row:
                        result_rows.append(row)
    for donor_set in ("all_donor_ids", "exclude_pooled_ids"):
        for metric in ("profile_hellinger", "raw_count_euclidean"):
            for adjust_depth in (False, True):
                global_rows.append(global_regression(donors, metric, donor_set, adjust_depth))

    results = pd.DataFrame(result_rows)
    results["q_bh"] = np.nan
    for (_, _, _), indices in results.groupby(["donor_set", "metric", "model"]).groups.items():
        results.loc[indices, "q_bh"] = multipletests(
            results.loc[indices, "pvalue"], method="fdr_bh"
        )[1]
    results.to_csv(RESULTS / "h3-expression-vector.csv", index=False)
    pd.DataFrame(global_rows).to_csv(RESULTS / "h3-expression-vector-global.csv", index=False)
    primary = results[
        (results["donor_set"] == "all_donor_ids")
        & (results["metric"] == "profile_hellinger")
    ]
    print(primary[["cell_type", "model", "n_donors", "slope_per_month",
                   "ci95_low", "ci95_high", "pvalue", "q_bh"]].to_string(index=False), flush=True)
    print(pd.DataFrame(global_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
