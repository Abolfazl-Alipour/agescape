#!/usr/bin/env python3

from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_h3_transcriptional_disorder import (
    AGE_MAP, TARGET_CELL_TYPES, excluded_gene_mask, offprogram_markers, reference_stats, shannon_entropy,
)

RESULTS = REPO / "results"


def global_fit(donors, metric, depth, library_covariate):
    data = donors.dropna(subset=[metric]).copy()
    rows = []
    for cell_type, part in data.groupby("cell_type"):
        part = part.copy()
        sd = part[metric].std(ddof=1)
        if not np.isfinite(sd) or sd <= 0:
            continue
        part["outcome"] = (part[metric] - part[metric].mean()) / sd
        part["age15"] = (part["age_months"] - 3) / 15
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        part["depth_centered"] = part["mean_log_umi_actual"] - part["mean_log_umi_actual"].mean()
        rows.append(part)
    data = pd.concat(rows)
    dummies = pd.get_dummies(data["cell_type"], drop_first=True, dtype=float)
    predictors = ["age15", "sex_centered"] + (["depth_centered"] if library_covariate else [])
    X = sm.add_constant(pd.concat([data[predictors], dummies], axis=1), has_constant="add")
    fit = sm.OLS(data["outcome"], X).fit(cov_type="cluster", cov_kwds={"groups": data["donor"]}, use_t=True)
    ci = fit.conf_int().loc["age15"]
    return {
        "depth": depth, "metric": metric,
        "model": "age_sex_library" if library_covariate else "age_sex",
        "n_donors": data["donor"].nunique(), "n_cell_types": data["cell_type"].nunique(),
        "n_donor_units": len(data), "effect_15mo_sd": fit.params["age15"], "se": fit.bse["age15"],
        "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1], "pvalue": fit.pvalues["age15"],
    }


def cell_type_fits(donors, metric, depth, library_covariate):
    rows = []
    for cell_type, part in donors.dropna(subset=[metric]).groupby("cell_type"):
        if len(part) < 6:
            continue
        predictors = ["age_months", "sex_male"] + (["mean_log_umi_actual"] if library_covariate else [])
        X = sm.add_constant(part[predictors].astype(float), has_constant="add")
        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            predictors.remove("sex_male")
            X = sm.add_constant(part[predictors].astype(float), has_constant="add")
        fit = sm.OLS(part[metric], X).fit()
        ci = fit.conf_int().loc["age_months"]
        rows.append({
            "depth": depth, "metric": metric, "cell_type": cell_type,
            "model": "age_sex_library" if library_covariate else "age_sex",
            "n_donors": len(part), "slope_per_month": fit.params["age_months"], "se": fit.bse["age_months"],
            "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1], "pvalue": fit.pvalues["age_months"],
        })
    if rows:
        qvalues = multipletests([row["pvalue"] for row in rows], method="fdr_bh")[1]
        for row, qvalue in zip(rows, qvalues):
            row["q_bh"] = qvalue
    return rows


def main():
    atlas = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    mask = ((atlas.obs["tissue"] == "Lung") & ~atlas.obs["predicted_doublet"]
            & (atlas.obs["n_counts"] >= 500)
            & (atlas.obs["cell_ontology_class"].astype(str) != "unknown"))
    lung = atlas[mask].to_memory()
    lung.obs["age_months"] = lung.obs["age"].astype(str).map(AGE_MAP)
    lung.obs["sex_male"] = (lung.obs["sex"].astype(str).str.lower() == "male").astype(float)
    lung.obs["mouse.id"] = lung.obs["mouse.id"].astype(str)
    lung.obs["actual_n_counts"] = np.asarray(lung.X.sum(axis=1)).ravel()
    young = lung[lung.obs["age_months"].isin([1, 3]) & (lung.obs["actual_n_counts"] >= 1500)].copy()
    stats = reference_stats(young.X.tocsr(), young.obs)
    gene_names = lung.var_names.astype(str).to_numpy()
    technical, _ = excluded_gene_mask(gene_names)
    marker_cache = {}
    for cell_type in TARGET_CELL_TYPES:
        for donor in lung.obs["mouse.id"].unique():
            marker_cache[(cell_type, donor)], _ = offprogram_markers(
                stats, cell_type, donor, gene_names, stringent=False
            )
    all_donors = []
    for depth, minimum in ((500, 500), (1500, 1500), (3000, 3000), ("full", 1500)):
        selected = lung[lung.obs["actual_n_counts"] >= minimum].copy()
        if depth == "full":
            matrix = selected.X.tocsr()
        else:
            sc.pp.downsample_counts(selected, counts_per_cell=depth, replace=False, random_state=42)
            matrix = selected.X.tocsr()
            assert np.all(np.asarray(matrix.sum(axis=1)).ravel() == depth)
        selected.obs["entropy"] = shannon_entropy(matrix[:, ~technical])
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        for cell_type in TARGET_CELL_TYPES:
            ct_mask = (selected.obs["cell_ontology_class"].astype(str) == cell_type).to_numpy()
            ct_obs = selected.obs[ct_mask]
            ct_matrix = matrix[ct_mask]
            for donor in ct_obs["mouse.id"].unique():
                donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
                if donor_mask.sum() < 20:
                    continue
                donor_obs = ct_obs.loc[donor_mask]
                donor_matrix = ct_matrix[donor_mask]
                markers = marker_cache[(cell_type, donor)]
                offprogram = (np.asarray(donor_matrix[:, markers].sum(axis=1)).ravel()
                              / totals[ct_mask][donor_mask] if len(markers)
                              else np.full(donor_matrix.shape[0], np.nan))
                all_donors.append({
                    "depth": depth, "cell_type": cell_type, "donor": donor,
                    "age_months": float(donor_obs["age_months"].iloc[0]),
                    "sex_male": float(donor_obs["sex_male"].iloc[0]),
                    "n_cells": int(donor_mask.sum()),
                    "mean_log_umi_actual": float(np.log(donor_obs["actual_n_counts"]).mean()),
                    "entropy": float(donor_obs["entropy"].mean()),
                    "offprogram_fraction": float(np.nanmean(offprogram)) if len(markers) else np.nan,
                })
        print(f"Completed depth {depth}", flush=True)
    donors = pd.DataFrame(all_donors)
    donors.to_csv(RESULTS / "reviewer-phase2-depth-donor-results.csv", index=False)
    global_rows = []
    cell_rows = []
    for depth in donors["depth"].unique():
        part = donors[donors["depth"] == depth]
        for library_covariate in (False, True):
            for metric in ("entropy", "offprogram_fraction"):
                global_rows.append(global_fit(part, metric, depth, library_covariate))
                cell_rows.extend(cell_type_fits(part, metric, depth, library_covariate))
    pd.DataFrame(global_rows).to_csv(RESULTS / "reviewer-phase2-depth-global.csv", index=False)
    pd.DataFrame(cell_rows).to_csv(RESULTS / "reviewer-phase2-depth-cell-types.csv", index=False)


if __name__ == "__main__":
    main()
