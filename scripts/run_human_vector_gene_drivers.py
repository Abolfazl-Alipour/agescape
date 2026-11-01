#!/usr/bin/env python3
"""Explore donor-pseudobulk gene drivers of human vector-distance associations."""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import AGE_RE, technical_gene_mask
from scripts.run_m9_pipeline import COHORTS

RESULTS = REPO / "results"
UNITS = [
    ("MSSM_DLPFC", "GABAergic neuron"),
    ("MSSM_DLPFC", "VIP GABAergic cortical interneuron"),
    ("MSSM_DLPFC", "pvalb GABAergic cortical interneuron"),
    ("MSSM_DLPFC", "sst GABAergic cortical interneuron"),
    ("Lung_Integrated", "acinar cell"),
    ("Lung_Integrated", "nasal mucosa goblet cell"),
]


def pseudobulk(adata):
    adata.obs["donor_id"] = adata.obs["donor_id"].astype(str)
    adata.obs["age"] = adata.obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
    adata = adata[adata.obs["age"].between(20, 120)].copy()
    adata.var["symbol"] = adata.var["feature_name"].astype(str)
    adata.var_names = adata.var["symbol"].values
    adata.var_names_make_unique()
    keep = ~technical_gene_mask(adata.var_names)
    adata = adata[:, keep].copy()
    matrix = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    count_rows = []
    metadata_rows = []
    for donor in adata.obs["donor_id"].unique():
        mask = (adata.obs["donor_id"] == donor).to_numpy()
        if mask.sum() < 20:
            continue
        donor_matrix = matrix[mask]
        count_rows.append(np.asarray(donor_matrix.sum(axis=0)).ravel())
        metadata_rows.append({
            "donor": donor,
            "age": float(adata.obs.loc[mask, "age"].iloc[0]),
            "sex_male": float(adata.obs.loc[mask, "sex"].astype(str).str.lower().iloc[0] == "male"),
            "n_cells": int(mask.sum()),
            "mean_log_umi": float(np.log(np.asarray(donor_matrix.sum(axis=1)).ravel()).mean()),
        })
    counts = np.vstack(count_rows)
    metadata = pd.DataFrame(metadata_rows)
    libraries = counts.sum(axis=1)
    log_cp10k = np.log1p(counts / libraries[:, None] * 1e4)
    return log_cp10k, metadata, adata.var_names.to_numpy()


def gene_regressions(expression, metadata, gene_names, cohort, cell_type):
    from scipy import stats

    rows = []
    for model_name, predictors in (
        ("age_plus_sex", ["age", "sex_male"]),
        ("age_plus_sex_depth", ["age", "sex_male", "mean_log_umi"]),
    ):
        X = sm.add_constant(metadata[predictors].astype(float), has_constant="add").to_numpy()
        if np.linalg.matrix_rank(X) < X.shape[1]:
            predictors.remove("sex_male")
            X = sm.add_constant(metadata[predictors].astype(float), has_constant="add").to_numpy()
        age_index = predictors.index("age") + 1
        coefficients = np.linalg.pinv(X) @ expression
        residuals = expression - X @ coefficients
        degrees_freedom = X.shape[0] - np.linalg.matrix_rank(X)
        residual_variance = np.sum(residuals ** 2, axis=0) / degrees_freedom
        age_variance = np.linalg.pinv(X.T @ X)[age_index, age_index]
        standard_errors = np.sqrt(residual_variance * age_variance)
        slopes = coefficients[age_index]
        t_statistics = np.divide(
            slopes, standard_errors, out=np.zeros_like(slopes), where=standard_errors > 0
        )
        pvalues = 2 * stats.t.sf(np.abs(t_statistics), degrees_freedom)
        critical = stats.t.ppf(0.975, degrees_freedom)
        qvalues = multipletests(pvalues, method="fdr_bh")[1]
        for gene_index, gene in enumerate(gene_names):
            rows.append({
                "cohort": cohort, "cell_type": cell_type, "model": model_name,
                "gene": gene, "n_donors": len(metadata),
                "slope_log_cp10k_per_year": slopes[gene_index],
                "se": standard_errors[gene_index],
                "ci95_low": slopes[gene_index] - critical * standard_errors[gene_index],
                "ci95_high": slopes[gene_index] + critical * standard_errors[gene_index],
                "pvalue": pvalues[gene_index], "q_bh": qvalues[gene_index],
            })
    return rows


def main():
    import cellxgene_census

    census = cellxgene_census.open_soma(census_version="2025-11-08")
    dataset_lookup = {label: dataset_id for dataset_id, (label, _) in COHORTS.items()}
    all_rows = []
    coverage_rows = []
    for cohort, cell_type in UNITS:
        dataset_id = dataset_lookup[cohort]
        adata = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=(
                f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                f"and raw_sum >= 1500 and cell_type == '{cell_type}'"
            ),
            obs_column_names=["donor_id", "development_stage", "sex"],
            var_column_names=["feature_name"],
        )
        expression, metadata, gene_names = pseudobulk(adata)
        print(f"{cohort} | {cell_type}: {len(metadata)} donors, {len(gene_names)} genes", flush=True)
        coverage_rows.append({"cohort": cohort, "cell_type": cell_type,
                              "n_donors": len(metadata), "n_genes": len(gene_names)})
        all_rows.extend(gene_regressions(expression, metadata, gene_names, cohort, cell_type))
    results = pd.DataFrame(all_rows)
    results.to_csv(RESULTS / "h3-human-vector-gene-drivers-all.csv", index=False)
    top = (results.sort_values(["cohort", "cell_type", "model", "pvalue"])
           .groupby(["cohort", "cell_type", "model"], group_keys=False).head(50))
    top.to_csv(RESULTS / "h3-human-vector-gene-drivers-top50.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(RESULTS / "h3-human-vector-gene-driver-coverage.csv", index=False)
    summary = (results.groupby(["cohort", "cell_type", "model"])
               .agg(n_genes=("gene", "size"), n_fdr=("q_bh", lambda values: int((values < 0.05).sum())),
                    min_q=("q_bh", "min")).reset_index())
    summary.to_csv(RESULTS / "h3-human-vector-gene-driver-summary.csv", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
