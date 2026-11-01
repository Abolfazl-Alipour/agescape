#!/usr/bin/env python3
"""Check assay and suspension confounding in human lung vector associations."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import AGE_RE
from scripts.run_m9_pipeline import COHORTS

RESULTS = REPO / "results"
UNITS = [("Lung_Integrated", "acinar cell"), ("Lung_Integrated", "nasal mucosa goblet cell")]


def main():
    import cellxgene_census

    dataset_id = next(key for key, value in COHORTS.items() if value[0] == "Lung_Integrated")
    census = cellxgene_census.open_soma(census_version="2025-11-08")
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    rows = []
    metadata_rows = []
    for cohort, cell_type in UNITS:
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter=(
                f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                f"and raw_sum >= 1500 and cell_type == '{cell_type}'"
            ),
            column_names=["donor_id", "assay", "suspension_type", "development_stage"],
        )
        obs["donor_id"] = obs["donor_id"].astype(str)
        metadata = (obs.groupby("donor_id", observed=True)
                    .agg(assay=("assay", lambda values: values.astype(str).mode().iloc[0]),
                         suspension=("suspension_type", lambda values: values.astype(str).mode().iloc[0]),
                         n_assays=("assay", lambda values: values.astype(str).nunique()))
                    .reset_index())
        metadata["cohort"] = cohort
        metadata["cell_type"] = cell_type
        metadata_rows.append(metadata)
        frame = donors[(donors["cohort"] == cohort) & (donors["cell_type"] == cell_type)].merge(
            metadata, left_on="donor", right_on="donor_id"
        )
        for model_name, include_assay in (("depth_adjusted", False), ("depth_and_assay_adjusted", True)):
            predictors = frame[["age", "sex_male", "mean_log_umi"]].astype(float)
            if include_assay:
                assay_dummies = pd.get_dummies(frame["assay"].astype(str), prefix="assay", drop_first=True, dtype=float)
                suspension_dummies = pd.get_dummies(frame["suspension"].astype(str), prefix="suspension", drop_first=True, dtype=float)
                predictors = pd.concat([predictors, assay_dummies, suspension_dummies], axis=1)
            predictors = predictors.loc[:, predictors.var() > 1e-12]
            X = sm.add_constant(predictors, has_constant="add")
            identifiable = np.linalg.matrix_rank(X.to_numpy()) == X.shape[1]
            if not identifiable:
                raise ValueError(f"Non-identifiable assay model: {cell_type}, {model_name}")
            fit = sm.OLS(frame["profile_hellinger"], X).fit()
            rows.append({
                "cohort": cohort, "cell_type": cell_type, "model": model_name,
                "n_donors": len(frame), "n_assays": frame["assay"].nunique(),
                "n_suspension_types": frame["suspension"].nunique(), "full_rank": identifiable,
                "age_slope": fit.params["age"], "age_pvalue": fit.pvalues["age"],
                "depth_slope": fit.params["mean_log_umi"], "depth_pvalue": fit.pvalues["mean_log_umi"],
            })
    pd.DataFrame(rows).to_csv(RESULTS / "h3-human-lung-assay-audit.csv", index=False)
    pd.concat(metadata_rows, ignore_index=True).to_csv(RESULTS / "h3-human-lung-donor-assays.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
