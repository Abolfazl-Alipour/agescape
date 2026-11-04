#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.proliferation import POST_MITOTIC, PROLIFERATIVE, classify_cell_type

RESULTS = REPO / "results"
AGE_SCALE = {
    "Mouse": 0.4,
    "Human": 12.25,
    "Fruit fly": 40.0,
    "C. elegans": 14.0,
    "Zebrafish": 22.0,
    "Naked mole-rat": 3.1,
}


def add_classification(frame):
    labels = frame[["species", "cell_type"]].drop_duplicates()
    curation = pd.DataFrame([
        classify_cell_type(row.cell_type, row.species) for row in labels.itertuples()
    ])
    return frame.merge(curation, on=["species", "cell_type"], how="left", validate="many_to_one")


def core_composition(sensitivity):
    from scripts.run_organ_composition import AGE_MAP, clr_table, count_table
    import anndata as ad

    frames = []
    adata = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    obs = adata.obs[~adata.obs["predicted_doublet"]].copy()
    obs["age"] = obs["age"].astype(str).map(AGE_MAP)
    obs["sex_male"] = (obs["sex"].astype(str).str.lower() == "male").astype(float)
    for tissue in sorted(obs["tissue"].unique()):
        part = obs[(obs["tissue"] == tissue) & (obs["cell_ontology_class"].astype(str) != "unknown")]
        counts = count_table(part, "mouse.id", "cell_ontology_class")
        metadata = part.groupby("mouse.id").agg(age=("age", "first"), sex_male=("sex_male", "first"))
        clr = clr_table(counts)
        for cell_type in counts.columns:
            donors = counts.index
            if sensitivity == "observed_records_only":
                donors = donors[counts[cell_type] > 0]
            if len(donors) < 6:
                continue
            frames.append(pd.DataFrame({
                "species": "Mouse", "tissue": tissue, "cohort": "TMS", "cell_type": cell_type,
                "donor": donors.astype(str), "age": metadata.loc[donors, "age"].to_numpy(),
                "sex_male": metadata.loc[donors, "sex_male"].to_numpy(),
                "outcome": clr.loc[donors, cell_type].to_numpy(),
            }))
    source = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    for cohort in sorted(source["cohort"].unique()):
        part = source[(source["cohort"] == cohort) & (source["age"] >= 20)].copy()
        tissue = part["group"].iloc[0]
        counts = part.set_index(["donor_id", "cell_type"])["n_cells"].unstack(fill_value=0).sort_index()
        counts = counts.loc[:, (counts > 0).sum(axis=0) >= 6]
        metadata = part.groupby("donor_id").agg(age=("age", "first"), sex_male=("sex", "first"))
        metadata["sex_male"] = (metadata["sex_male"].astype(str).str.lower() == "male").astype(float)
        clr = clr_table(counts)
        for cell_type in counts.columns:
            donors = counts.index
            if sensitivity == "observed_records_only":
                donors = donors[counts[cell_type] > 0]
            if len(donors) < 6:
                continue
            frames.append(pd.DataFrame({
                "species": "Human", "tissue": tissue, "cohort": cohort, "cell_type": cell_type,
                "donor": donors.astype(str), "age": metadata.loc[donors, "age"].to_numpy(),
                "sex_male": metadata.loc[donors, "sex_male"].to_numpy(),
                "outcome": clr.loc[donors, cell_type].to_numpy(),
            }))
    return pd.concat(frames, ignore_index=True)


def extended_composition():
    fly = pd.read_csv(RESULTS / "extended-afca-composition-model-input.csv")
    fly["outcome"] = fly["clr"]
    worm = pd.read_csv(RESULTS / "extended-worm-composition-model-input.csv")
    worm["outcome"] = worm["clr"]
    import anndata as ad
    from scripts.extended_atlas_utils import clr_table, composition_counts
    from scripts.run_zebrafish_testis_aging import standard_obs

    atlas = ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad", backed="r")
    fish_obs = standard_obs(atlas, "5-22mo")
    fish_obs = fish_obs[fish_obs["fine_cell_type"] != "unknown"]
    counts = composition_counts(fish_obs, "donor", "fine_cell_type")
    counts = counts.loc[:, counts.sum(axis=0) >= 50]
    metadata = fish_obs.groupby("donor").agg(age=("age", "first"), sex_male=("sex_male", "first"))
    clr = clr_table(counts)
    fish_frames = []
    for cell_type in counts.columns:
        fish_frames.append(pd.DataFrame({
            "species": "Zebrafish", "tissue": "Testis fine cell types", "cohort": "GSE275361",
            "cell_type": cell_type, "donor": counts.index.astype(str),
            "age": metadata.loc[counts.index, "age"].to_numpy(),
            "sex_male": metadata.loc[counts.index, "sex_male"].to_numpy(),
            "outcome": clr[cell_type].to_numpy(),
        }))
    fish = pd.concat(fish_frames, ignore_index=True)
    return pd.concat([fly, worm, fish], ignore_index=True)


def vector_donors(reference_age=6, depth_adjusted=True):
    core = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    core = core[core["species"] != "Naked mole-rat"].copy()
    fly = pd.read_csv(RESULTS / "extended-afca-vector-donors.csv")
    worm = pd.read_csv(RESULTS / f"extended-worm-vector-donors-reference{reference_age}.csv")
    fish = pd.read_csv(RESULTS / "extended-zebrafish-vector-donors-fine.csv")
    data = pd.concat([core, fly, worm, fish], ignore_index=True)
    data["outcome"] = data["profile_hellinger"]
    data["depth_adjusted"] = depth_adjusted
    return data


def unit_effects(data, analysis, sensitivity, depth_adjusted=False):
    data = add_classification(data)
    data["unit"] = (data["species"] + " | " + data["tissue"].astype(str) + " | "
                    + data["cohort"] + " | " + data["cell_type"])
    rows = []
    for unit, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 6 or part["outcome"].std(ddof=1) <= 0:
            continue
        predictors = ["age"]
        for column in ["sex_male", "platform_nextseq"]:
            if column in part and part[column].nunique() > 1:
                predictors.append(column)
        if depth_adjusted:
            predictors.append("mean_log_umi")
        X = sm.add_constant(part[predictors].astype(float), has_constant="add")
        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            continue
        fit = sm.OLS(part["outcome"].astype(float), X).fit()
        outcome_sd = part["outcome"].std(ddof=1)
        species = part["species"].iloc[0]
        scale = AGE_SCALE[species]
        rows.append({
            "species": species, "tissue": part["tissue"].iloc[0], "cohort": part["cohort"].iloc[0],
            "cell_type": part["cell_type"].iloc[0], "unit": unit,
            "proliferation_class": part["proliferation_class"].iloc[0],
            "classification_basis": part["classification_basis"].iloc[0],
            "analysis": analysis, "sensitivity": sensitivity, "n_donors": len(part),
            "n_cells": int(part["n_cells"].sum()) if "n_cells" in part else np.nan,
            "age_slope": fit.params["age"], "age_slope_se": fit.bse["age"],
            "outcome_sd": outcome_sd, "age_scale": scale,
            "effect_sd": fit.params["age"] * scale / outcome_sd,
            "se_sd": fit.bse["age"] * scale / outcome_sd,
            "pvalue": fit.pvalues["age"],
        })
    return pd.DataFrame(rows)


def design_matrix(part, response):
    selected = part[part["proliferation_class"].isin([PROLIFERATIVE, POST_MITOTIC])].copy()
    if selected.empty:
        return None
    selected["post_mitotic"] = (selected["proliferation_class"] == POST_MITOTIC).astype(float)
    selected["stratum"] = (selected["species"] + " | " + selected["tissue"].astype(str)
                           + " | " + selected["cohort"])
    strata = selected.groupby("stratum")["post_mitotic"].nunique()
    selected = selected[selected["stratum"].isin(strata[strata == 2].index)].copy()
    if selected.empty or selected["post_mitotic"].nunique() < 2:
        return None
    dummies = pd.get_dummies(selected["stratum"], drop_first=True, dtype=float)
    X = sm.add_constant(pd.concat([selected[["post_mitotic"]], dummies], axis=1), has_constant="add")
    return selected, selected[response].astype(float), X


def moderator_test(part, response, analysis, sensitivity, scope):
    designed = design_matrix(part, response)
    if designed is None:
        return None
    selected, y, X = designed
    fit = sm.OLS(y, X).fit(cov_type="HC3")
    ci = fit.conf_int().loc["post_mitotic"]
    return {
        "analysis": analysis, "sensitivity": sensitivity, "scope": scope,
        "response": response, "n_units": len(selected),
        "n_proliferative": int((selected["post_mitotic"] == 0).sum()),
        "n_post_mitotic": int((selected["post_mitotic"] == 1).sum()),
        "n_strata": selected["stratum"].nunique(),
        "post_minus_proliferative": fit.params["post_mitotic"],
        "se": fit.bse["post_mitotic"], "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1],
        "pvalue": fit.pvalues["post_mitotic"],
    }


def permutation_pvalue(part, response, permutations=2000, seed=20261104):
    designed = design_matrix(part, response)
    if designed is None:
        return np.nan
    selected, y, X = designed
    observed = float(sm.OLS(y, X).fit().params["post_mitotic"])
    rng = np.random.default_rng(seed + sum(map(ord, response)))
    labels = selected["post_mitotic"].to_numpy().copy()
    stratum_values = selected["stratum"].to_numpy()
    x_values = X.to_numpy()
    post_column = list(X.columns).index("post_mitotic")
    y_values = y.to_numpy()
    exceed = 0
    for _ in range(permutations):
        shuffled = labels.copy()
        for stratum in np.unique(stratum_values):
            indices = np.flatnonzero(stratum_values == stratum)
            shuffled[indices] = rng.permutation(shuffled[indices])
        x_values[:, post_column] = shuffled
        coefficient = np.linalg.lstsq(x_values, y_values, rcond=None)[0][post_column]
        exceed += abs(coefficient) >= abs(observed) - 1e-12
    return (1 + exceed) / (permutations + 1)


def group_summary(units):
    selected = units[units["proliferation_class"].isin([PROLIFERATIVE, POST_MITOTIC])].copy()
    selected["absolute_effect_sd"] = selected["effect_sd"].abs()
    rows = []
    keys = ["analysis", "sensitivity", "species", "proliferation_class"]
    for key, part in selected.groupby(keys):
        rows.append(dict(zip(keys, key)) | {
            "n_units": len(part), "mean_effect_sd": part["effect_sd"].mean(),
            "median_effect_sd": part["effect_sd"].median(),
            "mean_absolute_effect_sd": part["absolute_effect_sd"].mean(),
            "median_absolute_effect_sd": part["absolute_effect_sd"].median(),
            "n_p_lt_0.05": int((part["pvalue"] < 0.05).sum()),
        })
    return pd.DataFrame(rows)


def run_analysis(data, analysis, sensitivity, depth_adjusted, permutations):
    units = unit_effects(data, analysis, sensitivity, depth_adjusted)
    rows = []
    scopes = [(species, part) for species, part in units.groupby("species")]
    mammals = units[units["species"].isin(["Mouse", "Human"])]
    scopes.append(("Mammals combined", mammals))
    if sensitivity != "observed_records_only":
        scopes.append(("All species combined", units))
    for scope, part in scopes:
        for response in ["effect_sd", "absolute_effect_sd"]:
            part = part.copy()
            part["absolute_effect_sd"] = part["effect_sd"].abs()
            row = moderator_test(part, response, analysis, sensitivity, scope)
            if row:
                row["permutation_p"] = permutation_pvalue(part, response, permutations=permutations)
                rows.append(row)
    return units, pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=2000)
    parser.add_argument("--worm-reference-age", type=int, default=6)
    args = parser.parse_args()
    all_units = []
    all_tests = []
    composition = pd.concat(
        [core_composition("all_donors"), extended_composition()], ignore_index=True
    )
    specifications = [
        (vector_donors(args.worm_reference_age, True), "expression_vector", "depth_adjusted", True),
        (vector_donors(args.worm_reference_age, False), "expression_vector", "unadjusted", False),
        (composition, "cell_fraction", "all_donors", False),
        (core_composition("observed_records_only"), "cell_fraction", "observed_records_only", False),
    ]
    for data, analysis, sensitivity, depth in specifications:
        units, tests = run_analysis(data, analysis, sensitivity, depth, args.permutations)
        all_units.append(units)
        all_tests.append(tests)
    from statsmodels.stats.multitest import multipletests

    units = pd.concat(all_units, ignore_index=True)
    tests = pd.concat(all_tests, ignore_index=True)
    tests["q_bh_across_species"] = np.nan
    species_rows = ~tests["scope"].str.endswith("combined")
    for _, indices in tests[species_rows].groupby(["analysis", "sensitivity", "response"]).groups.items():
        tests.loc[indices, "q_bh_across_species"] = multipletests(
            tests.loc[indices, "pvalue"], method="fdr_bh"
        )[1]
    units.to_csv(RESULTS / "proliferation-unit-effects.csv", index=False)
    tests.to_csv(RESULTS / "proliferation-moderator-tests.csv", index=False)
    group_summary(units).to_csv(RESULTS / "proliferation-group-summary.csv", index=False)
    curation = units[["species", "cell_type", "proliferation_class", "classification_basis"]].drop_duplicates()
    curation.sort_values(["species", "cell_type"]).to_csv(
        RESULTS / "proliferation-cell-type-classification.csv", index=False
    )
    print(tests.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
