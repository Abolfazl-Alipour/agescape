#!/usr/bin/env python3
"""Bayesian hierarchical analysis of lifespan-scaled expression-vector displacement."""

import argparse
import json
from pathlib import Path
import sys

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_h3_bayesian_disorder import run_diagnostics, save_idata

RESULTS = REPO / "results"
ARTIFACTS = RESULTS / "h3-cross-species-bayesian"
SEED = sum(map(ord, "agescape-cross-species-expression-vector-bayes-v1"))


def prepare_data(donors, species, depth_adjusted):
    part = donors[donors["species"] == species].copy()
    part["unit"] = part["tissue"] + " | " + part["cohort"] + " | " + part["cell_type"]
    part["group"] = part["cohort"] if species == "Human" else part["tissue"]
    part["donor_cluster"] = part["cohort"] + " | " + part["donor"]
    minimum = 5 if species == "Naked mole-rat" else 6
    part = part.groupby("unit").filter(lambda frame: len(frame) >= minimum)
    units = sorted(part["unit"].unique())
    groups = sorted(part["group"].unique())
    donor_levels = sorted(part["donor_cluster"].unique())
    unit_map = {value: index for index, value in enumerate(units)}
    group_map = {value: index for index, value in enumerate(groups)}
    donor_map = {value: index for index, value in enumerate(donor_levels)}
    age_scale = {"Mouse": 0.4, "Human": 12.25, "Naked mole-rat": 3.1}[species]
    reference_age = {"Mouse": 0.0, "Human": 40.0, "Naked mole-rat": 0.0}[species]
    part["age_scaled"] = (part["age"] - reference_age) / age_scale
    part["sex_centered"] = part["sex_male"] - part.groupby("unit")["sex_male"].transform("mean")
    part["depth_centered"] = part["mean_log_umi"] - part.groupby("unit")["mean_log_umi"].transform("mean")
    means = part.groupby("unit")["profile_hellinger"].transform("mean")
    sds = part.groupby("unit")["profile_hellinger"].transform("std")
    part["y"] = (part["profile_hellinger"] - means) / sds
    part = part.sort_values(["unit", "donor_cluster"]).reset_index(drop=True)
    coords = {"unit": units, "group": groups, "donor": donor_levels, "obs": part.index.to_numpy()}
    unit_groups = part.groupby("unit")["group"].first().loc[units]
    data = {
        "age": part["age_scaled"].to_numpy(),
        "sex": part["sex_centered"].to_numpy(),
        "depth": part["depth_centered"].to_numpy(),
        "unit_index": part["unit"].map(unit_map).to_numpy(),
        "unit_group_index": unit_groups.map(group_map).to_numpy(),
        "donor_index": part["donor_cluster"].map(donor_map).to_numpy(),
        "y": part["y"].to_numpy(),
    }
    return part, coords, data


def build_model(coords, data, prior_sd, depth_adjusted, fixed_group_effects=False):
    with pm.Model(coords=coords) as model:
        age = pm.Data("age", data["age"], dims="obs")
        sex = pm.Data("sex", data["sex"], dims="obs")
        depth = pm.Data("depth", data["depth"], dims="obs")
        unit_index = pm.Data("unit_index", data["unit_index"], dims="obs")
        unit_group_index = pm.Data("unit_group_index", data["unit_group_index"], dims="unit")
        donor_index = pm.Data("donor_index", data["donor_index"], dims="obs")
        global_intercept = pm.Normal("global_intercept", 0, 0.5)
        global_age = pm.Normal("global_age", 0, prior_sd)
        global_sex = pm.Normal("global_sex", 0, 0.5)
        tau_intercept = pm.Gamma("tau_intercept", 2, 2)
        tau_unit_age = pm.Gamma("tau_unit_age", 2, 2)
        tau_donor = pm.Gamma("tau_donor", 2, 2)
        z_intercept = pm.Normal("z_intercept", 0, 1, dims="unit")
        z_unit_age = pm.Normal("z_unit_age", 0, 1, dims="unit")
        z_donor = pm.Normal("z_donor", 0, 1, dims="donor")
        unit_intercept = pm.Deterministic(
            "unit_intercept", global_intercept + tau_intercept * z_intercept, dims="unit"
        )
        if len(coords["group"]) > 1:
            if fixed_group_effects:
                group_offset = pm.Normal("group_offset", 0, 0.25, dims="group")
                group_age = pm.Deterministic("group_age", global_age + group_offset, dims="group")
            else:
                tau_group_age = pm.Gamma("tau_group_age", 2, 2)
                z_group_age = pm.Normal("z_group_age", 0, 1, dims="group")
                group_age = pm.Deterministic(
                    "group_age", global_age + tau_group_age * z_group_age, dims="group"
                )
            unit_age = pm.Deterministic(
                "unit_age", group_age[unit_group_index] + tau_unit_age * z_unit_age, dims="unit"
            )
        else:
            unit_age = pm.Deterministic(
                "unit_age", global_age + tau_unit_age * z_unit_age, dims="unit"
            )
        donor_intercept = pm.Deterministic("donor_intercept", tau_donor * z_donor, dims="donor")
        sigma = pm.Gamma("sigma", 2, 2, dims="unit")
        nu = pm.Gamma("nu", 2, 0.1)
        mu = (
            unit_intercept[unit_index]
            + unit_age[unit_index] * age
            + global_sex * sex
            + donor_intercept[donor_index]
        )
        if depth_adjusted:
            global_depth = pm.Normal("global_depth", 0, 0.5)
            mu = mu + global_depth * depth
        pm.StudentT(
            "y_obs", nu=nu, mu=mu, sigma=sigma[unit_index], observed=data["y"], dims="obs"
        )
    return model


def posterior_rows(idata, part, coords, species, model_name, prior_sd):
    hdi = az.hdi(idata, hdi_prob=0.94)
    posterior = idata.posterior
    rows = []

    def row(scope, label, values, hdi_values=None):
        return {
            "species": species, "model": model_name, "prior_sd": prior_sd,
            "scope": scope, "unit": label, "n_donors": part["donor_cluster"].nunique(),
            "posterior_mean": float(values.mean()), "posterior_sd": float(values.std()),
            "hdi94_low": float(hdi_values.sel(hdi="lower")) if hdi_values is not None else np.nan,
            "hdi94_high": float(hdi_values.sel(hdi="higher")) if hdi_values is not None else np.nan,
            "prob_positive": float(np.mean(values > 0)),
        }

    global_values = posterior["global_age"].to_numpy().ravel()
    rows.append(row("species", "ALL", global_values, hdi["global_age"]))
    if "group_age" in posterior:
        for group in coords["group"]:
            values = posterior["group_age"].sel(group=group).to_numpy().ravel()
            rows.append(row("tissue_or_cohort", group, values, hdi["group_age"].sel(group=group)))
    for unit in coords["unit"]:
        values = posterior["unit_age"].sel(unit=unit).to_numpy().ravel()
        rows.append(row("cell_type_unit", unit, values, hdi["unit_age"].sel(unit=unit)))
    return rows


def fit_model(donors, species, depth_adjusted, prior_sd, diagnostics=True, draws=2000,
              fixed_group_effects=False):
    model_name = "depth_adjusted" if depth_adjusted else "unadjusted"
    if fixed_group_effects:
        model_name += "_fixed_group"
    label = f"{species.lower().replace(' ', '-')}-{model_name}-prior-{prior_sd:.2f}"
    slug = ARTIFACTS / label
    slug.mkdir(parents=True, exist_ok=True)
    part, coords, data = prepare_data(donors, species, depth_adjusted)
    rng = np.random.default_rng(SEED + sum(map(ord, label)))
    model = build_model(coords, data, prior_sd, depth_adjusted, fixed_group_effects)
    with model:
        prior = pm.sample_prior_predictive(draws=500, random_seed=rng)
    prior_values = prior.prior_predictive["y_obs"].to_numpy().ravel()
    observed = data["y"]
    prior_summary = {
        "species": species, "model": model_name, "prior_sd": prior_sd,
        "n_obs": len(observed), "n_units": len(coords["unit"]),
        "observed_min": float(observed.min()), "observed_max": float(observed.max()),
        "prior_q001": float(np.quantile(prior_values, 0.001)),
        "prior_q999": float(np.quantile(prior_values, 0.999)),
        "observed_within_prior_998": bool(
            observed.min() >= np.quantile(prior_values, 0.001)
            and observed.max() <= np.quantile(prior_values, 0.999)
        ),
    }
    with model:
        idata = pm.sample(
            draws=draws, tune=draws, nuts_sampler="nutpie", target_accept=0.99,
            init="adapt_diag_grad", random_seed=rng, progressbar=False,
        )
        save_idata(idata, slug / "inference_data.nc")
        idata.extend(prior)
        idata.extend(pm.sample_posterior_predictive(idata, random_seed=rng, progressbar=False))
        pm.compute_log_likelihood(idata, model=model)
        pm.compute_log_prior(idata, model=model)
        save_idata(idata, slug / "inference_data.nc")
    az.summary(idata, hdi_prob=0.94).to_csv(slug / "summary.csv")
    if diagnostics:
        run_diagnostics(slug)
    rows = posterior_rows(idata, part, coords, species, model_name, prior_sd)
    pd.DataFrame(rows).to_csv(slug / "posterior_effects.csv", index=False)
    diagnostic_summary = {"species": species, "model": model_name, "prior_sd": prior_sd}
    diagnostics_path = slug / "diagnostics.json"
    if diagnostics_path.exists():
        diagnostic_summary.update(json.loads(diagnostics_path.read_text()))
    return rows, prior_summary, diagnostic_summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-diagnostics", action="store_true")
    parser.add_argument("--primary-only", action="store_true")
    args = parser.parse_args()
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    specifications = [
        ("Mouse", False, 0.5), ("Mouse", True, 0.5),
        ("Human", False, 0.5),
        ("Naked mole-rat", False, 0.5),
    ]
    all_rows = []
    prior_rows = []
    diagnostic_rows = []
    for species, depth_adjusted, prior_sd in specifications:
        print(f"Fitting {species}, depth_adjusted={depth_adjusted}, prior={prior_sd}", flush=True)
        rows, prior, diagnostics = fit_model(
            donors, species, depth_adjusted, prior_sd, diagnostics=not args.skip_diagnostics,
            draws=4000 if species == "Human" else 2000,
            fixed_group_effects=species == "Human",
        )
        all_rows.extend(rows)
        prior_rows.append(prior)
        diagnostic_rows.append(diagnostics)
    posterior = pd.DataFrame(all_rows)
    posterior.to_csv(RESULTS / "h3-cross-species-bayesian-posterior.csv", index=False)
    pd.DataFrame(prior_rows).to_csv(RESULTS / "h3-cross-species-bayesian-prior-predictive.csv", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(RESULTS / "h3-cross-species-bayesian-diagnostics.csv", index=False)
    print(posterior[posterior["scope"] == "species"].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
