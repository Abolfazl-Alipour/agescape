#!/usr/bin/env python3
"""Stable global Bayesian models for human expression-vector displacement."""

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

from scripts.run_cross_species_bayesian_vector import ARTIFACTS, RESULTS, SEED, prepare_data
from scripts.run_h3_bayesian_disorder import run_diagnostics, save_idata


def build_global_model(coords, data, prior_sd, depth_adjusted):
    with pm.Model(coords=coords) as model:
        age = pm.Data("age", data["age"], dims="obs")
        sex = pm.Data("sex", data["sex"], dims="obs")
        depth = pm.Data("depth", data["depth"], dims="obs")
        unit_index = pm.Data("unit_index", data["unit_index"], dims="obs")
        donor_index = pm.Data("donor_index", data["donor_index"], dims="obs")
        global_intercept = pm.Normal("global_intercept", 0, 0.5)
        global_age = pm.Normal("global_age", 0, prior_sd)
        global_sex = pm.Normal("global_sex", 0, 0.5)
        tau_donor = pm.Gamma("tau_donor", 2, 2)
        z_donor = pm.Normal("z_donor", 0, 1, dims="donor")
        donor_intercept = pm.Deterministic("donor_intercept", tau_donor * z_donor, dims="donor")
        sigma = pm.Gamma("sigma", 2, 2, dims="unit")
        nu = pm.Gamma("nu", 2, 0.1)
        mu = (
            global_intercept
            + global_age * age
            + global_sex * sex
            + donor_intercept[donor_index]
        )
        if depth_adjusted:
            global_depth = pm.Normal("global_depth", 0, 0.5)
            mu = mu + global_depth * depth
        pm.StudentT("y_obs", nu=nu, mu=mu, sigma=sigma[unit_index], observed=data["y"], dims="obs")
    return model


def fit_global(donors, depth_adjusted, prior_sd):
    model_name = "depth_adjusted" if depth_adjusted else "unadjusted"
    label = f"human-global-within-unit-{model_name}-prior-{prior_sd:.2f}"
    slug = ARTIFACTS / label
    slug.mkdir(parents=True, exist_ok=True)
    part, coords, data = prepare_data(donors, "Human", depth_adjusted)
    data["age"] = (
        (part["age"] - part.groupby("unit")["age"].transform("mean")) / 12.25
    ).to_numpy()
    rng = np.random.default_rng(SEED + sum(map(ord, label)))
    model = build_global_model(coords, data, prior_sd, depth_adjusted)
    with model:
        prior = pm.sample_prior_predictive(draws=500, random_seed=rng)
        idata = pm.sample(
            draws=4000, tune=4000, nuts_sampler="nutpie", target_accept=0.99,
            init="adapt_diag_grad", random_seed=rng, progressbar=False,
        )
        save_idata(idata, slug / "inference_data.nc")
        idata.extend(prior)
        idata.extend(pm.sample_posterior_predictive(idata, random_seed=rng, progressbar=False))
        pm.compute_log_likelihood(idata, model=model)
        pm.compute_log_prior(idata, model=model)
        save_idata(idata, slug / "inference_data.nc")
    az.summary(idata, hdi_prob=0.94).to_csv(slug / "summary.csv")
    run_diagnostics(slug)
    hdi = az.hdi(idata, hdi_prob=0.94)["global_age"]
    values = idata.posterior["global_age"].to_numpy().ravel()
    result = {
        "species": "Human", "model": f"global_{model_name}", "prior_sd": prior_sd,
        "scope": "species", "unit": "ALL", "n_donors": part["donor_cluster"].nunique(),
        "posterior_mean": float(values.mean()), "posterior_sd": float(values.std()),
        "hdi94_low": float(hdi.sel(hdi="lower")), "hdi94_high": float(hdi.sel(hdi="higher")),
        "prob_positive": float(np.mean(values > 0)),
    }
    pd.DataFrame([result]).to_csv(slug / "posterior_effects.csv", index=False)
    diagnostics = {"species": "Human", "model": f"global_{model_name}", "prior_sd": prior_sd}
    diagnostics.update(json.loads((slug / "diagnostics.json").read_text()))
    prior_values = prior.prior_predictive["y_obs"].to_numpy().ravel()
    prior_summary = {
        "species": "Human", "model": f"global_{model_name}", "prior_sd": prior_sd,
        "n_obs": len(data["y"]), "n_units": len(coords["unit"]),
        "observed_min": float(data["y"].min()), "observed_max": float(data["y"].max()),
        "prior_q001": float(np.quantile(prior_values, 0.001)),
        "prior_q999": float(np.quantile(prior_values, 0.999)),
    }
    return result, diagnostics, prior_summary


def main():
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    posterior_rows = []
    diagnostic_rows = []
    prior_rows = []
    for depth_adjusted, prior_sd in ((False, 0.5), (True, 0.5), (True, 0.25), (True, 1.0)):
        print(f"Global human depth={depth_adjusted} prior={prior_sd}", flush=True)
        posterior, diagnostics, prior = fit_global(donors, depth_adjusted, prior_sd)
        posterior_rows.append(posterior)
        diagnostic_rows.append(diagnostics)
        prior_rows.append(prior)
    pd.DataFrame(posterior_rows).to_csv(RESULTS / "h3-human-bayesian-global-posterior.csv", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(RESULTS / "h3-human-bayesian-global-diagnostics.csv", index=False)
    pd.DataFrame(prior_rows).to_csv(RESULTS / "h3-human-bayesian-global-prior-predictive.csv", index=False)
    print(pd.DataFrame(posterior_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
