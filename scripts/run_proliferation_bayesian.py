#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

RESULTS = REPO / "results"
SEED = sum(map(ord, "agescape-proliferation-moderator-bayes-v1"))


def prepare(part):
    selected = part[part["proliferation_class"].isin(["proliferative_capable", "post_mitotic"])].copy()
    if selected.empty:
        return None
    selected["post_mitotic"] = (selected["proliferation_class"] == "post_mitotic").astype(float)
    selected["stratum"] = (selected["species"] + " | " + selected["tissue"].astype(str)
                           + " | " + selected["cohort"])
    strata = selected.groupby("stratum")["post_mitotic"].nunique()
    selected = selected[selected["stratum"].isin(strata[strata == 2].index)].copy()
    if selected.empty or selected["post_mitotic"].nunique() < 2:
        return None
    strata = sorted(selected["stratum"].unique())
    stratum_map = {value: index for index, value in enumerate(strata)}
    coords = {"stratum": strata, "unit": np.arange(len(selected))}
    data = {
        "effect": selected["effect_sd"].to_numpy(),
        "se": np.maximum(selected["se_sd"].to_numpy(), 1e-8),
        "post": selected["post_mitotic"].to_numpy(),
        "stratum_index": selected["stratum"].map(stratum_map).to_numpy(),
    }
    return selected, coords, data


def fit_scope(part, analysis, sensitivity, scope, prior_sd, draws):
    prepared = prepare(part)
    if prepared is None:
        return None
    selected, coords, data = prepared
    with pm.Model(coords=coords) as model:
        stratum_index = pm.Data("stratum_index", data["stratum_index"], dims="unit")
        post = pm.Data("post", data["post"], dims="unit")
        stratum_effect = pm.Normal("stratum_effect", 0, 1, dims="stratum")
        contrast = pm.Normal("post_minus_proliferative", 0, prior_sd)
        tau = pm.HalfNormal("tau", 0.5)
        z = pm.Normal("z", 0, 1, dims="unit")
        mu = stratum_effect[stratum_index] + contrast * post + tau * z
        pm.Normal("effect_obs", mu=mu, sigma=data["se"], observed=data["effect"], dims="unit")
        idata = pm.sample(
            draws=draws, tune=draws, chains=4, nuts_sampler="nutpie", target_accept=0.95,
            random_seed=np.random.default_rng(SEED + sum(map(ord, scope + analysis + str(prior_sd)))),
            progressbar=False,
        )
    values = idata.posterior["post_minus_proliferative"].to_numpy().ravel()
    hdi = az.hdi(idata, hdi_prob=0.94)["post_minus_proliferative"].to_numpy()
    summary = az.summary(idata, var_names=["post_minus_proliferative", "tau"], hdi_prob=0.94)
    return {
        "analysis": analysis, "sensitivity": sensitivity, "scope": scope, "prior_sd": prior_sd,
        "n_units": len(selected), "n_proliferative": int((selected["post_mitotic"] == 0).sum()),
        "n_post_mitotic": int((selected["post_mitotic"] == 1).sum()),
        "n_strata": selected["stratum"].nunique(), "posterior_mean": values.mean(),
        "posterior_sd": values.std(), "hdi94_low": hdi[0], "hdi94_high": hdi[1],
        "prob_post_greater": np.mean(values > 0),
        "contrast_rhat": summary.loc["post_minus_proliferative", "r_hat"],
        "contrast_ess_bulk": summary.loc["post_minus_proliferative", "ess_bulk"],
        "tau_rhat": summary.loc["tau", "r_hat"],
        "n_divergences": int(idata.sample_stats["diverging"].sum()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--skip-prior-sensitivity", action="store_true")
    args = parser.parse_args()
    units = pd.read_csv(RESULTS / "proliferation-unit-effects.csv")
    specifications = [
        ("expression_vector", "depth_adjusted"),
        ("expression_vector", "unadjusted"),
        ("cell_fraction", "all_donors"),
        ("cell_fraction", "observed_records_only"),
    ]
    priors = [0.5] if args.skip_prior_sensitivity else [0.25, 0.5, 1.0]
    rows = []
    for analysis, sensitivity in specifications:
        selected = units[(units["analysis"] == analysis) & (units["sensitivity"] == sensitivity)]
        scopes = [(species, part) for species, part in selected.groupby("species")]
        mammals = selected[selected["species"].isin(["Mouse", "Human"])]
        scopes.append(("Mammals combined", mammals))
        if sensitivity != "observed_records_only":
            scopes.append(("All species combined", selected))
        for scope, part in scopes:
            for prior_sd in priors:
                print(f"{analysis}, {sensitivity}, {scope}, prior={prior_sd}", flush=True)
                row = fit_scope(part, analysis, sensitivity, scope, prior_sd, args.draws)
                if row:
                    rows.append(row)
    posterior = pd.DataFrame(rows)
    posterior.to_csv(RESULTS / "proliferation-bayesian-posterior.csv", index=False)
    print(posterior[posterior["prior_sd"] == 0.5].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
