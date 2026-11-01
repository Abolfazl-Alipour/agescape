#!/usr/bin/env python3
"""Hierarchical Bayesian analysis of age-related transcriptional disorder."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import arviz as az
import arviz_plots as azp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
ARTIFACTS = RESULTS / "h3-bayesian"
WORKFLOW = Path("/tmp/baygent-skills/bayesian-workflow")
WORKFLOW_COMMIT = "18ab0de9b63b5bb03caa96146235e8f1534d45d4"
RANDOM_SEED = sum(map(ord, "agescape-transcriptional-disorder-hierarchical-v1"))


def ensure_workflow():
    if (WORKFLOW / "scripts" / "diagnose_model.py").exists():
        return
    root = WORKFLOW.parent
    subprocess.run([
        "git", "clone", "https://github.com/Learning-Bayesian-Statistics/baygent-skills.git",
        str(root),
    ], check=True)
    subprocess.run(["git", "-C", str(root), "checkout", WORKFLOW_COMMIT], check=True)


def save_idata(idata, path):
    idata._attrs = {
        key: json.dumps(value) if isinstance(value, dict) else value
        for key, value in idata._attrs.items()
    }
    for group in idata.groups():
        dataset = getattr(idata, group)
        dataset.attrs = {
            key: json.dumps(value) if isinstance(value, dict) else value
            for key, value in dataset.attrs.items()
        }
    idata.to_netcdf(path)


def prepare_data(donors, metric, analysis):
    part = donors[donors["sampling"] == analysis["sampling"]].copy()
    if analysis["donor_set"] == "exclude_pooled_ids":
        pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                     .query("is_pooled_id")["mouse.id"].astype(str))
        part = part[~part["donor"].isin(pooled)]
    part = part.dropna(subset=[metric])
    part = part.groupby("cell_type").filter(lambda frame: len(frame) >= 6)
    cell_types = sorted(part["cell_type"].unique())
    part["cell_type"] = pd.Categorical(part["cell_type"], categories=cell_types, ordered=True)
    part = part.sort_values(["cell_type", "donor"]).reset_index(drop=True)
    part["age15"] = (part["age_months"] - 3) / 15
    part["sex_centered"] = part["sex_male"] - part.groupby("cell_type")["sex_male"].transform("mean")
    means = part.groupby("cell_type", observed=True)[metric].transform("mean")
    sds = part.groupby("cell_type", observed=True)[metric].transform("std")
    part["y"] = (part[metric] - means) / sds
    part["y_sd"] = sds
    donor_levels = sorted(part["donor"].unique())
    part["donor"] = pd.Categorical(part["donor"], categories=donor_levels, ordered=True)
    coords = {
        "cell_type": cell_types, "donor": donor_levels, "obs": part.index.to_numpy()
    }
    model_data = {
        "age15": part["age15"].to_numpy(),
        "sex_centered": part["sex_centered"].to_numpy(),
        "cell_index": part["cell_type"].cat.codes.to_numpy(),
        "donor_index": part["donor"].cat.codes.to_numpy(),
        "y": part["y"].to_numpy(),
    }
    return part, cell_types, coords, model_data


def build_model(coords, model_data, prior_sd, likelihood):
    with pm.Model(coords=coords) as model:
        age15 = pm.Data("age15", model_data["age15"], dims="obs")
        sex_centered = pm.Data("sex_centered", model_data["sex_centered"], dims="obs")
        cell_index = pm.Data("cell_index", model_data["cell_index"], dims="obs")
        donor_index = pm.Data("donor_index", model_data["donor_index"], dims="obs")
        global_intercept = pm.Normal("global_intercept", 0, 0.5)
        global_age = pm.Normal("global_age", 0, prior_sd)
        global_sex = pm.Normal("global_sex", 0, 0.5)
        tau_intercept = pm.Gamma("tau_intercept", alpha=2, beta=2)
        tau_age = pm.Gamma("tau_age", alpha=2, beta=2)
        tau_donor = pm.Gamma("tau_donor", alpha=2, beta=2)
        z_intercept = pm.Normal("z_intercept", 0, 1, dims="cell_type")
        z_age = pm.Normal("z_age", 0, 1, dims="cell_type")
        z_donor = pm.Normal("z_donor", 0, 1, dims="donor")
        ct_intercept = pm.Deterministic(
            "ct_intercept", global_intercept + tau_intercept * z_intercept, dims="cell_type"
        )
        ct_age = pm.Deterministic("ct_age", global_age + tau_age * z_age, dims="cell_type")
        donor_intercept = pm.Deterministic("donor_intercept", tau_donor * z_donor, dims="donor")
        sigma = pm.Gamma("sigma", alpha=2, beta=2, dims="cell_type")
        mu = (
            ct_intercept[cell_index] + ct_age[cell_index] * age15
            + global_sex * sex_centered + donor_intercept[donor_index]
        )
        if likelihood == "student_t":
            nu = pm.Gamma("nu", alpha=2, beta=0.1)
            pm.StudentT("y_obs", nu=nu, mu=mu, sigma=sigma[cell_index], observed=model_data["y"], dims="obs")
        else:
            pm.Normal("y_obs", mu=mu, sigma=sigma[cell_index], observed=model_data["y"], dims="obs")
    return model


def frequentist_cell_type(part, metric, cell_type):
    frame = part[part["cell_type"] == cell_type]
    X = sm.add_constant(frame[["age15", "sex_centered"]], has_constant="add")
    fit = sm.OLS(frame["y"], X).fit()
    ci = fit.conf_int().loc["age15"]
    sd = frame["y_sd"].iloc[0]
    return {
        "frequentist_effect_15mo_sd": fit.params["age15"],
        "frequentist_ci95_low_sd": ci.iloc[0],
        "frequentist_ci95_high_sd": ci.iloc[1],
        "frequentist_pvalue": fit.pvalues["age15"],
        "frequentist_slope_per_month": fit.params["age15"] * sd / 15,
    }


def frequentist_global(donors, metric, analysis):
    part, cell_types, coords, model_data = prepare_data(donors, metric, analysis)
    dummies = pd.get_dummies(part["cell_type"], dtype=float).reset_index(drop=True)
    X = pd.concat([dummies, part[["age15", "sex_centered"]].reset_index(drop=True)], axis=1)
    fit = sm.OLS(part["y"].to_numpy(), X).fit(
        cov_type="cluster", cov_kwds={"groups": part["donor"]}
    )
    ci = fit.conf_int().loc["age15"]
    return {
        "metric": metric,
        "analysis": analysis["name"],
        "n_observations": len(part),
        "n_donors": part["donor"].nunique(),
        "n_cell_types": len(cell_types),
        "effect_15mo_sd": fit.params["age15"],
        "ci95_low_sd15": ci.iloc[0],
        "ci95_high_sd15": ci.iloc[1],
        "pvalue": fit.pvalues["age15"],
    }


def run_diagnostics(slug):
    if not (WORKFLOW / "scripts" / "diagnose_model.py").exists():
        return
    subprocess.run([
        sys.executable, str(WORKFLOW / "scripts" / "diagnose_model.py"),
        "--idata", str(slug / "inference_data.nc"),
        "--output", str(slug / "diagnostics.json"),
    ], check=True)
    subprocess.run([
        sys.executable, str(WORKFLOW / "scripts" / "calibration_check.py"),
        "--idata", str(slug / "inference_data.nc"),
        "--output", str(slug / "calibration.json"),
        "--save-plots", "--plot-dir", str(slug),
    ], check=True)
    subprocess.run([
        sys.executable, str(WORKFLOW / "scripts" / "check_diagnostics.py"),
        "--diagnostics", str(slug / "diagnostics.json"),
        "--calibration", str(slug / "calibration.json"),
        "--output", str(slug / "check_report.json"),
    ], check=True)


def fit_model(donors, metric, analysis, prior_sd, likelihood, diagnostics=True):
    name = analysis["name"]
    slug = ARTIFACTS / f"{metric}-{name}-prior-{prior_sd:.2f}-{likelihood}"
    slug.mkdir(parents=True, exist_ok=True)
    part, cell_types, coords, model_data = prepare_data(donors, metric, analysis)
    rng = np.random.default_rng(RANDOM_SEED + sum(map(ord, f"{metric}-{name}-{prior_sd}-{likelihood}")))
    model = build_model(coords, model_data, prior_sd, likelihood)
    with model:
        prior = pm.sample_prior_predictive(draws=500, random_seed=rng)
    prior_values = prior.prior_predictive["y_obs"].to_numpy().ravel()
    observed = model_data["y"]
    prior_summary = {
        "metric": metric,
        "analysis": name,
        "prior_sd": prior_sd,
        "likelihood": likelihood,
        "n_obs": len(observed),
        "n_cell_types": len(cell_types),
        "observed_min": float(observed.min()),
        "observed_max": float(observed.max()),
        "prior_q001": float(np.quantile(prior_values, 0.001)),
        "prior_q999": float(np.quantile(prior_values, 0.999)),
        "prior_abs_gt5_fraction": float(np.mean(np.abs(prior_values) > 5)),
        "observed_within_prior_998": bool(
            observed.min() >= np.quantile(prior_values, 0.001)
            and observed.max() <= np.quantile(prior_values, 0.999)
        ),
    }
    pd.DataFrame([prior_summary]).to_csv(slug / "prior_predictive_summary.csv", index=False)
    graph = pm.model_to_graphviz(model)
    graph.render(directory=slug, filename="model_graph", format="png", cleanup=True)

    with model:
        idata = pm.sample(
            draws=2000, tune=2000, nuts_sampler="nutpie", target_accept=0.99,
            init="adapt_diag_grad", random_seed=rng, progressbar=False,
        )
        save_idata(idata, slug / "inference_data.nc")
        idata.extend(prior)
        idata.extend(pm.sample_posterior_predictive(idata, random_seed=rng, progressbar=False))
        pm.compute_log_likelihood(idata, model=model)
        pm.compute_log_prior(idata, model=model)
        save_idata(idata, slug / "inference_data.nc")

    az.plot_trace(idata, var_names=["global_age", "global_sex", "tau_age", "sigma"])
    plt.gcf().savefig(slug / "trace.png", dpi=130, bbox_inches="tight")
    plt.close()
    az.plot_forest(idata, var_names=["global_age", "ct_age"], combined=True, hdi_prob=0.94)
    plt.gcf().savefig(slug / "forest.png", dpi=130, bbox_inches="tight")
    plt.close()
    azp.plot_ppc_dist(idata).savefig(slug / "posterior_predictive.png")
    azp.plot_ppc_dist(prior, group="prior_predictive").savefig(slug / "prior_predictive.png")
    plt.close("all")
    az.summary(idata, hdi_prob=0.94).to_csv(slug / "summary.csv")
    if diagnostics:
        run_diagnostics(slug)

    hdi = az.hdi(idata, hdi_prob=0.94)
    posterior = idata.posterior
    rows = []
    global_values = posterior["global_age"].to_numpy().ravel()
    rows.append({
        "metric": metric, "analysis": name, "prior_sd": prior_sd, "likelihood": likelihood,
        "scope": "global", "cell_type": "ALL", "n_donors": len(part),
        "posterior_mean_sd15": global_values.mean(), "posterior_sd": global_values.std(),
        "hdi94_low_sd15": float(hdi["global_age"].sel(hdi="lower")),
        "hdi94_high_sd15": float(hdi["global_age"].sel(hdi="higher")),
        "prob_positive": float(np.mean(global_values > 0)),
    })
    for index, cell_type in enumerate(cell_types):
        values = posterior["ct_age"].sel(cell_type=cell_type).to_numpy().ravel()
        local_hdi = hdi["ct_age"].sel(cell_type=cell_type)
        sd = part.loc[part["cell_type"] == cell_type, "y_sd"].iloc[0]
        row = {
            "metric": metric, "analysis": name, "prior_sd": prior_sd, "likelihood": likelihood,
            "scope": "cell_type", "cell_type": cell_type,
            "n_donors": int((part["cell_type"] == cell_type).sum()),
            "posterior_mean_sd15": values.mean(), "posterior_sd": values.std(),
            "hdi94_low_sd15": float(local_hdi.sel(hdi="lower")),
            "hdi94_high_sd15": float(local_hdi.sel(hdi="higher")),
            "prob_positive": float(np.mean(values > 0)),
            "posterior_slope_per_month": values.mean() * sd / 15,
            "hdi94_low_per_month": float(local_hdi.sel(hdi="lower")) * sd / 15,
            "hdi94_high_per_month": float(local_hdi.sel(hdi="higher")) * sd / 15,
        }
        row.update(frequentist_cell_type(part, metric, cell_type))
        rows.append(row)
    pd.DataFrame(rows).to_csv(slug / "posterior_effects.csv", index=False)
    diagnostic_summary = {
        "metric": metric, "analysis": name, "prior_sd": prior_sd, "likelihood": likelihood
    }
    diagnostics_path = slug / "diagnostics.json"
    if diagnostics_path.exists():
        diagnostic_summary.update(json.loads(diagnostics_path.read_text()))
    pd.DataFrame([diagnostic_summary]).to_csv(slug / "diagnostic_summary.csv", index=False)
    return rows, prior_summary, diagnostic_summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--robust-only", action="store_true")
    parser.add_argument("--normal-likelihood", action="store_true")
    parser.add_argument("--skip-diagnostics", action="store_true")
    args = parser.parse_args()
    if not args.skip_diagnostics:
        ensure_workflow()
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    donors = pd.read_csv(RESULTS / "h3-transcriptional-disorder-donors.csv")
    analyses = [{"name": "all_cells_all_donor_ids", "sampling": "all_cells", "donor_set": "all_donor_ids"}]
    primary_only = args.primary_only or args.robust_only
    if not primary_only:
        analyses.extend([
            {"name": "all_cells_exclude_pooled_ids", "sampling": "all_cells", "donor_set": "exclude_pooled_ids"},
            {"name": "balanced_20_cells_all_donor_ids", "sampling": "balanced_20_cells", "donor_set": "all_donor_ids"},
        ])
    prior_sds = [0.5] if primary_only else [0.25, 0.5, 1.0]
    likelihood = "normal" if args.normal_likelihood else "student_t"
    all_rows = []
    prior_rows = []
    diagnostic_rows = []
    for metric in ("entropy", "offprogram_fraction"):
        for analysis in analyses:
            for prior_sd in prior_sds:
                print(f"Fitting {metric}, {analysis['name']}, prior {prior_sd}, {likelihood}", flush=True)
                rows, prior, diagnostics = fit_model(
                    donors, metric, analysis, prior_sd, likelihood,
                    diagnostics=not args.skip_diagnostics
                )
                all_rows.extend(rows)
                prior_rows.append(prior)
                diagnostic_rows.append(diagnostics)
    posterior = pd.DataFrame(all_rows)
    posterior.to_csv(RESULTS / "h3-bayesian-posterior.csv", index=False)
    pd.DataFrame(prior_rows).to_csv(RESULTS / "h3-bayesian-prior-predictive.csv", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(RESULTS / "h3-bayesian-diagnostics.csv", index=False)
    global_frequentist = [
        frequentist_global(donors, metric, analysis)
        for metric in ("entropy", "offprogram_fraction")
        for analysis in analyses
    ]
    pd.DataFrame(global_frequentist).to_csv(RESULTS / "h3-frequentist-global.csv", index=False)
    primary = posterior[
        (posterior["analysis"] == "all_cells_all_donor_ids")
        & (posterior["prior_sd"] == 0.5)
        & (posterior["scope"] == "cell_type")
    ].copy()
    frequentist = pd.read_csv(RESULTS / "h3-transcriptional-disorder.csv")
    frequentist = frequentist[frequentist["analysis"] == "all_cells_all_donor_ids"]
    primary = primary.merge(
        frequentist[["cell_type", "metric", "q_bh"]], on=["cell_type", "metric"], how="left"
    )
    primary.to_csv(RESULTS / "h3-bayesian-frequentist-comparison.csv", index=False)
    print(posterior[(posterior["analysis"] == "all_cells_all_donor_ids")
                    & (posterior["prior_sd"] == 0.5)][[
        "metric", "scope", "cell_type", "posterior_mean_sd15",
        "hdi94_low_sd15", "hdi94_high_sd15", "prob_positive"
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
