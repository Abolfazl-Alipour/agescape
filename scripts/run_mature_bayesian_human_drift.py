#!/usr/bin/env python3
"""Item 5 (human extension): Bayesian hierarchical whole-vector drift on
mature human donors (age >= 30, 30-39 reference).

Mirrors scripts/run_rev3_bayesian_vector_mature_only.py and its shared
model scripts/bayesian_diagnostics_helper.py::build_hierarchical_student_t
(Student-t likelihood, donor intercepts, unit residual scales, tissue-level
age deviations), with human scaling:
  - outcome profile_hellinger standardized WITHIN each tissue|cohort|cell-type
    unit (as in reports/followup-mature-drift.md);
  - age centered at 30 and scaled per 12.25 years (10% of the 122.5-year
    human maximum lifespan), matching the frequentist human slopes;
  - sex (sex_male) and mean log UMI centered within unit;
  - units with >= 6 rows kept (same floor as the mouse refit).
Priors, fixed before fitting: Normal(0, 0.5) on the global age slope
(primary), sensitivity Normal(0, 0.25) and Normal(0, 1.00) on the
depth-adjusted model; Normal(0, 0.5) on sex/depth coefficients;
Gamma(2, 2) on group/unit/donor scales and unit sigma; Gamma(2, 0.1) on nu.
Sampler: NUTS via nutpie, 4 chains, 1000 tuning draws, 2000 kept draws per
chain, target_accept 0.95. Seed fixed; the primary depth-adjusted spec is
rerun with a second seed and the headline move is reported.

Runs on beast in the agescape-bayesian conda env (environment-bayesian.yml).
Inputs (copied to beast): results/followup-mature-drift-human-donors.csv,
scripts/bayesian_diagnostics_helper.py (imported from the same directory).
Outputs: results/followup-mature-bayesian-human-drift.csv (posterior rows),
results/followup-mature-bayesian-human-drift-diagnostics.csv,
results/followup-mature-bayesian-human-drift-frequentist.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from bayesian_diagnostics_helper import (  # noqa: E402
    build_hierarchical_student_t,
    summarize_diagnostics,
)

SEED = sum(map(ord, "agescape-mature-human-drift-bayes-v1"))
SEED2 = SEED + 1
DATA = HERE / "followup-mature-drift-human-donors.csv"
OUT = HERE  # outputs land beside the script on beast; copied back to results/


def prepare_human_data(donors: pd.DataFrame) -> tuple[pd.DataFrame, dict, dict, bool]:
    part = donors.copy()
    part["unit"] = part["tissue"] + " | " + part["cohort"] + " | " + part["cell_type"]
    part["group"] = part["tissue"]
    part["donor_cluster"] = part["cohort"] + " | " + part["donor"].astype(str)
    part = part.groupby("unit").filter(lambda f: len(f) >= 6).copy()
    means = part.groupby("unit")["profile_hellinger"].transform("mean")
    sds = part.groupby("unit")["profile_hellinger"].transform("std")
    part["y"] = (part["profile_hellinger"] - means) / sds
    part = part[np.isfinite(part["y"])].copy()
    part = part.groupby("unit").filter(lambda f: len(f) >= 6).copy()
    # pooled within-unit SD of the raw outcome, for raw-scale reporting
    pooled_sd = float(np.sqrt(np.mean(part.groupby("unit")["profile_hellinger"].var(ddof=1))))
    units = sorted(part["unit"].unique())
    groups = sorted(part["group"].unique())
    donor_levels = sorted(part["donor_cluster"].unique())
    unit_map = {v: i for i, v in enumerate(units)}
    group_map = {v: i for i, v in enumerate(groups)}
    donor_map = {v: i for i, v in enumerate(donor_levels)}
    part["age_scaled"] = (part["age"] - 30.0) / 12.25
    part["sex_centered"] = part["sex_male"] - part.groupby("unit")["sex_male"].transform("mean")
    part["depth_centered"] = part["mean_log_umi"] - part.groupby("unit")["mean_log_umi"].transform("mean")
    part = part.sort_values(["unit", "donor_cluster"]).reset_index(drop=True)
    coords = {"unit": units, "group": groups, "donor": donor_levels,
              "obs": part.index.to_numpy()}
    unit_groups = part.groupby("unit")["group"].first().loc[units]
    data = {
        "age": part["age_scaled"].to_numpy(dtype=float),
        "depth": part["depth_centered"].to_numpy(dtype=float),
        "unit_index": part["unit"].map(unit_map).to_numpy(dtype=int),
        "unit_group_index": unit_groups.map(group_map).to_numpy(dtype=int),
        "donor_index": part["donor_cluster"].map(donor_map).to_numpy(dtype=int),
        "y": part["y"].to_numpy(dtype=float),
        "sex": part["sex_centered"].to_numpy(dtype=float),
    }
    return part, coords, data, pooled_sd


def frequentist_check(part, depth_adjusted, spec_id):
    dummies = pd.get_dummies(part["unit"], drop_first=True, dtype=float)
    preds = ["age_scaled", "sex_centered"]
    if depth_adjusted:
        preds.append("depth_centered")
    X = sm.add_constant(pd.concat([part[preds], dummies], axis=1), has_constant="add")
    fit = sm.OLS(part["y"], X).fit(
        cov_type="cluster", cov_kwds={"groups": part["donor_cluster"]}, use_t=True)
    ci = fit.conf_int().loc["age_scaled"]
    return {
        "spec_id": spec_id,
        "n_donors": int(part["donor_cluster"].nunique()),
        "n_tissues": int(part["group"].nunique()),
        "n_units": int(part["unit"].nunique()),
        "n_obs": int(len(part)),
        "freq_effect_per_12_25y": float(fit.params["age_scaled"]),
        "freq_se": float(fit.bse["age_scaled"]),
        "freq_ci95_low": float(ci.iloc[0]),
        "freq_ci95_high": float(ci.iloc[1]),
        "freq_pvalue": float(fit.pvalues["age_scaled"]),
    }


def extract_rows(idata, part, coords, spec_id, depth_adjusted, prior_sd, pooled_sd):
    hdi = az.hdi(idata, hdi_prob=0.94)
    post = idata.posterior
    rows = []

    def add_row(scope, label, vals, hdi_da):
        vals = np.asarray(vals, dtype=float).ravel()
        lo95, hi95 = np.quantile(vals, [0.025, 0.975])
        rows.append({
            "spec_id": spec_id,
            "model": "depth_adjusted" if depth_adjusted else "unadjusted",
            "prior_sd": float(prior_sd),
            "scope": scope,
            "unit": label,
            "n_donors": int(part["donor_cluster"].nunique()),
            "n_tissues": int(part["group"].nunique()),
            "n_units": int(len(coords["unit"])),
            "n_obs": int(len(part)),
            "posterior_median": float(np.median(vals)),
            "posterior_mean": float(vals.mean()),
            "posterior_sd": float(vals.std()),
            "ci95_low": float(lo95),
            "ci95_high": float(hi95),
            "hdi94_low": float(hdi_da.sel(hdi="lower")),
            "hdi94_high": float(hdi_da.sel(hdi="higher")),
            "prob_positive": float(np.mean(vals > 0.0)),
            "prob_negative": float(np.mean(vals < 0.0)),
            "prob_rope_0.05": float(np.mean(np.abs(vals) < 0.05)),
            "prob_rope_0.025": float(np.mean(np.abs(vals) < 0.025)),
            "prob_rope_0.10": float(np.mean(np.abs(vals) < 0.10)),
            # raw Hellinger per 12.25 years: multiply by pooled within-unit SD
            "posterior_median_raw": float(np.median(vals) * pooled_sd),
            "ci95_low_raw": float(lo95 * pooled_sd),
            "ci95_high_raw": float(hi95 * pooled_sd),
            "pooled_within_unit_sd": pooled_sd,
        })

    add_row("species", "ALL", post["global_age"].to_numpy(), hdi["global_age"])
    if "global_sex" in post:
        add_row("covariate", "global_sex", post["global_sex"].to_numpy(), hdi["global_sex"])
    if "global_depth" in post:
        add_row("covariate", "global_depth", post["global_depth"].to_numpy(), hdi["global_depth"])
    if "group_age" in post:
        for grp in coords["group"]:
            vals = post["group_age"].sel(group=grp).to_numpy().ravel()
            add_row("tissue", grp, vals, hdi["group_age"].sel(group=grp))
    for u in coords["unit"]:
        vals = post["unit_age"].sel(unit=u).to_numpy().ravel()
        add_row("cell_type_unit", u, vals, hdi["unit_age"].sel(unit=u))
    return rows


def run_one(part, coords, data, depth_adjusted, prior_sd, seed, spec_id):
    rng = np.random.default_rng(seed)
    model = build_hierarchical_student_t(
        coords, data, prior_sd=prior_sd, depth_adjusted=depth_adjusted,
        include_sex=True)
    with model:
        prior = pm.sample_prior_predictive(draws=500, random_seed=rng)
        idata = pm.sample(
            draws=2000, tune=1000, chains=4, nuts_sampler="nutpie",
            target_accept=0.95, random_seed=rng, progressbar=False)
        idata.extend(prior)
        idata.extend(pm.sample_posterior_predictive(idata, random_seed=rng, progressbar=False))
        pm.compute_log_likelihood(idata, model=model)
    diag, _ = summarize_diagnostics(
        idata, data["y"], model_id=spec_id, family="human_drift_mature",
        prior_sd=prior_sd)
    return idata, diag


def main():
    donors = pd.read_csv(DATA)
    part, coords, data, pooled_sd = prepare_human_data(donors)
    print(f"[human] rows={len(part)} donors={part['donor_cluster'].nunique()} "
          f"units={len(coords['unit'])} tissues={len(coords['group'])} "
          f"pooled_within_unit_sd={pooled_sd:.4f}", flush=True)
    specs = [
        ("mature_human-profile_hellinger-unadjusted-prior_0.50", False, 0.50, SEED),
        ("mature_human-profile_hellinger-depth_adjusted-prior_0.50", True, 0.50, SEED),
        ("mature_human-profile_hellinger-depth_adjusted-prior_0.25", True, 0.25, SEED),
        ("mature_human-profile_hellinger-depth_adjusted-prior_1.00", True, 1.00, SEED),
        # second-seed rerun of the primary depth-adjusted spec
        ("mature_human-profile_hellinger-depth_adjusted-prior_0.50-seed2", True, 0.50, SEED2),
    ]
    post_rows, diag_rows, freq_rows = [], [], []
    for spec_id, depth_adj, prior_sd, seed in specs:
        print(f"[human] Fitting {spec_id}", flush=True)
        freq_rows.append(frequentist_check(part, depth_adj, spec_id))
        idata, diag = run_one(part, coords, data, depth_adj, prior_sd, seed, spec_id)
        diag_rows.append(diag)
        rows_now = extract_rows(idata, part, coords, spec_id, depth_adj, prior_sd, pooled_sd)
        post_rows.extend(rows_now)
        sp = [r for r in rows_now if r["scope"] == "species"][0]
        print(f"[human] Done {spec_id}: median={sp['posterior_median']:+.4f} "
              f"95% CI [{sp['ci95_low']:+.4f}, {sp['ci95_high']:+.4f}] "
              f"P(>0)={sp['prob_positive']:.4f} max_rhat={diag['max_rhat']:.4f} "
              f"div={diag['divergences']} min_bulk_ESS={diag.get('min_ess_bulk', float('nan')):.0f}",
              flush=True)
    pd.DataFrame(post_rows).to_csv(OUT / "followup-mature-bayesian-human-drift.csv", index=False)
    pd.DataFrame(diag_rows).to_csv(OUT / "followup-mature-bayesian-human-drift-diagnostics.csv", index=False)
    pd.DataFrame(freq_rows).to_csv(OUT / "followup-mature-bayesian-human-drift-frequentist.csv", index=False)
    print("[human] Wrote followup-mature-bayesian-human-drift*.csv", flush=True)


if __name__ == "__main__":
    main()
