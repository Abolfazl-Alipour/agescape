#!/usr/bin/env python3
"""Shared ArviZ convergence, LOO, prior-predictive, and PPC diagnostics."""

from __future__ import annotations

import warnings
import arviz as az
import numpy as np
import pymc as pm


def build_hierarchical_student_t(
    coords: dict,
    data: dict,
    prior_sd: float,
    depth_adjusted: bool = False,
    fixed_group_effects: bool = False,
    include_sex: bool = True,
    include_donor: bool = True,
) -> pm.Model:
    """Hierarchical Student-t age-slope model shared across vector and full-paper Bayes."""
    with pm.Model(coords=coords) as model:
        age = pm.Data("age", data["age"], dims="obs")
        unit_index = pm.Data("unit_index", data["unit_index"], dims="obs")
        global_intercept = pm.Normal("global_intercept", 0.0, 0.5)
        global_age = pm.Normal("global_age", 0.0, float(prior_sd))
        tau_intercept = pm.Gamma("tau_intercept", 2.0, 2.0)
        tau_unit_age = pm.Gamma("tau_unit_age", 2.0, 2.0)
        z_intercept = pm.Normal("z_intercept", 0.0, 1.0, dims="unit")
        z_unit_age = pm.Normal("z_unit_age", 0.0, 1.0, dims="unit")
        unit_intercept = pm.Deterministic(
            "unit_intercept", global_intercept + tau_intercept * z_intercept, dims="unit"
        )
        if "group" in coords and len(coords["group"]) > 1 and "unit_group_index" in data:
            unit_group_index = pm.Data("unit_group_index", data["unit_group_index"], dims="unit")
            if fixed_group_effects:
                group_offset = pm.Normal("group_offset", 0.0, 0.25, dims="group")
                group_age = pm.Deterministic("group_age", global_age + group_offset, dims="group")
            else:
                tau_group_age = pm.Gamma("tau_group_age", 2.0, 2.0)
                z_group_age = pm.Normal("z_group_age", 0.0, 1.0, dims="group")
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

        mu = unit_intercept[unit_index] + unit_age[unit_index] * age
        if include_donor and "donor" in coords and "donor_index" in data:
            donor_index = pm.Data("donor_index", data["donor_index"], dims="obs")
            tau_donor = pm.Gamma("tau_donor", 2.0, 2.0)
            z_donor = pm.Normal("z_donor", 0.0, 1.0, dims="donor")
            donor_intercept = pm.Deterministic("donor_intercept", tau_donor * z_donor, dims="donor")
            mu = mu + donor_intercept[donor_index]
        if include_sex and "sex" in data:
            sex = pm.Data("sex", data["sex"], dims="obs")
            global_sex = pm.Normal("global_sex", 0.0, 0.5)
            mu = mu + global_sex * sex
        if depth_adjusted and "depth" in data:
            depth = pm.Data("depth", data["depth"], dims="obs")
            global_depth = pm.Normal("global_depth", 0.0, 0.5)
            mu = mu + global_depth * depth

        sigma = pm.Gamma("sigma", 2.0, 2.0, dims="unit")
        nu = pm.Gamma("nu", 2.0, 0.1)
        pm.StudentT(
            "y_obs", nu=nu, mu=mu, sigma=sigma[unit_index], observed=data["y"], dims="obs"
        )
    return model


def _var_stat(stat_ds) -> tuple[float, float]:
    vals = []
    for var in stat_ds.data_vars:
        arr = np.asarray(stat_ds[var].to_numpy(), dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size:
            vals.append((float(arr.min()), float(arr.max())))
    if not vals:
        return np.nan, np.nan
    return min(v[0] for v in vals), max(v[1] for v in vals)


def summarize_diagnostics(
    idata,
    observed_y: np.ndarray,
    model_id: str,
    family: str,
    prior_sd: float,
    obs_var: str = "y_obs",
) -> tuple[dict, dict]:
    """Compute MCMC convergence, Pareto-k LOO, PPC 90% coverage, and prior coverage."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rhat_ds = az.rhat(idata)
        ess_bulk_ds = az.ess(idata, method="bulk")
        ess_tail_ds = az.ess(idata, method="tail")
        _, max_rhat = _var_stat(rhat_ds)
        min_ess_bulk, _ = _var_stat(ess_bulk_ds)
        min_ess_tail, _ = _var_stat(ess_tail_ds)
        divergences = int(np.asarray(idata.sample_stats["diverging"]).sum())
        bfmi_vals = np.asarray(az.bfmi(idata), dtype=float)
        min_bfmi = float(np.nanmin(bfmi_vals)) if bfmi_vals.size else np.nan

        loo_res = az.loo(idata, pointwise=True)
        pareto_k = np.asarray(loo_res.pareto_k, dtype=float)
        pareto_k_max = float(np.nanmax(pareto_k)) if pareto_k.size else np.nan
        pareto_k_n_bad = int(np.sum(pareto_k > 0.7))

        ppc_vals = np.asarray(idata.posterior_predictive[obs_var].to_numpy(), dtype=float)
        ppc_flat = ppc_vals.reshape(-1, len(observed_y))
        q05 = np.quantile(ppc_flat, 0.05, axis=0)
        q95 = np.quantile(ppc_flat, 0.95, axis=0)
        ppc_90_cov = float(np.mean((observed_y >= q05) & (observed_y <= q95)))

        prior_vals = np.asarray(idata.prior_predictive[obs_var].to_numpy(), dtype=float).ravel()
        prior_q001 = float(np.quantile(prior_vals, 0.001))
        prior_q999 = float(np.quantile(prior_vals, 0.999))

    convergence_ok = bool(
        np.isfinite(max_rhat)
        and max_rhat <= 1.01
        and np.isfinite(min_ess_bulk)
        and min_ess_bulk >= 400
        and np.isfinite(min_ess_tail)
        and min_ess_tail >= 400
        and divergences == 0
    )

    diag = {
        "family": family,
        "model_id": model_id,
        "prior_sd": float(prior_sd),
        "n_obs": int(len(observed_y)),
        "max_rhat": float(max_rhat),
        "min_ess_bulk": float(min_ess_bulk),
        "min_ess_tail": float(min_ess_tail),
        "divergences": divergences,
        "min_bfmi": min_bfmi,
        "elpd_loo": float(loo_res.elpd_loo),
        "p_loo": float(loo_res.p_loo),
        "pareto_k_max": pareto_k_max,
        "pareto_k_n_bad": pareto_k_n_bad,
        "ppc_90_coverage": ppc_90_cov,
        "convergence_ok": convergence_ok,
    }
    prior = {
        "family": family,
        "model_id": model_id,
        "prior_sd": float(prior_sd),
        "n_obs": int(len(observed_y)),
        "observed_min": float(np.min(observed_y)),
        "observed_max": float(np.max(observed_y)),
        "prior_q001": prior_q001,
        "prior_q999": prior_q999,
        "observed_within_prior_998": bool(
            np.min(observed_y) >= prior_q001 and np.max(observed_y) <= prior_q999
        ),
    }
    return diag, prior
