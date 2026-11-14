#!/usr/bin/env python3
"""Item 7: Bayesian hierarchical synthesis across all five Agescape analytical families."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import warnings

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.bayesian_diagnostics_helper import (
    build_hierarchical_student_t,
    summarize_diagnostics,
)

RESULTS = REPO / "results"
SEED = sum(map(ord, "agescape-rev3-full-paper-bayes-v1"))


def _sample_and_diagnose(
    model: pm.Model,
    observed_y: np.ndarray,
    model_id: str,
    family: str,
    prior_sd: float,
    seed_offset: int,
    draws: int = 2000,
    tune: int = 2000,
    target_accept: float = 0.96,
):
    rng = np.random.default_rng(SEED + seed_offset)
    with model:
        prior = pm.sample_prior_predictive(draws=500, random_seed=rng)
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=4,
            nuts_sampler="nutpie",
            target_accept=target_accept,
            random_seed=rng,
            progressbar=False,
        )
        idata.extend(prior)
        idata.extend(pm.sample_posterior_predictive(idata, random_seed=rng, progressbar=False))
        pm.compute_log_likelihood(idata, model=model)
    diag, prior_sum = summarize_diagnostics(
        idata, observed_y, model_id=model_id, family=family, prior_sd=prior_sd
    )
    print(
        f"[item7] Done {model_id}: max_rhat={diag['max_rhat']:.4f} "
        f"ess={diag['min_ess_bulk']:.0f} div={diag['divergences']} ok={diag['convergence_ok']}",
        flush=True,
    )
    return idata, diag, prior_sum


def _summarize_array(
    family: str,
    model_id: str,
    prior_sd: float,
    scope: str,
    subgroup: str,
    parameter: str,
    vals: np.ndarray,
    n_units: int,
    n_obs: int,
    extra: dict | None = None,
) -> dict:
    vals = np.asarray(vals, dtype=float).ravel()
    low, high = az.hdi(vals, hdi_prob=0.94)
    p_pos = float(np.mean(vals > 0.0))
    lfsr = float(min(p_pos, 1.0 - p_pos))
    row = {
        "family": family,
        "model_id": model_id,
        "prior_sd": float(prior_sd),
        "scope": scope,
        "subgroup": subgroup,
        "parameter": parameter,
        "n_units": int(n_units),
        "n_obs": int(n_obs),
        "posterior_mean": float(np.mean(vals)),
        "posterior_sd": float(np.std(vals)),
        "hdi94_low": float(low),
        "hdi94_high": float(high),
        "prob_positive": p_pos,
        "lfsr": lfsr,
        "prob_abs_gt_0_10": float(np.mean(np.abs(vals) > 0.10)),
        "prob_abs_gt_0_20": float(np.mean(np.abs(vals) > 0.20)),
    }
    if extra:
        row.update(extra)
    return row


def run_family_composition(summary_rows, diag_rows, prior_rows):
    """Family (a): 620 CLR tests across 5 species (+ NMR note) and mature-only mouse/human."""
    comp = pd.read_csv(RESULTS / "rev3-composition-global-bh-624.csv")
    clr = comp[comp["se"].notna() & (comp["se"] > 0)].copy()
    clr["d_obs"] = (clr["effect"] / clr["se"]) / np.sqrt(clr["n_donors"])
    clr["d_se"] = 1.0 / np.sqrt(clr["n_donors"])
    species_list = sorted(clr["species"].unique())
    sp_map = {s: i for i, s in enumerate(species_list)}
    clr["sp_idx"] = clr["species"].map(sp_map)

    for prior_sd in [0.25, 0.50]:
        model_id = f"composition_global_620_std-prior_{prior_sd:.2f}"
        print(f"[item7] Fitting {model_id}", flush=True)
        coords = {"species": species_list, "unit": clr.index.to_numpy()}
        with pm.Model(coords=coords) as model:
            sp_idx = pm.Data("sp_idx", clr["sp_idx"].to_numpy(dtype=int), dims="unit")
            d_se = pm.Data("d_se", clr["d_se"].to_numpy(dtype=float), dims="unit")
            mu_sp = pm.Normal("mu_species", 0.0, float(prior_sd), dims="species")
            tau_sp = pm.HalfNormal("tau_species", float(prior_sd), dims="species")
            z_raw = pm.Normal("z_raw", 0.0, 1.0, dims="unit")
            theta_unit = pm.Deterministic("theta_unit", mu_sp[sp_idx] + tau_sp[sp_idx] * z_raw, dims="unit")
            pm.Normal("y_obs", mu=theta_unit, sigma=d_se, observed=clr["d_obs"].to_numpy(dtype=float), dims="unit")
        idata, diag, psum = _sample_and_diagnose(
            model, clr["d_obs"].to_numpy(dtype=float), model_id, "composition", prior_sd, seed_offset=10 + int(prior_sd * 100), draws=2500, tune=2500, target_accept=0.99
        )
        diag_rows.append(diag)
        prior_rows.append(psum)
        post = idata.posterior
        for sp in species_list:
            sub = clr[clr["species"] == sp]
            idx_mask = (clr["species"] == sp).to_numpy()
            mu_vals = post["mu_species"].sel(species=sp).to_numpy().ravel()
            tau_vals = post["tau_species"].sel(species=sp).to_numpy().ravel()
            theta_sp = post["theta_unit"].to_numpy()[:, :, idx_mask]
            # Per-unit local false sign rate across units in this species
            p_pos_unit = (theta_sp > 0).mean(axis=(0, 1))
            lfsr_unit = np.minimum(p_pos_unit, 1.0 - p_pos_unit)
            abs_mean_unit = np.abs(theta_sp).mean(axis=(0, 1))
            extra = {
                "mean_unit_abs_theta": float(np.mean(abs_mean_unit)),
                "frac_units_lfsr_lt_05": float(np.mean(lfsr_unit < 0.05)),
                "n_units_lfsr_lt_05": int(np.sum(lfsr_unit < 0.05)),
            }
            summary_rows.append(_summarize_array("composition", model_id, prior_sd, "species", sp, "mu_species_std", mu_vals, len(sub), int(sub["n_donors"].max()), extra))
            summary_rows.append(_summarize_array("composition", model_id, prior_sd, "species", sp, "tau_species_std", tau_vals, len(sub), int(sub["n_donors"].max()), extra))

    # Full vs mature-only CLR slope per 10% lifespan in Mouse (4y) and Human (122.5y)
    arms = [
        ("mouse_full_136", clr[clr["species"] == "Mouse"].copy(), 0.4),
        ("mouse_mature_18m_114", pd.read_csv(RESULTS / "mature-only-composition-mouse.csv"), 0.4),
        ("human_full_291", clr[clr["species"] == "Human"].copy(), 12.25),
        ("human_mature_30y_200", pd.read_csv(RESULTS / "mature-only-composition-human.csv"), 12.25),
    ]
    for idx_arm, (arm_name, df, span) in enumerate(arms):
        df = df[df["se"].notna() & (df["se"] > 0)].copy()
        y_obs = (df["effect"] * span).to_numpy(dtype=float)
        y_se = (df["se"] * span).to_numpy(dtype=float)
        for prior_sd in [0.25, 0.50]:
            model_id = f"composition_{arm_name}-prior_{prior_sd:.2f}"
            print(f"[item7] Fitting {model_id}", flush=True)
            coords = {"unit": np.arange(len(df))}
            with pm.Model(coords=coords) as model:
                mu = pm.Normal("mu", 0.0, float(prior_sd))
                tau = pm.HalfNormal("tau", float(prior_sd))
                nu = pm.Gamma("nu", 3.0, 0.2)
                t_raw = pm.StudentT("t_raw", nu=nu, mu=0.0, sigma=1.0, dims="unit")
                theta = pm.Deterministic("theta", mu + tau * t_raw, dims="unit")
                pm.Normal("y_obs", mu=theta, sigma=y_se, observed=y_obs, dims="unit")
            idata, diag, psum = _sample_and_diagnose(
                model, y_obs, model_id, "composition", prior_sd, seed_offset=100 + idx_arm * 10 + int(prior_sd * 100)
            )
            diag_rows.append(diag)
            prior_rows.append(psum)
            post = idata.posterior
            theta_draws = post["theta"].to_numpy()
            lfsr_u = np.minimum((theta_draws > 0).mean(axis=(0, 1)), (theta_draws < 0).mean(axis=(0, 1)))
            extra = {
                "mean_unit_abs_theta": float(np.abs(theta_draws).mean()),
                "frac_units_lfsr_lt_05": float(np.mean(lfsr_u < 0.05)),
                "n_units_lfsr_lt_05": int(np.sum(lfsr_u < 0.05)),
            }
            summary_rows.append(_summarize_array("composition", model_id, prior_sd, "arm", arm_name, "mu_clr_10pct_lifespan", post["mu"].to_numpy(), len(df), int(df["n_donors"].max()), extra))
            summary_rows.append(_summarize_array("composition", model_id, prior_sd, "arm", arm_name, "tau_clr_10pct_lifespan", post["tau"].to_numpy(), len(df), int(df["n_donors"].max()), extra))


def run_family_transport(summary_rows, diag_rows, prior_rows):
    """Family (b): 252 human P1_Clean transport slopes + random-gene-set excess in PsychAD & MSSM."""
    m9 = pd.read_csv(RESULTS / "m9-regressions.csv")[["cohort", "cell_type", "n_donors", "slope_per_year", "se"]].rename(columns={"slope_per_year": "slope"})
    m10 = pd.read_csv(RESULTS / "m10-celltype-regressions.csv")[["cell_type", "n_donors", "slope", "se"]].copy()
    m10.insert(0, "cohort", "PsychAD")
    all252 = pd.concat([m9, m10], ignore_index=True)
    all252["d_obs"] = (all252["slope"] / all252["se"]) / np.sqrt(all252["n_donors"])
    all252["d_se"] = 1.0 / np.sqrt(all252["n_donors"])
    cohorts = sorted(all252["cohort"].unique())
    cmap = {c: i for i, c in enumerate(cohorts)}
    c_idx = all252["cohort"].map(cmap).to_numpy(dtype=int)

    for prior_sd in [0.25, 0.50]:
        model_id = f"transport_252_std-prior_{prior_sd:.2f}"
        print(f"[item7] Fitting {model_id}", flush=True)
        coords = {"cohort": cohorts, "unit": np.arange(len(all252))}
        with pm.Model(coords=coords) as model:
            c_idx_d = pm.Data("c_idx", c_idx, dims="unit")
            d_se_d = pm.Data("d_se", all252["d_se"].to_numpy(dtype=float), dims="unit")
            mu_global = pm.Normal("mu_global", 0.0, float(prior_sd))
            tau_cohort = pm.HalfNormal("tau_cohort", float(prior_sd))
            z_cohort = pm.Normal("z_cohort", 0.0, 1.0, dims="cohort")
            mu_cohort = pm.Deterministic("mu_cohort", mu_global + tau_cohort * z_cohort, dims="cohort")
            tau_unit = pm.HalfNormal("tau_unit", float(prior_sd))
            z_raw = pm.Normal("z_raw", 0.0, 1.0, dims="unit")
            theta = pm.Deterministic("theta", mu_cohort[c_idx_d] + tau_unit * z_raw, dims="unit")
            pm.Normal("y_obs", mu=theta, sigma=d_se_d, observed=all252["d_obs"].to_numpy(dtype=float), dims="unit")
        idata, diag, psum = _sample_and_diagnose(
            model, all252["d_obs"].to_numpy(dtype=float), model_id, "transport", prior_sd, seed_offset=200 + int(prior_sd * 100), draws=2500, tune=2500, target_accept=0.99
        )
        diag_rows.append(diag)
        prior_rows.append(psum)
        summary_rows.append(_summarize_array("transport", model_id, prior_sd, "global", "ALL_252", "mu_global_std", idata.posterior["mu_global"].to_numpy(), len(all252), len(all252)))
        summary_rows.append(_summarize_array("transport", model_id, prior_sd, "global", "ALL_252", "tau_cohort_std", idata.posterior["tau_cohort"].to_numpy(), len(cohorts), len(all252)))
        for ch in cohorts:
            sub = all252[all252["cohort"] == ch]
            summary_rows.append(_summarize_array("transport", model_id, prior_sd, "cohort", ch, "mu_cohort_std", idata.posterior["mu_cohort"].sel(cohort=ch).to_numpy(), len(sub), int(sub["n_donors"].max())))

    # Random-gene-set calibration excess (actual_slope - random_mean) / random_sd scaled by 1/sqrt(n_donors) vs raw z
    rn_p = pd.read_csv(RESULTS / "reviewer-phase1-signature-random-null.csv")
    rn_p["cohort"] = "PsychAD"
    rn_p["z_excess"] = (rn_p["actual_slope"] - rn_p["random_mean_slope"]) / rn_p["random_sd"]
    rn_m = pd.read_csv(RESULTS / "reviewer-rev2-transport-null-mssm.csv")
    rn_m["z_excess"] = (rn_m["observed_slope"] - rn_m["random_mean_slope"]) / rn_m["random_sd"]
    rn = pd.concat([rn_p[["cohort", "cell_type", "n_donors", "z_excess"]], rn_m[["cohort", "cell_type", "n_donors", "z_excess"]]], ignore_index=True)
    # Fisher-like normalized excess per random SD unit (null mean=0, null SD=1 across random gene sets)
    y_z = rn["z_excess"].to_numpy(dtype=float)
    r_cohorts = sorted(rn["cohort"].unique())
    rmap = {c: i for i, c in enumerate(r_cohorts)}
    r_idx = rn["cohort"].map(rmap).to_numpy(dtype=int)

    for prior_sd in [0.50, 1.00]:
        model_id = f"transport_random_null_z_excess-prior_{prior_sd:.2f}"
        print(f"[item7] Fitting {model_id}", flush=True)
        coords = {"cohort": r_cohorts, "unit": np.arange(len(rn))}
        with pm.Model(coords=coords) as model:
            c_d = pm.Data("c_d", r_idx, dims="unit")
            mu_global = pm.Normal("mu_global", 0.0, float(prior_sd))
            tau_c = pm.HalfNormal("tau_cohort", 0.5)
            z_c = pm.Normal("z_cohort", 0.0, 1.0, dims="cohort")
            mu_c = pm.Deterministic("mu_cohort", mu_global + tau_c * z_c, dims="cohort")
            sigma = pm.HalfNormal("sigma", 1.5, dims="cohort")
            nu = pm.Gamma("nu", 3.0, 0.2)
            pm.StudentT("y_obs", nu=nu, mu=mu_c[c_d], sigma=sigma[c_d], observed=y_z, dims="unit")
        idata, diag, psum = _sample_and_diagnose(
            model, y_z, model_id, "transport", prior_sd, seed_offset=250 + int(prior_sd * 100), draws=2500, tune=2500, target_accept=0.99
        )
        diag_rows.append(diag)
        prior_rows.append(psum)
        summary_rows.append(_summarize_array("transport", model_id, prior_sd, "random_null", "ALL_48", "mu_global_z_excess", idata.posterior["mu_global"].to_numpy(), len(rn), len(rn)))
        for ch in r_cohorts:
            sub = rn[rn["cohort"] == ch]
            summary_rows.append(_summarize_array("transport", model_id, prior_sd, "random_null", ch, "mu_cohort_z_excess", idata.posterior["mu_cohort"].sel(cohort=ch).to_numpy(), len(sub), len(sub)))


def run_family_concordance(summary_rows, diag_rows, prior_rows):
    """Family (c): Cross-species within-tissue concordance on Fisher-z permutation excess."""
    rows = []
    # Mouse 11 tissues
    ms = pd.read_csv(RESULTS / "reviewer-rev2-concordance-strata.csv")
    for _, r in ms.iterrows():
        rows.append({"species": "Mouse", "stratum": r["tissue"], "obs": r["observed_mean_pearson_primary"], "null_mean": r["joint_null_mean_primary"], "null_sd": r["joint_null_sd_primary"], "n_units": r["n_units"], "n_pairs": r["n_pairs_within"]})
    # Fly, Worm, Human, Zebrafish
    for fn in [
        "reviewer-rev2-joint-concordance-fly-strata.csv",
        "reviewer-rev2-joint-concordance-worm-strata.csv",
        "reviewer-rev2-joint-concordance-human-strata.csv",
        "reviewer-rev2-joint-concordance-zebrafish-strata.csv",
    ]:
        df = pd.read_csv(RESULTS / fn)
        for _, r in df.iterrows():
            rows.append({"species": r["species"], "stratum": f"{r['tissue']} ({r['cohort']})", "obs": r["observed_mean_pearson"], "null_mean": r["joint_null_mean"], "null_sd": r["joint_null_sd"], "n_units": r["n_units"], "n_pairs": r["n_pairs"]})
    # NMR from main table
    cmain = pd.read_csv(RESULTS / "reviewer-rev2-concordance.csv")
    nmr = cmain[cmain["species"] == "Naked mole-rat"].iloc[0]
    rows.append({"species": "Naked mole-rat", "stratum": "Skin (GSE152886)", "obs": nmr["observed_mean_pearson"], "null_mean": nmr["joint_null_mean"], "null_sd": nmr["joint_null_sd"], "n_units": nmr["n_units"], "n_pairs": nmr["n_pairs"]})

    cdf = pd.DataFrame(rows)
    cdf["y_fisher_excess"] = np.arctanh(np.clip(cdf["obs"], -0.99, 0.99)) - np.arctanh(np.clip(cdf["null_mean"], -0.99, 0.99))
    lo_z = np.arctanh(np.clip(cdf["null_mean"] - cdf["null_sd"], -0.99, 0.99))
    hi_z = np.arctanh(np.clip(cdf["null_mean"] + cdf["null_sd"], -0.99, 0.99))
    cdf["s_fisher"] = 0.5 * (hi_z - lo_z)
    cdf["z_perm"] = (cdf["obs"] - cdf["null_mean"]) / cdf["null_sd"]

    species_list = sorted(cdf["species"].unique())
    sp_map = {s: i for i, s in enumerate(species_list)}
    sp_idx = cdf["species"].map(sp_map).to_numpy(dtype=int)

    for prior_sd in [0.25, 0.50]:
        model_id = f"concordance_fisher_excess_6species-prior_{prior_sd:.2f}"
        print(f"[item7] Fitting {model_id}", flush=True)
        coords = {"species": species_list, "stratum": np.arange(len(cdf))}
        with pm.Model(coords=coords) as model:
            s_idx = pm.Data("s_idx", sp_idx, dims="stratum")
            s_se = pm.Data("s_se", cdf["s_fisher"].to_numpy(dtype=float), dims="stratum")
            mu_all = pm.Normal("mu_global", 0.0, float(prior_sd))
            tau_between_sp = pm.HalfNormal("tau_between_species", float(prior_sd))
            z_sp = pm.Normal("z_sp", 0.0, 1.0, dims="species")
            mu_sp = pm.Deterministic("mu_species", mu_all + tau_between_sp * z_sp, dims="species")
            tau_stratum = pm.HalfNormal("tau_stratum", 0.25)
            z_st = pm.Normal("z_st", 0.0, 1.0, dims="stratum")
            theta_st = pm.Deterministic("theta_stratum", mu_sp[s_idx] + tau_stratum * z_st, dims="stratum")
            pm.Normal("y_obs", mu=theta_st, sigma=s_se, observed=cdf["y_fisher_excess"].to_numpy(dtype=float), dims="stratum")
        idata, diag, psum = _sample_and_diagnose(
            model, cdf["y_fisher_excess"].to_numpy(dtype=float), model_id, "concordance", prior_sd, seed_offset=300 + int(prior_sd * 100), draws=3000, tune=3000, target_accept=0.99
        )
        diag_rows.append(diag)
        prior_rows.append(psum)
        post = idata.posterior
        summary_rows.append(_summarize_array("concordance", model_id, prior_sd, "global", "ALL_6_SPECIES", "mu_global_fisher_excess", post["mu_global"].to_numpy(), len(cdf), len(cdf)))
        summary_rows.append(_summarize_array("concordance", model_id, prior_sd, "global", "ALL_6_SPECIES", "tau_between_species", post["tau_between_species"].to_numpy(), len(species_list), len(cdf)))
        for sp in species_list:
            sub = cdf[cdf["species"] == sp]
            vals = post["mu_species"].sel(species=sp).to_numpy()
            summary_rows.append(_summarize_array("concordance", model_id, prior_sd, "species", sp, "mu_species_fisher_excess", vals, int(sub["n_units"].sum()), len(sub), {"mean_raw_z_perm": float(sub["z_perm"].mean()), "n_strata": len(sub)}))
        # Species contrasts vs Human
        h_vals = post["mu_species"].sel(species="Human").to_numpy()
        for sp in ["Fruit fly", "C. elegans", "Mouse"]:
            diff = post["mu_species"].sel(species=sp).to_numpy() - h_vals
            summary_rows.append(_summarize_array("concordance", model_id, prior_sd, "contrast", f"{sp}_minus_Human", "fisher_excess_diff", diff, len(cdf), len(cdf)))


def _fit_donor_endpoint(df: pd.DataFrame, endpoint_col: str, age_col: str, family_label: str, spec_prefix: str, summary_rows, diag_rows, prior_rows, seed_base: int):
    part = df.copy()
    part["unit"] = part["cell_type"].astype(str)
    part["group"] = "Lung"
    part["donor_cluster"] = part["donor"].astype(str)
    part = part.groupby("unit").filter(lambda f: len(f) >= 6).copy()
    means = part.groupby("unit")[endpoint_col].transform("mean")
    sds = part.groupby("unit")[endpoint_col].transform("std")
    part["y"] = (part[endpoint_col] - means) / sds
    part = part[np.isfinite(part["y"])].copy()
    part = part.groupby("unit").filter(lambda f: len(f) >= 6).copy()
    units = sorted(part["unit"].unique())
    donors = sorted(part["donor_cluster"].unique())
    umap = {v: i for i, v in enumerate(units)}
    dmap = {v: i for i, v in enumerate(donors)}
    age_years = part[age_col].to_numpy(dtype=float) / 12.0
    part["age_scaled"] = (age_years - float(np.mean(age_years))) / 0.4
    part["sex_centered"] = part["sex_male"] - part.groupby("unit")["sex_male"].transform("mean")
    part["depth_centered"] = part["mean_log_umi_actual"] - part.groupby("unit")["mean_log_umi_actual"].transform("mean")
    part = part.sort_values(["unit", "donor_cluster"]).reset_index(drop=True)

    coords = {"unit": units, "donor": donors, "obs": part.index.to_numpy()}
    data = {
        "age": part["age_scaled"].to_numpy(dtype=float),
        "sex": part["sex_centered"].to_numpy(dtype=float),
        "depth": part["depth_centered"].to_numpy(dtype=float),
        "unit_index": part["unit"].map(umap).to_numpy(dtype=int),
        "donor_index": part["donor_cluster"].map(dmap).to_numpy(dtype=int),
        "y": part["y"].to_numpy(dtype=float),
    }
    for depth_adj in [False, True]:
        for prior_sd in [0.25, 0.50]:
            mod_name = "depth_adjusted" if depth_adj else "unadjusted"
            model_id = f"{spec_prefix}-{mod_name}-prior_{prior_sd:.2f}"
            print(f"[item7] Fitting {model_id}", flush=True)
            model = build_hierarchical_student_t(coords, data, prior_sd=prior_sd, depth_adjusted=depth_adj)
            idata, diag, psum = _sample_and_diagnose(
                model, data["y"], model_id, family_label, prior_sd, seed_offset=seed_base + (50 if depth_adj else 0) + int(prior_sd * 100), draws=2500, tune=2500, target_accept=0.98
            )
            diag_rows.append(diag)
            prior_rows.append(psum)
            summary_rows.append(
                _summarize_array(
                    family_label, model_id, prior_sd, "global", spec_prefix, f"global_age_10pct_lifespan_{mod_name}", idata.posterior["global_age"].to_numpy(), len(units), len(part)
                )
            )


def run_family_disorder_entropy(summary_rows, diag_rows, prior_rows):
    """Family (d): Mouse lung full-range (1-30m) off-program & entropy + mature-only (>=18m) entropy & D_KL."""
    p2 = pd.read_csv(RESULTS / "reviewer-phase2-depth-donor-results.csv")
    mat = pd.read_csv(RESULTS / "rev3-entropy-mature-only-mouse-lung-donors.csv")
    mat_t = mat[mat["is_manuscript_target"] == True].copy()
    jobs = [
        ("lung_full_1500_offprogram", p2[p2["depth"] == "1500"], "offprogram_fraction", "age_months", 400),
        ("lung_full_full_offprogram", p2[p2["depth"] == "full"], "offprogram_fraction", "age_months", 500),
        ("lung_full_1500_entropy", p2[p2["depth"] == "1500"], "entropy", "age_months", 600),
        ("lung_mature18m_1500_entropy", mat_t[mat_t["depth"] == "1500"], "shannon_entropy", "age_months", 700),
        ("lung_mature18m_1500_deman_kl", mat_t[mat_t["depth"] == "1500"], "deman_kl_methods", "age_months", 800),
    ]
    for name, df, col, age_col, sb in jobs:
        _fit_donor_endpoint(df, col, age_col, "disorder_entropy", name, summary_rows, diag_rows, prior_rows, sb)


def run_family_interneurons(summary_rows, diag_rows, prior_rows, only_eiv: bool = False):
    """Family (e): Human cortical interneurons EIV HBCC-vs-MSSM replication & 6 shared genes meta-model."""
    inm = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes-moderated.csv")
    # 1. Bayesian EIV replication slope regression per subtype across all MSSM-FDR genes (controls MSSM winner's curse)
    for idx_st, st in enumerate(["VIP", "pvalb", "sst"]):
        sub = inm[inm["short_type"] == st].copy()
        # Scale slopes by 100 (per century / 100y -> SD ~0.1-0.4) to match prior_sd [0.25, 0.50] cleanly
        scale = 100.0
        x_obs = (sub["slope_mssm"] * scale).to_numpy(dtype=float)
        x_se = (sub["se_mssm_mod"] * scale).to_numpy(dtype=float)
        y_obs = (sub["slope_hbcc"] * scale).to_numpy(dtype=float)
        y_se = (sub["se_hbcc_mod"] * scale).to_numpy(dtype=float)
        for prior_sd in [1.0, 2.0]:
            model_id = f"interneuron_eiv_{st}-prior_{prior_sd:.2f}"
            print(f"[item7] Fitting {model_id}", flush=True)
            coords = {"gene": np.arange(len(sub))}
            with pm.Model(coords=coords) as model:
                alpha = pm.Normal("alpha_rep", 0.0, 0.1)
                beta_rep = pm.Normal("beta_rep", 0.0, float(prior_sd))
                sigma_true = pm.HalfNormal("sigma_true_mssm", 0.5)
                mu_true = pm.Normal("mu_true_mssm", 0.0, 0.5)
                # Exact Normal-Normal conjugate factorization of p(x_true, x_obs | mu_true, sigma_true)
                pm.Normal("x_obs_lik", mu=mu_true, sigma=pm.math.sqrt(sigma_true ** 2 + x_se ** 2), observed=x_obs, dims="gene")
                var_cond = (sigma_true ** 2 * x_se ** 2) / (sigma_true ** 2 + x_se ** 2)
                mu_cond = var_cond * (mu_true / (sigma_true ** 2) + x_obs / (x_se ** 2))
                z_x = pm.Normal("z_x", 0.0, 1.0, dims="gene")
                x_true = pm.Deterministic("x_true", mu_cond + pm.math.sqrt(var_cond) * z_x, dims="gene")
                sigma_rep = pm.HalfNormal("sigma_rep_resid", 0.5)
                nu_rep = pm.Gamma("nu_rep", 3.0, 0.2)
                sigma_tot = pm.math.sqrt(sigma_rep ** 2 + y_se ** 2)
                pm.StudentT("y_obs", nu=nu_rep, mu=alpha + beta_rep * x_true, sigma=sigma_tot, observed=y_obs, dims="gene")
            idata, diag, psum = _sample_and_diagnose(
                model, y_obs, model_id, "interneurons", prior_sd, seed_offset=900 + idx_st * 20 + int(prior_sd * 10), draws=3000, tune=3000, target_accept=0.99
            )
            diag_rows.append(diag)
            prior_rows.append(psum)
            summary_rows.append(_summarize_array("interneurons", model_id, prior_sd, "subtype_eiv", st, "beta_hbcc_on_mssm", idata.posterior["beta_rep"].to_numpy(), len(sub), len(sub)))

        if only_eiv:
            continue

        # Independent HBCC-only replication magnitude on the 'both_FDR' genes (units: per 50y -> ~0.15)
        both = sub[sub["tier"] == "both_FDR"].copy()
        y50 = (both["slope_hbcc"] * 50.0).to_numpy(dtype=float)
        se50 = (both["se_hbcc_mod"] * 50.0).to_numpy(dtype=float)
        for prior_sd in [0.25, 0.50]:
            model_id = f"interneuron_hbcc_abs_{st}-prior_{prior_sd:.2f}"
            print(f"[item7] Fitting {model_id}", flush=True)
            coords = {"gene": np.arange(len(both))}
            with pm.Model(coords=coords) as model:
                mu = pm.Normal("mu_abs_hbcc", 0.0, float(prior_sd))
                tau = pm.HalfNormal("tau_abs_hbcc", float(prior_sd))
                nu = pm.Gamma("nu", 3.0, 0.2)
                t_raw = pm.StudentT("t_raw", nu=nu, mu=0.0, sigma=1.0, dims="gene")
                theta = pm.Deterministic("theta", mu + tau * t_raw, dims="gene")
                pm.Normal("y_obs", mu=theta, sigma=se50, observed=np.abs(y50), dims="gene")
            idata, diag, psum = _sample_and_diagnose(
                model, np.abs(y50), model_id, "interneurons", prior_sd, seed_offset=970 + idx_st * 20 + int(prior_sd * 100), target_accept=0.98
            )
            diag_rows.append(diag)
            prior_rows.append(psum)
            summary_rows.append(_summarize_array("interneurons", model_id, prior_sd, "both_fdr_hbcc", st, "mu_abs_50y_hbcc", idata.posterior["mu_abs_hbcc"].to_numpy(), len(both), 80))

    if only_eiv:
        return

    # 2. Hierarchical Bayesian pooling of the 6 moderated shared genes across VIP, pvalb, sst and MSSM+HBCC
    shared_genes = ["FBN1", "GABRB1", "IGF1R", "PHF3", "TRMT11", "ZBTB16"]
    sh = inm[(inm["short_type"].isin(["VIP", "pvalb", "sst"])) & (inm["gene"].isin(shared_genes))].copy()
    # Stack MSSM (unbiased check alongside HBCC) and HBCC observations per (gene, subtype) per 50 years
    obs_rows = []
    for _, r in sh.iterrows():
        obs_rows.append({"gene": r["gene"], "subtype": r["short_type"], "cohort": "MSSM", "y50": r["slope_mssm"] * 50.0, "se50": r["se_mssm_mod"] * 50.0})
        obs_rows.append({"gene": r["gene"], "subtype": r["short_type"], "cohort": "HBCC", "y50": r["slope_hbcc"] * 50.0, "se50": r["se_hbcc_mod"] * 50.0})
    sh_long = pd.DataFrame(obs_rows)
    g_map = {g: i for i, g in enumerate(shared_genes)}
    g_idx = sh_long["gene"].map(g_map).to_numpy(dtype=int)

    for prior_sd in [0.25, 0.50]:
        model_id = f"interneuron_6shared_twocohort-prior_{prior_sd:.2f}"
        print(f"[item7] Fitting {model_id}", flush=True)
        coords = {"gene": shared_genes, "obs": np.arange(len(sh_long))}
        with pm.Model(coords=coords) as model:
            g_d = pm.Data("g_idx", g_idx, dims="obs")
            se_d = pm.Data("se50", sh_long["se50"].to_numpy(dtype=float), dims="obs")
            mu_shared = pm.Normal("mu_shared", 0.0, float(prior_sd))
            tau_gene = pm.HalfNormal("tau_gene", float(prior_sd))
            z_g = pm.Normal("z_g", 0.0, 1.0, dims="gene")
            gene_mu = pm.Deterministic("gene_mu", mu_shared + tau_gene * z_g, dims="gene")
            tau_within = pm.HalfNormal("tau_within", 0.15)
            z_o = pm.Normal("z_o", 0.0, 1.0, dims="obs")
            theta_obs = pm.Deterministic("theta_obs", gene_mu[g_d] + tau_within * z_o, dims="obs")
            pm.Normal("y_obs", mu=theta_obs, sigma=se_d, observed=sh_long["y50"].to_numpy(dtype=float), dims="obs")
        idata, diag, psum = _sample_and_diagnose(
            model, sh_long["y50"].to_numpy(dtype=float), model_id, "interneurons", prior_sd, seed_offset=1050 + int(prior_sd * 100), target_accept=0.98
        )
        diag_rows.append(diag)
        prior_rows.append(psum)
        summary_rows.append(_summarize_array("interneurons", model_id, prior_sd, "shared6", "ALL_6_GENES", "mu_shared_50y", idata.posterior["mu_shared"].to_numpy(), len(shared_genes), len(sh_long)))
        for g in shared_genes:
            summary_rows.append(_summarize_array("interneurons", model_id, prior_sd, "shared6_gene", g, "gene_mu_50y", idata.posterior["gene_mu"].sel(gene=g).to_numpy(), 3, 6))


def _replace_model_rows(existing_path: Path, new_df: pd.DataFrame, model_ids: set[str]) -> pd.DataFrame:
    old = pd.read_csv(existing_path)
    kept = old[~old["model_id"].isin(model_ids)]
    return pd.concat([kept, new_df], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Item 7 Bayesian synthesis across five analytical families.")
    parser.add_argument(
        "--only-interneuron-eiv",
        action="store_true",
        help="Re-sample only the interneuron EIV models and update existing rev3-bayesian-full CSVs in place.",
    )
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    summary_rows: list[dict] = []
    diag_rows: list[dict] = []
    prior_rows: list[dict] = []

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if args.only_interneuron_eiv:
            run_family_interneurons(summary_rows, diag_rows, prior_rows, only_eiv=True)
        else:
            run_family_composition(summary_rows, diag_rows, prior_rows)
            run_family_transport(summary_rows, diag_rows, prior_rows)
            run_family_concordance(summary_rows, diag_rows, prior_rows)
            run_family_disorder_entropy(summary_rows, diag_rows, prior_rows)
            run_family_interneurons(summary_rows, diag_rows, prior_rows, only_eiv=False)

    sum_df = pd.DataFrame(summary_rows)
    diag_df = pd.DataFrame(diag_rows)
    prior_df = pd.DataFrame(prior_rows)
    if args.only_interneuron_eiv:
        eiv_ids = set(diag_df["model_id"])
        sum_df = _replace_model_rows(RESULTS / "rev3-bayesian-full-summary.csv", sum_df, eiv_ids)
        diag_df = _replace_model_rows(RESULTS / "rev3-bayesian-full-diagnostics.csv", diag_df, eiv_ids)
        prior_df = _replace_model_rows(RESULTS / "rev3-bayesian-full-prior-predictive.csv", prior_df, eiv_ids)

    sum_df.to_csv(RESULTS / "rev3-bayesian-full-summary.csv", index=False)
    diag_df.to_csv(RESULTS / "rev3-bayesian-full-diagnostics.csv", index=False)
    prior_df.to_csv(RESULTS / "rev3-bayesian-full-prior-predictive.csv", index=False)

    print("\n=== ITEM 7 DIAGNOSTICS SUMMARY ===", flush=True)
    print(f"Total models: {len(diag_df)}, converged: {int(diag_df['convergence_ok'].sum())}", flush=True)
    bad = diag_df[~diag_df["convergence_ok"]]
    if len(bad):
        print("Non-converged:", bad[["model_id", "max_rhat", "min_ess_bulk", "divergences"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
