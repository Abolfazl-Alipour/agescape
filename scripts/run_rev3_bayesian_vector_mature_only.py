#!/usr/bin/env python3
"""Item 6: Bayesian hierarchical whole-vector drift on mature mice (>= 18m vs 18m ref)."""

from __future__ import annotations

from pathlib import Path
import sys

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.bayesian_diagnostics_helper import (
    build_hierarchical_student_t,
    summarize_diagnostics,
)

RESULTS = REPO / "results"
SEED = sum(map(ord, "agescape-rev3-mature-18m-vector-bayes-v1"))


def prepare_mature_mouse_data(
    donors: pd.DataFrame,
    min_age_years: float,
    strict_gt: bool,
    endpoint_col: str,
    sex_filter: float | None = None,
) -> tuple[pd.DataFrame, dict, dict, bool]:
    part = donors[donors["table"] == "mouse_ref_18m"].copy()
    if strict_gt:
        part = part[part["age"] > min_age_years + 1e-6].copy()
    else:
        part = part[part["age"] >= min_age_years - 1e-6].copy()
    if sex_filter is not None:
        part = part[np.isclose(part["sex_male"], sex_filter)].copy()

    part["unit"] = part["tissue"] + " | " + part["cohort"] + " | " + part["cell_type"]
    part["group"] = part["tissue"]
    part["donor_cluster"] = part["cohort"] + " | " + part["donor"].astype(str)
    part = part.groupby("unit").filter(lambda frame: len(frame) >= 6).copy()

    means = part.groupby("unit")[endpoint_col].transform("mean")
    sds = part.groupby("unit")[endpoint_col].transform("std")
    part["y"] = (part[endpoint_col] - means) / sds
    part = part[np.isfinite(part["y"])].copy()
    part = part.groupby("unit").filter(lambda frame: len(frame) >= 6).copy()

    units = sorted(part["unit"].unique())
    groups = sorted(part["group"].unique())
    donor_levels = sorted(part["donor_cluster"].unique())
    unit_map = {v: i for i, v in enumerate(units)}
    group_map = {v: i for i, v in enumerate(groups)}
    donor_map = {v: i for i, v in enumerate(donor_levels)}

    # Scale in units of 10% of mouse 4-year max lifespan (0.4 y = 4.8 months), centered at min_age
    part["age_scaled"] = (part["age"] - float(min_age_years)) / 0.4
    part["sex_centered"] = part["sex_male"] - part.groupby("unit")["sex_male"].transform("mean")
    part["depth_centered"] = part["mean_log_umi"] - part.groupby("unit")["mean_log_umi"].transform("mean")
    part = part.sort_values(["unit", "donor_cluster"]).reset_index(drop=True)
    has_sex_var = bool(part["sex_centered"].abs().max() > 1e-9)

    coords = {"unit": units, "group": groups, "donor": donor_levels, "obs": part.index.to_numpy()}
    unit_groups = part.groupby("unit")["group"].first().loc[units]
    data = {
        "age": part["age_scaled"].to_numpy(dtype=float),
        "depth": part["depth_centered"].to_numpy(dtype=float),
        "unit_index": part["unit"].map(unit_map).to_numpy(dtype=int),
        "unit_group_index": unit_groups.map(group_map).to_numpy(dtype=int),
        "donor_index": part["donor_cluster"].map(donor_map).to_numpy(dtype=int),
        "y": part["y"].to_numpy(dtype=float),
    }
    if has_sex_var:
        data["sex"] = part["sex_centered"].to_numpy(dtype=float)
    return part, coords, data, has_sex_var


def frequentist_check(
    part: pd.DataFrame, depth_adjusted: bool, has_sex_var: bool, spec_id: str
) -> dict:
    dummies = pd.get_dummies(part["unit"], drop_first=True, dtype=float)
    preds = ["age_scaled"]
    if has_sex_var:
        preds.append("sex_centered")
    if depth_adjusted:
        preds.append("depth_centered")
    X = sm.add_constant(pd.concat([part[preds], dummies], axis=1), has_constant="add")
    fit = sm.OLS(part["y"], X).fit(
        cov_type="cluster", cov_kwds={"groups": part["donor_cluster"]}, use_t=True
    )
    ci = fit.conf_int().loc["age_scaled"]
    return {
        "spec_id": spec_id,
        "n_donors": int(part["donor_cluster"].nunique()),
        "n_tissues": int(part["group"].nunique()),
        "n_units": int(part["unit"].nunique()),
        "n_obs": int(len(part)),
        "freq_effect_10pct_lifespan": float(fit.params["age_scaled"]),
        "freq_se": float(fit.bse["age_scaled"]),
        "freq_ci95_low": float(ci.iloc[0]),
        "freq_ci95_high": float(ci.iloc[1]),
        "freq_pvalue": float(fit.pvalues["age_scaled"]),
    }


def extract_posterior_rows(
    idata, part: pd.DataFrame, coords: dict, spec_id: str, endpoint_col: str, subset_name: str, depth_adjusted: bool, prior_sd: float
) -> list[dict]:
    hdi = az.hdi(idata, hdi_prob=0.94)
    post = idata.posterior
    rows: list[dict] = []

    def add_row(scope: str, label: str, vals: np.ndarray, hdi_da) -> None:
        rows.append(
            {
                "spec_id": spec_id,
                "subset": subset_name,
                "endpoint": endpoint_col,
                "model": "depth_adjusted" if depth_adjusted else "unadjusted",
                "prior_sd": float(prior_sd),
                "scope": scope,
                "unit": label,
                "n_donors": int(part["donor_cluster"].nunique()),
                "n_tissues": int(part["group"].nunique()),
                "n_units": int(len(coords["unit"])),
                "n_obs": int(len(part)),
                "posterior_mean": float(vals.mean()),
                "posterior_sd": float(vals.std()),
                "hdi94_low": float(hdi_da.sel(hdi="lower")),
                "hdi94_high": float(hdi_da.sel(hdi="higher")),
                "prob_positive": float(np.mean(vals > 0.0)),
                "prob_gt_0_05": float(np.mean(vals > 0.05)),
                "prob_gt_0_10": float(np.mean(vals > 0.10)),
            }
        )

    add_row("species", "ALL", post["global_age"].to_numpy().ravel(), hdi["global_age"])
    if "global_sex" in post:
        add_row("covariate", "global_sex", post["global_sex"].to_numpy().ravel(), hdi["global_sex"])
    if "global_depth" in post:
        add_row("covariate", "global_depth", post["global_depth"].to_numpy().ravel(), hdi["global_depth"])
    if "group_age" in post:
        for grp in coords["group"]:
            vals = post["group_age"].sel(group=grp).to_numpy().ravel()
            add_row("tissue", grp, vals, hdi["group_age"].sel(group=grp))
    for u in coords["unit"]:
        vals = post["unit_age"].sel(unit=u).to_numpy().ravel()
        add_row("cell_type_unit", u, vals, hdi["unit_age"].sel(unit=u))
    return rows


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    donors = pd.read_csv(RESULTS / "reviewer-rev2-hellinger-donors.csv")

    specs = [
        ("mature_ge_18m", 1.5, False, None, "profile_hellinger_centroid_mean", False, 0.25),
        ("mature_ge_18m", 1.5, False, None, "profile_hellinger_centroid_mean", False, 0.50),
        ("mature_ge_18m", 1.5, False, None, "profile_hellinger_centroid_mean", True, 0.25),
        ("mature_ge_18m", 1.5, False, None, "profile_hellinger_centroid_mean", True, 0.50),
        ("post_ref_gt_18m", 1.5, True, None, "profile_hellinger_centroid_mean", False, 0.50),
        ("post_ref_gt_18m", 1.5, True, None, "profile_hellinger_centroid_mean", True, 0.50),
        ("male_ge_24m", 2.0, False, 1.0, "profile_hellinger_centroid_mean", False, 0.50),
        ("male_ge_24m", 2.0, False, 1.0, "profile_hellinger_centroid_mean", True, 0.50),
        ("mature_ge_18m", 1.5, False, None, "pseudobulk_hellinger", False, 0.50),
        ("mature_ge_18m", 1.5, False, None, "pseudobulk_hellinger", True, 0.50),
    ]

    post_rows: list[dict] = []
    diag_rows: list[dict] = []
    prior_rows: list[dict] = []
    freq_rows: list[dict] = []

    for subset_name, min_age, strict_gt, sex_filter, endpoint_col, depth_adj, prior_sd in specs:
        mod_label = "depth_adjusted" if depth_adj else "unadjusted"
        spec_id = f"{subset_name}-{endpoint_col}-{mod_label}-prior_{prior_sd:.2f}"
        print(f"[item6] Fitting {spec_id}", flush=True)
        part, coords, data, has_sex_var = prepare_mature_mouse_data(
            donors, min_age, strict_gt, endpoint_col, sex_filter=sex_filter
        )
        freq_rows.append(frequentist_check(part, depth_adj, has_sex_var, spec_id))

        rng = np.random.default_rng(SEED + sum(map(ord, spec_id)))
        model = build_hierarchical_student_t(
            coords,
            data,
            prior_sd=prior_sd,
            depth_adjusted=depth_adj,
            include_sex=has_sex_var,
        )
        with model:
            prior = pm.sample_prior_predictive(draws=500, random_seed=rng)
            idata = pm.sample(
                draws=2500,
                tune=2500,
                chains=4,
                nuts_sampler="nutpie",
                target_accept=0.99,
                random_seed=rng,
                progressbar=False,
            )
            idata.extend(prior)
            idata.extend(pm.sample_posterior_predictive(idata, random_seed=rng, progressbar=False))
            pm.compute_log_likelihood(idata, model=model)

        diag, prior_sum = summarize_diagnostics(
            idata, data["y"], model_id=spec_id, family="vector_mature_only", prior_sd=prior_sd
        )
        diag_rows.append(diag)
        prior_rows.append(prior_sum)
        rows_now = extract_posterior_rows(
            idata, part, coords, spec_id, endpoint_col, subset_name, depth_adj, prior_sd
        )
        post_rows.extend(rows_now)
        sp_row = [r for r in rows_now if r["scope"] == "species"][0]
        print(
            f"[item6] Done {spec_id}: mean={sp_row['posterior_mean']:+.4f} "
            f"[{sp_row['hdi94_low']:+.4f}, {sp_row['hdi94_high']:+.4f}] "
            f"P(>0)={sp_row['prob_positive']:.4f} | max_rhat={diag['max_rhat']:.4f} div={diag['divergences']}",
            flush=True,
        )

    post_df = pd.DataFrame(post_rows)
    freq_df = pd.DataFrame(freq_rows)
    post_df.to_csv(RESULTS / "rev3-bayesian-vector-mature-only-posterior.csv", index=False)
    pd.DataFrame(diag_rows).to_csv(RESULTS / "rev3-bayesian-vector-mature-only-diagnostics.csv", index=False)
    pd.DataFrame(prior_rows).to_csv(RESULTS / "rev3-bayesian-vector-mature-only-prior-predictive.csv", index=False)
    freq_df.to_csv(RESULTS / "rev3-bayesian-vector-mature-only-frequentist.csv", index=False)
    print("\n=== SPECIES SUMMARY ===", flush=True)
    print(post_df[post_df["scope"] == "species"].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
