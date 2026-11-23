#!/usr/bin/env python3
"""Bayesian refits for the mature Agescape analyses (followup-mature-bayesian).

Method (see ~/workspace/custom-skills/bayesian-testing/SKILL.md and the tex
Methods subsection "Bayesian refits"):
  - One Bayesian model per result, donor as the unit, never the cell.
  - Closed-form conjugate Bayesian linear regression (Normal-inverse-gamma)
    for every design-matrix model; exact, no MCMC. Summary-statistic models
    use a Normal likelihood with known SE plus a Normal prior (also exact).
  - Priors, fixed before fitting. Primary: Normal(0, 0.5) on the
    standardized age slope. Sensitivity: Normal(0, 1.0) and Normal(0, 0.25).
  - Standardization. Design-matrix models: the outcome is divided by its own
    SD and the age predictor by its SD, so the prior sits on the fully
    standardized slope (outcome SDs per predictor SD); posteriors are
    converted back to the raw slope (per year or per decade). Summary models
    have no design matrix, so the prior sits on the slope in outcome SDs per
    predictor unit (per decade); the source CSVs supply the outcome SD
    (Shannon rows carry outcome_sd = 1.0; the naive CD4 row carries
    outcome_sd_ref). The one exception is PsychAD Micro (item 4), which has
    no committed per-donor scores and therefore no outcome SD: its scale is
    borrowed from the MSSM microglia P1m score SD (same score family, same
    tissue, same cell type), with a half/double sensitivity reported.
  - Nuisance coefficients (intercept, sex, depth, cell-type dummies, batch
    fixed effects) carry a near-flat Normal(0, 1e4 * sigma2) prior; the
    residual variance carries the Jeffreys prior p(sigma2) ~ 1/sigma2.
  - Reporting: posterior median, 95% equal-tailed credible interval,
    probability in the claimed direction, ROPE probability. ROPE half-width
    is the slope that moves the outcome by less than 0.05 outcome SDs over
    the paper's reporting span for that slope: one decade for items 1-4
    (item 1 slopes are per year, so the ROPE is 0.05 * sd_y / 10 per year),
    0.4 years (mouse) / 12.25 years (human) for item 5, and 0.05 correlation
    points for item 6. Defined before fitting; ROPE sensitivity uses half
    and double width.
  - Seed 20261122; every headline number is rerun with seed 20261123 and the
    largest move is reported as a diagnostic.

Inputs: only committed CSVs and committed reports under results/ and reports/.
Outputs: per-result CSVs under results/ and reports/followup-mature-bayesian.md
(written by hand from the CSVs; this script prints the headline numbers).

DO NOT edit the tex, figures, or Q&A. DO NOT change any number in the paper.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
REPORTS = REPO / "reports"

SEED = 20261122
SEED2 = 20261123
N_DRAWS = 200_000
PRIORS = {"primary": 0.5, "wide": 1.0, "narrow": 0.25}
ROPE_HALF_SD = 0.05
ROPE_SENS = (0.025, 0.10)

rng_main = np.random.default_rng(SEED)


# ---------------------------------------------------------------- engine
def conjugate_lm_draws(X, y, slope_idx, tau, seed, n_draws=N_DRAWS):
    """Exact Normal-inverse-gamma posterior draws for a linear model.

    y = X beta + eps, eps ~ N(0, sigma2 I).
    beta | sigma2 ~ N(0, sigma2 * V0), V0 = diag(1e4) except V0[slope, slope]
    = tau^2. p(sigma2) ~ 1/sigma2. Returns (n_draws, k) beta draws.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, k = X.shape
    v0_inv = np.full(k, 1.0 / 1e4)
    v0_inv[slope_idx] = 1.0 / tau**2
    XtX = X.T @ X
    Vn = np.linalg.inv(XtX + np.diag(v0_inv))
    bn = Vn @ (X.T @ y)
    S = float(y @ y - bn @ np.linalg.inv(Vn) @ bn)
    assert S > 0, "non-positive posterior scale S"
    rng = np.random.default_rng(seed)
    sigma2 = 1.0 / rng.gamma(shape=n / 2.0, scale=2.0 / S, size=n_draws)
    L = np.linalg.cholesky(Vn)
    z = rng.standard_normal((n_draws, k))
    beta = bn[None, :] + np.sqrt(sigma2)[:, None] * (z @ L.T)
    return beta, bn, Vn, S


def summary_posterior_draws(est_std, se_std, tau, seed, n_draws=N_DRAWS):
    """Normal likelihood with known SE plus Normal(0, tau^2) prior."""
    post_var = 1.0 / (1.0 / tau**2 + 1.0 / se_std**2)
    post_mean = post_var * est_std / se_std**2
    rng = np.random.default_rng(seed)
    return rng.normal(post_mean, math.sqrt(post_var), size=n_draws)


def summarize(raw_draws, direction, rope_half):
    """raw_draws: posterior draws on the raw slope scale."""
    med = float(np.median(raw_draws))
    lo, hi = np.quantile(raw_draws, [0.025, 0.975])
    if direction == "negative":
        p_dir = float(np.mean(raw_draws < 0))
    elif direction == "positive":
        p_dir = float(np.mean(raw_draws > 0))
    else:
        raise ValueError(direction)
    p_rope = float(np.mean(np.abs(raw_draws) < rope_half))
    return {
        "median": med,
        "ci95_low": float(lo),
        "ci95_high": float(hi),
        "prob_direction": p_dir,
        "prob_rope": p_rope,
    }


def ols_slope(X, y, slope_idx):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n, k = X.shape
    s2 = float(resid @ resid) / (n - k)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = math.sqrt(cov[slope_idx, slope_idx])
    return float(beta[slope_idx]), se, float(np.sqrt(s2))


def cluster_se(X, y, groups, slope_idx):
    """Donor-clustered (CR1) standard error for one coefficient."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    XtX_inv = np.linalg.inv(X.T @ X)
    G = len(np.unique(groups))
    n, k = X.shape
    meat = np.zeros((k, k))
    for g in np.unique(groups):
        m = groups == g
        s = X[m].T @ resid[m]
        meat += np.outer(s, s)
    cov = XtX_inv @ meat @ XtX_inv * (G / (G - 1)) * ((n - 1) / (n - k))
    return float(beta[slope_idx]), math.sqrt(cov[slope_idx, slope_idx])


def design_fit(y, x_age, extra_cols, slope_name="age"):
    """Standardized conjugate fit. Returns dict with raw-scale draws per prior.

    y: raw outcome; x_age: raw age predictor (years or decades).
    extra_cols: list of (name, raw values) nuisance regressors.
    """
    sd_y = float(np.std(y, ddof=1))
    sd_x = float(np.std(x_age, ddof=1))
    assert sd_y > 0 and sd_x > 0
    y_s = np.asarray(y, dtype=float) / sd_y
    x_s = np.asarray(x_age, dtype=float) / sd_x
    cols = [np.ones(len(y)), x_s]
    names = ["intercept", slope_name]
    for nm, vals in extra_cols:
        v = np.asarray(vals, dtype=float)
        if np.unique(v).size >= 2:
            cols.append(v)
            names.append(nm)
    X = np.column_stack(cols)
    # drop exactly collinear columns (keep intercept and slope)
    keep = [0, 1]
    for j in range(2, X.shape[1]):
        if np.linalg.matrix_rank(X[:, keep + [j]]) == len(keep) + 1:
            keep.append(j)
    X = X[:, keep]
    names = [names[j] for j in keep]
    slope_idx = 1
    ols_b, ols_se, resid_sd = ols_slope(X, y_s, slope_idx)
    out = {}
    for label, tau in PRIORS.items():
        beta, _, _, _ = conjugate_lm_draws(X, y_s, slope_idx, tau, SEED)
        raw = beta[:, slope_idx] * sd_y / sd_x
        out[label] = raw
        beta2, _, _, _ = conjugate_lm_draws(X, y_s, slope_idx, tau, SEED2)
        out[label + "_seed2_median_raw"] = float(np.median(beta2[:, slope_idx] * sd_y / sd_x))
    out["sd_y"] = sd_y
    out["sd_x"] = sd_x
    out["ols_std"] = ols_b
    out["ols_se_std"] = ols_se
    out["ols_raw"] = ols_b * sd_y / sd_x
    out["ols_se_raw"] = ols_se * sd_y / sd_x
    out["resid_sd_std"] = resid_sd
    out["n"] = len(y)
    out["k"] = X.shape[1]
    out["design_names"] = names
    return out


def row_from_draws(result, unit, cohort, n, freq_est, freq_se, freq_p,
                   draws_by_prior, direction, rope_half, model, data_file,
                   extra=None):
    """Assemble one CSV row from raw-scale draws under the three priors."""
    r = {
        "result": result,
        "unit": unit,
        "cohort": cohort,
        "n": n,
        "model": model,
        "data_file": data_file,
        "freq_estimate": freq_est,
        "freq_se": freq_se,
        "freq_pvalue": freq_p,
        "direction": direction,
        "rope_half_raw": rope_half,
        "ci_type": "95% equal-tailed",
        "prior_sd_primary": 0.5,
    }
    for label, key in (("primary", "primary"), ("wide", "wide"), ("narrow", "narrow")):
        s = summarize(draws_by_prior[label], direction, rope_half)
        pfx = "" if label == "primary" else f"sens_{label}_"
        r[pfx + "posterior_median"] = s["median"]
        r[pfx + "ci95_low"] = s["ci95_low"]
        r[pfx + "ci95_high"] = s["ci95_high"]
        r[pfx + "prob_direction"] = s["prob_direction"]
        r[pfx + "prob_rope"] = s["prob_rope"]
    # ROPE-width sensitivity on the primary prior: half and double width.
    r["prob_rope_half_width"] = float(np.mean(np.abs(draws_by_prior["primary"]) < rope_half / 2.0))
    r["prob_rope_double_width"] = float(np.mean(np.abs(draws_by_prior["primary"]) < rope_half * 2.0))
    if extra:
        r.update(extra)
    return r


# ---------------------------------------------------------------- item 1
def item1_composition():
    print("=== Item 1: headline fraction shifts (CLR on age + sex) ===")
    df = pd.read_csv(RESULTS / "followup-mature-composition-donor-clr.csv")
    targets = [
        ("Human", "lung", "Lung_Integrated", "club cell", "negative",
         -0.08040728205, 0.01177686833),
        ("Human", "lung", "Lung_Integrated", "multiciliated epithelial cell", "negative",
         -0.04700185429, 0.02785398695),
        ("Human", "brain", "MSSM_DLPFC", "L6b glutamatergic neuron of the primary motor cortex",
         "positive", 0.01976532187, 0.01177686833),
        ("Human", "brain", "MSSM_DLPFC", "oligodendrocyte", "positive",
         0.01726917662, 0.01177686833),
        ("Human", "brain", "MSSM_DLPFC", "oligodendrocyte precursor cell", "negative",
         -0.01189556644, 0.01177686833),
    ]
    rows = []
    for species, organ, cohort, ct, direction, freq_slope, freq_q in targets:
        sub = df[(df.species == species) & (df.organ == organ)
                 & (df.cohort == cohort) & (df.cell_type == ct)].copy()
        assert len(sub) in (91, 209), (ct, len(sub))
        assert sub.donor_id.nunique() == len(sub), f"duplicate donors in {ct}"
        assert sub.sex_male.nunique() == 2, f"sex constant in {ct}"
        fit = design_fit(sub.clr_value.to_numpy(), sub.age.to_numpy(),
                         [("sex_male", sub.sex_male.to_numpy())])
        rope = ROPE_HALF_SD * fit["sd_y"] / 10.0  # 0.05 SD per decade, slope is per year
        draws = {k: fit[k] for k in ("primary", "wide", "narrow")}
        seed2_move = max(abs(fit["primary_seed2_median_raw"] - np.median(draws["primary"])),
                         abs(fit["wide_seed2_median_raw"] - np.median(draws["wide"])),
                         abs(fit["narrow_seed2_median_raw"] - np.median(draws["narrow"])))
        r = row_from_draws(
            "composition_clr_slope_per_year", ct, cohort, len(sub),
            freq_slope, None, freq_q, draws, direction, rope,
            "CLR ~ age_years + sex_male; conjugate Normal-inverse-gamma",
            "results/followup-mature-composition-donor-clr.csv",
            extra={"freq_q_bh": freq_q, "outcome_sd": fit["sd_y"],
                   "age_sd_years": fit["sd_x"], "ols_slope_raw": fit["ols_raw"],
                   "ols_se_raw": fit["ols_se_raw"],
                   "seed2_max_median_move_raw": seed2_move,
                   "n_design_cols": fit["k"]})
        rows.append(r)
        print(f"  {ct[:40]:42s} n={len(sub)} OLS={fit['ols_raw']:+.4f} "
              f"Bayes med={r['posterior_median']:+.4f} "
              f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] "
              f"P(dir)={r['prob_direction']:.4f} freq={freq_slope:+.4f} q={freq_q:.4g}")
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "followup-mature-bayesian-composition.csv", index=False)
    return out


# ---------------------------------------------------------------- item 2
def _stratified_scatter_fit(dataset):
    """Bayesian version of the frequentist stratified-global scatter model.

    Replicates scripts/run_followup_mature_fulldepth.py::fit_stratified_global_
    with_batch (stratified_global_primary): outcome scatter standardized
    WITHIN each cell type by that cell type's SD; sex and depth centered
    within cell type; cell-type dummies; study/pool fixed effects (batches
    with >= 2 donors; reference = ENCODE v4 (Snyder) for Gao, largest batch
    for Yazar); rank pruning identical in spirit to the original.
    Limitation, stated in the report: the conjugate fit treats donor x
    cell-type rows as independent given the design, while the frequentist
    fit clusters standard errors by donor. The clustered OLS SE is reported
    beside the posterior SD so the reader can see the gap.
    """
    import statsmodels.api as sm_mod  # noqa: F401  (kept for parity; not used)
    part1 = pd.read_csv(RESULTS / f"followup-mature-fulldepth-{dataset}-donor-part1.csv")
    d = part1[~part1["is_pooled"]].copy()
    d = d.dropna(subset=["scatter", "age_per_decade", "sex_male", "_depth",
                         "study_or_pool", "donor_id", "group"])
    # Eligible cell types: replicate the original rule exactly. A group is
    # eligible only if its per-cell-type PRIMARY scatter model
    # (age_per_decade + sex_male + _depth + study/pool FE, min residual df 5
    # Gao / 10 Yazar) returns yazar_pool_fe_status == "evaluated", i.e. the
    # batch fixed-effect block survived rank pruning with at least one dummy
    # kept. This is scripts/run_followup_between_subject_gao_yazar.py::
    # _prune_covariate_matrix + fit_ols_model_row, reimplemented here (the
    # module itself is not imported: it pulls h5py and the full pipeline).
    # Batches are counted on rows (>= 2 rows kept); rows in smaller batches
    # are dropped; reference is ENCODE v4 (Snyder) for Gao, largest batch
    # for Yazar. A no-batch fallback exists for Yazar in the original but
    # its status is never "evaluated", so it never confers eligibility.
    is_gao = dataset == "gao"
    min_df = 5 if is_gao else 10
    ref_lv = "ENCODE v4 (Snyder)" if is_gao else None
    elig = []
    for g, cpart in d.groupby("group"):
        if len(cpart) < 6 or cpart["age_per_decade"].nunique() < 2:
            continue
        work = cpart.copy()
        vc = work["study_or_pool"].value_counts()
        ge_b = vc[vc >= 2].index.tolist()
        work = work[work["study_or_pool"].isin(ge_b)].copy()
        vc2 = work["study_or_pool"].value_counts()
        if len(vc2) < 2:
            continue
        lv_sorted = [k for k, _ in sorted(vc2.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))]
        ref = ref_lv if (ref_lv is not None and ref_lv in lv_sorted) else lv_sorted[0]
        design = [np.ones(len(work)), work["age_per_decade"].to_numpy(dtype=float)]
        for cov in ("sex_male", "_depth"):
            vals = work[cov].to_numpy(dtype=float)
            if np.unique(vals).size < 2:
                continue
            if np.linalg.matrix_rank(np.column_stack([*design, vals])) == len(design) + 1:
                design.append(vals)
        kept_b = 0
        for lv in lv_sorted:
            if lv == ref:
                continue
            vals = (work["study_or_pool"].astype(str) == str(lv)).astype(float).to_numpy()
            if np.linalg.matrix_rank(np.column_stack([*design, vals])) == len(design) + 1:
                design.append(vals)
                kept_b += 1
        if kept_b == 0:
            continue
        if (len(work) - len(design)) < min_df:
            continue
        elig.append(g)
    d = d[d["group"].isin(elig)].copy()
    parts = []
    for g, cpart in d.groupby("group"):
        sd_y = float(cpart["scatter"].std(ddof=1))
        if not np.isfinite(sd_y) or sd_y <= 0:
            continue
        cp = cpart.copy()
        cp["y_std"] = (cp["scatter"] - cp["scatter"].mean()) / sd_y
        cp["sex_c"] = cp["sex_male"] - cp["sex_male"].mean()
        cp["depth_c"] = cp["_depth"] - cp["_depth"].mean()
        parts.append(cp)
    df_st = pd.concat(parts, ignore_index=True)
    # batch FE: batches with >= 2 unique donors; reference level
    d_batches = df_st[["donor_id", "study_or_pool"]].drop_duplicates()
    vc = d_batches["study_or_pool"].value_counts()
    ge_b = vc[vc >= 2].index.tolist()
    df_st = df_st[df_st["study_or_pool"].isin(ge_b)].copy()
    vc2 = df_st[["donor_id", "study_or_pool"]].drop_duplicates()["study_or_pool"].value_counts()
    lv_sorted = [k for k, _ in sorted(vc2.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))]
    ref_lv = "ENCODE v4 (Snyder)" if (dataset == "gao" and "ENCODE v4 (Snyder)" in lv_sorted) else lv_sorted[0]
    batch_cols = []
    for lv in lv_sorted:
        if lv != ref_lv:
            cname = f"batch_{str(lv).replace(' ', '_')}"
            df_st[cname] = (df_st["study_or_pool"].astype(str) == str(lv)).astype(float)
            batch_cols.append(cname)
    g_dum = pd.get_dummies(df_st["group"], prefix="g", drop_first=True, dtype=float)
    y = df_st["y_std"].to_numpy(dtype=float)
    sd_x = float(np.std(df_st["age_per_decade"].to_numpy(dtype=float), ddof=1))
    x_s = df_st["age_per_decade"].to_numpy(dtype=float) / sd_x
    cols = [np.ones(len(df_st)), x_s, df_st["sex_c"].to_numpy(dtype=float),
            df_st["depth_c"].to_numpy(dtype=float)]
    names = ["intercept", "age_per_decade", "sex_c", "depth_c"]
    for c in g_dum.columns:
        cols.append(g_dum[c].to_numpy(dtype=float))
        names.append(c)
    for c in batch_cols:
        cols.append(df_st[c].to_numpy(dtype=float))
        names.append(c)
    X = np.column_stack(cols)
    keep = list(range(X.shape[1]))
    # rank pruning, keeping intercept and slope
    pruned = [0, 1]
    for j in range(2, X.shape[1]):
        if np.linalg.matrix_rank(X[:, pruned + [j]]) == len(pruned) + 1:
            pruned.append(j)
    X = X[:, pruned]
    names = [names[j] for j in pruned]
    slope_idx = 1
    ols_b, ols_se, resid_sd = ols_slope(X, y, slope_idx)
    _, cl_se = cluster_se(X, y, df_st["donor_id"].to_numpy(), slope_idx)
    draws = {}
    seed2_move = 0.0
    for label, tau in PRIORS.items():
        beta, _, _, _ = conjugate_lm_draws(X, y, slope_idx, tau, SEED)
        draws[label] = beta[:, slope_idx] / sd_x  # y_std units per decade
        beta2, _, _, _ = conjugate_lm_draws(X, y, slope_idx, tau, SEED2)
        seed2_move = max(seed2_move, abs(float(np.median(beta2[:, slope_idx] / sd_x))
                                            - float(np.median(draws[label]))))
    return {
        "draws": draws, "sd_x_decades": sd_x, "ols_std": ols_b / sd_x,
        "ols_se_std": ols_se / sd_x, "cluster_se_std": cl_se / sd_x,
        "n_rows": len(df_st), "n_donors": int(df_st["donor_id"].nunique()),
        "n_cell_types": len(elig), "n_batches_kept": len(batch_cols),
        "ref_batch": ref_lv, "seed2_max_median_move": seed2_move,
        "design_names": names, "resid_sd_std": resid_sd,
    }


def item2_narrowing():
    print("=== Item 2: narrowing (scatter + Shannon) ===")
    rows = []
    refs = {"gao": (-0.0722825, 0.04103), "yazar": (-0.0196652, 6.685e-05)}
    for dataset in ("gao", "yazar"):
        fit = _stratified_scatter_fit(dataset)
        freq_slope, freq_p = refs[dataset]
        rope = ROPE_HALF_SD  # y_std units per decade
        r = row_from_draws(
            "scatter_slope_per_decade", f"{dataset} stratified scatter",
            dataset, fit["n_rows"], freq_slope, fit["cluster_se_std"], freq_p,
            fit["draws"], "negative", rope,
            "scatter_std(within cell type) ~ age_per_decade + sex_c + depth_c "
            "+ cell-type dummies + study/pool FE; conjugate Normal-inverse-gamma",
            f"results/followup-mature-fulldepth-{dataset}-donor-part1.csv",
            extra={"n_donors": fit["n_donors"], "n_cell_types": fit["n_cell_types"],
                   "n_batches_kept": fit["n_batches_kept"],
                   "ref_batch": fit["ref_batch"],
                   "age_sd_decades": fit["sd_x_decades"],
                   "ols_slope": fit["ols_std"], "ols_se": fit["ols_se_std"],
                   "ols_cluster_se": fit["cluster_se_std"],
                   "seed2_max_median_move": fit["seed2_max_median_move"]})
        rows.append(r)
        print(f"  {dataset} scatter n_rows={fit['n_rows']} n_donors={fit['n_donors']} "
              f"OLS={fit['ols_std']:+.4f} (cluster SE {fit['cluster_se_std']:.4f}) "
              f"Bayes med={r['posterior_median']:+.4f} "
              f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] "
              f"P(<0)={r['prob_direction']:.4f} freq={freq_slope:+.4f} p={freq_p:.3g}")
    # Shannon: summary-statistic model on the Stratified_CellType_Global rows.
    # model_spec in the CSVs is 'stratified_global_primary' (the task text
    # says 'primary'); full_depth is the primary depth mode. outcome_sd = 1.0
    # in these rows, so the slope is already in outcome-SD units per decade.
    shannon_specs = [
        ("gao", "stratified_global_primary", "full_depth", -0.100142, 0.051363, 0.053904),
        ("gao", "stratified_global_depth_x_ct", "full_depth", -0.101305, 0.047263, 0.034410),
        ("yazar", "stratified_global_primary", "full_depth", -0.029773, 0.010216, 0.003649),
        ("yazar", "stratified_global_depth_x_ct", "full_depth", -0.027407, 0.010141, 0.007006),
    ]
    spread = {}
    for dataset in ("gao", "yazar"):
        sp = pd.read_csv(RESULTS / f"followup-mature-fulldepth-{dataset}-spread.csv")
        spread[dataset] = sp
    for dataset, spec, depth_mode, slope, se, p in shannon_specs:
        sp = spread[dataset]
        match = sp[(sp.metric == "shannon") & (sp.model_spec == spec)
                   & (sp.depth_mode == depth_mode)
                   & (sp.group == "Stratified_CellType_Global")]
        assert len(match) == 1, (dataset, spec, depth_mode, len(match))
        m = match.iloc[0]
        assert abs(float(m.slope_raw_per_decade) - slope) < 1e-6
        assert abs(float(m.se_raw) - se) < 1e-6
        assert abs(float(m.pvalue) - p) < 1e-6
        assert float(m.outcome_sd) == 1.0
        draws = {}
        seed2_move = 0.0
        for label, tau in PRIORS.items():
            draws[label] = summary_posterior_draws(slope, se, tau, SEED)
            d2 = summary_posterior_draws(slope, se, tau, SEED2)
            seed2_move = max(seed2_move, abs(float(np.median(d2)) - float(np.median(draws[label]))))
        r = row_from_draws(
            "shannon_slope_per_decade", f"{dataset} Shannon {spec}",
            dataset, int(m.n_donors), slope, se, p, draws, "negative",
            ROPE_HALF_SD,
            "Normal likelihood with known SE + Normal prior on the slope",
            f"results/followup-mature-fulldepth-{dataset}-spread.csv",
            extra={"model_spec": spec, "depth_mode": depth_mode,
                   "seed2_max_median_move": seed2_move})
        rows.append(r)
        print(f"  {dataset} Shannon {spec}: Bayes med={r['posterior_median']:+.4f} "
              f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] "
              f"P(<0)={r['prob_direction']:.4f} freq={slope:+.4f} p={p:.3g}")
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "followup-mature-bayesian-narrowing.csv", index=False)
    return out


# ---------------------------------------------------------------- item 3
def item3_naive_cd4():
    print("=== Item 3: naive CD4 dispersion ===")
    df = pd.read_csv(RESULTS / "followup-mature-block1-residual-dispersion.csv")
    r0 = df[(df.dataset == "yazar")
            & (df.cell_type == "naive thymus-derived CD4-positive, alpha-beta T cell")
            & (df.arm == "primary") & (df.is_pooled == False)]
    assert len(r0) == 1
    m = r0.iloc[0]
    slope, se = float(m.slope_per_decade), float(m.se)
    p, q = float(m.pvalue), float(m.q_bh)
    sd_ref = float(m.outcome_sd_ref)
    assert abs(slope - 0.000714) < 1e-6 and abs(p - 5.174179e-08) < 1e-12
    est_std, se_std = slope / sd_ref, se / sd_ref
    rows = []
    draws = {}
    seed2_move = 0.0
    for label, tau in PRIORS.items():
        d_std = summary_posterior_draws(est_std, se_std, tau, SEED)
        draws[label] = d_std * sd_ref
        d2 = summary_posterior_draws(est_std, se_std, tau, SEED2)
        seed2_move = max(seed2_move, abs(float(np.median(d2 * sd_ref)) - float(np.median(draws[label]))))
    rows.append(row_from_draws(
        "naive_cd4_dispersion_slope_per_decade", "yazar naive CD4 T (primary)",
        "yazar", int(m.n_donors), slope, se, p, draws, "positive",
        ROPE_HALF_SD * sd_ref,
        "Normal likelihood with known SE + Normal prior on the slope",
        "results/followup-mature-block1-residual-dispersion.csv",
        extra={"freq_q_bh": q, "outcome_sd_ref": sd_ref,
               "slope_sd_per_decade": float(m.slope_sd_per_decade),
               "seed2_max_median_move": seed2_move}))
    # Debiased slope: exists ONLY in reports/followup-mature-block1.md §5.7
    # (0.000313, p = 0.000638). No SE is committed anywhere, so the SE is
    # derived here as slope / z, z = Phi^-1(1 - p/2) for the two-sided p.
    deb_slope, deb_p = 0.000313, 0.000638
    deb_z = float(norm.ppf(1.0 - deb_p / 2.0))
    deb_se = deb_slope / deb_z
    deb_est_std, deb_se_std = deb_slope / sd_ref, deb_se / sd_ref
    draws_d = {}
    seed2_move_d = 0.0
    for label, tau in PRIORS.items():
        d_std = summary_posterior_draws(deb_est_std, deb_se_std, tau, SEED)
        draws_d[label] = d_std * sd_ref
        d2 = summary_posterior_draws(deb_est_std, deb_se_std, tau, SEED2)
        seed2_move_d = max(seed2_move_d, abs(float(np.median(d2 * sd_ref)) - float(np.median(draws_d[label]))))
    rows.append(row_from_draws(
        "naive_cd4_dispersion_slope_per_decade", "yazar naive CD4 T (debiased)",
        "yazar", int(m.n_donors), deb_slope, deb_se, deb_p, draws_d, "positive",
        ROPE_HALF_SD * sd_ref,
        "Normal likelihood with SE derived as slope/z from the two-sided p "
        "(no SE committed); Normal prior on the slope",
        "reports/followup-mature-block1.md §5.7 (slope and p only; SE derived)",
        extra={"freq_q_bh": None, "outcome_sd_ref": sd_ref,
               "derived_z": deb_z, "derived_se_note": "SE = slope / z, z = Phi^-1(1 - p/2)",
               "seed2_max_median_move": seed2_move_d}))
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "followup-mature-bayesian-naive-cd4.csv", index=False)
    for r in rows:
        print(f"  {r['unit']}: Bayes med={r['posterior_median']:.6f} "
              f"[{r['ci95_low']:.6f}, {r['ci95_high']:.6f}] "
              f"P(>0)={r['prob_direction']:.4f} freq={r['freq_estimate']:.6f} "
              f"p={r['freq_pvalue']:.3g}")
    return out


# ---------------------------------------------------------------- item 4
def item4_transport():
    print("=== Item 4: transport null (P1m composite scores) ===")
    # NOTE on inputs. The task names results/followup-mature-block8-donor-scores.csv,
    # but that file holds ONLY gao and yazar rows: no PsychAD (DLPFC_Lifespan_Atlas)
    # and no MSSM rows exist in it, so the two microglia passes cannot be fit
    # from it. The per-donor P1m scores for the core cohorts live in
    # results/followup-mature-p1m-m9-donor-scores.csv (MSSM, HBCC, kidney, lung,
    # BAL). PsychAD per-donor scores are committed nowhere: the PsychAD Micro
    # pass is fit below as a summary-statistic model from its committed slope
    # and p value, with the SE derived as slope/z. The block8 file is still
    # used for two of the ten null units (one gao, one yazar).
    m9 = pd.read_csv(RESULTS / "followup-mature-p1m-m9-donor-scores.csv")
    m9["sex_male"] = (m9["sex"] == "male").astype(float)
    m9["age_per_decade"] = m9["age"] / 10.0
    transport = pd.read_csv(RESULTS / "followup-mature-p1m-transport.csv")
    p1m = transport[transport.gene_set == "P1m"].copy()

    def freq_ref(cohort, unit):
        r = p1m[(p1m.cohort == cohort) & (p1m.unit == unit)]
        assert len(r) == 1, (cohort, unit, len(r))
        m = r.iloc[0]
        return (float(m.composite_slope) * 10.0, float(m.composite_pvalue),
                float(m.composite_q), int(m.n_donors))

    rows = []
    # --- the two microglia passes (negative direction) ---
    # MSSM microglial cell: per-donor fit.
    sub = m9[(m9.cohort == "MSSM_DLPFC") & (m9.cell_type == "microglial cell")].copy()
    assert len(sub) == 196, len(sub)
    fit = design_fit(sub.score_P1m.to_numpy(), sub.age_per_decade.to_numpy(),
                     [("sex_male", sub.sex_male.to_numpy())])
    f_slope, f_p, f_q, _ = freq_ref("MSSM_DLPFC", "microglial cell")
    rope = ROPE_HALF_SD * fit["sd_y"]
    draws = {k: fit[k] for k in ("primary", "wide", "narrow")}
    seed2_move = max(abs(fit[f"{k}_seed2_median_raw"] - float(np.median(draws[k])))
                     for k in ("primary", "wide", "narrow"))
    rows.append(row_from_draws(
        "transport_composite_slope_per_decade", "MSSM_DLPFC microglial cell",
        "MSSM_DLPFC", len(sub), f_slope, None, f_p, draws, "negative", rope,
        "score_P1m ~ age_per_decade + sex_male; conjugate Normal-inverse-gamma",
        "results/followup-mature-p1m-m9-donor-scores.csv",
        extra={"freq_q_bh": f_q, "score_sd": fit["sd_y"],
               "age_sd_decades": fit["sd_x"], "ols_slope_raw": fit["ols_raw"],
               "ols_se_raw": fit["ols_se_raw"],
               "seed2_max_median_move_raw": seed2_move}))
    print(f"  MSSM microglia n={len(sub)} OLS={fit['ols_raw']:+.5f}/dec "
          f"Bayes med={rows[-1]['posterior_median']:+.5f} "
          f"[{rows[-1]['ci95_low']:+.5f}, {rows[-1]['ci95_high']:+.5f}] "
          f"P(<0)={rows[-1]['prob_direction']:.4f} freq={f_slope:+.5f} p={f_p:.3g}")
    # PsychAD Micro: summary-statistic model. Committed slope -0.001433/yr
    # (fully_adjusted model), p = 8.06e-05, n = 201. No SE is committed, so
    # SE = |slope| / z with z = Phi^-1(1 - p/2). No per-donor PsychAD scores
    # are committed, so the score SD is borrowed from the MSSM microglia
    # P1m scores (same score family, same tissue, same cell type); the
    # half/double-scale sensitivity is reported in the extra columns.
    s0 = float(sub.score_P1m.std(ddof=1))
    psy_slope_yr, psy_p = -0.001433, 8.06e-05
    psy_z = float(norm.ppf(1.0 - psy_p / 2.0))
    psy_se_yr = abs(psy_slope_yr) / psy_z
    psy_slope, psy_se = psy_slope_yr * 10.0, psy_se_yr * 10.0
    est_std, se_std = psy_slope / s0, psy_se / s0
    draws_p = {}
    seed2_move_p = 0.0
    for label, tau in PRIORS.items():
        d_std = summary_posterior_draws(est_std, se_std, tau, SEED)
        draws_p[label] = d_std * s0
        d2 = summary_posterior_draws(est_std, se_std, tau, SEED2)
        seed2_move_p = max(seed2_move_p, abs(float(np.median(d2 * s0)) - float(np.median(draws_p[label]))))
    # borrowed-scale sensitivity: P(<0) under s0/2 and 2*s0, primary prior
    sens_scale = {}
    for mult in (0.5, 2.0):
        d_std = summary_posterior_draws(psy_slope / (s0 * mult), psy_se / (s0 * mult), 0.5, SEED)
        sens_scale[f"prob_neg_scale_x{mult}"] = float(np.mean(d_std * (s0 * mult) < 0))
    rows.append(row_from_draws(
        "transport_composite_slope_per_decade", "DLPFC_Lifespan_Atlas Micro (PsychAD)",
        "DLPFC_Lifespan_Atlas", 201, psy_slope, psy_se, psy_p, draws_p, "negative",
        ROPE_HALF_SD * s0,
        "Normal likelihood with SE derived as |slope|/z from the two-sided p; "
        "score SD borrowed from MSSM microglia P1m scores (no per-donor "
        "PsychAD scores committed); Normal prior on the slope",
        "results/followup-mature-p1m-transport.csv (slope, p); "
        "results/followup-mature-p1m-m9-donor-scores.csv (borrowed scale)",
        extra={"freq_q_bh": 0.004568, "score_sd_borrowed": s0,
               "derived_z": psy_z,
               "derived_se_note": "SE = |slope| / z, z = Phi^-1(1 - p/2)",
               "seed2_max_median_move_raw": seed2_move_p, **sens_scale}))
    print(f"  PsychAD Micro (summary) Bayes med={rows[-1]['posterior_median']:+.5f}/dec "
          f"[{rows[-1]['ci95_low']:+.5f}, {rows[-1]['ci95_high']:+.5f}] "
          f"P(<0)={rows[-1]['prob_direction']:.4f} freq={psy_slope:+.5f} p={psy_p:.3g}")

    # --- ten null units: 8 from the P1m score family + 2 from block8 ---
    # Ten null units: 8 from the P1m score family + 2 from block8. Only
    # MSSM, HBCC and Lung_Integrated have units with >= 40 donors in the
    # per-donor score file (kidney cohorts max 19 donors, Lung_CellRef max
    # 10), so the spread across cohorts is 3 MSSM + 3 HBCC + 2 Lung.
    cand = p1m[(p1m.model == "age_plus_sex") & (p1m.composite_pvalue > 0.20)
               & (p1m.n_donors >= 40)].copy()
    cohort_order = ["MSSM_DLPFC", "HBCC_DLPFC", "Lung_Integrated",
                    "Kidney_SCR_v15", "Kidney_LivingDonor", "Lung_CellRef"]
    picked = []
    for cohort in cohort_order:
        if len(picked) >= 8:
            break
        taken_here = 0
        cc = cand[cand.cohort == cohort].sort_values("unit")
        for _, r0 in cc.iterrows():
            if len(picked) >= 8 or taken_here >= 3:
                break
            n_m9 = int(((m9.cohort == cohort) & (m9.cell_type == r0.unit)).sum())
            if n_m9 >= 40:
                picked.append((cohort, str(r0.unit)))
                taken_here += 1
    assert len(picked) == 8, picked
    for cohort, unit in picked:
        sub = m9[(m9.cohort == cohort) & (m9.cell_type == unit)].copy()
        fit = design_fit(sub.score_P1m.to_numpy(), sub.age_per_decade.to_numpy(),
                         [("sex_male", sub.sex_male.to_numpy())])
        f_slope, f_p, f_q, _ = freq_ref(cohort, unit)
        direction = "negative" if fit["ols_raw"] < 0 else "positive"
        rope = ROPE_HALF_SD * fit["sd_y"]
        draws = {k: fit[k] for k in ("primary", "wide", "narrow")}
        seed2_move = max(abs(fit[f"{k}_seed2_median_raw"] - float(np.median(draws[k])))
                         for k in ("primary", "wide", "narrow"))
        rows.append(row_from_draws(
            "transport_composite_slope_per_decade", f"{cohort} {unit}",
            cohort, len(sub), f_slope, None, f_p, draws, direction, rope,
            "score_P1m ~ age_per_decade + sex_male; conjugate Normal-inverse-gamma",
            "results/followup-mature-p1m-m9-donor-scores.csv",
            extra={"freq_q_bh": f_q, "score_sd": fit["sd_y"],
                   "age_sd_decades": fit["sd_x"], "ols_slope_raw": fit["ols_raw"],
                   "ols_se_raw": fit["ols_se_raw"],
                   "seed2_max_median_move_raw": seed2_move}))
        print(f"  NULL {cohort} {unit[:38]:38s} n={len(sub)} OLS={fit['ols_raw']:+.5f} "
              f"Bayes med={rows[-1]['posterior_median']:+.5f} "
              f"[{rows[-1]['ci95_low']:+.5f}, {rows[-1]['ci95_high']:+.5f}] "
              f"P(dir)={rows[-1]['prob_direction']:.3f} P(ROPE)={rows[-1]['prob_rope']:.3f} "
              f"freq p={f_p:.2f}")
    # two block8 nulls (gao + yazar), selected by in-script OLS p > 0.20, n >= 40
    b8 = pd.read_csv(RESULTS / "followup-mature-block8-donor-scores.csv")
    for dataset in ("gao", "yazar"):
        d8 = b8[b8.dataset == dataset]
        quals = []
        for ct, grp in d8.groupby("cell_type"):
            if len(grp) < 40:
                continue
            X = sm.add_constant(pd.DataFrame(
                {"age_per_decade": grp.age_per_decade.to_numpy(),
                 "sex_male": grp.sex_male.to_numpy()}), has_constant="add")
            fit8 = sm.OLS(grp.composite_score.to_numpy(), X).fit()
            if float(fit8.pvalues["age_per_decade"]) > 0.20:
                quals.append((ct, float(fit8.pvalues["age_per_decade"]), len(grp)))
        quals.sort()
        assert quals, dataset
        ct = quals[0][0]
        sub = d8[d8.cell_type == ct].copy()
        fit = design_fit(sub.composite_score.to_numpy(), sub.age_per_decade.to_numpy(),
                         [("sex_male", sub.sex_male.to_numpy())])
        direction = "negative" if fit["ols_raw"] < 0 else "positive"
        rope = ROPE_HALF_SD * fit["sd_y"]
        draws = {k: fit[k] for k in ("primary", "wide", "narrow")}
        seed2_move = max(abs(fit[f"{k}_seed2_median_raw"] - float(np.median(draws[k])))
                         for k in ("primary", "wide", "narrow"))
        rows.append(row_from_draws(
            "transport_composite_slope_per_decade", f"block8 {dataset} {ct}",
            dataset, len(sub), fit["ols_raw"], fit["ols_se_raw"], quals[0][1],
            draws, direction, rope,
            "composite_score ~ age_per_decade + sex_male; conjugate Normal-inverse-gamma",
            "results/followup-mature-block8-donor-scores.csv",
            extra={"freq_q_bh": None, "score_sd": fit["sd_y"],
                   "age_sd_decades": fit["sd_x"], "ols_slope_raw": fit["ols_raw"],
                   "ols_se_raw": fit["ols_se_raw"],
                   "seed2_max_median_move_raw": seed2_move,
                   "freq_note": "no slope reference committed for block8 units; "
                                "frequentist columns are the in-script OLS"}))
        print(f"  NULL block8 {dataset} {ct[:38]:38s} n={len(sub)} OLS={fit['ols_raw']:+.6f} "
              f"Bayes med={rows[-1]['posterior_median']:+.6f} "
              f"P(dir)={rows[-1]['prob_direction']:.3f} P(ROPE)={rows[-1]['prob_rope']:.3f}")
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "followup-mature-bayesian-transport.csv", index=False)
    return out


# ---------------------------------------------------------------- item 5
def item5_drift_mouse():
    print("=== Item 5: drift nulls (mouse carry-over; human = PyMC on beast) ===")
    post = pd.read_csv(RESULTS / "rev3-bayesian-vector-mature-only-posterior.csv")
    row = post[(post.spec_id == "mature_ge_18m-profile_hellinger_centroid_mean-depth_adjusted-prior_0.50")
               & (post.scope == "species")]
    assert len(row) == 1
    m = row.iloc[0]
    assert abs(float(m.posterior_mean) - 0.038944) < 1e-6
    assert abs(float(m.hdi94_low) - (-0.099779)) < 1e-6
    assert abs(float(m.hdi94_high) - 0.177334) < 1e-6
    assert abs(float(m.prob_positive) - 0.7068) < 1e-4
    mean, sd = float(m.posterior_mean), float(m.posterior_sd)
    # prior-0.25 sensitivity row (committed); no prior-1.00 mature row exists.
    row25 = post[(post.spec_id == "mature_ge_18m-profile_hellinger_centroid_mean-depth_adjusted-prior_0.25")
                 & (post.scope == "species")]
    assert len(row25) == 1
    m25 = row25.iloc[0]
    rope = ROPE_HALF_SD  # standardized Hellinger SDs per 0.4 years
    p_rope = float(norm.cdf((rope - mean) / sd) - norm.cdf((-rope - mean) / sd))
    p_rope_half = float(norm.cdf((rope / 2 - mean) / sd) - norm.cdf((-rope / 2 - mean) / sd))
    p_rope_double = float(norm.cdf((rope * 2 - mean) / sd) - norm.cdf((-rope * 2 - mean) / sd))
    p_rope25 = float(norm.cdf((rope - float(m25.posterior_mean)) / float(m25.posterior_sd))
                     - norm.cdf((-rope - float(m25.posterior_mean)) / float(m25.posterior_sd)))
    r = {
        "result": "drift_slope_per_0.4y", "unit": "mouse global (mature >= 18m)",
        "cohort": "TMS", "n": int(m.n_donors), "model": "hierarchical Student-t (committed rev3 refit)",
        "data_file": "results/rev3-bayesian-vector-mature-only-posterior.csv",
        "freq_estimate": -0.013715, "freq_se": 0.057853, "freq_pvalue": 0.816303,
        "direction": "positive", "rope_half_raw": rope,
        "ci_type": "94% HDI (as committed; 95% equal-tailed normal approx also given)",
        "prior_sd_primary": 0.5,
        "posterior_median": mean, "ci95_low": mean - 1.96 * sd,
        "ci95_high": mean + 1.96 * sd, "prob_direction": float(m.prob_positive),
        "prob_rope": p_rope,
        "prob_rope_half_width": p_rope_half,
        "prob_rope_double_width": p_rope_double,
        "hdi94_low_committed": float(m.hdi94_low),
        "hdi94_high_committed": float(m.hdi94_high),
        "sens_narrow_posterior_median": float(m25.posterior_mean),
        "sens_narrow_ci95_low": float(m25.posterior_mean) - 1.96 * float(m25.posterior_sd),
        "sens_narrow_ci95_high": float(m25.posterior_mean) + 1.96 * float(m25.posterior_sd),
        "sens_narrow_prob_direction": float(m25.prob_positive),
        "sens_narrow_prob_rope": p_rope25,
        "sens_wide_posterior_median": None, "sens_wide_ci95_low": None,
        "sens_wide_ci95_high": None, "sens_wide_prob_direction": None,
        "sens_wide_prob_rope": None,
        "sens_wide_note": "no prior-1.00 mature mouse refit is committed",
        "n_units": int(m.n_units), "n_obs": int(m.n_obs),
    }
    out = pd.DataFrame([r])
    out.to_csv(RESULTS / "followup-mature-bayesian-drift.csv", index=False)
    print(f"  mouse depth-adjusted: mean={mean:+.4f} 94% HDI "
          f"[{float(m.hdi94_low):+.4f}, {float(m.hdi94_high):+.4f}] "
          f"P(>0)={float(m.prob_positive):.4f} P(ROPE)={p_rope:.3f} "
          f"(freq -0.0137, p=0.816)")
    print(f"  prior-0.25: mean={float(m25.posterior_mean):+.4f} "
          f"P(>0)={float(m25.prob_positive):.4f} P(ROPE)={p_rope25:.3f}")
    return out


# ---------------------------------------------------------------- item 6
def item6_concordance():
    print("=== Item 6: mouse within-tissue concordance ===")
    df = pd.read_csv(RESULTS / "followup-mature-concordance-mouse.csv")
    row = df[(df.species == "Mouse") & (df.scope == "within_tissue")
             & (df.gene_set == "primary_detection_filtered")]
    assert len(row) == 1
    m = row.iloc[0]
    r_obs = float(m.observed_mean_pearson)
    r_null, sd_null = float(m.joint_null_mean), float(m.joint_null_sd)
    n_pairs, n_perm = int(m.n_pairs), int(m.n_permutations_joint)
    p_freq = float(m.joint_p_greater)
    assert abs(r_obs - 0.21047) < 1e-4 and abs(r_null - 0.14646) < 1e-4
    assert abs(sd_null - 0.02724) < 1e-4 and n_pairs == 260 and n_perm == 1000
    # Construction, exactly. Fisher z of the observed correlation carries
    # sampling error 1/sqrt(n_pairs - 3) under a flat prior on z. The null
    # correlation is Normal(joint_null_mean, joint_null_sd) from the 1000
    # joint permutations; on the z scale its SE is sd_null / (1 - r_null^2)
    # (delta method). The estimand is delta_z = z_obs - z_null, with a
    # Normal(0, tau^2) prior on delta_z. delta_r is then recovered by
    # drawing z_null from its Normal and back-transforming:
    # delta_r = tanh(z_null + delta_z) - tanh(z_null).
    z_obs = float(np.arctanh(r_obs))
    z_null = float(np.arctanh(r_null))
    se_obs_z = 1.0 / math.sqrt(n_pairs - 3)
    se_null_z = sd_null / (1.0 - r_null**2)
    d_hat = z_obs - z_null
    se_d = math.sqrt(se_obs_z**2 + se_null_z**2)
    rng = np.random.default_rng(SEED)
    z_null_draws = rng.normal(z_null, se_null_z, size=N_DRAWS)
    rows = []
    draws_by_prior = {}
    for label, tau in PRIORS.items():
        d_draws = summary_posterior_draws(d_hat, se_d, tau, SEED)
        draws_by_prior[label] = np.tanh(z_null_draws + d_draws) - np.tanh(z_null_draws)
    # direction: positive difference (observed above the joint null)
    r = row_from_draws(
        "concordance_delta_r", "Mouse within_tissue primary_detection_filtered",
        "Mouse", n_pairs, r_obs - r_null, se_d, p_freq, draws_by_prior,
        "positive", 0.05,
        "Fisher-z difference: z_obs ~ N(atanh(r_obs), 1/(n_pairs-3)), "
        "z_null ~ N(atanh(joint_null_mean), (sd_null/(1-r_null^2))^2); "
        "Normal prior on delta_z; delta_r by back-transform",
        "results/followup-mature-concordance-mouse.csv",
        extra={"r_observed": r_obs, "r_null_mean": r_null, "r_null_sd": sd_null,
               "n_permutations_joint": n_perm,
               "delta_z_hat": d_hat, "se_delta_z": se_d,
               "freq_p_greater": p_freq, "freq_q_bh_two_test": 0.0440})
    # also report the delta_z summaries in the extra columns
    for label, tau in PRIORS.items():
        d_draws = summary_posterior_draws(d_hat, se_d, tau, SEED)
        s = summarize(d_draws, "positive", 0.05)
        pfx = "" if label == "primary" else f"sens_{label}_"
        r[pfx + "delta_z_median"] = s["median"]
        r[pfx + "delta_z_ci95_low"] = s["ci95_low"]
        r[pfx + "delta_z_ci95_high"] = s["ci95_high"]
        r[pfx + "delta_z_prob_positive"] = s["prob_direction"]
    out = pd.DataFrame([r])
    out.to_csv(RESULTS / "followup-mature-bayesian-concordance.csv", index=False)
    print(f"  delta_r med={r['posterior_median']:+.4f} "
          f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] "
          f"P(>0)={r['prob_direction']:.4f} P(ROPE)={r['prob_rope']:.3f} "
          f"delta_z med={r['delta_z_median']:+.4f} freq p={p_freq:.4f}")
    return out


def main():
    item1_composition()
    item2_narrowing()
    item3_naive_cd4()
    item4_transport()
    item5_drift_mouse()
    item6_concordance()
    print("\nAll six items done. Human PyMC extension runs separately on beast "
          "(scripts/run_mature_bayesian_human_drift.py).")


if __name__ == "__main__":
    main()
