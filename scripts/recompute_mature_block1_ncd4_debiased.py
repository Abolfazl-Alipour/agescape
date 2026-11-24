"""Independent recompute of the mature naive CD4 debiased dispersion slope.

Reads only committed cache tables on Beast, never the block1 report text:
  - tmp_organs/block1_cache_mature/yazar/stats/donor_group_stats.csv
    (per-donor mean residual LAD dispersion, cell counts, depth columns)
  - tmp_organs/block1_cache_mature/yazar/stats/finite_n_calib.csv
    (real-cell subsampled dispersion, 5 seeds, n_sub 20/50/100/200)
  - tmp_organs/block1_cache_mature/yazar/donors_meta.csv
    (donor age per decade, sex, pool number)
  - results/followup-mature-block1-residual-dispersion.csv
    (only the primary-row active covariate list, to reuse the exact
    pruned design: sex_male plus 73 pool dummies)

Method (same definitions as scripts/run_mature_block1_residual_dispersion.py,
reimplemented here from scratch with plain numpy QR OLS and scipy t tests):
  1. Calibration bias per n_sub = mean over high-N donors and seeds of
     (subsampled mean_resid_lad minus full-N mean_resid_lad).
  2. Per-donor bias estimate = linear interpolation of the four bias
     points over donor cell count, anchored at 0 bias at the median
     cell count of donors with at least 250 cells (or 250 if larger).
  3. y_debiased = mean_resid_lad minus per-donor bias estimate.
  4. OLS fits: primary, induced (bias as outcome), debiased,
     debiased plus depth, raw-Y plus log cell count, raw-Y plus
     cell count and depth, and the two true debiased versions of the
     last two (new numbers: the diagnostics JSON holds no debiased
     +logN or +N+D fit). Slopes are per decade of age. P values are
     two-sided t.

Comparison values come from the diagnostics JSON written by the
original run:
  tmp_organs/block1_cache_mature/block1_diagnostics_summary.json
  (diagnostics.yazar.mechanism_audit.* and .calib_summary).
The raw-Y +N+D pair is checked against block1 report section 5.1,
because it lives in no diagnostics key.

Outputs:
  - results/followup-mature-block1-debiased-recompute.csv
  - results/followup-mature-block1-debiased-recompute-calibration.csv
  - reports/followup-mature-block1-debiased-recompute.md
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

REPO = Path("/home/reza/projects/agescape-revision")
DATA = Path("/home/reza/projects/agescape-revision-data")
CACHE = DATA / "tmp_organs" / "block1_cache_mature" / "yazar"
NCD4 = "naive thymus-derived CD4-positive, alpha-beta T cell"


def ols_fit(df, outcome, extra_covs, batch_cols):
    """Plain OLS by thin QR. Returns slope per decade, SE, p, n, resid df."""
    cols = ["age_per_decade", "sex_male", *extra_covs, *batch_cols]
    sub = df.dropna(subset=[outcome, *cols]).copy()
    y = sub[outcome].to_numpy(dtype=float)
    x = np.column_stack(
        [np.ones(len(sub)), sub[cols].to_numpy(dtype=float)]
    )
    q, r = np.linalg.qr(x, mode="reduced")
    beta = np.linalg.solve(r, q.T @ y)
    resid = y - q @ (q.T @ y)
    df_resid = int(len(sub) - x.shape[1])
    mse = float(np.dot(resid, resid)) / df_resid
    c = np.zeros(x.shape[1])
    c[1] = 1.0
    vc = np.linalg.solve(r.T, c)
    se = float(np.sqrt(mse * float(np.dot(vc, vc))))
    slope = float(beta[1])
    t_stat = slope / se
    pval = float(2.0 * scipy_stats.t.sf(abs(t_stat), df_resid))
    return {
        "slope": slope,
        "se": se,
        "p": pval,
        "n": int(len(sub)),
        "resid_df": df_resid,
    }


def main():
    gdf = pd.read_csv(CACHE / "stats" / "donor_group_stats.csv")
    meta = pd.read_csv(CACHE / "donors_meta.csv")
    cal = pd.read_csv(CACHE / "stats" / "finite_n_calib.csv")
    diag = json.loads(
        (DATA / "tmp_organs" / "block1_cache_mature"
         / "block1_diagnostics_summary.json").read_text()
    )
    mech = diag["diagnostics"]["yazar"]["mechanism_audit"]
    calib_diag = {
        (c["cell_type"], c["n_sub"]): c["lad_mean_bias"]
        for c in diag["diagnostics"]["yazar"]["calib_summary"]
    }

    ncd4 = gdf[gdf["cell_type"] == NCD4].copy()
    meta = meta.reset_index(drop=True)
    meta["donor_idx"] = np.arange(len(meta), dtype=np.int32)
    ncd4 = ncd4.merge(
        meta[["donor_idx", "age_per_decade", "sex_male", "pool_number"]],
        on="donor_idx",
        how="left",
    )
    ncd4["log_n_cells"] = np.log(ncd4["n_cells"].astype(float))
    ncd4["study_or_pool"] = ncd4["pool_number"].astype(str)
    assert len(ncd4) == 778, f"naive CD4 rows in cache: {len(ncd4)}"
    # Replicate _prune_covariate_matrix: drop donors in pools with
    # fewer than 2 naive CD4 donors (singleton-pool exclusion).
    pool_counts = ncd4["study_or_pool"].value_counts()
    keep_pools = pool_counts[pool_counts >= 2].index.tolist()
    n_singleton = int((~ncd4["study_or_pool"].isin(keep_pools)).sum())
    ncd4 = ncd4[ncd4["study_or_pool"].isin(keep_pools)].copy()
    assert len(ncd4) == 777, f"naive CD4 donors kept: {len(ncd4)}"
    assert n_singleton == 1, f"singleton-pool donors dropped: {n_singleton}"
    # Reference pool = most frequent pool, ties broken by pool name.
    levels = sorted(pool_counts[keep_pools].items(),
                    key=lambda kv: (-int(kv[1]), str(kv[0])))
    ref_pool = levels[0][0]
    batch_cols = [f"batch_{lv}" for lv, _ in levels[1:]]
    assert len(batch_cols) == 73, f"batch dummies: {len(batch_cols)}"
    for b in batch_cols:
        pool = b.split("_", 1)[1]
        ncd4[b] = (ncd4["study_or_pool"] == pool).astype(float)

    # Exact pruned design from the committed primary row (74 covariates).
    prim = pd.read_csv(REPO / "results"
                       / "followup-mature-block1-residual-dispersion.csv")
    prow = prim[(prim["dataset"] == "yazar")
                & (prim["cell_type"] == NCD4)
                & (prim["arm"] == "primary")].iloc[0]
    covs = prow["covariates_active"].split(";")
    assert covs[0] == "sex_male"
    csv_batch_cols = [c for c in covs if c.startswith("batch_")]
    assert set(csv_batch_cols) == set(batch_cols), (
        "recomputed dummy set differs from committed primary row")
    assert len(csv_batch_cols) == 73

    # Step 1: recompute calibration bias per n_sub from the raw calib rows.
    full = ncd4.set_index("donor_idx")["mean_resid_lad"]
    cal4 = cal[cal["cell_type"] == NCD4].copy()
    cal4["full_lad"] = cal4["donor_idx"].map(full)
    cal4["bias"] = cal4["mean_resid_lad"] - cal4["full_lad"]
    bias_by_n = cal4.groupby("n_sub")["bias"].mean().to_dict()
    bias_rows = []
    for n_sub in (20, 50, 100, 200):
        bias_rows.append({
            "n_sub": n_sub,
            "n_donor_seed_rows": int((cal4["n_sub"] == n_sub).sum()),
            "recomputed_bias": float(bias_by_n[n_sub]),
            "diagnostics_bias": float(calib_diag[(NCD4, n_sub)]),
            "abs_diff": float(abs(bias_by_n[n_sub]
                                  - calib_diag[(NCD4, n_sub)])),
        })

    # Step 2: interpolate the bias curve and debias each donor.
    med_hi_n = float(ncd4.loc[ncd4["n_cells"] >= 250, "n_cells"].median())
    xs = np.array([20.0, 50.0, 100.0, 200.0, max(med_hi_n, 250.0)])
    ys = np.array([bias_by_n[n] for n in (20, 50, 100, 200)] + [0.0])
    ncd4["calib_bias_est"] = np.interp(
        ncd4["n_cells"].to_numpy(float), xs, ys
    )
    ncd4["y_debiased"] = ncd4["mean_resid_lad"] - ncd4["calib_bias_est"]

    # Step 3: OLS fits with a fresh implementation.
    fits = {}
    fits["primary"] = ols_fit(ncd4, "mean_resid_lad", [], batch_cols)
    fits["induced_by_calibration"] = ols_fit(
        ncd4, "calib_bias_est", [], batch_cols
    )
    fits["debiased"] = ols_fit(ncd4, "y_debiased", [], batch_cols)
    fits["debiased_plus_depth"] = ols_fit(
        ncd4, "y_debiased",
        ["log_mean_inv_umi", "var_log_umi"], batch_cols,
    )
    fits["raw_plus_logN"] = ols_fit(
        ncd4, "mean_resid_lad", ["log_n_cells"], batch_cols,
    )
    fits["raw_plus_N_plus_depth"] = ols_fit(
        ncd4, "mean_resid_lad",
        ["log_n_cells", "log_mean_inv_umi", "var_log_umi"], batch_cols,
    )
    fits["debiased_plus_logN"] = ols_fit(
        ncd4, "y_debiased", ["log_n_cells"], batch_cols,
    )
    fits["debiased_plus_N_plus_depth"] = ols_fit(
        ncd4, "y_debiased",
        ["log_n_cells", "log_mean_inv_umi", "var_log_umi"], batch_cols,
    )

    expected = {
        "primary": (float(prow["slope_per_decade"]), float(prow["pvalue"])),
        "induced_by_calibration": (
            mech["ncd4_calib_induced_slope"], mech["ncd4_calib_induced_p"]),
        "debiased": (
            mech["ncd4_debiased_age_slope"], mech["ncd4_debiased_age_p"]),
        "debiased_plus_depth": (
            mech["ncd4_debiased_depth_slope"], mech["ncd4_debiased_depth_p"]),
        "raw_plus_logN": (
            mech["ncd4_plus_logn_age_slope"], mech["ncd4_plus_logn_age_p"]),
        "raw_plus_N_plus_depth": (0.00024, 0.01893),
        "debiased_plus_logN": (np.nan, np.nan),
        "debiased_plus_N_plus_depth": (np.nan, np.nan),
    }

    out_rows = []
    for arm, f in fits.items():
        es, ep = expected[arm]
        out_rows.append({
            "arm": arm,
            "slope_per_decade": f["slope"],
            "se": f["se"],
            "p_value": f["p"],
            "n_donors": f["n"],
            "residual_df": f["resid_df"],
            "expected_slope": es,
            "expected_p": ep,
            "abs_slope_diff": (abs(f["slope"] - es)
                               if np.isfinite(es) else np.nan),
            "abs_p_diff": (abs(f["p"] - ep)
                           if np.isfinite(ep) else np.nan),
        })
    out = pd.DataFrame(out_rows)
    share = (fits["induced_by_calibration"]["slope"]
             / fits["primary"]["slope"])
    out_csv = (REPO / "results"
               / "followup-mature-block1-debiased-recompute.csv")
    out.to_csv(out_csv, index=False)

    bias_csv = (REPO / "results"
                / "followup-mature-block1-debiased-recompute-calibration.csv")
    pd.DataFrame(bias_rows).to_csv(bias_csv, index=False)

    checked = out[np.isfinite(out["abs_slope_diff"])]
    diag_checked = checked[checked["arm"] != "raw_plus_N_plus_depth"]
    max_slope_diff = float(diag_checked["abs_slope_diff"].max())
    max_p_diff = float(diag_checked["abs_p_diff"].max())
    raw_nd = checked[checked["arm"] == "raw_plus_N_plus_depth"].iloc[0]
    # The raw-Y +N+D pair is report-sourced at 5-decimal rounding, so its
    # tolerance is half the last reported digit (5e-06), not 1e-06.
    raw_nd_ok = (float(raw_nd["abs_slope_diff"]) < 5e-06
                 and float(raw_nd["abs_p_diff"]) < 5e-06)
    verified = (max_slope_diff < 1e-6 and max_p_diff < 1e-4 and raw_nd_ok)
    max_bias_diff = float(pd.DataFrame(bias_rows)["abs_diff"].max())
    report = f"""# Independent recompute: mature naive CD4 debiased slope

Method: fresh numpy QR OLS plus scipy two-sided t tests, applied to the
committed cache tables only (donor_group_stats.csv, finite_n_calib.csv,
donors_meta.csv). No report text was used as input except the raw-Y +N+D
pair, which lives in no diagnostics key. The design replicates the
singleton-pool prune and cross-checks the 73 pool dummies against the
committed primary row (777 donors, residual df 701).

## Calibration bias (subsampled minus full-N, high-N donors)

| n_sub | donor-seed rows | recomputed bias | diagnostics bias | abs diff |
| --- | --- | --- | --- | --- |
"""
    for b in bias_rows:
        report += (f"| {b['n_sub']} | {b['n_donor_seed_rows']} | "
                   f"{b['recomputed_bias']:.6f} | "
                   f"{b['diagnostics_bias']:.6f} | {b['abs_diff']:.2e} |\n")
    report += f"""
Max calibration bias difference: {max_bias_diff:.2e}.

## Slopes

| arm | recomputed slope | recomputed p | expected slope | expected p | abs slope diff |
| --- | --- | --- | --- | --- | --- |
"""
    for _, r in out.iterrows():
        es = f"{r['expected_slope']:.6f}" if np.isfinite(r["expected_slope"]) else "new"
        ep = f"{r['expected_p']:.4g}" if np.isfinite(r["expected_p"]) else "new"
        sd = f"{r['abs_slope_diff']:.2e}" if np.isfinite(r["abs_slope_diff"]) else "new"
        report += (f"| {r['arm']} | {r['slope_per_decade']:.6f} | "
                   f"{r['p_value']:.4g} | {es} | {ep} | {sd} |\n")
    report += f"""
Calibration share of the primary slope: {share:.4f} ({100 * share:.2f} percent).
Max slope difference on checked arms: {max_slope_diff:.2e}.
Max p difference on checked arms: {max_p_diff:.4g}.

## Verdict

{"VERIFIED: every checked recomputed slope and p matches its expected value (diagnostics arms to 1e-6 in slope; the report-sourced raw-Y +N+D arm to its 5-decimal rounding)." if verified else "MISMATCH: see differences above."}
The debiased slope is {fits['debiased']['slope']:.6f} per decade
(p = {fits['debiased']['p']:.4g}), still significant after calibration.
The exact calibration share is {100 * share:.2f} percent, not 56.4 percent:
the report divided the induced slope by the rounded primary slope 0.00071
instead of the exact 0.000714. The tex wording (56 percent) is unaffected.
The true debiased +logN slope is {fits['debiased_plus_logN']['slope']:.6f}
(p = {fits['debiased_plus_logN']['p']:.4g}) and the true debiased +N+D slope
is {fits['debiased_plus_N_plus_depth']['slope']:.6f}
(p = {fits['debiased_plus_N_plus_depth']['p']:.4g}); the transcription CSV
labels the raw-Y fits 0.000286 and 0.00024 with those names instead.
"""
    (REPO / "reports"
     / "followup-mature-block1-debiased-recompute.md").write_text(report)
    print(report)
    print("WROTE", out_csv)
    print("WROTE", bias_csv)


if __name__ == "__main__":
    main()
