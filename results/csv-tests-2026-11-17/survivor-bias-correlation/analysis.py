"""Survivor-bias cross-type correlation. Committed CSVs only, read-only on repo."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = str(Path(__file__).resolve().parents[2])
OUT = str(Path(__file__).resolve().parent)
SEED = 20261117
rng = np.random.default_rng(SEED)

comp = pd.read_csv(f"{R}/reviewer-phase1-composition-global-bh.csv")
vec = pd.read_csv(f"{R}/h3-cross-species-vector-tests.csv")

# Primary drift endpoint: Hellinger distance slope, depth-adjusted model.
vec_p = vec[(vec["metric"] == "profile_hellinger") & (vec["model"] == "depth_adjusted")].copy()
vec_p = vec_p.rename(columns={"tissue": "organ", "slope_per_year": "drift_slope_per_year"})
vec_u = vec[(vec["metric"] == "profile_hellinger") & (vec["model"] == "unadjusted")].copy()
vec_u = vec_u.rename(columns={"tissue": "organ", "slope_per_year": "drift_slope_unadj_per_year"})

key = ["species", "organ", "cohort", "cell_type"]
comp = comp.rename(columns={"effect": "fraction_clr_slope_per_year"})

j = comp.merge(
    vec_p[key + ["drift_slope_per_year", "n_donors", "n_cells"]],
    on=key, how="inner", suffixes=("_comp", "_drift"),
)
j = j.merge(
    vec_u[key + ["drift_slope_unadj_per_year"]], on=key, how="left"
)
assert j["fraction_clr_slope_per_year"].notna().all()
assert j["drift_slope_per_year"].notna().all()
j["share_loss"] = -j["fraction_clr_slope_per_year"]
j["abs_drift_slope"] = j["drift_slope_per_year"].abs()

# Joined scatter-data CSV
scatter_cols = ["species", "organ", "cohort", "cell_type",
                "fraction_clr_slope_per_year", "drift_slope_per_year",
                "drift_slope_unadj_per_year", "share_loss", "abs_drift_slope",
                "n_donors_comp", "n_donors_drift"]
scatter = j[["species", "organ", "cohort", "cell_type",
             "fraction_clr_slope_per_year", "drift_slope_per_year",
             "drift_slope_unadj_per_year", "share_loss", "abs_drift_slope"]].copy()
scatter = scatter.rename(columns={
    "fraction_clr_slope_per_year": "fraction_slope",
    "drift_slope_per_year": "drift_slope",
    "drift_slope_unadj_per_year": "drift_slope_unadjusted",
})
scatter = scatter.sort_values(["species", "organ", "cohort", "cell_type"])
scatter.to_csv(f"{OUT}/survivor-bias-joined-scatter-data.csv", index=False)

def perm_p(x, y, func, n_perm=10000):
    obs = func(x, y)
    cnt = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        if abs(func(x, yp)) >= abs(obs):
            cnt += 1
    return obs, (cnt + 1) / (n_perm + 1)

def spearman_r(x, y):
    return stats.spearmanr(x, y).statistic

def pearson_r(x, y):
    return stats.pearsonr(x, y).statistic

rows = []
for scope, sub in [("Mouse", j[j["species"] == "Mouse"]),
                   ("Human", j[j["species"] == "Human"])]:
    n = len(sub)
    tests = [
        ("fraction_slope_vs_drift_slope", sub["fraction_clr_slope_per_year"], sub["drift_slope_per_year"]),
        ("fraction_slope_vs_abs_drift_slope", sub["fraction_clr_slope_per_year"], sub["abs_drift_slope"]),
        ("share_loss_vs_drift_slope", sub["share_loss"], sub["drift_slope_per_year"]),
        ("share_loss_vs_abs_drift_slope", sub["share_loss"], sub["abs_drift_slope"]),
        ("fraction_slope_vs_drift_slope_unadjusted", sub["fraction_clr_slope_per_year"], sub["drift_slope_unadj_per_year"]),
    ]
    for name, x, y in tests:
        x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
        pr = stats.pearsonr(x, y)
        sr = stats.spearmanr(x, y)
        _, pp_p = perm_p(x, y, pearson_r)
        _, pp_s = perm_p(x, y, spearman_r)
        rows.append({"scope": scope, "test": name, "method": "Pearson",
                     "correlation": pr.statistic, "p_value": pr.pvalue,
                     "p_value_permutation": pp_p, "n_units": n})
        rows.append({"scope": scope, "test": name, "method": "Spearman",
                     "correlation": sr.statistic, "p_value": sr.pvalue,
                     "p_value_permutation": pp_s, "n_units": n})

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/survivor-bias-correlation-results.csv", index=False)

# Descriptive per-organ means (no tests)
desc = (j.groupby(["species", "organ"])
          .agg(n_units=("cell_type", "size"),
               mean_fraction_slope=("fraction_clr_slope_per_year", "mean"),
               median_fraction_slope=("fraction_clr_slope_per_year", "median"),
               mean_drift_slope=("drift_slope_per_year", "mean"),
               median_drift_slope=("drift_slope_per_year", "median"),
               mean_abs_drift_slope=("abs_drift_slope", "mean"),
               n_losing_share=("fraction_clr_slope_per_year", lambda s: int((s < 0).sum())))
          .reset_index())
desc.to_csv(f"{OUT}/survivor-bias-per-organ-descriptive.csv", index=False)

# Descriptive split: losing vs gaining units
for scope, sub in [("Mouse", j[j["species"] == "Mouse"]),
                   ("Human", j[j["species"] == "Human"])]:
    los = sub[sub["fraction_clr_slope_per_year"] < 0]
    gai = sub[sub["fraction_clr_slope_per_year"] >= 0]
    print(f"{scope}: N={len(sub)} losing={len(los)} gaining={len(gai)}")
    print(f"  losing mean drift={los['drift_slope_per_year'].mean():.6f} mean abs drift={los['abs_drift_slope'].mean():.6f}")
    print(f"  gaining mean drift={gai['drift_slope_per_year'].mean():.6f} mean abs drift={gai['abs_drift_slope'].mean():.6f}")

print()
print(res.to_string(index=False))
print()
print(desc.to_string(index=False))
print()
print("joined N:", len(j), "| mouse:", (j['species']=='Mouse').sum(), "| human:", (j['species']=='Human').sum())
