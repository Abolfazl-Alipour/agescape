"""Donor-spread variance test: distance-variance version of donor-cloud spread.

Reads ONLY the committed donor tables. Donor is the unit. Pooled-group rows
only (Gao All_LV_Nuclei, Yazar All_Cells), one row per donor.
Read-only on the repo; writes only to the output directory.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)
SEED = 7
N_BOOT = 2000
rng = np.random.default_rng(SEED)

PRIMARY = "a1_dist_full_1500"          # Hellinger distance to young reference
MATURE = "a1_dist_mature_ref_mid_loo"  # Hellinger distance to mature (MID) reference


def descriptives(s):
    s = s.dropna()
    q1, med, q3 = s.quantile(0.25), s.quantile(0.50), s.quantile(0.75)
    return {
        "n": int(s.count()),
        "mean": float(s.mean()),
        "sd": float(s.std(ddof=1)),
        "q1": float(q1),
        "median": float(med),
        "q3": float(q3),
        "iqr": float(q3 - q1),
    }


def compare(a, b):
    """Variance comparison of group a vs group b (arrays, no NaN).
    Returns Levene (median-centered) stat/p, variances, ratio a/b,
    bootstrap percentile CI for the ratio."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    stat, p = stats.levene(a, b, center="median")
    va = float(np.var(a, ddof=1))
    vb = float(np.var(b, ddof=1))
    ratio = va / vb
    boots = np.empty(N_BOOT)
    for i in range(N_BOOT):
        ra = rng.choice(a, size=len(a), replace=True)
        rb = rng.choice(b, size=len(b), replace=True)
        boots[i] = np.var(ra, ddof=1) / np.var(rb, ddof=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {
        "levene_stat": float(stat),
        "levene_p": float(p),
        "var_a": va,
        "var_b": vb,
        "var_ratio": float(ratio),
        "ci_low": float(lo),
        "ci_high": float(hi),
    }


def get_bin(df, col, b):
    return df.loc[df["bin"] == b, col].dropna().to_numpy(dtype=float)


results_rows = []
desc_rows = []
notes = []

DATASETS = [
    ("Gao heart", "followup-between-subject-gao-donors.csv", "All_LV_Nuclei"),
    ("Yazar blood", "followup-between-subject-yazar-donors.csv", "All_Cells"),
]

for label, fname, pooled_group in DATASETS:
    df = pd.read_csv(REPO / fname)
    p = df[df["group"] == pooled_group].copy()
    assert p["donor_id"].nunique() == len(p), "pooled rows are not one per donor"
    notes.append(f"{label}: pooled rows {len(p)}, unique donors {p['donor_id'].nunique()}")

    # Per-bin descriptives, primary distance
    for b in ["YOUNG", "MID", "OLD", "LT20"]:
        s = p.loc[p["bin"] == b, PRIMARY]
        if len(s.dropna()) == 0:
            continue
        d = descriptives(s)
        desc_rows.append({"dataset": label, "distance": "young_ref_full_1500",
                          "variable": PRIMARY, "bin": b, **d})

    # Primary comparison: OLD vs YOUNG on young-ref distance
    old = get_bin(p, PRIMARY, "OLD")
    young = get_bin(p, PRIMARY, "YOUNG")
    mid = get_bin(p, PRIMARY, "MID")
    for comp, a, b, an, bn in [
        ("OLD_vs_YOUNG", old, young, "OLD", "YOUNG"),
        ("MID_vs_YOUNG", mid, young, "MID", "YOUNG"),
        ("OLD_vs_MID", old, mid, "OLD", "MID"),
    ]:
        c = compare(a, b)
        results_rows.append({
            "dataset": label, "distance": "young_ref_full_1500", "variable": PRIMARY,
            "comparison": comp, "group_a": an, "group_b": bn,
            "n_a": len(a), "n_b": len(b),
            "sd_a": float(np.std(a, ddof=1)), "sd_b": float(np.std(b, ddof=1)),
            "method": "Levene median-centered (Brown-Forsythe); ratio bootstrap percentile CI, 2000 resamples, seed 7",
            **c,
        })

    # Sensitivity: mature-ref distance. YOUNG has no mature-ref distance by
    # construction (reference is built from MID donors), so compare OLD vs MID.
    old_m = get_bin(p, MATURE, "OLD")
    mid_m = get_bin(p, MATURE, "MID")
    young_m = get_bin(p, MATURE, "YOUNG")
    notes.append(f"{label}: mature-ref distance present: OLD {len(old_m)}, MID {len(mid_m)}, YOUNG {len(young_m)}")
    for b in ["MID", "OLD"]:
        d = descriptives(p.loc[p["bin"] == b, MATURE])
        desc_rows.append({"dataset": label, "distance": "mature_ref_mid_loo",
                          "variable": MATURE, "bin": b, **d})
    c = compare(old_m, mid_m)
    results_rows.append({
        "dataset": label, "distance": "mature_ref_mid_loo", "variable": MATURE,
        "comparison": "OLD_vs_MID", "group_a": "OLD", "group_b": "MID",
        "n_a": len(old_m), "n_b": len(mid_m),
        "sd_a": float(np.std(old_m, ddof=1)), "sd_b": float(np.std(mid_m, ddof=1)),
        "method": "Levene median-centered (Brown-Forsythe); ratio bootstrap percentile CI, 2000 resamples, seed 7",
        **c,
    })

    # Yazar-only pool-residualized version
    if label == "Yazar blood":
        pool_sizes = p.groupby("study_or_pool").size()
        singletons = pool_sizes[pool_sizes == 1].index.tolist()
        notes.append(f"Yazar: {p['study_or_pool'].nunique()} pools, sizes {pool_sizes.min()} to {pool_sizes.max()}, "
                     f"singleton pools: {len(singletons)}")
        pr = p[~p["study_or_pool"].isin(singletons)].copy()
        grand = pr[PRIMARY].mean()
        pool_mean = pr.groupby("study_or_pool")[PRIMARY].transform("mean")
        pr["dist_pool_resid"] = pr[PRIMARY] - pool_mean + grand
        notes.append(f"Yazar pool-residualized: excluded {len(p) - len(pr)} donors in singleton pools, kept {len(pr)}")
        for b in ["YOUNG", "MID", "OLD", "LT20"]:
            s = pr.loc[pr["bin"] == b, "dist_pool_resid"]
            if len(s.dropna()) == 0:
                continue
            d = descriptives(s)
            desc_rows.append({"dataset": label, "distance": "young_ref_pool_residualized",
                              "variable": "a1_dist_full_1500 minus pool mean plus grand mean", "bin": b, **d})
        old_r = get_bin(pr, "dist_pool_resid", "OLD")
        young_r = get_bin(pr, "dist_pool_resid", "YOUNG")
        mid_r = get_bin(pr, "dist_pool_resid", "MID")
        for comp, a, b, an, bn in [
            ("OLD_vs_YOUNG", old_r, young_r, "OLD", "YOUNG"),
            ("OLD_vs_MID", old_r, mid_r, "OLD", "MID"),
        ]:
            c = compare(a, b)
            results_rows.append({
                "dataset": label, "distance": "young_ref_pool_residualized",
                "variable": "a1_dist_full_1500 minus pool mean plus grand mean",
                "comparison": comp, "group_a": an, "group_b": bn,
                "n_a": len(a), "n_b": len(b),
                "sd_a": float(np.std(a, ddof=1)), "sd_b": float(np.std(b, ddof=1)),
                "method": "Levene median-centered (Brown-Forsythe); ratio bootstrap percentile CI, 2000 resamples, seed 7; singleton-pool donors excluded",
                **c,
            })

# Per-cell-type supplement: SD ratio OLD vs YOUNG within each group,
# groups with at least 20 donors with non-missing distance in each bin.
supp_rows = []
for label, fname, pooled_group in DATASETS:
    df = pd.read_csv(REPO / fname)
    for grp, g in df[df["group"] != pooled_group].groupby("group"):
        old = g.loc[g["bin"] == "OLD", PRIMARY].dropna().to_numpy(dtype=float)
        young = g.loc[g["bin"] == "YOUNG", PRIMARY].dropna().to_numpy(dtype=float)
        if len(old) < 20 or len(young) < 20:
            continue
        stat, pval = stats.levene(old, young, center="median")
        supp_rows.append({
            "dataset": label, "cell_type_group": grp,
            "n_old": len(old), "n_young": len(young),
            "sd_old": float(np.std(old, ddof=1)),
            "sd_young": float(np.std(young, ddof=1)),
            "var_ratio_old_over_young": float(np.var(old, ddof=1) / np.var(young, ddof=1)),
            "levene_stat": float(stat), "levene_p": float(pval),
        })

res_df = pd.DataFrame(results_rows)
res_df.to_csv(OUT / "donor-spread-variance-results.csv", index=False)
desc_df = pd.DataFrame(desc_rows)
desc_df.to_csv(OUT / "donor-spread-variance-descriptives.csv", index=False)
supp_df = pd.DataFrame(supp_rows)
supp_df.to_csv(OUT / "donor-spread-variance-by-cell-type.csv", index=False)

print("\n".join(notes))
print("\n=== DESCRIPTIVES ===")
print(desc_df.to_string(index=False))
print("\n=== COMPARISONS ===")
print(res_df[["dataset", "distance", "comparison", "n_a", "n_b", "sd_a", "sd_b",
              "var_ratio", "ci_low", "ci_high", "levene_stat", "levene_p"]].to_string(index=False))
print("\n=== PER-GROUP SUPPLEMENT ===")
print(supp_df.to_string(index=False))
