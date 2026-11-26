#!/usr/bin/env python3
"""Study-aware evidence and closest analogs for the fly/worm fraction rows.

Study-aware check: IMPOSSIBLE for both species. The mature fly data are a
single study (AFCA) by construction, and the mature worm data are a single
dataset (Gao GSE229022) with no study column. This script recomputes that
evidence from the loaded data and saves it, so the claim is traceable.

Analogs (clearly labeled as analogs, NOT the study-aware check):
(a) worm platform-aware: refit the ARSC slope with platform_nextseq as a
    covariate (9 pools; platform is balanced 1 vs 2 within each age), plus
    per-platform slopes (platform 0: 3 pools, 1 residual df, fragile;
    platform 1: 6 pools).
(b) fly sex-stratified: refit each of the 25 fly rows within male pools only
    (18 pools) and within female pools only (18 pools).

All refits use the exact scan model (CLR with 0.5 pseudocount, OLS on age;
sex dropped when constant). Writes results/followup-mature-fly-worm-analogs.csv.
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_afca_aging import afca_obs
from scripts.run_gao_worm_aging import normal_obs as worm_normal_obs
from scripts.run_mature_all_species_fractions import (
    MIN_PRESENT, PSEUDOCOUNT, fit_cell_type, load_strata,
)
from scripts.run_mature_composition_complete import clr_table

SCAN = REPO / "results" / "mature-all-species-composition.csv"
OUTPUT = REPO / "results" / "followup-mature-fly-worm-analogs.csv"


def ols_slope(y, design):
    """OLS of y on design (DataFrame, first column is age). Returns effect, se, pvalue, df_resid."""
    X = sm.add_constant(design.astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or len(y) <= X.shape[1]:
        return None
    fit = sm.OLS(y, X).fit()
    if not np.isfinite(fit.pvalues["age"]):
        return None
    return {"effect": float(fit.params["age"]), "se": float(fit.bse["age"]),
            "pvalue": float(fit.pvalues["age"]), "df_resid": int(fit.df_resid)}


def study_evidence():
    """Recompute the single-study evidence from the loaded data."""
    rows = []
    for organ in ("head", "body"):
        obs = afca_obs(ad.read_h5ad(REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad", backed="r"))
        raw_counts = obs["dataset"].value_counts().to_dict()
        mature = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated") & (obs["age"] >= 30)]
        rows.append({
            "check": "study_evidence", "species": "Fruit fly", "organ": organ, "cohort": "AFCA",
            "cell_type": "", "subset": "all", "n_pools": int(mature["sample"].nunique()),
            "effect": np.nan, "se": np.nan, "pvalue": np.nan, "sign_matches_scan": np.nan,
            "note": (f"raw file dataset counts {raw_counts}; FCA cells are all age 5 and "
                     f"excluded by the mature filter; after dataset==AFCA and age>=30, "
                     f"{int(mature['sample'].nunique())} pools, every cell dataset AFCA; "
                     f"no study/batch column exists in the AFCA obs"),
        })
    obs = worm_normal_obs(ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r"))
    mature = obs[obs["age"] >= 6]
    study_cols = [c for c in obs.columns if "study" in c.lower() or "dataset" in c.lower() or "batch" in c.lower()]
    rows.append({
        "check": "study_evidence", "species": "C. elegans", "organ": "Whole animal",
        "cohort": "Gao_GSE229022", "cell_type": "", "subset": "all",
        "n_pools": int(mature["donor"].nunique()),
        "effect": np.nan, "se": np.nan, "pvalue": np.nan, "sign_matches_scan": np.nan,
        "note": (f"single integrated file from one GEO series (GSE229022); no study/dataset/batch "
                 f"column in obs (checked: {study_cols}); mature pools are wild-type N2 only "
                 f"(genotypes in mature set: {sorted(mature['genotype'].unique())}); "
                 f"orig.ident values are sample IDs inside the one study"),
    })
    return rows


def main():
    scan = pd.read_csv(SCAN)
    passing = scan[(scan["significant_pooled_05"]) & (scan["species"].isin(["Fruit fly", "C. elegans"]))].copy()
    assert len(passing) == 26, f"expected 26 passing fly/worm rows, got {len(passing)}"
    strata = {(s, o, c): (counts, meta) for s, o, c, counts, meta in load_strata()
              if s in ("Fruit fly", "C. elegans")}
    records = study_evidence()

    # (b) fly sex-stratified refits
    for row in passing[passing["species"] == "Fruit fly"].itertuples():
        key = (row.species, row.organ, row.cohort)
        counts_all, metadata = strata[key]
        present = (counts_all > 0).sum(axis=0)
        counts = counts_all.loc[:, present >= MIN_PRESENT].copy()
        clr = clr_table(counts, pseudocount=PSEUDOCOUNT)
        for label, sex_value in (("male_only", 1.0), ("female_only", 0.0)):
            pools = metadata.index[metadata["sex"] == sex_value]
            refit, reason = fit_cell_type(clr.loc[pools], counts.loc[pools], metadata,
                                          row.cell_type, *key)
            if refit is None:
                records.append({
                    "check": "sex_analog", "species": row.species, "organ": row.organ,
                    "cohort": row.cohort, "cell_type": row.cell_type, "subset": label,
                    "n_pools": int(len(pools)), "effect": np.nan, "se": np.nan,
                    "pvalue": np.nan, "sign_matches_scan": np.nan,
                    "note": f"refit impossible: {reason}",
                })
                continue
            records.append({
                "check": "sex_analog", "species": row.species, "organ": row.organ,
                "cohort": row.cohort, "cell_type": row.cell_type, "subset": label,
                "n_pools": int(len(pools)), "effect": refit["effect"], "se": refit["se"],
                "pvalue": refit["pvalue"],
                "sign_matches_scan": bool(np.sign(refit["effect"]) == np.sign(row.effect)),
                "note": f"scan effect {row.effect:+.4f} per day; {label} refit on {len(pools)} pools",
            })

    # (a) worm platform-aware refits for ARSC
    wrow = passing[passing["species"] == "C. elegans"].iloc[0]
    key = (wrow.species, wrow.organ, wrow.cohort)
    counts_all, metadata = strata[key]
    present = (counts_all > 0).sum(axis=0)
    counts = counts_all.loc[:, present >= MIN_PRESENT].copy()
    clr = clr_table(counts, pseudocount=PSEUDOCOUNT)
    obs = worm_normal_obs(ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r"))
    plat = obs[obs["age"] >= 6].groupby("donor")["platform_nextseq"].first().loc[counts.index].astype(float)
    y = clr[wrow.cell_type]
    design = pd.DataFrame({"age": metadata.loc[counts.index, "age"], "platform_nextseq": plat})
    adj = ols_slope(y, design)
    assert adj is not None, "platform-adjusted worm fit failed"
    records.append({
        "check": "platform_analog", "species": wrow.species, "organ": wrow.organ,
        "cohort": wrow.cohort, "cell_type": wrow.cell_type, "subset": "platform_adjusted",
        "n_pools": int(len(counts)), "effect": adj["effect"], "se": adj["se"],
        "pvalue": adj["pvalue"],
        "sign_matches_scan": bool(np.sign(adj["effect"]) == np.sign(wrow.effect)),
        "note": (f"scan effect {wrow.effect:+.4f} per day; platform balanced 1 vs 2 pools "
                 f"within each age (6, 12, 14 days); residual df {adj['df_resid']}"),
    })
    for label, value in (("platform_0", 0.0), ("platform_1", 1.0)):
        pools = plat.index[plat == value]
        sub = ols_slope(y.loc[pools], pd.DataFrame({"age": metadata.loc[pools, "age"]}))
        if sub is None:
            records.append({
                "check": "platform_analog", "species": wrow.species, "organ": wrow.organ,
                "cohort": wrow.cohort, "cell_type": wrow.cell_type, "subset": label,
                "n_pools": int(len(pools)), "effect": np.nan, "se": np.nan,
                "pvalue": np.nan, "sign_matches_scan": np.nan,
                "note": f"refit impossible on {len(pools)} pools",
            })
            continue
        records.append({
            "check": "platform_analog", "species": wrow.species, "organ": wrow.organ,
            "cohort": wrow.cohort, "cell_type": wrow.cell_type, "subset": label,
            "n_pools": int(len(pools)), "effect": sub["effect"], "se": sub["se"],
            "pvalue": sub["pvalue"],
            "sign_matches_scan": bool(np.sign(sub["effect"]) == np.sign(wrow.effect)),
            "note": (f"scan effect {wrow.effect:+.4f} per day; {len(pools)} pools, "
                     f"residual df {sub['df_resid']}" +
                     ("; FRAGILE: 3 pools, 1 residual df" if label == "platform_0" else "")),
        })

    out = pd.DataFrame(records)
    out.to_csv(OUTPUT, index=False)
    print(f"wrote {len(out)} rows to {OUTPUT}")
    sex_rows = out[out["check"] == "sex_analog"]
    print(f"sex analog: sign matches scan in {int(sex_rows['sign_matches_scan'].sum())} of {len(sex_rows)} refits")


if __name__ == "__main__":
    main()
