#!/usr/bin/env python3
"""Item 1: Fresh CLR and propeller/limma cell-type composition analysis for Fly, Worm, and Zebrafish.

Matches the exact primary pipeline used for Mouse and Human in scripts/run_organ_composition.py:
  - Donor/pool eligibility: present (count > 0) in >= 6 pools/fish (no total_cells >= 50 exclusion in primary).
  - Primary transformation: centered log-ratio (CLR) with pseudocount = 0.5.
  - Sensitivity transformations: CLR with pseudocount = 0.1 and 1.0.
  - Primary model: OLS of CLR(cell_type) ~ 1 + age + sex_male (dropping sex_male if constant/rank-deficient).
  - Additional sensitivities (matching extended scripts & reviewer phase 1):
      * observed_records_only (samples with count > 0, pseudocount 0.5)
      * total_cells_adjusted (+ log_total_cells)
      * platform_adjusted / platform_and_total_adjusted (+ platform_nextseq for Worm)
      * total_cells_ge_50 (excluding cell types with < 50 total cells)
      * somatic_only (excluding reproductive compartments in Worm)
      * zebrafish_10fish_broad (pooling 5-22mo and 27mo cohorts on 3 broad compartments for reconciliation)
  - Propeller + limma-style empirical Bayes check:
      * Logit-transformed pseudocount proportions log(p / (1 - p)) where p = (count + 0.5) / rowSum(counts + 0.5).
      * OLS fit per cell type followed by Smyth (2004) empirical-Bayes variance moderation (fitFDist).
      * Cross-checked against the committed R limma/speckle output in results/reviewer-phase1-composition-propeller.csv.
  - Multiple testing:
      * Within-stratum BH correction (q_within_bh) per (species, organ, cohort, analysis, sensitivity, pseudocount).
      * One combined global BH family (q_global_bh) uniting the 193 new primary tests with the 431 existing
        Mouse, Human, and NMR primary tests in results/reviewer-phase1-composition-global-bh.csv (total m = 624).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import special, stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_afca_aging import afca_obs
from scripts.run_gao_worm_aging import REPRODUCTIVE as WORM_REPRODUCTIVE, normal_obs as worm_normal_obs
from scripts.run_zebrafish_testis_aging import standard_obs as zfish_standard_obs

RESULTS = REPO / "results"
AGESCAPE_DATA = Path("/home/reza/projects/agescape/data/interim")

# Reference maximum-lifespan 10%-increment scales used across the Agescape codebase
LIFESPAN_10PCT_SCALES = {
    "Fruit fly": ("day", 4.0),       # max lifespan 40 days -> 4.0 days
    "C. elegans": ("day", 1.4),      # max lifespan 14 days -> 1.4 days
    "Zebrafish": ("month", 2.2),     # cohort span/reference scale 22 months -> 2.2 months
}


def clr_table(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    adjusted = counts.astype(float) + float(pseudocount)
    proportions = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_proportions = np.log(proportions)
    return log_proportions.sub(log_proportions.mean(axis=1), axis=0)


def filter_present_ge(counts: pd.DataFrame, min_present: int = 6) -> pd.DataFrame:
    return counts.loc[:, (counts > 0).sum(axis=0) >= min_present].copy()


def _trigamma_inverse(x: float) -> float:
    if x <= 0:
        return float(np.inf)
    if x > 1e7:
        return float(1.0 / np.sqrt(x))
    if x < 1e-6:
        return float(1.0 / x)
    y = 0.5 + 1.0 / x
    for _ in range(60):
        tri = float(special.polygamma(1, y))
        dif = tri * (1.0 - tri / x) / float(special.polygamma(2, y))
        y_new = max(y + dif, 1e-12)
        if abs(dif / y_new) < 1e-11:
            return float(y_new)
        y = y_new
    return float(y)


def fit_f_dist_smyth(s2: np.ndarray, df1: float) -> tuple[float, float]:
    """Smyth (2004) method-of-moments empirical Bayes F-prior estimation."""
    s2_clean = np.maximum(np.asarray(s2, dtype=float), 1e-15)
    e = np.log(s2_clean) - float(special.digamma(df1 / 2.0)) + np.log(df1 / 2.0)
    emean = float(np.mean(e))
    if len(e) < 2:
        return float(np.exp(emean)), 0.0
    evar = float(np.var(e, ddof=1))
    diff = evar - float(special.polygamma(1, df1 / 2.0))
    if diff <= 0:
        return float(np.exp(emean)), float(np.inf)
    df0 = 2.0 * _trigamma_inverse(diff)
    s02 = float(np.exp(emean + float(special.digamma(df0 / 2.0)) - np.log(df0 / 2.0)))
    return s02, float(df0)


def _design_matrix(
    metadata: pd.DataFrame,
    donors: pd.Index,
    counts: pd.DataFrame,
    include_log_total: bool = False,
    include_platform: bool = False,
    include_spermatogenesis: bool = False,
) -> tuple[pd.DataFrame, str]:
    sub = metadata.loc[donors].copy()
    cols = ["age"]
    tags = ["age"]
    if "sex_male" in sub.columns and sub["sex_male"].nunique() > 1:
        cols.append("sex_male")
        tags.append("sex")
    else:
        tags.append("rank_deficient_sex_dropped")
    if include_log_total:
        log_tot = np.log(counts.loc[donors].sum(axis=1).astype(float))
        sub["log_total_cells"] = log_tot - log_tot.mean()
        cols.append("log_total_cells")
        tags.append("log_total_cells")
    if include_platform and "platform_nextseq" in sub.columns and sub["platform_nextseq"].nunique() > 1:
        cols.append("platform_nextseq")
        tags.append("platform_nextseq")
    if (
        include_spermatogenesis
        and "incomplete_spermatogenesis" in sub.columns
        and sub["incomplete_spermatogenesis"].nunique() > 1
    ):
        cols.append("incomplete_spermatogenesis")
        tags.append("incomplete_spermatogenesis")
    X = sm.add_constant(sub[cols].astype(float), has_constant="add")
    return X, "+".join(tags)


def run_clr_regressions(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    species: str,
    organ: str,
    cohort: str,
    analysis: str,
    sensitivity: str,
    pseudocount: float = 0.5,
    observed_only: bool = False,
    include_log_total: bool = False,
    include_platform: bool = False,
    include_spermatogenesis: bool = False,
    min_donors: int = 6,
) -> list[dict]:
    clr = clr_table(counts, pseudocount=pseudocount)
    age_unit, scale_10pct = LIFESPAN_10PCT_SCALES[species]
    rows: list[dict] = []
    for cell_type in counts.columns:
        donors = counts.index[counts[cell_type] > 0] if observed_only else counts.index
        if len(donors) < min_donors:
            continue
        X, model_tag = _design_matrix(
            metadata,
            donors,
            counts,
            include_log_total=include_log_total,
            include_platform=include_platform,
            include_spermatogenesis=include_spermatogenesis,
        )
        if len(donors) <= X.shape[1] + 1 or np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            continue
        fit = sm.OLS(clr.loc[donors, cell_type], X).fit()
        ci = fit.conf_int().loc["age"]
        est = float(fit.params["age"])
        se = float(fit.bse["age"])
        rows.append(
            {
                "species": species,
                "organ": organ,
                "cohort": cohort,
                "analysis": analysis,
                "cell_type": str(cell_type),
                "method": "CLR_OLS",
                "sensitivity": sensitivity,
                "pseudocount": float(pseudocount),
                "model": model_tag,
                "effect_type": f"CLR_slope_per_{age_unit}",
                "age_unit": age_unit,
                "effect": est,
                "se": se,
                "ci95_low": float(ci.iloc[0]),
                "ci95_high": float(ci.iloc[1]),
                "effect_per_10pct_lifespan": est * scale_10pct,
                "se_per_10pct_lifespan": se * scale_10pct,
                "tvalue": float(fit.tvalues["age"]),
                "df_residual": float(fit.df_resid),
                "pvalue": float(fit.pvalues["age"]),
                "n_donors": int(len(donors)),
                "n_present": int((counts.loc[donors, cell_type] > 0).sum()),
                "total_cells": int(counts.loc[donors, cell_type].sum()),
            }
        )
    if rows:
        qvals = multipletests([r["pvalue"] for r in rows], method="fdr_bh")[1]
        for r, q in zip(rows, qvals):
            r["q_within_bh"] = float(q)
    return rows


def run_propeller_limma(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    species: str,
    organ: str,
    cohort: str,
    analysis: str,
    sensitivity: str,
    pseudocount: float = 0.5,
    include_log_total: bool = False,
    include_platform: bool = False,
) -> list[dict]:
    """Logit-proportion + Smyth (2004) empirical Bayes t-test mirroring scripts/run_reviewer_phase1_propeller.R."""
    props = (counts.astype(float) + pseudocount).div((counts.astype(float) + pseudocount).sum(axis=1), axis=0)
    trans = np.log(props / (1.0 - props))
    donors = counts.index
    X, model_tag = _design_matrix(
        metadata,
        donors,
        counts,
        include_log_total=include_log_total,
        include_platform=include_platform,
        include_spermatogenesis=False,
    )
    if len(donors) <= X.shape[1] or np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return []
    df_res = float(X.shape[0] - X.shape[1])
    age_unit, scale_10pct = LIFESPAN_10PCT_SCALES[species]
    coefs = []
    s2 = []
    stdev_unscaled = []
    cell_types = list(counts.columns)
    for ct in cell_types:
        fit = sm.OLS(trans[ct], X).fit()
        coefs.append(float(fit.params["age"]))
        mse = float(fit.mse_resid)
        s2.append(mse)
        stdev_unscaled.append(float(fit.bse["age"] / np.sqrt(mse)) if mse > 0 else np.nan)
    coefs_arr = np.asarray(coefs, dtype=float)
    s2_arr = np.asarray(s2, dtype=float)
    su_arr = np.asarray(stdev_unscaled, dtype=float)
    s02, df0 = fit_f_dist_smyth(s2_arr, df_res)
    if np.isinf(df0):
        s2_post = np.full_like(s2_arr, s02)
        df_total = np.full_like(s2_arr, np.inf)
    else:
        s2_post = (df_res * s2_arr + df0 * s02) / (df_res + df0)
        df_total = np.full_like(s2_arr, df_res + df0)
    se_post = np.sqrt(s2_post) * su_arr
    t_post = coefs_arr / se_post
    pvals = np.where(
        np.isinf(df_total),
        2.0 * stats.norm.sf(np.abs(t_post)),
        2.0 * stats.t.sf(np.abs(t_post), df=df_total),
    )
    crit = np.where(
        np.isinf(df_total),
        stats.norm.ppf(0.975),
        stats.t.ppf(0.975, df=df_total),
    )
    qvals = multipletests(pvals, method="fdr_bh")[1] if len(pvals) else np.array([])
    rows: list[dict] = []
    for i, ct in enumerate(cell_types):
        est = float(coefs_arr[i])
        se = float(se_post[i])
        rows.append(
            {
                "species": species,
                "organ": organ,
                "cohort": cohort,
                "analysis": analysis,
                "cell_type": str(ct),
                "method": "propeller_logit_eBayes",
                "sensitivity": sensitivity,
                "pseudocount": float(pseudocount),
                "model": model_tag,
                "effect_type": f"logit_slope_per_{age_unit}",
                "age_unit": age_unit,
                "effect": est,
                "se": se,
                "ci95_low": est - float(crit[i]) * se,
                "ci95_high": est + float(crit[i]) * se,
                "effect_per_10pct_lifespan": est * scale_10pct,
                "se_per_10pct_lifespan": se * scale_10pct,
                "tvalue": float(t_post[i]),
                "df_residual": float(df_total[i]),
                "pvalue": float(pvals[i]),
                "n_donors": int(len(donors)),
                "n_present": int((counts[ct] > 0).sum()),
                "total_cells": int(counts[ct].sum()),
                "q_within_bh": float(qvals[i]),
            }
        )
    return rows


def load_stratum_counts() -> tuple[list[dict], dict]:
    """Load donor/pool counts and metadata directly from raw/interim H5AD .obs."""
    strata = []
    summary_meta = {}

    # 1. AFCA head and body (72 pools total: 36 head, 36 body)
    for organ in ("head", "body"):
        adata = ad.read_h5ad(REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad", backed="r")
        obs = afca_obs(adata)
        adata.file.close()
        obs = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated")].copy()
        full_counts = obs.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0).sort_index()
        meta = obs.groupby("sample").agg(age=("age", "first"), sex_male=("sex_male", "first")).loc[full_counts.index]
        counts_p6 = filter_present_ge(full_counts, 6)
        strata.append(
            {
                "species": "Fruit fly",
                "organ": organ,
                "cohort": "AFCA",
                "analysis": "all_annotated_cells",
                "counts_p6": counts_p6,
                "counts_ge50": counts_p6.loc[:, counts_p6.sum(axis=0) >= 50].copy(),
                "metadata": meta,
                "has_platform": False,
                "is_primary": True,
            }
        )
        summary_meta[f"afca_{organ}"] = {
            "n_pools": int(counts_p6.shape[0]),
            "n_total_cells": int(full_counts.to_numpy().sum()),
            "n_cell_types_all": int(full_counts.shape[1]),
            "n_cell_types_present6": int(counts_p6.shape[1]),
            "n_cell_types_present6_sum50": int((counts_p6.sum(axis=0) >= 50).sum()),
        }

    # 2. Gao worm (12 normal N2 pools)
    wad = ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r")
    wobs = worm_normal_obs(wad)
    wad.file.close()
    wfull = wobs.groupby(["donor", "cell_type"], observed=True).size().unstack(fill_value=0).sort_index()
    wmeta = (
        wobs.groupby("donor")
        .agg(age=("age", "first"), sex_male=("sex_male", "first"), platform_nextseq=("platform_nextseq", "first"))
        .loc[wfull.index]
    )
    w_p6 = filter_present_ge(wfull, 6)
    strata.append(
        {
            "species": "C. elegans",
            "organ": "Whole animal",
            "cohort": "Gao_GSE229022",
            "analysis": "all_annotated_cells",
            "counts_p6": w_p6,
            "counts_ge50": w_p6.loc[:, w_p6.sum(axis=0) >= 50].copy(),
            "metadata": wmeta,
            "has_platform": True,
            "is_primary": True,
        }
    )
    wobs_som = wobs[~wobs["cell_type"].isin(WORM_REPRODUCTIVE)].copy()
    wfull_som = wobs_som.groupby(["donor", "cell_type"], observed=True).size().unstack(fill_value=0).sort_index()
    w_som_p6 = filter_present_ge(wfull_som, 6)
    strata.append(
        {
            "species": "C. elegans",
            "organ": "Whole animal",
            "cohort": "Gao_GSE229022",
            "analysis": "somatic_only",
            "counts_p6": w_som_p6,
            "counts_ge50": w_som_p6.loc[:, w_som_p6.sum(axis=0) >= 50].copy(),
            "metadata": wmeta.loc[w_som_p6.index],
            "has_platform": True,
            "is_primary": False,
        }
    )
    summary_meta["worm"] = {
        "n_pools": int(w_p6.shape[0]),
        "n_total_cells": int(wfull.to_numpy().sum()),
        "n_cell_types_all": int(wfull.shape[1]),
        "n_cell_types_present6": int(w_p6.shape[1]),
        "n_cell_types_present6_sum50": int((w_p6.sum(axis=0) >= 50).sum()),
        "n_somatic_cell_types": int(w_som_p6.shape[1]),
        "reproductive_cells": int(wfull.loc[:, wfull.columns.isin(WORM_REPRODUCTIVE)].to_numpy().sum()),
    }

    # 3. Zebrafish testis: Primary 6-fish fine cell types (5 to 22 months)
    zad6 = ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad", backed="r")
    zobs6 = zfish_standard_obs(zad6, "5-22mo")
    zad6.file.close()
    zobs6_fine = zobs6[zobs6["fine_cell_type"] != "unknown"].copy()
    z6_full = zobs6_fine.groupby(["donor", "fine_cell_type"], observed=True).size().unstack(fill_value=0).sort_index()
    z6_meta = (
        zobs6_fine.groupby("donor")
        .agg(
            age=("age", "first"),
            sex_male=("sex_male", "first"),
            incomplete_spermatogenesis=("incomplete_spermatogenesis", "first"),
        )
        .loc[z6_full.index]
    )
    z6_p6 = filter_present_ge(z6_full, 6)
    strata.append(
        {
            "species": "Zebrafish",
            "organ": "Testis",
            "cohort": "GSE275361",
            "analysis": "fine_5_to_22_months",
            "counts_p6": z6_p6,
            "counts_ge50": z6_p6.loc[:, z6_p6.sum(axis=0) >= 50].copy(),
            "metadata": z6_meta,
            "has_platform": False,
            "is_primary": True,
        }
    )

    # 3b. Zebrafish reconciliation: 6-fish broad and 10-fish broad (5-22mo + 27mo)
    zobs6_broad = zobs6[zobs6["broad_cell_type"] != "unknown"].copy()
    z6_broad = zobs6_broad.groupby(["donor", "broad_cell_type"], observed=True).size().unstack(fill_value=0).sort_index()
    z6_broad_p6 = filter_present_ge(z6_broad, 6)
    strata.append(
        {
            "species": "Zebrafish",
            "organ": "Testis",
            "cohort": "GSE275361",
            "analysis": "broad_6fish_5_to_22_months",
            "counts_p6": z6_broad_p6,
            "counts_ge50": z6_broad_p6.loc[:, z6_broad_p6.sum(axis=0) >= 50].copy(),
            "metadata": z6_meta.loc[z6_broad_p6.index],
            "has_platform": False,
            "is_primary": False,
        }
    )

    z27_path = AGESCAPE_DATA / "zebrafish-testis-27mo.h5ad"
    if z27_path.exists():
        zad27 = ad.read_h5ad(z27_path, backed="r")
        zobs27 = zfish_standard_obs(zad27, "27mo")
        zad27.file.close()
        zobs10 = pd.concat([zobs6, zobs27], ignore_index=False)
        zobs10_broad = zobs10[zobs10["broad_cell_type"] != "unknown"].copy()
        z10_broad = (
            zobs10_broad.groupby(["donor", "broad_cell_type"], observed=True).size().unstack(fill_value=0).sort_index()
        )
        z10_meta = (
            zobs10_broad.groupby("donor")
            .agg(
                age=("age", "first"),
                sex_male=("sex_male", "first"),
                incomplete_spermatogenesis=("incomplete_spermatogenesis", "first"),
            )
            .loc[z10_broad.index]
        )
        z10_p6 = filter_present_ge(z10_broad, 6)
        strata.append(
            {
                "species": "Zebrafish",
                "organ": "Testis",
                "cohort": "GSE275361",
                "analysis": "broad_10fish_all_ages_reconciliation",
                "counts_p6": z10_p6,
                "counts_ge50": z10_p6.loc[:, z10_p6.sum(axis=0) >= 50].copy(),
                "metadata": z10_meta,
                "has_platform": False,
                "is_primary": False,
            }
        )
        summary_meta["zebrafish_10fish"] = {
            "n_fish": int(z10_p6.shape[0]),
            "n_total_cells": int(z10_broad.to_numpy().sum()),
            "germ_cells": int(z10_broad["germ cells"].sum()),
            "germ_fraction": float(z10_broad["germ cells"].sum() / z10_broad.to_numpy().sum()),
        }

    summary_meta["zebrafish_6fish"] = {
        "n_fish": int(z6_p6.shape[0]),
        "n_total_cells": int(z6_full.to_numpy().sum()),
        "n_cell_types_all": int(z6_full.shape[1]),
        "n_cell_types_present6": int(z6_p6.shape[1]),
        "n_cell_types_present6_sum50": int((z6_p6.sum(axis=0) >= 50).sum()),
        "excluded_cell_types": {
            c: {"n_present": int((z6_full[c] > 0).sum()), "total_cells": int(z6_full[c].sum())}
            for c in sorted(set(z6_full.columns) - set(z6_p6.columns))
        },
        "germ_cells": int(z6_broad["germ cells"].sum()),
        "germ_fraction_of_known": float(z6_broad["germ cells"].sum() / z6_broad.to_numpy().sum()),
    }
    return strata, summary_meta


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    strata, summary_meta = load_stratum_counts()
    all_rows: list[dict] = []

    for st in strata:
        sp_name = st["species"]
        organ = st["organ"]
        cohort = st["cohort"]
        analysis = st["analysis"]
        counts_p6 = st["counts_p6"]
        counts_ge50 = st["counts_ge50"]
        meta = st["metadata"]
        has_platform = st["has_platform"]

        # 1. Primary CLR (pseudocount = 0.5) and pseudocount sensitivities (0.1, 1.0)
        for pc in (0.5, 0.1, 1.0):
            sens_label = "primary_clr_0.5" if pc == 0.5 else f"sensitivity_clr_{pc:.1f}"
            all_rows.extend(
                run_clr_regressions(
                    counts_p6, meta, sp_name, organ, cohort, analysis, sens_label, pseudocount=pc
                )
            )
            if counts_ge50.shape[1] != counts_p6.shape[1]:
                all_rows.extend(
                    run_clr_regressions(
                        counts_ge50,
                        meta,
                        sp_name,
                        organ,
                        cohort,
                        analysis,
                        f"{sens_label}_total_cells_ge_50",
                        pseudocount=pc,
                    )
                )

        # 2. Additional regression sensitivities matching Mouse/Human/extended pipelines
        all_rows.extend(
            run_clr_regressions(
                counts_p6,
                meta,
                sp_name,
                organ,
                cohort,
                analysis,
                "observed_records_only",
                pseudocount=0.5,
                observed_only=True,
            )
        )
        all_rows.extend(
            run_clr_regressions(
                counts_p6,
                meta,
                sp_name,
                organ,
                cohort,
                analysis,
                "total_cells_adjusted",
                pseudocount=0.5,
                include_log_total=True,
            )
        )
        if has_platform:
            all_rows.extend(
                run_clr_regressions(
                    counts_p6,
                    meta,
                    sp_name,
                    organ,
                    cohort,
                    analysis,
                    "platform_adjusted",
                    pseudocount=0.5,
                    include_platform=True,
                )
            )
            all_rows.extend(
                run_clr_regressions(
                    counts_p6,
                    meta,
                    sp_name,
                    organ,
                    cohort,
                    analysis,
                    "platform_and_total_adjusted",
                    pseudocount=0.5,
                    include_log_total=True,
                    include_platform=True,
                )
            )
        if analysis == "broad_10fish_all_ages_reconciliation":
            all_rows.extend(
                run_clr_regressions(
                    counts_p6,
                    meta,
                    sp_name,
                    organ,
                    cohort,
                    analysis,
                    "spermatogenesis_adjusted",
                    pseudocount=0.5,
                    include_spermatogenesis=True,
                )
            )

        # 3. Propeller + limma empirical Bayes checks (primary design AND total-adjusted design)
        all_rows.extend(
            run_propeller_limma(
                counts_p6,
                meta,
                sp_name,
                organ,
                cohort,
                analysis,
                "propeller_limma_primary_design",
                pseudocount=0.5,
                include_log_total=False,
                include_platform=False,
            )
        )
        all_rows.extend(
            run_propeller_limma(
                counts_p6,
                meta,
                sp_name,
                organ,
                cohort,
                analysis,
                "propeller_limma_total_adjusted",
                pseudocount=0.5,
                include_log_total=True,
                include_platform=has_platform,
            )
        )

    df = pd.DataFrame(all_rows)
    df["is_primary_stratum"] = df["analysis"].isin(["all_annotated_cells", "fine_5_to_22_months"])
    df["is_primary_test"] = df["is_primary_stratum"] & (df["sensitivity"] == "primary_clr_0.5")

    # Cross-check committed R propeller output if available
    prop_r = pd.read_csv(RESULTS / "reviewer-phase1-composition-propeller.csv")
    prop_r_sub = prop_r[prop_r["species"].isin(["Fruit fly", "C. elegans", "Zebrafish"])].copy()
    prop_r_sub["analysis"] = np.where(
        prop_r_sub["species"] == "Zebrafish",
        "broad_10fish_all_ages_reconciliation",
        "all_annotated_cells",
    )
    prop_r_sub = prop_r_sub.rename(
        columns={
            "effect_logit": "r_propeller_effect_logit",
            "se": "r_propeller_se",
            "pvalue": "r_propeller_pvalue",
            "q_bh": "r_propeller_q_within_bh",
            "df_total": "r_propeller_df_total",
        }
    )[
        [
            "species",
            "organ",
            "cohort",
            "analysis",
            "cell_type",
            "r_propeller_effect_logit",
            "r_propeller_se",
            "r_propeller_pvalue",
            "r_propeller_q_within_bh",
            "r_propeller_df_total",
        ]
    ]
    df = df.merge(prop_r_sub, on=["species", "organ", "cohort", "analysis", "cell_type"], how="left")

    # Global BH family combining the 431 existing Mouse + Human + NMR primary tests with the 193 new primary tests (m = 624)
    existing_431 = pd.read_csv(RESULTS / "reviewer-phase1-composition-global-bh.csv")
    assert len(existing_431) == 431, f"Expected 431 existing rows, got {len(existing_431)}"
    new_primary = df[df["is_primary_test"]].copy()
    assert len(new_primary) == 193, f"Expected 193 primary new tests, got {len(new_primary)}"

    combined_pvals = np.concatenate(
        [existing_431["pvalue"].to_numpy(dtype=float), new_primary["pvalue"].to_numpy(dtype=float)]
    )
    combined_q = multipletests(combined_pvals, method="fdr_bh")[1]
    n_existing = len(existing_431)

    df["q_global_bh_624"] = np.nan
    df.loc[new_primary.index, "q_global_bh_624"] = combined_q[n_existing:]
    df["significant_within_05"] = df["q_within_bh"] < 0.05
    df["significant_global_624_05"] = df["q_global_bh_624"] < 0.05

    out_csv = RESULTS / "rev3-composition-fly-worm-zebrafish.csv"
    df.to_csv(out_csv, index=False)

    # Also save the full 624-row global BH table so every existing + new global q-value traces directly to disk
    existing_rows = existing_431.copy()
    existing_rows["source_family"] = "existing_431_mouse_human_nmr"
    existing_rows["q_within_bh"] = existing_rows["q_bh"]
    existing_rows["q_global_bh_431"] = existing_rows["q_global_bh"]
    existing_rows["q_global_bh_624"] = combined_q[:n_existing]
    existing_rows["significant_global_624_05"] = existing_rows["q_global_bh_624"] < 0.05

    new_rows_624 = df.loc[
        new_primary.index,
        [
            "species",
            "organ",
            "cohort",
            "cell_type",
            "sensitivity",
            "model",
            "effect_type",
            "effect",
            "se",
            "pvalue",
            "n_donors",
            "n_present",
            "total_cells",
            "q_within_bh",
            "q_global_bh_624",
            "significant_global_624_05",
        ],
    ].copy()
    new_rows_624["q_bh"] = new_rows_624["q_within_bh"]
    new_rows_624["source_family"] = "new_193_fly_worm_zebrafish"
    new_rows_624["q_global_bh_431"] = np.nan
    new_rows_624["significant_global_05"] = np.nan

    global_624 = pd.concat([existing_rows, new_rows_624], ignore_index=True)
    global_624.to_csv(RESULTS / "rev3-composition-global-bh-624.csv", index=False)

    # Print structured verification summary
    summary = {
        "metadata": summary_meta,
        "total_rows_written": int(len(df)),
        "primary_tests_count": int(len(new_primary)),
        "global_family_count": int(len(global_624)),
        "primary_by_stratum": [],
        "existing_431_under_624": {
            "sig_under_431": int((existing_rows["q_global_bh_431"] < 0.05).sum()),
            "sig_under_624": int(existing_rows["significant_global_624_05"].sum()),
        },
        "total_global_624_sig": int(global_624["significant_global_624_05"].sum()),
    }
    for (sp_name, organ, analysis), grp in df[df["is_primary_test"]].groupby(
        ["species", "organ", "analysis"], sort=False
    ):
        summary["primary_by_stratum"].append(
            {
                "species": sp_name,
                "organ": organ,
                "analysis": analysis,
                "n_tests": int(len(grp)),
                "sig_within_05": int((grp["q_within_bh"] < 0.05).sum()),
                "sig_global_624_05": int((grp["q_global_bh_624"] < 0.05).sum()),
                "min_pvalue": float(grp["pvalue"].min()),
                "min_q_within": float(grp["q_within_bh"].min()),
                "min_q_global_624": float(grp["q_global_bh_624"].min()),
            }
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
