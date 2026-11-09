#!/usr/bin/env python3
"""Rev2 fix 7: study-aware composition analysis for Agescape.

Human lung (Lung_Integrated, 119 donors age >= 20, 44 cell types) plus
Lung_CellRef replication analogues and a Tabula Muris Senis (mouse) design audit.

CPU-only reanalysis of files already in the repo. No new data, no new
experiments, no census re-query, no .tex or manuscript edits.

What this script does
---------------------
1. Derives contributing-study proxy groups from donor-ID prefixes (the
   Sikkema et al. integrated atlas carries no study column in-repo; the
   prefix is the only study signal, same proxy as
   scripts/run_d3_human_composition_sensitivity.py). Fits the club-cell and
   multiciliated-cell CLR slopes WITHIN each study group, then combines
   them by fixed-effect (inverse-variance) and DerSimonian-Laird
   random-effects meta-analysis, with Cochran Q and I^2.
2. Pooled models with study fixed effects, assay fixed effects (assay
   coverage is partial: 66/119 donors, from
   results/h3-human-lung-donor-assays.csv), log total recovered cells, and
   a sampling-site proxy (airway- vs alveolar-epithelial balance).
3. Re-analyses the Phase 2 age-aggregated QC-exclusion / unassigned-cell
   rates (results/reviewer-phase2-annotation-exclusion-by-age.csv) with
   cell-count weighting and influence checks, maps the age-group excluded
   fraction onto donors as an ecological covariate, and flags the
   club-cell marker bug in results/reviewer-phase2-marker-stability.csv
   (CD79A, a B-cell marker, was assigned to club cells because "club cell"
   contains the substring "b cell"; the correct SCGB1A1 test was never run).
4. Zero-sum checks: CLR slope sum, which types gain share (raw fractions,
   young vs old), leave-one-type-out CLR refits, compositional balances,
   and raw-proportion (non-CLR) regressions.
5. Mouse: design audit of TMS lung mice (sequencing wave x age x sex x
   pooling) from results/m1-animal-tissue-design.csv, a donor-level
   fraction comparison from results/m2-composition-vs-age.csv (only 3m vs
   18m lung mice exist there), and citation of the published summary
   slopes. No donor-level mouse count refit is possible in-repo
   (data/interim/m1-qc-cleaned.h5ad is absent from this checkout).
6. Lung_CellRef: club and multiciliated cells are absent (0 donors), so
   only the nearest analogues (ciliated cell, secretory cells) are tested.

All findings are reported strictly as associations among fractions of
RECOVERED ANNOTATED cells, never as absolute tissue cell numbers.

Output: results/reviewer-rev2-composition.csv (long format; one row per
estimate). Columns: analysis, species, cohort, cell_type, group, model,
n_donors, n_studies, age_mean, age_sd, estimate, se, ci95_low, ci95_high,
pvalue, q_bh, note. For descriptive/design rows, estimate/se/ci carry the
described quantity named in `note`.
"""

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm, spearmanr
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.run_h2_composition import compute_clr_table, load_human_counts

RESULTS = REPO / "results"
OUT = RESULTS / "reviewer-rev2-composition.csv"
PREFIX = "homosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747"
FOCAL = ["club cell", "multiciliated epithelial cell",
         "multiciliated columnar cell of tracheobronchial tree",
         "pulmonary alveolar type 2 cell"]
# Airway vs alveolar epithelial sets (sampling-site proxy, derived from counts only).
AIRWAY_EPI = ["club cell", "multiciliated epithelial cell",
              "multiciliated columnar cell of tracheobronchial tree",
              "respiratory basal cell", "bronchial goblet cell",
              "tracheobronchial goblet cell", "respiratory tract hillock cell",
              "brush cell of tracheobronchial tree", "ionocyte",
              "pulmonary neuroendocrine cell", "tracheobronchial serous cell",
              "serous secreting cell", "mucus secreting cell"]
# Same proxy with the three focal types removed, so its age trend cannot be
# driven by the focal declines themselves.
AIRWAY_EPI_EXCL_FOCAL = ["respiratory basal cell", "bronchial goblet cell",
                         "tracheobronchial goblet cell", "respiratory tract hillock cell",
                         "brush cell of tracheobronchial tree", "ionocyte",
                         "pulmonary neuroendocrine cell", "tracheobronchial serous cell",
                         "serous secreting cell", "mucus secreting cell"]
ALVEOLAR_EPI = ["pulmonary alveolar type 1 cell", "pulmonary alveolar type 2 cell",
                "pulmonary alveolar epithelial cell"]

ROWS = []


def add(analysis, species="", cohort="", cell_type="", group="", model="", n_donors=np.nan,
        n_studies=np.nan, age_mean=np.nan, age_sd=np.nan, estimate=np.nan, se=np.nan,
        ci_low=np.nan, ci_high=np.nan, pvalue=np.nan, q_bh=np.nan, note=""):
    ROWS.append({"analysis": analysis, "species": species, "cohort": cohort,
                 "cell_type": cell_type, "group": group, "model": model,
                 "n_donors": n_donors, "n_studies": n_studies, "age_mean": age_mean,
                 "age_sd": age_sd, "estimate": estimate, "se": se, "ci95_low": ci_low,
                 "ci95_high": ci_high, "pvalue": pvalue, "q_bh": q_bh, "note": note})


def ols_slope(y, X_df):
    """OLS slope on 'age' with 95% CI; returns dict or None if rank-deficient."""
    X = sm.add_constant(X_df.astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or len(y) <= X.shape[1] + 1:
        return None
    fit = sm.OLS(y.astype(float), X).fit()
    ci = fit.conf_int().loc["age"]
    return {"slope": float(fit.params["age"]), "se": float(fit.bse["age"]),
            "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
            "p": float(fit.pvalues["age"]), "n": int(len(y))}


def study_of(donor_id):
    raw = str(donor_id).replace(PREFIX, "")
    m = re.match(r"^([A-Za-z]+|\d+-\d+-[A-Za-z]+)", raw)
    return m.group(1) if m else "other"


def dl_meta(ests):
    """Inverse-variance fixed effect + DerSimonian-Laird random effects."""
    ests = [e for e in ests if np.isfinite(e["slope"]) and np.isfinite(e["se"]) and e["se"] > 0]
    k = len(ests)
    if k < 2:
        return None
    b = np.array([e["slope"] for e in ests])
    v = np.array([e["se"] ** 2 for e in ests])
    w = 1.0 / v
    fe = float((w * b).sum() / w.sum())
    fe_se = float(np.sqrt(1.0 / w.sum()))
    q = float((w * (b - fe) ** 2).sum())
    df = k - 1
    i2 = float(max(0.0, (q - df) / q)) if q > 0 else 0.0
    tau2 = float(max(0.0, (q - df) / (w.sum() - (w ** 2).sum() / w.sum())))
    wr = 1.0 / (v + tau2)
    re_ = float((wr * b).sum() / wr.sum())
    re_se = float(np.sqrt(1.0 / wr.sum()))
    return {"k": k, "fe": fe, "fe_se": fe_se, "q": q, "df": df, "i2": i2,
            "tau2": tau2, "re": re_, "re_se": re_se}


def main():
    counts, meta = load_human_counts(RESULTS)
    meta["study"] = [study_of(d) for d in meta.index]
    meta["total_cells"] = counts.sum(axis=1)
    clr = compute_clr_table(counts)
    props = counts.div(counts.sum(axis=1), axis=0)
    assay = (pd.read_csv(RESULTS / "h3-human-lung-donor-assays.csv")
             [["donor_id", "assay"]].drop_duplicates().set_index("donor_id")["assay"])
    meta["assay"] = meta.index.map(assay)

    # ------------------------------------------------------------------ 1. design
    for study, part in meta.groupby("study"):
        assays = part["assay"].value_counts(dropna=True)
        dom = f"{assays.index[0]} ({assays.iloc[0]}/{len(part)})" if len(assays) else "no assay in-repo"
        add("study_design", "Human", "Lung_Integrated", group=study, model="donor_id_prefix_group",
            n_donors=len(part), age_mean=part["age"].mean(), age_sd=part["age"].std(),
            estimate=part["age"].min(), se=part["age"].max(),
            ci_low=part["total_cells"].median(), ci_high=(part["sex"] == 1).mean(),
            note=(f"age range {part['age'].min():.0f} to {part['age'].max():.0f}; "
                  f"median recovered annotated cells/donor {part['total_cells'].median():.0f}; "
                  f"dominant assay in-repo: {dom}; "
                  f"estimate=age_min, se=age_max, ci95_low=median total cells, ci95_high=male fraction"))
    add("study_design", "Human", "Lung_Integrated", model="assay_coverage",
        n_donors=int(meta["assay"].notna().sum()),
        note=("Assay known for 66 of 119 donors (10x 3' v2: 19; 10x 5' v1: 25; "
              "10x 3' v3: 15; 10x 5' v2: 4; Seq-Well: 3). The other 53 donors have no "
              "assay value in any in-repo file; assay-adjusted models below use n=66. "
              "No smoking-history or sampling-site (airway vs parenchyma) variable "
              "exists in any in-repo file; site is proxied only by the airway/alveolar "
              "epithelial balance derived from the counts themselves."))

    # per-study mean recovered fractions and presence (the key confounding display)
    for ct in ["club cell", "multiciliated epithelial cell",
               "multiciliated columnar cell of tracheobronchial tree",
               "pulmonary alveolar type 2 cell", "nasal mucosa goblet cell",
               "respiratory basal cell"]:
        for study, part in meta.groupby("study"):
            frac = props.loc[part.index, ct]
            add("study_mean_fractions", "Human", "Lung_Integrated", cell_type=ct, group=study,
                model="mean_observed_fraction_and_presence", n_donors=len(part),
                age_mean=part["age"].mean(), age_sd=part["age"].std(),
                estimate=float(frac.mean()), se=float(frac.std()),
                ci_low=float((counts.loc[part.index, ct] > 0).mean()),
                note=(f"mean recovered fraction {frac.mean():.4f}; present in "
                      f"{int((counts.loc[part.index, ct] > 0).sum())}/{len(part)} donors "
                      f"(ci95_low column = presence fraction); se column = SD of the fraction"))
        vals = [props.loc[meta[meta['study'] == s].index, ct].mean()
                for s in meta["study"].unique() if (meta["study"] == s).sum() >= 5]
        pos = [v for v in vals if v > 0]
        n_zero = sum(1 for v in vals if v == 0)
        fold = max(pos) / min(pos) if pos else float("nan")
        add("study_mean_fractions", "Human", "Lung_Integrated", cell_type=ct,
            model="fold_range_across_studies_n_ge_5", n_donors=119,
            estimate=float(fold),
            note=(f"mean recovered fraction ranges {min(pos):.5f} to {max(pos):.4f} across "
                  f"study groups with n>=5 and nonzero mean, a {fold:.0f}-fold range; "
                  f"{n_zero} further group(s) with n>=5 have exactly zero recovered cells of "
                  f"this type; estimate is that fold range"))

    # ------------------------------------------------------- 2. per-study slopes
    major = meta["study"].value_counts()
    major_studies = major[major >= 5].index.tolist()
    per_study = {ct: [] for ct in FOCAL}
    for ct in FOCAL:
        for study, part in meta.groupby("study"):
            if len(part) < 5 or part["age"].std() < 3:
                add("per_study_slope", "Human", "Lung_Integrated", cell_type=ct, group=study,
                    model="within_study_ols_age_sex", n_donors=len(part),
                    age_mean=part["age"].mean(), age_sd=part["age"].std(),
                    note="not estimated: n<5 or within-study age SD<3y (no usable age contrast)")
                continue
            cols = ["age"] + (["sex"] if part["sex"].nunique() > 1 else [])
            r = ols_slope(clr.loc[part.index, ct], part[cols])
            if r is None:
                add("per_study_slope", "Human", "Lung_Integrated", cell_type=ct, group=study,
                    model="within_study_ols_age_sex", n_donors=len(part),
                    age_mean=part["age"].mean(), age_sd=part["age"].std(),
                    note="not estimated: rank-deficient design")
                continue
            r["study"] = study
            per_study[ct].append(r)
            add("per_study_slope", "Human", "Lung_Integrated", cell_type=ct, group=study,
                model="within_study_ols_age_sex" if "sex" in cols else "within_study_ols_age_only",
                n_donors=r["n"], age_mean=part["age"].mean(), age_sd=part["age"].std(),
                estimate=r["slope"], se=r["se"], ci_low=r["ci_low"], ci_high=r["ci_high"],
                pvalue=r["p"],
                note=(f"age {part['age'].min():.0f} to {part['age'].max():.0f}; "
                      f"slope sign {'negative' if r['slope'] < 0 else 'positive'}"))
        # BH across the estimated per-study p-values for this cell type
        est = per_study[ct]
        if est:
            q = multipletests([e["p"] for e in est], method="fdr_bh")[1]
            for e, qq in zip(est, q):
                for row in ROWS:
                    if (row["analysis"] == "per_study_slope" and row["cell_type"] == ct
                            and row["group"] == e["study"] and pd.notna(row["estimate"])):
                        row["q_bh"] = float(qq)
            n_neg = sum(e["slope"] < 0 for e in est)
            n_sig = sum(e["p"] < 0.05 for e in est)
            n_sig_neg = sum(e["p"] < 0.05 and e["slope"] < 0 for e in est)
            add("per_study_summary", "Human", "Lung_Integrated", cell_type=ct,
                model="per_study_count", n_donors=int(sum(e["n"] for e in est)),
                n_studies=len(est),
                note=(f"{n_neg}/{len(est)} estimable studies negative; {n_sig_neg}/{len(est)} "
                      f"negative and p<0.05; {n_sig}/{len(est)} p<0.05 either direction"))
            m = dl_meta(est)
            if m:
                add("meta_combine", "Human", "Lung_Integrated", cell_type=ct,
                    model="fixed_effect_ivw", n_donors=int(sum(e["n"] for e in est)),
                    n_studies=m["k"], estimate=m["fe"], se=m["fe_se"],
                    ci_low=m["fe"] - 1.96 * m["fe_se"], ci_high=m["fe"] + 1.96 * m["fe_se"],
                    pvalue=float(2 * norm.sf(abs(m["fe"] / m["fe_se"]))),
                    note=(f"Cochran Q={m['q']:.2f} on {m['df']} df; I2={m['i2']:.2f}; "
                          f"tau2={m['tau2']:.5f}"))
                add("meta_combine", "Human", "Lung_Integrated", cell_type=ct,
                    model="random_effects_dl", n_donors=int(sum(e["n"] for e in est)),
                    n_studies=m["k"], estimate=m["re"], se=m["re_se"],
                    ci_low=m["re"] - 1.96 * m["re_se"], ci_high=m["re"] + 1.96 * m["re_se"],
                    pvalue=float(2 * norm.sf(abs(m["re"] / m["re_se"]))),
                    note=(f"Cochran Q={m['q']:.2f} on {m['df']} df; I2={m['i2']:.2f}; "
                          f"tau2={m['tau2']:.5f}"))

    # ------------------------------------------------------- 3. pooled models
    base = pd.read_csv(RESULTS / "h2-composition-human-lung.csv").set_index("cell_type")
    meta["study_grouped"] = meta["study"].apply(lambda s: s if s in major_studies else "other")
    meta["log_total"] = np.log(meta["total_cells"])
    # sampling-site proxy: log airway / alveolar epithelial counts
    airway = counts[[c for c in AIRWAY_EPI if c in counts.columns]].sum(axis=1)
    alveolar = counts[[c for c in ALVEOLAR_EPI if c in counts.columns]].sum(axis=1)
    meta["log_airway_alveolar"] = np.log((airway + 0.5) / (alveolar + 0.5))
    meta["airway_share"] = airway / meta["total_cells"]
    airway_xf = counts[[c for c in AIRWAY_EPI_EXCL_FOCAL if c in counts.columns]].sum(axis=1)
    meta["log_airway_alveolar_excl_focal"] = np.log((airway_xf + 0.5) / (alveolar + 0.5))
    for ct in FOCAL:
        y = clr[ct].loc[meta.index]
        r = ols_slope(y, meta[["age", "sex"]])
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct, model="baseline_age_sex_reproduced",
            n_donors=119, estimate=r["slope"], se=r["se"], ci_low=r["ci_low"], ci_high=r["ci_high"],
            pvalue=r["p"], note=f"primary file value slope {base.loc[ct, 'slope_per_year_or_month']:.6f}, "
                                 f"q {base.loc[ct, 'q_bh']:.6f}")
        dummies = pd.get_dummies(meta["study_grouped"], drop_first=True, dtype=float)
        r = ols_slope(y, pd.concat([meta[["age", "sex"]], dummies], axis=1))
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct, model="study_fixed_effects",
            n_donors=119, n_studies=len(major_studies) + 1, estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note="study dummies for groups with n>=5, smaller groups pooled as 'other' (same grouping as D3)")
        r = ols_slope(y, pd.concat([meta[["age", "sex", "log_total"]], dummies], axis=1))
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct,
            model="study_fixed_effects_plus_log_total_cells",
            n_donors=119, n_studies=len(major_studies) + 1, estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"], note="")
        r = ols_slope(y, pd.concat([meta[["age", "sex", "log_airway_alveolar"]], dummies], axis=1))
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct,
            model="study_fixed_effects_plus_airway_alveolar_balance",
            n_donors=119, n_studies=len(major_studies) + 1, estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note="airway/alveolar balance is a counts-derived sampling-site proxy, not a measured site; includes the focal types")
        r = ols_slope(y, pd.concat([meta[["age", "sex", "log_airway_alveolar_excl_focal"]], dummies], axis=1))
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct,
            model="study_fixed_effects_plus_airway_alveolar_balance_excl_focal",
            n_donors=119, n_studies=len(major_studies) + 1, estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note="same proxy with the three focal types removed from the airway set, so the covariate cannot be driven by the focal declines")
        sub = meta[meta["assay"].notna()]
        rb = ols_slope(y.loc[sub.index], sub[["age", "sex"]])
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct,
            model="baseline_age_sex_assay_known_subset", n_donors=len(sub),
            estimate=rb["slope"], se=rb["se"], ci_low=rb["ci_low"], ci_high=rb["ci_high"],
            pvalue=rb["p"], note="same 66 donors as the assay-adjusted model, for comparability")
        ad = pd.get_dummies(sub["assay"], drop_first=True, dtype=float)
        ra = ols_slope(y.loc[sub.index], pd.concat([sub[["age", "sex"]], ad], axis=1))
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct,
            model="assay_fixed_effects", n_donors=len(sub), n_studies=int(sub["assay"].nunique()),
            estimate=ra["slope"], se=ra["se"], ci_low=ra["ci_low"], ci_high=ra["ci_high"],
            pvalue=ra["p"], note=f"assay dummies, {sub['assay'].nunique()} assays; n=66 subset only")
        # between-study (ecological) slope: study means
        means = pd.DataFrame({"y": clr[ct].groupby(meta["study"]).mean(),
                              "age": meta.groupby("study")["age"].mean(),
                              "n": meta.groupby("study").size()})
        means = means[means["n"] >= 5]
        rb2 = ols_slope(means["y"], means[["age"]])
        add("pooled_model", "Human", "Lung_Integrated", cell_type=ct,
            model="between_study_ecological_means", n_donors=int(means["n"].sum()),
            n_studies=len(means), estimate=rb2["slope"], se=rb2["se"],
            ci_low=rb2["ci_low"], ci_high=rb2["ci_high"], pvalue=rb2["p"],
            note="OLS across study-group means (weighted equally); shows how much of the pooled signal is between-study")
    # does the sampling-site proxy itself change with age / differ by study?
    for outcome, col in (("log airway/alveolar balance", "log_airway_alveolar"),
                         ("log airway/alveolar balance excluding focal types",
                          "log_airway_alveolar_excl_focal"),
                         ("airway epithelial share of recovered cells", "airway_share")):
        r = ols_slope(meta[col], meta[["age", "sex"]])
        add("sampling_site_proxy", "Human", "Lung_Integrated", cell_type="(all epithelial)",
            model="proxy_vs_age", n_donors=119, estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note=f"outcome is {outcome}; a strong age trend here would mean sampling site is age-confounded")
    add("sampling_site_proxy", "Human", "Lung_Integrated", cell_type="(all epithelial)",
        model="proxy_by_study", n_donors=119,
        note=("median log airway/alveolar balance by study: " +
              "; ".join(f"{s}: {v:.2f}" for s, v in
                        meta.groupby("study")["log_airway_alveolar"].median().sort_values().items()) +
              ". Excluding focal types: " +
              "; ".join(f"{s}: {v:.2f}" for s, v in
                        meta.groupby("study")["log_airway_alveolar_excl_focal"].median().sort_values().items())))

    # ------------------------------------------------------- 4. exclusion / unassigned
    excl = pd.read_csv(RESULTS / "reviewer-phase2-annotation-exclusion-by-age.csv")
    lung_ex = excl[(excl["cohort"] == "Lung_Integrated")].copy()
    lung_ex["age"] = lung_ex["age"].astype(float)
    for metric in ("fraction_excluded_qc", "fraction_unassigned_raw"):
        rho, rp = spearmanr(lung_ex["age"], lung_ex[metric])
        add("exclusion_model", "Human", "Lung_Integrated", cell_type="(all cells)",
            model=f"spearman_age_vs_{metric}", n_donors=int(lung_ex["n_raw_cells"].sum()),
            estimate=float(rho), pvalue=float(rp),
            note=f"unweighted Spearman across {len(lung_ex)} age-group points")
        X = sm.add_constant(lung_ex[["age"]].astype(float), has_constant="add")
        wfit = sm.WLS(lung_ex[metric], X, weights=lung_ex["n_raw_cells"]).fit()
        ci = wfit.conf_int().loc["age"]
        add("exclusion_model", "Human", "Lung_Integrated", cell_type="(all cells)",
            model=f"wls_age_vs_{metric}_weighted_by_n_raw", n_donors=int(lung_ex["n_raw_cells"].sum()),
            estimate=float(wfit.params["age"]), se=float(wfit.bse["age"]),
            ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), pvalue=float(wfit.pvalues["age"]),
            note="weights = n_raw_cells per age group; estimate is fraction points per year of age")
        ufit = sm.OLS(lung_ex[metric], X).fit()
        ci = ufit.conf_int().loc["age"]
        add("exclusion_model", "Human", "Lung_Integrated", cell_type="(all cells)",
            model=f"ols_age_vs_{metric}_unweighted", n_donors=int(lung_ex["n_raw_cells"].sum()),
            estimate=float(ufit.params["age"]), se=float(ufit.bse["age"]),
            ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), pvalue=float(ufit.pvalues["age"]),
            note="unweighted OLS across age-group points, for comparison with the weighted fit")
        keep = lung_ex[lung_ex["n_raw_cells"] >= 3000]
        wfit2 = sm.WLS(keep[metric], sm.add_constant(keep[["age"]].astype(float), has_constant="add"),
                       weights=keep["n_raw_cells"]).fit()
        ci = wfit2.conf_int().loc["age"]
        add("exclusion_model", "Human", "Lung_Integrated", cell_type="(all cells)",
            model=f"wls_age_vs_{metric}_drop_small_n_ages", n_donors=int(keep["n_raw_cells"].sum()),
            estimate=float(wfit2.params["age"]), se=float(wfit2.bse["age"]),
            ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), pvalue=float(wfit2.pvalues["age"]),
            note=f"drops {len(lung_ex) - len(keep)} age groups with n_raw<3000 (incl. the extreme ages 60/73/75)")
    # young vs old weighted mean fractions
    for label, mask in (("age 20-39", lung_ex["age"] <= 39), ("age 60-81", lung_ex["age"] >= 60)):
        part = lung_ex[mask]
        for metric in ("fraction_excluded_qc", "fraction_unassigned_raw"):
            wmean = float((part[metric] * part["n_raw_cells"]).sum() / part["n_raw_cells"].sum())
            add("exclusion_model", "Human", "Lung_Integrated", cell_type="(all cells)",
                model=f"weighted_mean_{metric}", group=label, n_donors=int(part["n_raw_cells"].sum()),
                estimate=wmean, note=f"n_raw_cells={int(part['n_raw_cells'].sum())}; estimate is the cell-weighted mean fraction")
    # annotated share of raw cells by age band (m9 annotated totals vs phase2 raw totals)
    ann_by_age = counts.sum(axis=1).groupby(meta["age"]).sum()
    raw_by_age = lung_ex.set_index("age")["n_raw_cells"]
    common = ann_by_age.index.intersection(raw_by_age.index)
    share = (ann_by_age.loc[common] / raw_by_age.loc[common])
    for label, ages in (("age 20-39", [a for a in common if a <= 39]),
                        ("age 40-59", [a for a in common if 40 <= a <= 59]),
                        ("age 60-81", [a for a in common if a >= 60])):
        num = ann_by_age.loc[ages].sum()
        den = raw_by_age.loc[ages].sum()
        add("exclusion_model", "Human", "Lung_Integrated", cell_type="(all cells)",
            model="annotated_share_of_raw_cells", group=label,
            n_donors=int(den), estimate=float(num / den),
            note=(f"annotated recovered cells {int(num)} / raw cells {int(den)} in this age band; "
                  "remainder = QC-excluded + unassigned + overlap between the two flags (overlap not separable in-repo)"))
    # ecological adjustment: map age-group excluded fraction onto donors
    excl_map = lung_ex.set_index("age")["fraction_excluded_qc"]
    ages_sorted = np.array(sorted(excl_map.index))
    frac_sorted = np.array([excl_map.loc[a] for a in ages_sorted])
    meta["excl_frac_agegroup"] = np.interp(meta["age"], ages_sorted, frac_sorted)
    meta["unassigned_frac_agegroup"] = np.interp(
        meta["age"], ages_sorted,
        np.array([lung_ex.set_index("age").loc[a, "fraction_unassigned_raw"] for a in ages_sorted]))
    for ct in FOCAL[:2]:
        y = clr[ct].loc[meta.index]
        r = ols_slope(y, meta[["age", "sex", "excl_frac_agegroup"]])
        add("exclusion_model", "Human", "Lung_Integrated", cell_type=ct,
            model="clr_age_sex_plus_agegroup_excluded_fraction", n_donors=119,
            estimate=r["slope"], se=r["se"], ci_low=r["ci_low"], ci_high=r["ci_high"],
            pvalue=r["p"],
            note="ECOLOGICAL covariate: every donor gets its age-group's excluded fraction; not a donor-level measurement")
        r = ols_slope(y, meta[["age", "sex", "unassigned_frac_agegroup"]])
        add("exclusion_model", "Human", "Lung_Integrated", cell_type=ct,
            model="clr_age_sex_plus_agegroup_unassigned_fraction", n_donors=119,
            estimate=r["slope"], se=r["se"], ci_low=r["ci_low"], ci_high=r["ci_high"],
            pvalue=r["p"], note="ECOLOGICAL covariate, as above")
    # presence with study fixed effects. The logit is separation-prone (some
    # study groups have all-present or all-absent donors), so convergence is
    # checked and a linear-probability fallback is reported instead when it fails.
    for ct in FOCAL[:2]:
        present = (counts[ct] > 0).astype(int).loc[meta.index]
        dummies = pd.get_dummies(meta["study_grouped"], drop_first=True, dtype=float)
        X = sm.add_constant(pd.concat([meta[["age", "sex"]], dummies], axis=1).astype(float), has_constant="add")
        converged = False
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = sm.Logit(present, X).fit(disp=0, maxiter=200)
            converged = bool(fit.mle_retvals.get("converged", False))
        except Exception:  # noqa: BLE001 - separation is a real outcome here
            fit = None
        if converged:
            add("exclusion_model", "Human", "Lung_Integrated", cell_type=ct,
                model="presence_logit_age_sex_study_fe", n_donors=119,
                estimate=float(fit.params["age"]), se=float(fit.bse["age"]),
                ci_low=float(fit.conf_int().loc["age"].iloc[0]),
                ci_high=float(fit.conf_int().loc["age"].iloc[1]),
                pvalue=float(fit.pvalues["age"]),
                note="log-odds of the type being observed at all in a donor, per year of age, study FE included")
        else:
            add("exclusion_model", "Human", "Lung_Integrated", cell_type=ct,
                model="presence_logit_age_sex_study_fe", n_donors=119,
                note=("logit did NOT converge (separation: some study groups are all-present or "
                      "all-absent for this type, e.g. multiciliated epithelial absent in 7/7 Patient "
                      "and 10/11 muc donors but present in 11/11 D and 13/13 T donors); no logit "
                      "estimate is reported"))
            lp = sm.OLS(present, X).fit()
            ci = lp.conf_int().loc["age"]
            add("exclusion_model", "Human", "Lung_Integrated", cell_type=ct,
                model="presence_linear_probability_age_sex_study_fe", n_donors=119,
                estimate=float(lp.params["age"]), se=float(lp.bse["age"]),
                ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]),
                pvalue=float(lp.pvalues["age"]),
                note="linear-probability fallback (OLS on 0/1 presence) with study FE; estimate is change in presence probability per year of age")
    # marker audit flag (values cited from the Phase 2 file; bug verified in the script source)
    add("marker_audit", "Human", "Lung_Integrated", cell_type="club cell",
        model="phase2_marker_stability_INVALID", n_donors=43,
        estimate=0.000094, pvalue=0.219078,
        note=("Phase 2 assigned marker CD79A to club cells (positive fraction 0.0034) because "
              "scripts/run_reviewer_phase2_annotation.py tests the ('b cell',) rule before the "
              "('club cell',) rule and 'club cell' contains 'b cell'. CD79A is a B-cell marker, "
              "so this row is not a club-cell identity test at all. The correct SCGB1A1 test was "
              "never run and cannot be run in-repo (no expression data in this checkout). "
              "Identity loss for club cells is therefore UNTESTED, not excluded."))
    add("marker_audit", "Human", "Lung_Integrated", cell_type="multiciliated epithelial cell",
        model="phase2_marker_stability_FOXJ1_valid", n_donors=18,
        estimate=0.000595, se=0.005577, ci_low=-0.011367, ci_high=0.012557, pvalue=0.916582,
        note=("FOXJ1 log-CP10K per year, valid marker (positive fraction 0.81). No age trend: "
              "per-cell FOXJ1 expression is stable, which argues against loss of ciliated "
              "identity within annotated cells, but does not rule out cells dropping out of "
              "annotation entirely. Source: results/reviewer-phase2-marker-stability.csv:159."))
    add("marker_audit", "Human", "Lung_Integrated",
        cell_type="multiciliated columnar cell of tracheobronchial tree",
        model="phase2_marker_stability_FOXJ1_valid", n_donors=80,
        estimate=-0.002713, se=0.002639, ci_low=-0.007968, ci_high=0.002542, pvalue=0.307107,
        note=("FOXJ1 log-CP10K per year, positive fraction 0.79. No age trend. "
              "Source: results/reviewer-phase2-marker-stability.csv:150."))

    # ------------------------------------------------------- 5. zero-sum checks
    slopes = base["slope_per_year_or_month"]
    add("zerosum", "Human", "Lung_Integrated", cell_type="(all 44 types)",
        model="clr_slope_sum_check", n_donors=119, estimate=float(slopes.sum()),
        note=f"CLR slopes must sum to ~0 by construction; observed sum {slopes.sum():.2e} across 44 types")
    pos = slopes[slopes > 0]
    neg = slopes[slopes < 0]
    add("zerosum", "Human", "Lung_Integrated", cell_type="(all 44 types)",
        model="positive_vs_negative_slope_mass", n_donors=119,
        estimate=float(pos.sum()), se=float(neg.sum()),
        note=(f"sum of positive CLR slopes {pos.sum():.4f} (estimate) vs sum of negative slopes "
              f"{neg.sum():.4f} (se column); the two focal declines are part of the negative mass"))
    for ct in base.sort_values("slope_per_year_or_month", ascending=False).index:
        if base.loc[ct, "slope_per_year_or_month"] > 0:
            add("zerosum_gainer", "Human", "Lung_Integrated", cell_type=ct,
                model="primary_clr_slope_positive", n_donors=119,
                estimate=float(base.loc[ct, "slope_per_year_or_month"]),
                se=float(base.loc[ct, "se"]), pvalue=float(base.loc[ct, "pvalue"]),
                q_bh=float(base.loc[ct, "q_bh"]),
                note="types gaining CLR share; the largest gainer is the mirror image of the focal declines")
    # young vs old raw fractions
    young = meta["age"] <= 39
    old = meta["age"] >= 60
    for ct in props.columns:
        yf = props.loc[young, ct].mean()
        of = props.loc[old, ct].mean()
        add("zerosum_raw_fractions", "Human", "Lung_Integrated", cell_type=ct,
            model="mean_observed_fraction_old_minus_young", n_donors=int(young.sum() + old.sum()),
            estimate=float(of - yf), se=float(yf), ci_low=float(of),
            note=(f"mean observed fraction young (20-39y, n={int(young.sum())}) {yf:.4f} (se column) "
                  f"vs old (60-81y, n={int(old.sum())}) {of:.4f} (ci95_low column); "
                  "estimate is old minus young in fraction points"))
    # leave-one-type-out CLR refits
    for drop, label in ((["club cell"], "drop_club"),
                        (["multiciliated epithelial cell"], "drop_multiciliated_epithelial"),
                        (["club cell", "multiciliated epithelial cell"], "drop_both_focal"),
                        (["pulmonary alveolar type 2 cell"], "drop_AT2"),
                        (["club cell", "multiciliated epithelial cell",
                          "multiciliated columnar cell of tracheobronchial tree"],
                         "drop_all_airway_ciliated_secretory_focal")):
        sub_counts = counts.drop(columns=[c for c in drop if c in counts.columns])
        sub_clr = compute_clr_table(sub_counts)
        for ct in FOCAL:
            if ct in drop:
                continue
            r = ols_slope(sub_clr[ct].loc[meta.index], meta[["age", "sex"]])
            add("zerosum_loto", "Human", "Lung_Integrated", cell_type=ct, model=f"clr_{label}",
                n_donors=119, estimate=r["slope"], se=r["se"], ci_low=r["ci_low"],
                ci_high=r["ci_high"], pvalue=r["p"],
                note=f"CLR recomputed without {', '.join(drop)}; primary slope was "
                     f"{base.loc[ct, 'slope_per_year_or_month']:.4f}")
    # compositional balances (log-ratios, pseudocount 0.5)
    for label, num, den in (
            ("club / AT2", ["club cell"], ["pulmonary alveolar type 2 cell"]),
            ("multiciliated epithelial / AT2", ["multiciliated epithelial cell"],
             ["pulmonary alveolar type 2 cell"]),
            ("club+multiciliated epithelial / AT2",
             ["club cell", "multiciliated epithelial cell"], ["pulmonary alveolar type 2 cell"]),
            ("club / multiciliated columnar",
             ["club cell"], ["multiciliated columnar cell of tracheobronchial tree"]),
            ("(club+multiciliated epithelial) / rest of recovered cells",
             ["club cell", "multiciliated epithelial cell"], None)):
        n_ = counts[num].sum(axis=1) + 0.5
        d_ = counts[den].sum(axis=1) + 0.5 if den else (meta["total_cells"] - counts[num].sum(axis=1)) + 0.5
        bal = np.log(n_ / d_)
        r = ols_slope(bal, meta[["age", "sex"]])
        add("zerosum_balance", "Human", "Lung_Integrated", cell_type=label,
            model="log_ratio_vs_age_sex", n_donors=119, estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note="log-ratio slope per year; a negative balance slope means the numerator loses share against that denominator specifically")
    # raw proportions (non-CLR confirmation)
    for ct in FOCAL[:2]:
        r = ols_slope(props[ct].loc[meta.index], meta[["age", "sex"]])
        add("zerosum_raw_proportions", "Human", "Lung_Integrated", cell_type=ct,
            model="observed_fraction_vs_age_sex", n_donors=119, estimate=r["slope"],
            se=r["se"], ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note="OLS on the raw observed fraction (not CLR); slope is fraction points per year")
        rho, rp = spearmanr(meta["age"], props[ct])
        add("zerosum_raw_proportions", "Human", "Lung_Integrated", cell_type=ct,
            model="observed_fraction_spearman_vs_age", n_donors=119,
            estimate=float(rho), pvalue=float(rp), note="")

    # ------------------------------------------------------- 6. mouse
    design = pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
    lung_design = design[design["tissue"] == "Lung"]
    add("mouse_design", "Mouse", "Lung", model="tms_lung_design_audit",
        n_donors=len(lung_design),
        note=("TMS lung mice by age: " +
              "; ".join(f"{a}: n={len(p)}, waves={sorted(p['sequencing_wave'].unique())}, "
                        f"sexes={sorted(p['sex'].unique())}, pooled={sorted(p['is_pooled_id'].unique())}, "
                        f"immune-endo-depleted={sorted(p['is_immune_endo_depleted'].unique())}"
                        for a, p in lung_design.groupby("age")) +
              ". Age is near-perfectly confounded with sequencing wave (1m/3m=10X_P, 18m/21m=MACA, "
              "30m=10X_P) and with sex at 21m (all female) and 30m (all male in the design file); "
              "the two 3m male lung samples are pooled IDs and immune/endothelial-depleted. "
              "No 24m lung mice appear in the design file."))
    m2 = pd.read_csv(RESULTS / "m2-composition-vs-age.csv")
    m2lung = m2[m2["tissue"] == "Lung"]
    for ct in ("club cell of bronchiole", "ciliated columnar cell of tracheobronchial tree"):
        sub = m2lung[m2lung["cell_ontology_class"] == ct]
        young_f = sub[sub["age"] == "3m"]["fraction"]
        old_f = sub[sub["age"] == "18m"]["fraction"]
        add("mouse_refit", "Mouse", "Lung", cell_type=ct,
            model="m2_fractions_3m_vs_18m_only", n_donors=int(m2lung["mouse.id"].nunique()),
            estimate=float(old_f.mean() - young_f.mean()) if len(young_f) and len(old_f) else np.nan,
            note=(f"donor-level fractions exist in-repo only for 3m vs 18m lung mice "
                  f"(n={int(m2lung['mouse.id'].nunique())} mice); present in {len(young_f)}/4 3m mice "
                  f"(mean fraction {young_f.mean():.4f} where present) vs {len(old_f)}/4 18m mice "
                  f"(mean {old_f.mean():.4f} where present). The 3m vs 18m contrast is exactly "
                  f"confounded with sequencing wave (10X_P vs MACA), so it cannot separate age from wave."))
    organ = pd.read_csv(RESULTS / "organ-composition-age-tests.csv")
    for ct, sens in (("club cell of bronchiole", "all_donor_ids"),
                     ("club cell of bronchiole", "exclude_pooled_ids"),
                     ("ciliated columnar cell of tracheobronchial tree", "all_donor_ids"),
                     ("ciliated columnar cell of tracheobronchial tree", "exclude_pooled_ids")):
        row = organ[(organ["organ"] == "Lung") & (organ["cell_type"] == ct)
                    & (organ["sensitivity"] == sens)]
        if row.empty:
            continue
        row = row.iloc[0]
        total_note = ("mouse club cells total only 15 cells across 16 mice" if "club" in ct
                      else "mouse ciliated cells total only 56 cells across 16 mice (26 excluding pooled)")
        add("mouse_refit", "Mouse", "Lung", cell_type=ct, model=f"published_summary_{sens}",
            n_donors=int(row["n_donors"]), estimate=float(row["effect"]), se=float(row["se"]),
            pvalue=float(row["pvalue"]), q_bh=float(row["q_bh"]),
            note=(f"cited from results/organ-composition-age-tests.csv (CLR slope per month, "
                  f"n_present={int(row['n_present'])}, total cells={int(row['total_cells'])}); "
                  f"{total_note}"))
    add("mouse_refit", "Mouse", "Lung", model="within_wave_refit_NOT_POSSIBLE",
        note=("No donor-level mouse lung counts exist in-repo for the 1m/21m/24m/30m mice "
              "(data/interim/m1-qc-cleaned.h5ad is absent from this checkout), so a "
              "within-sequencing-wave age slope cannot be fitted. The design audit above is "
              "as far as study/batch-aware mouse modelling can go with the data in hand."))

    # ------------------------------------------------------- 7. Lung_CellRef analogues
    src = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    cr = src[(src["cohort"] == "Lung_CellRef") & (src["age"] >= 20)].copy()
    cr_counts = cr.set_index(["donor_id", "cell_type"])["n_cells"].unstack(fill_value=0).sort_index()
    cr_meta = cr.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"))
    cr_meta["sex"] = (cr_meta["sex"].astype(str).str.lower() == "male").astype(float)
    add("cellref_replication", "Human", "Lung_CellRef", cell_type="club cell",
        model="absent_from_cohort", n_donors=len(cr_meta),
        note="club cell occurs in 0 Lung_CellRef donors; direct replication impossible (different annotation ontology)")
    add("cellref_replication", "Human", "Lung_CellRef", cell_type="multiciliated epithelial cell",
        model="absent_from_cohort", n_donors=len(cr_meta),
        note="multiciliated epithelial cell occurs in 0 Lung_CellRef donors; direct replication impossible")
    cr_clr = compute_clr_table(cr_counts)
    cr_total = cr_counts.sum(axis=1)
    for ct in ("ciliated cell", "secretory cell", "lung secretory cell", "goblet cell"):
        if ct not in cr_clr.columns:
            continue
        r = ols_slope(cr_clr[ct].loc[cr_meta.index], cr_meta[["age", "sex"]])
        add("cellref_replication", "Human", "Lung_CellRef", cell_type=ct,
            model="clr_age_sex", n_donors=len(cr_meta), age_mean=cr_meta["age"].mean(),
            age_sd=cr_meta["age"].std(), estimate=r["slope"], se=r["se"],
            ci_low=r["ci_low"], ci_high=r["ci_high"], pvalue=r["p"],
            note=("nearest available analogue of the focal airway cell types in this independent "
                  "lung cohort (n=17 donors, ages 20-76); a positive slope here is the OPPOSITE "
                  "sign to the Lung_Integrated focal declines"))
        X = cr_meta[["age", "sex"]].copy()
        X["log_total"] = np.log(cr_total)
        r2 = ols_slope(cr_clr[ct].loc[cr_meta.index], X)
        add("cellref_replication", "Human", "Lung_CellRef", cell_type=ct,
            model="clr_age_sex_log_total_cells", n_donors=len(cr_meta),
            age_mean=cr_meta["age"].mean(), age_sd=cr_meta["age"].std(),
            estimate=r2["slope"], se=r2["se"], ci_low=r2["ci_low"], ci_high=r2["ci_high"],
            pvalue=r2["p"],
            note=("same model as Phase 1 (results/reviewer-phase1-composition-total-adjusted.csv), "
                  "which reports ciliated +0.0425 (p 0.014), lung secretory +0.0435 (p 0.007), "
                  "secretory +0.0290 (p 0.113); small differences vs Phase 1 come from the "
                  "cell-type filter (Phase 1 dropped types present in <6 donors)"))
    bal = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    bal_types = sorted(bal[bal["cohort"] == "ImmuneAging_BAL"]["cell_type"].unique())
    add("cellref_replication", "Human", "ImmuneAging_BAL", cell_type="(epithelial cells)",
        model="absent_from_cohort", n_donors=int(bal[bal["cohort"] == "ImmuneAging_BAL"]["donor_id"].nunique()),
        note="BAL cohort contains immune cells only (no club, ciliated, secretory, or any epithelial type); no replication possible there")

    out = pd.DataFrame(ROWS)
    out.to_csv(OUT, index=False)
    print(f"Wrote {OUT} with {len(out)} rows")
    print(out.groupby("analysis").size().to_string())


if __name__ == "__main__":
    main()
