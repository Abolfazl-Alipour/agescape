#!/usr/bin/env python3
"""Repeat the D3 and Rev2 study-aware lung composition arms at age 30-plus."""

import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_h2_composition import compute_clr_table, load_human_counts
from scripts.run_rev2_composition import dl_meta, study_of
from scripts.run_mature_composition_complete import MODEL_USED, RUNTIME_FILE, utc_now, write_report

RESULTS = REPO / "results"
OUTPUT = RESULTS / "followup-mature-composition-lung-study.csv"
CUTOFF = 30.0


def ols_slope(y, predictors):
    X = sm.add_constant(predictors.astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or len(y) <= X.shape[1] + 1:
        return None
    fit = sm.OLS(y.astype(float), X).fit()
    ci = fit.conf_int().loc["age"]
    return {
        "effect": float(fit.params["age"]), "se": float(fit.bse["age"]),
        "ci95_low": float(ci.iloc[0]), "ci95_high": float(ci.iloc[1]),
        "pvalue": float(fit.pvalues["age"]), "n_donors": int(len(y)),
    }


def record(rows, analysis, cell_type, model, group="", result=None, status="TESTED",
           note="", n_donors=np.nan, n_studies=np.nan, metadata=None):
    row = {
        "analysis": analysis, "species": "Human", "organ": "lung",
        "cohort": "Lung_Integrated", "cell_type": cell_type, "group": group,
        "model": model, "status": status, "n_donors": n_donors, "n_studies": n_studies,
        "age_min": np.nan, "age_max": np.nan, "age_mean": np.nan, "age_sd": np.nan,
        "effect": np.nan, "se": np.nan, "ci95_low": np.nan, "ci95_high": np.nan,
        "pvalue": np.nan, "q_bh": np.nan, "primary_pooled_q_431_family": np.nan,
        "primary_pooled_q_refit_family": np.nan, "note": note,
    }
    if metadata is not None and len(metadata):
        row.update({
            "n_donors": len(metadata), "age_min": metadata["age"].min(),
            "age_max": metadata["age"].max(), "age_mean": metadata["age"].mean(),
            "age_sd": metadata["age"].std(),
        })
    if result is not None:
        row.update(result)
    rows.append(row)


def append_runtime(pipeline, started):
    ended = utc_now()
    duration = (datetime.fromisoformat(ended) - datetime.fromisoformat(started)).total_seconds()
    row = pd.DataFrame([{
        "pipeline": pipeline, "start_utc": started, "end_utc": ended,
        "duration_seconds": duration, "model_used": MODEL_USED, "status": "completed",
    }])
    if RUNTIME_FILE.exists():
        previous = pd.read_csv(RUNTIME_FILE)
        previous = previous[previous["pipeline"] != pipeline]
        row = pd.concat([previous, row], ignore_index=True)
    row.to_csv(RUNTIME_FILE, index=False)


def main():
    started = utc_now()
    print(f"START {started} MODEL {MODEL_USED}", flush=True)
    counts, meta = load_human_counts(RESULTS)
    keep = meta["age"] >= CUTOFF
    meta = meta.loc[keep].copy()
    counts = counts.loc[meta.index].copy()
    if meta.empty or meta["age"].min() < CUTOFF:
        raise RuntimeError("Mature lung donor restriction failed")
    meta["study"] = [study_of(donor) for donor in meta.index]
    clr = compute_clr_table(counts)
    study_counts = meta["study"].value_counts()
    major_studies = study_counts[study_counts >= 5].index.tolist()
    meta["study_grouped"] = meta["study"].apply(lambda s: s if s in major_studies else "other")
    primary = pd.read_csv(RESULTS / "followup-mature-composition-human.csv")
    primary = primary[primary["cohort"] == "Lung_Integrated"].set_index("cell_type")
    rows = []
    for cell_type in clr.columns:
        y = clr[cell_type].loc[meta.index]
        result = ols_slope(y, meta[["age", "sex"]])
        note = ("D3 adult-eligible cell universe with mature donors only; primary q values are "
                "references from the complete mature model and use its 6-donor mature cell universe")
        record(rows, "pooled_model", cell_type, "baseline_age_sex", result=result,
               status="TESTED" if result else "UNTESTABLE", note=note, metadata=meta)
        dummies = pd.get_dummies(meta["study_grouped"], drop_first=True, dtype=float)
        adjusted_predictors = pd.concat([meta[["age", "sex"]], dummies], axis=1)
        adjusted = ols_slope(y, adjusted_predictors)
        record(rows, "pooled_model", cell_type, "study_fixed_effects", result=adjusted,
               status="TESTED" if adjusted else "UNTESTABLE", n_studies=len(major_studies) + 1,
               note="study groups with at least 5 mature donors; smaller groups pooled as other",
               metadata=meta)
        for row in rows[-2:]:
            if cell_type in primary.index:
                row["primary_pooled_q_431_family"] = primary.loc[cell_type, "q_431_family"]
                row["primary_pooled_q_refit_family"] = primary.loc[cell_type, "q_refit_family"]
        per_study = []
        for study, part in meta.groupby("study"):
            if len(part) < 5 or part["age"].std() < 3:
                record(rows, "per_study_slope", cell_type, "within_study_ols_age_sex", group=study,
                       result=None, status="UNTESTABLE", metadata=part,
                       note="not estimated: fewer than 5 donors or within-study age SD below 3 years")
                continue
            columns = ["age"] + (["sex"] if part["sex"].nunique() > 1 else [])
            result = ols_slope(y.loc[part.index], part[columns])
            model = "within_study_ols_age_sex" if "sex" in columns else "within_study_ols_age_only"
            if result is None:
                record(rows, "per_study_slope", cell_type, model, group=study,
                       status="UNTESTABLE", metadata=part, note="rank-deficient or insufficient residual degrees of freedom")
                continue
            result["study"] = study
            result["slope"] = result["effect"]
            per_study.append(result)
            record(rows, "per_study_slope", cell_type, model, group=study, result=result,
                   metadata=part, note="mature donors only")
        if per_study:
            qvalues = multipletests([r["pvalue"] for r in per_study], method="fdr_bh")[1]
            for estimate, qvalue in zip(per_study, qvalues):
                for row in rows:
                    if (row["analysis"] == "per_study_slope" and row["cell_type"] == cell_type
                            and row["group"] == estimate["study"] and row["status"] == "TESTED"):
                        row["q_bh"] = float(qvalue)
            negative = sum(r["effect"] < 0 for r in per_study)
            negative_p = sum(r["effect"] < 0 and r["pvalue"] < 0.05 for r in per_study)
            record(rows, "per_study_summary", cell_type, "per_study_count",
                   n_donors=sum(r["n_donors"] for r in per_study), n_studies=len(per_study),
                   note=(f"{negative}/{len(per_study)} estimable studies negative; "
                         f"{negative_p}/{len(per_study)} negative and p<0.05"))
            combined = dl_meta(per_study)
            if combined is not None:
                common = {"n_donors": sum(r["n_donors"] for r in per_study), "n_studies": combined["k"]}
                for model, effect_key, se_key in (("fixed_effect_ivw", "fe", "fe_se"),
                                                  ("random_effects_dl", "re", "re_se")):
                    effect = combined[effect_key]
                    se = combined[se_key]
                    result = {
                        "effect": effect, "se": se, "ci95_low": effect - 1.96 * se,
                        "ci95_high": effect + 1.96 * se,
                        "pvalue": float(2 * norm.sf(abs(effect / se))),
                    }
                    record(rows, "meta_combine", cell_type, model, result=result, **common,
                           note=(f"Cochran Q={combined['q']:.8g}; df={combined['df']}; "
                                 f"I2={combined['i2']:.8g}; tau2={combined['tau2']:.8g}"))
        for study in major_studies:
            retained = meta["study"] != study
            if retained.sum() < 20:
                record(rows, "loso", cell_type, "leave_one_study_out_age_sex", group=study,
                       status="UNTESTABLE", n_donors=int(retained.sum()),
                       note="fewer than 20 remaining mature donors")
                continue
            result = ols_slope(y.loc[retained], meta.loc[retained, ["age", "sex"]])
            record(rows, "loso", cell_type, "leave_one_study_out_age_sex", group=study,
                   result=result, status="TESTED" if result else "UNTESTABLE",
                   n_donors=int(retained.sum()), metadata=meta.loc[retained],
                   note=f"dropped study {study}; all remaining donors are at least 30 years")
    output = pd.DataFrame(rows)
    output.to_csv(OUTPUT, index=False)
    append_runtime("mature_lung_study_sensitivity", started)
    write_report()
    print(f"END {utc_now()} rows={len(output)} donors={len(meta)} major_studies={len(major_studies)}", flush=True)
    print(output[(output["cell_type"].isin(["club cell", "multiciliated epithelial cell"]))
                 & (output["analysis"].isin(["pooled_model", "meta_combine", "per_study_summary"]))]
          [["cell_type", "analysis", "model", "effect", "se", "pvalue", "n_donors", "n_studies", "note"]]
          .to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
