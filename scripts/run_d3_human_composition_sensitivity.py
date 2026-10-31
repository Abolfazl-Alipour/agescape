#!/usr/bin/env python3
"""Assess study-of-origin sensitivity for every tested human lung cell type."""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.run_h2_composition import compute_clr_table, load_human_counts


def extract_study(donor_id: str) -> str:
    prefix = "homosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747"
    raw = str(donor_id).replace(prefix, "")
    match = re.match(r"^([A-Za-z]+|\d+-\d+-[A-Za-z]+)", raw)
    return match.group(1) if match else "other"


def main():
    counts, meta = load_human_counts(REPO / "results")
    meta["study"] = [extract_study(donor) for donor in meta.index]
    clr = compute_clr_table(counts)
    study_counts = meta["study"].value_counts()
    major_studies = study_counts[study_counts >= 5].index.tolist()
    meta["study_grouped"] = meta["study"].apply(lambda s: s if s in major_studies else "other")

    study_records = []
    loso_records = []
    for cell_type in clr.columns:
        y = clr[cell_type].loc[meta.index]
        base = sm.OLS(y, sm.add_constant(meta[["age", "sex"]])).fit()
        dummies = pd.get_dummies(meta["study_grouped"], drop_first=True, dtype=float)
        adjusted = sm.OLS(y, sm.add_constant(pd.concat([meta[["age", "sex"]], dummies], axis=1))).fit()
        study_records.append({
            "cell_type": cell_type,
            "baseline_slope": base.params["age"], "baseline_se": base.bse["age"],
            "baseline_p": base.pvalues["age"], "study_adj_slope": adjusted.params["age"],
            "study_adj_se": adjusted.bse["age"], "study_adj_p": adjusted.pvalues["age"],
            "n_donors": len(meta),
        })
        for study in major_studies:
            keep = meta["study"] != study
            if keep.sum() < 20:
                continue
            fit = sm.OLS(y.loc[keep], sm.add_constant(meta.loc[keep, ["age", "sex"]])).fit()
            loso_records.append({
                "cell_type": cell_type, "dropped_study": study,
                "n_remaining_donors": int(keep.sum()), "slope": fit.params["age"],
                "se": fit.bse["age"], "p": fit.pvalues["age"],
            })

    study_df = pd.DataFrame(study_records)
    loso_df = pd.DataFrame(loso_records)
    study_df.to_csv(REPO / "results" / "d3-human-composition-study-sensitivity.csv", index=False)
    loso_df.to_csv(REPO / "results" / "d3-human-composition-loso.csv", index=False)
    print(study_df.sort_values("study_adj_p").head(10).to_string(index=False))


if __name__ == "__main__":
    main()
