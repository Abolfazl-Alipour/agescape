#!/usr/bin/env python3
"""
scripts/run_d3_human_composition_sensitivity.py

D3 Sensitivity Analysis: Human lung compositional remodeling across study-of-origin batches.
Evaluates the 4 FDR-significant adult cell types (club, vein endothelial, AT2, mast cell)
in CELLxGENE Lung_Integrated (Sikkema et al. 2023, N=117 adult donors >= 20y).

Tests:
1. Study-covariate adjusted OLS: CLR ~ Age + Sex + Study
2. Leave-One-Study-Out (LOSO) regression: re-fitting Age + Sex dropping each study center

Outputs:
  - results/d3-human-composition-study-sensitivity.csv
  - results/d3-human-composition-loso.csv
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent

def compute_clr_table(count_matrix: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    props = (count_matrix + pseudocount).div((count_matrix + pseudocount).sum(axis=1), axis=0)
    log_props = np.log(props)
    geom_mean_log = log_props.mean(axis=1)
    clr = log_props.sub(geom_mean_log, axis=0)
    return clr

def extract_study(donor_id: str) -> str:
    prefix = "homosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747"
    raw = donor_id.replace(prefix, "")
    m = re.match(r"^([A-Za-z]+|\d+-\d+-[A-Za-z]+)", raw)
    return m.group(1) if m else "other"

def main():
    print("[D3 Composition Sensitivity] Loading results/m9-donor-scores.csv ...")
    m9_df = pd.read_csv(REPO / "results" / "m9-donor-scores.csv")
    human_lung = m9_df[m9_df["cohort"] == "Lung_Integrated"].copy()

    # Adult donors (age >= 20.0, N=117)
    adult_lung = human_lung[human_lung["age"] >= 20.0].copy()
    h_ct_counts = adult_lung["cell_type"].value_counts()
    valid_cts = h_ct_counts[h_ct_counts >= 20].index.tolist()

    counts = (
        adult_lung[adult_lung["cell_type"].isin(valid_cts)]
        .set_index(["donor_id", "cell_type"])["n_cells"]
        .unstack(fill_value=0)
    )

    meta = (
        adult_lung.groupby("donor_id")
        .agg({"age": "first", "sex": "first"})
    )
    meta["sex_bin"] = (meta["sex"].astype(str).str.lower() == "male").astype(float)
    meta["study"] = [extract_study(d) for d in meta.index]

    clr = compute_clr_table(counts)

    target_cts = [
        "club cell",
        "vein endothelial cell",
        "pulmonary alveolar type 2 cell",
        "mast cell"
    ]

    # Map to major studies (studies with at least 5 donors; others to 'other')
    study_counts = meta["study"].value_counts()
    major_studies = study_counts[study_counts >= 5].index.tolist()
    meta["study_grouped"] = meta["study"].apply(lambda s: s if s in major_studies else "other")

    study_records = []
    loso_records = []

    for ct in target_cts:
        if ct not in clr.columns:
            continue
        y = clr[ct].loc[meta.index]
        
        # 1. Baseline model: Age + Sex
        X_base = sm.add_constant(meta[["age", "sex_bin"]])
        res_base = sm.OLS(y, X_base).fit()

        # 2. Study-adjusted model: Age + Sex + C(study_grouped)
        dummies = pd.get_dummies(meta["study_grouped"], drop_first=True, dtype=float)
        X_adj = sm.add_constant(pd.concat([meta[["age", "sex_bin"]], dummies], axis=1))
        res_adj = sm.OLS(y, X_adj).fit()

        study_records.append({
            "cell_type": ct,
            "baseline_slope": res_base.params["age"],
            "baseline_se": res_base.bse["age"],
            "baseline_p": res_base.pvalues["age"],
            "study_adj_slope": res_adj.params["age"],
            "study_adj_se": res_adj.bse["age"],
            "study_adj_p": res_adj.pvalues["age"]
        })

        # 3. Leave-one-study-out
        for s in major_studies:
            s_mask = meta["study"] != s
            if s_mask.sum() >= 20:
                X_loso = sm.add_constant(meta.loc[s_mask, ["age", "sex_bin"]])
                y_loso = y.loc[s_mask]
                res_l = sm.OLS(y_loso, X_loso).fit()
                loso_records.append({
                    "cell_type": ct,
                    "dropped_study": s,
                    "n_remaining_donors": int(s_mask.sum()),
                    "slope": res_l.params["age"],
                    "se": res_l.bse["age"],
                    "p": res_l.pvalues["age"]
                })

    res_df = pd.DataFrame(study_records)
    res_df.to_csv(REPO / "results" / "d3-human-composition-study-sensitivity.csv", index=False)

    loso_df = pd.DataFrame(loso_records)
    loso_df.to_csv(REPO / "results" / "d3-human-composition-loso.csv", index=False)

    print("Saved D3 study-adjusted results to results/d3-human-composition-study-sensitivity.csv:")
    print(res_df.to_string())
    print("\nSummary of LOSO slope ranges:")
    for ct in target_cts:
        sub = loso_df[loso_df["cell_type"] == ct]
        print(f"  {ct}: slope range [{sub['slope'].min():+.4f}, {sub['slope'].max():+.4f}], p range [{sub['p'].min():.4f}, {sub['p'].max():.4f}]")

if __name__ == "__main__":
    main()
