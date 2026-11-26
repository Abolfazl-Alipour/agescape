#!/usr/bin/env python3
"""Join the fly/worm check outputs into one per-row summary table.

Reads the scan CSV (26 passing fly/worm rows), the propeller CSV, the LODO
CSV, and the analogs CSV. Adds a pooled-across-the-26-rows BH q for the
propeller p-values (written back into the propeller CSV as q_bh_pooled_26,
NA for rows outside the 26) and writes
results/followup-mature-fly-worm-checks-summary.csv, one row per passing row.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

SCAN = REPO / "results" / "mature-all-species-composition.csv"
PROP = REPO / "results" / "followup-mature-fly-worm-propeller.csv"
LODO = REPO / "results" / "followup-mature-fly-worm-lodo.csv"
ANALOGS = REPO / "results" / "followup-mature-fly-worm-analogs.csv"
OUTPUT = REPO / "results" / "followup-mature-fly-worm-checks-summary.csv"

KEY = ["species", "organ", "cohort", "cell_type"]


def main():
    scan = pd.read_csv(SCAN)
    passing = scan[(scan["significant_pooled_05"]) & (scan["species"].isin(["Fruit fly", "C. elegans"]))].copy()
    assert len(passing) == 26, f"expected 26 passing rows, got {len(passing)}"

    prop = pd.read_csv(PROP)
    merged = passing.merge(prop, on=KEY, how="left", suffixes=("", "_prop"))
    assert merged["pvalue_prop"].notna().all(), "propeller rows missing for passing rows"
    merged = merged.rename(columns={"effect": "scan_effect", "q_pooled_mature_bh": "scan_q_pooled",
                                  "q_bh": "scan_q_within_stratum"})
    pooled_q = multipletests(merged["pvalue_prop"], method="fdr_bh")[1]
    merged["propeller_q_pooled_26"] = pooled_q
    merged["propeller_pass_pooled_26"] = pooled_q < 0.05
    # write the pooled q back into the propeller CSV (NA outside the 26)
    prop = prop.drop(columns=["propeller_q_pooled_26"], errors="ignore")
    prop = prop.merge(merged[KEY + ["propeller_q_pooled_26"]], on=KEY, how="left")
    prop.to_csv(PROP, index=False)

    lodo = pd.read_csv(LODO)
    merged = merged.merge(lodo, on=KEY + ["scan_effect", "scan_q_pooled"], how="left")
    assert merged["lodo_slope_min"].notna().all(), "LODO rows missing for passing rows"

    analogs = pd.read_csv(ANALOGS)
    sex = analogs[analogs["check"] == "sex_analog"].copy()
    for label, prefix in (("male_only", "male"), ("female_only", "female")):
        sub = sex[sex["subset"] == label][KEY + ["effect", "pvalue", "sign_matches_scan"]]
        sub = sub.rename(columns={"effect": f"{prefix}_effect", "pvalue": f"{prefix}_pvalue",
                                  "sign_matches_scan": f"{prefix}_sign_matches"})
        merged = merged.merge(sub, on=KEY, how="left")
    plat = analogs[analogs["check"] == "platform_analog"].copy()
    for label, prefix in (("platform_adjusted", "platform_adjusted"),
                          ("platform_0", "platform_0"), ("platform_1", "platform_1")):
        sub = plat[plat["subset"] == label][KEY + ["effect", "pvalue", "sign_matches_scan"]]
        sub = sub.rename(columns={"effect": f"{prefix}_effect", "pvalue": f"{prefix}_pvalue",
                                  "sign_matches_scan": f"{prefix}_sign_matches"})
        merged = merged.merge(sub, on=KEY, how="left")

    merged["sex_sign_holds_both"] = (merged["male_sign_matches"] == True) & (merged["female_sign_matches"] == True)
    cols = KEY + [
        "scan_effect", "scan_q_pooled",
        "propeller_effect_logit", "propeller_pvalue", "propeller_q_within_cohort",
        "propeller_q_pooled_26", "propeller_pass_pooled_26",
        "lodo_slope_min", "lodo_slope_max", "lodo_same_sign",
        "lodo_max_abs_slope_change", "lodo_max_change_pool", "lodo_pvalue_max",
        "male_effect", "male_pvalue", "male_sign_matches",
        "female_effect", "female_pvalue", "female_sign_matches", "sex_sign_holds_both",
        "platform_adjusted_effect", "platform_adjusted_pvalue", "platform_adjusted_sign_matches",
        "platform_0_effect", "platform_0_pvalue", "platform_0_sign_matches",
        "platform_1_effect", "platform_1_pvalue", "platform_1_sign_matches",
    ]
    out = merged.rename(columns={"effect_logit": "propeller_effect_logit",
                                 "pvalue_prop": "propeller_pvalue",
                                 "q_bh_prop": "propeller_q_within_cohort"})
    out = out[[c for c in cols if c in out.columns]]
    out.to_csv(OUTPUT, index=False)
    print(f"wrote {len(out)} rows to {OUTPUT}")
    print(f"propeller pooled q < 0.05: {int(out['propeller_pass_pooled_26'].sum())} of 26")
    print(f"LODO sign holds: {int(out['lodo_same_sign'].sum())} of 26")
    fly = out[out["species"] == "Fruit fly"]
    print(f"sex sign holds in both sexes: {int(fly['sex_sign_holds_both'].sum())} of {len(fly)} fly rows")
    print(out[["species", "organ", "cell_type", "scan_effect", "propeller_q_pooled_26",
               "lodo_same_sign", "sex_sign_holds_both"]].to_string(index=False))


if __name__ == "__main__":
    main()
