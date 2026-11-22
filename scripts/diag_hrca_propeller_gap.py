#!/usr/bin/env python3
"""Diagnostic: why does HRCA propeller disagree with the batch-B CLR result?

Refits the AC fraction slope in Sn_Adult_Unenriched under different
transforms and covariate sets, all OLS, to locate the source of the gap.
No files are written.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def ols(y, X, names):
    Xd = np.column_stack([np.ones(len(y))] + [np.asarray(v, dtype=float) for v in X])
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    dof = len(y) - Xd.shape[1]
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(Xd.T @ Xd)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return {n: (float(b), float(s), float(pp)) for n, b, s, pp in zip(["const"] + names, beta, se, p)}


def main():
    from scripts.run_followup_mature_organ_transports_b_hrca import SUITES, build_donor_lists
    from scripts.run_followup_organ_10_eye_hrca import CLASSES, SC_PATH, SN_PATH, _aggregate, _donor_meta

    sn = _aggregate(SN_PATH, "sn")
    sc = _aggregate(SC_PATH, "sc")
    sck = sc["keys"]
    sck = sck[~((sck["study"] == "Sanes_GSE148077") & (sck["regiongroup"] != "macular"))].copy()
    donor_lists = build_donor_lists(sn["keys"], sck)
    suite = [s for s in SUITES if s["name"] == "Sn_Adult_Unenriched"][0]
    keys = sn["keys"]
    donors = donor_lists[suite["donors_key"]]
    k = keys[keys["donor_id"].isin(donors) & keys["regiongroup"].isin(suite["regions"])].copy()
    meta = _donor_meta(k, donors, str(suite["studies_ref"]))
    cc = k.pivot_table(index="donor_id", columns="majorclass", values="n_cells",
                       aggfunc="sum", fill_value=0)
    cc = cc.reindex(index=meta.index, columns=[c for c in CLASSES if c in cc.columns], fill_value=0)
    counts = cc.values.astype(float)
    shares = (counts + 0.5) / (counts + 0.5).sum(axis=1, keepdims=True)
    clr = np.log(shares) - np.log(shares).mean(axis=1, keepdims=True)
    j = list(cc.columns).index("AC")
    y_clr = clr[:, j]
    p_ac = shares[:, j]
    y_logit = np.log(p_ac / (1 - p_ac))
    age_dec = meta["age"].to_numpy(dtype=float) / 10.0
    age_yr = meta["age"].to_numpy(dtype=float)
    sex = (meta["sex"] == "male").astype(float).to_numpy()
    log_total = np.log(counts.sum(axis=1))
    log_lib = np.log(counts.sum(axis=1))
    study_cols = [c for c in meta.columns if c.startswith("study_")]
    frac_per = meta["frac_peripheral"].to_numpy(dtype=float)

    print("n donors:", len(meta))
    print("CLR(AC) ~ age_decade : slope/decade, se, p =",
          ols(y_clr, [age_dec], ["age_dec"])["age_dec"])
    print("CLR(AC) ~ age_decade + sex + log_total :",
          ols(y_clr, [age_dec, sex, log_total], ["age_dec", "sex", "log_total"])["age_dec"])
    print("logit(AC) ~ age_years + sex + log_total : per-year slope =",
          ols(y_logit, [age_yr, sex, log_total], ["age_yr", "sex", "log_total"])["age_yr"])
    full = ols(y_clr, [age_dec] + [meta[c].to_numpy(dtype=float) for c in study_cols]
               + [sex, frac_per, log_lib],
               ["age_dec"] + study_cols + ["sex", "frac_per", "log_lib"])
    print("CLR(AC) ~ age_decade + study + sex + frac_per + log_lib (batch-B-like) :",
          full["age_dec"])
    nostudy = ols(y_clr, [age_dec, sex, frac_per, log_lib],
                  ["age_dec", "sex", "frac_per", "log_lib"])
    print("CLR(AC) ~ age_decade + sex + frac_per + log_lib (no study) :",
          nostudy["age_dec"])
    studonly = ols(y_clr, [age_dec] + [meta[c].to_numpy(dtype=float) for c in study_cols] + [sex],
                   ["age_dec"] + study_cols + ["sex"])
    print("CLR(AC) ~ age_decade + study + sex :",
          studonly["age_dec"])
    print("study columns:", study_cols)


if __name__ == "__main__":
    main()
