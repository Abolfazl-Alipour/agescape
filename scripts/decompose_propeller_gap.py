#!/usr/bin/env python3
"""Decompose scan-vs-propeller discrepancies for extreme entered rows.

For each target row, refit plain OLS under four models on the SAME
propeller inputs (data/interim/mature-propeller-fly-worm) and compare the
age slope and p-value:
  (a) CLR ~ age + sex            (scan design)
  (b) logit ~ age + sex          (transform change only)
  (c) logit ~ age + sex + log_total   (propeller covariates, no eBayes)
  (d) logit ~ age + sex + log_total + platform  (worm only)
The committed propeller p-value (limma eBayes on model c/d) is printed
alongside for reference. If (b) already shows the signal, the transform
(CLR centering vs own-share logit) explains it; if only (c)/(d) does, the
covariates explain it; if none does, limma variance shrinkage does.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
ROOT = REPO / "data" / "interim" / "mature-propeller-fly-worm"

TARGETS = [
    ("Fruit fly", "body", "AFCA", "gustatory receptor neuron"),
    ("Fruit fly", "body", "AFCA", "intestinal stem cell"),
    ("Fruit fly", "body", "AFCA", "adult reticular neuropil associated glial cell_body"),
    ("Fruit fly", "body", "AFCA", "adult midgut-hindgut hybrid zone"),
    ("C. elegans", "Whole animal", "Gao_GSE229022", "Coelomocyte"),
    ("C. elegans", "Whole animal", "Gao_GSE229022", "Neuron"),
    ("C. elegans", "Whole animal", "Gao_GSE229022", "ARSC"),
]


def ols_age_p(y, X, names):
    Xd = np.column_stack([np.ones(len(y))] + [np.asarray(v, dtype=float) for v in X])
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    dof = len(y) - Xd.shape[1]
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(Xd.T @ Xd)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    i = names.index("age") + 1
    return float(beta[i]), float(se[i]), float(p[i])


def main():
    manifest = pd.read_csv(ROOT / "manifest.csv")
    full = pd.read_csv(REPO / "results" / "followup-mature-propeller-full-577.csv")
    full_index = {(r.species, r.organ, r.cohort, r.cell_type): r for r in full.itertuples()}

    for species, organ, cohort, cell_type in TARGETS:
        m = manifest[(manifest["species"] == species) & (manifest["organ"] == organ)
                     & (manifest["cohort"] == cohort)]
        assert len(m) == 1, f"manifest match {species} {organ} {cohort}: {len(m)}"
        key = m.iloc[0]["key"]
        counts = pd.read_csv(ROOT / key / "counts.csv", index_col=0)
        meta = pd.read_csv(ROOT / key / "metadata.csv", index_col=0)
        counts = counts.loc[meta.index]
        assert cell_type in counts.columns, f"{cell_type} not in counts"
        totals = counts.sum(axis=1).to_numpy(dtype=float)
        props = ((counts + 0.5) / (totals[:, None] + 0.5 * counts.shape[1])).to_numpy(dtype=float)
        j = list(counts.columns).index(cell_type)
        shares = props[:, j]
        clr_all = np.log(props) - np.log(props).mean(axis=1, keepdims=True)
        y_clr = clr_all[:, j]
        y_logit = np.log(shares / (1 - shares))
        age = meta["age"].to_numpy(dtype=float)
        sex = meta["sex_male"].to_numpy(dtype=float) if "sex_male" in meta else None
        log_total = np.log(totals)
        has_sex = sex is not None and len(np.unique(sex)) > 1
        plat = meta["platform_nextseq"].to_numpy(dtype=float) if "platform_nextseq" in meta else None
        has_plat = plat is not None and len(np.unique(plat)) > 1

        def fit(y, extra):
            cols, names = [age], ["age"]
            if has_sex:
                cols.append(sex)
                names.append("sex")
            for v, n in extra:
                cols.append(v)
                names.append(n)
            return ols_age_p(y, cols, names)

        a = fit(y_clr, [])
        b = fit(y_logit, [])
        c = fit(y_logit, [(log_total, "log_total")])
        d = fit(y_logit, [(log_total, "log_total"), (plat, "platform")]) if has_plat else None
        ref = full_index[(species, organ, cohort, cell_type)]
        print(f"{cell_type} [{species} {organ}]", flush=True)
        print(f"  (a) CLR age+sex: slope {a[0]:+.4f} p {a[2]:.3g}", flush=True)
        print(f"  (b) logit age+sex: slope {b[0]:+.4f} p {b[2]:.3g}", flush=True)
        print(f"  (c) logit age+sex+log_total: slope {c[0]:+.4f} p {c[2]:.3g}", flush=True)
        if d is not None:
            print(f"  (d) logit +platform: slope {d[0]:+.4f} p {d[2]:.3g}", flush=True)
        print(f"  committed propeller (limma eBayes): effect {float(ref.propeller_effect_logit):+.4f} "
              f"p {float(ref.propeller_pvalue):.3g} q577 {float(ref.q_propeller_pooled_577):.3g} | "
              f"scan CLR effect {float(ref.scan_effect):+.4f} p {float(ref.scan_pvalue):.3g}",
              flush=True)


if __name__ == "__main__":
    main()
