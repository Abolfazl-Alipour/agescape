#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

RESULTS = Path(__file__).resolve().parent.parent / "results"
MOUSE_AGES = {"1m": 1, "3m": 3, "18m": 18, "21m": 21, "24m": 24, "30m": 30}


def main():
    excluded = pd.read_csv(RESULTS / "reviewer-phase2-annotation-exclusion-by-age.csv")
    excluded["age_numeric"] = pd.to_numeric(excluded["age"].map(MOUSE_AGES).fillna(excluded["age"]))
    rows = []
    for (species, cohort, tissue), part in excluded.groupby(["species", "cohort", "tissue"]):
        for metric in ("fraction_excluded_qc", "fraction_unassigned_raw"):
            if len(part) < 4 or part[metric].nunique() < 2:
                continue
            result = spearmanr(part["age_numeric"], part[metric])
            rows.append({"species": species, "cohort": cohort, "tissue": tissue, "metric": metric,
                         "n_age_groups": len(part), "spearman_rho": result.statistic, "pvalue": result.pvalue})
    tests = pd.DataFrame(rows)
    if not tests.empty:
        tests["q_bh"] = multipletests(tests["pvalue"], method="fdr_bh")[1]
    tests.to_csv(RESULTS / "reviewer-phase2-annotation-exclusion-age-tests.csv", index=False)
    markers = pd.read_csv(RESULTS / "reviewer-phase2-marker-stability.csv")
    summary = markers.groupby("species").agg(n_tests=("pvalue", "size"),
        n_fdr=("q_bh", lambda x: int((x < 0.05).sum())), min_q=("q_bh", "min")).reset_index()
    summary.to_csv(RESULTS / "reviewer-phase2-marker-stability-summary.csv", index=False)


if __name__ == "__main__":
    main()
