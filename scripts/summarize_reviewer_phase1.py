#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    pairs = pd.read_csv(RESULTS / "reviewer-phase1-fly-hit-gene-overlap.csv")
    summary = pd.DataFrame([{
        "n_vector_hit_units": 38, "n_pairs": len(pairs),
        "n_units_with_gene_hits": int((pd.concat([
            pairs[["unit_1", "n_1"]].rename(columns={"unit_1": "unit", "n_1": "n"}),
            pairs[["unit_2", "n_2"]].rename(columns={"unit_2": "unit", "n_2": "n"}),
        ]).drop_duplicates()["n"] > 0).sum()),
        "n_pairs_with_overlap": int((pairs["overlap"] > 0).sum()),
        "n_pairs_q_below_05": int((pairs["q_bh"] < 0.05).sum()),
        "median_jaccard": pairs["jaccard"].median(), "max_gene_overlap": pairs["overlap"].max(),
    }])
    summary.to_csv(RESULTS / "reviewer-phase1-fly-hit-overlap-summary.csv", index=False)
    rows = []
    for analysis, filename in [("interneuron", "reviewer-phase1-interneuron-go-enrichment.csv"),
                               ("fly", "reviewer-phase1-fly-go-enrichment.csv")]:
        result = pd.read_csv(RESULTS / filename)
        rows.append({"analysis": analysis, "n_tested_terms_with_overlap": len(result),
                     "n_terms_q_below_05": int((result["q_bh"] < 0.05).sum()),
                     "min_q": result["q_bh"].min(),
                     "top_term": result.sort_values("q_bh").iloc[0]["term"],
                     "top_list": result.sort_values("q_bh").iloc[0]["list_name"]})
    pd.DataFrame(rows).to_csv(RESULTS / "reviewer-phase1-go-summary.csv", index=False)
    signature = pd.read_csv(RESULTS / "reviewer-phase1-signature-random-null.csv")
    pd.DataFrame([{
        "n_cell_types": len(signature), "n_q_below_05": int((signature["q_bh"] < 0.05).sum()),
        "min_q": signature["q_bh"].min(), "n_random_per_cell_type": signature["n_random_signatures"].min(),
        "coverage_min": signature["p1_coverage"].min(), "coverage_max": signature["p1_coverage"].max(),
        "mean_actual_depth_r": signature["actual_score_depth_r"].mean(),
        "mean_random_depth_r": signature["random_score_depth_r_mean"].mean(),
    }]).to_csv(RESULTS / "reviewer-phase1-signature-null-summary.csv", index=False)
    composition = pd.read_csv(RESULTS / "reviewer-phase1-composition-total-adjusted.csv")
    primary = composition[composition["pseudocount"] == 0.5]
    rows = []
    for species, part in primary.groupby("species"):
        rows.append({"species": species, "n_tests": len(part),
                     "n_within_stratum_fdr": int((part["q_bh"] < 0.05).sum()),
                     "n_global_fdr": int((part["q_global_bh"] < 0.05).sum()),
                     "n_strata_with_fdr": part.groupby(["organ", "cohort"])["q_bh"].apply(lambda x: (x < 0.05).any()).sum(),
                     "n_strata_with_global_fdr": part.groupby(["organ", "cohort"])["q_global_bh"].apply(lambda x: (x < 0.05).any()).sum()})
    pd.DataFrame(rows).to_csv(RESULTS / "reviewer-phase1-composition-total-summary.csv", index=False)


if __name__ == "__main__":
    main()
