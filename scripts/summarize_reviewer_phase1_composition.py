#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    data = pd.read_csv(RESULTS / "reviewer-phase1-composition-total-adjusted.csv")
    rows = []
    for species, part in data.groupby("species"):
        wide = part.pivot(index=["organ", "cohort", "cell_type"], columns="pseudocount", values="effect")
        qwide = part.pivot(index=["organ", "cohort", "cell_type"], columns="pseudocount", values="q_bh")
        rows.append({
            "species": species, "n_units": len(wide),
            "fraction_same_sign_all_pseudocounts": ((wide.min(axis=1) > 0) | (wide.max(axis=1) < 0)).mean(),
            "fraction_fdr_all_pseudocounts": (qwide < 0.05).all(axis=1).mean(),
            "fraction_fdr_any_pseudocount": (qwide < 0.05).any(axis=1).mean(),
            "effect_pearson_0.1_vs_1.0": pearsonr(wide[0.1], wide[1.0]).statistic,
            "effect_spearman_0.1_vs_1.0": spearmanr(wide[0.1], wide[1.0]).statistic,
        })
    pd.DataFrame(rows).to_csv(RESULTS / "reviewer-phase1-composition-pseudocount-summary.csv", index=False)
    propeller = pd.read_csv(RESULTS / "reviewer-phase1-composition-propeller.csv")
    primary = data[data["pseudocount"] == 0.5]
    merged = primary.merge(propeller, on=["species", "organ", "cohort", "cell_type"], suffixes=("_clr", "_propeller"))
    rows = []
    for species, part in merged.groupby("species"):
        rows.append({
            "species": species, "n_units": len(part),
            "sign_concordance": np.mean(part["effect"] * part["effect_logit"] > 0),
            "effect_pearson": pearsonr(part["effect"], part["effect_logit"]).statistic,
            "effect_spearman": spearmanr(part["effect"], part["effect_logit"]).statistic,
            "n_clr_fdr": int((part["q_bh_clr"] < 0.05).sum()),
            "n_propeller_fdr": int((part["q_bh_propeller"] < 0.05).sum()),
            "n_both_fdr": int(((part["q_bh_clr"] < 0.05) & (part["q_bh_propeller"] < 0.05)).sum()),
        })
    pd.DataFrame(rows).to_csv(RESULTS / "reviewer-phase1-composition-propeller-summary.csv", index=False)


if __name__ == "__main__":
    main()
