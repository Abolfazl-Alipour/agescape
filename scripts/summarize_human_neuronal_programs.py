#!/usr/bin/env python3
"""Summarize shared depth-adjusted neuronal gene associations in MSSM."""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
SUBTYPES = [
    "VIP GABAergic cortical interneuron",
    "pvalb GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron",
]


def main():
    results = pd.read_csv(RESULTS / "h3-human-vector-gene-drivers-all.csv")
    results = results[
        (results["model"] == "age_plus_sex_depth")
        & (results["cell_type"].isin(SUBTYPES))
    ].copy()
    results["sign"] = results["slope_log_cp10k_per_year"].apply(lambda value: 1 if value > 0 else -1)
    rows = []
    for gene, part in results.groupby("gene"):
        discovery = part[part["q_bh"] < 0.05]
        if len(discovery) < 2:
            continue
        directions = set(discovery["sign"])
        rows.append({
            "gene": gene,
            "n_fdr_subtypes": len(discovery),
            "consistent_direction": len(directions) == 1,
            "direction": discovery["sign"].iloc[0] if len(directions) == 1 else 0,
            "cell_types": ";".join(discovery["cell_type"]),
            "min_q": discovery["q_bh"].min(),
        })
    shared = pd.DataFrame(rows).sort_values(["n_fdr_subtypes", "min_q"], ascending=[False, True])
    shared.to_csv(RESULTS / "h3-human-mssm-shared-neuronal-genes.csv", index=False)
    detail = results[results["gene"].isin(shared["gene"])].copy()
    detail.to_csv(RESULTS / "h3-human-mssm-shared-neuronal-gene-details.csv", index=False)
    print(shared.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
