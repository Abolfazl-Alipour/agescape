#!/usr/bin/env python3
"""Stage 1 validation: how many of the 28 name-matched mouse-human pairs
(results/followup-mature-matched-p1m.csv) does SAMap recover?

A pair counts as recovered only if the SAMap pair for exactly that
(mouse tissue + cell type, human cohort + cell type) passes the full
pre-stated match rule (alignment >= 0.30, reciprocal best hit,
>= 40 enriched gene pairs, top 1% of pairs in the comparison).

Gate (pre-stated): recover >= 14 of 28 (half). Below that, STOP and report;
do not proceed to distant arms.

For every pair this also reports the per-rule diagnostics (alignment score,
reciprocal-best flag, gene-pair count, top-1% flag) so a gate failure can be
attributed to a specific sub-rule instead of guessed at.

Usage: python scripts/validate_samap_hums.py
Writes results/samap-hums/validation.csv and prints the recovery count.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
GATE = 14

TISSUE_OF_COHORT = {
    "Kidney_LivingDonor": "kidney",
    "Kidney_SCR_v15": "kidney",
    "Lung_CellRef": "lung",
    "Lung_Integrated": "lung",
    "ImmuneAging_BAL": "lung",
    "Yazar_Blood": "blood",
}


def main():
    matches = pd.read_csv(RESULTS / "samap-hums" / "matches.csv")
    matches["ms_cluster"] = matches.apply(
        lambda r: r["cluster1"] if r["sid1"] == "ms" else r["cluster2"], axis=1)
    matches["hu_cluster"] = matches.apply(
        lambda r: r["cluster1"] if r["sid1"] == "hu" else r["cluster2"], axis=1)
    by_pair = matches.set_index(["ms_cluster", "hu_cluster"])

    val = pd.read_csv(RESULTS / "followup-mature-matched-p1m.csv")
    pairs = val[["mouse_tissue", "mouse_unit", "human_cohort",
                 "human_unit"]].drop_duplicates()
    assert len(pairs) == 28, f"expected 28 name-matched pairs, got {len(pairs)}"

    rows = []
    for r in pairs.itertuples():
        cell = r.mouse_unit.split("::", 1)[1]
        ms_cl = f"{r.mouse_tissue} | TMS | {cell}"
        hu_cl = (f"{TISSUE_OF_COHORT[r.human_cohort]} | "
                 f"{r.human_cohort} | {r.human_unit}")
        try:
            row = by_pair.loc[(ms_cl, hu_cl)]
            rows.append({
                "mouse_tissue": r.mouse_tissue, "mouse_unit": cell,
                "human_cohort": r.human_cohort, "human_unit": r.human_unit,
                "in_matches_table": True,
                "score": float(row["score"]),
                "reciprocal_best": bool(row["reciprocal_best"]),
                "n_gene_pairs": int(row["n_gene_pairs"]),
                "top1pct": bool(row["top1pct"]),
                "passes": bool(row["passes"]),
            })
        except KeyError:
            rows.append({
                "mouse_tissue": r.mouse_tissue, "mouse_unit": cell,
                "human_cohort": r.human_cohort, "human_unit": r.human_unit,
                "in_matches_table": False,
                "score": None, "reciprocal_best": None,
                "n_gene_pairs": None, "top1pct": None, "passes": False,
            })
    out = pd.DataFrame(rows)
    out_path = RESULTS / "samap-hums" / "validation.csv"
    out.to_csv(out_path, index=False)
    n_recovered = int(out["passes"].sum())
    print(f"Recovered {n_recovered} of 28 name-matched pairs (gate: >= {GATE})")
    print(f"GATE: {'PASS' if n_recovered >= GATE else 'FAIL'}")
    print(out.to_string(index=False))
    # failure attribution: which sub-rule fails most often (among table rows)
    tab = out[out["in_matches_table"]]
    print("\nAmong the pairs present in the matches table:")
    for rule in ("score", "reciprocal_best", "n_gene_pairs", "top1pct"):
        if rule == "score":
            print(f"  score >= 0.30: {int((tab['score'] >= 0.30).sum())} of {len(tab)}")
        elif rule == "n_gene_pairs":
            print(f"  gene pairs >= 40: {int((tab['n_gene_pairs'] >= 40).sum())} of {len(tab)}")
        else:
            print(f"  {rule}: {int(tab[rule].sum())} of {len(tab)}")


if __name__ == "__main__":
    main()
