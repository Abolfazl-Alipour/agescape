#!/usr/bin/env python3
"""Variant 3b (supplement, not in the original brief): donor split that keeps
marrow eligible. The specified split (hold out 18-M-53 and 30-M-5) drops all
13 marrow units because marrow has only 3 mice at 30m (30-M-2, 30-M-3,
30-M-5); holding out 30-M-5 leaves 2. Here we hold out 18-M-53 and 30-M-4
instead: 30-M-4 has no marrow cells, so every tissue keeps exactly 3 mice
per side and the frozen >=3 rule still holds. Same frozen derivation rules,
same direction check for both held-out mice.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import scripts.run_mature_heldout_validation as hv
import scripts.run_mature_p1m_derivation as deriv


def main():
    for tissue in hv.TISSUES:
        hv.EXPR[tissue] = hv.build_expression_matrix(tissue)
    held_out = ["18-M-53", "30-M-4"]
    _, source_mice, _ = hv.get_eligibility()
    all_mice = sorted(source_mice["mouse.id"].astype(str).unique())
    allow = [m for m in all_mice if m not in held_out]
    print(f"held out: {held_out}; derive mice: {allow}", flush=True)
    genes, ranking = hv.derive_signature(hv.TISSUES, allow_mice=allow, tag="donor-split-3b")
    genes.to_csv(hv.RESULTS / "followup-mature-heldout-donor-split-3b-genes.csv", index=False)
    ranking.to_csv(hv.RESULTS / "followup-mature-heldout-donor-split-3b-unit-ranking.csv", index=False)
    print("3b signature:", len(genes), "genes;",
          int((genes["direction"] == "UP").sum()), "up,",
          int((genes["direction"] == "DOWN").sum()), "down", flush=True)
    print("3b eligible units:", len(ranking), flush=True)
    print(ranking.sort_values("n_de_uniform_before_clean", ascending=False).head(12)[
        ["unit", "n_de_uniform_before_clean"]].to_string(index=False), flush=True)
    up_genes = genes.loc[genes["direction"] == "UP", "mouse_gene_symbol"].tolist()
    down_genes = genes.loc[genes["direction"] == "DOWN", "mouse_gene_symbol"].tolist()
    rows = []
    for tissue in hv.TISSUES:
        genes_arr, pairs, matrix = hv.EXPR[tissue]
        units = ranking[ranking["tissue"] == tissue]["cell_ontology_class"].tolist()
        for ctype in units:
            sub = pairs[(pairs["cell_type"] == ctype)
                        & (pairs["age"].isin(deriv.PRIMARY_AGES))
                        & (pairs["n_cells"] >= hv.MIN_CELLS_PER_DONOR_CT)]
            if len(sub) == 0:
                continue
            scores, _, _ = hv.score_from_matrix(
                matrix[:, sub.index.to_numpy()], genes_arr, up_genes, down_genes)
            if scores is None:
                continue
            frame = sub.copy()
            frame["score"] = scores
            for held in held_out:
                hrow = frame[frame["mouse_id"] == held]
                if len(hrow) == 0:
                    continue
                held_age = hrow["age"].iloc[0]
                held_score = float(hrow["score"].iloc[0])
                derive = frame[frame["mouse_id"] != held]
                mean_18 = float(derive.loc[derive["age"] == "18m", "score"].mean())
                mean_30 = float(derive.loc[derive["age"] == "30m", "score"].mean())
                match = (held_score > mean_18) if held_age == "30m" else (held_score < mean_30)
                rows.append({"unit": f"{tissue}::{ctype}", "tissue": tissue,
                             "cell_type": ctype, "held_out_mouse": held,
                             "held_out_age": held_age, "held_out_score": held_score,
                             "derive_mean_18m": mean_18, "derive_mean_30m": mean_30,
                             "direction_match": bool(match)})
    out = pd.DataFrame(rows)
    out.to_csv(hv.RESULTS / "followup-mature-heldout-donor-split-3b.csv", index=False)
    print(f"3b: {int(out['direction_match'].sum())} of {len(out)} direction matches", flush=True)
    print(out.groupby("tissue")["direction_match"].agg(["sum", "count"]).to_string(), flush=True)


if __name__ == "__main__":
    main()
