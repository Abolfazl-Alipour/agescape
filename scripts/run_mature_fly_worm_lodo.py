#!/usr/bin/env python3
"""Leave-one-pool-out (LODO) for the 26 non-mammalian passing rows of the pooled mature fraction scan.

Refits the exact scan model (CLR with 0.5 pseudocount, OLS on age plus sex,
sex dropped when rank-deficient; slopes per day) with each pool removed in
turn. Mirrors the mammalian LODO in scripts/run_mature_composition_complete.py
(secondary_checks): the CLR table is computed once on the full stratum and
restricted to the kept pools; the eligible cell-type set is not recomputed.
The inferential unit is the pool: fly 72 pools (36 head, 36 body), worm 9 pools.

Writes results/followup-mature-fly-worm-lodo.csv. Verifies every full fit
against the committed scan effects, and hand-recomputes one row from scratch
(full fit and one leave-one-out fit) with numpy least squares on a CLR table
computed independently from the raw counts.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_mature_all_species_fractions import (
    MIN_PRESENT, PSEUDOCOUNT, fit_cell_type, load_strata,
)
from scripts.run_mature_composition_complete import clr_table

SCAN = REPO / "results" / "mature-all-species-composition.csv"
OUTPUT = REPO / "results" / "followup-mature-fly-worm-lodo.csv"


def hand_clr_slope(counts, metadata, cell_type, keep=None):
    """Independent slope from scratch: CLR from raw counts, OLS via lstsq."""
    sub = counts if keep is None else counts.loc[keep]
    meta = metadata.loc[sub.index]
    adj = sub.astype(float) + 0.5
    props = adj.div(adj.sum(axis=1), axis=0)
    logp = np.log(props)
    clr = logp.sub(logp.mean(axis=1), axis=0)
    y = clr[cell_type].to_numpy(float)
    age = meta["age"].to_numpy(float)
    sex = meta["sex"].to_numpy(float)
    if len(np.unique(sex)) > 1:
        X = np.column_stack([np.ones(len(y)), age, sex])
    else:
        X = np.column_stack([np.ones(len(y)), age])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta[1]


def main():
    scan = pd.read_csv(SCAN)
    passing = scan[(scan["significant_pooled_05"]) & (scan["species"].isin(["Fruit fly", "C. elegans"]))].copy()
    assert len(passing) == 26, f"expected 26 passing fly/worm rows, got {len(passing)}"
    strata = {(s, o, c): (counts, meta) for s, o, c, counts, meta in load_strata()
              if s in ("Fruit fly", "C. elegans")}
    records = []
    hand_done = False
    max_full_diff = 0.0
    for row in passing.itertuples():
        key = (row.species, row.organ, row.cohort)
        counts_all, metadata = strata[key]
        present = (counts_all > 0).sum(axis=0)
        counts = counts_all.loc[:, present >= MIN_PRESENT].copy()
        clr = clr_table(counts, pseudocount=PSEUDOCOUNT)
        full, reason = fit_cell_type(clr, counts, metadata, row.cell_type, *key)
        assert full is not None, f"full fit failed for {key} {row.cell_type}: {reason}"
        max_full_diff = max(max_full_diff, abs(full["effect"] - row.effect))
        if not hand_done:
            hand_full = hand_clr_slope(counts, metadata, row.cell_type)
            first_pool = counts.index[0]
            hand_lodo = hand_clr_slope(counts, metadata, row.cell_type, keep=counts.index != first_pool)
            refit, _ = fit_cell_type(clr.loc[counts.index != first_pool],
                                     counts.loc[counts.index != first_pool],
                                     metadata, row.cell_type, *key)
            print(f"hand check {key} {row.cell_type}: scan={row.effect:.12g} "
                  f"fit={full['effect']:.12g} hand={hand_full:.12g} "
                  f"lodo_drop_{first_pool}: fit={refit['effect']:.12g} hand={hand_lodo:.12g}")
            assert abs(hand_full - full["effect"]) < 1e-9, "hand full fit disagrees"
            assert abs(hand_lodo - refit["effect"]) < 1e-9, "hand LODO fit disagrees"
            hand_done = True
        slopes = []
        for pool in counts.index:
            keep = counts.index != pool
            refit, _ = fit_cell_type(clr.loc[keep], counts.loc[keep], metadata, row.cell_type, *key)
            if refit is None:
                continue
            slopes.append({"pool": pool, "effect": refit["effect"], "pvalue": refit["pvalue"]})
        sdf = pd.DataFrame(slopes)
        assert len(sdf) == len(counts), f"LODO fits missing for {key} {row.cell_type}"
        changes = (sdf["effect"] - full["effect"]).abs()
        records.append({
            "species": row.species, "organ": row.organ, "cohort": row.cohort,
            "cell_type": row.cell_type,
            "scan_effect": row.effect, "scan_q_pooled": row.q_pooled_mature_bh,
            "n_pools": int(len(counts)), "lodo_n_fits": int(len(sdf)),
            "lodo_slope_min": float(sdf["effect"].min()),
            "lodo_slope_max": float(sdf["effect"].max()),
            "lodo_same_sign": bool((np.sign(sdf["effect"]) == np.sign(full["effect"])).all()),
            "lodo_max_abs_slope_change": float(changes.max()),
            "lodo_max_change_pool": str(sdf.loc[changes.idxmax(), "pool"]),
            "lodo_pvalue_max": float(sdf["pvalue"].max()),
        })
    out = pd.DataFrame(records)
    out.to_csv(OUTPUT, index=False)
    print(f"max |full-fit effect - scan effect| = {max_full_diff:.3g}")
    print(f"wrote {len(out)} rows to {OUTPUT}")
    print(f"sign holds in every refit: {int(out['lodo_same_sign'].sum())} of {len(out)} rows")


if __name__ == "__main__":
    main()
