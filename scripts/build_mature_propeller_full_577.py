#!/usr/bin/env python3
"""Assemble propeller p-values for the ENTIRE pooled mature fraction family
(577 tests: mouse 114, human 285, fly 162, worm 16) and apply ONE
Benjamini-Hochberg correction across all 577 propeller p-values.

Sources (both mature-only, both produced by the unchanged
scripts/run_reviewer_phase1_propeller.R: limma regression on
logit-transformed proportions with a 0.5 pseudocount, design age plus log
total cells plus sex when variable, empirical Bayes on the age coefficient,
donor as the unit for mouse/human, pool as the unit for fly/worm):
  - results/followup-mature-propeller.csv: 399 mammalian scan rows
    (mature inputs from scripts/prepare_mature_propeller.py plus the 9 HRCA
    retina suites from scripts/prepare_mature_propeller_hrca.py).
  - results/followup-mature-fly-worm-propeller.csv: 178 fly/worm scan rows
    (mature inputs from scripts/prepare_mature_fly_worm_propeller.py).
    Worm rows carry the platform term the R script adds when it varies;
    fly rows have no platform column. This matches the committed fly/worm
    run exactly, so the p-values are reused as-is.

NOT reused: results/reviewer-phase1-composition-propeller.csv (612 rows).
That file was built from scripts/prepare_reviewer_phase1_propeller.py,
which loads core_inputs/extended_inputs from
scripts/run_reviewer_phase1_composition (full-range young-to-old inputs;
it includes zebrafish, naked mole rat, and mouse tissues the mature run
skips). Empirically its p-values differ from the mature runs on 398 of
the 577 scan rows, so it is a different design and is excluded.

Outputs: results/followup-mature-propeller-full-577.csv, one row per scan
row, with q_propeller_pooled_577 (single BH across all 577 propeller
p-values) and a status column (entered / dropped / unchanged) versus the
scan's q_pooled_mature_bh pass at 0.05.
"""

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"

SCAN = RESULTS / "mature-all-species-composition.csv"
MAT = RESULTS / "followup-mature-propeller.csv"
FLYWORM = RESULTS / "followup-mature-fly-worm-propeller.csv"
OUT = RESULTS / "followup-mature-propeller-full-577.csv"


def bh_qvalues(pvals):
    """Benjamini-Hochberg adjusted p-values, pure Python.

    q[i] = min over j >= rank(i) of n / j * p_sorted[j], monotone enforced.
    """
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running_min = min(running_min, pvals[i] * n / rank)
        q[i] = running_min
    return q


def main():
    scan = pd.read_csv(SCAN)
    mat = pd.read_csv(MAT)
    fw = pd.read_csv(FLYWORM)
    print(f"scan rows: {len(scan)}", flush=True)
    print(f"mature mammalian propeller rows: {len(mat)}", flush=True)
    print(f"fly/worm propeller rows: {len(fw)}", flush=True)
    assert len(scan) == 577, f"scan is not 577 rows: {len(scan)}"
    assert len(scan.drop_duplicates(
        subset=["species", "organ", "cohort", "cell_type"])) == 577, "scan keys not unique"

    key_cols = ["species", "organ", "cohort", "cell_type"]
    mat_index = {(r.species, r.organ, r.cohort, r.cell_type): r for r in mat.itertuples()}
    fw_index = {(r.species, r.organ, r.cohort, r.cell_type): r for r in fw.itertuples()}
    assert len(mat_index) == len(mat), "duplicate keys in mature propeller file"
    assert len(fw_index) == len(fw), "duplicate keys in fly/worm propeller file"
    assert not set(mat_index) & set(fw_index), "mammalian and fly/worm propeller overlap"

    rows = []
    missing = []
    for s in scan.itertuples():
        k = (s.species, s.organ, s.cohort, s.cell_type)
        prow = mat_index.get(k) or fw_index.get(k)
        if prow is None:
            missing.append(k)
            continue
        rows.append({
            "species": s.species, "organ": s.organ, "cohort": s.cohort,
            "cell_type": s.cell_type,
            "scan_effect": float(s.effect),
            "scan_pvalue": float(s.pvalue),
            "scan_q_pooled_mature_bh": float(s.q_pooled_mature_bh),
            "scan_pass": bool(float(s.q_pooled_mature_bh) < 0.05),
            "propeller_effect_logit": float(prow.effect_logit),
            "propeller_se": float(prow.se),
            "propeller_pvalue": float(prow.pvalue),
            "propeller_q_within_cohort": float(prow.q_bh),
            "propeller_source": ("followup-mature-propeller.csv"
                                 if k in mat_index
                                 else "followup-mature-fly-worm-propeller.csv"),
        })
    assert not missing, f"{len(missing)} scan rows lack a propeller match: {missing[:5]}"
    assert len(rows) == 577, f"expected 577 assembled rows, got {len(rows)}"
    frame = pd.DataFrame(rows)
    assert frame["propeller_pvalue"].notna().all(), "missing propeller p-value"
    assert ((frame["propeller_pvalue"] > 0) & (frame["propeller_pvalue"] <= 1)).all(), \
        "propeller p-value out of range"

    # Single BH correction across all 577 propeller p-values.
    q577 = bh_qvalues(frame["propeller_pvalue"].tolist())
    frame["q_propeller_pooled_577"] = q577
    frame["propeller_pass"] = frame["q_propeller_pooled_577"] < 0.05
    frame["sign_agree"] = (frame["scan_effect"] > 0) == (frame["propeller_effect_logit"] > 0)
    frame["status"] = frame.apply(
        lambda r: ("entered" if (not r["scan_pass"] and r["propeller_pass"]) else
                   ("dropped" if (r["scan_pass"] and not r["propeller_pass"]) else
                    "unchanged")), axis=1)

    # Independent cross-check of the BH values with statsmodels when available.
    try:
        from statsmodels.stats.multitest import multipletests
        expected = multipletests(frame["propeller_pvalue"].to_numpy(), method="fdr_bh")[1]
        diff = abs(expected - frame["q_propeller_pooled_577"].to_numpy()).max()
        print(f"BH cross-check vs statsmodels max abs diff: {diff:.3e}", flush=True)
        assert diff <= 1e-12, "BH cross-check mismatch"
    except ImportError:
        print("statsmodels unavailable; BH cross-check skipped", flush=True)

    frame.to_csv(OUT, index=False)
    print(f"wrote {len(frame)} rows to {OUT}", flush=True)
    print(f"propeller passes at q<0.05: {int(frame['propeller_pass'].sum())} of 577", flush=True)
    print("per-species passes:", flush=True)
    for sp, g in frame.groupby("species"):
        print(f"  {sp}: {int(g['propeller_pass'].sum())} of {len(g)}", flush=True)
    print(f"sign agreement: {int(frame['sign_agree'].sum())} of 577 "
          f"({frame['sign_agree'].mean():.4f})", flush=True)
    print(f"entered: {int((frame['status'] == 'entered').sum())}, "
          f"dropped: {int((frame['status'] == 'dropped').sum())}, "
          f"unchanged: {int((frame['status'] == 'unchanged').sum())}", flush=True)
    for _, r in frame[frame["status"] != "unchanged"].iterrows():
        print(f"  {r['status'].upper()}: {r['species']} {r['organ']} {r['cohort']} "
              f"{r['cell_type']} | scan q {r['scan_q_pooled_mature_bh']:.4g} -> "
              f"prop q577 {r['q_propeller_pooled_577']:.4g}", flush=True)


if __name__ == "__main__":
    main()
