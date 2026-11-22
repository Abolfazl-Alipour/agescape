"""Tests for the mature M10 PsychAD pipeline at 30-plus.

Checks row counts, hand-recomputed BH q values (independent implementation),
the exact cross-check of fully adjusted slopes against the committed Item 2
transport CSV, the 200-draw empirical p floor, and the verdict-change counts.
Run from the repo root: python3 -m pytest tests/test_followup_mature_psychad_m10.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"


def bh_manual(pvals):
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def test_row_counts():
    assert len(pd.read_csv(RESULTS / "followup-mature-psychad-m10-regressions.csv")) == 54
    assert len(pd.read_csv(RESULTS / "followup-mature-psychad-m10-confounder-audit.csv")) == 54
    assert len(pd.read_csv(RESULTS / "followup-mature-psychad-m10-cohort-summary.csv")) == 2
    assert len(pd.read_csv(RESULTS / "followup-mature-psychad-m10-null200.csv")) == 54


def test_bh_values_match_hand_recomputation():
    reg = pd.read_csv(RESULTS / "followup-mature-psychad-m10-regressions.csv")
    for pcol, qcol in (("p", "q_bh"), ("p_covar_adj", "q_bh_covar_adj"), ("p_raw_counts", "q_bh_raw")):
        for tier in ("P0m", "P1m"):
            sub = reg[reg.program_tier == tier]
            rec = bh_manual(sub[pcol])
            assert np.max(np.abs(rec - sub[qcol].to_numpy(float))) < 1e-12
    aud = pd.read_csv(RESULTS / "followup-mature-psychad-m10-confounder-audit.csv")
    for pcol in ("p_all", "p_covar_adj", "p_age_ge40", "p_hbcc", "p_mssm"):
        qcol = pcol.replace("p_", "q_")
        for tier in ("P0m", "P1m"):
            sub = aud[aud.program_tier == tier]
            valid = sub[pcol].notna()
            rec = bh_manual(sub.loc[valid, pcol])
            assert np.max(np.abs(rec - sub.loc[valid, qcol].to_numpy(float))) < 1e-12
    nul = pd.read_csv(RESULTS / "followup-mature-psychad-m10-null200.csv")
    for tier in ("P0m", "P1m"):
        sub = nul[nul.gene_set == tier]
        rec = bh_manual(sub["empirical_pvalue_200"])
        assert np.max(np.abs(rec - sub["empirical_q_200"].to_numpy(float))) < 1e-12


def test_covar_slopes_match_item2_transport():
    reg = pd.read_csv(RESULTS / "followup-mature-psychad-m10-regressions.csv")
    tr = pd.read_csv(RESULTS / "followup-mature-p1m-transport.csv")
    tr = tr[tr.cohort == "DLPFC_Lifespan_Atlas"].set_index(["unit", "gene_set"])
    diffs = [
        abs(r.slope_covar_adj - tr.loc[(r.cell_type, r.program_tier), "composite_slope"])
        for r in reg.itertuples()
    ]
    assert max(diffs) < 1e-12
    qdiffs = [
        abs(r.q_composite_mature170 - tr.loc[(r.cell_type, r.program_tier), "composite_q"])
        for r in reg.itertuples()
    ]
    assert max(qdiffs) < 1e-12


def test_null200_floor_and_verdict_changes():
    nul = pd.read_csv(RESULTS / "followup-mature-psychad-m10-null200.csv")
    assert abs(nul["empirical_pvalue_200"].min() - 1 / 201) < 1e-12
    assert nul.groupby("gene_set")["empirical_pass_200"].sum().to_dict() == {"P0m": 0, "P1m": 0}
    assert nul.groupby("gene_set")["empirical_pass_2000"].sum().to_dict() == {"P0m": 3, "P1m": 1}
    assert nul["verdict_change"].value_counts().to_dict() == {"SAME_FAIL": 50, "PASS_2000_ONLY": 4}


def test_pooled_summary_donors():
    summ = pd.read_csv(RESULTS / "followup-mature-psychad-m10-cohort-summary.csv")
    assert set(summ["n_donors"]) == {201}
    assert set(summ["age_min"]) == {30.0}
    assert set(summ["n_cell_types"]) == {27}
