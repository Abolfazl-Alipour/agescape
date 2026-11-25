"""Tests for the SAMap match rule (scripts/run_samap_pair.py:apply_match_rule).

The rule is pre-stated and fixed before seeing results; a pair passes only if
ALL four hold:
  1. alignment score >= 0.30
  2. reciprocal best hit at cluster level (each cluster's top partner is the other)
  3. >= 40 enriched gene pairs
  4. score in the top 1% of all pairs in that comparison

Run: python -m pytest tests/test_samap_match_rule.py -v
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from run_samap_pair import apply_match_rule  # noqa: E402


def make_scores(cells):
    """cells: dict {(sid, cluster): {(sid, cluster): score}} -> wide DataFrame."""
    labels = sorted(cells)
    frame = pd.DataFrame(0.0, index=[f"{s}_{c}" for s, c in labels],
                         columns=[f"{s}_{c}" for s, c in labels])
    for (s1, c1), partners in cells.items():
        for (s2, c2), score in partners.items():
            frame.loc[f"{s1}_{c1}", f"{s2}_{c2}"] = score
    return frame


def test_passing_pair():
    # Two species, two clusters each. A1<->B1 is the reciprocal best pair with
    # a high score and 50 gene pairs; A2<->B2 is a weaker reciprocal pair.
    scores = make_scores({
        ("hu", "A1"): {("ms", "B1"): 0.50, ("ms", "B2"): 0.10},
        ("hu", "A2"): {("ms", "B1"): 0.05, ("ms", "B2"): 0.40},
        ("ms", "B1"): {("hu", "A1"): 0.50, ("hu", "A2"): 0.05},
        ("ms", "B2"): {("hu", "A1"): 0.10, ("hu", "A2"): 0.40},
    })
    counts = {"hu_A1;ms_B1": 50, "hu_A2;ms_B2": 50}
    out = apply_match_rule(scores, counts)
    assert len(out) == 4  # one row per unordered pair
    row = out[(out["cluster1"] == "A1")].iloc[0]
    assert bool(row["passes"])
    assert bool(row["reciprocal_best"])
    assert row["n_gene_pairs"] == 50


def test_below_alignment_threshold_fails():
    scores = make_scores({
        ("hu", "A1"): {("ms", "B1"): 0.29, ("ms", "B2"): 0.10},
        ("hu", "A2"): {("ms", "B1"): 0.05, ("ms", "B2"): 0.28},
        ("ms", "B1"): {("hu", "A1"): 0.29, ("hu", "A2"): 0.05},
        ("ms", "B2"): {("hu", "A1"): 0.10, ("hu", "A2"): 0.28},
    })
    out = apply_match_rule(scores, {"hu_A1;ms_B1": 100, "hu_A2;ms_B2": 100})
    assert not out["passes"].any()


def test_non_reciprocal_best_fails():
    # A1's best is B1, but B1's best is A2, so A1-B1 is not reciprocal.
    scores = make_scores({
        ("hu", "A1"): {("ms", "B1"): 0.60, ("ms", "B2"): 0.10},
        ("hu", "A2"): {("ms", "B1"): 0.70, ("ms", "B2"): 0.20},
        ("ms", "B1"): {("hu", "A1"): 0.60, ("hu", "A2"): 0.70},
        ("ms", "B2"): {("hu", "A1"): 0.10, ("hu", "A2"): 0.20},
    })
    out = apply_match_rule(scores, {"hu_A1;ms_B1": 100, "hu_A2;ms_B1": 100})
    row = out[(out["cluster1"] == "A1")].iloc[0]
    assert not bool(row["reciprocal_best"])
    assert not bool(row["passes"])


def test_fewer_than_40_gene_pairs_fails():
    scores = make_scores({
        ("hu", "A1"): {("ms", "B1"): 0.80, ("ms", "B2"): 0.10},
        ("hu", "A2"): {("ms", "B1"): 0.05, ("ms", "B2"): 0.70},
        ("ms", "B1"): {("hu", "A1"): 0.80, ("hu", "A2"): 0.05},
        ("ms", "B2"): {("hu", "A1"): 0.10, ("hu", "A2"): 0.70},
    })
    out = apply_match_rule(scores, {"hu_A1;ms_B1": 39, "hu_A2;ms_B2": 100})
    row = out[(out["cluster1"] == "A1")].iloc[0]
    assert row["n_gene_pairs"] == 39
    assert not bool(row["passes"])


def test_top_1pct_gate():
    # 150 pairs: one strong pair (0.9) and 149 weak pairs (0.05). Top 1% of
    # 150 is the top 1.5, so only the 0.9 pair clears the top-1% gate.
    cells = {}
    for i in range(1, 76):
        cells[("hu", f"A{i}")] = {("ms", f"B{i}"): 0.05, ("ms", "B1"): 0.05}
        cells[("ms", f"B{i}")] = {("hu", f"A{i}"): 0.05, ("hu", "A1"): 0.05}
    cells[("hu", "A1")][("ms", "B1")] = 0.9
    cells[("ms", "B1")][("hu", "A1")] = 0.9
    scores = make_scores(cells)
    counts = {f"hu_A{i};ms_B{i}": 100 for i in range(1, 76)}
    out = apply_match_rule(scores, counts)
    passing = out[out["passes"]]
    assert len(passing) == 1
    assert passing.iloc[0]["cluster1"] == "A1"


def test_empty_scores():
    out = apply_match_rule(pd.DataFrame(), {})
    assert out.empty
