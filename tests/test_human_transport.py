from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent


def test_primary_signature_tiers_and_global_multiplicity():
    m9 = pd.read_csv(REPO / "results" / "m9-regressions.csv")
    m10 = pd.read_csv(REPO / "results" / "m10-celltype-regressions.csv")
    family = pd.read_csv(REPO / "results" / "human-transport-global-multiplicity.csv")
    assert set(m9["program_tier"]) == {"P1_Clean"}
    assert set(m10["program_tier"]) == {"P1_Clean"}
    assert len(family) == len(m9) + len(m10)
    assert family["test_id"].is_unique
    assert family["q_global_bh"].between(0, 1).all()
    survivors = pd.read_csv(REPO / "results" / "human-transport-global-survivors.csv")
    assert (survivors["q_global_bh"] < 0.05).all()
    assert len(survivors) == int((family["q_global_bh"] < 0.05).sum())


def test_sensitivity_arms_use_distinct_tiers_and_normalization():
    m9_p0 = pd.read_csv(REPO / "results" / "m9-regressions-p0.csv")
    m9_raw = pd.read_csv(REPO / "results" / "m9-regressions-raw-counts.csv")
    m10_p0 = pd.read_csv(REPO / "results" / "m10-celltype-regressions-p0.csv")
    assert set(m9_p0["program_tier"]) == {"P0_Raw257"}
    assert set(m10_p0["program_tier"]) == {"P0_Raw257"}
    assert set(m9_raw["normalization"]) == {"log1p_raw_counts"}
