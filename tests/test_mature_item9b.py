from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.run_mature_between_subject import CONTRASTS
from scripts.run_mature_pbmc import DESIGN, select_donors

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def test_mature_pbmc_donor_selection():
    cache = {
        "donor": np.array(["a", "b", "c", "d", "e", "f"]),
        "age": np.array([29.0, 30.0, 39.0, 40.0, 60.0, 85.0]),
        "sex": np.array(["female"] * 6),
        "cmv": np.array([""] * 6),
        "batch": np.array([""] * 6),
        "pool": np.array([""] * 6),
    }
    selected = select_donors(cache)
    assert selected["donor"].tolist() == ["b", "c", "e", "f"]
    assert selected["old"].tolist() == [0.0, 0.0, 1.0, 1.0]
    assert selected["age"].min() >= 30
    assert DESIGN == "mature_30_39_vs_60_plus"


def test_mature_contrasts_keep_old_mid_separate():
    assert CONTRASTS[0][0] == "OLD_vs_MATURE_30_39"
    assert CONTRASTS[0][1:3] == ("YOUNG", "OLD")
    assert CONTRASTS[-1][0] == "OLD_vs_MID"
    assert CONTRASTS[-1][1:3] == ("MID", "OLD")


def test_mature_part_a_output_contract_when_present():
    for dataset in ("gao", "yazar"):
        path = RESULTS / f"followup-mature-between-subject-{dataset}-donors.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        assert frame["age_years"].min() >= 30
        assert set(frame["mature_bin"]) <= {"MATURE_30_39", "MID_40_59", "OLD_60_PLUS"}
        valid = frame["a1_dist_full_1500"].notna()
        assert frame.loc[valid, "reference_age_min"].min() >= 30
        assert frame.loc[valid, "reference_age_max"].max() <= 39
        assert frame.loc[valid, "a1_full_1500_source"].eq(
            "results/followup-mature-nov17-donor-spread-distances.csv"
        ).all()


def test_mature_pbmc_p1m_contract_when_present():
    paths = list(RESULTS.glob("followup-mature-pbmc-GSE*-p1m-donor-scores.csv"))
    if not paths:
        pytest.skip("Mature PBMC outputs have not been generated")
    for path in paths:
        frame = pd.read_csv(path)
        assert frame["age"].min() >= 30
        assert ((frame["age"] <= 39) | (frame["age"] >= 60)).all()
    for path in RESULTS.glob("followup-mature-pbmc-GSE*-p1m-models.csv"):
        models = pd.read_csv(path)
        assert models["p1m_source"].eq("results/followup-mature-p1m-genes.csv").all()
        assert models["endpoint"].eq("P1m_transport").all()
