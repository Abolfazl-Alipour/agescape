import numpy as np
import pandas as pd

from scripts import run_rev2_joint_concordance as joint


def test_detection_mask_uses_donor_and_unit_floors():
    units = [
        {"expression": np.array([[1, 1, 0], [1, 0, 0], [1, 0, 0]])},
        {"expression": np.array([[1, 0, 0], [1, 0, 0], [1, 0, 0]])},
    ]
    mask, _ = joint.detection_mask(units, 3)
    assert mask.tolist() == [True, False, False]


def test_joint_age_permutation_preserves_cohort_and_sex():
    donors = ["A|d1", "A|d2", "B|d1", "B|d2"]
    ages = np.array([20.0, 80.0, 30.0, 90.0])
    strata = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    joint.CURRENT_DONORS = donors
    permuted = joint.permute_donor_ages(ages, strata, np.random.default_rng(7))
    assert sorted(permuted[:2]) == [20.0, 80.0]
    assert sorted(permuted[2:]) == [30.0, 90.0]


def test_species_analysis_uses_one_row_per_donor_per_unit():
    metadata = pd.DataFrame({
        "donor": ["d1", "d2", "d3", "d4", "d5", "d6"],
        "age": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        "sex_male": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        "mean_log_umi": np.linspace(7.0, 8.0, 6),
    })
    rng = np.random.default_rng(11)
    units = []
    for offset in (0.0, 1.0):
        expression = rng.normal(size=(6, 20)) + metadata["age"].to_numpy()[:, None] + offset
        units.append({"metadata": metadata.copy(), "expression": expression})
    strata = [("Tissue", "Cohort", np.array([f"g{i}" for i in range(20)]), units)]
    summary, rows = joint.analyze_species("Human", strata, permutations=20, seed=13)
    assert summary["n_donors"] == 6
    assert summary["n_units"] == 2
    assert summary["n_pairs"] == 1
    assert summary["n_permutations_joint"] == 20
    assert len(rows) == 1
    assert np.isfinite(summary["joint_null_mean"])
