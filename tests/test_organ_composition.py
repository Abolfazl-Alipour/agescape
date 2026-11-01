from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent


def test_all_species_and_organs_are_tested():
    results = pd.read_csv(REPO / "results" / "organ-composition-age-tests.csv")
    assert set(results["species"]) == {"Mouse", "Human", "Naked mole-rat"}
    summary = pd.read_csv(REPO / "results" / "organ-composition-summary.csv")
    assert summary.loc[summary["species"] == "Mouse", "organ"].nunique() == 16
    assert summary.loc[summary["species"] == "Human", "cohort"].nunique() == 11
    assert results.loc[results["species"] == "Naked mole-rat", "cell_type"].nunique() == 4


def test_nmr_fractions_sum_to_one_and_permutation_floor_is_explicit():
    fractions = pd.read_csv(REPO / "results" / "nmr-skin-cell-fractions-resolution-0.4.csv", index_col=0)
    np.testing.assert_allclose(fractions.sum(axis=1), 1.0)
    results = pd.read_csv(REPO / "results" / "organ-composition-age-tests.csv")
    nmr = results[results["species"] == "Naked mole-rat"]
    assert nmr["pvalue"].min() >= 0.1
    assert set(nmr["cell_type"]) == {"basal", "intermediate", "corneum", "immune"}


def test_summary_uses_primary_models_only():
    summary = pd.read_csv(REPO / "results" / "organ-composition-summary.csv")
    assert len(summary) == 28
    assert set(summary["answer"]) <= {"age-associated fraction change detected", "no fraction change detected at FDR 0.05", "insufficient donor coverage"}
