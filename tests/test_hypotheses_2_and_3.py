"""
tests/test_hypotheses_2_and_3.py

Verifies mathematical invariants, data integrity, and reproducibility
for Hypothesis 2 (Composition Shifts) and Hypothesis 3 (Cellular Entropy).
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


def test_h2_composition_outputs_exist():
    res_dir = REPO_ROOT / "results"
    assert (res_dir / "h2-composition-mouse-lung.csv").exists()
    assert (res_dir / "h2-composition-human-lung.csv").exists()
    assert (res_dir / "h2-composition-cross-species.csv").exists()


def test_h2_composition_simplex_invariant():
    """
    Simplex geometry invariant: Centred Log-Ratio (CLR) coordinates
    must sum to zero across all cell types for every single donor.
    """
    from scripts.run_h2_composition import compute_clr_table
    dummy_counts = pd.DataFrame({
        "cell_a": [10, 50, 100],
        "cell_b": [20, 25, 5],
        "cell_c": [5, 10, 20]
    }, index=["donor1", "donor2", "donor3"])
    clr = compute_clr_table(dummy_counts)
    np.testing.assert_allclose(clr.sum(axis=1).values, 0.0, atol=1e-10)


def test_h2_composition_results_content():
    comp = pd.read_csv(REPO_ROOT / "results" / "h2-composition-cross-species.csv")
    assert len(comp) == 11, "Expected 11 homologous cell-type pairs"
    assert not comp["mouse_slope_per_month"].isna().any()
    assert not comp["human_slope_per_year"].isna().any()
    
    # Verify that sign concordance was properly calculated
    concordant = comp["sign_concordance"].sum()
    assert concordant == 3, f"Expected 3 concordant pairs, got {concordant}"


def test_h3_entropy_outputs_exist_and_valid():
    res_dir = REPO_ROOT / "results"
    out_csv = res_dir / "h3-entropy-mouse-lung.csv"
    assert out_csv.exists(), "Mouse lung entropy CSV missing!"
    df = pd.read_csv(out_csv)
    assert len(df) == 7, f"Expected 7 tested cell types, got {len(df)}"
    assert not df["dm_slope_per_month"].isna().any()
    assert not df["dm_pval"].isna().any()

    # Fibroblasts must show significant positive noise increase (q < 0.05)
    fibro = df[df["cell_type"] == "fibroblast of lung"].iloc[0]
    assert fibro["dm_slope_per_month"] > 0
    assert fibro["dm_q_bh"] < 0.05

    # Classical monocytes must show flat noise (p > 0.5)
    mono = df[df["cell_type"] == "classical monocyte"].iloc[0]
    assert mono["dm_pval"] > 0.5


def test_h2_and_h3_figures_exist():
    fig_dir = REPO_ROOT / "figures"
    assert (fig_dir / "h2-composition-shifts-mouse-vs-human.png").exists()
    assert (fig_dir / "h3-transcriptional-noise-vs-age.png").exists()
