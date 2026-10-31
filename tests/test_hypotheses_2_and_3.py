"""
tests/test_hypotheses_2_and_3.py

Verifies mathematical invariants, data integrity, and reproducibility
for Hypothesis 1 (Epigenetic Clock), Hypothesis 2 (Composition Shifts),
Hypothesis 3 (Cellular Entropy), and Fibroblast Marker Audits.
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
    assert (res_dir / "h2-composition-cross-species-homologs.csv").exists()


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
    comp = pd.read_csv(REPO_ROOT / "results" / "h2-composition-cross-species-homologs.csv")
    assert len(comp) >= 7, "Expected at least 7 homologous cell-type pairs"
    assert not comp["mouse_slope_per_month"].isna().any()
    assert not comp["human_slope_per_year"].isna().any()
    assert "sign_concordance" in comp.columns


def test_h3_entropy_outputs_exist_and_valid():
    res_dir = REPO_ROOT / "results"
    out_csv = res_dir / "h3-entropy-mouse-lung.csv"
    assert out_csv.exists(), "Mouse lung entropy CSV missing!"
    df = pd.read_csv(out_csv)
    assert len(df) >= 7, f"Expected at least 7 tested cell types, got {len(df)}"
    assert not df["dm_slope_per_month"].isna().any()
    assert not df["dm_pval"].isna().any()
    assert not df["dispersion_slope_per_month"].isna().any()

    # Fibroblasts: positive unadjusted DM slope
    fibro = df[df["cell_type"] == "fibroblast of lung"].iloc[0]
    assert fibro["dm_slope_per_month"] > 0
    assert fibro["dm_pval"] < 0.05

    # Classical monocytes: non-significant DM slope
    mono = df[df["cell_type"] == "classical monocyte"].iloc[0]
    assert mono["dm_pval"] > 0.5


def test_fibroblast_markers_audit():
    res_dir = REPO_ROOT / "results"
    out_csv = res_dir / "fibroblast-markers-audit.csv"
    assert out_csv.exists(), "fibroblast-markers-audit.csv missing!"
    df = pd.read_csv(out_csv).set_index("gene")
    assert "Cdkn1a" in df.index
    assert "Sirt6" in df.index
    # Cdkn1a (p21) upregulation survives FDR
    assert df.loc["Cdkn1a", "slope_per_month"] > 0
    assert df.loc["Cdkn1a", "q_value_bh"] < 0.05
    # Sirt6 downregulation
    assert df.loc["Sirt6", "slope_per_month"] < 0
    assert df.loc["Sirt6", "p_value"] < 0.05


def test_h2_and_h3_figures_exist():
    fig_dir = REPO_ROOT / "figures"
    assert (fig_dir / "h2-composition-shifts-mouse-vs-human.png").exists()
    assert (fig_dir / "h3-transcriptional-noise-vs-age.png").exists()
    assert (fig_dir / "fibroblast-candidate-markers.png").exists()


def test_h1_epigenetic_outputs_exist_and_accurate():
    res_dir = REPO_ROOT / "results"
    fig_dir = REPO_ROOT / "figures"
    assert (res_dir / "h1-epigenetic-clock-human.csv").exists()
    assert (res_dir / "h1-epigenetic-clock-mouse.csv").exists()
    assert (res_dir / "h1-epigenetic-cross-species-summary.csv").exists()
    assert (fig_dir / "h1-pan-mammalian-epigenetic-clock.png").exists()

    # Verify Mouse Liver with corrected link function
    m_df = pd.read_csv(res_dir / "h1-epigenetic-clock-mouse.csv")
    r_mouse = float(m_df.iloc[0]["pearson_r"])
    mae_mouse = float(m_df.iloc[0]["mae_years"])
    assert r_mouse > 0.90, f"Expected mouse r > 0.90, got {r_mouse}"
    assert mae_mouse < 0.50, f"Expected mouse MAE < 0.50 years (with link function), got {mae_mouse}"
    assert int(m_df.iloc[0]["n_samples"]) == 339

    # Verify Human Liver and All Tissues with corrected link function
    h_df = pd.read_csv(res_dir / "h1-epigenetic-clock-human.csv").set_index("tissue")
    r_h_liver = float(h_df.loc["Liver", "pearson_r"])
    r_h_all = float(h_df.loc["All Tissues", "pearson_r"])
    mae_h_all = float(h_df.loc["All Tissues", "mae_years"])
    assert r_h_liver > 0.80, f"Expected human liver r > 0.80, got {r_h_liver}"
    assert r_h_all > 0.75, f"Expected human all tissues r > 0.75, got {r_h_all}"
    assert mae_h_all < 10.0, f"Expected human MAE < 10.0 years (with link function), got {mae_h_all}"
    assert int(h_df.loc["All Tissues", "n_samples"]) == 661


def test_h2_prethreshold_counts_and_presence_outputs():
    res_dir = REPO_ROOT / "results"
    counts = pd.read_csv(res_dir / "m9-donor-cell-counts-prethreshold.csv")
    lung = counts[(counts["cohort"] == "Lung_Integrated") & (counts["age"] >= 20)]
    assert ((lung["n_cells"] > 0) & (lung["n_cells"] < 20)).any()
    primary = pd.read_csv(res_dir / "h2-composition-human-lung.csv")
    observed = pd.read_csv(res_dir / "h2-composition-human-lung-observed-only.csv")
    presence = pd.read_csv(res_dir / "h2-composition-human-lung-presence.csv")
    assert len(primary) == len(observed) == len(presence)
    assert primary["n_donors"].nunique() == 1
    assert (observed["n_donors"] <= primary["n_donors"].iloc[0]).all()
    assert (presence["n_present"] + presence["n_absent"] == primary["n_donors"].iloc[0]).all()


def test_h1_predictive_r2_recomputes_from_per_sample_predictions():
    for species, age_column in (("human", "age"), ("mouse", "age_years")):
        predictions = pd.read_csv(REPO_ROOT / "results" / f"h1-epigenetic-clock-{species}-predictions.csv")
        observed = predictions[age_column].to_numpy()
        predicted = predictions["pred_age_years"].to_numpy()
        expected = 1 - np.sum((observed - predicted) ** 2) / np.sum((observed - observed.mean()) ** 2)
        summary = pd.read_csv(REPO_ROOT / "results" / f"h1-epigenetic-clock-{species}.csv")
        if species == "human":
            summary = summary[summary["tissue"] == "All Tissues"]
        assert np.isclose(summary.iloc[0]["predictive_r2"], expected)
        assert "r2" not in summary.columns
