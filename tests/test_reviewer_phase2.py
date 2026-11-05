import numpy as np
import pandas as pd

from scripts.run_reviewer_phase2_annotation import marker_for_cell_type


def test_depth_sweep_reproduces_original_point_estimates():
    depth = pd.read_csv("results/reviewer-phase2-depth-global.csv")
    original = pd.read_csv("results/h3-frequentist-global.csv")
    for metric, column in [("entropy", "entropy"), ("offprogram_fraction", "offprogram_fraction")]:
        current = depth[(depth["depth"].astype(str) == "1500") & (depth["metric"] == metric)
                        & (depth["model"] == "age_sex")].iloc[0]
        reference = original[(original["metric"] == column)
                             & (original["analysis"] == "all_cells_all_donor_ids")].iloc[0]
        assert np.isclose(current["effect_15mo_sd"], reference["effect_15mo_sd"])


def test_full_depth_and_adult_outputs_have_expected_scope():
    depth = pd.read_csv("results/reviewer-phase2-depth-global.csv")
    full = depth[(depth["depth"] == "full") & (depth["model"] == "age_sex")]
    assert set(full["metric"]) == {"entropy", "offprogram_fraction"}
    ages = pd.read_csv("results/reviewer-phase2-human-age-distributions.csv")
    assert ages["all_adult"].all()
    assert ages["age_min"].min() >= 20


def test_clock_overlap_remains_unresolved():
    audit = pd.read_csv("results/reviewer-phase2-clock-overlap-audit.csv").iloc[0]
    assert "unresolved" in audit["training_overlap_status"]
    assert audit["accession_mentions"] == 0
    assert audit["tested_sample_id_mentions"] == 0


def test_marker_panel_distinguishes_precursors_from_mature_oligodendrocytes():
    assert marker_for_cell_type("oligodendrocyte precursor cell") == "PDGFRA"
    assert marker_for_cell_type("oligodendrocyte") == "MBP"
    assert marker_for_cell_type("VIP GABAergic cortical interneuron") == "VIP"
    assert marker_for_cell_type("tracheobronchial goblet cell") == "MUC5AC"
