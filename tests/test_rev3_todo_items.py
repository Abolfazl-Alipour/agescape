"""Tests for Rev3 overnight to-do list scripts and deliverables."""

from pathlib import Path

import numpy as np
import pandas as pd

from scripts.run_rev3_composition_fly_worm_zebrafish import (
    clr_table,
    fit_f_dist_smyth,
    run_clr_regressions,
    run_propeller_limma,
)

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def test_clr_zero_sum_and_propeller_synthetic():
    counts = pd.DataFrame(
        {
            "cell_a": [100, 120, 150, 200, 220, 250, 300, 320],
            "cell_b": [400, 380, 350, 300, 280, 250, 200, 180],
            "cell_c": [50, 52, 48, 51, 49, 50, 52, 48],
        },
        index=[f"d{i}" for i in range(8)],
    )
    meta = pd.DataFrame(
        {
            "age": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            "sex_male": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
        },
        index=counts.index,
    )
    clr = clr_table(counts, pseudocount=0.5)
    assert np.allclose(clr.sum(axis=1).to_numpy(), 0.0, atol=1e-12)

    s02, df0 = fit_f_dist_smyth(np.array([0.05, 0.06, 0.07]), df1=5.0)
    assert s02 > 0
    assert df0 > 0

    clr_rows = run_clr_regressions(counts, meta, "Fruit fly", "head", "AFCA", "all", "primary", 0.5)
    assert len(clr_rows) == 3
    assert np.isclose(sum(r["effect"] for r in clr_rows), 0.0, atol=1e-12)

    prop_rows = run_propeller_limma(counts, meta, "Fruit fly", "head", "AFCA", "all", "propeller", 0.5)
    assert len(prop_rows) == 3
    cell_a = [r for r in prop_rows if r["cell_type"] == "cell_a"][0]
    assert cell_a["effect"] > 0
    assert cell_a["q_within_bh"] < 0.05


def test_item1_outputs_and_global_bh_family():
    comp = pd.read_csv(RESULTS / "rev3-composition-fly-worm-zebrafish.csv")
    primary = comp[comp["is_primary_test"]]
    assert len(primary) == 193
    assert set(comp.loc[comp["is_primary_stratum"], "pseudocount"].dropna().unique()) >= {0.1, 0.5, 1.0}
    assert (comp["sensitivity"] == "propeller_limma_primary_design").sum() >= 193

    g624 = pd.read_csv(RESULTS / "rev3-composition-global-bh-624.csv")
    assert len(g624) == 624
    assert int((g624["source_family"] == "existing_431_mouse_human_nmr").sum()) == 431
    assert int((g624["source_family"] == "new_193_fly_worm_zebrafish").sum()) == 193
