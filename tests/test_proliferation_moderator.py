import numpy as np
import pandas as pd

from scripts.run_proliferation_moderator import design_matrix, permutation_pvalue, unit_effects
from src.proliferation import MEIOTIC_OR_MIXED, POST_MITOTIC, PROLIFERATIVE, UNCERTAIN, classify_cell_type


def test_biological_classification_distinguishes_mature_cells_from_progenitors():
    assert classify_cell_type("GABAergic neuron", "Human")["proliferation_class"] == POST_MITOTIC
    assert classify_cell_type("oligodendrocyte precursor cell", "Human")["proliferation_class"] == PROLIFERATIVE
    assert classify_cell_type("skeletal muscle satellite cell", "Mouse")["proliferation_class"] == PROLIFERATIVE
    assert classify_cell_type("skeletal muscle cell", "Mouse")["proliferation_class"] == POST_MITOTIC
    assert classify_cell_type("spermatocytes 1", "Zebrafish")["proliferation_class"] == MEIOTIC_OR_MIXED
    assert classify_cell_type("bladder cell", "Mouse")["proliferation_class"] == UNCERTAIN


def test_unit_effects_use_donor_rows_and_lifespan_scaling():
    rows = []
    for cell_type in ("neuron", "fibroblast"):
        for index in range(8):
            rows.append({
                "species": "Mouse", "tissue": "Lung", "cohort": "TMS", "cell_type": cell_type,
                "donor": str(index), "age": index / 7, "sex_male": 0,
                "mean_log_umi": 7.0, "n_cells": 20,
                "outcome": 1 + index * (0.1 if cell_type == "neuron" else 0.2),
            })
    result = unit_effects(pd.DataFrame(rows), "expression_vector", "test")
    expected_sd = np.std(np.arange(8), ddof=1)
    for row in result.itertuples():
        assert np.isclose(row.effect_sd, 0.4 * 7 / expected_sd)
        assert row.n_donors == 8


def test_strata_with_only_one_class_do_not_contribute_to_contrast():
    frame = pd.DataFrame({
        "species": ["Mouse"] * 4,
        "tissue": ["brain", "brain", "lung", "lung"],
        "cohort": ["TMS"] * 4,
        "cell_type": ["neuron 1", "neuron 2", "fibroblast 1", "fibroblast 2"],
        "proliferation_class": [POST_MITOTIC, POST_MITOTIC, PROLIFERATIVE, PROLIFERATIVE],
        "effect_sd": [1.0, 2.0, 3.0, 4.0],
    })
    assert design_matrix(frame, "effect_sd") is None
    assert np.isnan(permutation_pvalue(frame, "effect_sd", permutations=10))


def test_committed_classification_covers_every_tested_unit():
    units = pd.read_csv("results/proliferation-unit-effects.csv")
    curation = pd.read_csv("results/proliferation-cell-type-classification.csv")
    merged = units.merge(curation, on=["species", "cell_type"], suffixes=("", "_curated"))
    assert len(merged) == len(units)
    assert merged["proliferation_class_curated"].notna().all()
    assert (merged["proliferation_class"] == merged["proliferation_class_curated"]).all()
