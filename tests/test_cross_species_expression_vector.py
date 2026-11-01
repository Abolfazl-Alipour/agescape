import numpy as np
import pandas as pd

from scripts.run_cross_species_expression_vector import (
    build_reference,
    nmr_tests,
    reference_without_donor,
    technical_gene_mask,
)


def test_technical_gene_filter_across_species_naming():
    genes = ["MT-CO1", "RPLP0", "Mrpl12", "Xist", "DDX3Y", "ACTB"]
    np.testing.assert_array_equal(
        technical_gene_mask(genes), [True, True, True, True, True, False]
    )


def test_reference_excludes_evaluated_animal():
    matrix = np.array([[1, 0], [3, 0], [0, 4]])
    obs = pd.DataFrame({"cell_type": ["a", "a", "a"], "donor": ["young1", "young2", "old"]})
    reference = build_reference(matrix[:2], obs.iloc[:2], "cell_type", "donor")
    excluded = reference_without_donor(reference["a"], "young1")
    np.testing.assert_array_equal(excluded["counts"], [3, 0])
    assert excluded["n_cells"] == 1


def test_nmr_global_permutation_uses_animal_labels_not_row_order():
    donors = pd.DataFrame({
        "species": ["Naked mole-rat"] * 10,
        "tissue": ["Skin"] * 10,
        "cohort": ["GSE152886"] * 10,
        "cell_type": ["a"] * 5 + ["b"] * 5,
        "donor": ["1J", "5J", "6J", "3V", "4V"] * 2,
        "group": ["young", "young", "young", "aged", "aged"] * 2,
        "n_cells": [20] * 10,
        "profile_hellinger": [0.1, 0.2, 0.3, 0.4, 0.5, 0.2, 0.3, 0.4, 0.5, 0.6],
        "raw_count_euclidean": [10, 20, 30, 40, 50, 20, 30, 40, 50, 60],
    })
    tests, global_tests = nmr_tests(donors)
    assert len(tests) == 4
    assert set(global_tests["pvalue"]) == {0.2}
    assert (global_tests["effect_sd"] > 0).all()


def test_cross_species_outputs_trace_to_public_donor_rows():
    donors = pd.read_csv("results/h3-cross-species-vector-donors.csv")
    tests = pd.read_csv("results/h3-cross-species-vector-tests.csv")
    global_tests = pd.read_csv("results/h3-cross-species-vector-global.csv")
    assert set(donors["species"]) == {"Mouse", "Human", "Naked mole-rat"}
    assert set(global_tests["species"]) == {"Mouse", "Human", "Naked mole-rat"}
    nmr = tests[tests["species"] == "Naked mole-rat"]
    assert (nmr["pvalue"] >= 0.1).all()
    for row in tests.itertuples():
        expected = donors[
            (donors["species"] == row.species)
            & (donors["tissue"] == row.tissue)
            & (donors["cohort"] == row.cohort)
            & (donors["cell_type"] == row.cell_type)
        ]
        assert len(expected) == row.n_donors
        assert expected["n_cells"].sum() == row.n_cells


def test_global_model_keeps_tissues_as_independent_units():
    from scripts.run_cross_species_expression_vector import fit_global

    rows = []
    for tissue in ("lung", "skin"):
        for index in range(8):
            rows.append({
                "species": "Mouse", "tissue": tissue, "cohort": "TMS", "cell_type": "fibroblast",
                "donor": f"{tissue}-{index}", "age": 0.25 + index * 0.25,
                "sex_male": index % 2, "mean_log_umi": 7.0,
                "profile_hellinger": 0.1 + index * 0.01 + (tissue == "skin") * 0.2,
            })
    result = fit_global(pd.DataFrame(rows), "profile_hellinger", "Mouse", "All", "All", "unadjusted", 0.4)
    assert result["n_units"] == 2
    assert result["n_donors"] == 16
