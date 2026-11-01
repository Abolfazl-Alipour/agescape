import numpy as np
import pandas as pd
import scipy.sparse as sp

from scripts.run_h3_expression_vector import hellinger_distance, raw_euclidean_distance


def test_hellinger_distance_known_values():
    matrix = sp.csr_matrix(np.array([[4, 0], [1, 1], [0, 4]], dtype=float))
    reference = np.array([0.5, 0.5])
    expected = [
        np.sqrt(1 - np.sqrt(0.5)),
        0.0,
        np.sqrt(1 - np.sqrt(0.5)),
    ]
    np.testing.assert_allclose(hellinger_distance(matrix, reference), expected)


def test_raw_euclidean_distance_known_values():
    matrix = sp.csr_matrix(np.array([[4, 0], [1, 1]], dtype=float))
    np.testing.assert_allclose(raw_euclidean_distance(matrix, [1, 1]), [np.sqrt(10), 0.0])


def test_expression_vector_results_trace_to_animal_rows():
    results = pd.read_csv("results/h3-expression-vector.csv")
    donors = pd.read_csv("results/h3-expression-vector-donors.csv")
    primary = results[
        (results["donor_set"] == "all_donor_ids")
        & (results["model"] == "age_plus_sex")
    ]
    assert primary[primary["metric"] == "profile_hellinger"]["cell_type"].nunique() == 8
    assert "vein endothelial cell" not in set(primary["cell_type"])
    for row in primary.itertuples():
        expected = donors[donors["cell_type"] == row.cell_type]
        assert len(expected) == row.n_donors
        assert expected["n_cells"].sum() == row.n_cells
