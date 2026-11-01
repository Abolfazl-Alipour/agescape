import numpy as np
import pandas as pd
import scipy.sparse as sp

from scripts.run_h3_transcriptional_disorder import shannon_entropy


def test_shannon_entropy_known_values():
    matrix = sp.csr_matrix(np.array([[1, 1, 1, 1], [4, 0, 0, 0]], dtype=float))
    np.testing.assert_allclose(shannon_entropy(matrix), [np.log(4), 0.0])


def test_primary_transcriptional_disorder_results_trace_to_donors():
    results = pd.read_csv("results/h3-transcriptional-disorder.csv")
    donors = pd.read_csv("results/h3-transcriptional-disorder-donors.csv")
    primary = results[results["analysis"] == "all_cells_all_donor_ids"]
    assert primary[primary["metric"] == "entropy"]["cell_type"].nunique() == 9
    assert (primary["q_bh"] >= 0.05).all()
    for row in primary.itertuples():
        expected = donors[
            (donors["cell_type"] == row.cell_type)
            & (donors["sampling"] == "all_cells")
        ].dropna(subset=[row.metric])
        assert len(expected) == row.n_donors
        assert expected["n_cells"].sum() == row.n_cells


def test_offprogram_marker_provenance_matches_donor_marker_counts():
    donors = pd.read_csv("results/h3-transcriptional-disorder-donors.csv")
    donors = donors[donors["sampling"] == "all_cells"]
    provenance = pd.read_csv("results/h3-offprogram-marker-provenance.csv")
    expected = provenance.groupby(["cell_type", "excluded_donor"])["n_marker_genes"].sum()
    for row in donors.itertuples():
        if row.n_offprogram_genes:
            assert expected.loc[(row.cell_type, row.donor)] == row.n_offprogram_genes
