import numpy as np
import scipy.sparse as sp

from scripts.run_followup_de_man_entropy import (
    binomial_downsample, kl_uniform, multinomial_downsample, shannon)


def test_entropy_functions_match_known_values():
    uniform = sp.csr_matrix(np.array([[1, 1], [1, 0]], dtype=float))
    assert np.isclose(shannon(uniform)[0], np.log(2))
    assert np.isclose(kl_uniform(uniform)[0], 0)
    assert np.isclose(shannon(uniform)[1], 0)
    assert np.isclose(kl_uniform(uniform)[1], np.log(2))


def test_downsampling_preserves_expected_library_totals():
    matrix = sp.csr_matrix(np.array([[4, 4], [2, 2], [1, 1]], dtype=np.int64))
    thinned = binomial_downsample(matrix, np.random.default_rng(4))
    totals = np.asarray(thinned.sum(axis=1)).ravel()
    assert np.all(totals >= 2)
    assert np.any(totals == 2)
    fixed = multinomial_downsample(matrix, 2, np.random.default_rng(4))
    assert np.all(np.asarray(fixed.sum(axis=1)).ravel() == 2)
