import numpy as np
import pandas as pd
import statsmodels.api as sm

from src.reviewer_stats import pairwise_pearson, random_effects_meta, vectorized_ols


def test_vectorized_ols_matches_statsmodels():
    rng = np.random.default_rng(42)
    X = np.column_stack([np.ones(12), np.arange(12), rng.normal(size=12)])
    Y = rng.normal(size=(12, 5))
    slopes, se, pvalues = vectorized_ols(Y, X)
    for index in range(Y.shape[1]):
        fit = sm.OLS(Y[:, index], X).fit()
        assert np.isclose(slopes[index], fit.params[1])
        assert np.isclose(se[index], fit.bse[1])
        assert np.isclose(pvalues[index], fit.pvalues[1])


def test_pairwise_pearson_excludes_diagonal():
    effects = np.array([[1, 2, 3], [1, 2, 3], [3, 2, 1]])
    np.testing.assert_allclose(pairwise_pearson(effects), [1, -1, -1])


def test_random_effects_handles_zero_heterogeneity():
    result = random_effects_meta([0.1, 0.2, 0.15], [0.05, 0.05, 0.05])
    assert result["tau2"] == 0
    assert np.isclose(result["pooled_effect"], 0.15)


def test_composition_global_correction_has_expected_scope():
    result = pd.read_csv("results/reviewer-phase1-composition-global-bh.csv")
    assert result["q_global_bh"].notna().all()
    summary = pd.read_csv("results/reviewer-phase1-composition-global-bh-by-organ.csv")
    mouse = summary[summary["species"] == "Mouse"]
    human = summary[summary["species"] == "Human"]
    assert int((mouse["n_significant_global"] > 0).sum()) == 7
    assert int((human["n_significant_global"] > 0).sum()) == 2


def test_random_signature_score_subtracts_down_genes():
    from scripts.run_reviewer_phase1_signature_null import random_signature_scores

    matrix = np.array([[1.0, 1.0], [3.0, 3.0], [5.0, 5.0], [7.0, 7.0]])
    scores = random_signature_scores(
        {"matrix": matrix, "genes": np.array(["a", "b", "c", "d"])},
        np.array([0]), np.array([1]), 1, np.random.default_rng(42), np.ones(4, dtype=bool),
    )
    np.testing.assert_allclose(scores, [[-2.0, -2.0]])


def test_age_permutation_preserves_sex_strata():
    from scripts.run_reviewer_phase1_gene_concordance import permuted_ages

    metadata = pd.DataFrame({"age": [1, 2, 3, 4, 5, 6, 7, 8],
                             "sex_male": [0, 0, 0, 0, 1, 1, 1, 1]})
    permuted = permuted_ages(metadata, np.random.default_rng(42))
    for sex in (0, 1):
        expected = metadata.loc[metadata["sex_male"] == sex, "age"]
        observed = permuted[metadata["sex_male"].to_numpy() == sex]
        assert set(observed) == set(expected)
    assert not np.array_equal(permuted, metadata["age"].to_numpy())


def test_signature_null_preserves_observed_slopes():
    null = pd.read_csv("results/reviewer-phase1-signature-random-null.csv")
    observed = pd.read_csv("results/m10-celltype-regressions.csv")
    merged = null.merge(observed[["cell_type", "slope_covar_adj"]], on="cell_type")
    np.testing.assert_allclose(merged["actual_slope"], merged["slope_covar_adj"])
