import numpy as np
import pandas as pd
import statsmodels.api as sm

from scripts.run_followup_ts2_blood_nk import (
    DONORS, SEED, design_matrix, fit_age_models, prepare_p1, run_p1,
    young_reference_floors,
)


def test_vectorized_age_models_match_statsmodels():
    rng = np.random.default_rng(19)
    donors = pd.DataFrame({
        "age_years": [59, 61, 69, 33, 59, 45, 60],
        "male": [0, 0, 0, 1, 1, 0, 1],
        "assay_5prime_present": [0, 1, 0, 0, 1, 1, 1],
        "mean_log_umi": [8.94, 8.78, 9.07, 8.90, 8.79, 8.76, 8.41],
    })
    design = design_matrix(donors)
    outcomes = rng.normal(size=(7, 2))
    result = fit_age_models(outcomes, design)
    fit = sm.OLS(outcomes[:, 0], design).fit()
    np.testing.assert_allclose(result["effect_age_10y"][0], fit.params[1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result["se"][0], fit.bse[1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result["pvalue"][0], fit.pvalues[1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result["ci95_low"][0], fit.conf_int()[1, 0], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result["ci95_high"][0], fit.conf_int()[1, 1], rtol=1e-12, atol=1e-12)
    assert result["residual_df"][0] == 2


def test_young_reference_requires_two_donors_and_crossfit_floor():
    nk = pd.DataFrame({
        "donor": DONORS,
        "age": [59, 61, 69, 33, 59, 45, 60],
    })
    assert not young_reference_floors(nk)["passes"]
    expanded = nk.loc[nk.index.repeat(60)].reset_index(drop=True)
    expanded.loc[expanded["donor"] == "TSP14", "age"] = 35
    floors = young_reference_floors(expanded)
    assert floors["n_young_donors"] == 2
    assert floors["minimum_remaining_cells"] == 60
    assert floors["passes"]


def synthetic_expression():
    from src.aging_programs import load_program_tiers
    from pathlib import Path
    program = load_program_tiers(Path(__file__).resolve().parent.parent, "human")["P1_Clean"]
    genes = program["up"] + program["down"] + [f"NULL{i}" for i in range(2000)]
    rng = np.random.default_rng(72)
    matrix = rng.gamma(shape=1.2, scale=1.0, size=(7, len(genes)))
    return {
        "donors": np.asarray(DONORS),
        "gene_ids": np.asarray(genes),
        "gene_names": np.asarray(genes),
        "pseudobulk": matrix,
    }


def test_p1_null_is_reproducible_and_uses_two_residual_df():
    expression = synthetic_expression()
    setup = prepare_p1(expression)
    assert setup["valid"]
    assert len(setup["detected"]["up"]) == 76
    assert len(setup["detected"]["down"]) == 42
    donors = pd.DataFrame({
        "age_years": [59, 61, 69, 33, 59, 45, 60],
        "male": [0, 0, 0, 1, 1, 0, 1],
        "assay_5prime_present": [0, 1, 0, 0, 1, 1, 1],
        "mean_log_umi": [8.94, 8.78, 9.07, 8.90, 8.79, 8.76, 8.41],
    })
    design = design_matrix(donors)
    observed1, model1, null1 = run_p1(expression, design, setup)
    observed2, model2, null2 = run_p1(expression, design, setup)
    assert observed1.shape == (7,)
    assert model1["residual_df"] == 2
    assert model1["null_seed"] == SEED
    assert 1 / 201 <= model1["empirical_pvalue"] <= 1
    np.testing.assert_array_equal(null1["effect_age_10y"], null2["effect_age_10y"])
    assert model1["effect_age_10y"] == model2["effect_age_10y"]


def test_committed_outputs_preserve_result_family_and_p1_null():
    genes = pd.read_csv("results/followup-ts2-blood-nk-gene-models.csv")
    assert len(genes) == 32830
    assert genes["ensembl_id"].is_unique
    assert genes["q_bh"].notna().all()
    assert not (genes["q_bh"] < 0.05).any()
    model = pd.read_csv("results/followup-ts2-blood-nk-p1-model.csv").iloc[0]
    null = pd.read_csv("results/followup-ts2-blood-nk-p1-null.csv")
    assert len(null) == 200
    assert np.isclose(model["effect_age_10y"], 0.03172447652517948)
    assert np.isclose(model["ci95_low"], -0.13387231240683528)
    assert np.isclose(model["ci95_high"], 0.19732126545719425)
    centered = null["effect_age_10y"] - null["effect_age_10y"].mean()
    observed = model["effect_age_10y"] - null["effect_age_10y"].mean()
    empirical = (1 + (centered.abs() >= abs(observed)).sum()) / 201
    assert np.isclose(empirical, model["empirical_pvalue"])
    endpoints = pd.read_csv("results/followup-ts2-blood-nk-endpoint-status.csv").set_index("endpoint")
    assert endpoints.loc["donor_pseudobulk_gene_models", "status"] == "completed"
    assert endpoints.loc["P1_transport", "status"] == "completed"
    assert endpoints.loc["off_program_expression", "status"] == "untestable"
    assert endpoints.loc["Hellinger_displacement", "status"] == "untestable"
    assert not endpoints.loc["off_program_expression", "outcomes_examined"]
    assert not endpoints.loc["Hellinger_displacement", "outcomes_examined"]
