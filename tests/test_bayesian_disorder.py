import ast

import pandas as pd


def test_bayesian_primary_posteriors_are_complete():
    posterior = pd.read_csv("results/h3-bayesian-posterior.csv")
    primary = posterior[
        (posterior["analysis"] == "all_cells_all_donor_ids")
        & (posterior["prior_sd"] == 0.5)
        & (posterior["likelihood"] == "student_t")
    ]
    assert len(primary[primary["scope"] == "global"]) == 2
    assert len(primary[primary["scope"] == "cell_type"]) == 17
    assert primary["prob_positive"].between(0, 1).all()
    assert (primary["hdi94_low_sd15"] <= primary["hdi94_high_sd15"]).all()


def test_bayesian_prior_and_sampling_sensitivities_exist():
    posterior = pd.read_csv("results/h3-bayesian-posterior.csv")
    assert set(posterior["prior_sd"]) == {0.25, 0.5, 1.0}
    assert set(posterior["analysis"]) == {
        "all_cells_all_donor_ids",
        "all_cells_exclude_pooled_ids",
        "balanced_20_cells_all_donor_ids",
    }
    prior = pd.read_csv("results/h3-bayesian-prior-predictive.csv")
    assert prior["observed_within_prior_998"].all()


def test_bayesian_convergence_and_diagnostics_are_recorded():
    diagnostics = pd.read_csv("results/h3-bayesian-diagnostics.csv")
    convergence = diagnostics["convergence"].apply(ast.literal_eval)
    assert convergence.apply(lambda value: value["all_ok"]).all()
    assert (diagnostics["likelihood"] == "student_t").all()
