import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from scripts.run_followup_lake_kidney import (
    MIN_CELLS,
    MIN_DONORS,
    NULL_ITERATIONS,
    NULL_SEED,
    clr_transform,
    eligibility_table,
    expected_ages,
    expected_healthy,
    fit_donor_model,
    hellinger_distances,
    parse_exact_age,
    random_signature_scores,
)

REPO = Path(__file__).resolve().parent.parent


def test_frozen_donor_sets_and_primary_age_rule():
    assert len(expected_healthy("scCv3")) == 18
    assert len(expected_ages("scCv3")) == 18
    assert len(expected_healthy("snCv3")) == 13
    assert len(expected_ages("snCv3")) == 10
    assert expected_ages("snCv3")["3490"] == 14
    primary_sn = [age for age in expected_ages("snCv3").values() if 20 <= age <= 120]
    assert len(primary_sn) == 9
    assert 14 not in primary_sn
    plan = json.loads((REPO / "data/manifest/lake_kidney_analysis_plan.json").read_text())
    assert plan["endpoints"]["P1_transport"]["matched_random_null"].startswith("200 draws")
    assert plan["endpoints"]["donor_pseudobulk_gene_models"]["status"] == "not_predeclared_not_performed"


def test_exact_age_parser_never_uses_decade_midpoints():
    assert parse_exact_age("45-year-old stage") == 45
    assert np.isnan(parse_exact_age("fifth decade stage"))
    assert np.isnan(parse_exact_age("14-year-old"))
    assert np.isnan(parse_exact_age("45"))


def test_clr_zero_handling_and_closure():
    counts = pd.DataFrame({"a": [0, 2, 0], "b": [3, 0, 0], "c": [0, 0, 4]})
    transformed = clr_transform(counts)
    assert np.isfinite(transformed.to_numpy()).all()
    np.testing.assert_allclose(transformed.sum(axis=1), 0.0, atol=1e-12)
    assert transformed.loc[0, "a"] < 0
    assert transformed.loc[2, "c"] > 0


def test_cell_type_eligibility_requires_six_donors_with_twenty_cells():
    rows = []
    for donor in range(6):
        for cell_type in ["pass", "fail"]:
            for _ in range(MIN_CELLS if cell_type == "pass" or donor < 5 else MIN_CELLS - 1):
                rows.append((f"d{donor}", cell_type))
    obs = pd.DataFrame(rows, columns=["donor_id", "cell_type"])
    obs["disease"] = "normal"
    obs["age"] = 50
    eligible = eligibility_table(obs).set_index("cell_type")["eligible"]
    assert bool(eligible["pass"])
    assert not bool(eligible["fail"])
    assert MIN_DONORS == 6


def test_donor_model_is_full_rank_and_reports_effect_interval():
    donors = pd.DataFrame({
        "age": np.arange(30, 60, 5), "sex": ["female", "male"] * 3,
        "n_cells": [20] * 6, "mean_log_umi": np.linspace(7, 8, 6),
    })
    donors["score"] = 0.01 * donors["age"] + 0.1 * (donors["sex"] == "male")
    fit = fit_donor_model(donors, "score")
    assert fit["n_donors"] == 6
    assert np.isclose(fit["slope_per_year"], 0.01)
    assert fit["ci95_low"] <= fit["slope_per_year"] <= fit["ci95_high"]
    assert np.isclose(fit["effect_per_30y"], 0.3)


def test_matched_random_null_is_reproducible():
    rng = np.random.default_rng(4)
    matrix = rng.normal(size=(60, 8))
    up = np.array([0, 1, 2])
    down = np.array([3, 4])
    first = random_signature_scores(matrix, up, down, 5, np.random.default_rng(NULL_SEED))
    second = random_signature_scores(matrix, up, down, 5, np.random.default_rng(NULL_SEED))
    np.testing.assert_allclose(first, second)
    assert first.shape == (5, 8)
    assert NULL_ITERATIONS == 200


def test_hellinger_distance_zero_for_identical_profiles():
    counts = sp.csr_matrix([[1, 1, 2], [2, 2, 4]], dtype=float)
    reference = np.array([0.25, 0.25, 0.5])
    totals = np.asarray(counts.sum(axis=1)).ravel()
    distances = hellinger_distances(counts, reference, totals)
    np.testing.assert_allclose(distances, 0.0, atol=1e-12)


def test_integer_aggregate_preserves_raw_count_dtype():
    from scripts.run_followup_lake_kidney import aggregate_block

    matrix = sp.csr_matrix([[1, 0, 2], [3, 4, 0]], dtype=np.float64)
    output = np.zeros((1, 3), dtype=np.int64)
    aggregate_block(matrix, np.array([0, 0]), output)
    np.testing.assert_array_equal(output, [[4, 4, 2]])
    assert output.dtype == np.int64


def test_delivered_csvs_preserve_assays_donors_and_multiplicity_families():
    for assay, healthy, exact, adults in [("scCv3", 18, 18, 18), ("snCv3", 13, 10, 9)]:
        verification = json.loads((REPO / f"results/followup-lake-kidney-{assay}-verification.json").read_text())
        assert verification["raw_X_all_stored_values_finite_nonnegative_integer"]
        assert verification["healthy_donors"] == healthy
        assert verification["exact_age_healthy_donors"] == exact
        assert verification["primary_adult_donors"] == adults
        assert verification["n_obs"] > 0

    donors = pd.read_csv(REPO / "results/followup-lake-kidney-donor-crosswalk.csv")
    assert not donors.duplicated(["assay_release", "cellxgene_donor_id"]).any()
    primary_sn = donors[(donors["assay_release"] == "snCv3") & donors["primary_adult"]]
    assert len(primary_sn) == 9
    assert primary_sn["age"].min() >= 20
    assert "3490" not in set(primary_sn["cellxgene_donor_id"].astype(str))

    p1 = pd.read_csv(REPO / "results/followup-lake-kidney-p1-models.csv")
    assert (p1["n_donors"] >= MIN_DONORS).all()
    primary = p1[p1["model"] == "primary_age_sex"]
    assert primary.groupby(["assay_release", "donor_set", "cell_type"]).size().eq(1).all()
    null = primary[primary["donor_set"] == "primary_adult"]
    assert null["random_iterations"].eq(NULL_ITERATIONS).all()
    assert null["random_seed"].eq(NULL_SEED).all()
    assert null["random_pvalue"].between(0, 1).all()
    assert null["random_q_bh"].between(0, 1).all()
    assert (null["random_pvalue"] >= 1 / (NULL_ITERATIONS + 1) - 1e-12).all()

    composition = pd.read_csv(REPO / "results/followup-lake-kidney-composition-models.csv")
    assert (composition["n_donors"] >= MIN_DONORS).all()
    assert composition["q_bh"].between(0, 1).all()
    counts = pd.read_csv(REPO / "results/followup-lake-kidney-composition-donor-counts.csv")
    assert not counts.duplicated(["assay_release", "donor_id", "cell_type"]).any()
    assert (counts["n_cells"] >= 0).all()
    assert (counts["n_cells"] == 0).any()


def test_hellinger_outputs_respect_floors_and_untestable_labels():
    eligibility = pd.read_csv(REPO / "results/followup-lake-kidney-hellinger-eligibility.csv")
    sn = eligibility[eligibility["assay_release"] == "snCv3"]
    assert not sn["hellinger_eligible"].any()
    model_path = REPO / "results/followup-lake-kidney-hellinger-models.csv"
    if model_path.exists():
        models = pd.read_csv(model_path)
        assert (models["n_donors"] >= MIN_DONORS).all()
        assert models["q_bh"].between(0, 1).all()
        donors = pd.read_csv(REPO / "results/followup-lake-kidney-hellinger-donors.csv")
        assert (donors["n_cells"] >= MIN_CELLS).all()
        assert (donors["cross_fitted_reference_cells"] >= 50).all()
        assert not donors.duplicated(["assay_release", "cell_type", "donor_id"]).any()
