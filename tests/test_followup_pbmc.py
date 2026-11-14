import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm

from src.followup_pbmc import (
    clr_transform, common_cell_type, fit_vectorized, hellinger_distances,
    random_signature_scores,
)

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
PREFIX = "followup-pbmc"


def test_clr_and_released_label_harmonization():
    counts = pd.DataFrame([[1, 2, 0], [3, 4, 5]], columns=["a", "b", "c"])
    transformed = clr_transform(counts)
    np.testing.assert_allclose(transformed.sum(axis=1), 0, atol=1e-12)
    assert common_cell_type("CD8-Naive") == "CD8 naive T"
    assert common_cell_type("CD8 Naive") == "CD8 naive T"
    assert common_cell_type("CD16-MC") == "non-classical monocytes"
    assert common_cell_type("CD16+ Monocytes") == "non-classical monocytes"


def test_vectorized_donor_model_matches_ols():
    frame = pd.DataFrame({
        "old": [0, 0, 0, 1, 1, 1], "age": [25, 30, 35, 61, 70, 80],
        "sex": ["male", "female"] * 3, "mean_log_umi": [1.0] * 6,
        "qc_mean_log_umi": [1.0] * 6, "n_cells": [20] * 6,
    })
    outcomes = np.array([[1, 2, 4, 3, 5, 8], [0, 1, 1, 2, 2, 4]], dtype=float)
    fitted = fit_vectorized(outcomes, frame, "strict_group", "primary")
    X = sm.add_constant(frame.assign(male=(frame["sex"] == "male").astype(float))[["old", "male"]].astype(float))
    expected = sm.OLS(outcomes[0], X).fit()
    np.testing.assert_allclose(fitted["effect"][0], expected.params.iloc[1])
    np.testing.assert_allclose(fitted["se"][0], expected.bse.iloc[1])
    np.testing.assert_allclose(fitted["pvalue"][0], expected.pvalues.iloc[1])


def test_matched_null_and_hellinger_are_deterministic():
    matrix = np.arange(1, 601, dtype=float).reshape(100, 6)
    up, down = np.array([0, 1]), np.array([2, 3])
    first = random_signature_scores(matrix, up, down, 5, np.random.default_rng(42))
    second = random_signature_scores(matrix, up, down, 5, np.random.default_rng(42))
    np.testing.assert_allclose(first, second)
    counts = sp.csr_matrix(np.array([[4, 0], [1, 1], [0, 4]], dtype=float))
    observed = hellinger_distances(counts, np.array([0.5, 0.5]), np.array([4, 2, 4], dtype=float))
    np.testing.assert_allclose(observed, [np.sqrt(1 - np.sqrt(0.5)), 0, np.sqrt(1 - np.sqrt(0.5))])


def test_raw_validation_and_genome_restriction():
    zhu = json.loads((RESULTS / f"{PREFIX}-GSE213516-validation.json").read_text())
    thomson = json.loads((RESULTS / f"{PREFIX}-GSE214546-validation.json").read_text())
    assert len(zhu["donors"]) == 17
    assert len(thomson["donors"]) == 32
    for verification in [zhu, thomson]:
        for donor in verification["donors"]:
            assert donor["all_stored_counts_finite_nonnegative_integer"]
            assert donor["barcode_join_one_to_one"]
            assert donor["annotation_rows_matched"] > 0
    assert sum(d["unmatched_annotation_rows_excluded"] for d in zhu["donors"]) == 199
    assert zhu["genome_reference_adjudication"]["verdict"] == "unresolved_from_released_features"
    assert "symbol-level only" in zhu["genome_reference_adjudication"]["restriction"]


def test_strict_and_sensitivity_result_contract():
    p1 = pd.read_csv(RESULTS / f"{PREFIX}-p1-models.csv")
    strict_p1 = p1[(p1["design"] == "strict_group") & (p1["model"] == "primary")]
    assert strict_p1.groupby("cohort").size().to_dict() == {"GSE213516": 23, "GSE214546": 13}
    targeted = pd.read_csv(RESULTS / f"{PREFIX}-pseudobulk-targeted.csv")
    strict_targeted = targeted[(targeted["design"] == "strict_group") & (targeted["model"] == "primary")]
    assert len(strict_targeted) == 12
    assert strict_targeted.groupby("cohort")["n_donors"].first().to_dict() == {"GSE213516": 10, "GSE214546": 23}
    atlas = targeted[targeted["design"] == "atlas_group"]
    assert atlas.groupby("cohort")["n_donors"].first().to_dict() == {"GSE213516": 12, "GSE214546": 25}
    exact = targeted[targeted["design"] == "exact_age"]
    assert exact.groupby("cohort")["n_donors"].first().to_dict() == {"GSE213516": 17, "GSE214546": 32}
    composition = pd.read_csv(RESULTS / f"{PREFIX}-clr-targeted.csv")
    strict_composition = composition[(composition["design"] == "strict_group") & (composition["model"] == "primary")]
    cd8 = strict_composition[strict_composition["common_cell_type"] == "CD8 naive T"]
    assert len(cd8) == 2 and (cd8["effect"] < 0).all() and (cd8["q_bh"] < 0.05).all()
    mait = strict_composition[strict_composition["common_cell_type"] == "MAIT"]
    assert len(mait) == 2 and mait["status"].eq("untestable").all()
    hellinger = pd.read_csv(RESULTS / f"{PREFIX}-hellinger-models.csv")
    strict_h = hellinger[(hellinger["design"] == "strict_group") & (hellinger["model"] == "primary")]
    assert strict_h.groupby("cohort").size().to_dict() == {"GSE213516": 22, "GSE214546": 13}
    fitted = strict_h[strict_h["status"] == "fit"]
    for column in ["effect", "se", "ci95_low", "ci95_high", "pvalue", "q_bh"]:
        assert np.isfinite(fitted[column]).all()


def test_exploratory_families_resources_and_report():
    expected = {"GSE213516": 344557, "GSE214546": 232716}
    for cohort, rows in expected.items():
        frame = pd.read_csv(RESULTS / f"{PREFIX}-{cohort}-pseudobulk-exploratory.csv.gz",
                            usecols=["design", "model", "q_bh"])
        strict = frame[(frame["design"] == "strict_group") & (frame["model"] == "primary")]
        assert len(strict) == rows
        assert strict["q_bh"].between(0, 1).all()
    resources = json.loads((RESULTS / f"{PREFIX}-resource-use.json").read_text())
    for cohort in expected:
        assert resources[cohort]["guard_memory_gb"] == 8
        assert resources[cohort]["guard_min_free_gb"] == 30
        assert resources[cohort]["parallel_heavy_jobs"] == 1
        assert resources[cohort]["max_rss_kib"] < 8 * 1024**2
    for name in ["targeted-forest", "p1-forest", "hellinger-forest"]:
        assert (REPO / "figures" / f"{PREFIX}-{name}.png").stat().st_size > 10000
    report = (REPO / "reports" / "followup-pbmc-analysis.md").read_text()
    assert "| nan |" not in report
    assert "untestable" in report and "symbol-level only" in report



def test_thomson_robustness_is_explicit_in_report():
    p1 = pd.read_csv(RESULTS / f"{PREFIX}-GSE214546-p1-models.csv")
    strict = p1[p1["design"] == "strict_group"]
    survivors = strict.groupby("model")["q_bh"].agg(lambda x: int((x < 0.05).sum())).to_dict()
    assert survivors == {"primary": 11, "depth": 8, "cmv": 11, "pool": 9, "batch": 0}
    assert np.isclose(strict.loc[strict["model"] == "batch", "q_bh"].min(), 0.1094438352051024)
    clr = pd.read_csv(RESULTS / f"{PREFIX}-GSE214546-clr-models.csv")
    batch = clr[(clr["design"] == "strict_group") & (clr["model"] == "batch")
                & clr["cell_type"].isin(["CD8 Naive", "CD16+ Monocytes"])]
    assert len(batch) == 2 and np.allclose(batch["q_bh"], 0.139044)
    hell = pd.read_csv(RESULTS / f"{PREFIX}-GSE214546-hellinger-models.csv")
    strict_h = hell[hell["design"] == "strict_group"]
    assert (strict_h.loc[strict_h["model"] == "batch_depth", "q_bh"] < 0.05).sum() == 2
    assert (strict_h.loc[strict_h["model"] == "cmv_depth", "q_bh"] < 0.05).sum() == 1
    assert (strict_h.loc[strict_h["model"] == "primary", "q_bh"] < 0.05).sum() == 0
    report = (REPO / "reports" / "followup-pbmc-analysis.md").read_text()
    assert "minimum batch q=0.109444" in report
    assert "Broad P1 is not robust to batch" in report
    assert "-1.396839" in report and "+1.040266" in report
    assert "exploratory all-label batch-family q=0.139044" in report
    assert "exploratory batch-depth family has 2 survivors" in report
    assert "exploratory CMV-depth family has 1" in report
