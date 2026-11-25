"""Tests for the mature concordance hierarchy driver.

Covers the safe (QR / Frisch-Waugh-Lovell) OLS, the jitter-UMI clean design,
within-unit BH counting, cross-pair construction and null, the 311-row NA
parity for unidentifiable units, and the human decomposition helpers.
All data here is synthetic; no real biological data is used.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import run_mature_concordance_hierarchy as hier  # noqa: E402
from scripts import run_rev2_concordance as rev2  # noqa: E402


def _metadata(n=12, ages=None, sex=None, umi_ptp=1.0):
    rng = np.random.default_rng(7)
    if ages is None:
        ages = np.linspace(30, 80, n)
    if sex is None:
        sex = np.array([0, 1] * (n // 2) + [0] * (n % 2), dtype=float)
    umi = 7.0 + umi_ptp * rng.random(n)
    return pd.DataFrame({
        "donor": [f"D{i}" for i in range(n)],
        "age": np.asarray(ages, dtype=float),
        "sex_male": np.asarray(sex, dtype=float),
        "mean_log_umi": umi,
        "n_cells": np.full(n, 100),
    })


def _expression(n=12, n_genes=50, seed=3):
    rng = np.random.default_rng(seed)
    return rng.random((n, n_genes)) * 5.0 + 0.1


def test_safe_ols_matches_reference_on_well_conditioned_design():
    metadata = _metadata()
    design = hier.clean_unit_design(metadata)
    expression = _expression()
    slopes, se, pvalues, ok = hier.safe_vectorized_ols(expression, design)
    assert ok
    # Reference: statsmodels OLS per gene for the first five genes.
    import statsmodels.api as sm

    for gene in range(5):
        fit = sm.OLS(expression[:, gene], design).fit()
        assert slopes[gene] == pytest.approx(fit.params[1], abs=1e-8)
        assert se[gene] == pytest.approx(fit.bse[1], abs=1e-8)
        assert pvalues[gene] == pytest.approx(fit.pvalues[1], abs=1e-8)
    # And it must agree with the pipeline pinv OLS on a well-conditioned design.
    ref_slopes, ref_se, _ = rev2.vectorized_ols(expression, design)
    assert np.allclose(slopes, ref_slopes, rtol=1e-8, atol=1e-10)
    assert np.allclose(se, ref_se, rtol=1e-6, atol=1e-10)


def test_safe_ols_single_age_is_unidentifiable():
    metadata = _metadata(ages=np.full(12, 50.0))
    design = hier.clean_unit_design(metadata)
    slopes, se, pvalues, ok = hier.safe_vectorized_ols(_expression(), design)
    assert not ok
    assert np.isnan(slopes).all()
    assert (pvalues == 1.0).all()


def test_safe_ols_age_sex_collinear_is_unidentifiable():
    # Age perfectly collinear with sex: old donors are all male, young female.
    ages = np.array([30, 30, 30, 30, 30, 30, 80, 80, 80, 80, 80, 80], dtype=float)
    sex = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1], dtype=float)
    metadata = _metadata(ages=ages, sex=sex, umi_ptp=0.0)
    design = hier.clean_unit_design(metadata)
    _, _, _, ok = hier.safe_vectorized_ols(_expression(), design)
    assert not ok


def test_safe_ols_zero_residual_and_constant_genes_are_not_hits():
    metadata = _metadata(n=8)
    design = hier.clean_unit_design(metadata)
    expression = _expression(n=8, n_genes=6, seed=11)
    # Gene 0: an exact linear combination of the design columns (perfect fit).
    expression[:, 0] = 2.0 * design[:, 0] + 0.5 * design[:, 1] - design[:, 2]
    # Gene 1: constant across donors.
    expression[:, 1] = 3.0
    _, se, pvalues, ok = hier.safe_vectorized_ols(expression, design)
    assert ok
    assert se[0] == 0.0 and pvalues[0] == 1.0
    assert se[1] == 0.0 and pvalues[1] == 1.0
    assert np.isfinite(pvalues).all()


def test_clean_design_drops_jitter_umi_keeps_genuine_umi():
    jitter = hier.clean_unit_design(_metadata(umi_ptp=1e-6))
    assert jitter.shape[1] == 3  # const, age, sex
    genuine = hier.clean_unit_design(_metadata(umi_ptp=1.0))
    assert genuine.shape[1] == 4  # const, age, sex, umi
    # Column order and values match the pipeline design when UMI is genuine.
    pipeline = rev2.unit_design(_metadata(umi_ptp=1.0))
    assert np.array_equal(genuine, pipeline)


def test_bh_count_is_within_unit_not_pooled():
    # Unit A: one p=0.01 among 99 p=0.5. Within-unit BH kills it
    # (q = 0.01*100/1 = 1.0, then monotonicity caps it at 0.5): 0 hits.
    count_a, _ = hier.bh_count(np.array([0.01] + [0.5] * 99))
    assert count_a == 0
    # Unit B: one hundred p=0.001. Within-unit BH: 100 hits.
    count_b, _ = hier.bh_count(np.array([0.001] * 100))
    assert count_b == 100
    # Pooled over all 200, the 0.01 survives (q = 0.01*200/101 = 0.0198):
    # 101 hits. So within-unit counting is not pooled counting.
    pooled, _ = hier.bh_count(np.array([0.01] + [0.5] * 99 + [0.001] * 100))
    assert pooled == 101
    assert pooled != count_a + count_b


def test_cross_pair_mask():
    mask = hier.cross_pair_mask(["A", "A", "B", "B"])
    # triu order: (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
    assert mask.tolist() == [False, True, True, True, True, False]
    assert int(mask.sum()) == 4


def test_fly_pair_id():
    assert hier.fly_pair_id("S1_head_A") == "S1_A"
    assert hier.fly_pair_id("S1_body_A") == "S1_A"
    with pytest.raises(ValueError):
        hier.fly_pair_id("S1_A")


def test_technical_gene_names():
    genes = ["RpL3", "mt:CoI", "MT-ND1", "MRPL12", "GAPDH", "Rps6", "Act5C"]
    assert hier.technical_gene_names(genes).tolist() == [True, True, True, True, False, True, False]


def _cross_records():
    rng = np.random.default_rng(5)
    records = []
    for cohort in ("C1", "C2"):
        for unit_no in range(2):
            metadata = pd.DataFrame({
                "donor": [f"{cohort}_D{i}" for i in range(6)],
                "age": [30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
                "sex_male": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                "mean_log_umi": 7.0 + rng.random(6),
                "n_cells": [100] * 6,
            })
            records.append({
                "species": "Human", "tissue": "T", "cohort": cohort,
                "cell_type": f"cell_{cohort}_{unit_no}", "metadata": metadata,
                "expression": rng.random((6, 20)) * 5.0 + 0.1,
            })
    return records


def test_cross_null_synthetic_pair_construction_and_null():
    records = _cross_records()
    labels = np.array([r["cohort"] for r in records])
    cross = hier.cross_pair_mask(labels)
    assert int(cross.sum()) == 4
    masks = {"primary": np.ones(20, dtype=bool)}
    result = hier.cross_null(records, masks, cross, "independent_by_cohort",
                             permutations=10, seed=1)
    assert np.isfinite(result["observed"]["primary"])
    assert len(result["nulls"]["primary"]) == 10
    summary = hier.summarize_null(result["observed"]["primary"], result["nulls"]["primary"])
    assert 0.0 <= summary["joint_p_greater"] <= 1.0
    assert result["null_pairs"].shape == (10, 6)  # 4 units -> 6 triu pairs
    assert result["masked_iterations"] == 0
    assert len(result["observed_values"]["primary"]) == 4


def test_level0_table_keeps_unidentifiable_units_as_na_rows():
    good_md = _metadata(n=8)
    bad_md = _metadata(n=8, ages=np.full(8, 50.0))
    genes = np.array([f"g{i}" for i in range(10)])
    good = {"cell_type": "good_cell", "metadata": good_md,
            "expression": _expression(n=8, n_genes=10, seed=21)}
    bad = {"cell_type": "bad_cell", "metadata": bad_md,
           "expression": _expression(n=8, n_genes=10, seed=22)}
    units = [good, bad]
    prepared = [{
        "tissue": "T", "cohort": "C", "genes": genes, "units": units,
        "fits": [hier.unit_ols(u) for u in units],
        "identifiable": [hier.age_identifiable(u) for u in units],
        "stratum_primary": hier.primary_mask(units),
    }]
    table, excluded = hier.level0_table("Human", prepared, None)
    assert set(table["gene_set"].unique()) <= {"stratum_primary", "unfiltered", "species_primary"}
    primary = table[table["pipeline_gene_set"]]
    assert len(primary) == 2
    bad_row = primary[primary["cell_type"] == "bad_cell"].iloc[0]
    assert pd.isna(bad_row["n_significant_q05"])
    assert pd.isna(bad_row["frac_significant"])
    assert bad_row["n_genes_tested"] == int(prepared[0]["stratum_primary"].sum())
    assert len(excluded) == 1
    assert excluded.iloc[0]["cell_type"] == "bad_cell"


def test_human_decomposition_helpers():
    observed = np.array([0.1, 0.2, 0.3, 0.4])
    rng = np.random.default_rng(0)
    null = rng.normal(0.0, 0.1, size=(100, 4))
    subset = hier.subset_null_summary(observed, null, np.array([True, True, False, False]))
    assert subset["observed_mean_pearson"] == pytest.approx(0.15)
    assert 0.0 <= subset["joint_p_greater"] <= 1.0
    grouped = hier.group_weighted_null_summary(observed, null, np.array(["A", "A", "B", "B"]))
    # Equal weighting: (mean(A) + mean(B)) / 2 = (0.15 + 0.35) / 2 = 0.25.
    assert grouped["observed_mean_pearson"] == pytest.approx(0.25)


def test_group_weighted_null_summary_accepts_tuple_labels():
    # Regression: the real human decompositions pass cohort-pair tuples
    # (e.g. ("HBCC_DLPFC", "MSSM_DLPFC")); pd.Series misreads each tuple
    # as a 2D row and raises ValueError. Tuple and string labels must agree.
    observed = np.array([0.1, 0.2, 0.3, 0.4])
    rng = np.random.default_rng(0)
    null = rng.normal(0.0, 0.1, size=(100, 4))
    tuples = np.array([("A", "B"), ("A", "B"), ("A", "C"), ("A", "C")], dtype=object)
    strings = np.array(["A|B", "A|B", "A|C", "A|C"])
    by_tuple = hier.group_weighted_null_summary(observed, null, tuples)
    by_string = hier.group_weighted_null_summary(observed, null, strings)
    assert by_tuple["observed_mean_pearson"] == pytest.approx(
        by_string["observed_mean_pearson"])
    assert by_tuple["observed_mean_pearson"] == pytest.approx(0.25)
    assert 0.0 <= by_tuple["joint_p_greater"] <= 1.0


def test_human_level2_decompositions_end_to_end():
    # Regression: the real human decompositions failed twice, first on tuple
    # group labels (pd.Series read each tuple as a 2D row) and then on
    # np.array of tuples becoming a 2D array (pair_key.tolist() gave lists,
    # unhashable in set()). This drives the real function on synthetic data.
    rng = np.random.default_rng(2)
    pairs = [
        ("HBCC_DLPFC", "MSSM_DLPFC", "neuron A", "neuron A", 0.20),
        ("HBCC_DLPFC", "Lung_CellRef", "neuron A", "macrophage", 0.10),
        ("HBCC_DLPFC", "ImmuneAging_BAL", "neuron A", "T cell", 0.05),
        ("MSSM_DLPFC", "Lung_CellRef", "neuron A", "macrophage", 0.12),
        ("MSSM_DLPFC", "ImmuneAging_BAL", "neuron A", "T cell", 0.03),
        ("Lung_CellRef", "ImmuneAging_BAL", "macrophage", "T cell", 0.08),
    ]
    pair_frame = pd.DataFrame(pairs, columns=[
        "cohort_1", "cohort_2", "cell_type_1", "cell_type_2", "pearson"])
    null_pairs = rng.normal(0.0, 0.1, size=(10, 6))
    detail_rows, cohort_pair_rows = hier.human_level2_decompositions(
        pair_frame, null_pairs, np.arange(6))
    # 1 equal-weighted + 1 strict + 1 same-organ + 1 same-type +
    # 1 different-type + 1 BAL + 4 leave-one-out = 10 detail rows.
    assert len(detail_rows) == 10
    assert len(cohort_pair_rows) == 6
    for row in cohort_pair_rows:
        # Full cohort names, not single characters from indexing a string.
        assert row["cohort_1"] in {"HBCC_DLPFC", "MSSM_DLPFC", "Lung_CellRef",
                                   "ImmuneAging_BAL"}
        assert row["cohort_2"] in {"HBCC_DLPFC", "MSSM_DLPFC", "Lung_CellRef",
                                   "ImmuneAging_BAL"}
        assert row["n_pairs"] == 1
        assert 0.0 <= row["joint_p_greater"] <= 1.0
    scopes = {r["scope"] for r in detail_rows}
    assert "cross_cohort_equal_cohort_pair_weighted" in scopes
    assert "cross_cohort_bal_involved" in scopes
    assert "cross_cohort_same_organ_same_cell_type" in scopes


def test_mouse_mean_fidelity_tolerance_covers_rebuild_offset():
    # Documented rebuild offset: the from-raw mouse rebuild shifts the
    # within-tissue mean by 0.000667 (Fat/Kidney scrublet doublet calls
    # differ from the original). That must pass the mouse mean tolerance
    # and fail the strict 1e-6 tolerance, while a wrong-gene-set error
    # (primary 0.2105 vs lenient 0.2039, diff 0.0066) must fail both.
    assert hier.FIDELITY_TOLERANCE == 1e-6
    assert hier.MOUSE_MEAN_FIDELITY_TOLERANCE == 1e-3
    assert 0.0006673714130476793 <= hier.MOUSE_MEAN_FIDELITY_TOLERANCE
    assert 0.0006673714130476793 > hier.FIDELITY_TOLERANCE
    assert 0.0066 > hier.MOUSE_MEAN_FIDELITY_TOLERANCE


def test_single_stratum_rows_read_committed_csvs():
    rows = hier.single_stratum_rows()
    assert [r["species"] for r in rows] == ["C. elegans", "Zebrafish", "Naked mole-rat"]
    for row in rows:
        assert np.isfinite(row["observed_mean_pearson"])
        assert 0.0 <= row["joint_p_greater"] <= 1.0
