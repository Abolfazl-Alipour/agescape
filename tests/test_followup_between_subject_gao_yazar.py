"""Tests for Part A between-subject drift/variance/dispersion and Part B spread in Gao Heart & Yazar Blood."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse as sp

from scripts.run_cross_species_expression_vector import reference_without_donor
from scripts.run_h3_expression_vector import hellinger_distance
from scripts.run_followup_between_subject_gao_yazar import (
    _compute_a3_loo_distances_for_matrix,
    _fit_stratified_global_row,
    _prune_covariate_matrix,
    _run_a1_a3_models_from_donors,
    _summarize_a2_perm_stats,
    apply_slice_bh,
    build_human_gene_projection,
    compute_a2_variance_metrics,
    evaluate_verdicts,
    fit_ols_model_row,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"


def test_gene_filtering_and_duplicate_symbol_merging():
    raw_syms = [
        "MT-CO1",
        "RPL10",
        "RPS3",
        "MRPL12",
        "XIST",
        "ERCC-0001",
        "MYH6",
        "MYH6",
        "GAPDH",
        "ACTB",
    ]
    info = build_human_gene_projection(raw_syms)
    retained = info["audit"]["retained_symbols"]
    assert retained == ["MYH6", "GAPDH", "ACTB"]
    assert info["audit"]["n_genes_raw"] == len(raw_syms)
    assert info["audit"]["n_genes_retained"] == 3
    assert info["audit"]["n_duplicate_collapsed"] == 1

    counts = sp.csr_matrix(
        np.array(
            [
                [100, 50, 20, 10, 5, 8, 11, 29, 7, 3],
                [0, 1, 2, 3, 4, 5, 40, 0, 12, 8],
            ],
            dtype=np.float64,
        )
    )
    projected = (counts @ info["proj"]).toarray()
    assert projected.shape == (2, 3)
    # Columns preserved in first-seen unique symbol order: ["MYH6", "GAPDH", "ACTB"]
    assert np.allclose(projected[0], [40.0, 7.0, 3.0])
    assert np.allclose(projected[1], [40.0, 12.0, 8.0])


def test_a1_hellinger_loo_vs_nonloo_bias():
    ref_entry = {
        "counts": np.array([100.0, 200.0, 300.0], dtype=np.float64),
        "n_cells": 60,
        "donors": {
            "Y1": {"counts": np.array([80.0, 20.0, 20.0]), "n_cells": 30},
            "Y2": {"counts": np.array([20.0, 180.0, 280.0]), "n_cells": 30},
        },
    }
    cells_y1 = sp.csr_matrix(np.array([[4.0, 1.0, 1.0], [3.0, 1.0, 0.0]], dtype=np.float64))
    r_loo = reference_without_donor(ref_entry, "Y1")
    assert r_loo["n_cells"] == 30
    hd_loo = hellinger_distance(cells_y1, r_loo["counts"])
    hd_self = hellinger_distance(cells_y1, ref_entry["counts"])
    # Including Y1 inside its own reference pulls its reference toward Y1, biasing self distance lower than LOO
    assert np.all(hd_self < hd_loo)


def test_a3_loo_hellinger_exact_equivalence():
    rng = np.random.default_rng(123)
    counts = rng.integers(0, 20, size=(35, 8)).astype(np.float64)
    counts[:, 0] += 5.0
    donor_idx = np.array([0] * 15 + [1] * 20, dtype=np.int64)
    lookup = np.array([0, 1], dtype=np.int64)
    n_cells_vec = np.array([15, 20], dtype=np.int64)

    l1_mat = sp.csr_matrix(counts / counts.sum(axis=1, keepdims=True))
    # Pass probability SUMS across all n_r cells (as accumulated in Pass 1 / Pass 2b)
    centroid_sums = np.vstack([l1_mat[:15].sum(axis=0).A1, l1_mat[15:].sum(axis=0).A1])

    out = _compute_a3_loo_distances_for_matrix(
        sp.csr_matrix(counts),
        donor_idx,
        centroid_sums,
        n_cells_vec,
        lookup,
        is_l1_already=False,
    )
    # Exact brute-force LOO Hellinger distance: sqrt(1 - sum_g sqrt(p_ig * p_{-i,g}))
    l1_dense = l1_mat.toarray()
    for d_id, sl in ((0, slice(0, 15)), (1, slice(15, 35))):
        sub = l1_dense[sl]
        exact = []
        for i in range(len(sub)):
            others = np.delete(sub, i, axis=0).mean(axis=0)
            bc = float(np.clip(np.sum(np.sqrt(sub[i] * others)), 0.0, 1.0))
            exact.append(float(np.sqrt(max(0.0, 1.0 - bc))))
        sm, cnt = out[d_id]
        assert cnt == len(exact)
        assert np.isclose(sm / cnt, float(np.mean(exact)), atol=1e-12)


def test_a2_residualized_variance_calibration_and_signal_separation():
    rng = np.random.default_rng(2026)
    n_y, n_c, n_g, n_b = 18, 36, 120, 4
    is_comp = np.array([False] * n_y + [True] * n_c)
    # Create batch-correlated covariate Z with unbalanced old/young fractions
    batches = rng.choice(n_b, size=n_y + n_c, p=[0.4, 0.3, 0.2, 0.1])
    # Slightly skew batches across age
    batches[:10] = 0
    batches[-12:] = 2
    b_dum = pd.get_dummies(pd.Series(batches), drop_first=True, dtype=float).to_numpy()
    depth_z = rng.normal(0.0, 1.0, size=n_y + n_c)
    z_mat = np.column_stack([np.ones(n_y + n_c), depth_z, b_dum])

    # 1. Pure homoscedastic null with strong batch mean offsets -> mean LVR near 0, no negative df bias
    batch_effects = np.outer(batches - 1.0, rng.normal(0.0, 1.2, size=n_g))
    y_null = batch_effects + 0.5 * depth_z[:, None] + rng.normal(0.0, 1.0, size=(n_y + n_c, n_g))
    res_null = compute_a2_variance_metrics(y_null, is_comp, cov_mat_with_intercept=z_mat, n_perm=400, seed=7, eps=1e-8)
    assert abs(res_null["mean_log_var_ratio"]) < 0.08
    assert res_null["p_lvr_one_sided_comp_gt_young"] > 0.05
    assert res_null["p_eta2_perm"] > 0.05

    # 2. Pure DEG mean shift with equal residual variance -> eta2 significant, LVR remains null
    y_deg = y_null.copy()
    y_deg[is_comp, :40] += 1.0
    res_deg = compute_a2_variance_metrics(y_deg, is_comp, cov_mat_with_intercept=z_mat, n_perm=400, seed=8, eps=1e-8)
    assert res_deg["p_eta2_perm"] < 0.01
    assert abs(res_deg["mean_log_var_ratio"]) < 0.12
    assert res_deg["p_lvr_one_sided_comp_gt_young"] > 0.05

    # 3. True excess variance in OLD -> positive LVR & significant p_lvr
    y_var = y_null.copy()
    y_var[is_comp] += rng.normal(0.0, 1.25, size=(n_c, n_g))
    res_var = compute_a2_variance_metrics(y_var, is_comp, cov_mat_with_intercept=z_mat, n_perm=400, seed=9, eps=1e-8)
    assert res_var["mean_log_var_ratio"] > 0.35
    assert res_var["p_lvr_one_sided_comp_gt_young"] < 0.01


def test_atomic_batch_pruning_and_singleton_exclusion():
    df = pd.DataFrame(
        {
            "donor_id": [f"D{i}" for i in range(16)],
            "age_per_decade": np.linspace(2.0, 8.0, 16),
            "sex_male": [0, 1] * 8,
            "mean_log_umi_all": np.linspace(7.0, 8.5, 16),
            "batch": ["B1"] * 7 + ["B2"] * 8 + ["SINGLETON"],
        }
    )
    sub, active, omitted, status = _prune_covariate_matrix(
        df,
        nonbatch_covs=["sex_male", "mean_log_umi_all"],
        batch_col="batch",
        min_batch_donors=2,
        min_residual_df=5,
        ref_level="B1",
        require_is_comp_collinear_check="age_per_decade",
        atomic_batch_block=True,
    )
    assert status == "evaluated"
    assert len(sub) == 15  # singleton excluded before fitting
    assert "batch_B2" in active
    assert "SINGLETON" not in sub["batch"].values

    # Too few residual df forces atomic drop of whole batch block while keeping non-collinear nonbatch covariates
    small = df.iloc[[0, 1, 7, 8, 10, 11, 12, 13]].copy()
    _, active_s, omitted_s, status_s = _prune_covariate_matrix(
        small,
        nonbatch_covs=["sex_male", "mean_log_umi_all"],
        batch_col="batch",
        min_batch_donors=2,
        min_residual_df=10,
        ref_level="B1",
        require_is_comp_collinear_check="age_per_decade",
        atomic_batch_block=True,
    )
    assert status_s == "not_evaluable_within_batch_df"
    assert [c for c in active_s if c.startswith("batch_")] == []
    assert "batch_block(df_limit)" in omitted_s


def test_stratified_global_excluded_from_bh_family_and_shared_helper():
    rng = np.random.default_rng(11)
    rows = []
    for ct in ("CT_A", "CT_B", "CT_C"):
        for d in range(25):
            age = 2.0 + 0.25 * d
            rows.append(
                {
                    "donor_id": f"D{d}",
                    "age_per_decade": age,
                    "sex_male": d % 2,
                    "n_cells_ge1500": 50,
                    "mean_log_umi_ge1500": 7.5 + 0.01 * d + float(rng.normal(0.0, 0.2)),
                    "group": ct,
                    "val": 0.1 * age + rng.normal(0.0, 0.2),
                }
            )
    df = pd.DataFrame(rows)
    srow = _fit_stratified_global_row(
        df,
        elig_cts=["CT_A", "CT_B", "CT_C"],
        out_col="val",
        cell_col="n_cells_ge1500",
        depth_col="mean_log_umi_ge1500",
        ds_name="gao",
        part="B",
        metric="shannon",
        depth="1500",
        ref_type="none",
        include_part_a_cols=False,
    )
    assert srow is not None
    assert srow["is_stratified_global"] is True
    assert srow["group"] == "Stratified_CellType_Global"

    # Build slice table with 3 CTs + 1 pooled + 1 stratified
    slice_df = pd.DataFrame(
        [
            {"dataset": "gao", "part": "B", "metric": "shannon", "depth": "1500", "reference_type": "none", "model_spec": "primary", "group": "CT_A", "is_pooled": False, "is_stratified_global": False, "pvalue": 0.01},
            {"dataset": "gao", "part": "B", "metric": "shannon", "depth": "1500", "reference_type": "none", "model_spec": "primary", "group": "CT_B", "is_pooled": False, "is_stratified_global": False, "pvalue": 0.20},
            {"dataset": "gao", "part": "B", "metric": "shannon", "depth": "1500", "reference_type": "none", "model_spec": "primary", "group": "CT_C", "is_pooled": False, "is_stratified_global": False, "pvalue": 0.50},
            {"dataset": "gao", "part": "B", "metric": "shannon", "depth": "1500", "reference_type": "none", "model_spec": "primary", "group": "All_LV_Nuclei", "is_pooled": True, "is_stratified_global": False, "pvalue": 0.03},
            {"dataset": "gao", "part": "B", "metric": "shannon", "depth": "1500", "reference_type": "none", "model_spec": "stratified_global_primary", "group": "Stratified_CellType_Global", "is_pooled": True, "is_stratified_global": True, "pvalue": srow["pvalue"]},
        ]
    )
    bh = apply_slice_bh(slice_df, ["dataset", "part", "metric", "depth", "reference_type", "model_spec"])
    ct_m = bh.loc[~bh["is_pooled"], "family_size_m"].unique().tolist()
    assert ct_m == [3.0]
    assert bh.loc[bh["group"] == "All_LV_Nuclei", "family_size_m"].iloc[0] == 1.0
    assert pd.isna(bh.loc[bh["group"] == "Stratified_CellType_Global", "q_bh_within_slice"].iloc[0])


def test_verdict_engine_strict_boundaries_and_sign_guards():
    def _mk_pa(a1_f_s, a1_f_p, a1_f_old, a1_f_sh, a1_f_mde, a1_15_s, a1_15_p, a3_15_s, a3_15_p, a3_f_s, a3_15_mde):
        return pd.DataFrame(
            [
                {"part": "A1", "depth": "full_1500", "reference_type": "young_loo", "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": a1_f_s, "pvalue": a1_f_p, "old_mean": a1_f_old, "split_half_young_baseline": a1_f_sh, "mde80_sd_per_decade": a1_f_mde},
                {"part": "A1", "depth": "1500", "reference_type": "young_loo", "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": a1_15_s, "pvalue": a1_15_p, "old_mean": a1_f_old, "split_half_young_baseline": a1_f_sh, "mde80_sd_per_decade": a1_f_mde},
                # Ensure Stratified_CellType_Global cannot spoof All_Cells
                {"part": "A1", "depth": "full_1500", "reference_type": "young_loo", "model_spec": "primary", "group": "All_Cells", "is_stratified_global": True, "slope_raw_per_decade": 9.9, "pvalue": 1e-9, "old_mean": 9.9, "split_half_young_baseline": 0.1, "mde80_sd_per_decade": 0.1},
                {"part": "A3", "depth": "1500", "reference_type": "none", "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": a3_15_s, "pvalue": a3_15_p, "mde80_sd_per_decade": a3_15_mde},
                {"part": "A3", "depth": "full_1500", "reference_type": "none", "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": a3_f_s, "pvalue": a3_15_p, "mde80_sd_per_decade": a3_15_mde},
            ]
        )

    var_df = pd.DataFrame(
        [
            {"contrast": "OLD_vs_YOUNG", "arm": "batch_residualized", "eps_used": 1e-8, "group": "All_Cells", "is_stratified_global": False, "mean_log_var_ratio": -0.02, "p_eta2_perm": 0.30, "p_lvr_one_sided_comp_gt_young": 0.70},
            {"contrast": "OLD_vs_YOUNG", "arm": "depth_only_residualized", "eps_used": 1e-8, "group": "All_Cells", "is_stratified_global": False, "mean_log_var_ratio": 0.01, "p_eta2_perm": 0.25, "p_lvr_one_sided_comp_gt_young": 0.40},
        ]
    )
    spread_rows = []
    rob_rows = []
    # Shannon decreases significantly across all 3 depths -> must be DOES NOT PROVE with OPPOSITE note!
    for d in ("500", "1500", "3000"):
        spread_rows.append({"metric": "shannon", "depth": d, "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": -0.05, "pvalue": 0.001, "mde80_sd_per_decade": 0.20})
        spread_rows.append({"metric": "deman_kl_methods", "depth": d, "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": -0.05, "pvalue": 0.001, "mde80_sd_per_decade": 0.20})
        spread_rows.append({"metric": "deman_kl_sqrt", "depth": d, "model_spec": "primary", "group": "All_Cells", "is_stratified_global": False, "slope_raw_per_decade": 0.01, "pvalue": 0.70, "mde80_sd_per_decade": 0.75})
    for m, s_pri in (("shannon", -0.05), ("deman_kl_methods", -0.05), ("deman_kl_sqrt", 0.01)):
        rob_rows.append({"metric": m, "depth": 1500, "model_spec": "primary", "group": "All_Cells", "n_seeds": 20 if m != "shannon" else 5, "n_same_sign": 20 if m != "shannon" else 5, "passes_robustness": True})

    # Case 1: A1 full_1500 significant but OLD <= split_half baseline -> strictly DOES NOT PROVE!
    pa_fail_sh = _mk_pa(0.05, 0.001, 0.25, 0.26, 0.25, 0.05, 0.001, -0.05, 0.001, -0.04, 0.25)
    v1 = evaluate_verdicts(pa_fail_sh, var_df, pd.DataFrame(spread_rows), pd.DataFrame(rob_rows), 1500, "All_Cells")
    assert v1["A1"]["verdict"] == "DOES NOT PROVE"
    assert v1["A2"]["verdict"] == "DOES NOT PROVE"
    assert v1["A3"]["verdict"] == "DOES NOT PROVE"  # Negative slope never proves increase
    assert v1["Part_B"]["shannon"]["verdict"] == "DOES NOT PROVE"
    assert "OPPOSITE" in v1["Part_B"]["shannon"]["reason"]
    assert v1["Part_B"]["deman_kl_methods"]["verdict"] == "PROVES"  # KL lower => more uniform
    assert v1["Part_B"]["deman_kl_sqrt"]["verdict"] == "INCONCLUSIVE"  # MDE80 > 0.5

    # Case 2: A1 full & 1500 significant & OLD > split_half -> PROVES
    pa_pass = _mk_pa(0.05, 0.001, 0.30, 0.26, 0.25, 0.05, 0.001, 0.05, 0.001, 0.04, 0.25)
    v2 = evaluate_verdicts(pa_pass, var_df, pd.DataFrame(spread_rows), pd.DataFrame(rob_rows), 1500, "All_Cells")
    assert v2["A1"]["verdict"] == "PROVES"
    assert v2["A3"]["verdict"] == "PROVES"


def test_split_half_fallback_propagates_to_a1_models():
    rows = []
    for i in range(12):
        age = 22.0 + 6.0 * i
        b = "YOUNG" if age < 40 else ("MID" if age < 60 else "OLD")
        rows.append(
            {
                "donor_id": f"D{i}",
                "age_years": age,
                "age_per_decade": age / 10.0,
                "bin": b,
                "sex_male": i % 2,
                "study_or_pool": "P1" if i < 6 else "P2",
                "technology": "10x",
                "split_half_group": 1,  # Deliberate global H1-only corner case!
                "group": "RareCT",
                "is_pooled": False,
                "n_cells_all": 30,
                "n_cells_ge500": 30,
                "n_cells_ge1500": 30,
                "n_cells_ge3000": 30,
                "mean_log_umi_all": 7.6,
                "mean_log_umi_ge500": 7.6,
                "mean_log_umi_ge1500": 7.6,
                "mean_log_umi_ge3000": 7.6,
                "a1_dist_full_1500": 0.20 + 0.01 * i,
                "a1_dist_1500": 0.21 + 0.01 * i,
                "a1_dist_500": 0.22 + 0.01 * i,
                "a1_dist_3000": 0.19 + 0.01 * i,
                "a1_dist_mature_ref_mid_loo": 0.20 + 0.01 * i if age >= 40 else np.nan,
                "a1_dist_mature_ref_1500": 0.21 + 0.01 * i if age >= 40 else np.nan,
                "a1_self_non_loo": 0.15 if b == "YOUNG" else np.nan,
                "a1_split_half_h2_to_h1": 0.19 if i == 1 else np.nan,
                "a3_disp_1500": 0.30,
                "a3_disp_full_1500": 0.29,
                "eligible_a1_1500": True,
                "eligible_a1_3000": True,
                "eligible_a3_1500": True,
                "yazar_pool_fe_all_evaluable": True,
            }
        )
    donors_long = pd.DataFrame(rows)
    refs_summary = {
        "ct_full_1500": {
            "RareCT": {
                "half1": {
                    "h1_n_donors": 2,
                    "h1_n_cells": 60,
                    "h2_donors": {"D1"},
                    "fallback_used": True,
                }
            }
        }
    }
    out_rows = []
    _run_a1_a3_models_from_donors("yazar", donors_long, refs_summary, out_rows)
    m_df = pd.DataFrame(out_rows)
    pri = m_df[(m_df["part"] == "A1") & (m_df["depth"] == "full_1500") & (m_df["reference_type"] == "young_loo") & (m_df["model_spec"] == "primary")].iloc[0]
    assert bool(pri["split_half_fallback_flag"]) is True
    assert int(pri["split_half_n_donors_h1"]) == 2
    assert int(pri["split_half_n_donors_h2"]) == 1
    assert np.isclose(pri["split_half_young_baseline"], 0.19)


def test_generated_deliverables_hygiene_and_invariants_when_present():
    expected_csvs = [
        RESULTS_DIR / "followup-between-subject-gao-donors.csv",
        RESULTS_DIR / "followup-between-subject-gao-models.csv",
        RESULTS_DIR / "followup-between-subject-variance-gao.csv",
        RESULTS_DIR / "followup-between-subject-yazar-donors.csv",
        RESULTS_DIR / "followup-between-subject-yazar-models.csv",
        RESULTS_DIR / "followup-between-subject-variance-yazar.csv",
        RESULTS_DIR / "followup-spread-gao-models.csv",
        RESULTS_DIR / "followup-spread-yazar-models.csv",
        RESULTS_DIR / "followup-between-subject-robustness.csv",
    ]
    for path in expected_csvs:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        assert "/home/reza" not in text
        assert "\u2013" not in text and "\u2014" not in text
        df = pd.read_csv(path)
        assert not df.empty
        if "q_bh_within_slice" in df.columns and "pvalue" in df.columns:
            m = df["q_bh_within_slice"].notna() & df["pvalue"].notna()
            assert (df.loc[m, "q_bh_within_slice"] >= df.loc[m, "pvalue"] - 1e-12).all()

    qc_path = RESULTS_DIR / "followup-between-subject-gao-yazar-qc.json"
    if qc_path.exists():
        qc_txt = qc_path.read_text(encoding="utf-8")
        assert "/home/reza" not in qc_txt
        assert "\u2013" not in qc_txt and "\u2014" not in qc_txt
        qc = json.loads(qc_txt)
        assert set(qc["verdicts"].keys()) == {"gao", "yazar"}

    rep_path = REPORTS_DIR / "followup-between-subject-and-spread-gao-yazar.md"
    if rep_path.exists():
        rep_txt = rep_path.read_text(encoding="utf-8")
        assert "/home/reza" not in rep_txt
        assert "\u2013" not in rep_txt and "\u2014" not in rep_txt
