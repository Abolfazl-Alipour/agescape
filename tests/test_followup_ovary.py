import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from src.followup_ovary import build_jin_design, build_postmenopausal_design

REPO = Path(__file__).resolve().parent.parent


def test_frozen_ovary_plans_preserve_sequence_units_and_blocked_arms():
    jin = json.loads((REPO / "data/manifest/ovary_jin_analysis_plan.json").read_text())
    post = json.loads((REPO / "data/manifest/ovary_postmenopausal_analysis_plan.json").read_text())
    assert jin["models_frozen_before_expression_download_or_outcomes"]
    assert post["models_frozen_before_expression_download_or_outcomes"]
    assert jin["current_plan_frozen_before_outcomes"]
    assert post["current_plan_frozen_before_outcomes"]
    assert (jin["execution_order"], post["execution_order"]) == (1, 2)
    assert jin["resources"]["analysis_memory_gb"] == 8
    assert post["resources"]["maximum_compressed_atlas_download_gb"] == 7
    assert jin["resources"]["minimum_free_gb"] == post["resources"]["minimum_free_gb"] == 30
    assert jin["donors"]["continuous_age_model"] == "blocked_not_planned"
    assert "UNRESOLVED_47" in jin["donors"]["ovary5_exact_age"]
    assert not jin["common"]["cross_cohort_pooling"]
    assert not post["common"]["cross_cohort_pooling"]
    assert post["endpoints"]["Hellinger_displacement"]["status"] == "untestable_no_20_to_40_year_reference"
    assert "not separately released" in post["count_processing"]["uncorrected_sensitivity"]
    assert jin["endpoints"]["same_nucleus_RNA_ATAC_coupling"]["status"].startswith("blocked")


def test_jin_design_requires_exact_four_versus_four():
    frame = pd.DataFrame({"old": [0, 0, 0, 0, 1, 1, 1, 1], "mean_log_umi": np.linspace(7, 8, 8)})
    design, names = build_jin_design(frame, "primary_group_depth_adjusted")
    assert names == ["intercept", "reproductively_aged", "mean_log_umi"]
    assert design.shape == (8, 3)
    unadjusted, _ = build_jin_design(frame, "unadjusted_group_sensitivity")
    assert unadjusted.shape == (8, 2)
    try:
        build_jin_design(frame.iloc[:7], "primary_group_depth_adjusted")
    except ValueError as exc:
        assert "four young" in str(exc)
    else:
        raise AssertionError("Seven donors passed the Jin design gate")


def post_frame(n=24):
    rng = np.random.default_rng(29)
    return pd.DataFrame({
        "age_years": np.linspace(50, 84, n),
        "mean_log_umi": rng.normal(6, 0.5, n),
        "plex_fraction_1": rng.uniform(0.05, 0.3, n),
        "plex_fraction_2": rng.uniform(0.05, 0.3, n),
        "plex_fraction_3": rng.uniform(0.05, 0.3, n),
        "BMI": rng.normal(30, 5, n),
        "non_ovarian_abnormal": np.tile([0, 0, 1], n // 3),
        "parity_live_births": np.tile([0, 1, 2, 3], n // 4),
    })


def test_postmenopausal_design_is_per_decade_and_matches_statsmodels():
    frame = post_frame()
    design, names = build_postmenopausal_design(frame, "primary_depth_plex_clinical")
    assert names[1] == "age_per_10y"
    assert design.shape == (24, 8)
    outcome = (0.4 * design[:, 1] - 0.2 * design[:, 2] + np.linspace(-0.1, 0.1, 24)
               + np.tile([0.01, -0.02, 0.015, -0.005], 6))
    from src.followup_muscle import fit_ols
    fitted = fit_ols(outcome, design)
    reference = sm.OLS(outcome, design).fit()
    np.testing.assert_allclose(fitted["effect"][0], reference.params[1])
    np.testing.assert_allclose(fitted["se"][0], reference.bse[1])
    np.testing.assert_allclose(fitted["ci95_low"][0], reference.conf_int()[1, 0])
    np.testing.assert_allclose(fitted["ci95_high"][0], reference.conf_int()[1, 1])
    technical, _ = build_postmenopausal_design(frame, "technical_depth_plex_sensitivity")
    assert technical.shape == (24, 6)
    unadjusted, _ = build_postmenopausal_design(frame, "unadjusted_sensitivity")
    assert unadjusted.shape == (24, 2)
    parity, _ = build_postmenopausal_design(frame, "secondary_primary_plus_parity")
    assert parity.shape == (24, 9)


def test_postmenopausal_design_rejects_rank_failure_and_small_donor_sets():
    frame = post_frame()
    frame["BMI"] = frame["age_years"]
    try:
        build_postmenopausal_design(frame, "primary_depth_plex_clinical")
    except ValueError as exc:
        assert "full rank" in str(exc)
    else:
        raise AssertionError("Aliased BMI passed the frozen design")
    try:
        build_postmenopausal_design(post_frame(24).iloc[:19], "unadjusted_sensitivity")
    except ValueError as exc:
        assert "20" in str(exc)
    else:
        raise AssertionError("Nineteen donors passed the postmenopausal gate")


def test_jin_raw_gates_and_delivered_untestable_arms():
    verification = json.loads((REPO / "results/followup-ovary-jin-verification.json").read_text())
    assert verification["donors"] == 8
    assert verification["annotated_nuclei"] == 42568
    assert verification["raw_h5_barcodes"] == 43494
    assert verification["stored_count_values_checked"] == 69548923
    assert verification["all_counts_finite_nonnegative_integer_valued"]
    assert verification["all_annotations_matched_one_to_one"]
    assert not verification["ovary5_numeric_age_used"]
    assert verification["peak_rss_bytes"] < 8 * 1024**3
    eligible = pd.read_csv(REPO / "results/followup-ovary-jin-celltype-eligibility.csv")
    assert set(eligible.loc[~eligible["eligible"], "cell_type"]) == {"LEC", "TC"}
    assert eligible.loc[eligible["eligible"], "n_young_ge20"].eq(4).all()
    assert eligible.loc[eligible["eligible"], "n_reproductively_aged_ge20"].eq(4).all()
    status = pd.read_csv(REPO / "results/followup-ovary-jin-model-status.csv")
    blocked = status[status["status"].ne("completed")]
    assert not blocked["outcomes_examined"].any()
    for endpoint in ["continuous_chronological_age", "same_nucleus_RNA_ATAC_coupling", "donor_matched_RNA_ATAC_linkage"]:
        assert endpoint in set(status["endpoint"])


def test_jin_effect_intervals_and_explicit_bh_families():
    for name in ["p1-models", "composition-models", "hellinger-models", "gene-models"]:
        frame = pd.read_csv(REPO / f"results/followup-ovary-jin-{name}.csv")
        assert (frame["ci95_low"] <= frame["effect"]).all()
        assert (frame["effect"] <= frame["ci95_high"]).all()
        assert frame["q_bh"].between(0, 1).all()
        assert frame["n_donors"].eq(8).all()
        for _, part in frame.groupby("model"):
            expected = multipletests(part["pvalue"], method="fdr_bh")[1]
            np.testing.assert_allclose(part["q_bh"], expected)
    p1 = pd.read_csv(REPO / "results/followup-ovary-jin-p1-models.csv")
    primary = p1[p1["model"].eq("primary_group_depth_adjusted")]
    expected = multipletests(primary["random_pvalue"], method="fdr_bh")[1]
    np.testing.assert_allclose(primary["random_q_bh"], expected)
    assert primary["random_iterations"].eq(200).all()
    assert primary["random_seed"].eq(20261105).all()


def test_jin_no_continuous_or_cross_cohort_output():
    results = list((REPO / "results").glob("followup-ovary-jin-*.csv"))
    for path in results:
        frame = pd.read_csv(path, nrows=5)
        assert "age_per_10y" not in frame.columns
    p1 = pd.read_csv(REPO / "results/followup-ovary-jin-p1-models.csv")
    assert set(p1["effect_definition"]) == {"reproductively_aged_minus_young"}
    h = pd.read_csv(REPO / "results/followup-ovary-jin-hellinger-eligibility.csv")
    assert set(h.loc[h["hellinger_eligible"], "cell_type"]) == {"BEC", "GC", "IC", "SC", "SMC"}


def test_postmenopausal_pipeline_synthetic_atlas_end_to_end(tmp_path, monkeypatch):
    import anndata as ad
    import scipy.sparse as sp
    from scripts import run_followup_ovary_postmenopausal as post
    from src.aging_programs import load_program_tiers

    manifest = json.loads((REPO / "data/manifest/ovary_cohorts.json").read_text())["GSE329417"]
    cell_types = [f"CT{i:02d}" for i in range(21)]
    donor_map = {donor["donor"]: donor for donor in manifest["donors"]}
    obs_parts = []
    for sample_index, sample in enumerate(manifest["samples"]):
        repeats = 20 + (np.arange(21) + sample_index) % 5
        donor = donor_map[sample["donor"]]
        frame = pd.DataFrame({
            "SenNet_ID": sample["donor"],
            "unique_sample_ID": post.atlas_sample_identifier(sample),
            "cell_type": np.repeat(cell_types, repeats),
            "Age": donor["exact_age_years"],
            "BMI": donor["BMI"],
            "Para": donor["parity_live_births"],
            "NonOvarian_Cancer_State": donor["non_ovarian_cancer_state"],
            "plex": sample["plex"],
        })
        obs_parts.append(frame)
    obs = pd.concat(obs_parts, ignore_index=True)
    obs.index = [f"nucleus{i}" for i in range(len(obs))]
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    genes = np.asarray(program["up"] + program["down"] + [f"NULL{i}" for i in range(11)])
    rng = np.random.default_rng(20261114)
    counts = rng.poisson(0.5, size=(len(obs), len(genes))).astype(np.float32)
    counts[:, 0] += 1
    adata = ad.AnnData(
        X=sp.csr_matrix(counts), obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name=None)),
    )
    adata.layers["raw_counts"] = adata.X.copy()
    adata.layers["log1p_norm"] = adata.X.copy()
    path = tmp_path / "synthetic_atlas.h5ad"
    adata.write_h5ad(path)
    monkeypatch.setattr(post, "ATLAS", path)
    monkeypatch.setattr(post, "CACHE", tmp_path / "cache.npz")
    monkeypatch.setattr(post, "RESULTS", tmp_path)
    monkeypatch.setattr(post, "EXPECTED_NUCLEI", len(obs))
    post.build_cache()
    verification = json.loads((tmp_path / "followup-ovary-postmenopausal-verification.json").read_text())
    assert verification["all_counts_finite_nonnegative_integer_valued"]
    assert verification["x_identical_to_log1p_norm_in_every_chunk"]
    assert verification["join"]["all_64_samples_present"]
    assert verification["join"]["technical_replicates_combined_within_donor"] == ["3_4", "16_4"]
    eligibility = pd.read_csv(tmp_path / "followup-ovary-postmenopausal-celltype-eligibility.csv")
    assert eligibility["eligible"].all()
    assert eligibility["n_donors_ge20"].eq(28).all()
    post.analyze()
    p1 = pd.read_csv(tmp_path / "followup-ovary-postmenopausal-p1-models.csv")
    assert p1["n_donors"].eq(28).all()
    assert p1["effect_definition"].eq("within_postmenopausal_per_10_exact_years").all()
    assert set(p1["model"]) == set(post.POST_MODELS)
    composition = pd.read_csv(tmp_path / "followup-ovary-postmenopausal-composition-models.csv")
    assert composition["q_bh"].between(0, 1).all()
    for model in post.POST_MODELS:
        genes_out = pd.read_csv(tmp_path / f"followup-ovary-postmenopausal-gene-models-{model}.csv")
        assert genes_out["q_bh"].between(0, 1).all()
        assert genes_out["n_donors"].eq(28).all()
    status = pd.read_csv(tmp_path / "followup-ovary-postmenopausal-model-status.csv")
    blocked = status[status["status"].ne("completed")]
    assert not blocked["outcomes_examined"].any()
    assert status.loc[status["endpoint"].eq("Hellinger_displacement"), "status"].eq("untestable").all()
    assert status.loc[status["endpoint"].eq("uncorrected_Cell_Ranger_sensitivity"), "status"].eq("blocked").all()


def test_postmenopausal_raw_gates_and_blocked_arms_delivered():
    verification = json.loads((REPO / "results/followup-ovary-postmenopausal-verification.json").read_text())
    assert verification["atlas_nuclei"] == 439011
    assert verification["donors"] == 28
    assert verification["samples"] == 64
    assert verification["coarse_cell_types"] == 21
    assert verification["all_counts_finite_nonnegative_integer_valued"]
    assert verification["x_identical_to_log1p_norm_in_every_chunk"]
    assert not verification["uncorrected_cell_ranger_counts_available"]
    assert not verification["uncorrected_sensitivity_performed"]
    assert verification["join"]["all_observations_exactly_mapped"]
    assert verification["peak_rss_bytes"] < 8 * 1024**3
    status = pd.read_csv(REPO / "results/followup-ovary-postmenopausal-model-status.csv")
    hellinger = status[status["endpoint"].eq("Hellinger_displacement")]
    assert len(hellinger) == 21
    assert hellinger["status"].eq("untestable").all()
    assert not hellinger["outcomes_examined"].any()
    uncorrected = status[status["endpoint"].eq("uncorrected_Cell_Ranger_sensitivity")]
    assert uncorrected["status"].eq("blocked").all()
    assert not uncorrected["outcomes_examined"].any()
    assert status.loc[status["status"].ne("completed"), "outcomes_examined"].eq(False).all()


def test_postmenopausal_effect_intervals_donor_units_and_bh_families():
    for name in ["p1-models", "composition-models"]:
        frame = pd.read_csv(REPO / f"results/followup-ovary-postmenopausal-{name}.csv")
        assert frame["n_donors"].ge(20).all()
        assert frame["n_donors"].le(28).all()
        assert frame["effect_definition"].eq("within_postmenopausal_per_10_exact_years").all()
        assert (frame["ci95_low"] <= frame["effect"]).all()
        assert (frame["effect"] <= frame["ci95_high"]).all()
        for _, part in frame.groupby("model"):
            expected = multipletests(part["pvalue"], method="fdr_bh")[1]
            np.testing.assert_allclose(part["q_bh"], expected)
    models = [
        "primary_depth_plex_clinical", "technical_depth_plex_sensitivity",
        "unadjusted_sensitivity", "secondary_primary_plus_parity",
    ]
    for model in models:
        frame = pd.read_csv(REPO / f"results/followup-ovary-postmenopausal-gene-models-{model}.csv")
        assert frame["n_donors"].ge(20).all()
        assert frame["effect_definition"].eq("within_postmenopausal_per_10_exact_years").all()
        expected = multipletests(frame["pvalue"], method="fdr_bh")[1]
        np.testing.assert_allclose(frame["q_bh"], expected)
        if model == "primary_depth_plex_clinical":
            assert frame["design_columns"].str.contains("BMI").all()
            assert frame["design_columns"].str.contains("non_ovarian_abnormal").all()
        if model == "secondary_primary_plus_parity":
            assert frame["design_columns"].str.contains("parity_live_births").all()
