import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm

from src.followup_muscle import (
    MIN_CELLS,
    MIN_DONORS,
    MIN_PER_GROUP,
    MIN_REFERENCE_CELLS,
    NULL_ITERATIONS,
    NULL_SEED,
    assert_count_gate,
    build_design,
    clr_transform,
    fit_ols,
    hellinger_distances,
    matched_null_effect_pvalue,
    random_signature_scores,
)

REPO = Path(__file__).resolve().parent.parent


def test_frozen_plan_preserves_donor_and_age_rules():
    plan = json.loads((REPO / "data/manifest/muscle_analysis_plan.json").read_text())
    assert plan["frozen_before_expression_download_or_outcomes"]
    assert plan["execution_order"] == ["Kedlian", "Perez"]
    assert plan["resource_guard"]["minimum_free_gb"] == 30
    assert plan["resource_guard"]["analysis_memory_gb"] == 6
    assert plan["resource_guard"]["kedlian_maximum_download_gb"] == 5
    assert "470BR" in plan["Kedlian"]["primary_donor_set"]
    assert "470BR" in plan["Kedlian"]["sensitivity_donor_set"]
    assert plan["Kedlian"]["endpoints"]["Hellinger_displacement"]["status"] == "nonnumeric_age_group_sensitivity_only"
    assert plan["Perez"]["donor_groups"] == {"Young": 6, "Old": 11}
    assert plan["Perez"]["bulk_72_man_arm"].startswith("excluded")
    assert not plan["common_design"]["studies_pooled"]


def test_count_gate_rejects_nonfinite_negative_and_fractional_values():
    assert_count_gate(np.array([0, 1, 2**20], dtype=np.float32))
    for bad in ([np.nan], [np.inf], [-1], [1.5]):
        try:
            assert_count_gate(np.asarray(bad, dtype=np.float32))
        except ValueError as exc:
            assert "Raw-count gate failed" in str(exc)
        else:
            raise AssertionError("Invalid raw counts passed the gate")


def test_clr_zero_handling_and_closure():
    counts = pd.DataFrame({"a": [0, 2, 0], "b": [3, 0, 0], "c": [0, 0, 4]})
    transformed = clr_transform(counts)
    assert np.isfinite(transformed.to_numpy()).all()
    np.testing.assert_allclose(transformed.sum(axis=1), 0.0, atol=1e-12)
    assert transformed.loc[0, "a"] < 0
    assert transformed.loc[2, "c"] > 0


def test_group_design_requires_group_donor_floors_and_full_rank():
    donors = pd.DataFrame({
        "old": [0, 0, 0, 0, 1, 1, 1, 1],
        "male": [0, 1, 0, 1, 1, 0, 1, 0],
        "dcd": [0, 0, 1, 1, 1, 1, 0, 0],
        "mean_log_umi": np.linspace(7, 8, 8),
    })
    design, names, omitted = build_design(donors)
    assert names == ["intercept", "old", "male", "dcd", "mean_log_umi"]
    assert omitted == []
    assert np.linalg.matrix_rank(design) == design.shape[1]
    donors["dcd"] = donors["old"]
    _, names, omitted = build_design(donors)
    assert "dcd" in omitted and "dcd" not in names
    donors.loc[[0, 1, 2], "old"] = 1
    try:
        build_design(donors)
    except ValueError as exc:
        assert "at least two donors" in str(exc)
    else:
        raise AssertionError("A one-donor group passed")
    assert (MIN_DONORS, MIN_PER_GROUP, MIN_CELLS, MIN_REFERENCE_CELLS) == (6, 2, 20, 50)


def test_vectorized_ols_matches_statsmodels_effect_interval():
    donors = pd.DataFrame({
        "old": [0, 0, 0, 0, 1, 1, 1, 1],
        "male": [0, 1, 0, 1, 1, 0, 1, 0],
        "dcd": [0, 0, 1, 1, 1, 1, 0, 0],
        "mean_log_umi": np.linspace(7, 8, 8),
    })
    design, _, _ = build_design(donors)
    outcomes = np.column_stack([
        0.2 * donors["old"] + 0.01 * donors["mean_log_umi"]
        + np.array([0.01, -0.02, 0.03, 0.01, -0.01, 0.02, -0.03, 0.0]),
        np.linspace(-1, 1, 8) + np.array([0.1, 0, -0.1, 0.2, 0, -0.2, 0.1, 0]),
    ])
    fitted = fit_ols(outcomes, design)
    reference = sm.OLS(outcomes[:, 0], design).fit()
    np.testing.assert_allclose(fitted["effect"][0], reference.params[1])
    np.testing.assert_allclose(fitted["se"][0], reference.bse[1])
    np.testing.assert_allclose(fitted["pvalue"][0], reference.pvalues[1])
    np.testing.assert_allclose(fitted["ci95_low"][0], reference.conf_int()[1, 0])
    np.testing.assert_allclose(fitted["ci95_high"][0], reference.conf_int()[1, 1])


def test_expression_matched_null_is_reproducible_and_centered():
    rng = np.random.default_rng(4)
    matrix = rng.normal(size=(60, 8))
    up = np.array([0, 1, 2])
    down = np.array([3, 4])
    first = random_signature_scores(matrix, up, down)
    second = random_signature_scores(matrix, up, down)
    np.testing.assert_allclose(first, second)
    assert first.shape == (NULL_ITERATIONS, 8)
    assert NULL_SEED == 20261105
    effects = first.mean(axis=1)
    observed = effects.mean() + effects.std()
    expected = (1 + (np.abs(effects - effects.mean()) >= abs(observed - effects.mean())).sum()) / 201
    assert matched_null_effect_pvalue(observed, effects) == expected


def test_hellinger_distance_zero_for_identical_profiles():
    counts = sp.csr_matrix([[1, 1, 2], [2, 2, 4]], dtype=float)
    reference = np.array([0.25, 0.25, 0.5])
    totals = np.asarray(counts.sum(axis=1)).ravel()
    distances = hellinger_distances(counts, reference, totals)
    np.testing.assert_allclose(distances, 0.0, atol=1e-12)


def test_kedlian_pipeline_synthetic_end_to_end(tmp_path, monkeypatch):
    import anndata as ad
    from scripts import run_followup_kedlian as kedlian
    from src.aging_programs import load_program_tiers

    selected = ["583B", "411C", "464C", "591C", "621B", "582C", "362C", "343B"]
    manifest = json.loads((REPO / "data/manifest/E-MTAB-13874.json").read_text())
    samples = {
        donor["id"]: donor["scRNA_sample_ids"][0]
        for donor in manifest["donors"] if donor["id"] in selected
    }
    obs = pd.DataFrame({
        "sample": np.repeat([samples[donor] for donor in selected], 20),
        "cell": "typeA",
        "nuclei": "nan",
    })
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    genes = np.asarray(program["up"] + program["down"] + [f"NULL{i}" for i in range(500)])
    rng = np.random.default_rng(81)
    counts = rng.poisson(1, size=(len(obs), len(genes))) + 3
    adata = ad.AnnData(
        X=sp.csr_matrix(counts.astype(np.float32)), obs=obs,
        var=pd.DataFrame({"feature_name": genes}, index=[f"ENSG{i}" for i in range(len(genes))]),
    )
    adata.layers["counts"] = adata.X.copy()
    path = tmp_path / "synthetic.h5ad"
    adata.write_h5ad(path)
    monkeypatch.setattr(kedlian, "H5AD", path)
    monkeypatch.setattr(kedlian, "CACHE", tmp_path / "cache.npz")
    monkeypatch.setattr(kedlian, "RESULTS", tmp_path)
    kedlian.build_cache()
    verification = json.loads((tmp_path / "followup-muscle-kedlian-verification.json").read_text())
    assert verification["all_counts_finite_nonnegative_integer_valued"]
    kedlian.analyze()
    p1 = pd.read_csv(tmp_path / "followup-muscle-kedlian-p1-models.csv")
    assert set(p1["assay"]) == {"scRNA"}
    assert set(p1["cell_type"]) == {"typeA"}
    assert (p1["n_donors"] == 8).all()
    assert (p1["residual_df"] >= 2).all()
    genes_out = pd.read_csv(tmp_path / "followup-muscle-kedlian-scRNA-gene-models.csv")
    assert genes_out["q_bh"].between(0, 1).all()
    h_eligible = pd.read_csv(tmp_path / "followup-muscle-kedlian-hellinger-eligibility.csv")
    assert h_eligible["hellinger_eligible"].all()
    h = pd.read_csv(tmp_path / "followup-muscle-kedlian-hellinger-models.csv")
    assert h["analysis_label"].eq("NONNUMERIC_AGE_GROUP_SENSITIVITY").all()


def test_perez_whole_nucleus_pipeline_synthetic(tmp_path, monkeypatch):
    import h5py
    from scripts import run_followup_perez as perez

    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads((REPO / "data/manifest/GSE167186.json").read_text())
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    rng = np.random.default_rng(93)
    for donor in manifest["donors"]:
        donor["filtered_nuclei"] = 5
        matrix = sp.csc_matrix(rng.poisson(1, size=(30, 5)).astype(np.int32))
        path = data_dir / f"{donor['gsm']}_{donor['hm_id']}_filtered_feature_bc_matrix.h5"
        with h5py.File(path, "w") as handle:
            group = handle.create_group("matrix")
            group.create_dataset("data", data=matrix.data)
            group.create_dataset("indices", data=matrix.indices)
            group.create_dataset("indptr", data=matrix.indptr)
            group.create_dataset("shape", data=matrix.shape)
            group.create_dataset("barcodes", data=np.asarray([f"cell{i}" for i in range(5)], dtype="S"))
            features = group.create_group("features")
            features.create_dataset("id", data=np.asarray([f"ENSG{i}" for i in range(30)], dtype="S"))
            features.create_dataset("name", data=np.asarray([f"GENE{i}" for i in range(30)], dtype="S"))
        donor["h5_size_bytes"] = path.stat().st_size
    manifest_path.write_text(json.dumps(manifest))
    monkeypatch.setattr(perez, "DATA", data_dir)
    monkeypatch.setattr(perez, "RESULTS", tmp_path)
    monkeypatch.setattr(perez, "MANIFEST", manifest_path)
    monkeypatch.setattr(perez, "EXPECTED_NUCLEI", 85)
    perez.run()
    donors = pd.read_csv(tmp_path / "followup-muscle-perez-donor-summary.csv")
    assert len(donors) == 17
    assert donors["n_nuclei"].sum() == 85
    genes = pd.read_csv(tmp_path / "followup-muscle-perez-gene-models.csv")
    assert set(genes["model"]) == {"primary_group_depth_adjusted", "exact_age_10y_sensitivity"}
    assert genes["composition_confounded"].all()
    assert genes["q_bh"].between(0, 1).all()
    endpoints = pd.read_csv(tmp_path / "followup-muscle-perez-endpoint-status.csv").set_index("endpoint")
    assert not endpoints.loc["P1_transport", "outcomes_examined"]
    assert not endpoints.loc["donor_CLR_composition", "outcomes_examined"]
    assert not endpoints.loc["Hellinger_displacement", "outcomes_examined"]


def test_delivered_perez_outputs_and_multiplicity():
    verification = json.loads((REPO / "results/followup-muscle-perez-verification.json").read_text())
    assert verification["released_filtered_nuclei"] == 134486
    assert verification["stored_count_values"] == 58309907
    assert verification["all_counts_finite_nonnegative_integer_valued"]
    assert verification["peak_rss_bytes"] < 6 * 1024**3
    donors = pd.read_csv(REPO / "results/followup-muscle-perez-donor-summary.csv")
    assert len(donors) == 17
    assert donors["n_nuclei"].sum() == 134486
    assert donors.groupby("group")["gsm"].nunique().to_dict() == {"Old": 11, "Young": 6}
    assert donors["h5_root_groups"].eq("matrix").all()
    assert not donors["released_nucleus_annotations_present"].any()
    genes = pd.read_csv(REPO / "results/followup-muscle-perez-gene-models.csv")
    assert len(genes) == 40686
    assert genes.groupby("model").size().eq(20343).all()
    assert genes["q_bh"].between(0, 1).all()
    assert (genes["ci95_low"] <= genes["effect"]).all()
    assert (genes["effect"] <= genes["ci95_high"]).all()
    assert genes["composition_confounded"].all()
    assert (genes["residual_df"] == 14).all()
    from statsmodels.stats.multitest import multipletests
    for _, part in genes.groupby("model"):
        expected = multipletests(part["pvalue"], method="fdr_bh")[1]
        np.testing.assert_allclose(part["q_bh"], expected)
    endpoints = pd.read_csv(REPO / "results/followup-muscle-perez-endpoint-status.csv").set_index("endpoint")
    for endpoint in ["P1_transport", "cell_type_models", "donor_CLR_composition", "Hellinger_displacement"]:
        assert endpoints.loc[endpoint, "status"].startswith("blocked")
        assert not endpoints.loc[endpoint, "outcomes_examined"]
    assert endpoints.loc["bulk_72_man_arm", "status"] == "excluded"


def test_kedlian_incomplete_download_has_no_outcome_claims():
    resource = json.loads((REPO / "results/followup-muscle-kedlian-download-resource.json").read_text())
    assert resource["partial_bytes"] < resource["expected_bytes"]
    assert not resource["expression_values_read"]
    assert not resource["count_validation_performed"]
    assert not resource["biological_outcomes_examined"]
    assert not resource["alternate_cellxgene_object_substituted"]
    manifest = json.loads((REPO / "data/manifest/E-MTAB-13874.json").read_text())
    assert manifest["analysis"]["status"] == "not_run_incomplete_required_raw_download"
    assert not manifest["verdict"]["expression_analysis_performed"]
    assert not (REPO / "results/followup-muscle-kedlian-p1-models.csv").exists()
    report = (REPO / "reports/followup-muscle-analysis.md").read_text()
    assert "Step 4 is incomplete" in report
    assert "no Kedlian expression values" in report
