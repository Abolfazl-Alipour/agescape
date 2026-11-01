import numpy as np
import pandas as pd

from scripts.run_human_vector_gene_drivers import gene_regressions


def test_gene_regressions_recover_known_age_effects():
    metadata = pd.DataFrame({
        "age": [20, 30, 40, 50, 60, 70, 80, 90],
        "sex_male": [0, 1, 0, 1, 0, 1, 0, 1],
        "mean_log_umi": [8, 8.4, 8.1, 8.7, 8.2, 8.5, 8.3, 8.8],
    })
    expression = np.column_stack([
        0.1 * metadata["age"] + metadata["sex_male"],
        -0.1 * metadata["age"] + 1 - metadata["sex_male"],
    ])
    results = pd.DataFrame(gene_regressions(expression, metadata, ["A", "B"], "c", "ct"))
    primary = results[results["model"] == "age_plus_sex"]
    np.testing.assert_allclose(primary["slope_log_cp10k_per_year"], [0.1, -0.1])
    adjusted = results[results["model"] == "age_plus_sex_depth"]
    np.testing.assert_allclose(adjusted["slope_log_cp10k_per_year"], [0.1, -0.1])


def test_human_primary_hits_do_not_survive_depth_adjustment():
    audit = pd.read_csv("results/h3-human-vector-design-audit.csv")
    hits = audit[audit["primary_fdr_hit"]]
    assert len(hits) == 6
    assert not hits["depth_adjusted_fdr_hit"].any()
    neuronal = hits[hits["cohort"] == "MSSM_DLPFC"]
    epithelial = hits[hits["cohort"] == "Lung_Integrated"]
    assert (neuronal["depth_age_slope_per_year"] < 0).all()
    assert (epithelial["depth_age_slope_per_year"] > 0).all()
    assert (hits["distance_depth_slope"] < 0).all()


def test_bayesian_hierarchy_maps_groups_to_units_before_observations():
    import pytest
    pytest.importorskip("pymc")
    from scripts.run_cross_species_bayesian_vector import prepare_data

    rows = []
    for cohort in ("a", "b"):
        for cell_type in ("neuron", "glia"):
            for index in range(6):
                rows.append({
                    "species": "Human", "tissue": "brain", "cohort": cohort,
                    "cell_type": cell_type, "donor": f"{cohort}-{index}",
                    "age": 20 + index * 10, "sex_male": index % 2,
                    "mean_log_umi": 8 + index * 0.1, "profile_hellinger": 0.1 + index * 0.01,
                })
    part, coords, data = prepare_data(pd.DataFrame(rows), "Human", True)
    assert len(data["unit_group_index"]) == len(coords["unit"])
    assert len(data["unit_index"]) == len(part)
    assert set(data["unit_group_index"]) == {0, 1}


def test_gene_level_depth_adjusted_results_match_report():
    summary = pd.read_csv("results/h3-human-vector-gene-driver-summary.csv")
    depth = summary[summary["model"] == "age_plus_sex_depth"].set_index("cell_type")
    assert depth.loc["acinar cell", "n_fdr"] == 0
    assert depth.loc["nasal mucosa goblet cell", "n_fdr"] == 0
    assert depth.loc["GABAergic neuron", "n_fdr"] == 5
    assert depth.loc["VIP GABAergic cortical interneuron", "n_fdr"] == 294
    assert depth.loc["pvalb GABAergic cortical interneuron", "n_fdr"] == 333
    assert depth.loc["sst GABAergic cortical interneuron", "n_fdr"] == 123
    shared = pd.read_csv("results/h3-human-mssm-shared-neuronal-genes.csv")
    all_three = shared[shared["n_fdr_subtypes"] == 3]
    assert len(all_three) == 7
    assert all_three["consistent_direction"].all()


def test_lung_assay_models_are_identifiable_and_null():
    assay = pd.read_csv("results/h3-human-lung-assay-audit.csv")
    assert assay["full_rank"].all()
    adjusted = assay[assay["model"] == "depth_and_assay_adjusted"]
    assert (adjusted["age_pvalue"] > 0.05).all()


def test_neuronal_gene_replication_matches_report():
    summary = pd.read_csv("results/h3-human-neuronal-gene-replication-summary.csv").set_index("cell_type")
    assert summary.loc["VIP GABAergic cortical interneuron", "same_direction_in_hbcc"] == 272
    assert summary.loc["pvalb GABAergic cortical interneuron", "same_direction_in_hbcc"] == 295
    assert summary.loc["sst GABAergic cortical interneuron", "same_direction_in_hbcc"] == 118
    assert summary.loc["VIP GABAergic cortical interneuron", "both_fdr"] == 24
    assert summary.loc["pvalb GABAergic cortical interneuron", "both_fdr"] == 20
    assert summary.loc["sst GABAergic cortical interneuron", "both_fdr"] == 62
    details = pd.read_csv("results/h3-human-neuronal-gene-replication-significant.csv")
    seven = ["GABRB1", "RNLS", "STAG1", "TTC28", "CAMK2N1", "OLFM1", "RNF152"]
    shared = details[details["gene"].isin(seven)]
    assert len(shared) == 21
    assert shared["same_direction"].all()
    rnls = shared[shared["gene"] == "RNLS"]
    assert (rnls["q_bh_hbcc"] < 0.05).all()


def test_bayesian_vector_summary_preserves_uncertainty():
    posterior = pd.read_csv("results/h3-cross-species-bayesian-posterior.csv")
    species = posterior[posterior["scope"] == "species"]
    mouse = species[(species["species"] == "Mouse") & (species["model"] == "depth_adjusted")].iloc[0]
    assert mouse["hdi94_low"] > 0
    human = species[species["model"] == "global_depth_adjusted"]
    assert len(human) == 3
    assert ((human["hdi94_low"] < 0) & (human["hdi94_high"] > 0)).all()
    nmr = species[species["species"] == "Naked mole-rat"].iloc[0]
    assert nmr["hdi94_low"] < 0 < nmr["hdi94_high"]
    diagnostics = pd.read_csv("results/h3-cross-species-bayesian-diagnostics.csv")
    assert (diagnostics["divergences"] == 0).all()
    assert (diagnostics["max_rhat"] < 1.02).all()
