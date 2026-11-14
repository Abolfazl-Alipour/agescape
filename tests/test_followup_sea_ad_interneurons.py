import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.run_followup_sea_ad_interneurons import (
    FAMILY_ORDER,
    REPO,
    checked_design,
    load_plan,
    load_targets,
    permute_within_blocks,
    sparse_counts,
)


def test_frozen_plan_and_source_gene_sets():
    plan = load_plan()
    assert plan["frozen_before_expression_outcomes"] is True
    assert plan["raw_layer"] == "UMIs"
    assert plan["donor_gate"]["minimum_donors_per_subtype"] == 10
    assert "Hellinger displacement" in plan["forbidden"]
    for short, expected in (("VIP", 76), ("pvalb", 66), ("sst", 249)):
        _, both, discovery, _, p1 = load_targets(short)
        assert len(both) == expected
        assert len(discovery) == {"VIP": 202, "pvalb": 143, "sst": 375}[short]
        assert both < discovery
        assert len(p1) == 118
    assert len(FAMILY_ORDER) == 6


def test_raw_count_validation_rejects_normalized_or_negative_values():
    matrix = np.array([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]])
    converted = sparse_counts(matrix)
    assert converted.sum() == 6
    with pytest.raises(AssertionError):
        sparse_counts(np.array([[0.1, 2.0]]))
    with pytest.raises(AssertionError):
        sparse_counts(np.array([[-1.0, 2.0]]))


def test_block_permutation_preserves_each_block_multiset():
    values = np.arange(12)
    blocks = np.repeat(np.arange(3), 4)
    permuted = permute_within_blocks(values, blocks, np.random.default_rng(42))
    for block in np.unique(blocks):
        assert np.array_equal(
            np.sort(values[blocks == block]), np.sort(permuted[blocks == block])
        )


def test_design_floor_and_rank_guards():
    metadata = pd.DataFrame({
        "age": np.linspace(70, 100, 15),
        "sex_male": np.tile([0, 1], 8)[:15],
        "mean_log_umi": np.linspace(8, 9, 15) + np.sin(np.arange(15)) * 0.01,
        "multiome_fraction": np.repeat([0, 1], [10, 5]),
        "adnc_low": np.repeat([0, 1], [8, 7]),
    })
    active, dof = checked_design(metadata, list(metadata.columns))
    assert len(active) == 5 and dof == 9
    with pytest.raises(ValueError, match="rank-deficient"):
        checked_design(metadata.assign(age=metadata["adnc_low"]), list(metadata.columns))
    with pytest.raises(ValueError, match="floor"):
        checked_design(metadata.iloc[:9], list(metadata.columns))


def test_public_outputs_are_aggregate_and_complete():
    eligibility = pd.read_csv(REPO / "results" / "followup-sea-ad-eligibility.csv")
    assert eligibility["eligible_donors_20_nuclei"].eq(15).all()
    assert eligibility["not_ad_donors"].eq(8).all()
    assert eligibility["raw_count_layer"].eq("UMIs").all()
    assert "donor" not in eligibility.columns

    genome = pd.read_csv(REPO / "results" / "followup-sea-ad-genome-wide.csv")
    assert genome["n_donors"].eq(15).all()
    assert genome["n_shared_genes"].gt(1000).all()
    assert genome["block_perm_p_greater"].between(0, 1).all()

    targeted = pd.read_csv(REPO / "results" / "followup-sea-ad-replication-summary.csv")
    both = targeted[targeted["family"] == "psychad_hbcc_both_fdr"].set_index("subtype")
    assert both["n_source_genes"].to_dict() == {"VIP": 76, "pvalb": 66, "sst": 249}
    assert targeted["n_tested"].le(targeted["n_source_genes"]).all()
    assert targeted["same_direction_fraction"].dropna().between(0, 1).all()

    genes = pd.read_csv(REPO / "results" / "followup-sea-ad-gene-effects.csv")
    assert len(genes) > 5000
    assert genes["q_bh"].between(0, 1).all()
    assert "donor" not in genes.columns
    sensitivities = pd.read_csv(REPO / "results" / "followup-sea-ad-sensitivities.csv")
    not_ad = sensitivities[sensitivities["arm"] == "not_ad_only_wls"]
    assert not_ad["status"].eq("untestable").all()
    assert not_ad["n_donors"].eq(8).all()


def test_manifest_has_verified_inputs_and_no_public_donor_expression():
    manifest = json.loads((REPO / "data" / "manifest" / "sea-ad-mtg-2024-02-13.json").read_text())
    assert len(manifest["sources"]) == 16
    assert all(len(source["sha256"]) == 64 for source in manifest["sources"])
    assert manifest["raw_count_layer"] == "UMIs"
    assert manifest["normalized_X_used"] is False
    assert manifest["public_outputs_contain_donor_expression"] is False
