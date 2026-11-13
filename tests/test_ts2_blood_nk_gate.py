import copy

import pandas as pd
import pytest

from scripts.prepare_ts2_blood_nk_gate import (
    COLLECTION_ID,
    DATASET_ID,
    DATASET_VERSION_ID,
    EXPECTED_CANDIDATES,
    GIB,
    build_donor_eligibility,
    build_endpoint_status,
    resolve_blood_asset,
    within_file_cap,
)


def collection_fixture():
    return {
        "collection_id": COLLECTION_ID,
        "collection_version_id": "collection-version",
        "revised_at": "2026-06-11T16:54:33+00:00",
        "datasets": [{
            "dataset_id": DATASET_ID,
            "dataset_version_id": DATASET_VERSION_ID,
            "title": "Tabula Sapiens - Blood",
            "tissue": [{"label": "blood"}],
            "assets": [{
                "filetype": "H5AD",
                "filesize": 2568875081,
                "url": f"https://datasets.cellxgene.cziscience.com/{DATASET_VERSION_ID}.h5ad",
            }],
        }],
    }


def test_exact_blood_asset_is_resolved_and_blocked():
    collection, dataset, asset, size = resolve_blood_asset(collection_fixture())
    assert collection["collection_id"] == COLLECTION_ID
    assert dataset["dataset_id"] == DATASET_ID
    assert asset["filetype"] == "H5AD"
    assert size == 2568875081
    assert size > 2 * GIB
    assert not within_file_cap(size)


def test_size_cap_boundary_is_inclusive():
    assert within_file_cap(2 * GIB)
    assert not within_file_cap(2 * GIB + 1)


def test_changed_dataset_version_fails_closed():
    collection = collection_fixture()
    collection["datasets"][0]["dataset_version_id"] = "new-version"
    with pytest.raises(AssertionError, match="version changed"):
        resolve_blood_asset(collection)


def test_duplicate_blood_dataset_fails_closed():
    collection = collection_fixture()
    collection["datasets"].append(copy.deepcopy(collection["datasets"][0]))
    with pytest.raises(AssertionError, match="exactly one"):
        resolve_blood_asset(collection)


def test_supplementary_donor_freeze_has_seven_candidates_and_no_expression_eligibility():
    coverage = pd.read_csv("results/followup-ts2-celltype-coverage.csv")
    donors = build_donor_eligibility(coverage)
    candidates = donors[donors["supplementary_candidate"]]
    assert len(donors) == 9
    assert set(candidates["donor_id"]) == EXPECTED_CANDIDATES
    assert candidates["age_years"].std(ddof=1) == pytest.approx(12.061351104921474)
    assert not donors["expression_eligible"].any()
    assert not donors["raw_counts_verified"].any()
    assert not donors["primary_data_deduplication_verified"].any()
    excluded = donors[~donors["supplementary_candidate"]]
    assert set(excluded["donor_id"]) == {"TSP8", "TSP27"}
    assert set(excluded["status"]) == {"excluded_supplementary_20_cell_floor"}


def test_all_expression_endpoints_are_untestable_without_outcomes():
    endpoints = build_endpoint_status()
    assert set(endpoints["endpoint"]) == {
        "donor_pseudobulk_gene_models", "P1_transport", "off_program_expression",
        "Hellinger_displacement",
    }
    assert set(endpoints["status"]) == {"untestable"}
    assert not endpoints["outcomes_examined"].any()
    assert not {"effect", "pvalue", "qvalue"}.intersection(endpoints.columns)
