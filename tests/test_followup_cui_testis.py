import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from src.followup_cui_testis import decade, summarize_metric, validate_raw_counts

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
PREFIX = "followup-cui-testis"


def test_count_validation_and_descriptive_summaries():
    matrix = sp.csr_matrix(np.array([[0, 1, 2], [3, 0, 4]], dtype=np.int64))
    checked = validate_raw_counts(matrix)
    assert checked["raw_umi_sum"] == 10
    assert checked["stored_nonzero_entries"] == 4
    with pytest.raises(ValueError, match="Non-integer"):
        validate_raw_counts(sp.csr_matrix(np.array([[1.5]])))
    with pytest.raises(ValueError, match="Negative"):
        validate_raw_counts(sp.csr_matrix(np.array([[-1]])))
    summary = summarize_metric(np.array([1, 2, 3, 4]), "umi")
    assert summary["umi_median"] == 2.5
    assert summary["umi_min"] == 1 and summary["umi_max"] == 4
    assert decade(39) == "30s" and decade(40) == "40s"


def test_25_geo_donor_and_35_atlas_mapping_is_explicit():
    frame = pd.read_csv(RESULTS / f"{PREFIX}-donor-map.csv")
    included = frame[frame["included_in_inferential_n"] == 1]
    assert len(frame) == 35 and len(included) == 25
    assert included["gsm"].nunique() == included["donor"].nunique() == 25
    assert included["library"].nunique() == 25
    external = frame[frame["scope"] == "external_atlas_only"]["donor"].tolist()
    assert external == [
        "Donor1", "Donor2", "Donor3", "Donor26", "Donor27", "Donor28",
        "Donor29", "Donor31", "Donor32", "Donor33",
    ]
    donor8 = included[included["donor"] == "Donor8"].iloc[0]
    assert donor8["primary_age"] == 40 and donor8["geo_age"] == 39
    assert included.groupby("primary_decade").size().to_dict() == {
        "30s": 4, "40s": 10, "50s": 8, "60s": 3,
    }
    manifest = json.loads((REPO / "data" / "manifest" / "GSE254315.json").read_text())
    assert manifest["mapping_resolution"]["maximum_independent_n"] == 25
    assert not manifest["mapping_resolution"]["pooled_or_demultiplexed_libraries"]


def test_frozen_plan_and_all_cell_resolved_endpoints_are_blocked():
    plan_path = REPO / "data" / "manifest" / "cui_testis_analysis_plan.json"
    plan = json.loads(plan_path.read_text())
    assert plan["frozen_after_sample_donor_mapping_and_before_expression_download_or_outcomes"]
    assert plan["scope"]["maximum_independent_n"] == 25
    assert plan["contrasts_if_annotations_become_available"]["primary"]["name"] == "ordinal_decade_depth"
    sensitivity_names = {
        row["name"] for row in plan["contrasts_if_annotations_become_available"]["sensitivities"]
    }
    assert {"old_young_depth", "exact_age_depth", "geo_donor8_39", "spermatogenesis_adjusted"} <= sensitivity_names
    status = pd.read_csv(RESULTS / f"{PREFIX}-endpoint-status.csv")
    assert len(status) == 56
    assert status["status"].eq("blocked").all()
    assert status["reason"].eq("no_released_donor_barcode_cell_type_crosswalk").all()
    assert status[["effect", "ci95_low", "ci95_high", "pvalue", "q_bh"]].isna().all().all()
    assert status["caveat"].str.contains("germ-cell-stage composition").all()


def test_public_inventory_qc_and_resource_contract():
    inventory = pd.read_csv(RESULTS / f"{PREFIX}-file-inventory.csv")
    assert len(inventory) == 75
    assert inventory["download_verified"].all()
    assert inventory["sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert inventory.groupby("donor")["kind"].nunique().eq(3).all()
    qc = pd.read_csv(RESULTS / f"{PREFIX}-qc.csv")
    assert len(qc) == 25 and qc["donor"].nunique() == 25
    assert qc["all_stored_counts_finite_nonnegative_integer"].all()
    assert not qc["annotation_crosswalk_available"].any()
    assert qc["analysis"].eq("descriptive_raw_count_qc_only_no_age_outcome_model").all()
    assert (qc["paper_qc_pass_cells"] + qc["paper_qc_fail_cells"] == qc["raw_barcodes"]).all()
    validation = json.loads((RESULTS / f"{PREFIX}-validation.json").read_text())
    assert validation["raw_barcodes_total"] == qc["raw_barcodes"].sum()
    assert not validation["raw_barcode_count_matches_geo_summary"]
    assert not validation["published_cell_count_reconciled"]
    assert validation["geo_summary_minus_paper_qc_pass"] == 99130 - qc["paper_qc_pass_cells"].sum()
    assert not validation["expression_outcomes_computed"]
    assert not validation["age_associations_computed"]
    resources = json.loads((RESULTS / f"{PREFIX}-resource-use.json").read_text())
    assert resources["guard_memory_gb"] == 8 and resources["guard_min_free_gb"] == 30
    assert resources["max_rss_kib"] < 8 * 1024**2


def test_spermatogenesis_information_is_not_imputed():
    frame = pd.read_csv(RESULTS / f"{PREFIX}-spermatogenesis-metadata.csv")
    assert len(frame) == 25
    assert frame["incomplete_spermatogenesis_status"].eq("not_released").all()
    assert frame["donor_pas_sperm_count_per_tubule"].eq("not_released").all()
    assert frame["released_freq_status"].value_counts().to_dict() == {
        "released": 24, "missing_from_Fig6a": 1,
    }
    missing = frame[frame["released_freq_status"] == "missing_from_Fig6a"]
    assert missing["donor"].tolist() == ["Donor34"]
    assert missing["released_elongating_spermatid_freq"].isna().all()
    report = (REPO / "reports" / "followup-cui-testis.md").read_text()
    assert "no public donor-and-barcode cell-type crosswalk" in report.lower()
    assert "P1" in report and "Hellinger" in report and "not run" in report
    assert "whole-body aging" in report
