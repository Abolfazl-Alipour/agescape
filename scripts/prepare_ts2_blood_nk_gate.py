#!/usr/bin/env python3
"""Freeze the TS2 blood NK plan and record the CELLxGENE access gate.

This script queries collection metadata only. It never downloads or opens an
H5AD and never computes expression outcomes. The expected blood asset is
larger than the 2 GiB task ceiling, so the expression arm remains blocked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
MANIFEST = REPO / "data" / "manifest" / "tabula_sapiens_2_0_blood.json"
COVERAGE = RESULTS / "followup-ts2-celltype-coverage.csv"
COLLECTION_ID = "e5f58829-1a66-40b5-a624-9046778e74f5"
COLLECTION_URL = (
    "https://api.cellxgene.cziscience.com/curation/v1/collections/"
    + COLLECTION_ID
)
DATASET_ID = "983d5ec9-40e8-4512-9e65-a572a9c486cb"
DATASET_VERSION_ID = "1bf55482-ad49-446a-87c6-f3442bebcef1"
GIB = 1024 ** 3
MAX_FILE_GIB = 2
MIN_FREE_GB = 30
EXPECTED_DONORS = {
    "TSP1", "TSP2", "TSP7", "TSP8", "TSP10", "TSP14", "TSP21", "TSP25", "TSP27",
}
EXPECTED_CANDIDATES = {"TSP1", "TSP2", "TSP7", "TSP10", "TSP14", "TSP21", "TSP25"}
BLOCKER = "blood_H5AD_exceeds_2_GiB_download_cap"
H5AD_VERIFICATION = "not_run_blood_H5AD_exceeds_2_GiB_download_cap"


def fetch_collection(url=COLLECTION_URL):
    request = urllib.request.Request(url, headers={"User-Agent": "agescape-ts2-gate/1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def labels(records):
    return ";".join(sorted(record["label"] for record in records or []))


def resolve_blood_asset(collection):
    assert collection["collection_id"] == COLLECTION_ID
    matches = [
        dataset for dataset in collection["datasets"]
        if dataset["dataset_id"] == DATASET_ID and dataset["title"] == "Tabula Sapiens - Blood"
    ]
    assert len(matches) == 1, "expected exactly one authenticated blood dataset"
    dataset = matches[0]
    assert dataset["dataset_version_id"] == DATASET_VERSION_ID, (
        "blood dataset version changed; review the new asset instead of reusing this gate"
    )
    assert [term["label"] for term in dataset["tissue"]] == ["blood"]
    assets = [asset for asset in dataset["assets"] if asset["filetype"] == "H5AD"]
    assert len(assets) == 1, "expected exactly one blood H5AD asset"
    asset = assets[0]
    assert asset["url"].endswith(f"/{DATASET_VERSION_ID}.h5ad")
    size = int(asset["filesize"])
    assert size > 0
    return collection, dataset, asset, size


def within_file_cap(size_bytes):
    return int(size_bytes) <= MAX_FILE_GIB * GIB


def build_asset_row(collection, dataset, asset, size):
    return {
        "collection_id": collection["collection_id"],
        "collection_version_id": collection["collection_version_id"],
        "collection_revised_at": collection["revised_at"],
        "dataset_id": dataset["dataset_id"],
        "dataset_version_id": dataset["dataset_version_id"],
        "dataset_title": dataset["title"],
        "dataset_published_at": dataset["published_at"],
        "dataset_revised_at": dataset["revised_at"],
        "cell_count_catalog": dataset["cell_count"],
        "asset_filetype": asset["filetype"],
        "asset_url": asset["url"],
        "asset_size_bytes": size,
        "asset_size_gib": size / GIB,
        "maximum_file_gib": MAX_FILE_GIB,
        "minimum_free_gb": MIN_FREE_GB,
        "within_file_cap": within_file_cap(size),
        "catalog_donor_ids": ";".join(sorted(dataset["donor_id"], key=donor_sort_key)),
        "catalog_assays": labels(dataset.get("assay")),
        "catalog_diseases": labels(dataset.get("disease")),
        "catalog_development_stages": labels(dataset.get("development_stage")),
        "catalog_suspension_types": ";".join(sorted(dataset.get("suspension_type", []))),
        "H5AD_verification": H5AD_VERIFICATION,
    }


def donor_sort_key(donor_id):
    return int(donor_id.removeprefix("TSP"))


def build_donor_eligibility(coverage):
    selected = coverage[
        coverage["tissue"].eq("Blood")
        & coverage["compartment"].eq("Immune")
        & coverage["cell_type"].eq("natural killer cell")
        & coverage["method"].eq("10X")
    ].copy()
    assert set(selected["donor_id"]) == EXPECTED_DONORS
    assert not selected["donor_id"].duplicated().any()
    floor = selected["meets_20_cell_metadata_floor"].astype(str).str.lower().eq("true")
    assert floor.equals(selected["n_cells_supplementary"].ge(20))
    selected["supplementary_candidate"] = floor
    candidates = selected[selected["supplementary_candidate"]]
    assert set(candidates["donor_id"]) == EXPECTED_CANDIDATES
    assert candidates["age_years"].notna().all()
    assert candidates["sex"].isin(["F", "M"]).all()
    assert candidates["age_years"].std(ddof=1) >= 12
    selected["H5AD_metadata_verification"] = H5AD_VERIFICATION
    selected["raw_counts_verified"] = False
    selected["primary_data_deduplication_verified"] = False
    selected["expression_eligible"] = False
    selected["status"] = "blocked_before_H5AD_verification"
    selected.loc[~selected["supplementary_candidate"], "status"] = (
        "excluded_supplementary_20_cell_floor"
    )
    selected["blocker"] = BLOCKER
    columns = list(coverage.columns) + [
        "supplementary_candidate", "H5AD_metadata_verification", "raw_counts_verified",
        "primary_data_deduplication_verified", "expression_eligible", "status", "blocker",
    ]
    return selected[columns].sort_values(
        "donor_id", key=lambda series: series.map(donor_sort_key)
    ).reset_index(drop=True)


def build_endpoint_status():
    common_reason = (
        "The sole blood H5AD exceeds the 2 GiB ceiling; raw counts, primary-data "
        "deduplication, per-cell donor/age/sex/assay/disease fields, and endpoint floors "
        "cannot be verified. No outcomes were examined."
    )
    designs = {
        "donor_pseudobulk_gene_models": (
            "Exploratory seven-donor mean cell-level log1p(CPM) per gene; continuous-age "
            "OLS with effect per 10 years and 95% CI; sex, mean log actual UMI, and "
            "nonconstant assay terms only when the design is full rank with positive "
            "residual df; all variable genes with finite donor values; BH across all "
            "successfully fitted genes; no cell-level inference."
        ),
        "P1_transport": (
            "Existing P1_Clean up/down signature and expression-matched random-null "
            "procedure only after nonempty detected up/down sets and valid null pools "
            "are verified; donor-level continuous-age model only if covariates are identifiable."
        ),
        "off_program_expression": (
            "Existing off-program definition and technical-gene exclusions only after "
            "its gene/signature floors are verified; donor-level continuous-age model "
            "only if covariates are identifiable."
        ),
        "Hellinger_displacement": (
            "Existing young-reference Hellinger endpoint with leave-one-donor-out "
            "reference, at least 50 remaining reference cells and at least 20 cells per "
            "donor; donor-level inference only if all reference and covariate floors pass."
        ),
    }
    return pd.DataFrame([
        {
            "endpoint": endpoint,
            "inferential_unit": "donor",
            "donor_family": "TSP1;TSP2;TSP7;TSP10;TSP14;TSP21;TSP25",
            "donor_family_status": "seven_donor_exploratory_supplementary_screen_only",
            "frozen_design": design,
            "status": "untestable",
            "reason": common_reason,
            "outcomes_examined": False,
        }
        for endpoint, design in designs.items()
    ])


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_manifest(collection, dataset, asset, size, donors):
    candidates = donors.loc[donors["supplementary_candidate"], "donor_id"].tolist()
    return {
        "study": "Tabula Sapiens 2.0 blood natural-killer-cell arm",
        "collection_api_url": COLLECTION_URL,
        "collection_id": collection["collection_id"],
        "collection_version_id": collection["collection_version_id"],
        "collection_revised_at": collection["revised_at"],
        "dataset": {
            "dataset_id": dataset["dataset_id"],
            "dataset_version_id": dataset["dataset_version_id"],
            "title": dataset["title"],
            "published_at": dataset["published_at"],
            "revised_at": dataset["revised_at"],
            "cell_count_catalog": dataset["cell_count"],
            "donor_ids_catalog": dataset["donor_id"],
            "assay_catalog": dataset.get("assay", []),
            "disease_catalog": dataset.get("disease", []),
            "development_stage_catalog": dataset.get("development_stage", []),
            "tissue_catalog": dataset["tissue"],
            "assets": [{
                "filetype": asset["filetype"],
                "url": asset["url"],
                "filesize_bytes": size,
                "filesize_gib": size / GIB,
            }],
            "catalog_metadata_is_not_H5AD_verification": True,
        },
        "access_gate": {
            "maximum_file_gib": MAX_FILE_GIB,
            "minimum_free_gb": MIN_FREE_GB,
            "within_file_cap": within_file_cap(size),
            "download_permitted": False,
            "local_H5AD": None,
            "download_command_if_policy_changes": (
                "python scripts/resource_guard.py download --url " + asset["url"]
                + " --output /home/reza/projects/agescape-revision-data/ts2/blood.h5ad"
                + " --min-free-gb 30 --max-file-gb 2"
            ),
        },
        "supplementary_coverage": {
            "path": str(COVERAGE.relative_to(REPO)),
            "sha256": sha256(COVERAGE),
        },
        "donor_gate": {
            "method": "10X",
            "tissue": "Blood",
            "cell_type": "natural killer cell",
            "minimum_cells_per_donor": 20,
            "candidate_donors": candidates,
            "n_candidate_donors": len(candidates),
            "analysis_label": "exploratory",
            "candidate_ages_are_H5AD_verified": False,
        },
        "verification": {
            "raw_counts": False,
            "primary_data_and_deduplication": False,
            "per_cell_donor_age_sex_assay_disease": False,
            "P1_signature_floors": False,
            "off_program_floors": False,
            "Hellinger_reference_floors": False,
        },
        "expression_analysis": {
            "status": "blocked_before_H5AD_access",
            "outcomes_examined": False,
            "reason": BLOCKER,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection-json", type=Path, help="offline collection API response")
    args = parser.parse_args()
    if args.collection_json:
        collection = json.loads(args.collection_json.read_text())
    else:
        collection = fetch_collection()
    collection, dataset, asset, size = resolve_blood_asset(collection)
    assert not within_file_cap(size), (
        "blood asset now meets the size cap; stop and review a separate H5AD analysis plan"
    )
    donors = build_donor_eligibility(pd.read_csv(COVERAGE))
    endpoints = build_endpoint_status()
    asset_row = build_asset_row(collection, dataset, asset, size)

    outputs = {
        RESULTS / "followup-ts2-blood-nk-asset.csv": pd.DataFrame([asset_row]),
        RESULTS / "followup-ts2-blood-nk-donor-eligibility.csv": donors,
        RESULTS / "followup-ts2-blood-nk-endpoint-status.csv": endpoints,
    }
    for path, frame in outputs.items():
        frame.to_csv(path, index=False)
        print(f"wrote {path.relative_to(REPO)}: {len(frame)} rows")
    MANIFEST.write_text(json.dumps(
        build_manifest(collection, dataset, asset, size, donors), indent=2, sort_keys=True
    ) + "\n")
    print(f"wrote {MANIFEST.relative_to(REPO)}")
    print(f"BLOCKED: {size} bytes ({size / GIB:.3f} GiB) exceeds {MAX_FILE_GIB} GiB; no download or outcomes")


if __name__ == "__main__":
    main()
