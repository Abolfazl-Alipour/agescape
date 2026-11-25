#!/usr/bin/env python3
"""Audit source sufficiency and exact-name matches for other-species P1m tests."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_afca_aging import afca_obs
from scripts.run_gao_worm_aging import normal_obs
from scripts.run_mature_matched_p1m import human_roster, normalize_name
from scripts.run_zebrafish_testis_aging import standard_obs

RESULTS = REPO / "results"
DATA_ROOT = REPO.parent / "agescape-revision-data"
CELL_TYPE_CSV = RESULTS / "followup-mature-other-species-pergene-cell-types.csv"
SUMMARY_CSV = RESULTS / "followup-mature-other-species-pergene-sufficiency.csv"


def reference_rosters() -> tuple[pd.DataFrame, pd.DataFrame]:
    mouse = pd.read_csv(RESULTS / "followup-mature-p1m-eligibility.csv")
    mouse = mouse[mouse["eligible"].astype(str).eq("yes")].copy()
    mouse["target_tissue"] = mouse["tissue"].replace({"Marrow": "Blood"})
    mouse["normalized"] = mouse["cell_ontology_class"].map(normalize_name)
    assert len(mouse) == 36
    human = human_roster()
    human = human[human["source_status"].eq("present")].copy()
    return mouse, human


def matches(frame: pd.DataFrame, tissue: str, normalized: str) -> tuple[str, str]:
    same_tissue = frame[(frame["target_tissue"] == tissue) & (frame["normalized"] == normalized)]
    any_tissue = frame[frame["normalized"] == normalized]
    return ";".join(sorted(same_tissue["match_label"].unique())), ";".join(
        sorted(any_tissue["match_label"].unique())
    )


def audit_counts(
    rows: list[dict],
    obs: pd.DataFrame,
    *,
    species: str,
    cohort: str,
    tissue: str,
    scope: str,
    unit_column: str,
    cell_column: str,
    source: str,
    raw_status: str,
    mouse: pd.DataFrame,
    human: pd.DataFrame,
) -> None:
    required = [unit_column, cell_column, "age"]
    if obs.empty or any(column not in obs for column in required):
        raise ValueError(f"Incomplete source metadata: {species} {scope}")
    selected = obs.dropna(subset=required).copy()
    selected[unit_column] = selected[unit_column].astype(str)
    selected[cell_column] = selected[cell_column].astype(str)
    if "n_cells" in selected:
        selected = selected[selected["n_cells"] > 0].copy()
    unit_ages = selected.groupby(unit_column, observed=True)["age"].agg(
        lambda values: sorted(set(values.astype(float)))
    )
    assert unit_ages.map(len).eq(1).all(), f"Non-unique unit age: {species} {scope}"
    for cell_type, part in selected.groupby(cell_column, observed=True, sort=True):
        grouped = part.groupby(unit_column, observed=True)
        per_unit = (grouped["n_cells"].sum() if "n_cells" in part else grouped.size()).sort_index()
        ages = sorted(part.groupby(unit_column, observed=True)["age"].first().unique())
        normalized = normalize_name(cell_type)
        mouse_same, mouse_any = matches(mouse, tissue, normalized)
        human_same, human_any = matches(human, tissue, normalized)
        rows.append({
            "species": species,
            "cohort": cohort,
            "tissue": tissue,
            "scope": scope,
            "cell_type": str(cell_type),
            "normalized_cell_type": normalized,
            "inferential_unit": unit_column,
            "n_units": len(per_unit),
            "unit_ages": ";".join(f"{age:g}" for age in ages),
            "n_age_levels": len(ages),
            "n_cells": int(per_unit.sum()),
            "n_units_ge20_cells": int((per_unit >= 20).sum()),
            "min_cells_per_unit": int(per_unit.min()),
            "max_cells_per_unit": int(per_unit.max()),
            "slope_fittable_min2_units": len(per_unit) >= 2 and len(ages) >= 2,
            "passes_matched_floor_min3_units_50cells": len(per_unit) >= 3 and per_unit.sum() >= 50,
            "mouse_exact_same_tissue_matches": mouse_same,
            "human_exact_same_tissue_matches": human_same,
            "mouse_exact_name_any_tissue_matches": mouse_any,
            "human_exact_name_any_tissue_matches": human_any,
            "mouse_pergene_decision": "eligible_match" if mouse_same else "blocked_no_exact_name_match",
            "human_pergene_decision": "eligible_match" if human_same else "blocked_no_exact_name_match",
            "source": source,
            "raw_status": raw_status,
        })


def summarize(audit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["species", "cohort", "tissue", "scope", "inferential_unit", "source", "raw_status"]
    for key, part in audit.groupby(keys, sort=True, dropna=False):
        record = dict(zip(keys, key))
        ages = sorted({float(age) for text in part["unit_ages"] for age in text.split(";")})
        record.update({
            "n_cell_types": len(part),
            "n_slope_fittable_cell_types": int(part["slope_fittable_min2_units"].sum()),
            "n_floored_cell_types": int(part["passes_matched_floor_min3_units_50cells"].sum()),
            "min_units_per_cell_type": int(part["n_units"].min()),
            "max_units_per_cell_type": int(part["n_units"].max()),
            "total_cell_type_assignments": int(part["n_cells"].sum()),
            "unit_ages": ";".join(f"{age:g}" for age in ages),
            "n_mouse_same_tissue_matches": int(part["mouse_exact_same_tissue_matches"].ne("").sum()),
            "n_human_same_tissue_matches": int(part["human_exact_same_tissue_matches"].ne("").sum()),
            "n_mouse_name_matches_any_tissue": int(part["mouse_exact_name_any_tissue_matches"].ne("").sum()),
            "n_human_name_matches_any_tissue": int(part["human_exact_name_any_tissue_matches"].ne("").sum()),
            "mouse_test_status": "NOT_RUN_NO_MATCHED_PAIRS",
            "human_test_status": "NOT_RUN_NO_MATCHED_PAIRS",
        })
        rows.append(record)
    return pd.DataFrame(rows)


def main() -> None:
    mouse, human = reference_rosters()
    mouse["match_label"] = mouse["unit"]
    human["match_label"] = human["human_cohort"] + "::" + human["human_unit"]
    human = human.rename(columns={"human_tissue": "target_tissue"})
    rows = []

    for organ in ("head", "body"):
        path = REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad"
        atlas = ad.read_h5ad(path, backed="r")
        try:
            obs = afca_obs(atlas)
            obs = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated") & (obs["age"] >= 30)]
            audit_counts(
                rows, obs, species="Fruit fly", cohort="AFCA_GSE218661", tissue=organ,
                scope="mature_30plus_days", unit_column="sample", cell_column="cell_type",
                source=str(path), raw_status=f"present_resolved={path.resolve().exists()}", mouse=mouse, human=human,
            )
        finally:
            atlas.file.close()

    path = REPO / "data" / "interim" / "worm-gao-integrated.h5ad"
    atlas = ad.read_h5ad(path, backed="r")
    try:
        obs = normal_obs(atlas)
        for scope, selected in (("full_1plus_days", obs), ("mature_6plus_days", obs[obs["age"] >= 6])):
            audit_counts(
                rows, selected, species="C. elegans", cohort="Gao_GSE229022", tissue="Whole animal",
                scope=scope, unit_column="donor", cell_column="cell_type", source=str(path),
                raw_status=f"present_resolved={path.resolve().exists()}", mouse=mouse, human=human,
            )
    finally:
        atlas.file.close()

    path = REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad"
    atlas = ad.read_h5ad(path, backed="r")
    try:
        obs = standard_obs(atlas, "5-22mo")
        obs = obs[obs["fine_cell_type"] != "unknown"]
        for scope, selected in (("full_5_to_22_months", obs), ("mature_12plus_months", obs[obs["age"] >= 12])):
            audit_counts(
                rows, selected, species="Zebrafish", cohort="GSE275361", tissue="Testis",
                scope=scope, unit_column="donor", cell_column="fine_cell_type", source=str(path),
                raw_status=f"present_resolved={path.resolve().exists()}", mouse=mouse, human=human,
            )
    finally:
        atlas.file.close()

    nmr = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    nmr = nmr[nmr["species"] == "Naked mole-rat"].rename(columns={"donor": "animal"})
    nmr_skin_raw = REPO.parent / "agescape" / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz"
    audit_counts(
        rows, nmr, species="Naked mole-rat", cohort="GSE152886", tissue="Skin",
        scope="all_1.1_and_11.3_years_no_defensible_mature_window", unit_column="animal",
        cell_column="cell_type", source=str(nmr_skin_raw),
        raw_status=f"present_elsewhere_on_beast={nmr_skin_raw.exists()};committed_donor_table_used_for_counts",
        mouse=mouse, human=human,
    )

    donor_path = RESULTS / "followup-nmr-spleen-t-gse214390-donors.csv"
    donors = pd.read_csv(donor_path)
    spleen_raw = DATA_ROOT / "tmp_organs" / "nmr_gse214390"
    for cell_type in ("CD4.T.cell", "CD8.T.cell", "non.cytotoxic.gdT.cell", "cytotoxic.gdT.cell", "ILC"):
        obs = donors[["sample_id", "age_years"]].rename(columns={"sample_id": "animal", "age_years": "age"})
        obs["cell_type"] = cell_type
        obs["n_cells"] = donors[f"cells_{cell_type}"].astype(int)
        audit_counts(
            rows, obs, species="Naked mole-rat", cohort="GSE214390", tissue="Spleen",
            scope="all_2_and_24_years_no_defensible_mature_window", unit_column="animal",
            cell_column="cell_type", source=str(spleen_raw),
            raw_status=f"raw_directory_absent={not spleen_raw.exists()};committed_donor_table_used_for_counts",
            mouse=mouse, human=human,
        )

    audit = pd.DataFrame(rows)
    summary = summarize(audit)
    assert not audit["mouse_exact_same_tissue_matches"].ne("").any()
    assert not audit["human_exact_same_tissue_matches"].ne("").any()
    assert not audit["mouse_exact_name_any_tissue_matches"].ne("").any()
    assert not audit["human_exact_name_any_tissue_matches"].ne("").any()
    audit.to_csv(CELL_TYPE_CSV, index=False)
    summary.to_csv(SUMMARY_CSV, index=False)
    reread = pd.read_csv(CELL_TYPE_CSV).fillna("")
    assert len(reread) == len(audit)
    assert not reread["mouse_exact_same_tissue_matches"].ne("").any()
    assert not reread["human_exact_same_tissue_matches"].ne("").any()
    print(summary[[
        "species", "scope", "n_cell_types", "n_slope_fittable_cell_types", "n_floored_cell_types",
        "min_units_per_cell_type", "max_units_per_cell_type", "unit_ages", "raw_status",
    ]].to_string(index=False))
    print(f"Wrote {CELL_TYPE_CSV} ({len(audit)} rows)")
    print(f"Wrote {SUMMARY_CSV} ({len(summary)} rows)")


if __name__ == "__main__":
    main()
