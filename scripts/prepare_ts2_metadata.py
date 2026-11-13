#!/usr/bin/env python3
"""Build the TS2 supplementary-table metadata gate without reading H5AD files."""

from itertools import combinations
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / f"{REPO.name}-data" / "ts2"
RESULTS = REPO / "results"
COUNTS_SOURCE = DATA / "media-18.xlsx"
DONOR_SOURCE = DATA / "media-24.xlsx"
COUNT_COLUMNS = ["tissue", "compartment", "cell_ontology_class", "method", "donor"]
VERIFICATION = "supplementary_tables_only_no_H5AD_verification"


def donor_sort_key(donor_id):
    return int(donor_id.removeprefix("TSP"))


def read_counts():
    sheets = pd.read_excel(
        COUNTS_SOURCE,
        sheet_name=["all", "tissue_donor_cell-type", "tissue", "donor"],
    )
    counts = sheets["all"]
    counts[COUNT_COLUMNS[:-1]] = counts[COUNT_COLUMNS[:-1]].ffill()
    assert counts[COUNT_COLUMNS + ["number of cells"]].notna().all().all()
    assert not counts.duplicated(COUNT_COLUMNS).any()
    assert set(counts["method"]) == {"10X", "smartseq"}
    assert counts["number of cells"].gt(0).all()
    assert counts["number of cells"].mod(1).eq(0).all()
    assert counts["donor"].nunique() == 24
    assert counts["tissue"].nunique() == 28
    counts["number of cells"] = counts["number of cells"].astype(int)

    donor_cell = sheets["tissue_donor_cell-type"]
    donor_cell[["tissue", "cell_ontology_class"]] = donor_cell[
        ["tissue", "cell_ontology_class"]
    ].ffill()
    donor_cell_keys = ["tissue", "cell_ontology_class", "donor"]
    assert donor_cell[donor_cell_keys + ["number of cells"]].notna().all().all()
    assert not donor_cell.duplicated(donor_cell_keys).any()
    expected = counts.groupby(donor_cell_keys)["number of cells"].sum().sort_index()
    observed = donor_cell.set_index(donor_cell_keys)["number of cells"].sort_index()
    pd.testing.assert_series_equal(expected, observed, check_dtype=False)

    for sheet_name, key in (("tissue", "tissue"), ("donor", "donor")):
        aggregate = sheets[sheet_name]
        assert aggregate[key].notna().all()
        assert not aggregate[key].duplicated().any()
        expected = counts.groupby(key)["number of cells"].sum().sort_index()
        observed = aggregate.set_index(key)["number of cells"].sort_index()
        pd.testing.assert_series_equal(expected, observed, check_dtype=False)

    return counts.rename(
        columns={
            "cell_ontology_class": "cell_type",
            "donor": "donor_id",
            "number of cells": "n_cells_supplementary",
        }
    )


def read_donor_metadata(count_donors):
    manual = pd.read_excel(DONOR_SOURCE, sheet_name="Manually_Curated_metadata")
    manual = manual[manual["Unnamed: 0"].astype(str).str.fullmatch(r"Donor \d+")].copy()
    manual["donor_id"] = "TSP" + manual["Unnamed: 0"].str.extract(r"(\d+)", expand=False).str.lstrip("0")
    manual = manual.rename(columns={"Age": "age_years", "Sex": "sex"})
    manual["age_years"] = pd.to_numeric(manual["age_years"], errors="raise").astype(int)
    assert manual["donor_id"].is_unique
    assert manual[["donor_id", "age_years", "sex"]].notna().all().all()
    assert set(manual["sex"]) == {"F", "M"}

    chat = pd.read_excel(DONOR_SOURCE, sheet_name="ChatTS_metadata")
    chat = chat[chat["donor_id"].notna()].copy()
    chat["donor_id"] = chat["donor_id"].str.upper()
    chat = chat.rename(columns={"Age": "chatts_age_years", "Sex": "chatts_sex"})
    chat["chatts_age_years"] = pd.to_numeric(chat["chatts_age_years"], errors="raise").astype(int)
    chat["chatts_sex"] = chat["chatts_sex"].str.upper().map({"FEMALE": "F", "MALE": "M"})
    assert chat["donor_id"].is_unique

    donors = manual[["donor_id", "age_years", "sex"]].merge(
        chat[["donor_id", "chatts_age_years", "chatts_sex"]],
        on="donor_id",
        how="outer",
        validate="one_to_one",
    )
    assert set(donors["donor_id"]) <= set(count_donors)
    donors["sex_crosssheet_agrees"] = donors["sex"].eq(donors["chatts_sex"]).astype("object")
    donors.loc[donors["sex"].isna() | donors["chatts_sex"].isna(), "sex_crosssheet_agrees"] = pd.NA
    donors["age_reconciliation"] = "matches_ChatTS_metadata"
    donors.loc[donors["age_years"].isna() | donors["chatts_age_years"].isna(), "age_reconciliation"] = "missing_from_one_media-24_sheet"
    discrepancy = donors["age_years"].notna() & donors["chatts_age_years"].notna() & donors["age_years"].ne(donors["chatts_age_years"])
    donors.loc[discrepancy, "age_reconciliation"] = "manual_curated_value_retained_ChatTS_conflicts"

    all_donors = pd.DataFrame({"donor_id": sorted(count_donors, key=donor_sort_key)})
    donors = all_donors.merge(donors, on="donor_id", how="left", validate="one_to_one")
    donors["age_reconciliation"] = donors["age_reconciliation"].fillna("absent_from_media-24")
    donors[["age_years", "chatts_age_years"]] = donors[["age_years", "chatts_age_years"]].astype("Int64")
    donors["metadata_verification"] = VERIFICATION
    return donors


def method_summary(methods):
    return ";".join(sorted(set(methods)))


def build_donor_summary(counts, donors):
    summary = counts.groupby("donor_id").agg(
        n_tissues_supplementary=("tissue", "nunique"),
        n_cells_supplementary=("n_cells_supplementary", "sum"),
        n_methods_supplementary=("method", "nunique"),
        methods_supplementary=("method", method_summary),
    ).reset_index()
    summary = donors.merge(summary, on="donor_id", how="left", validate="one_to_one")
    return summary[[
        "donor_id", "age_years", "sex", "chatts_age_years", "age_reconciliation",
        "sex_crosssheet_agrees", "n_tissues_supplementary", "n_cells_supplementary",
        "n_methods_supplementary", "methods_supplementary", "metadata_verification",
    ]]


def build_crosswalk(counts, donors):
    crosswalk = counts.groupby(["donor_id", "tissue"]).agg(
        n_cells_supplementary=("n_cells_supplementary", "sum"),
        n_cell_types_supplementary=("cell_type", "nunique"),
        methods_supplementary=("method", method_summary),
    ).reset_index()
    method_cells = counts.pivot_table(
        index=["donor_id", "tissue"], columns="method", values="n_cells_supplementary",
        aggfunc="sum", fill_value=0,
    ).reset_index().rename(columns={
        "10X": "n_cells_10x_supplementary",
        "smartseq": "n_cells_smartseq_supplementary",
    })
    crosswalk = crosswalk.merge(method_cells, on=["donor_id", "tissue"], validate="one_to_one")
    crosswalk = crosswalk.merge(donors[["donor_id", "age_years", "sex"]], on="donor_id", validate="many_to_one")
    crosswalk["metadata_verification"] = VERIFICATION
    return crosswalk.sort_values(["tissue", "donor_id"], key=lambda s: s.map(donor_sort_key) if s.name == "donor_id" else s).reset_index(drop=True)


def build_organ_summary_from_counts(counts, crosswalk):
    summary = counts.groupby("tissue").agg(
        n_cells_supplementary=("n_cells_supplementary", "sum"),
        n_cell_types_supplementary=("cell_type", "nunique"),
        methods_supplementary=("method", method_summary),
    ).reset_index()
    donor_summary = crosswalk.groupby("tissue").agg(
        n_donors_supplementary=("donor_id", "nunique"),
        n_donors_age_verified=("age_years", "count"),
        age_min=("age_years", "min"),
        age_max=("age_years", "max"),
        age_sd=("age_years", "std"),
    ).reset_index()
    summary = summary.merge(donor_summary, on="tissue", validate="one_to_one")
    summary["metadata_verification"] = VERIFICATION
    return summary[[
        "tissue", "n_donors_supplementary", "n_donors_age_verified",
        "n_cells_supplementary", "age_min", "age_max", "age_sd",
        "n_cell_types_supplementary", "methods_supplementary", "metadata_verification",
    ]].sort_values("tissue").reset_index(drop=True)


def build_coverage(counts, donors):
    coverage = counts.merge(donors[["donor_id", "age_years", "sex"]], on="donor_id", validate="many_to_one")
    coverage["meets_20_cell_metadata_floor"] = coverage["n_cells_supplementary"] >= 20
    coverage["metadata_verification"] = VERIFICATION
    return coverage[[
        "donor_id", "tissue", "compartment", "cell_type", "method", "age_years", "sex",
        "n_cells_supplementary", "meets_20_cell_metadata_floor", "metadata_verification",
    ]].sort_values(["tissue", "cell_type", "method", "donor_id"], key=lambda s: s.map(donor_sort_key) if s.name == "donor_id" else s).reset_index(drop=True)


def build_organ_pairs(counts, donors):
    ages = donors.dropna(subset=["age_years"]).set_index("donor_id")["age_years"]
    rows = []
    for method in sorted(counts["method"].unique()):
        method_counts = counts[counts["method"] == method]
        tissues = sorted(method_counts["tissue"].unique())
        donors_by_tissue = method_counts.groupby("tissue")["donor_id"].agg(set)
        for organ_1, organ_2 in combinations(tissues, 2):
            shared = donors_by_tissue[organ_1] & donors_by_tissue[organ_2]
            shared_age_verified = sorted(shared & set(ages.index), key=donor_sort_key)
            shared_ages = ages.loc[shared_age_verified]
            age_sd = shared_ages.std() if len(shared_ages) >= 2 else float("nan")
            six_donors = len(shared_age_verified) >= 6
            age_sd_eligible = pd.notna(age_sd) and age_sd >= 12
            rows.append({
                "method_supplementary": method,
                "organ_1": organ_1,
                "organ_2": organ_2,
                "n_shared_donors_supplementary": len(shared),
                "n_shared_donors_age_verified": len(shared_age_verified),
                "shared_donors_age_verified": ";".join(shared_age_verified),
                "age_sd_shared": age_sd,
                "eligible_six_age_verified_donors": six_donors,
                "eligible_age_sd_12": age_sd_eligible,
                "passes_organ_level_metadata_screen": six_donors and age_sd_eligible,
                "metadata_verification": VERIFICATION,
            })
    return pd.DataFrame(rows)


def main():
    counts = read_counts()
    donors = read_donor_metadata(set(counts["donor_id"]))
    donor_summary = build_donor_summary(counts, donors)
    crosswalk = build_crosswalk(counts, donors)
    organ_summary = build_organ_summary_from_counts(counts, crosswalk)
    coverage = build_coverage(counts, donors)
    pairs = build_organ_pairs(counts, donors)

    outputs = {
        "followup-ts2-donor-summary.csv": donor_summary,
        "followup-ts2-organ-summary.csv": organ_summary,
        "followup-ts2-donor-organ-crosswalk.csv": crosswalk,
        "followup-ts2-celltype-coverage.csv": coverage,
        "followup-ts2-organ-pair-eligibility.csv": pairs,
    }
    for filename, frame in outputs.items():
        frame.to_csv(RESULTS / filename, index=False)
        print(f"wrote {RESULTS / filename}: {len(frame)} rows")
    print(f"supplementary cells reconciled: {counts['n_cells_supplementary'].sum()}")
    print(f"donors with supplementary counts: {counts['donor_id'].nunique()}")
    print(f"donors with curated ages: {donors['age_years'].notna().sum()}")
    print(f"organ-level metadata screens passed: {pairs['passes_organ_level_metadata_screen'].sum()} of {len(pairs)}")


if __name__ == "__main__":
    main()
