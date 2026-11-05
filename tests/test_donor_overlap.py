import pandas as pd

from scripts.audit_m9_m10_donor_overlap import identifier_overlap


def test_identifier_overlap_requires_matching_age_and_sex():
    public = pd.DataFrame({
        "donor_id": ["a", "b", "c", "d"],
        "age": [50, 60, 70, 80],
        "sex": ["male", "female", "male", "female"],
    })
    controlled = pd.DataFrame({
        "anonymID": ["a", "b", "c", "z"],
        "age": [50, 65, 70, 80],
        "sex": ["male", "female", "female", "female"],
    })
    exact_strings, exact_concordant = identifier_overlap(public, controlled)
    assert exact_strings == 3
    assert exact_concordant == 1


def test_committed_overlap_counts_are_nonempty_and_bounded():
    result = pd.read_csv("results/m9-m10-donor-overlap-audit.csv").set_index("cohort")
    assert result.loc["MSSM_DLPFC", "exact_identifier_overlap"] == 103
    assert result.loc["HBCC_DLPFC", "exact_identifier_overlap"] == 131
    for row in result.itertuples():
        assert 0 <= row.exact_identifier_overlap <= row.exact_identifier_string_overlap
        assert row.exact_identifier_string_overlap <= min(row.m9_donors, row.m10_donors_in_bank)
