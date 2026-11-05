#!/usr/bin/env python3
"""Bound M9/M10 donor overlap without exporting controlled donor identifiers."""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent


def identifier_overlap(public, controlled):
    matched = public.merge(
        controlled, left_on="donor_id", right_on="anonymID", suffixes=("_public", "_controlled")
    )
    concordant = matched[
        (matched["age_public"] == matched["age_controlled"])
        & (matched["sex_public"] == matched["sex_controlled"])
    ]
    return len(matched), len(concordant)


def main():
    m9 = pd.read_csv(REPO / "results" / "m9-donor-cell-counts-prethreshold.csv")
    m9 = m9[m9["cohort"].isin(["MSSM_DLPFC", "HBCC_DLPFC", "RADC_DLPFC", "SEAAD_DLPFC"])]
    m9 = (m9.groupby(["cohort", "donor_id"], observed=True)
          .agg(age=("age", "first"), sex=("sex", "first")).reset_index())
    m9["sex"] = m9["sex"].astype(str).str.lower()
    m10 = pd.read_csv(REPO / "data" / "external" / "dlpfc-atlas" / "SupplementaryData1_Jan2026.csv")
    m10["sex"] = m10["Sex"].astype(str).str.lower()
    m10["age"] = m10["Age_rounded"].astype(float)

    rows = []
    for cohort in ["MSSM_DLPFC", "HBCC_DLPFC", "RADC_DLPFC", "SEAAD_DLPFC"]:
        public = m9[m9["cohort"] == cohort].copy()
        bank = "MSSM" if cohort == "MSSM_DLPFC" else "HBCC"
        comparable = cohort in ("MSSM_DLPFC", "HBCC_DLPFC")
        controlled = m10[m10["Brain_bank"] == bank] if comparable else m10.iloc[0:0]
        exact_strings, exact_concordant = identifier_overlap(public, controlled)
        candidate_counts = []
        unique_controlled_ids = set()
        for donor in public.itertuples():
            candidates = controlled[(controlled["age"] == donor.age) & (controlled["sex"] == donor.sex)]
            candidate_counts.append(len(candidates))
            if len(candidates) == 1:
                unique_controlled_ids.add(candidates.iloc[0]["anonymID"])
        counts = pd.Series(candidate_counts, dtype=int)
        rows.append({
            "cohort": cohort,
            "brain_bank_compared": bank if comparable else "not applicable",
            "m9_donors": len(public),
            "m10_donors_in_bank": len(controlled),
            "exact_identifier_overlap": exact_concordant,
            "exact_identifier_string_overlap": exact_strings,
            "m9_with_any_age_sex_bank_match": int((counts > 0).sum()),
            "m9_with_unique_age_sex_bank_match": int((counts == 1).sum()),
            "m9_with_ambiguous_age_sex_bank_match": int((counts > 1).sum()),
            "m9_without_age_sex_bank_match": int((counts == 0).sum()),
            "distinct_m10_ids_in_unique_matches": len(unique_controlled_ids),
            "overlap_status": (
                "exact identifiers with concordant age and sex confirm overlap"
                if comparable else "no corresponding M10 brain-bank arm"
            ),
        })
    result = pd.DataFrame(rows)
    result.to_csv(REPO / "results" / "m9-m10-donor-overlap-audit.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
