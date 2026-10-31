#!/usr/bin/env python3
"""Apply one BH family to primary M9 and M10 human cell-type transport tests."""

from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def main():
    m9 = pd.read_csv(RESULTS / "m9-regressions.csv")
    m10 = pd.read_csv(RESULTS / "m10-celltype-regressions.csv")
    family = pd.concat([
        m9[["cohort", "cell_type", "slope_per_year", "p"]].assign(source="M9", model="age_plus_sex"),
        m10[["cell_type", "slope_covar_adj", "p_covar_adj"]]
        .rename(columns={"slope_covar_adj": "slope_per_year", "p_covar_adj": "p"})
        .assign(cohort="DLPFC_Lifespan_Atlas", source="M10", model="fully_adjusted"),
    ], ignore_index=True)
    family.insert(0, "test_id", [f"transport_{i + 1:03d}" for i in range(len(family))])
    family["program_tier"] = "P1_Clean"
    family["multiplicity_family"] = "primary_M9_M10_cell_type_transport"
    family["q_global_bh"] = multipletests(family["p"], method="fdr_bh")[1]
    family["survives_global_bh_05"] = family["q_global_bh"] < 0.05
    family.to_csv(RESULTS / "human-transport-global-multiplicity.csv", index=False)
    family[family["survives_global_bh_05"]].to_csv(RESULTS / "human-transport-global-survivors.csv", index=False)

    m9 = m9.drop(columns=["q_global_bh"], errors="ignore")
    m9 = m9.merge(family.loc[family["source"] == "M9", ["cohort", "cell_type", "q_global_bh"]],
                  on=["cohort", "cell_type"], how="left", validate="one_to_one")
    m9.to_csv(RESULTS / "m9-regressions.csv", index=False)
    m10 = m10.drop(columns=["q_global_bh", "q_global_bh_x", "q_global_bh_y"], errors="ignore")
    m10 = m10.merge(family.loc[family["source"] == "M10", ["cell_type", "q_global_bh"]],
                    on="cell_type", how="left", validate="one_to_one")
    m10.to_csv(RESULTS / "m10-celltype-regressions.csv", index=False)
    print(f"Primary family: {len(family)} tests; survivors: {int(family['survives_global_bh_05'].sum())}")
    print(family[family["survives_global_bh_05"]].to_string(index=False))


if __name__ == "__main__":
    main()
