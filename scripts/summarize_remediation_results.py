#!/usr/bin/env python3
"""Write derived cohort and analysis counts used in the manuscript abstract."""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def main():
    m9 = pd.read_csv(RESULTS / "m9-cohort-summary.csv")
    m10 = pd.read_csv(RESULTS / "m10-cohort-summary.csv").iloc[0]
    transport = pd.read_csv(RESULTS / "human-transport-global-multiplicity.csv")
    human_lung = pd.read_csv(RESULTS / "h2-composition-human-lung.csv")
    mouse_lung = pd.read_csv(RESULTS / "h2-composition-mouse-lung.csv")
    comparison = pd.read_csv(RESULTS / "h2-composition-cross-species-test.csv").iloc[0]
    mouse_clock = pd.read_csv(RESULTS / "h1-epigenetic-clock-mouse.csv").iloc[0]
    human_clock = pd.read_csv(RESULTS / "h1-epigenetic-clock-human.csv")
    rows = [
        ("m9_cohorts", len(m9)),
        ("m9_scored_control_donors", int(m9["n_donors"].sum())),
        ("psychad_control_donors", int(m10["n_donors"])),
        ("methylation_array_profiles", int(mouse_clock["n_samples"] + human_clock.loc[human_clock["tissue"] == "All Tissues", "n_samples"].iloc[0])),
        ("primary_human_transport_tests", len(transport)),
        ("primary_human_transport_survivors", int(transport["survives_global_bh_05"].sum())),
        ("human_lung_adult_donors", int(human_lung["n_donors"].max())),
        ("human_lung_tested_cell_types", len(human_lung)),
        ("human_lung_fdr_significant_cell_types", int((human_lung["q_bh"] < 0.05).sum())),
        ("mouse_lung_donors", int(mouse_lung["n_donors"].max())),
        ("mouse_lung_tested_cell_types", len(mouse_lung)),
        ("mouse_lung_fdr_significant_cell_types", int((mouse_lung["q_bh"] < 0.05).sum())),
        ("cross_species_homologous_pairs", int(comparison["n_pairs"])),
        ("cross_species_sign_concordant_pairs", int(comparison["n_sign_concordant"])),
        ("mouse_clock_samples", int(mouse_clock["n_samples"])),
        ("human_clock_samples", int(human_clock.loc[human_clock["tissue"] == "All Tissues", "n_samples"].iloc[0])),
        ("mouse_clock_probes_used", 724),
        ("universal_clock_probes", 816),
    ]
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(RESULTS / "analysis-cohort-summary.csv", index=False)
    hierarchy = pd.read_csv(RESULTS / "m5-program-hierarchy.csv")
    tier_rows = []
    for tier in ("P0_Raw257", "P1_Clean"):
        subset = hierarchy[hierarchy["tier"] == tier]
        tier_rows.append({
            "tier": tier,
            "n_mouse_genes": len(subset),
            "n_mouse_up": int((subset["direction"] == "UP").sum()),
            "n_mouse_down": int((subset["direction"] == "DOWN").sum()),
            "n_human_up": 120 if tier == "P0_Raw257" else 76,
            "n_human_down": 68 if tier == "P0_Raw257" else 42,
        })
    pd.DataFrame(tier_rows).to_csv(RESULTS / "signature-tier-summary.csv", index=False)


if __name__ == "__main__":
    main()
