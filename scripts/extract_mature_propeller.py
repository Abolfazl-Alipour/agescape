#!/usr/bin/env python3
"""Extract mature propeller results for the headline and secondary rows.

Reads results/followup-mature-propeller.csv (produced by the unchanged
scripts/run_reviewer_phase1_propeller.R), verifies the within-cohort BH
q values by independent recomputation, and prints the propeller p and q
for every row the paper's fraction section claims or discusses.
"""

import sys
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"

# (species, organ, cohort, cell_type, label)
WANTED = [
    ("Human", "lung", "Lung_Integrated", "club cell", "HEADLINE lung club"),
    ("Human", "lung", "Lung_Integrated", "multiciliated epithelial cell", "HEADLINE lung multiciliated"),
    ("Human", "brain", "MSSM_DLPFC", "L6b glutamatergic neuron of the primary motor cortex", "HEADLINE MSSM L6b"),
    ("Human", "brain", "MSSM_DLPFC", "oligodendrocyte", "HEADLINE MSSM oligodendrocyte"),
    ("Human", "brain", "MSSM_DLPFC", "oligodendrocyte precursor cell", "HEADLINE MSSM OPC"),
    ("Mouse", "Limb_Muscle", "TMS", "macrophage", "SECONDARY mouse limb muscle macrophage"),
    ("Mouse", "Limb_Muscle", "TMS", "mesenchymal stem cell", "SECONDARY mouse limb muscle mesenchymal stem cell"),
    ("Human", "lung", "Lung_Integrated", "conventional dendritic cell", "SECONDARY lung conventional dendritic"),
    ("Human", "brain", "SEAAD_DLPFC", "L5 extratelencephalic projecting glutamatergic cortical neuron", "SECONDARY SEA-AD L5"),
    ("Human", "lung", "Lung_Integrated", "multiciliated columnar cell of tracheobronchial tree", "SECONDARY lung multiciliated columnar"),
    ("Human", "blood", "PBMC_GSE214546", "CD8 Naive", "PBMC Thomson CD8 naive"),
]


def main():
    prop = pd.read_csv(RESULTS / "followup-mature-propeller.csv")
    print(f"propeller rows: {len(prop)}, cohorts: {prop['cohort'].nunique()}", flush=True)

    # Independent BH verification within each (species, organ, cohort).
    worst = 0.0
    for _, group in prop.groupby(["species", "organ", "cohort"]):
        expected = multipletests(group["pvalue"].to_numpy(), method="fdr_bh")[1]
        diff = abs(expected - group["q_bh"].to_numpy()).max()
        worst = max(worst, float(diff))
    print(f"BH verification max abs diff: {worst:.3e}", flush=True)
    assert worst <= 1e-12, "BH mismatch"

    # CLR effects for sign agreement.
    mouse = pd.read_csv(RESULTS / "followup-mature-composition-mouse.csv")
    human = pd.read_csv(RESULTS / "followup-mature-composition-human.csv")
    clr = pd.concat([mouse, human], ignore_index=True)
    clr_index = {(r.species, r.organ, r.cohort, r.cell_type): r
                 for r in clr.itertuples()}
    pbmc = pd.read_csv(RESULTS / "followup-mature-pbmc-clr-models.csv")
    pbmc_index = {(r.cohort, r.cell_type): r for r in pbmc.itertuples()
                  if r.model == "primary" and r.endpoint == "CLR_composition"}

    out_rows = []
    for species, organ, cohort, cell_type, label in WANTED:
        match = prop[(prop["species"] == species) & (prop["organ"] == organ)
                     & (prop["cohort"] == cohort) & (prop["cell_type"] == cell_type)]
        if len(match) != 1:
            print(f"MISSING/AMBIGUOUS {label}: {len(match)} rows", flush=True)
            continue
        row = match.iloc[0]
        if species == "Human" and organ == "blood":
            crow = pbmc_index.get((cohort.replace("PBMC_", ""), cell_type))
            clr_effect = float(crow.effect) if crow is not None else float("nan")
            clr_q = float(crow.q_bh) if crow is not None else float("nan")
        else:
            crow = clr_index.get((species, organ, cohort, cell_type))
            clr_effect = float(crow.effect) if crow is not None else float("nan")
            clr_q = float(crow.q_431_family) if crow is not None else float("nan")
        same_sign = (row["effect_logit"] > 0) == (clr_effect > 0)
        verdict = "AGREE" if (row["q_bh"] < 0.05 and same_sign) else (
            "SIGN_CONFLICT" if not same_sign else "NOT_CONFIRMED")
        out_rows.append({
            "label": label, "species": species, "organ": organ, "cohort": cohort,
            "cell_type": cell_type, "clr_effect": clr_effect, "clr_q_431": clr_q,
            "prop_effect_logit": float(row["effect_logit"]), "prop_se": float(row["se"]),
            "prop_p": float(row["pvalue"]), "prop_q": float(row["q_bh"]),
            "prop_df": float(row["df_total"]), "same_sign": bool(same_sign),
            "verdict": verdict,
        })
        print(f"{label}: CLR effect {clr_effect:.4f} q431 {clr_q:.4g} | "
              f"prop effect {row['effect_logit']:.4f} p {row['pvalue']:.4g} "
              f"q {row['q_bh']:.4g} df {row['df_total']:.2f} | {verdict}", flush=True)

    out = pd.DataFrame(out_rows)
    out.to_csv(RESULTS / "followup-mature-propeller-key-rows.csv", index=False)
    print(f"wrote {len(out)} key rows to results/followup-mature-propeller-key-rows.csv", flush=True)


if __name__ == "__main__":
    main()
