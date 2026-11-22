#!/usr/bin/env python3
"""Append HRCA propeller results to the mature propeller outputs.

Reads results/followup-mature-propeller-hrca.csv (output of the unchanged
scripts/run_reviewer_phase1_propeller.R), appends all 88 rows to
results/followup-mature-propeller.csv, appends the 9 amacrine (AC) key rows
to results/followup-mature-propeller-key-rows.csv, and replaces the HRCA
placeholder line in reports/followup-mature-composition.md with the
per-suite table plus the diagnostic note. Verifies BH q values by
independent recomputation with statsmodels.

Notes on the HRCA key rows: the CLR effect is the batch-B mature effect,
which is per decade (batch B used age per decade), and the CLR q is the
batch-B family-A q within the suite, not a 431-family q. The report table
states both.
"""

import sys
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
REPORT = REPO / "reports" / "followup-mature-composition.md"

HRCA_CSV = RESULTS / "followup-mature-propeller-hrca.csv"
FULL_CSV = RESULTS / "followup-mature-propeller.csv"
KEY_CSV = RESULTS / "followup-mature-propeller-key-rows.csv"
BATCHB_CSV = RESULTS / "followup-mature-organ-transports-b-hrca.csv"

PLACEHOLDER = ("HRCA retina amacrine was not run: the HRCA aggregate files tripped "
               "an execution-risk check and need the user's explicit go-ahead before loading.")


def bh_check(df):
    worst = 0.0
    for _, g in df.groupby(["species", "organ", "cohort"]):
        exp = multipletests(g["pvalue"].to_numpy(), method="fdr_bh")[1]
        worst = max(worst, float(abs(exp - g["q_bh"].to_numpy()).max()))
    return worst


def main():
    hrca = pd.read_csv(HRCA_CSV)
    assert len(hrca) == 88 and hrca["cohort"].nunique() == 9, "unexpected HRCA output"
    full = pd.read_csv(FULL_CSV)
    assert not full["cohort"].str.startswith("HRCA_").any(), "HRCA rows already in full CSV"
    full2 = pd.concat([full, hrca], ignore_index=True)
    full2.to_csv(FULL_CSV, index=False)
    print(f"full propeller rows: {len(full)} -> {len(full2)}", flush=True)

    worst_hrca = bh_check(hrca)
    worst_all = bh_check(full2)
    print(f"BH max abs diff HRCA-only: {worst_hrca:.3e}; all rows: {worst_all:.3e}", flush=True)
    assert worst_all <= 1e-12, "BH mismatch"

    bb = pd.read_csv(BATCHB_CSV)
    frac = bb[(bb["endpoint_family"] == "fractions_clr") & (bb["cell_type"] == "AC")]
    assert len(frac) == 9, f"expected 9 AC fraction rows, got {len(frac)}"
    key = pd.read_csv(KEY_CSV)
    rows = []
    for _, f in frac.iterrows():
        cohort = str(f["cohort"])
        prow = hrca[(hrca["cohort"] == cohort) & (hrca["cell_type"] == "AC")]
        assert len(prow) == 1, f"{cohort}: expected 1 AC propeller row"
        p = prow.iloc[0]
        clr_sig = float(f["mature_q"]) < 0.05
        prop_sig = float(p["q_bh"]) < 0.05
        same = (float(p["effect_logit"]) > 0) == (float(f["mature_effect"]) > 0)
        if clr_sig and prop_sig and same:
            verdict = "AGREE"
        elif clr_sig:
            verdict = "NOT_CONFIRMED"
        elif prop_sig:
            verdict = "PROPELLER_ONLY"
        else:
            verdict = "NULL"
        rows.append({
            "label": f"HRCA {cohort.replace('HRCA_', '')} amacrine",
            "species": "Human", "organ": "retina", "cohort": cohort,
            "cell_type": "AC",
            "clr_effect": float(f["mature_effect"]), "clr_q_431": float(f["mature_q"]),
            "prop_effect_logit": float(p["effect_logit"]), "prop_se": float(p["se"]),
            "prop_p": float(p["pvalue"]), "prop_q": float(p["q_bh"]),
            "prop_df": float(p["df_total"]), "same_sign": bool(same), "verdict": verdict,
        })
    key2 = pd.concat([key, pd.DataFrame(rows)], ignore_index=True)
    key2.to_csv(KEY_CSV, index=False)
    for r in rows:
        print(f"{r['label']}: CLR {r['clr_effect']:.4f} q {r['clr_q_431']:.5g} | "
              f"prop p {r['prop_p']:.4g} q {r['prop_q']:.4g} | {r['verdict']}", flush=True)

    table_lines = ["HRCA retina amacrine is now run: 9 suites, 88 rows, all donors 30 years "
                   "and older, inputs built by scripts/prepare_mature_propeller_hrca.py. "
                   "Propeller confirms none of the 9 suites: every amacrine q is 0.43 or "
                   "higher. The three CLR-significant suites (Sn_Adult_Unenriched, "
                   "Sn_Adult_ExactAge, Sn_AllAges_Unenriched) are NOT_CONFIRMED. The other "
                   "six suites are NULL on both methods.",
                   "",
                   "For HRCA rows, the CLR effect is per decade (batch B used age per "
                   "decade) and the CLR q is the batch-B family-A q within the suite, not "
                   "a 431-family q. NULL means neither CLR nor propeller passes q below 0.05.",
                   "",
                   "| Row | CLR effect (per decade) | CLR q (family A) | Propeller effect (logit, per year) | Propeller p | Propeller q | Verdict |",
                   "| --- | ---: | ---: | ---: | ---: | ---: | --- |"]
    for r in rows:
        table_lines.append(f"| {r['label']} | {r['clr_effect']:.4f} | {r['clr_q_431']:.5g} | "
                           f"{r['prop_effect_logit']:.4f} | {r['prop_p']:.4g} | "
                           f"{r['prop_q']:.4g} | {r['verdict']} |")
    table_lines += ["",
                    "Why propeller disagrees: the propeller design has no study term, while "
                    "the batch-B CLR model adjusts for study, sex, peripheral fraction, and "
                    "log library size. A refit in Sn_Adult_Unenriched (scripts/diag_hrca_propeller_gap.py) "
                    "shows the amacrine CLR signal depends on that adjustment: unadjusted slope "
                    "plus 0.082 per decade (p 0.36); with study dummies plus sex, plus 0.218 "
                    "(p 0.0017); full batch-B model, plus 0.248 (p 0.00019). The logit refit "
                    "with propeller-like covariates reproduces the propeller estimate (plus "
                    "0.00245 per year), so the gap is the covariate set, not the transform "
                    "or the inputs."]
    new_text = "\n".join(table_lines)
    assert "\u2014" not in new_text and "\u2013" not in new_text, "forbidden dash in new text"

    text = REPORT.read_text(encoding="utf-8")
    assert PLACEHOLDER in text, "HRCA placeholder line not found in report"
    text = text.replace(PLACEHOLDER, new_text)
    old_count = ("23 cohorts ran; 6 mouse tissues were skipped for fewer than 5 mature "
                 "donors (Bladder, Large_Intestine, Mammary_Gland, Skin, Tongue, Trachea).")
    assert old_count in text, "cohort count sentence not found"
    text = text.replace(old_count, old_count + " 9 HRCA retina suites ran (88 rows).")
    old_verify = ("Mature propeller estimates are verified: independent BH recomputation "
                  "matches the R output to a maximum absolute difference of 2.665e-15.")
    assert old_verify in text, "BH verified line not found"
    text = text.replace(old_verify,
                        "Mature propeller estimates are verified: independent BH recomputation "
                        f"matches the R output to a maximum absolute difference of {worst_all:.3e} "
                        "across all 577 rows including the 88 HRCA rows.")
    REPORT.write_text(text, encoding="utf-8")
    print("report updated", flush=True)


if __name__ == "__main__":
    main()
