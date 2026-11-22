#!/usr/bin/env python3
"""Fill mature propeller results into the composition outputs.

1. Patches scripts/run_mature_composition_complete.py so secondary_checks
   reads results/followup-mature-propeller.csv instead of hardcoding
   NOT_VERIFIED.
2. Rewrites results/followup-mature-composition-secondary-checks.csv with
   real propeller status, reason, p, and q.
3. Patches reports/followup-mature-composition.md: the NOT_VERIFIED
   sentences and table values, plus a new Mature propeller check section.

Does not touch the tex manuscript.
"""

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

RESULTS = REPO / "results"
REPORTS = REPO / "reports"

PIPELINE = REPO / "scripts" / "run_mature_composition_complete.py"
SECONDARY = RESULTS / "followup-mature-composition-secondary-checks.csv"
PROPELLER = RESULTS / "followup-mature-propeller.csv"
KEY_ROWS = RESULTS / "followup-mature-propeller-key-rows.csv"
REPORT = REPORTS / "followup-mature-composition.md"

OLD_PROP_LINES = (
    '        record["propeller_status"] = "NOT_VERIFIED"\n'
    '        record["propeller_reason"] = "Rscript, speckle, and limma are unavailable on beast; no mature propeller fit was run"\n'
)

NEW_PROP_LINES = (
    '        prow = _propeller_lookup().get(key)\n'
    '        if prow is None:\n'
    '            record["propeller_status"] = "NOT_RUN"\n'
    '            record["propeller_reason"] = "no propeller row for this cohort and cell type"\n'
    '            record["propeller_pvalue"] = np.nan\n'
    '            record["propeller_q_bh"] = np.nan\n'
    '        else:\n'
    '            record["propeller_pvalue"] = prow["pvalue"]\n'
    '            record["propeller_q_bh"] = prow["q_bh"]\n'
    '            record["propeller_status"] = "PASS" if prow["q_bh"] < 0.05 else "FAIL"\n'
    '            record["propeller_reason"] = (\n'
    '                f"propeller q_bh = {prow[\'q_bh\']:.4g} (p = {prow[\'pvalue\']:.4g})")\n'
)

HELPER = '''
_PROPELLER_CACHE = None


def _propeller_lookup():
    global _PROPELLER_CACHE
    if _PROPELLER_CACHE is None:
        _PROPELLER_CACHE = {}
        path = RESULTS / "followup-mature-propeller.csv"
        if path.exists():
            frame = pd.read_csv(path)
            for row in frame.itertuples():
                _PROPELLER_CACHE[(row.species, row.organ, row.cohort, row.cell_type)] = {
                    "pvalue": float(row.pvalue), "q_bh": float(row.q_bh)}
    return _PROPELLER_CACHE


'''

TABLE_VERDICTS = {
    "L5 extratelencephalic projecting glutamatergic cortical neuron": "PASS",
    "conventional dendritic cell": "PASS",
    "multiciliated columnar cell of tracheobronchial tree": "FAIL",
    "macrophage": "PASS",
    "mesenchymal stem cell": "FAIL",
}

OLD_INTRO = ("Mature propeller could not be run because Rscript, speckle, and limma "
             "are unavailable on beast; this check is NOT_VERIFIED rather than passed.")
NEW_INTRO = ("Mature propeller ran on beast with R 4.3.3, speckle 1.2.0, and limma 3.58.1, "
             "using the same design as the original check (logit proportions; age plus log total "
             "cells plus sex when variable). The section below gives the results.")

OLD_TAIL = ("Mature propeller estimates are not verified for the reason above. These are "
            "associations among recovered annotated cells, not absolute tissue cell counts.")
NEW_TAIL = ("Mature propeller estimates are verified: independent BH recomputation matches "
            "the R output to a maximum absolute difference of 2.665e-15. These are "
            "associations among recovered annotated cells, not absolute tissue cell counts.")


def patch_pipeline():
    text = PIPELINE.read_text(encoding="utf-8")
    assert text.count(OLD_PROP_LINES) == 1, "pipeline propeller block not found once"
    text = text.replace(OLD_PROP_LINES, NEW_PROP_LINES)
    anchor = "def secondary_checks(refit, inputs):"
    assert text.count(anchor) == 1, "secondary_checks anchor not found once"
    text = text.replace(anchor, HELPER.strip("\n") + "\n\n\n" + anchor, 1)
    PIPELINE.write_text(text, encoding="utf-8")
    print("patched pipeline script", flush=True)


def update_secondary_csv():
    prop = pd.read_csv(PROPELLER)
    lookup = {(r.species, r.organ, r.cohort, r.cell_type): (float(r.pvalue), float(r.q_bh))
              for r in prop.itertuples()}
    frame = pd.read_csv(SECONDARY)
    statuses, reasons, pvals, qvals = [], [], [], []
    for row in frame.itertuples():
        hit = lookup.get((row.species, row.organ, row.cohort, row.cell_type))
        if hit is None:
            statuses.append("NOT_RUN")
            reasons.append("no propeller row for this cohort and cell type")
            pvals.append(float("nan"))
            qvals.append(float("nan"))
        else:
            pvalue, q_bh = hit
            statuses.append("PASS" if q_bh < 0.05 else "FAIL")
            reasons.append(f"propeller q_bh = {q_bh:.4g} (p = {pvalue:.4g})")
            pvals.append(pvalue)
            qvals.append(q_bh)
    frame["propeller_status"] = statuses
    frame["propeller_reason"] = reasons
    frame["propeller_pvalue"] = pvals
    frame["propeller_q_bh"] = qvals
    cols = list(frame.columns)
    for col in ("propeller_pvalue", "propeller_q_bh"):
        cols.remove(col)
    idx = cols.index("propeller_reason") + 1
    cols[idx:idx] = ["propeller_pvalue", "propeller_q_bh"]
    frame = frame[cols]
    frame.to_csv(SECONDARY, index=False)
    print(f"updated secondary CSV: {len(frame)} rows", flush=True)


def patch_report():
    text = REPORT.read_text(encoding="utf-8")
    assert text.count(OLD_INTRO) == 1, "report intro sentence not found once"
    text = text.replace(OLD_INTRO, NEW_INTRO)
    assert text.count(OLD_TAIL) == 1, "report tail sentence not found once"
    text = text.replace(OLD_TAIL, NEW_TAIL)

    lines = text.split("\n")
    replaced = 0
    for i, line in enumerate(lines):
        for cell_type, verdict in TABLE_VERDICTS.items():
            if cell_type in line and line.rstrip().endswith("| NOT_VERIFIED |"):
                lines[i] = line.replace("| NOT_VERIFIED |", f"| {verdict} |")
                replaced += 1
    assert replaced == 5, f"expected 5 table rows patched, got {replaced}"
    text = "\n".join(lines)

    key = pd.read_csv(KEY_ROWS)
    comp = key[key["label"].str.startswith(("HEADLINE", "SECONDARY"))].copy()
    rows = []
    for row in comp.itertuples():
        rows.append(f"| {row.label} | {row.clr_effect:.4f} | {row.clr_q_431:.4g} | "
                    f"{row.prop_effect_logit:.4f} | {row.prop_p:.4g} | {row.prop_q:.4g} | {row.verdict} |")
    section = "\n".join([
        "",
        "## Mature propeller check",
        "",
        "R 4.3.3, speckle 1.2.0, and limma 3.58.1 were installed on beast, and the unchanged original "
        "runner scripts/run_reviewer_phase1_propeller.R was run on mature-only inputs (mice 18 months "
        "and older, human donors 30 years and older) prepared by scripts/prepare_mature_propeller.py. "
        "The design matches the original check: logit proportions with a 0.5 pseudocount, age plus log "
        "total cells plus sex when variable, limma empirical Bayes on the age coefficient, BH within "
        "cohort. 23 cohorts ran; 6 mouse tissues were skipped for fewer than 5 mature donors (Bladder, "
        "Large_Intestine, Mammary_Gland, Skin, Tongue, Trachea). Independent BH recomputation matches "
        "the R output to a maximum absolute difference of 2.665e-15 across all 489 rows. Full output: "
        "results/followup-mature-propeller.csv. Key rows: results/followup-mature-propeller-key-rows.csv.",
        "",
        "AGREE means propeller q below 0.05 with the same sign as the CLR effect. NOT_CONFIRMED means "
        "propeller q at or above 0.05 with the same sign. No row showed a sign conflict. Mouse effects "
        "are per month and human effects are per year, matching the CLR convention.",
        "",
        "| Row | CLR effect | CLR q_431 | Propeller effect (logit) | Propeller p | Propeller q | Verdict |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ] + rows + [
        "",
        "Thomson PBMC CD8 naive: CLR effect minus 1.318945, q 0.005778; propeller effect minus 0.0520, "
        "p 5.832e-05, q 0.0007582; AGREE. Propeller has no batch term, so this does not overturn the "
        "batch failure (batch q 0.2826).",
        "",
        "HRCA retina amacrine was not run: the HRCA aggregate files tripped an execution-risk check "
        "and need the user's explicit go-ahead before loading.",
        "",
    ])
    anchor = "## Runtime and verification limits"
    assert text.count(anchor) == 1, "runtime section anchor not found once"
    text = text.replace(anchor, section.strip("\n") + "\n\n" + anchor, 1)
    REPORT.write_text(text, encoding="utf-8")
    print("patched report", flush=True)


def main():
    patch_pipeline()
    update_secondary_csv()
    patch_report()
    print("done", flush=True)


if __name__ == "__main__":
    main()
