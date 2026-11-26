#!/usr/bin/env python3
"""Incorporate the full-577 propeller correction into the manuscript.

Replaces the 26-row fly/worm propeller figures (feb9be2) with the single
BH correction across all 577 propeller p-values:
  - Results paragraph (fly/worm checks): rewritten with 68 of 577 passes,
    32 of 39 scan passes also passing propeller, 7 drops, 36 entered rows,
    and the reason propeller is a different test rather than a strict
    confirmation of the scan. Headline rows unchanged, so the abstract is
    NOT touched.
  - Methods sentence (Primary composition family paragraph): propeller now
    covers all 577 rows under one correction; points at the new report.
  - Limitations sentence (family-pooling paragraph): states the
    interpretation limit alongside the study-aware impossibility.
No .csv string is introduced (checker bans it); sentences point at reports.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TEX = REPO / "report" / "plos-one" / "final-plos-one.tex"

OLD_RESULTS = (
    "The 26 fly and worm passes were put through the same propeller and "
    "leave-one-pool-out checks as the mammalian rows, with the pool as the "
    "unit instead of the animal: 72 fly pools and 9 worm pools. An independent "
    "logit-proportion regression adjusting for total cells and sex confirms "
    "23 of 26 within each cohort and 24 of 26 under a pooled correction "
    "across the 26 rows, with sign agreement in 24 of 26. Leave-one-pool-out "
    "refits hold the sign for all 26 rows, and every refit stays significant. "
    "Four fly rows carry a flag. Hemocyte cells show nothing under the "
    "independent method (p = 0.960). Ejaculatory bulb cells are null "
    "(p = 0.301). Adult ventral nervous system cells point the other way "
    "under the independent method and hold only in male pools. Posterior "
    "midgut cells pass the independent method but hold only in female pools. "
    "The worm row, ARSC, passes every check, including a refit adjusted for "
    "sequencing platform that leaves its slope unchanged at plus 0.2275 per "
    "day. No study-aware check is possible for fly or worm because each comes "
    "from a single study."
)

NEW_RESULTS = (
    "The propeller check now covers all 577 rows. An independent "
    "logit-proportion regression adjusting for total cells and sex was run on "
    "every row, with the pool as the unit for fly and worm and the donor as "
    "the unit for mouse and human, and its p-values were corrected once "
    "across all 577 tests. Propeller passes 68 of 577 rows (5 mouse, 14 "
    "human, 44 fly, and 5 worm), with sign agreement against the scan in 478 "
    "of 577. Of the 39 scan passes, 32 also pass propeller, and 7 drop: mouse "
    "limb muscle mesenchymal stem cells, human kidney cortex collecting duct "
    "epithelial cells, human lung multiciliated columnar cells, fly hemocyte "
    "cells, fly ejaculatory bulb cells, fly adult ventral nervous system "
    "cells, and fly male accessory gland main cells. Propeller also passes 36 "
    "rows the scan does not pass (4 mouse, 5 human, 23 fly, and 4 worm), "
    "because it tests the own-share logit slope while the scan tests the slope "
    "relative to the cross-type mean and because it adjusts for total cells, "
    "so it is a different test of the same question rather than a strict "
    "confirmation of the scan. All four headline rows pass the pooled "
    "propeller correction: MSSM L6b neurons (q = 0.0144), oligodendrocytes "
    "(q = 0.000282), and oligodendrocyte precursor cells (q = 4.17e-06), and "
    "mouse limb muscle macrophages (q = 0.0486, the closest to the line). "
    "Leave-one-pool-out refits hold the sign for all 26 fly and worm scan "
    "passes, and every refit stays significant. Four fly rows carry a flag. "
    "Hemocyte cells show nothing under the independent method (p = 0.960). "
    "Ejaculatory bulb cells are null (p = 0.301). Adult ventral nervous "
    "system cells point the other way under the independent method and hold "
    "only in male pools. Posterior midgut cells pass the independent method "
    "but hold only in female pools. The worm row, ARSC, passes every check, "
    "including a refit adjusted for sequencing platform that leaves its slope "
    "unchanged at plus 0.2275 per day. No study-aware check is possible for "
    "fly or worm because each comes from a single study."
)

OLD_METHODS = (
    "Fruit fly and worm rows were checked with a propeller regression (a "
    "logit-proportion regression adjusting for total cells and sex) and with "
    "leave-one-pool-out refits, using the pool as the unit instead of the "
    "animal, as documented in reports/followup-mature-fly-worm-checks.md, "
    "which lists the per-row result tables."
)

NEW_METHODS = (
    "All 577 rows were checked with a propeller regression (a logit-proportion "
    "regression adjusting for total cells and sex) under one "
    "Benjamini-Hochberg correction across all 577 propeller p-values, as "
    "documented in reports/followup-mature-propeller-full-577.md, which lists "
    "the per-row result tables. Fruit fly and worm rows were additionally "
    "checked with leave-one-pool-out refits, using the pool as the unit "
    "instead of the animal, as documented in "
    "reports/followup-mature-fly-worm-checks.md, which lists the per-row "
    "result tables."
)

OLD_LIMITS = (
    "No study-aware check is possible for fly or worm because each comes from "
    "a single study; the worm row was checked against sequencing platform "
    "instead and the fly rows against sex."
)

NEW_LIMITS = (
    "Propeller was corrected once across all 577 tests, so it is not a strict "
    "subset confirmation of the scan: it passes 68 rows, including 36 the "
    "scan does not pass, because it tests the own-share logit slope while the "
    "scan tests the slope relative to the cross-type mean and because it "
    "adjusts for total cells. No study-aware check is possible for fly or "
    "worm because each comes from a single study; the worm row was checked "
    "against sequencing platform instead and the fly rows against sex."
)


def main():
    text = TEX.read_text(encoding="utf-8")
    for label, old, new in (("results", OLD_RESULTS, NEW_RESULTS),
                            ("methods", OLD_METHODS, NEW_METHODS),
                            ("limitations", OLD_LIMITS, NEW_LIMITS)):
        count = text.count(old)
        assert count == 1, f"{label}: old string occurs {count} times, expected 1"
        text = text.replace(old, new)
        print(f"replaced {label} passage", flush=True)
    assert ".csv" not in text, "forbidden .csv string in tex"
    assert "\u2014" not in text and "\u2013" not in text, "forbidden dash in tex"
    TEX.write_text(text, encoding="utf-8")
    print(f"old results words: {len(OLD_RESULTS.split())}, "
          f"new results words: {len(NEW_RESULTS.split())}", flush=True)
    print("tex updated", flush=True)


if __name__ == "__main__":
    main()
