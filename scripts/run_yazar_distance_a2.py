"""Yazar distance A2 + runtime rows for the parallel expression A2.

Run after scripts/run_mature_between_subject_yazar_a2_parallel.py finishes.
Calls the same run_distance_a2 used by the serial pipeline, writes the
Yazar distance CSV, and appends the two runtime rows (parallel expression
A2 completion, distance A2 completion). No model involved in either step.
"""
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import scripts.run_mature_between_subject as rmb


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def append_row(pipeline: str, start_utc: str, duration: float, model: str, status: str, details: str) -> None:
    row = pd.DataFrame([{
        "pipeline": pipeline, "start_utc": start_utc, "end_utc": utc_now(),
        "duration_seconds": duration, "model_used": model, "status": status,
        "max_rss_mb": "", "details": details,
    }])
    previous = pd.read_csv(rmb.RUNTIME)
    previous = previous[previous["pipeline"] != pipeline]
    pd.concat([previous, row], ignore_index=True).to_csv(rmb.RUNTIME, index=False)


def main() -> int:
    t0 = time.monotonic()
    start_utc = utc_now()
    donors = pd.read_csv(rmb.RESULTS / "followup-mature-between-subject-yazar-donors.csv")
    distance = rmb.run_distance_a2("yazar", donors)
    distance.to_csv(rmb.RESULTS / "followup-mature-between-subject-distance-variance-yazar.csv", index=False)
    print(f"distance rows={len(distance)} tested={int((distance['status'] == 'TESTED').sum())} "
          f"untestable={int((distance['status'] == 'UNTESTABLE').sum())} in {time.monotonic() - t0:.1f}s", flush=True)
    details_d = (f"donors={donors['donor_id'].nunique()}; "
                 f"distance_a2={len(distance)}; tested={int((distance['status'] == 'TESTED').sum())}; "
                 f"untestable={int((distance['status'] == 'UNTESTABLE').sum())}")
    append_row("mature_between_subject_yazar_distance_a2", start_utc,
               time.monotonic() - t0, "direct run, no model", "completed", details_d)
    # Parallel expression A2 runtime row, duration parsed from its log.
    variance = pd.read_csv(rmb.RESULTS / "followup-mature-between-subject-variance-yazar.csv")
    ineligible = pd.read_csv(rmb.RESULTS / "followup-mature-between-subject-variance-yazar-ineligible.csv")
    log_path = Path("/tmp/yazar_a2_full.log")
    log = log_path.read_text(errors="replace")
    m = re.search(r"wrote .* in (\d+)s", log)
    dur = float(m.group(1)) if m else float("nan")
    start_a2 = datetime.fromtimestamp(log_path.stat().st_mtime, tz=timezone.utc).isoformat()
    details_a = (f"donors={donors['donor_id'].nunique()}; "
                 f"baseline_donors={donors.loc[donors['bin'] == 'YOUNG', 'donor_id'].nunique()}; "
                 f"expression_a2={len(variance)}; ineligible={len(ineligible)}")
    append_row("mature_between_subject_yazar", start_a2,
               dur, "direct parallel run, no model (6 workers)", "completed", details_a)
    print("runtime rows appended", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
