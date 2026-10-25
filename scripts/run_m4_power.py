#!/usr/bin/env python3
"""
run_m4_power.py — Deterministic analytical power calculations for human aging cohorts.
Regenerates results/m4-human-power.csv bit-identically.
"""

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def main():
    donor_counts = [4, 6, 8, 10, 12, 16, 20, 30, 50]
    age_sds = [8.0, 12.0, 16.0, 20.0]
    dispersions = [0.1, 0.2, 0.4]

    power_records = []
    z_alpha = stats.norm.ppf(1 - 0.1 / 2)  # FDR 0.1 nominal two-sided (~1.64485)
    z_beta = stats.norm.ppf(0.80)         # 80% power (~0.84162)

    for n in donor_counts:
        for asd in age_sds:
            for phi in dispersions:
                se_slope = np.sqrt(phi / (n * (asd**2)))
                min_detectable_slope = (z_alpha + z_beta) * se_slope
                min_detectable_span_fc = min_detectable_slope * 30.0
                
                eligible = (n >= 6) and (asd >= 12.0)
                status = "Adequately Powered" if eligible else "Low-Powered (<6 donors or SD<12y)"
                
                power_records.append({
                    "n_donors": n,
                    "donor_age_sd_years": asd,
                    "biological_dispersion": phi,
                    "se_slope_per_year": round(se_slope, 5),
                    "min_detectable_slope_per_year": round(min_detectable_slope, 5),
                    "min_detectable_30y_log2fc": round(min_detectable_span_fc, 3),
                    "power_status": status,
                    "eligible_cohort": eligible
                })

    power_df = pd.DataFrame(power_records)
    power_path = RESULTS_DIR / "m4-human-power.csv"
    power_df.to_csv(power_path, index=False)
    print(f"Saved deterministic analytical power table ({len(power_df)} configurations) to {power_path}")


if __name__ == "__main__":
    main()
