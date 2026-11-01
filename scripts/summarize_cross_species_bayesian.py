#!/usr/bin/env python3
"""Combine committed summaries from the stable cross-species Bayesian models."""

import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
ARTIFACTS = RESULTS / "h3-cross-species-bayesian"
MODEL_DIRS = [
    "mouse-unadjusted-prior-0.50",
    "mouse-depth_adjusted-prior-0.50",
    "human-unadjusted_fixed_group-prior-0.50",
    "naked-mole-rat-unadjusted-prior-0.50",
    "human-global-within-unit-unadjusted-prior-0.50",
    "human-global-within-unit-depth_adjusted-prior-0.50",
    "human-global-within-unit-depth_adjusted-prior-0.25",
    "human-global-within-unit-depth_adjusted-prior-1.00",
]


def main():
    frames = []
    diagnostics = []
    for directory in MODEL_DIRS:
        path = ARTIFACTS / directory
        effect_path = path / "posterior_effects.csv"
        if effect_path.exists():
            frames.append(pd.read_csv(effect_path))
        diagnostic = json.loads((path / "diagnostics.json").read_text())
        report = json.loads((path / "check_report.json").read_text())
        summary_path = path / "summary.csv"
        max_rhat = (pd.read_csv(summary_path)["r_hat"].max() if summary_path.exists()
                    else diagnostic["convergence"]["rhat"]["max"])
        first_effect = pd.read_csv(effect_path).iloc[0] if effect_path.exists() else None
        diagnostics.append({
            "model_directory": directory,
            "species": first_effect["species"],
            "model": first_effect["model"],
            "prior_sd": first_effect["prior_sd"],
            "max_rhat": max_rhat,
            "divergences": diagnostic["convergence"]["divergences"]["count"],
            "convergence_rating": report["convergence"]["rating"],
            "pareto_k_max": diagnostic["loo"]["pareto_k"]["max"],
            "pareto_k_n_bad": diagnostic["loo"]["pareto_k"]["n_bad"],
            "calibration_rating": report["calibration"]["rating"],
            "calibration_diagnosis": report["calibration"]["diagnosis"],
        })
    global_effects = pd.read_csv(RESULTS / "h3-human-bayesian-global-posterior.csv")
    frames = [frame for frame in frames if "global_" not in frame["model"].iloc[0]]
    frames.append(global_effects)
    pd.concat(frames, ignore_index=True).to_csv(RESULTS / "h3-cross-species-bayesian-posterior.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(RESULTS / "h3-cross-species-bayesian-diagnostics.csv", index=False)


if __name__ == "__main__":
    main()
