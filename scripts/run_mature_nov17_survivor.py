#!/usr/bin/env python3
"""Rerun the November 17 survivor-bias tests with mature-only inputs."""

import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import fit_unit

RESULTS = REPO / "results"
SEED = 20261117
N_PERM = 10000
KEY = ["species", "organ", "cohort", "cell_type"]
DRIFT_KEY = ["species", "tissue", "cohort", "cell_type"]


def load_composition() -> pd.DataFrame:
    frames = []
    for species, filename, sensitivity, minimum_age in [
        ("Mouse", "followup-mature-composition-mouse.csv", "mature_only_18mo", 18.0),
        ("Human", "followup-mature-composition-human.csv", "mature_only_30y", 30.0),
    ]:
        frame = pd.read_csv(RESULTS / filename)
        frame = frame[
            (frame["species"] == species)
            & (frame["sensitivity"] == sensitivity)
            & (frame["model"] == "age_plus_sex")
            & (frame["effect_type"] == "CLR_slope_per_year")
        ].copy()
        assert frame["min_age"].min() >= minimum_age
        assert not frame.duplicated(KEY).any()
        frames.append(frame)
    return pd.concat(frames, ignore_index=True).rename(
        columns={"effect": "fraction_clr_slope_per_year"}
    )


def load_drift_donors() -> pd.DataFrame:
    mouse = pd.read_csv(RESULTS / "followup-mature-drift-mouse-donors.csv")
    mouse = mouse[(mouse["depth"] == "full") & (mouse["age"] >= 1.5)].copy()
    assert mouse["age"].min() >= 1.5
    assert (mouse["table"] == "mouse_ref_18m").all()
    mouse = mouse.rename(columns={"profile_hellinger_centroid_mean": "profile_hellinger"})

    human = pd.read_csv(RESULTS / "followup-mature-drift-human-donors.csv")
    human = human[human["age"] >= 30].copy()
    assert human["age"].min() >= 30
    assert (human["reference_window_min"] == 30).all()
    assert (human["reference_window_max"] == 39).all()

    donors = pd.concat([mouse, human], ignore_index=True)
    assert not donors.duplicated(DRIFT_KEY + ["donor"]).any()
    return donors


def fit_mature_units(donors: pd.DataFrame, model: str) -> pd.DataFrame:
    rows = []
    for (species, tissue, cohort, cell_type), part in donors.groupby(DRIFT_KEY):
        row = fit_unit(part, "profile_hellinger", species, tissue, cohort, cell_type, model)
        if row is not None:
            rows.append(row)
    result = pd.DataFrame(rows)
    assert not result.duplicated(DRIFT_KEY + ["model"]).any()
    return result


def pearson_r(x, y):
    return stats.pearsonr(x, y).statistic


def spearman_r(x, y):
    return stats.spearmanr(x, y).statistic


def perm_p(x, y, func, rng, n_perm=N_PERM):
    observed = func(x, y)
    count = 0
    for _ in range(n_perm):
        permuted = rng.permutation(y)
        if abs(func(x, permuted)) >= abs(observed):
            count += 1
    return observed, (count + 1) / (n_perm + 1)


def main() -> int:
    start_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    rng = np.random.default_rng(SEED)

    composition = load_composition()
    donors = load_drift_donors()
    primary = fit_mature_units(donors, "depth_adjusted").rename(
        columns={"tissue": "organ", "slope_per_year": "drift_slope_per_year"}
    )
    unadjusted = fit_mature_units(donors, "unadjusted").rename(
        columns={"tissue": "organ", "slope_per_year": "drift_slope_unadj_per_year"}
    )

    joined = composition.merge(
        primary[KEY + ["drift_slope_per_year", "n_donors", "n_cells"]],
        on=KEY,
        how="inner",
        suffixes=("_comp", "_drift"),
        validate="one_to_one",
    )
    joined = joined.merge(
        unadjusted[KEY + ["drift_slope_unadj_per_year"]],
        on=KEY,
        how="left",
        validate="one_to_one",
    )
    assert joined["fraction_clr_slope_per_year"].notna().all()
    assert joined["drift_slope_per_year"].notna().all()
    assert joined["drift_slope_unadj_per_year"].notna().all()
    joined["share_loss"] = -joined["fraction_clr_slope_per_year"]
    joined["abs_drift_slope"] = joined["drift_slope_per_year"].abs()

    scatter = joined[
        KEY
        + [
            "fraction_clr_slope_per_year",
            "drift_slope_per_year",
            "drift_slope_unadj_per_year",
            "share_loss",
            "abs_drift_slope",
        ]
    ].copy()
    scatter = scatter.rename(
        columns={
            "fraction_clr_slope_per_year": "fraction_slope",
            "drift_slope_per_year": "drift_slope",
            "drift_slope_unadj_per_year": "drift_slope_unadjusted",
        }
    ).sort_values(KEY)
    scatter.to_csv(RESULTS / "followup-mature-nov17-survivor-joined.csv", index=False)

    rows = []
    for scope, subset in [
        ("Mouse", joined[joined["species"] == "Mouse"]),
        ("Human", joined[joined["species"] == "Human"]),
    ]:
        n_units = len(subset)
        tests = [
            (
                "fraction_slope_vs_drift_slope",
                subset["fraction_clr_slope_per_year"],
                subset["drift_slope_per_year"],
            ),
            (
                "fraction_slope_vs_abs_drift_slope",
                subset["fraction_clr_slope_per_year"],
                subset["abs_drift_slope"],
            ),
            (
                "share_loss_vs_drift_slope",
                subset["share_loss"],
                subset["drift_slope_per_year"],
            ),
            (
                "share_loss_vs_abs_drift_slope",
                subset["share_loss"],
                subset["abs_drift_slope"],
            ),
            (
                "fraction_slope_vs_drift_slope_unadjusted",
                subset["fraction_clr_slope_per_year"],
                subset["drift_slope_unadj_per_year"],
            ),
        ]
        for name, x, y in tests:
            x = np.asarray(x, dtype=float)
            y = np.asarray(y, dtype=float)
            pearson = stats.pearsonr(x, y)
            spearman = stats.spearmanr(x, y)
            _, permutation_pearson = perm_p(x, y, pearson_r, rng)
            _, permutation_spearman = perm_p(x, y, spearman_r, rng)
            rows.append(
                {
                    "scope": scope,
                    "test": name,
                    "method": "Pearson",
                    "correlation": pearson.statistic,
                    "p_value": pearson.pvalue,
                    "p_value_permutation": permutation_pearson,
                    "n_units": n_units,
                }
            )
            rows.append(
                {
                    "scope": scope,
                    "test": name,
                    "method": "Spearman",
                    "correlation": spearman.statistic,
                    "p_value": spearman.pvalue,
                    "p_value_permutation": permutation_spearman,
                    "n_units": n_units,
                }
            )

    result = pd.DataFrame(rows)
    result.to_csv(RESULTS / "followup-mature-nov17-survivor.csv", index=False)

    for scope, subset in [
        ("Mouse", joined[joined["species"] == "Mouse"]),
        ("Human", joined[joined["species"] == "Human"]),
    ]:
        losing = subset[subset["fraction_clr_slope_per_year"] < 0]
        gaining = subset[subset["fraction_clr_slope_per_year"] >= 0]
        print(f"{scope}: N={len(subset)} losing={len(losing)} gaining={len(gaining)}")
        print(
            f"  losing mean drift={losing['drift_slope_per_year'].mean():.9f} "
            f"mean abs drift={losing['abs_drift_slope'].mean():.9f}"
        )
        print(
            f"  gaining mean drift={gaining['drift_slope_per_year'].mean():.9f} "
            f"mean abs drift={gaining['abs_drift_slope'].mean():.9f}"
        )

    print()
    print(result.to_string(index=False))
    print()
    print(
        f"joined N: {len(joined)} | mouse: {(joined['species'] == 'Mouse').sum()} "
        f"| human: {(joined['species'] == 'Human').sum()}"
    )
    end_utc = datetime.now(timezone.utc).isoformat()
    duration = time.perf_counter() - started
    print(
        f"RUNTIME start_utc={start_utc} end_utc={end_utc} "
        f"duration_seconds={duration:.6f} model=kimi-code/k3.1-pro"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
