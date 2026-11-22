#!/usr/bin/env python3
"""Rerun PBMC targeted genes, P1m, CLR, and Hellinger at 30-39 versus 60-plus."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import resource
import sys
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_pbmc import (
    MIN_CELLS,
    MIN_DONORS,
    MIN_REFERENCE_CELLS,
    MIN_UMI,
    NULL_ITERATIONS,
    NULL_SEED,
    TARGET_COMPOSITION,
    TARGET_GENES,
    base_row,
    build_cache,
    donor_table,
    load_thomson,
    load_zhu,
    pseudobulk_outcomes,
)
from src.followup_pbmc import (
    add_bh,
    clr_transform,
    common_cell_type,
    design_matrix,
    fast_effect,
    fit_donor_model,
    fit_vectorized,
    hellinger_distances,
    random_signature_scores,
    technical_gene_mask,
)

RESULTS = REPO / "results"
PREFIX = "followup-mature-pbmc"
DEFAULT_DATA = REPO.parent / "agescape-revision-data" / "pbmc"
DESIGN = "mature_30_39_vs_60_plus"
INTERNAL_DESIGN = "strict_group"
RUNTIME = RESULTS / f"{PREFIX}-runtime.csv"
COHORTS = ["GSE213516", "GSE214546"]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def append_runtime(cohort: str, start_utc: str, started: float, status: str, details: str = "") -> None:
    row = pd.DataFrame([{
        "pipeline": f"mature_pbmc_{cohort}", "start_utc": start_utc, "end_utc": utc_now(),
        "duration_seconds": time.monotonic() - started, "model_used": "Gemini 4 Neon XHigh",
        "status": status, "max_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        "details": details,
    }])
    if RUNTIME.exists():
        previous = pd.read_csv(RUNTIME)
        previous = previous[previous["pipeline"] != f"mature_pbmc_{cohort}"]
        row = pd.concat([previous, row], ignore_index=True)
    row.to_csv(RUNTIME, index=False)


def mature_row(cohort: str, endpoint: str, cell_type: str, model: str) -> dict:
    return {**base_row(cohort, endpoint, cell_type, DESIGN, model)}


def select_donors(cache: dict) -> pd.DataFrame:
    donors = donor_table(cache).copy()
    baseline = donors["age"].between(30, 39, inclusive="both")
    old = donors["age"] >= 60
    donors = donors[baseline | old].copy()
    donors["old"] = old[baseline | old].astype(float)
    donors["mature_group"] = np.where(donors["old"] == 0, "MATURE_30_39", "OLD_60_PLUS")
    if donors.empty or donors["age"].min() < 30:
        raise ValueError("Invalid mature donor selection")
    if ((donors["old"] == 0) & ~donors["age"].between(30, 39)).any():
        raise ValueError("Baseline donor outside 30-39")
    return donors


def unit_frame(cache: dict, cell_type: str, cell_floor: str = "n_cells") -> pd.DataFrame:
    mask = (cache["cell_type"] == cell_type) & (cache[cell_floor] >= MIN_CELLS)
    columns = [
        "donor", "cell_type", "age", "sex", "cmv", "batch", "pool", "n_cells",
        "mean_log_umi", "sum_umi", "qc_n", "qc_mean_log_umi", "qc_sum_umi",
    ]
    frame = pd.DataFrame({key: cache[key][mask] for key in columns})
    donors = select_donors(cache)
    return frame.merge(donors[["donor", "old", "mature_group"]], on="donor", how="inner", validate="many_to_one")


def fit(frame: pd.DataFrame, outcome: str, variant: str, endpoint: str) -> dict:
    result = fit_donor_model(frame, outcome, INTERNAL_DESIGN, variant, endpoint)
    if "n_young" in result:
        result["n_mature_30_39"] = result["n_young"]
    if result.get("status") == "fit":
        n_params = len(result["covariates"].split(";"))
        rdf = int(result["n_donors"]) - n_params
        result["residual_df"] = rdf
        if rdf >= 1 and np.isfinite(result.get("se", np.nan)):
            result["mde80"] = float((stats.t.ppf(0.975, rdf) + stats.t.ppf(0.80, rdf)) * result["se"])
        else:
            result["mde80"] = np.nan
    return result


def verify_bh(frame: pd.DataFrame, groups: list[str], p_column: str = "pvalue", q_column: str = "q_bh") -> None:
    for _, indices in frame.groupby(groups, dropna=False).groups.items():
        valid = frame.loc[indices, p_column].dropna().index
        if not len(valid):
            continue
        expected = multipletests(frame.loc[valid, p_column], method="fdr_bh")[1]
        if not np.allclose(frame.loc[valid, q_column], expected, rtol=1e-12, atol=1e-15):
            raise ValueError(f"BH verification failed for {groups}")


def donor_counts(cache: dict, cohort: str) -> dict:
    donors = select_donors(cache)
    n_baseline = int((donors["old"] == 0).sum())
    n_old = int((donors["old"] == 1).sum())
    testable = n_baseline >= 3 and n_old >= 3 and len(donors) >= MIN_DONORS
    return {
        "cohort": cohort, "design": DESIGN, "n_total_selected": len(donors),
        "n_mature_30_39": n_baseline, "n_old_60_plus": n_old,
        "baseline_ages": ";".join(str(int(age)) for age in sorted(donors.loc[donors["old"] == 0, "age"])),
        "old_ages": ";".join(str(int(age)) for age in sorted(donors.loc[donors["old"] == 1, "age"])),
        "status": "TESTABLE_AT_COHORT_LEVEL" if testable else "UNTESTABLE",
        "reason": "cell_type_and_covariate_floors_still_apply" if testable else "fewer_than_three_donors_in_one_group_or_fewer_than_six_total",
    }


def eligibility(cache: dict, cohort: str) -> pd.DataFrame:
    rows = []
    for cell_type in sorted(set(cache["cell_type"])):
        frame = unit_frame(cache, cell_type)
        n_baseline = int((frame["old"] == 0).sum())
        n_old = int((frame["old"] == 1).sum())
        eligible = len(frame) >= MIN_DONORS and min(n_baseline, n_old) >= 3
        rows.append({
            **mature_row(cohort, "eligibility", cell_type, "primary"),
            "n_donors": len(frame), "n_cells": int(frame["n_cells"].sum()),
            "n_mature_30_39": n_baseline, "n_old": n_old, "eligible": eligible,
            "reason": "pass" if eligible else "donor_or_group_floor",
        })
    return pd.DataFrame(rows)


def load_p1m() -> dict[str, list[str]]:
    frame = pd.read_csv(RESULTS / "followup-mature-p1m-genes.csv")
    frame = frame[frame["human_ortholog"].notna()]
    if len(frame) != 1724 or frame["human_ortholog"].nunique() != 1724:
        raise ValueError("Unexpected P1m human-mapped gene table")
    return {
        "up": frame.loc[frame["direction"] == "UP", "human_ortholog"].astype(str).tolist(),
        "down": frame.loc[frame["direction"] == "DOWN", "human_ortholog"].astype(str).tolist(),
    }


def run_p1(cohort: str, cache: dict, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    program = load_p1m()
    symbols = cache["symbol"].astype(str)
    unique = ~np.asarray(pd.Index(symbols).duplicated(keep=False))
    up = np.flatnonzero(np.isin(symbols, program["up"]) & unique)
    down = np.flatnonzero(np.isin(symbols, program["down"]) & unique)
    if not len(up) or not len(down):
        raise ValueError(f"P1m lacks unique detected genes for {cohort}")
    model_rows, donor_rows = [], []
    rng = np.random.default_rng(NULL_SEED)
    for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"].unique()):
        mask = (cache["cell_type"] == cell_type) & (cache["n_cells"] >= MIN_CELLS)
        indices = np.flatnonzero(mask)
        frame = unit_frame(cache, cell_type)
        lookup = pd.Series(indices, index=cache["donor"][mask])
        order = lookup.loc[frame["donor"]].to_numpy()
        matrix = (cache["log_sums"][order] / cache["n_cells"][order, None]).astype(float).T
        frame["score"] = matrix[up].mean(axis=0) - matrix[down].mean(axis=0)
        frame["cohort"] = cohort
        frame["design"] = DESIGN
        frame["common_cell_type"] = common_cell_type(cell_type)
        donor_rows.append(frame)
        variants = ["primary", "depth"]
        if cohort == "GSE214546":
            variants += ["cmv", "batch", "pool"]
        for variant in variants:
            fitted = fit(frame, "score", variant, "P1_transport")
            model_rows.append({
                **mature_row(cohort, "P1m_transport", cell_type, variant),
                "p1m_up_detected": len(up), "p1m_down_detected": len(down),
                "p1m_coverage": (len(up) + len(down)) / (len(program["up"]) + len(program["down"])),
                "p1m_source": "results/followup-mature-p1m-genes.csv", **fitted,
            })
        target = next(row for row in reversed(model_rows) if row["cell_type"] == cell_type and row["model"] == "primary")
        if target.get("status") == "fit" and np.isfinite(target.get("effect", np.nan)):
            observed = float(target["effect"])
            all_p1m = np.flatnonzero(np.isin(symbols, program["up"] + program["down"]))
            random_scores = random_signature_scores(matrix, up, down, NULL_ITERATIONS, rng, exclude=all_p1m)
            random_effects = np.array([fast_effect(score, frame, INTERNAL_DESIGN) for score in random_scores])
            null_mean = float(random_effects.mean())
            empirical = (1 + int(np.sum(np.abs(random_effects - null_mean) >= abs(observed - null_mean)))) / (NULL_ITERATIONS + 1)
            target.update({
                "random_iterations": NULL_ITERATIONS, "random_seed": NULL_SEED,
                "random_mean_effect": null_mean, "random_sd": float(random_effects.std(ddof=1)),
                "random_q025": float(np.quantile(random_effects, 0.025)),
                "random_q975": float(np.quantile(random_effects, 0.975)),
                "random_pvalue": float(empirical),
            })
            print(f"{cohort} P1m {cell_type}: n={len(frame)}, null p={empirical:.4g}", flush=True)
    models = add_bh(pd.DataFrame(model_rows), ["cohort", "endpoint", "design", "model"])
    models["random_q_bh"] = np.nan
    primary = models["model"].eq("primary") & models["status"].eq("fit") & models["random_pvalue"].notna()
    if primary.any():
        models.loc[primary, "random_q_bh"] = multipletests(models.loc[primary, "random_pvalue"], method="fdr_bh")[1]
    verify_bh(models, ["cohort", "endpoint", "design", "model"])
    random_valid = models.loc[primary, "random_pvalue"].dropna()
    if len(random_valid):
        expected_random = multipletests(random_valid, method="fdr_bh")[1]
        if not np.allclose(models.loc[random_valid.index, "random_q_bh"], expected_random, rtol=1e-12, atol=1e-15):
            raise ValueError("P1m matched-null BH verification failed")
    donors = pd.concat(donor_rows, ignore_index=True) if donor_rows else pd.DataFrame()
    return donors, models


def run_targeted(cohort: str, cache: dict, eligible: pd.DataFrame) -> pd.DataFrame:
    symbols = cache["symbol"].astype(str)
    unique = ~np.asarray(pd.Index(symbols).duplicated(keep=False))
    rows = []
    cell_types = sorted({label for label in cache["cell_type"] if common_cell_type(label) == "CD8 naive T"})
    for cell_type in cell_types:
        mask = (cache["cell_type"] == cell_type) & (cache["n_cells"] >= MIN_CELLS)
        indices = np.flatnonzero(mask)
        frame = unit_frame(cache, cell_type)
        lookup = pd.Series(indices, index=cache["donor"][mask])
        order = lookup.loc[frame["donor"]].to_numpy()
        counts = cache["counts"][order].astype(float)
        outcomes = pseudobulk_outcomes(counts)
        is_eligible = bool(eligible.query("cell_type == @cell_type")["eligible"].iloc[0])
        variants = ["primary", "depth"]
        if cohort == "GSE214546":
            variants += ["cmv", "batch", "pool"]
        for variant in variants:
            for gene in TARGET_GENES:
                base = {
                    **mature_row(cohort, "pseudobulk_targeted", cell_type, variant),
                    "gene": gene, "expected_direction": "increase",
                }
                if not is_eligible:
                    rows.append({**base, "status": "untestable", "reason": "donor_or_cell_floor", "n_donors": len(frame)})
                    continue
                selected = np.flatnonzero((symbols == gene) & unique)
                if not len(selected):
                    rows.append({**base, "status": "untestable", "reason": "gene_missing_or_ambiguous_symbol"})
                    continue
                fitted = fit_vectorized(outcomes[selected], frame, INTERNAL_DESIGN, variant)
                if fitted is None:
                    rows.append({
                        **base, "feature_id": cache["feature_id"][selected[0]],
                        "status": "untestable", "reason": "rank_or_group_floor", "n_donors": len(frame),
                    })
                    continue
                X_t, _ = design_matrix(frame, INTERNAL_DESIGN, variant, "pseudobulk")
                rdf = int(len(frame) - X_t.shape[1])
                se_val = float(fitted["se"][0])
                mde80_val = float((stats.t.ppf(0.975, rdf) + stats.t.ppf(0.80, rdf)) * se_val) if rdf >= 1 else np.nan
                rows.append({
                    **base, "feature_id": cache["feature_id"][selected[0]], "status": "fit",
                    "n_donors": len(frame), "n_cells": int(frame["n_cells"].sum()),
                    "n_mature_30_39": int((frame["old"] == 0).sum()), "n_old": int((frame["old"] == 1).sum()),
                    "residual_df": rdf, "mde80": mde80_val,
                    **{key: float(fitted[key][0]) for key in ["effect", "se", "ci95_low", "ci95_high", "pvalue"]},
                })
    result = add_bh(pd.DataFrame(rows), ["cohort", "endpoint", "design", "model"])
    verify_bh(result, ["cohort", "endpoint", "design", "model"])
    return result


def composition_counts(cache: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    donors = select_donors(cache)
    labels = sorted(set(cache["cell_type"]))
    units = pd.DataFrame({"donor": cache["donor"], "cell_type": cache["cell_type"], "n_cells": cache["n_cells"]})
    wide = units.pivot(index="donor", columns="cell_type", values="n_cells").fillna(0)
    wide = wide.reindex(index=donors["donor"], columns=labels, fill_value=0).astype(int)
    meta = donors.set_index("donor", drop=False)
    meta["total_cells"] = wide.sum(axis=1)
    return wide, meta


def run_composition(cohort: str, cache: dict, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    exploratory_rows, targeted_rows = [], []
    wide, meta = composition_counts(cache)
    clr = clr_transform(wide)
    donor_long = wide.stack().rename("n_cells").reset_index()
    donor_long["cohort"] = cohort
    donor_long["design"] = DESIGN
    donor_long["total_cells"] = donor_long["donor"].map(meta["total_cells"])
    for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"].unique()):
        frame = meta.copy().reset_index(drop=True)
        frame["n_cells"] = wide[cell_type].to_numpy()
        frame["clr"] = clr[cell_type].to_numpy()
        frame["fraction"] = frame["n_cells"] / frame["total_cells"]
        variants = ["primary", "log_total_cells"]
        if cohort == "GSE214546":
            variants += ["cmv", "batch", "pool"]
        for variant in variants:
            fitted = fit(frame, "clr", variant, "CLR_composition")
            exploratory_rows.append({**mature_row(cohort, "CLR_composition", cell_type, variant), **fitted})
        fitted = fit(frame, "fraction", "primary", "CLR_composition")
        exploratory_rows.append({**mature_row(cohort, "raw_fraction_sensitivity", cell_type, "primary"), **fitted})
    for common, direction in TARGET_COMPOSITION.items():
        matching = [label for label in wide.columns if common_cell_type(label) == common]
        base = mature_row(cohort, "CLR_targeted", common, "primary")
        base["expected_direction"] = direction
        if not matching:
            targeted_rows.append({
                **base, "cell_type": common, "common_cell_type": common,
                "status": "untestable", "reason": "released_label_absent_or_not_harmonizable",
            })
            continue
        label = matching[0]
        if not bool(eligible.query("cell_type == @label")["eligible"].iloc[0]):
            targeted_rows.append({**base, "cell_type": label, "status": "untestable", "reason": "donor_or_cell_floor"})
            continue
        frame = meta.copy().reset_index(drop=True)
        frame["n_cells"] = wide[label].to_numpy()
        frame["clr"] = clr[label].to_numpy()
        fitted = fit(frame, "clr", "primary", "CLR_composition")
        targeted_rows.append({**base, "cell_type": label, **fitted})
    exploratory = add_bh(pd.DataFrame(exploratory_rows), ["cohort", "endpoint", "design", "model"])
    targeted = add_bh(pd.DataFrame(targeted_rows), ["cohort", "endpoint", "design", "model"])
    verify_bh(exploratory, ["cohort", "endpoint", "design", "model"])
    verify_bh(targeted, ["cohort", "endpoint", "design", "model"])
    return exploratory, targeted, donor_long


def run_hellinger(cohort: str, cache: dict, data_root: Path, manifest: dict,
                  annotations: pd.DataFrame | None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    technical = technical_gene_mask(cache["symbol"].astype(str))
    keep_genes = ~technical
    eligibility_rows, donor_rows, model_rows = [], [], []
    references = {}
    reference_cells = {}
    for cell_type in sorted(set(cache["cell_type"])):
        baseline_indices = np.flatnonzero(
            (cache["cell_type"] == cell_type) & (cache["age"] >= 30) & (cache["age"] <= 39) & (cache["qc_n"] >= MIN_CELLS)
        )
        total_reference = cache["qc_counts"][baseline_indices].sum(axis=0).astype(np.int64)
        total_reference[technical] = 0
        total_cells = int(cache["qc_n"][baseline_indices].sum())
        references[cell_type] = {}
        reference_cells[cell_type] = {}
        frame = unit_frame(cache, cell_type, cell_floor="qc_n")
        minimum_reference = []
        for row in frame.itertuples():
            index = np.flatnonzero((cache["cell_type"] == cell_type) & (cache["donor"] == row.donor))[0]
            is_baseline = 30 <= cache["age"][index] <= 39
            own = cache["qc_counts"][index].astype(np.int64) if is_baseline else np.zeros(len(technical), dtype=np.int64)
            own[technical] = 0
            reference = total_reference - own
            own_cells = int(cache["qc_n"][index]) if is_baseline else 0
            reference_n = total_cells - own_cells
            if len(baseline_indices) < 2 or reference_n < MIN_REFERENCE_CELLS or reference.sum() <= 0:
                continue
            references[cell_type][row.donor] = reference / reference.sum()
            reference_cells[cell_type][row.donor] = reference_n
            minimum_reference.append(reference_n)
        eligible = len(references[cell_type]) >= MIN_DONORS
        eligibility_rows.append({
            **mature_row(cohort, "Hellinger_eligibility", cell_type, "primary"),
            "n_donors_qc20": len(frame), "n_mature_reference_donors": len(baseline_indices),
            "n_scored_donors": len(references[cell_type]),
            "min_cross_fitted_reference_cells": min(minimum_reference) if minimum_reference else 0,
            "eligible": eligible, "reason": "pass" if eligible else "mature_reference_or_scored_donor_floor",
        })
    needed = {(cell_type, donor) for cell_type, mapping in references.items() for donor in mapping}
    loader = load_zhu if cohort == "GSE213516" else load_thomson
    donor_lookup = {record["donor"]: record for record in manifest["donors"]}
    for donor in sorted({donor for _, donor in needed}):
        donor_meta = donor_lookup[donor]
        matrix, obs, features, _ = loader(data_root, donor_meta, annotations)
        libraries = np.asarray(matrix.sum(axis=1)).ravel()
        nontechnical = np.asarray(matrix[:, keep_genes].sum(axis=1)).ravel()
        for cell_type in sorted({cell_type for cell_type, selected_donor in needed if selected_donor == donor}):
            mask = (obs["cell_type"].to_numpy() == cell_type) & (libraries >= MIN_UMI) & (nontechnical > 0)
            if not mask.any():
                continue
            distances = hellinger_distances(matrix[mask], references[cell_type][donor], nontechnical[mask])
            donor_rows.append({
                "cohort": cohort, "cell_type": cell_type, "common_cell_type": common_cell_type(cell_type),
                "design": DESIGN, "donor": donor, "age": float(donor_meta["age"]), "sex": donor_meta["sex"],
                "cmv": donor_meta.get("cmv", ""), "n_cells": int(mask.sum()),
                "mean_log_umi": float(np.log(libraries[mask]).mean()),
                "qc_mean_log_umi": float(np.log(libraries[mask]).mean()),
                "mean_umi": float(libraries[mask].mean()), "hellinger": float(distances.mean()),
                "cross_fitted_reference_cells": int(reference_cells[cell_type][donor]),
                "reference_age_min": 30, "reference_age_max": 39,
            })
        del matrix, obs, features
        gc.collect()
    donors = pd.DataFrame(donor_rows)
    eligibility_frame = pd.DataFrame(eligibility_rows)
    if not donors.empty:
        cache_donors = donor_table(cache)
        for cell_type, part in donors.groupby("cell_type"):
            eligibility_row = eligibility_frame[eligibility_frame["cell_type"] == cell_type]
            if eligibility_row.empty or not bool(eligibility_row["eligible"].iloc[0]):
                continue
            part = part.merge(cache_donors[["donor", "batch", "pool"]], on="donor", how="left", validate="many_to_one")
            part = part.merge(select_donors(cache)[["donor", "old"]], on="donor", how="left", validate="many_to_one")
            variants = ["primary", "unadjusted"]
            if cohort == "GSE214546":
                variants += ["cmv_depth", "batch_depth", "pool_depth"]
            for variant in variants:
                fitted = fit(part, "hellinger", variant, "Hellinger")
                model_rows.append({**mature_row(cohort, "Hellinger", cell_type, variant), **fitted})
    models = add_bh(pd.DataFrame(model_rows), ["cohort", "endpoint", "design", "model"])
    verify_bh(models, ["cohort", "endpoint", "design", "model"])
    if not donors.empty and donors["age"].min() < 30:
        raise ValueError("Under-30 Hellinger donor")
    return donors, models, eligibility_frame


def identifiability(cohort: str, cache: dict) -> pd.DataFrame:
    rows = []
    donors = select_donors(cache)
    for variant in ["primary", "depth", "log_total_cells", "cmv", "batch", "pool"]:
        if cohort != "GSE214546" and variant in {"cmv", "batch", "pool"}:
            continue
        frame = donors.copy()
        frame["n_cells"] = 1
        frame["mean_log_umi"] = frame["donor"].map(pd.DataFrame({"donor": cache["donor"], "value": cache["mean_log_umi"]}).groupby("donor")["value"].mean())
        frame["qc_mean_log_umi"] = frame["mean_log_umi"]
        frame["total_cells"] = frame["donor"].map(pd.DataFrame({"donor": cache["donor"], "value": cache["n_cells"]}).groupby("donor")["value"].sum())
        matrix, names = design_matrix(frame, INTERNAL_DESIGN, variant, "P1_transport")
        rows.append({
            "cohort": cohort, "design": DESIGN, "model": variant, "n_donors": len(frame),
            "n_parameters": matrix.shape[1], "rank": int(np.linalg.matrix_rank(matrix)),
            "residual_df": len(frame) - matrix.shape[1],
            "identifiable": np.linalg.matrix_rank(matrix) == matrix.shape[1] and len(frame) > matrix.shape[1] + 1,
            "covariates": ";".join(names), "sex_levels": ";".join(sorted(frame["sex"].unique())),
            "cmv_levels": ";".join(sorted(frame["cmv"].unique())),
            "batch_levels": ";".join(sorted(frame["batch"].unique())),
            "pool_levels": ";".join(sorted(frame["pool"].unique())),
        })
    return pd.DataFrame(rows)


def run_cohort(cohort: str, data_root: Path) -> None:
    start_utc = utc_now()
    started = time.monotonic()
    try:
        full_manifest = json.loads((REPO / "data" / "manifest" / "pbmc_geo_cohorts.json").read_text())
        metadata = dict(full_manifest[cohort])
        metadata["accession"] = cohort
        cache, features, verification = build_cache(cohort, data_root, full_manifest)
        cache["batch"] = cache["batch"].astype(str)
        cache["pool"] = cache["pool"].astype(str)
        counts = donor_counts(cache, cohort)
        pd.DataFrame([counts]).to_csv(RESULTS / f"{PREFIX}-{cohort}-donor-counts.csv", index=False)
        eligible = eligibility(cache, cohort)
        eligible.to_csv(RESULTS / f"{PREFIX}-{cohort}-eligibility.csv", index=False)
        identifiability(cohort, cache).to_csv(RESULTS / f"{PREFIX}-{cohort}-covariate-identifiability.csv", index=False)
        annotations = None
        if cohort == "GSE213516":
            annotations = pd.read_csv(data_root / "raw" / cohort / metadata["annotations"]["file"])
            annotations = annotations[annotations["Cohort"].eq("SE")].copy()
        p1_donors, p1_models = run_p1(cohort, cache, eligible)
        p1_donors.to_csv(RESULTS / f"{PREFIX}-{cohort}-p1m-donor-scores.csv", index=False)
        p1_models.to_csv(RESULTS / f"{PREFIX}-{cohort}-p1m-models.csv", index=False)
        targeted = run_targeted(cohort, cache, eligible)
        targeted.to_csv(RESULTS / f"{PREFIX}-{cohort}-pseudobulk-targeted.csv", index=False)
        composition, targeted_composition, composition_donors = run_composition(cohort, cache, eligible)
        composition.to_csv(RESULTS / f"{PREFIX}-{cohort}-clr-models.csv", index=False)
        targeted_composition.to_csv(RESULTS / f"{PREFIX}-{cohort}-clr-targeted.csv", index=False)
        composition_donors.to_csv(RESULTS / f"{PREFIX}-{cohort}-composition-donor-counts.csv", index=False)
        hellinger_donors, hellinger_models, hellinger_eligibility = run_hellinger(cohort, cache, data_root, metadata, annotations)
        hellinger_donors.to_csv(RESULTS / f"{PREFIX}-{cohort}-hellinger-donors.csv", index=False)
        hellinger_models.to_csv(RESULTS / f"{PREFIX}-{cohort}-hellinger-models.csv", index=False)
        hellinger_eligibility.to_csv(RESULTS / f"{PREFIX}-{cohort}-hellinger-eligibility.csv", index=False)
        details = (
            f"baseline={counts['n_mature_30_39']}; old={counts['n_old_60_plus']}; "
            f"eligible_types={int(eligible['eligible'].sum())}; p1m_rows={len(p1_models)}; "
            f"hellinger_rows={len(hellinger_models)}"
        )
        append_runtime(cohort, start_utc, started, "completed", details)
        print(f"Completed {cohort}: {details}", flush=True)
    except Exception as exc:
        append_runtime(cohort, start_utc, started, "FAILED", f"{type(exc).__name__}: {exc}")
        raise


def concatenate(suffix: str, output: str) -> None:
    frames = [pd.read_csv(RESULTS / f"{PREFIX}-{cohort}-{suffix}") for cohort in COHORTS]
    pd.concat(frames, ignore_index=True, sort=False).to_csv(RESULTS / f"{PREFIX}-{output}", index=False)


def aggregate() -> None:
    for suffix, output in (
        ("donor-counts.csv", "donor-counts.csv"),
        ("p1m-models.csv", "p1m-models.csv"),
        ("pseudobulk-targeted.csv", "pseudobulk-targeted.csv"),
        ("clr-targeted.csv", "clr-targeted.csv"),
        ("clr-models.csv", "clr-models.csv"),
        ("hellinger-models.csv", "hellinger-models.csv"),
        ("covariate-identifiability.csv", "covariate-identifiability.csv"),
    ):
        concatenate(suffix, output)
    eligibility_frames = []
    for cohort in COHORTS:
        eligibility_frames.append(pd.read_csv(RESULTS / f"{PREFIX}-{cohort}-eligibility.csv"))
        eligibility_frames.append(pd.read_csv(RESULTS / f"{PREFIX}-{cohort}-hellinger-eligibility.csv"))
    pd.concat(eligibility_frames, ignore_index=True, sort=False).to_csv(RESULTS / f"{PREFIX}-eligibility-all.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=[*COHORTS, "all"], default="all")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    for cohort in (COHORTS if args.cohort == "all" else [args.cohort]):
        run_cohort(cohort, args.data_root)
    aggregate()


if __name__ == "__main__":
    main()
