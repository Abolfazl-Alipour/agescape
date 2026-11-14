from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from src.followup_muscle import fit_ols

JIN_MODELS = ["primary_group_depth_adjusted", "unadjusted_group_sensitivity"]
POST_MODELS = [
    "primary_depth_plex_clinical",
    "technical_depth_plex_sensitivity",
    "unadjusted_sensitivity",
    "secondary_primary_plus_parity",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False)
    print(f"wrote {path}: {len(frame)} rows", flush=True)


def build_jin_design(frame: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
    if len(frame) != 8 or frame["old"].value_counts().to_dict() != {0: 4, 1: 4}:
        raise ValueError("Jin requires all four young and four reproductively-aged donors")
    columns = [np.ones(8), frame["old"].to_numpy(float)]
    names = ["intercept", "reproductively_aged"]
    if model == "primary_group_depth_adjusted":
        columns.append(frame["mean_log_umi"].to_numpy(float))
        names.append("mean_log_umi")
    elif model != "unadjusted_group_sensitivity":
        raise ValueError(f"Unknown Jin model: {model}")
    design = np.column_stack(columns)
    if not np.isfinite(design).all() or np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Jin design is not finite and full rank")
    if len(frame) - design.shape[1] < 2:
        raise ValueError("Jin design lacks residual degrees of freedom")
    return design, names


def build_postmenopausal_design(frame: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
    if len(frame) < 20:
        raise ValueError("Postmenopausal models require at least 20 eligible donors")
    columns = [np.ones(len(frame)), frame["age_years"].to_numpy(float) / 10]
    names = ["intercept", "age_per_10y"]
    if model != "unadjusted_sensitivity":
        columns.append(frame["mean_log_umi"].to_numpy(float))
        names.append("mean_log_umi")
        for plex in [1, 2, 3]:
            name = f"plex_fraction_{plex}"
            columns.append(frame[name].to_numpy(float))
            names.append(name)
    if model in {"primary_depth_plex_clinical", "secondary_primary_plus_parity"}:
        for name in ["BMI", "non_ovarian_abnormal"]:
            columns.append(frame[name].to_numpy(float))
            names.append(name)
    if model == "secondary_primary_plus_parity":
        columns.append(frame["parity_live_births"].to_numpy(float))
        names.append("parity_live_births")
    elif model not in POST_MODELS:
        raise ValueError(f"Unknown postmenopausal model: {model}")
    design = np.column_stack(columns)
    if not np.isfinite(design).all():
        raise ValueError("Postmenopausal design contains missing or nonfinite values")
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Frozen postmenopausal covariates are not identifiable/full rank")
    if len(frame) - design.shape[1] < 2:
        raise ValueError("Postmenopausal design lacks residual degrees of freedom")
    return design, names


def build_design(frame: pd.DataFrame, model: str, cohort: str) -> tuple[np.ndarray, list[str]]:
    if cohort == "GSE202601":
        return build_jin_design(frame, model)
    if cohort == "GSE329417":
        return build_postmenopausal_design(frame, model)
    raise ValueError(f"Unknown cohort: {cohort}")


def fit_one(frame: pd.DataFrame, outcome: np.ndarray, model: str, cohort: str) -> tuple[dict, np.ndarray]:
    design, names = build_design(frame, model, cohort)
    fitted = fit_ols(np.asarray(outcome, dtype=float), design)
    result = {key: float(value[0]) for key, value in fitted.items()}
    result.update({
        "n_donors": len(frame),
        "n_cells": int(frame["n_cells"].sum()),
        "design_columns": ";".join(names),
    })
    if cohort == "GSE202601":
        result.update({
            "n_young": int((frame["old"] == 0).sum()),
            "n_reproductively_aged": int((frame["old"] == 1).sum()),
            "effect_definition": "reproductively_aged_minus_young",
        })
    else:
        result.update({
            "age_min": float(frame["age_years"].min()),
            "age_max": float(frame["age_years"].max()),
            "effect_definition": "within_postmenopausal_per_10_exact_years",
        })
    return result, design


def fit_gene_models(frame: pd.DataFrame, counts: np.ndarray, expression: np.ndarray,
                    model: str, cohort: str) -> tuple[dict, np.ndarray]:
    design, names = build_design(frame, model, cohort)
    fitted = fit_ols(expression, design)
    metadata = {
        "n_donors": np.full(expression.shape[1], len(frame)),
        "n_cells": np.full(expression.shape[1], int(frame["n_cells"].sum())),
        "design_columns": np.full(expression.shape[1], ";".join(names)),
    }
    fitted.update(metadata)
    return fitted, design


def status_row(cohort: str, model: str, endpoint: str, cell_type: str,
               status: str, reason: str, outcomes_examined: bool | None = None) -> dict:
    if outcomes_examined is None:
        outcomes_examined = status == "completed"
    return {
        "cohort": cohort,
        "model": model,
        "endpoint": endpoint,
        "cell_type": cell_type,
        "status": status,
        "reason": reason,
        "outcomes_examined": outcomes_examined,
    }
