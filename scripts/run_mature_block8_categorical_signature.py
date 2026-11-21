#!/usr/bin/env python3
"""Mature Block 8: categorical and nonlinear P1m signature models at 30-plus."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import optimize
from scipy import sparse as sp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_between_subject_gao_yazar import (  # noqa: E402
    _prune_covariate_matrix,
    assign_age_bin,
)
from scripts.run_followup_11pt_block1_residual_dispersion import mde80_exact_nct  # noqa: E402
from scripts.run_mature_p1m_transport import load_programs  # noqa: E402
from src.followup_organs import sanitize_dashes, write_clean_csv, write_clean_markdown  # noqa: E402

DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
BLOCK1_CACHE = DATA_ROOT / "block1_cache"
CACHE_ROOT = DATA_ROOT / "block8_mature_categorical_signature_cache"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
RESULT_PATH = RESULTS_DIR / "followup-mature-block8-categorical-signature.csv"
DONOR_PATH = RESULTS_DIR / "followup-mature-block8-donor-scores.csv"
REPORT_PATH = REPORTS_DIR / "followup-mature-block8-categorical.md"
RUNTIME_PATH = RESULTS_DIR / "followup-mature-block8-runtime.csv"
MODEL = "muse-spark-direct"
FROZEN_PATH = CACHE_ROOT / "frozen_spec.md"

DATASETS = ("gao", "yazar")
BINS = ("MATURE", "MID", "OLD")
MIN_AGE = 30.0
MIN_CELLS = 20
MIN_DONORS = 6
MIN_RESID_DF = {"gao": 5, "yazar": 10}
BATCH_REF = {"gao": "ENCODE v4 (Snyder)", "yazar": None}
MIN_AVAILABLE_GIB = 12.0
MDE80_SMALL_SD = 0.50
PRIMARY_FAMILY_TESTS = ("MID_vs_MATURE", "OLD_vs_MID", "trend")


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def available_memory_gib() -> float:
    fields = {}
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        key, value = line.split(":", 1)
        fields[key] = float(value.strip().split()[0])
    return fields["MemAvailable"] / 1024.0 / 1024.0


def wait_for_memory(label: str) -> None:
    while True:
        available = available_memory_gib()
        if available >= MIN_AVAILABLE_GIB:
            log(f"Memory guard passes for {label}: {available:.2f} Gi available")
            return
        log(f"Memory guard waits for {label}: {available:.2f} Gi available, need {MIN_AVAILABLE_GIB:.2f} Gi")
        time.sleep(60)


def mature_bin(age_years: float) -> str:
    if age_years < MIN_AGE:
        return "UNDER30"
    if age_years < 40:
        return "MATURE"
    if age_years < 60:
        return "MID"
    return "OLD"


def verify_program() -> dict[str, list[str]]:
    frame = pd.read_csv(RESULTS_DIR / "followup-mature-p1m-genes.csv")
    if len(frame) != 1914 or frame["direction"].value_counts().to_dict() != {"UP": 1119, "DOWN": 795}:
        raise ValueError("P1m must contain 1914 mouse genes: 1119 UP and 795 DOWN")
    program = load_programs()["P1m"]
    if len(program["up"]) != 1013 or len(program["down"]) != 711:
        raise ValueError("Expected 1013 mapped human UP and 711 mapped human DOWN P1m genes")
    if set(program["up"]) & set(program["down"]):
        raise ValueError("P1m up and down human gene lists overlap")
    return program


def frozen_spec() -> str:
    return """# Block 8: Nonlinear / categorical age-group signature models

Status: specification frozen before running any Block 8 score or model.

## Goal

Test whether the signature effect is linear in age or concentrated in one age transition. All gene-level and signature tests so far use a linear age slope, which assumes linearity. The drift result itself is a nonlinearity finding (effect concentrated in young-to-mature, dead against the mature reference). This block asks the same question of the signature: is the effect a maturation contrast (MATURE to MID), an aging contrast (MID to OLD), or a steady trend?

## Inputs

- Gao and Yazar h5ad files (same paths as Block 1), same donor tables, same technical gene mask.
- P1m gene set (1,914 mouse genes, 1,119 up / 795 down; 1,724 human-mapped, 1,013 up / 711 down) from `results/followup-mature-p1m-genes.csv` (commit d3dc14c); verify counts before running. Do NOT re-derive it and never score the original young-derived P1 in this block.
- Scoring code: reuse the composite-score computation from the transport pipeline (`scripts/run_m9_pipeline.py` / follow-up transport scripts). Per-donor per-cell-type composite scores must be computed from full-depth within-cell proportions, the same way the transport models score donors.
- Age bins: donors under 30 are excluded from all models in this block. MATURE is 30 to 39 and replaces YOUNG (20 to 39) as the reference bin; MID stays 40 to 59 and OLD stays 60 and older.

## Exact analysis specification

- Per donor per cell type (eligibility floors per shared rule 6): P1m composite score = mean up-gene score minus mean down-gene score, donor pseudobulk, full depth, technical mask applied.
- Primary model, per dataset and cell type: OLS of composite score on categorical age bin (MATURE reference) + sex_male + batch (Gao study indicators; Yazar pool fixed effects). Report MID vs MATURE and OLD vs MATURE coefficients with 95 percent CIs and p-values, plus the OLD vs MID contrast (mature-to-old change, computed from the same fit, not a separate model).
- Trend test secondary: ordinal model with bin coded 1, 2, 3 (MATURE, MID, OLD), same covariates; report slope, CI, p-value. A significant categorical pattern with a null trend test means the effect is a step, not a slope; report that distinction explicitly.
- Unadjusted sensitivity: categorical model without batch covariates.
- Spline sensitivity (labeled, secondary): cubic B-spline on continuous age with 3 degrees of freedom, F test vs the linear model, mirroring the Phase 2 nonlinearity analysis in `reports/reviewer_response_phase2.md` section 2 and `scripts/run_reviewer_phase2_human_age.py`. Read that script and reuse its spline construction.
- Family (state before correction): BH q within each dataset across cell types, separately for the MID vs MATURE contrast, the OLD vs MID contrast, and the trend test. Three stated families per dataset; do not merge them.

## Fixed implementation details

- The P1m count was verified before model fitting from `results/followup-mature-p1m-genes.csv`: 1,914 unique mouse genes, 1,119 UP and 795 DOWN, with 1,724 human-mapped genes (1,013 UP and 711 DOWN) from `scripts/run_mature_p1m_transport.py` `load_programs()`. Mapped genes absent from a dataset are omitted from the corresponding directional denominator; detected mapped genes are not filtered by expression level.
- The existing Block 1 donor-sharded cache is reused. It contains raw integer counts after duplicate-symbol collapse and the human technical mask, restricted to post-mask UMI >= 1500 cells. The large h5ad matrices are not reread. The human technical gene mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. ERCC controls are also excluded by the cache construction.
- For each cell and direction, the gene score is the mean within-cell proportion over the mapped genes present in that dataset. For each donor and cell type, up and down scores are arithmetic means of these cell scores across eligible cells, and composite = up minus down. This is an equal-cell-weight donor pseudobulk of full-depth within-cell proportions, with no log transform and no downsampling. Cells below 1500 post-mask UMIs were excluded as an eligibility filter only.
- A donor-cell-type pair requires at least 20 eligible cells. Categorical group models require at least 6 eligible donors in each of MATURE, MID, and OLD, both before and after batch pruning. Donors under 30 are excluded before any floor is checked. At least 6 donors total are required for every fit. The repository covariate pruning requires at least 2 donors per retained batch, full rank, and residual df floors of 5 for Gao and 10 for Yazar in adjusted models and 5 in unadjusted models. Units below a floor are UNTESTABLE, never nulls, and no unadjusted fallback replaces a failed primary fit.
- Age bins are mature bins derived from `age_years`: UNDER30 under 30, MATURE 30 to 39, MID 40 to 59, OLD 60 and older. UNDER30 scores may be present in the donor-score CSV but are excluded from every model. The Block 1 cache bins are still asserted unchanged as a cache-integrity check. Gao uses intercept plus 6 indicators for its 7 studies with ENCODE v4 (Snyder) as reference when present. Yazar uses intercept plus 74 indicators for its 75 pools. Both include sex_male.
- The OLD vs MID contrast is beta_OLD minus beta_MID from the categorical fit, with its standard error from that fit's covariance matrix. The ordinal trend estimate is per one-bin increment. OLD vs MATURE is reported descriptively and is not one of the three correction families. BH q-values are assigned only to evaluable primary rows in the three stated families. Unadjusted and spline sensitivities and OLD vs MATURE have no q-value because the specification defines exactly three families per dataset.
- The spline uses `dmatrix("bs(age, df=3, degree=3, include_intercept=False)", part)` exactly as in `scripts/run_reviewer_phase2_human_age.py`, on the same pruned donors and with the same sex and batch covariates as the primary. The nested linear model uses continuous age in years plus those covariates. The F test compares spline to linear and has 2 numerator degrees of freedom when both designs are full rank.
- For categorical contrasts and trend, MDE80 is the exact two-sided noncentral-t 80 percent power effect, computed as the fitted-contrast standard error times the required noncentrality. It is reported in composite-score units and divided by the donor-score SD in the fitted cell type. A small MDE80 is predefined as at most 0.50 donor-score SD per contrast or ordinal-bin step, matching the repository convention. For the spline F test, MDE80 is the exact noncentral-F noncentrality giving 80 percent power at alpha 0.05; its partial-R2 equivalent is also reported.
- No split-half technical null is required for these donor-score regressions. No random procedure is used, so no random seed is needed. No new downloads are used. Available memory is checked from `/proc/meminfo` MemAvailable before each cache-shard stage; below 20 Gi the process waits and retries every 60 seconds.

## Outputs

- `results/followup-mature-block8-categorical-signature.csv`: one row per dataset, cell type, contrast or test, model (primary, unadjusted, spline), with estimate, CI, p, q, N donors per bin, MDE80, and an age_window column recording 30-plus.
- `results/followup-mature-block8-donor-scores.csv`: per-donor per-cell-type P1m composite scores with age bin, sex, and batch, so the categorical fits are reproducible from the CSV alone. UNDER30 rows are present in this CSV but excluded from every model.
- `reports/followup-mature-block8-categorical.md`: spec, per-cell-type classification (maturation-concentrated, aging-concentrated, steady trend, null), and verdicts.
- `scripts/run_mature_block8_categorical_signature.py`.

## Pre-stated classification and verdict rules

A cell type is classified maturation-concentrated only if MID vs MATURE passes BH correction and OLD vs MID does not; aging-concentrated only if OLD vs MID passes correction and MID vs MATURE does not; steady trend only if the trend test passes correction and both steps agree in direction. Transition-concentrated classifications take precedence if more than one rule is met. Both steps agree in direction only when both estimated steps are nonzero and have the same sign. All other evaluable cell types are null. Unevaluable cell types are UNTESTABLE.

Block-level verdict: PROVES (signature effect is nonlinear / transition-concentrated) only if at least one cell type is classified maturation-concentrated or aging-concentrated after correction in the batch-adjusted primary model. DOES NOT PROVE if all passing cell types are steady-trend or if nothing passes with small MDE80. INCONCLUSIVE if nothing passes with large MDE80.

For a no-transition dataset or overall summary, any steady-trend cell type gives DOES NOT PROVE. If no cell type passes any of the three primary families, DOES NOT PROVE requires every evaluated primary MDE80 in those families to be at most 0.50 donor-score SD; otherwise the summary is INCONCLUSIVE. A significant pattern that meets no classification rule is INCONCLUSIVE. The overall summary is PROVES if either dataset is PROVES, DOES NOT PROVE only if both datasets are DOES NOT PROVE, and otherwise INCONCLUSIVE.

The classification will be compared against mature linear-slope results from `results/followup-mature-block5-gao-yazar-transport.csv` (primary arm, gene set P1m, composite BH q). For each passing cell type, the report states whether the mature linear model had already flagged it at composite BH q < 0.05. These mature transport fits use their log1p(CPM) depth-adjusted implementation while the categorical fits use full-depth proportion scores, so the comparison concerns prior flagging, not numerical identity with the new full-depth proportion score.

Exactly one deep review round follows completion of the outputs, followed by fixes, verification, commit, and push.
"""


def freeze_spec() -> None:
    verify_program()
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    text = frozen_spec()
    if FROZEN_PATH.exists() and FROZEN_PATH.read_text(encoding="utf-8") != text:
        raise RuntimeError("A different Block 8 specification is already frozen")
    if REPORT_PATH.exists():
        report = REPORT_PATH.read_text(encoding="utf-8")
        if report != text and not report.startswith(text + "\n## Results\n"):
            raise RuntimeError("The frozen Block 8 report prefix differs from the specification")
    write_clean_markdown(text, FROZEN_PATH)
    write_clean_markdown(text, REPORT_PATH)
    log("Froze the exact Block 8 specification before scoring or model fitting")


def verify_frozen_spec() -> None:
    if not FROZEN_PATH.exists() or not REPORT_PATH.exists():
        raise FileNotFoundError("Run --freeze-spec before analysis")
    if FROZEN_PATH.read_text(encoding="utf-8") != frozen_spec():
        raise RuntimeError("Frozen specification differs from analysis script")
    report = REPORT_PATH.read_text(encoding="utf-8")
    if report != frozen_spec() and not report.startswith(frozen_spec() + "\n## Results\n"):
        raise RuntimeError("Frozen report prefix was edited")


def load_info(dataset: str) -> dict[str, Any]:
    path = BLOCK1_CACHE / dataset / "CACHE_DONE.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing required Block 1 cache: {path}")
    info = json.loads(path.read_text(encoding="utf-8"))
    if not (BLOCK1_CACHE / dataset / "donors_meta.csv").exists():
        raise FileNotFoundError(f"Missing donor metadata for {dataset}")
    return info


def donor_metadata(dataset: str, info: dict[str, Any]) -> pd.DataFrame:
    donors = pd.read_csv(BLOCK1_CACHE / dataset / "donors_meta.csv")
    donors["donor_id"] = donors["donor_id"].astype(str)
    donors["donor_idx"] = np.arange(len(donors), dtype=np.int32)
    donors = donors.drop(columns=["n_cells"], errors="ignore")
    donors["study_or_pool"] = (
        donors["study"].astype(str) if dataset == "gao" else donors["pool_number"].astype(int).astype(str)
    )
    donors["batch_label"] = "study" if dataset == "gao" else "pool_number"
    if donors["donor_id"].tolist() != info["donor_ids"]:
        raise AssertionError(f"Donor order mismatch in {dataset} cache")
    recomputed_bins = donors["age_years"].map(assign_age_bin)
    if not recomputed_bins.equals(donors["bin"]):
        raise AssertionError(f"Age-bin mismatch in {dataset} cache")
    donors["mature_bin"] = donors["age_years"].map(mature_bin)
    if int((donors["age_years"] >= MIN_AGE).sum()) == 0:
        raise ValueError(f"No mature donors in {dataset} cache")
    return donors


def run_scores(dataset: str, program: dict[str, list[str]]) -> None:
    info = load_info(dataset)
    out_dir = CACHE_ROOT / dataset / "scores"
    out_dir.mkdir(parents=True, exist_ok=True)
    symbols = list(info["retained_symbols"])
    symbol_to_idx = {gene: idx for idx, gene in enumerate(symbols)}
    up_genes = [gene for gene in program["up"] if gene in symbol_to_idx]
    down_genes = [gene for gene in program["down"] if gene in symbol_to_idx]
    if not up_genes or not down_genes:
        raise ValueError(f"No mapped P1 genes in {dataset}")
    up_idx = np.array([symbol_to_idx[gene] for gene in up_genes], dtype=np.int32)
    down_idx = np.array([symbol_to_idx[gene] for gene in down_genes], dtype=np.int32)
    n_donors = len(info["donor_ids"])
    n_cell_types = len(info["cell_types"])
    n_genes = len(symbols)
    for shard in range(int(info["n_shards"])):
        output = out_dir / f"scores_shard{shard:02d}.npz"
        if output.exists():
            continue
        wait_for_memory(f"{dataset} score shard {shard}")
        files = sorted((BLOCK1_CACHE / dataset).glob(f"shard{shard:02d}_chunk*.npz"))
        up_sum = np.zeros(n_donors * n_cell_types, dtype=np.float64)
        down_sum = np.zeros_like(up_sum)
        counts = np.zeros(n_donors * n_cell_types, dtype=np.int64)
        for file in files:
            with np.load(file) as z:
                n_rows = len(z["donor"])
                matrix = sp.csr_matrix(
                    (z["data"], z["indices"], z["indptr"]), shape=(n_rows, n_genes)
                )
                donor = z["donor"].astype(np.int64)
                cell_type = z["ct"].astype(np.int64)
                umi = z["umi"].astype(np.float64)
            up_cell = np.asarray(matrix[:, up_idx].sum(axis=1)).ravel() / (len(up_idx) * umi)
            down_cell = np.asarray(matrix[:, down_idx].sum(axis=1)).ravel() / (len(down_idx) * umi)
            group = donor * n_cell_types + cell_type
            up_sum += np.bincount(group, weights=up_cell, minlength=len(up_sum))
            down_sum += np.bincount(group, weights=down_cell, minlength=len(down_sum))
            counts += np.bincount(group, minlength=len(counts))
        selected = np.flatnonzero(counts)
        np.savez_compressed(
            output,
            donor_idx=(selected // n_cell_types).astype(np.int32),
            ct_idx=(selected % n_cell_types).astype(np.int32),
            n_cells=counts[selected].astype(np.int32),
            up_score=up_sum[selected] / counts[selected],
            down_score=down_sum[selected] / counts[selected],
        )
        log(f"{dataset} shard {shard + 1}/{info['n_shards']}: aggregated {len(selected)} donor-cell-type pairs")
    audit = {
        "dataset": dataset,
        "p1_mouse_genes": 1914,
        "p1_mouse_up": 1119,
        "p1_mouse_down": 795,
        "mapped_human_up": len(program["up"]),
        "mapped_human_down": len(program["down"]),
        "present_up": up_genes,
        "present_down": down_genes,
        "absent_up": sorted(set(program["up"]) - set(up_genes)),
        "absent_down": sorted(set(program["down"]) - set(down_genes)),
        "n_shards": int(info["n_shards"]),
    }
    write_clean_markdown(json.dumps(audit, indent=2) + "\n", out_dir / "SCORES_DONE.json")


def build_donor_scores() -> pd.DataFrame:
    pieces = []
    coverage = {}
    for dataset in DATASETS:
        info = load_info(dataset)
        donors = donor_metadata(dataset, info)
        score_dir = CACHE_ROOT / dataset / "scores"
        audit_path = score_dir / "SCORES_DONE.json"
        if not audit_path.exists():
            raise FileNotFoundError(f"Missing completed score cache: {audit_path}")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        coverage[dataset] = audit
        arrays = []
        for file in sorted(score_dir.glob("scores_shard*.npz")):
            with np.load(file) as z:
                arrays.append({key: z[key].copy() for key in z.files})
        if not arrays:
            raise ValueError(f"No scores found for {dataset}")
        merged = {key: np.concatenate([piece[key] for piece in arrays]) for key in arrays[0]}
        frame = pd.DataFrame(merged)
        frame["donor_idx"] = frame["donor_idx"].astype(int)
        frame["ct_idx"] = frame["ct_idx"].astype(int)
        frame["cell_type"] = frame["ct_idx"].map(dict(enumerate(info["cell_types"])))
        frame = frame.merge(donors, on="donor_idx", how="left", validate="many_to_one")
        frame["dataset"] = dataset
        frame["age_bin"] = frame["mature_bin"]
        frame["composite_score"] = frame["up_score"] - frame["down_score"]
        frame["n_up_genes"] = len(audit["present_up"])
        frame["n_down_genes"] = len(audit["present_down"])
        frame["n_p1_genes"] = frame["n_up_genes"] + frame["n_down_genes"]
        pieces.append(frame)
    scores = pd.concat(pieces, ignore_index=True)
    columns = [
        "dataset", "donor_id", "cell_type", "age_years", "age_per_decade", "age_bin",
        "sex", "sex_male", "batch_label", "study_or_pool", "n_cells", "up_score",
        "down_score", "composite_score", "n_up_genes", "n_down_genes", "n_p1_genes",
        "age_window",
    ]
    scores["age_window"] = "30-plus; MATURE 30 to 39 baseline"
    scores = scores[columns].sort_values(["dataset", "cell_type", "donor_id"]).reset_index(drop=True)
    if scores[["dataset", "donor_id", "cell_type"]].duplicated().any():
        raise AssertionError("Duplicate donor-cell-type score rows")
    if not np.isfinite(scores[["up_score", "down_score", "composite_score"]].to_numpy()).all():
        raise AssertionError("Non-finite donor scores")
    write_clean_csv(scores, DONOR_PATH)
    write_clean_markdown(json.dumps(coverage, indent=2) + "\n", CACHE_ROOT / "score_coverage.json")
    return scores


def base_row(dataset: str, cell_type: str, test: str, model: str) -> dict[str, Any]:
    return {
        "dataset": dataset,
        "cell_type": cell_type,
        "test": test,
        "model": model,
        "status": "UNTESTABLE",
        "estimate": np.nan,
        "se": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "pvalue": np.nan,
        "q_bh": np.nan,
        "n_donors": 0,
        "n_mature": 0,
        "n_mid": 0,
        "n_old": 0,
        "residual_df": np.nan,
        "outcome_sd": np.nan,
        "mde80": np.nan,
        "mde80_sd": np.nan,
        "mde80_scale": "composite_score",
        "spline_df_num": np.nan,
        "spline_df_den": np.nan,
        "mde80_partial_r2": np.nan,
        "covariates_active": "none",
        "covariates_omitted": "none",
    }


def bin_counts(frame: pd.DataFrame) -> dict[str, int]:
    counts = frame["age_bin"].value_counts()
    return {
        "n_donors": int(len(frame)),
        "n_mature": int(counts.get("MATURE", 0)),
        "n_mid": int(counts.get("MID", 0)),
        "n_old": int(counts.get("OLD", 0)),
    }


def prune_model_frame(dataset: str, frame: pd.DataFrame, adjusted: bool) -> tuple[pd.DataFrame, list[str], list[str], str]:
    work = frame.copy()
    work["ordinal"] = work["age_bin"].map({"MATURE": 1.0, "MID": 2.0, "OLD": 3.0})
    return _prune_covariate_matrix(
        work,
        nonbatch_covs=["sex_male"],
        batch_col="study_or_pool" if adjusted else None,
        min_batch_donors=2,
        min_residual_df=MIN_RESID_DF[dataset] if adjusted else 5,
        ref_level=BATCH_REF[dataset] if adjusted else None,
        require_is_comp_collinear_check="ordinal",
        atomic_batch_block=True,
    )


def coefficient_result(fit: Any, label: str, sd: float) -> dict[str, float]:
    estimate = float(fit.params[label])
    se = float(fit.bse[label])
    df = float(fit.df_resid)
    critical = float(stats.t.ppf(0.975, df))
    mde = float(mde80_exact_nct(se, df))
    return {
        "estimate": estimate,
        "se": se,
        "ci95_low": estimate - critical * se,
        "ci95_high": estimate + critical * se,
        "pvalue": float(fit.pvalues[label]),
        "residual_df": df,
        "outcome_sd": sd,
        "mde80": mde,
        "mde80_sd": mde / sd,
    }


def old_mid_result(fit: Any, sd: float) -> dict[str, float]:
    contrast = np.zeros(len(fit.params), dtype=float)
    contrast[list(fit.params.index).index("OLD")] = 1.0
    contrast[list(fit.params.index).index("MID")] = -1.0
    estimate = float(contrast @ fit.params.to_numpy())
    variance = float(contrast @ fit.cov_params().to_numpy() @ contrast)
    se = float(np.sqrt(max(variance, 0.0)))
    df = float(fit.df_resid)
    critical = float(stats.t.ppf(0.975, df))
    statistic = estimate / se if se > 0 else np.nan
    mde = float(mde80_exact_nct(se, df))
    return {
        "estimate": estimate,
        "se": se,
        "ci95_low": estimate - critical * se,
        "ci95_high": estimate + critical * se,
        "pvalue": float(2.0 * stats.t.sf(abs(statistic), df)),
        "residual_df": df,
        "outcome_sd": sd,
        "mde80": mde,
        "mde80_sd": mde / sd,
    }


def mde80_noncentral_f(df_num: int, df_den: int, alpha: float = 0.05, power: float = 0.80) -> float:
    critical = float(stats.f.ppf(1.0 - alpha, df_num, df_den))

    def objective(noncentrality: float) -> float:
        return float(stats.ncf.sf(critical, df_num, df_den, noncentrality)) - power

    upper = 10.0
    while objective(upper) < 0 and upper < 1e6:
        upper *= 2.0
    return float(optimize.brentq(objective, 0.0, upper, xtol=1e-10, rtol=1e-10))


def fit_cell_type(dataset: str, cell_type: str, frame: pd.DataFrame) -> list[dict[str, Any]]:
    rows = []
    for model, adjusted in (("primary", True), ("unadjusted", False)):
        work = frame.dropna(subset=["composite_score", "age_years", "sex_male", "study_or_pool"]).copy()
        initial_counts = bin_counts(work)
        test_names = ("MID_vs_MATURE", "OLD_vs_MATURE", "OLD_vs_MID", "trend") if adjusted else (
            "MID_vs_MATURE", "OLD_vs_MATURE", "OLD_vs_MID"
        )
        model_rows = [base_row(dataset, cell_type, test, model) for test in test_names]
        for row in model_rows:
            row.update(initial_counts)
        if len(work) < MIN_DONORS or min(initial_counts["n_mature"], initial_counts["n_mid"], initial_counts["n_old"]) < MIN_DONORS:
            for row in model_rows:
                row["status"] = "UNTESTABLE_DONOR_OR_BIN_FLOOR"
            rows.extend(model_rows)
            continue
        sub, active, omitted, status = prune_model_frame(dataset, work, adjusted)
        counts = bin_counts(sub) if status == "evaluated" else initial_counts
        for row in model_rows:
            row.update(counts)
            row["covariates_active"] = ";".join(active) if active else "none"
            row["covariates_omitted"] = ";".join(omitted) if omitted else "none"
        if status != "evaluated":
            for row in model_rows:
                row["status"] = f"UNTESTABLE_{status.upper()}"
            rows.extend(model_rows)
            continue
        if min(counts["n_mature"], counts["n_mid"], counts["n_old"]) < MIN_DONORS:
            for row in model_rows:
                row["status"] = "UNTESTABLE_POST_PRUNE_BIN_FLOOR"
            rows.extend(model_rows)
            continue
        y = sub["composite_score"].to_numpy(dtype=float)
        sd = float(np.std(y, ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            for row in model_rows:
                row["status"] = "UNTESTABLE_ZERO_VARIANCE"
            rows.extend(model_rows)
            continue
        design = pd.DataFrame({
            "intercept": np.ones(len(sub)),
            "MID": sub["age_bin"].eq("MID").to_numpy(dtype=float),
            "OLD": sub["age_bin"].eq("OLD").to_numpy(dtype=float),
        }, index=sub.index)
        for column in active:
            design[column] = sub[column].to_numpy(dtype=float)
        min_df = MIN_RESID_DF[dataset] if adjusted else 5
        if np.linalg.matrix_rank(design.to_numpy()) != design.shape[1] or len(sub) - design.shape[1] < min_df:
            for row in model_rows:
                row["status"] = "UNTESTABLE_CATEGORICAL_RANK_OR_DF"
            rows.extend(model_rows)
            continue
        fit = sm.OLS(y, design).fit()
        results = {
            "MID_vs_MATURE": coefficient_result(fit, "MID", sd),
            "OLD_vs_MATURE": coefficient_result(fit, "OLD", sd),
            "OLD_vs_MID": old_mid_result(fit, sd),
        }
        for row in model_rows[:3]:
            row.update(results[row["test"]])
            row["status"] = "EVALUATED"
        if adjusted:
            ordinal_design = pd.DataFrame({
                "intercept": np.ones(len(sub)),
                "ordinal": sub["ordinal"].to_numpy(dtype=float),
            }, index=sub.index)
            for column in active:
                ordinal_design[column] = sub[column].to_numpy(dtype=float)
            if np.linalg.matrix_rank(ordinal_design.to_numpy()) != ordinal_design.shape[1] or len(sub) - ordinal_design.shape[1] < min_df:
                model_rows[3]["status"] = "UNTESTABLE_TREND_RANK_OR_DF"
            else:
                ordinal_fit = sm.OLS(y, ordinal_design).fit()
                model_rows[3].update(coefficient_result(ordinal_fit, "ordinal", sd))
                model_rows[3]["status"] = "EVALUATED"
        rows.extend(model_rows)
        if adjusted:
            spline_row = base_row(dataset, cell_type, "spline_vs_linear", "spline")
            spline_row.update(counts)
            spline_row["covariates_active"] = ";".join(active) if active else "none"
            spline_row["covariates_omitted"] = ";".join(omitted) if omitted else "none"
            spline = np.asarray(dmatrix("bs(age, df=3, degree=3, include_intercept=False)", {"age": sub["age_years"].to_numpy()}))
            linear_columns = [sub["age_years"].to_numpy(dtype=float)]
            spline_columns = [spline[:, index] for index in range(spline.shape[1])]
            for columns in (linear_columns, spline_columns):
                for column in active:
                    columns.append(sub[column].to_numpy(dtype=float))
            linear_x = sm.add_constant(np.column_stack(linear_columns))
            spline_x = sm.add_constant(np.column_stack(spline_columns))
            if (
                np.linalg.matrix_rank(linear_x) != linear_x.shape[1]
                or np.linalg.matrix_rank(spline_x) != spline_x.shape[1]
                or spline_x.shape[1] - linear_x.shape[1] != 2
                or len(sub) - spline_x.shape[1] < min_df
            ):
                spline_row["status"] = "UNTESTABLE_SPLINE_RANK_OR_DF"
            else:
                linear_fit = sm.OLS(y, linear_x).fit()
                spline_fit = sm.OLS(y, spline_x).fit()
                fvalue, pvalue, _ = spline_fit.compare_f_test(linear_fit)
                df_num = int(spline_fit.df_model - linear_fit.df_model)
                df_den = int(spline_fit.df_resid)
                noncentrality = mde80_noncentral_f(df_num, df_den)
                spline_row.update({
                    "status": "EVALUATED",
                    "estimate": float(fvalue),
                    "pvalue": float(pvalue),
                    "residual_df": df_den,
                    "outcome_sd": sd,
                    "mde80": noncentrality,
                    "mde80_scale": "noncentral_f_noncentrality",
                    "spline_df_num": df_num,
                    "spline_df_den": df_den,
                    "mde80_partial_r2": noncentrality / (noncentrality + df_den),
                })
            rows.append(spline_row)
    return rows


def fit_models(scores: pd.DataFrame) -> pd.DataFrame:
    info_by_dataset = {dataset: load_info(dataset) for dataset in DATASETS}
    rows = []
    for dataset in DATASETS:
        for cell_type in info_by_dataset[dataset]["cell_types"]:
            frame = scores[
                scores["dataset"].eq(dataset)
                & scores["cell_type"].eq(cell_type)
                & scores["n_cells"].ge(MIN_CELLS)
                & scores["age_bin"].isin(BINS)
            ].copy()
            cell_rows = fit_cell_type(dataset, cell_type, frame)
            if not any(row["model"] == "spline" for row in cell_rows):
                primary = next(row for row in cell_rows if row["model"] == "primary" and row["test"] == "MID_vs_MATURE")
                spline_row = base_row(dataset, cell_type, "spline_vs_linear", "spline")
                for key in ("n_donors", "n_mature", "n_mid", "n_old", "covariates_active", "covariates_omitted"):
                    spline_row[key] = primary[key]
                spline_row["status"] = primary["status"]
                cell_rows.append(spline_row)
            rows.extend(cell_rows)
        log(f"{dataset}: categorical, trend, unadjusted, and spline models complete")
    result = pd.DataFrame(rows)
    result["family_size"] = np.nan
    for dataset in DATASETS:
        for test in PRIMARY_FAMILY_TESTS:
            mask = (
                result["dataset"].eq(dataset)
                & result["model"].eq("primary")
                & result["test"].eq(test)
                & result["status"].eq("EVALUATED")
            )
            if mask.any():
                result.loc[mask, "q_bh"] = multipletests(result.loc[mask, "pvalue"], method="fdr_bh")[1]
                result.loc[mask, "family_size"] = int(mask.sum())
    result = result.sort_values(["dataset", "cell_type", "model", "test"]).reset_index(drop=True)
    result["age_window"] = "30-plus; MATURE 30 to 39 baseline"
    write_clean_csv(result, RESULT_PATH)
    return result


def classification(result: pd.DataFrame, dataset: str, cell_type: str) -> str:
    sub = result[
        result["dataset"].eq(dataset)
        & result["cell_type"].eq(cell_type)
        & result["model"].eq("primary")
    ].set_index("test")
    if not all(test in sub.index and sub.loc[test, "status"] == "EVALUATED" for test in PRIMARY_FAMILY_TESTS):
        return "UNTESTABLE"
    mid = sub.loc["MID_vs_MATURE"]
    old_mid = sub.loc["OLD_vs_MID"]
    trend = sub.loc["trend"]
    if mid["q_bh"] < 0.05 and not old_mid["q_bh"] < 0.05:
        return "maturation-concentrated"
    if old_mid["q_bh"] < 0.05 and not mid["q_bh"] < 0.05:
        return "aging-concentrated"
    if trend["q_bh"] < 0.05 and mid["estimate"] * old_mid["estimate"] > 0:
        return "steady trend"
    return "null"


def dataset_verdict(result: pd.DataFrame, dataset: str) -> str:
    classes = [
        classification(result, dataset, cell_type)
        for cell_type in result.loc[result["dataset"].eq(dataset), "cell_type"].unique()
    ]
    if any(value in ("maturation-concentrated", "aging-concentrated") for value in classes):
        return "PROVES"
    if any(value == "steady trend" for value in classes):
        return "DOES NOT PROVE"
    primary = result[
        result["dataset"].eq(dataset)
        & result["model"].eq("primary")
        & result["test"].isin(PRIMARY_FAMILY_TESTS)
        & result["status"].eq("EVALUATED")
    ]
    if primary.empty:
        return "INCONCLUSIVE"
    if (primary["q_bh"] < 0.05).any():
        return "INCONCLUSIVE"
    if (primary["mde80_sd"] <= MDE80_SMALL_SD).all():
        return "DOES NOT PROVE"
    return "INCONCLUSIVE"


def linear_transport_results() -> pd.DataFrame:
    block5 = pd.read_csv(RESULTS_DIR / "followup-mature-block5-gao-yazar-transport.csv")
    block5 = block5[
        block5["arm"].eq("primary")
        & block5["gene_set"].eq("P1m")
        & block5["status"].eq("TESTED")
        & ~block5["unit"].isin(("All_LV_Nuclei", "All_Cells"))
    ].copy()
    if block5.empty:
        raise ValueError("No mature Block 5 primary P1m rows for the linear comparison")
    cohort_map = {"Gao_Heart": "gao", "Yazar_Blood": "yazar"}
    return pd.DataFrame({
        "dataset": block5["cohort"].map(cohort_map),
        "cell_type": block5["unit"],
        "linear_slope": block5["composite_slope"],
        "linear_p": block5["composite_pvalue"],
        "linear_q": block5["composite_q"],
        "linear_source": "followup-mature-block5-gao-yazar-transport.csv",
    })


def fmt(value: Any, digits: int = 5) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (float, np.floating)):
        if value != 0 and abs(value) < 0.0001:
            return f"{value:.3e}"
        return f"{value:.{digits}g}"
    return str(value)


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    def clean(value: Any) -> str:
        return sanitize_dashes(str(value)).replace("|", "/").replace("\n", " ")
    lines = ["| " + " | ".join(clean(value) for value in headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    lines.extend("| " + " | ".join(clean(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def result_row(result: pd.DataFrame, dataset: str, cell_type: str, test: str, model: str = "primary") -> pd.Series:
    selected = result[
        result["dataset"].eq(dataset)
        & result["cell_type"].eq(cell_type)
        & result["test"].eq(test)
        & result["model"].eq(model)
    ]
    if len(selected) != 1:
        raise AssertionError(f"Expected one row for {dataset}, {cell_type}, {test}, {model}")
    return selected.iloc[0]


def generate_report(result: pd.DataFrame, scores: pd.DataFrame) -> dict[str, str]:
    verify_frozen_spec()
    linear = linear_transport_results()
    verdicts = {dataset: dataset_verdict(result, dataset) for dataset in DATASETS}
    overall = "PROVES" if "PROVES" in verdicts.values() else (
        "DOES NOT PROVE" if all(value == "DOES NOT PROVE" for value in verdicts.values()) else "INCONCLUSIVE"
    )
    lines = ["## Results", ""]
    lines.append("All estimates are on the full-depth within-cell-proportion composite scale. MID vs MATURE and OLD vs MATURE are coefficients with MATURE (30 to 39) as reference. OLD vs MID is OLD minus MID from the same categorical fit. Trend is per ordinal-bin increment. NA denotes UNTESTABLE or not applicable. BH q-values exist only for the three pre-stated primary families. Every fitted donor is at least 30 years old; donors under 30 are excluded from every model.")
    lines.append("")
    for dataset in DATASETS:
        info = load_info(dataset)
        coverage = json.loads((CACHE_ROOT / dataset / "scores" / "SCORES_DONE.json").read_text(encoding="utf-8"))
        ds_scores = scores[scores["dataset"].eq(dataset)]
        lines.extend([
            f"### {dataset.capitalize()} data and family sizes",
            "",
            f"The cache supplied {len(ds_scores):,} scored donor-cell-type rows, including {int((ds_scores['age_bin'] == 'UNDER30').sum()):,} UNDER30 rows excluded from models. Scoring used {len(coverage['present_up'])} of {coverage['mapped_human_up']} mapped UP genes and {len(coverage['present_down'])} of {coverage['mapped_human_down']} mapped DOWN genes. Absent UP genes: {', '.join(coverage['absent_up']) or 'none'}. Absent DOWN genes: {', '.join(coverage['absent_down']) or 'none'}. The source cache has {info['n_eligible_cells_ge1500']:,} eligible full-depth cells after the technical mask and 1500-count eligibility filter.",
            "",
        ])
        family_rows = []
        for test in PRIMARY_FAMILY_TESTS:
            family = result[
                result["dataset"].eq(dataset) & result["model"].eq("primary") & result["test"].eq(test)
            ]
            evaluated = family[family["status"].eq("EVALUATED")]
            family_rows.append([test, len(evaluated), int((evaluated["q_bh"] < 0.05).sum())])
        lines.append(markdown_table(["Primary family", "Evaluated cell types", "BH q < 0.05"], family_rows))
        lines.append("")

    lines.extend(["### Primary categorical and trend results", ""])
    for dataset in DATASETS:
        lines.extend([f"#### {dataset.capitalize()}", ""])
        table_rows = []
        for cell_type in result.loc[result["dataset"].eq(dataset), "cell_type"].unique():
            mid = result_row(result, dataset, cell_type, "MID_vs_MATURE")
            old = result_row(result, dataset, cell_type, "OLD_vs_MATURE")
            old_mid = result_row(result, dataset, cell_type, "OLD_vs_MID")
            trend = result_row(result, dataset, cell_type, "trend")
            status = mid["status"] if mid["status"] != "EVALUATED" else trend["status"]
            table_rows.append([
                cell_type, status,
                f"{mid['n_mature']}/{mid['n_mid']}/{mid['n_old']}",
                f"{fmt(mid['estimate'])} ({fmt(mid['ci95_low'])} to {fmt(mid['ci95_high'])}); p={fmt(mid['pvalue'])}; q={fmt(mid['q_bh'])}",
                f"{fmt(old['estimate'])} ({fmt(old['ci95_low'])} to {fmt(old['ci95_high'])}); p={fmt(old['pvalue'])}",
                f"{fmt(old_mid['estimate'])} ({fmt(old_mid['ci95_low'])} to {fmt(old_mid['ci95_high'])}); p={fmt(old_mid['pvalue'])}; q={fmt(old_mid['q_bh'])}",
                f"{fmt(trend['estimate'])} ({fmt(trend['ci95_low'])} to {fmt(trend['ci95_high'])}); p={fmt(trend['pvalue'])}; q={fmt(trend['q_bh'])}",
                f"{fmt(mid['mde80_sd'])} / {fmt(old_mid['mde80_sd'])} / {fmt(trend['mde80_sd'])}",
            ])
        lines.append(markdown_table([
            "Cell type", "Status", "N MAT/MID/OLD", "MID vs MATURE estimate (95 percent CI); p; q",
            "OLD vs MATURE estimate (95 percent CI); p", "OLD vs MID estimate (95 percent CI); p; q",
            "Trend estimate (95 percent CI); p; q", "MDE80 SD MID / OLD-MID / trend",
        ], table_rows))
        lines.append("")

    lines.extend(["### Unadjusted categorical sensitivity", ""])
    for dataset in DATASETS:
        unadjusted = result[result["dataset"].eq(dataset) & result["model"].eq("unadjusted")]
        table_rows = []
        for cell_type in unadjusted["cell_type"].unique():
            mid = result_row(result, dataset, cell_type, "MID_vs_MATURE", "unadjusted")
            old = result_row(result, dataset, cell_type, "OLD_vs_MATURE", "unadjusted")
            old_mid = result_row(result, dataset, cell_type, "OLD_vs_MID", "unadjusted")
            table_rows.append([
                cell_type, mid["status"], f"{mid['n_mature']}/{mid['n_mid']}/{mid['n_old']}",
                f"{fmt(mid['estimate'])}; p={fmt(mid['pvalue'])}",
                f"{fmt(old['estimate'])}; p={fmt(old['pvalue'])}",
                f"{fmt(old_mid['estimate'])}; p={fmt(old_mid['pvalue'])}",
            ])
        lines.extend([f"#### {dataset.capitalize()}", ""])
        lines.append(markdown_table(["Cell type", "Status", "N MAT/MID/OLD", "MID vs MATURE", "OLD vs MATURE", "OLD vs MID"], table_rows))
        lines.append("")

    lines.extend(["### Cubic B-spline sensitivity", ""])
    spline_rows = []
    for row in result[result["model"].eq("spline")].itertuples(index=False):
        spline_rows.append([
            row.dataset, row.cell_type, row.status, fmt(row.estimate), fmt(row.pvalue),
            fmt(row.spline_df_num), fmt(row.spline_df_den), fmt(row.mde80), fmt(row.mde80_partial_r2),
        ])
    lines.append(markdown_table(["Dataset", "Cell type", "Status", "F", "p", "df num", "df den", "MDE80 noncentrality", "MDE80 partial R2"], spline_rows))
    lines.append("")

    lines.extend(["### Classification and comparison with prior linear models", ""])
    lines.append("Transition classifications require the stated BH pattern in the batch-adjusted primary model. A significant categorical step with a null trend is a step pattern, not a steady slope. Linear flagged means the mature Block 5 linear model had composite BH q < 0.05 in its primary P1m arm.")
    lines.append("")
    for dataset in DATASETS:
        table_rows = []
        for cell_type in result.loc[result["dataset"].eq(dataset), "cell_type"].unique():
            label = classification(result, dataset, cell_type)
            prior = linear[linear["dataset"].eq(dataset) & linear["cell_type"].eq(cell_type)]
            if prior.empty:
                slope = q = np.nan
                flagged = "not tested"
            else:
                slope = float(prior.iloc[0]["linear_slope"])
                q = float(prior.iloc[0]["linear_q"])
                flagged = "yes" if q < 0.05 else "no"
            table_rows.append([cell_type, label, fmt(slope), fmt(q), flagged])
        lines.extend([f"#### {dataset.capitalize()}", ""])
        lines.append(markdown_table(["Cell type", "Classification", "Prior linear slope per decade", "Prior linear q", "Linear flagged"], table_rows))
        lines.append("")

    for dataset in DATASETS:
        classes = [
            classification(result, dataset, cell_type)
            for cell_type in result.loc[result["dataset"].eq(dataset), "cell_type"].unique()
        ]
        counts = {label: classes.count(label) for label in ("maturation-concentrated", "aging-concentrated", "steady trend", "null", "UNTESTABLE")}
        lines.append(f"- {dataset.capitalize()} verdict: {verdicts[dataset]}. Classifications: {counts['maturation-concentrated']} maturation-concentrated, {counts['aging-concentrated']} aging-concentrated, {counts['steady trend']} steady trend, {counts['null']} null, {counts['UNTESTABLE']} UNTESTABLE.")
    lines.extend([
        "",
        f"Block-level verdict: {overall}. This verdict applies to the P1m score in Gao heart and Yazar blood under the pre-stated categorical and trend rules.",
        "",
        "The unadjusted and spline models are sensitivities and do not determine the primary classification. The spline F test addresses departure from continuous linearity but is not substituted for the three categorical BH families.",
        "",
        "## Reproducibility",
        "",
        "Run `~/miniconda3/envs/agescape/bin/python scripts/run_mature_block8_categorical_signature.py --freeze-spec`, then `~/miniconda3/envs/agescape/bin/python scripts/run_mature_block8_categorical_signature.py`. The second command requires the frozen report prefix, reuses the Block 1 cache, writes both CSVs, and regenerates this results section.",
        "",
        "One deep review round was performed after output completion. No other analysis block was started.",
    ])
    report = frozen_spec() + "\n" + "\n".join(lines) + "\n"
    write_clean_markdown(report, REPORT_PATH)
    return {"gao": verdicts["gao"], "yazar": verdicts["yazar"], "overall": overall}


def validate_outputs(result: pd.DataFrame, scores: pd.DataFrame) -> None:
    expected_cell_types = sum(len(load_info(dataset)["cell_types"]) for dataset in DATASETS)
    if len(result) != expected_cell_types * 8:
        raise AssertionError(f"Expected {expected_cell_types * 8} model rows, found {len(result)}")
    if scores.empty or result.empty:
        raise AssertionError("Empty output")
    if scores["n_cells"].isna().any() or (scores["n_cells"] < 1).any():
        raise AssertionError("Missing or invalid donor-cell-type cell counts")
    evaluated = result[result["status"].eq("EVALUATED")]
    if not np.isfinite(evaluated[["estimate", "pvalue", "mde80"]].to_numpy(dtype=float)).all():
        raise AssertionError("Non-finite evaluated result")
    if (evaluated["n_donors"] < MIN_DONORS).any():
        raise AssertionError("Evaluated model below donor floor")
    categorical = evaluated[evaluated["test"].isin(("MID_vs_MATURE", "OLD_vs_MATURE", "OLD_vs_MID"))]
    if (categorical[["n_mature", "n_mid", "n_old"]].min(axis=1) < MIN_DONORS).any():
        raise AssertionError("Evaluated categorical model below bin floor")
    for dataset in DATASETS:
        for test in PRIMARY_FAMILY_TESTS:
            family = result[
                result["dataset"].eq(dataset) & result["model"].eq("primary")
                & result["test"].eq(test) & result["status"].eq("EVALUATED")
            ]
            if len(family) and not np.isfinite(family["q_bh"]).all():
                raise AssertionError("Missing BH q in evaluated primary family")
    for path in (RESULT_PATH, DONOR_PATH, REPORT_PATH):
        text = path.read_text(encoding="utf-8")
        if "\u2013" in text or "\u2014" in text:
            raise AssertionError(f"Forbidden dash in {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-spec", action="store_true")
    parser.add_argument("--stage", choices=("all", "scores", "models"), default="all")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.freeze_spec:
        freeze_spec()
        return
    verify_frozen_spec()
    program = verify_program()
    runtimes: list[dict[str, Any]] = []
    overall_start = time.monotonic()
    overall_start_utc = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()) + "Z"
    if args.stage in ("all", "scores"):
        stage_start = time.monotonic()
        stage_start_utc = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()) + "Z"
        for dataset in DATASETS:
            run_scores(dataset, program)
        build_donor_scores()
        runtimes.append({"model": MODEL, "stage": "scores", "start_utc": stage_start_utc,
                         "end_utc": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()) + "Z",
                         "duration_seconds": time.monotonic() - stage_start})
    if args.stage in ("all", "models"):
        stage_start = time.monotonic()
        stage_start_utc = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()) + "Z"
        scores = pd.read_csv(DONOR_PATH, float_precision="round_trip")
        result = fit_models(scores)
        verdicts = generate_report(result, scores)
        validate_outputs(result, scores)
        runtimes.append({"model": MODEL, "stage": "models", "start_utc": stage_start_utc,
                         "end_utc": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()) + "Z",
                         "duration_seconds": time.monotonic() - stage_start})
        log(f"Block 8 verdicts: {verdicts}")
    runtimes.append({"model": MODEL, "stage": "total", "start_utc": overall_start_utc,
                     "end_utc": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()) + "Z",
                     "duration_seconds": time.monotonic() - overall_start})
    write_clean_csv(pd.DataFrame(runtimes), RUNTIME_PATH)


if __name__ == "__main__":
    main()
