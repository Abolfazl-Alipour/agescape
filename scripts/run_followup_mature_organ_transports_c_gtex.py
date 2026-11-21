#!/usr/bin/env python3
"""Phase 2 Audit Item 5 Batch C: GTEx v8 mature (30-plus bins) transport rerun.

Frozen protocol, same as Batch A (scripts/run_followup_mature_organ_transports_a.py):
  - Mature signatures from Item 1: P1m 1,914 mouse genes (1,119 up, 795 down),
    1,724 human-mapped (1,013 up, 711 down); P0m 2,598 mouse genes (1,474 up,
    1,124 down), 2,049 human-mapped (1,157 up, 892 down). Loaded and asserted
    via load_verified_mature_programs() imported from the Batch A script.
  - Seed 20261105. 2000-draw expression-decile-matched null as primary.
  - Donors 30-plus only. Donor is the inferential unit, never the cell.
  - Explicit multiple-testing families declared before BH correction; BH by
    the Batch A script's bh_checked (explicit BH cross-checked against
    statsmodels, assert diff <= 1e-12).

GTEx-specific adaptations (age is ordinal bins, not exact years):
  - GTEx v8 open-access AGE has only decade bins 20-29, 30-39, 40-49, 50-59,
    60-69, 70-79. Restrict to bins 30-39 through 70-79. Do NOT impute exact
    ages. The ordinal predictor age_bin_order (1.0 through 6.0 in the
    original) keeps its original values 2.0 through 6.0, so a slope stays
    per decade bin and stays comparable with the original.
  - Bulk RNA-seq only: gene-level age tests and signature tests ONLY. NO
    single-cell fraction/CLR models and NO cross-tissue concordance, exactly
    as the original GTEx report states, because bulk expression changes
    conflate within-cell aging with shifts in cell-type proportions.
  - Track 1: OLS log2(CPM+1) per-gene tests plus P1m/P0m signature transport
    in all 49 original groups (47 primary tissues plus 2 isolated cell
    lines) with mature n_complete >= 50, same covariate design as the
    original script.
  - Track 2: PyDESeq2 NB-GLM per-gene tests plus P1m/P0m signature transport
    in the 6 focal tissues (Whole Blood, Muscle - Skeletal, Lung,
    Kidney - Cortex, Ovary, Testis). Primary continuous arm only; the 5
    frozen sensitivity arms are out of scope for this item and are not rerun.
  - Family B/C functions below mirror the Batch A run_family_b_mature and
    run_family_c_mature logic (same fit_ols, same fast_random_signature_scores
    with parity verification, same detection filter) with the maturity gate
    expressed on age_bin_order (>= 2.0) instead of an exact-age column,
    because GTEx has no exact ages to assert on.

Multiple-testing families (declared before testing):
  1. T1B-primary: all mature gene OLS p-values across the 47 primary tissues.
  2. T1B-cellline: all mature gene OLS p-values across the 2 cell lines.
  3. T1C-primary-P1m / T1C-primary-P0m: 47 composite p-values and 47
     empirical p-values each (BH separately per p-value type).
  4. T1C-cellline-P1m / T1C-cellline-P0m: 2 composite and 2 empirical each.
  5. T2B-global: all mature NB-GLM gene p-values across the 6 tissues.
  6. T2C-P1m / T2C-P0m: 6 composite Wald p-values and 6 empirical p-values
     each.
  Within-tissue q_bh_within_group is also reported per tissue, as in the
  original, but verdicts use the family q values above.

Outputs (results/followup-mature-organ-transports-c-gtex*.csv, report
reports/followup-mature-organ-transports-batch-c-gtex.md). Nothing under
report/ is touched. Nothing is submitted anywhere.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
from typing import Any

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    import os

    os.environ[_var] = "1"

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_gtex_v8 import (  # noqa: E402
    AGE_BIN_TO_ORDER,
    EXPECTED_GCT_BYTES,
    EXPECTED_GCT_MD5_B64,
    ISOLATED_CELL_LINES,
    TRACK2_PRIMARY_TISSUES,
    build_clean_covariate_design,
    compute_p1_on_nb_results,
    fit_pydeseq2_continuous,
    load_gtex_sample_subject_metadata,
    stream_gct_to_memmap,
    verify_gct_file,
)
from scripts.run_followup_mature_organ_transports_a import (  # noqa: E402
    ALPHA,
    ORIG_NULL_SEED,
    P1_ITERATIONS,
    bh_checked,
    fast_random_signature_scores,
    load_verified_mature_programs,
    make_row,
    mde80_from_null_sd,
    mde80_from_se,
    order_cohort_df,
    transport_class,
    verdict_standard,
)
from src.followup_muscle import fit_ols, matched_null_effect_pvalue  # noqa: E402
from src.followup_organs import (  # noqa: E402
    DATA_TMP_ROOT,
    assert_no_forbidden_dashes,
    counts_to_log2cpm,
    extended_technical_mask,
    sanitize_dashes,
    select_design_columns,
    write_clean_csv,
    write_clean_markdown,
)

RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
GCT_URL = "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
GCT_PATH = DATA_ROOT / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
WORK_DIR = DATA_TMP_ROOT / "gtex_v8"

MATURE_BINS = ("30-39", "40-49", "50-59", "60-69", "70-79")
EXPECTED_SUBJECT_BIN_COUNTS = {"20-29": 84, "30-39": 78, "40-49": 153, "50-59": 315, "60-69": 317, "70-79": 33}
EXPECTED_GATE = {"n_genes": 56200, "n_samples": 17382, "nonfinite": 0, "negative": 0, "noninteger_lines": 0, "total_nonzero_entries": 515414926}
MIN_DONORS_TRACK1 = 50
MIN_CPM = 1.0
MIN_EXPR_FRAC = 0.5
MIN_RESIDUAL_DF = 10
N_CPUS_TRACK2 = 6
MODEL_USED = "Muse Spark direct subagent run (no CLI model)"
COHORT = "GTEx_v8"

PREFIX = "followup-mature-organ-transports-c-gtex"
MAIN_CSV = RESULTS_DIR / f"{PREFIX}.csv"
BH_CHECK_CSV = RESULTS_DIR / f"{PREFIX}-bh-check.csv"
RUNTIME_CSV = RESULTS_DIR / f"{PREFIX}-runtime.csv"
T1_GENES_CSV = RESULTS_DIR / f"{PREFIX}-track1-genes.csv"
T1_GENES_SUMMARY_CSV = RESULTS_DIR / f"{PREFIX}-track1-genes-summary.csv"
T1_P1_CSV = RESULTS_DIR / f"{PREFIX}-track1-p1.csv"
T2_GENES_CSV = RESULTS_DIR / f"{PREFIX}-track2-nbglm-genes.csv"
T2_SUMMARY_CSV = RESULTS_DIR / f"{PREFIX}-track2-nbglm-summary.csv"
T2_P1_CSV = RESULTS_DIR / f"{PREFIX}-track2-p1.csv"
REPORT_MD = REPORTS_DIR / "followup-mature-organ-transports-batch-c-gtex.md"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def peak_rss_mb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0


def append_runtime(step: str, start_utc: str, t0: float, status: str, notes: str = "") -> None:
    row = pd.DataFrame(
        [
            {
                "cohort": COHORT,
                "step": step,
                "start_utc": start_utc,
                "end_utc": utc_now(),
                "duration_seconds": round(time.monotonic() - t0, 3),
                "peak_rss_mb": round(peak_rss_mb(), 2),
                "model_used": MODEL_USED,
                "status": status,
                "notes": sanitize_dashes(notes),
            }
        ]
    )
    if RUNTIME_CSV.exists():
        prev = pd.read_csv(RUNTIME_CSV)
        prev = prev[~((prev["cohort"] == COHORT) & (prev["step"] == step))]
        out = pd.concat([prev, row], ignore_index=True)
    else:
        out = row
    write_clean_csv(out, RUNTIME_CSV)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


# =============================================================================
# GTEx Family B (per-gene OLS) and Family C (signature transport), mirroring
# Batch A run_family_b_mature / run_family_c_mature with a bin-based gate.
# =============================================================================

def gtex_family_b(
    log2cpm_df: pd.DataFrame,
    meta: pd.DataFrame,
    candidate_covariates: tuple[str, ...],
    min_donors: int = MIN_DONORS_TRACK1,
    min_cpm: float = MIN_CPM,
    min_expr_frac: float = MIN_EXPR_FRAC,
    min_residual_df: int = MIN_RESIDUAL_DF,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Donor-level per-gene OLS on age_bin_order for one tissue (mature only)."""
    assert float(meta["age_bin_order"].min()) >= 2.0, "Under-30 bin leaked into Family B"
    assert meta["age_bin_order"].nunique() >= 2, "Fewer than 2 mature age bins present"
    common_donors = [d for d in meta.index if d in log2cpm_df.columns]
    if len(common_donors) < min_donors:
        return pd.DataFrame(), {"status": f"untestable_low_donors({len(common_donors)}<{min_donors})", "n_donors": len(common_donors), "n_genes_tested": 0}
    mat_df = log2cpm_df[common_donors]
    m_meta = meta.loc[common_donors]
    log2_thresh = float(np.log2(min_cpm + 1.0))
    min_req = max(3, int(np.ceil(min_expr_frac * len(common_donors))))
    keep_mask = (mat_df.to_numpy(dtype=np.float64) >= log2_thresh).sum(axis=1) >= min_req
    if int(keep_mask.sum()) == 0:
        return pd.DataFrame(), {"status": "untestable_no_expressed_genes", "n_donors": len(common_donors), "n_genes_tested": 0}
    mat_sub = mat_df.loc[keep_mask]
    try:
        design, names, omitted = select_design_columns(
            m_meta, age_col="age_bin_order", candidate_covariates=candidate_covariates, min_residual_df=min_residual_df
        )
    except ValueError as exc:
        return pd.DataFrame(), {"status": f"untestable_design({exc})", "n_donors": len(common_donors), "n_genes_tested": 0}
    outcomes = mat_sub.to_numpy(dtype=np.float64).T
    ols = fit_ols(outcomes, design, effect_column=1)
    t_stats = np.where(ols["se"] > 0, ols["effect"] / ols["se"], np.nan)
    df_val = int(ols["residual_df"][0])
    mde_factor = float(scipy_stats.t.ppf(0.975, df_val) + scipy_stats.t.ppf(0.80, df_val))
    mde_vec = mde_factor * ols["se"]
    q_within, _, _ = bh_checked(ols["pvalue"])
    gene_df = pd.DataFrame(
        {
            "gene": mat_sub.index.astype(str),
            "n_donors": int(len(common_donors)),
            "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
            "omitted_covariates": ";".join(omitted) if omitted else "none",
            "mean_log2cpm": outcomes.mean(axis=0),
            "effect": ols["effect"],
            "se": ols["se"],
            "ci95_low": ols["ci95_low"],
            "ci95_high": ols["ci95_high"],
            "t_stat": t_stats,
            "pvalue": ols["pvalue"],
            "residual_df": ols["residual_df"].astype(int),
            "r2": ols["r2"],
            "q_bh_within_group": q_within,
            "mde80": mde_vec,
        }
    )
    min_idx = int(np.nanargmin(ols["pvalue"]))
    summary = {
        "status": "tested",
        "n_donors": int(len(common_donors)),
        "age_bins_present": int(m_meta["age_bin"].nunique()),
        "n_genes_tested": int(len(gene_df)),
        "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
        "omitted_covariates": ";".join(omitted) if omitted else "none",
        "residual_df": df_val,
        "median_se": float(np.nanmedian(ols["se"])),
        "median_mde80": float(np.nanmedian(mde_vec)),
        "top_gene": str(mat_sub.index[min_idx]),
        "top_gene_effect": float(ols["effect"][min_idx]),
        "top_gene_se": float(ols["se"][min_idx]),
        "min_pvalue": float(np.nanmin(ols["pvalue"])),
        "n_sig_within_q05": int((q_within < ALPHA).sum()),
    }
    return gene_df, summary


def gtex_family_c(
    log2cpm_df: pd.DataFrame,
    meta: pd.DataFrame,
    program: dict[str, list[str]],
    candidate_covariates: tuple[str, ...],
    min_donors: int = MIN_DONORS_TRACK1,
    min_cpm: float = MIN_CPM,
    min_expr_frac: float = MIN_EXPR_FRAC,
    min_residual_df: int = MIN_RESIDUAL_DF,
    iterations: int = P1_ITERATIONS,
    seed: int = ORIG_NULL_SEED,
    verify_parity: bool = False,
) -> tuple[dict[str, Any], str]:
    """Signature transport for one tissue and one mature program tier."""
    assert float(meta["age_bin_order"].min()) >= 2.0, "Under-30 bin leaked into Family C"
    up_set = set(program["up"])
    down_set = set(program["down"])
    log2_thresh = float(np.log2(min_cpm + 1.0))
    common_donors = [d for d in meta.index if d in log2cpm_df.columns]
    if len(common_donors) < min_donors:
        return {}, f"n_donors={len(common_donors)} < min_donors={min_donors}"
    mat_df = log2cpm_df[common_donors]
    m_meta = meta.loc[common_donors]
    min_req = max(3, int(np.ceil(min_expr_frac * len(common_donors))))
    keep_mask = (mat_df.to_numpy(dtype=np.float64) >= log2_thresh).sum(axis=1) >= min_req
    mat_sub = mat_df.loc[keep_mask]
    if mat_sub.shape[0] < 200:
        return {}, f"n_expressed_pool={mat_sub.shape[0]} < 200"
    genes = mat_sub.index.astype(str).tolist()
    up_idx = np.array([i for i, g in enumerate(genes) if g in up_set], dtype=int)
    down_idx = np.array([i for i, g in enumerate(genes) if g in down_set], dtype=int)
    if len(up_idx) == 0 or len(down_idx) == 0:
        return {}, f"detected up={len(up_idx)}, down={len(down_idx)}"
    try:
        design, names, omitted = select_design_columns(
            m_meta, age_col="age_bin_order", candidate_covariates=candidate_covariates, min_residual_df=min_residual_df
        )
    except ValueError as exc:
        return {}, f"design_error: {exc}"
    mat = mat_sub.to_numpy(dtype=np.float64)
    up_score = mat[up_idx, :].mean(axis=0)
    down_score = mat[down_idx, :].mean(axis=0)
    comp_score = up_score - down_score
    ols = fit_ols(np.column_stack([up_score, down_score, comp_score]), design, effect_column=1)
    try:
        null_scores, min_decile_margin = fast_random_signature_scores(
            mat, up_idx, down_idx, iterations=iterations, seed=seed, verify_parity=verify_parity
        )
    except ValueError as exc:
        return {}, f"matched_null_pool_error: {exc}"
    null_ols = fit_ols(null_scores.T, design, effect_column=1)
    null_effects = null_ols["effect"]
    p_emp = matched_null_effect_pvalue(float(ols["effect"][2]), null_effects)
    return {
        "n_donors": int(len(common_donors)),
        "n_expressed_pool": int(len(genes)),
        "n_up_detected": int(len(up_idx)),
        "n_down_detected": int(len(down_idx)),
        "n_total_detected": int(len(up_idx) + len(down_idx)),
        "min_decile_margin": int(min_decile_margin),
        "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
        "omitted_covariates": ";".join(omitted) if omitted else "none",
        "up_effect": float(ols["effect"][0]),
        "up_pvalue": float(ols["pvalue"][0]),
        "down_effect": float(ols["effect"][1]),
        "down_pvalue": float(ols["pvalue"][1]),
        "composite_effect": float(ols["effect"][2]),
        "composite_se": float(ols["se"][2]),
        "composite_ci95_low": float(ols["ci95_low"][2]),
        "composite_ci95_high": float(ols["ci95_high"][2]),
        "composite_pvalue": float(ols["pvalue"][2]),
        "null_mean_effect": float(null_effects.mean()),
        "null_sd_effect": float(null_effects.std(ddof=1)),
        "null_iterations": int(iterations),
        "empirical_null_pvalue": float(p_emp),
        "residual_df": int(ols["residual_df"][2]),
    }, ""


def nb_p1_mature(
    nb_res: pd.DataFrame,
    gene_desc_map: dict[str, str],
    program: dict[str, list[str]],
    iterations: int = P1_ITERATIONS,
) -> dict[str, Any]:
    """P1m/P0m effect, Stouffer z, and 2000-draw expression-matched null on NB-GLM results.

    Mirrors scripts/run_followup_gtex_v8.py compute_p1_on_nb_results exactly
    (same decile pools from log2(baseMean+1), same seed, same
    matched_null_effect_pvalue) and additionally returns the null mean/sd.
    """
    work = nb_res.copy()
    work["symbol"] = work["ensembl_id"].map(gene_desc_map).fillna("")
    work = work.loc[(work["symbol"] != "") & work["log2FoldChange"].notna() & work["lfcSE"].notna() & (work["lfcSE"] > 0)].copy()
    work = work.sort_values("baseMean", ascending=False).drop_duplicates(subset=["symbol"], keep="first")
    z_all = (work["log2FoldChange"] / work["lfcSE"]).to_numpy(dtype=float)
    b_all = work["log2FoldChange"].to_numpy(dtype=float)
    mean_expr = np.log2(np.maximum(work["baseMean"].to_numpy(dtype=float), 0.0) + 1.0)
    sym_to_idx = {s: i for i, s in enumerate(work["symbol"].tolist())}
    up_idx = np.asarray([sym_to_idx[g] for g in program["up"] if g in sym_to_idx], dtype=int)
    down_idx = np.asarray([sym_to_idx[g] for g in program["down"] if g in sym_to_idx], dtype=int)
    eff = float(np.mean(b_all[up_idx]) - np.mean(b_all[down_idx]))
    z_signed = np.concatenate([z_all[up_idx], -z_all[down_idx]])
    comp_z = float(z_signed.sum() / np.sqrt(len(z_signed)))
    comp_p = float(2.0 * scipy_stats.norm.sf(abs(comp_z)))
    rng = np.random.default_rng(ORIG_NULL_SEED)
    deciles = pd.Series(pd.qcut(mean_expr, 10, labels=False, duplicates="drop")).fillna(-1).to_numpy()
    available = np.ones(len(b_all), dtype=bool)
    available[np.concatenate([up_idx, down_idx])] = False
    up_dec_counts = pd.Series(deciles[up_idx]).value_counts().to_dict()
    down_dec_counts = pd.Series(deciles[down_idx]).value_counts().to_dict()
    pools_by_dec = {dec: np.flatnonzero((deciles == dec) & available) for dec in set(up_dec_counts) | set(down_dec_counts)}
    null_draws = np.empty(iterations, dtype=float)
    for it in range(iterations):
        r_up: list[np.ndarray] = []
        r_dn: list[np.ndarray] = []
        for dec, cnt in up_dec_counts.items():
            r_up.append(rng.choice(pools_by_dec[dec], size=int(cnt), replace=False))
        for dec, cnt in down_dec_counts.items():
            r_dn.append(rng.choice(pools_by_dec[dec], size=int(cnt), replace=False))
        null_draws[it] = float(b_all[np.concatenate(r_up)].mean() - b_all[np.concatenate(r_dn)].mean())
    emp_p = matched_null_effect_pvalue(eff, null_draws)
    return {
        "nb_p1_n_up": int(len(up_idx)),
        "nb_p1_n_down": int(len(down_idx)),
        "nb_p1_effect_lfc": eff,
        "nb_p1_composite_z": comp_z,
        "nb_p1_composite_p": comp_p,
        "nb_p1_empirical_p": emp_p,
        "nb_p1_null_mean": float(null_draws.mean()),
        "nb_p1_null_sd": float(null_draws.std(ddof=1)),
        "nb_p1_null_iterations": int(iterations),
    }

# =============================================================================
# MAIN
# =============================================================================

def load_original_tables() -> dict[str, pd.DataFrame]:
    orig: dict[str, pd.DataFrame] = {}
    orig["t1_audit"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-tissue-audit.csv").set_index("SMTSD")
    orig["t1_genes_sum"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-track1-genes-summary.csv").set_index("group")
    orig["t1_top"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-track1-top-genes.csv")
    orig["t1_p1"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-track1-p1-signature.csv").set_index("group")
    orig["t2_sum"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-track2-nbglm-summary.csv").set_index("tissue")
    orig["t2_top"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-track2-nbglm-top-genes.csv")
    orig["t2_p1"] = pd.read_csv(RESULTS_DIR / "followup-gtex-v8-track2-p1-nbglm.csv").set_index("tissue")
    # Original per-tissue minimum family q: BH q is monotone in p, so the
    # minimum-p gene carries the minimum family q.
    t1_minq = (
        orig["t1_top"].sort_values("pvalue").groupby("group", sort=False).first()["q_bh_family"].to_dict()
    )
    orig["t1_min_family_q"] = t1_minq
    t2_minq = (
        orig["t2_top"].sort_values("pvalue").groupby("tissue", sort=False).first()["q_bh_track2_global"].to_dict()
    )
    orig["t2_min_global_q"] = t2_minq
    return orig


def run_track1(
    mm: np.memmap,
    genes_df: pd.DataFrame,
    sample_id_to_col: dict[str, int],
    merged_meta: pd.DataFrame,
    programs: dict[str, dict[str, list[str]]],
    orig: dict[str, pd.DataFrame],
    groups_49: list[str],
) -> dict[str, Any]:
    """Track 1 OLS scan on mature bins for the 49 original groups."""
    from scipy import sparse as sp

    tech_mask = extended_technical_mask(genes_df["Description"].astype(str).to_numpy())
    nontech_idx = np.flatnonzero(~tech_mask)
    nontech_syms = genes_df["Description"].astype(str).to_numpy()[nontech_idx]
    uniq_syms, sym_inv = np.unique(nontech_syms, return_inverse=True)
    valid_sym = np.array([s != "" and s.lower() not in {"nan", "na"} for s in uniq_syms], dtype=bool)
    sym_agg = sp.csr_matrix(
        (np.ones(len(nontech_idx), dtype=np.float64), (sym_inv, np.arange(len(nontech_idx)))),
        shape=(len(uniq_syms), len(nontech_idx)),
    )[valid_sym, :]
    clean_symbols = uniq_syms[valid_sym]
    sym_index = pd.Index(clean_symbols, dtype=str)

    gene_frames: list[pd.DataFrame] = []
    gene_summaries: list[dict[str, Any]] = []
    p1_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    parity_verified = {"P1m": False, "P0m": False}

    for smtsd in groups_49:
        t0 = time.monotonic()
        s0 = utc_now()
        cohort_class = str(orig["t1_audit"].loc[smtsd, "cohort_class"])
        orig_n_samples = int(orig["t1_audit"].loc[smtsd, "n_samples"])
        orig_n_complete = int(orig["t1_audit"].loc[smtsd, "n_complete"])
        t_meta_full = merged_meta.loc[merged_meta["SMTSD"] == smtsd].copy()
        n_20_29 = int((t_meta_full["age_bin"] == "20-29").sum())
        t_meta = t_meta_full.loc[t_meta_full["age_bin"].isin(MATURE_BINS)].copy()
        assert len(t_meta) == orig_n_samples - n_20_29
        assert not (t_meta["SUBJID"].duplicated().any()), f"Duplicate donor in {smtsd}"
        assert set(t_meta["age_bin"].unique()) <= set(MATURE_BINS)

        sids = list(t_meta.index.astype(str))
        col_indices = [sample_id_to_col[sid] for sid in sids]
        # Column-first fancy indexing: same values as the original script's
        # row-first mm[nontech_idx, :][:, col_indices], but it materialises
        # only this tissue's columns instead of all 17,382 sample columns.
        sub_raw = np.asarray(mm[:, col_indices][nontech_idx, :], dtype=np.float64)
        sym_counts = sym_agg @ sub_raw
        counts_df = pd.DataFrame(sym_counts, index=sym_index, columns=sids)
        log2cpm, lib_sizes = counts_to_log2cpm(counts_df)
        t_meta["lib_size"] = lib_sizes.loc[sids].to_numpy(dtype=float)
        t_meta["log_lib"] = np.log1p(t_meta["lib_size"].to_numpy(dtype=float))

        cov_design = build_clean_covariate_design(
            t_meta,
            numeric_candidates=["sex_male", "log_lib", "SMRIN", "SMTSISCH_hr", "SMEXNCRT"],
            categorical_candidates=["SMCENTER", "SMNABTCHT"],
            age_col="age_bin_order",
        )
        m_comp = cov_design.meta
        e_comp = log2cpm[list(m_comp.index)]
        n_mature_complete = int(cov_design.n_samples_complete)
        n_excluded = orig_n_complete - n_mature_complete

        audit_rows.append(
            {
                "SMTSD": smtsd,
                "cohort_class": cohort_class,
                "n_samples_original": orig_n_samples,
                "n_samples_mature_bins": int(len(t_meta)),
                "n_20_29_excluded": n_20_29,
                "n_complete_mature": n_mature_complete,
                "n_excluded_vs_original_complete": n_excluded,
                "n_male_mature": int(m_comp["sex_male"].sum()),
                "n_female_mature": int((1 - m_comp["sex_male"]).sum()),
                "age_bins_present": int(m_comp["age_bin"].nunique()),
                "covariates_used": ";".join(cov_design.all_cov_cols),
                "included_in_analysis": bool(n_mature_complete >= MIN_DONORS_TRACK1),
            }
        )

        if n_mature_complete < MIN_DONORS_TRACK1:
            gene_summaries.append(
                {"group": smtsd, "status": f"untestable_lt50_mature(n={n_mature_complete})", "n_donors": n_mature_complete, "n_genes_tested": 0}
            )
            for tier in ("P1m", "P0m"):
                p1_rows.append({"group": smtsd, "tier": tier, "status": f"untestable_lt50_mature(n={n_mature_complete})", "n_donors": n_mature_complete})
            append_runtime(f"track1_{smtsd}", s0, t0, "ok", f"UNTESTABLE mature n_complete={n_mature_complete} < 50")
            continue

        cov_tuple = tuple(cov_design.all_cov_cols)
        g_df, g_sum = gtex_family_b(e_comp, m_comp, candidate_covariates=cov_tuple)
        g_df.insert(0, "group", smtsd)
        g_df["predictor"] = "age_bin_order"
        g_df["cohort_class"] = cohort_class
        if g_sum["status"] == "tested":
            gene_frames.append(g_df)
        g_sum["group"] = smtsd
        g_sum["cohort_class"] = cohort_class
        g_sum["n_excluded_vs_original_complete"] = n_excluded
        g_sum["n_20_29_excluded"] = n_20_29
        gene_summaries.append(g_sum)

        for tier in ("P1m", "P0m"):
            verify = not parity_verified[tier]
            prow, reason = gtex_family_c(
                e_comp, m_comp, program=programs[tier], candidate_covariates=cov_tuple, verify_parity=verify
            )
            if verify and prow:
                parity_verified[tier] = True
            if not prow:
                p1_rows.append({"group": smtsd, "tier": tier, "status": f"untestable_{reason}", "n_donors": n_mature_complete})
            else:
                prow["group"] = smtsd
                prow["tier"] = tier
                prow["status"] = "tested"
                prow["predictor"] = "age_bin_order"
                prow["cohort_class"] = cohort_class
                p1_rows.append(prow)
        append_runtime(
            f"track1_{smtsd}",
            s0,
            t0,
            "ok",
            f"n={n_mature_complete} excl={n_excluded} genes={g_sum.get('n_genes_tested', 0)} within_q05={g_sum.get('n_sig_within_q05', 0)} status={g_sum['status']}",
        )
        print(
            f"[Track1] {smtsd}: n={n_mature_complete}/{orig_n_complete}, genes={g_sum.get('n_genes_tested', 0)}, "
            f"within_q05={g_sum.get('n_sig_within_q05', 0)}, top={g_sum.get('top_gene', '')} p={g_sum.get('min_pvalue', float('nan')):.3g}",
            flush=True,
        )
        gc.collect()

    genes_all = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    p1_all = pd.DataFrame(p1_rows)

    # Family BH: primary 47 and cellline 2, separately, per p-value type.
    def _family_bh(df: pd.DataFrame, mask: pd.Series, pcol: str, qcol: str) -> None:
        q_exp, _, _ = bh_checked(df.loc[mask, pcol])
        df.loc[mask, qcol] = q_exp

    genes_all["bh_family_name"] = ""
    genes_all["bh_family_size"] = 0
    genes_all["q_bh_family"] = np.nan
    for cls, fam_name in (("primary_tissue", "track1_47_primary_tissues_mature"), ("isolated_cell_line", "track1_2_isolated_cell_lines_mature")):
        m = genes_all["cohort_class"] == cls
        genes_all.loc[m, "bh_family_name"] = fam_name
        genes_all.loc[m, "bh_family_size"] = int(m.sum())
        _family_bh(genes_all, m, "pvalue", "q_bh_family")

    p1_tested = p1_all[p1_all["status"] == "tested"].copy()
    p1_tested["bh_family_name"] = ""
    p1_tested["family_c_size"] = 0
    p1_tested["q_bh_composite"] = np.nan
    p1_tested["q_bh_empirical"] = np.nan
    for tier in ("P1m", "P0m"):
        for cls, fam_name in (("primary_tissue", f"track1_47_primary_tissues_mature_{tier}"), ("isolated_cell_line", f"track1_2_isolated_cell_lines_mature_{tier}")):
            m = (p1_tested["tier"] == tier) & (p1_tested["cohort_class"] == cls)
            p1_tested.loc[m, "bh_family_name"] = fam_name
            p1_tested.loc[m, "family_c_size"] = int(m.sum())
            _family_bh(p1_tested, m, "composite_pvalue", "q_bh_composite")
            _family_bh(p1_tested, m, "empirical_null_pvalue", "q_bh_empirical")

    sum_df = pd.DataFrame(gene_summaries)
    sig_fam = genes_all.groupby("group")["q_bh_family"].apply(lambda s: int((s < ALPHA).sum())).to_dict()
    minq_fam = genes_all.groupby("group")["q_bh_family"].min().to_dict()
    minp = genes_all.groupby("group")["pvalue"].min().to_dict()
    sum_df["n_sig_family_q05"] = sum_df["group"].map(sig_fam).fillna(0).astype(int)
    sum_df["min_q_bh_family"] = sum_df["group"].map(minq_fam)
    sum_df["min_pvalue_all"] = sum_df["group"].map(minp)
    # Merge the family q columns back into the full P1 table (tested and
    # untested rows) so the detail CSV carries them.
    p1_all = p1_all.merge(
        p1_tested[["group", "tier", "q_bh_composite", "q_bh_empirical", "family_c_size", "bh_family_name"]],
        on=["group", "tier"],
        how="left",
    )
    return {
        "genes_all": genes_all,
        "p1_all": p1_all,
        "p1_tested": p1_tested,
        "summary": sum_df,
        "audit": pd.DataFrame(audit_rows),
    }


def run_track2(
    mm: np.memmap,
    genes_df: pd.DataFrame,
    sample_id_to_col: dict[str, int],
    merged_meta: pd.DataFrame,
    programs: dict[str, dict[str, list[str]]],
    orig: dict[str, pd.DataFrame],
    n_cpus_track2: int = N_CPUS_TRACK2,
) -> dict[str, Any]:
    """Track 2 PyDESeq2 NB-GLM primary arm on mature bins for the 6 tissues."""
    tech_mask = extended_technical_mask(genes_df["Description"].astype(str).to_numpy())
    nontech_idx = np.flatnonzero(~tech_mask)
    nontech_ens = genes_df["ensembl_id"].astype(str).to_numpy()[nontech_idx]
    ens_to_sym = dict(zip(genes_df["ensembl_id"].astype(str), genes_df["Description"].astype(str)))

    def _prepare_counts(sids: list[str]) -> pd.DataFrame:
        cols = [sample_id_to_col[sid] for sid in sids]
        # Column-first fancy indexing: same values as row-first indexing,
        # but it materialises only this tissue's columns.
        mat = np.asarray(mm[:, cols][nontech_idx, :], dtype=np.int32)
        min_samp = max(5, int(np.ceil(0.20 * len(sids))))
        keep_g = (mat >= 10).sum(axis=1) >= min_samp
        return pd.DataFrame(mat[keep_g, :], index=pd.Index(nontech_ens[keep_g], dtype=str), columns=sids)

    gene_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    p1_rows: list[dict[str, Any]] = []
    parity_checked = False

    for tissue in TRACK2_PRIMARY_TISSUES:
        t0 = time.monotonic()
        s0 = utc_now()
        t_meta_full = merged_meta.loc[merged_meta["SMTSD"] == tissue].copy()
        n_20_29 = int((t_meta_full["age_bin"] == "20-29").sum())
        t_meta = t_meta_full.loc[t_meta_full["age_bin"].isin(MATURE_BINS)].copy()
        orig_n_complete = int(orig["t2_sum"].loc[tissue, "n_donors"])
        num_cands = ["sex_male", "SMRIN", "SMTSISCH_hr", "SMEXNCRT"]
        cat_cands = ["SMCENTER", "SMNABTCHT"]
        cov_d = build_clean_covariate_design(t_meta, num_cands, cat_cands, age_col="age_bin_order")
        c_mat = _prepare_counts(list(cov_d.meta.index))
        n_mature = int(cov_d.n_samples_complete)
        print(f"[Track2] Fitting PyDESeq2 {tissue}: n={n_mature} (orig {orig_n_complete}), genes={len(c_mat)}", flush=True)
        res_p = fit_pydeseq2_continuous(c_mat, cov_d.meta, cov_d, age_col="age_bin_order", n_cpus=n_cpus_track2)
        res_p.insert(0, "tissue", tissue)
        res_p["symbol"] = res_p["ensembl_id"].map(ens_to_sym)
        valid = res_p["pvalue"].notna().to_numpy()
        q_in = np.full(len(res_p), np.nan)
        if valid.any():
            q_in[valid], _, _ = bh_checked(res_p.loc[valid, "pvalue"])
        res_p["q_bh_within_tissue"] = q_in
        gene_frames.append(res_p)
        min_idx = int(res_p["pvalue"].idxmin()) if valid.any() else -1
        summary_rows.append(
            {
                "tissue": tissue,
                "arm": "primary_continuous",
                "n_donors": n_mature,
                "n_excluded_vs_original": orig_n_complete - n_mature,
                "n_20_29_excluded": n_20_29,
                "residual_df_approx": int(res_p["residual_df_approx"].iloc[0]),
                "genes_tested": int(valid.sum()),
                "sig_q05_within": int((res_p["q_bh_within_tissue"] < ALPHA).sum()),
                "n_up_q05_within": int(((res_p["q_bh_within_tissue"] < ALPHA) & (res_p["log2FoldChange"] > 0)).sum()),
                "n_down_q05_within": int(((res_p["q_bh_within_tissue"] < ALPHA) & (res_p["log2FoldChange"] < 0)).sum()),
                "min_pvalue": float(res_p["pvalue"].min()),
                "top_gene": str(res_p.loc[min_idx, "symbol"]) if min_idx >= 0 else "",
                "top_gene_lfc": float(res_p.loc[min_idx, "log2FoldChange"]) if min_idx >= 0 else np.nan,
                "top_gene_se": float(res_p.loc[min_idx, "lfcSE"]) if min_idx >= 0 else np.nan,
                "covariates": ";".join(cov_d.all_cov_cols),
                "dropped_collinear": ";".join(cov_d.dropped_collinear) if cov_d.dropped_collinear else "none",
            }
        )
        for tier in ("P1m", "P0m"):
            prow = nb_p1_mature(res_p, ens_to_sym, programs[tier], iterations=P1_ITERATIONS)
            prow["tissue"] = tissue
            prow["tier"] = tier
            prow["n_donors"] = n_mature
            prow["residual_df_approx"] = int(res_p["residual_df_approx"].iloc[0])
            p1_rows.append(prow)
            if not parity_checked and tier == "P1m":
                ref = compute_p1_on_nb_results(res_p, ens_to_sym, programs[tier], iterations=P1_ITERATIONS)
                for key in ("nb_p1_effect_lfc", "nb_p1_composite_z", "nb_p1_composite_p", "nb_p1_empirical_p"):
                    assert abs(float(prow[key]) - float(ref[key])) < 1e-12, f"NB P1 parity failed on {key}: {prow[key]} vs {ref[key]}"
                parity_checked = True
                print(f"[Track2] NB P1 parity vs original compute_p1_on_nb_results verified on {tissue} P1m", flush=True)
        append_runtime(
            f"track2_{tissue}", s0, t0, "ok",
            f"n={n_mature} excl={orig_n_complete - n_mature} genes={int(valid.sum())} within_q05={int((res_p['q_bh_within_tissue'] < ALPHA).sum())}",
        )
        print(
            f"[Track2] {tissue}: n={n_mature}, genes={int(valid.sum())}, within_q05={int((res_p['q_bh_within_tissue'] < ALPHA).sum())}, "
            f"P1m eff={p1_rows[-2]['nb_p1_effect_lfc']:+.4f} emp_p={p1_rows[-2]['nb_p1_empirical_p']:.4g}, "
            f"P0m eff={p1_rows[-1]['nb_p1_effect_lfc']:+.4f} emp_p={p1_rows[-1]['nb_p1_empirical_p']:.4g}",
            flush=True,
        )
        del res_p, c_mat
        gc.collect()

    all_genes = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    valid_g = all_genes["pvalue"].notna().to_numpy()
    all_genes["family_b_track2_size"] = int(valid_g.sum())
    all_genes["q_bh_track2_global"] = np.nan
    if valid_g.any():
        q_glob, _, _ = bh_checked(all_genes.loc[valid_g, "pvalue"])
        all_genes.loc[valid_g, "q_bh_track2_global"] = q_glob
    sig_glob = all_genes.groupby("tissue")["q_bh_track2_global"].apply(lambda s: int((s < ALPHA).sum())).to_dict()
    minq_glob = all_genes.groupby("tissue")["q_bh_track2_global"].min().to_dict()
    sum_df = pd.DataFrame(summary_rows)
    sum_df["sig_q05_global_6tissues"] = sum_df["tissue"].map(sig_glob).fillna(0).astype(int)
    sum_df["min_q_bh_track2_global"] = sum_df["tissue"].map(minq_glob)

    p1_df = pd.DataFrame(p1_rows)
    p1_df["family_size"] = int(len(p1_df) // 2)
    p1_df["q_bh_composite"] = np.nan
    p1_df["q_bh_empirical"] = np.nan
    for tier in ("P1m", "P0m"):
        m = p1_df["tier"] == tier
        q_c, _, _ = bh_checked(p1_df.loc[m, "nb_p1_composite_p"])
        q_e, _, _ = bh_checked(p1_df.loc[m, "nb_p1_empirical_p"])
        p1_df.loc[m, "q_bh_composite"] = q_c
        p1_df.loc[m, "q_bh_empirical"] = q_e
    return {"genes_all": all_genes, "summary": sum_df, "p1": p1_df}

def assemble_main_rows(
    t1: dict[str, Any],
    t2: dict[str, Any],
    orig: dict[str, pd.DataFrame],
    groups_49: list[str],
) -> pd.DataFrame:
    """Assemble the Batch A schema main table from mature and original values."""
    out_rows: list[dict[str, Any]] = []
    t1_sum = t1["summary"].set_index("group")
    t1_p1 = t1["p1_tested"].set_index(["group", "tier"])
    t1_p1_all = t1["p1_all"].set_index(["group", "tier"])
    t2_sum = t2["summary"].set_index("tissue")
    t2_p1 = t2["p1"].set_index(["tissue", "tier"])

    def _extra(**kw: Any) -> dict[str, Any]:
        base: dict[str, Any] = {
            "arm": "",
            "residual_df": np.nan,
            "original_effect": np.nan,
            "mature_effect": np.nan,
            "mature_se": np.nan,
            "null_sd": np.nan,
            "mde80": np.nan,
            "mde80_method": "",
            "covariates": "",
            "omitted_covariates": "",
            "n_genes_tested": np.nan,
            "n_sig_q05": np.nan,
            "original_empirical_p": np.nan,
            "original_empirical_q": np.nan,
            "orig_transport_class": "",
            "mature_transport_class": "",
            "n_up_detected": np.nan,
            "n_down_detected": np.nan,
        }
        base.update(kw)
        return base

    for grp in groups_49:
        o_gsum = orig["t1_genes_sum"].loc[grp]
        o_p1 = orig["t1_p1"].loc[grp]
        o_audit = orig["t1_audit"].loc[grp]
        cohort_class = str(o_audit["cohort_class"])
        orig_n_complete = int(o_gsum["n_donors"])
        orig_sig_genes = bool(int(o_gsum["sig_q_family_0_05"]) > 0)
        orig_minq_genes = float(orig["t1_min_family_q"][grp])
        orig_cls = transport_class(float(o_p1["q_bh_composite"]), float(o_p1["q_bh_empirical"]), True)
        orig_both = orig_cls == "BOTH_PASS"

        s = t1_sum.loc[grp]
        n_mature = int(s["n_donors"])
        n_ex = orig_n_complete - n_mature
        if str(s["status"]) != "tested":
            out_rows.append(
                make_row(
                    COHORT, grp, "genes_track1", n_mature, n_ex,
                    float(o_gsum["min_pvalue"]), orig_minq_genes,
                    None, None, None, None, "UNTESTABLE",
                    f"reason={s['status']};cohort_class={cohort_class}",
                    **_extra(arm="track1_ols"),
                )
            )
            for tier in ("P1m", "P0m"):
                out_rows.append(
                    make_row(
                        COHORT, grp, f"transport_{tier}_track1", n_mature, n_ex,
                        float(o_p1["composite_pvalue"]), float(o_p1["q_bh_composite"]),
                        None, None, None, None, "UNTESTABLE",
                        f"reason={s['status']};orig_class={orig_cls}",
                        **_extra(
                            arm="track1_ols", original_effect=float(o_p1["composite_effect"]),
                            original_empirical_p=float(o_p1["empirical_null_pvalue"]),
                            original_empirical_q=float(o_p1["q_bh_empirical"]),
                            orig_transport_class=orig_cls, mature_transport_class="UNTESTABLE",
                        ),
                    )
                )
            continue

        mat_sig_genes = bool(int(s["n_sig_family_q05"]) > 0)
        verd_g, _ = verdict_standard(orig_sig_genes, mat_sig_genes)
        out_rows.append(
            make_row(
                COHORT, grp, "genes_track1", n_mature, n_ex,
                float(o_gsum["min_pvalue"]), orig_minq_genes,
                float(s["min_pvalue_all"]), float(s["min_q_bh_family"]),
                None, None, verd_g,
                f"n_genes={int(s['n_genes_tested'])};orig_sig={int(o_gsum['sig_q_family_0_05'])};"
                f"mature_sig={int(s['n_sig_family_q05'])};top={s['top_gene']};"
                f"within_q05={int(s['n_sig_within_q05'])};cohort_class={cohort_class}",
                **_extra(
                    arm="track1_ols", residual_df=int(s["residual_df"]),
                    mature_effect=float(s["top_gene_effect"]), mature_se=float(s["median_se"]),
                    mde80=float(s["median_mde80"]), mde80_method="median_t_quantile_se",
                    covariates=str(s["covariates"]), omitted_covariates=str(s["omitted_covariates"]),
                    n_genes_tested=int(s["n_genes_tested"]), n_sig_q05=int(s["n_sig_family_q05"]),
                ),
            )
        )

        for tier in ("P1m", "P0m"):
            prow_all = t1_p1_all.loc[(grp, tier)]
            if str(prow_all["status"]) != "tested":
                out_rows.append(
                    make_row(
                        COHORT, grp, f"transport_{tier}_track1", n_mature, n_ex,
                        float(o_p1["composite_pvalue"]), float(o_p1["q_bh_composite"]),
                        None, None, None, None, "UNTESTABLE",
                        f"reason={prow_all['status']};orig_class={orig_cls}",
                        **_extra(
                            arm="track1_ols", original_effect=float(o_p1["composite_effect"]),
                            original_empirical_p=float(o_p1["empirical_null_pvalue"]),
                            original_empirical_q=float(o_p1["q_bh_empirical"]),
                            orig_transport_class=orig_cls, mature_transport_class="UNTESTABLE",
                        ),
                    )
                )
                continue
            r = t1_p1.loc[(grp, tier)]
            mat_cls = transport_class(float(r["q_bh_composite"]), float(r["q_bh_empirical"]), True)
            mat_both = mat_cls == "BOTH_PASS"
            verd, vnote = verdict_standard(orig_both, mat_both, float(o_p1["composite_effect"]), float(r["composite_effect"]))
            mde_null = mde80_from_null_sd(float(r["null_sd_effect"]))
            mde_se = mde80_from_se(float(r["composite_se"]), int(r["residual_df"]))
            note_parts = [
                f"orig_class={orig_cls}",
                f"mature_class={mat_cls}",
                f"up_eff={float(r['up_effect']):+.4f}(p={float(r['up_pvalue']):.3g})",
                f"dn_eff={float(r['down_effect']):+.4f}(p={float(r['down_pvalue']):.3g})",
                f"mde80_se={mde_se:.4f}",
                f"draws={P1_ITERATIONS}",
                f"seed={ORIG_NULL_SEED}",
                f"min_decile_margin={int(r['min_decile_margin'])}",
            ]
            if vnote:
                note_parts.append(vnote)
            out_rows.append(
                make_row(
                    COHORT, grp, f"transport_{tier}_track1", n_mature, n_ex,
                    float(o_p1["composite_pvalue"]), float(o_p1["q_bh_composite"]),
                    float(r["composite_pvalue"]), float(r["q_bh_composite"]),
                    float(r["empirical_null_pvalue"]), float(r["q_bh_empirical"]),
                    verd, ";".join(note_parts),
                    **_extra(
                        arm="track1_ols", residual_df=int(r["residual_df"]),
                        original_effect=float(o_p1["composite_effect"]),
                        mature_effect=float(r["composite_effect"]), mature_se=float(r["composite_se"]),
                        null_sd=float(r["null_sd_effect"]), mde80=mde_null, mde80_method="2.4865_times_null_sd",
                        covariates=str(r["covariates"]), omitted_covariates=str(r["omitted_covariates"]),
                        original_empirical_p=float(o_p1["empirical_null_pvalue"]),
                        original_empirical_q=float(o_p1["q_bh_empirical"]),
                        orig_transport_class=orig_cls, mature_transport_class=mat_cls,
                        n_up_detected=int(r["n_up_detected"]), n_down_detected=int(r["n_down_detected"]),
                    ),
                )
            )

    for tissue in TRACK2_PRIMARY_TISSUES:
        o_sum = orig["t2_sum"].loc[tissue]
        o_p1 = orig["t2_p1"].loc[tissue]
        orig_n = int(o_sum["n_donors"])
        orig_sig_genes = bool(int(o_sum["sig_q05_global_6tissues"]) > 0)
        orig_minq_genes = float(orig["t2_min_global_q"][tissue])
        orig_cls = transport_class(float(o_p1["q_bh_composite"]), float(o_p1["q_bh_empirical"]), True)
        orig_both = orig_cls == "BOTH_PASS"
        s = t2_sum.loc[tissue]
        n_mature = int(s["n_donors"])
        n_ex = orig_n - n_mature
        mat_sig_genes = bool(int(s["sig_q05_global_6tissues"]) > 0)
        verd_g, _ = verdict_standard(orig_sig_genes, mat_sig_genes)
        out_rows.append(
            make_row(
                COHORT, tissue, "genes_track2", n_mature, n_ex,
                float(o_sum["min_pvalue"]), orig_minq_genes,
                float(s["min_pvalue"]), float(s["min_q_bh_track2_global"]),
                None, None, verd_g,
                f"n_genes={int(s['genes_tested'])};orig_sig_global={int(o_sum['sig_q05_global_6tissues'])};"
                f"mature_sig_global={int(s['sig_q05_global_6tissues'])};top={s['top_gene']};"
                f"within_q05={int(s['sig_q05_within'])}",
                **_extra(
                    arm="track2_nbglm", residual_df=int(s["residual_df_approx"]),
                    mature_effect=float(s["top_gene_lfc"]), mature_se=float(s["top_gene_se"]),
                    mde80=np.nan, mde80_method="not_defined_nb_glm",
                    covariates=str(s["covariates"]), omitted_covariates=str(s["dropped_collinear"]),
                    n_genes_tested=int(s["genes_tested"]), n_sig_q05=int(s["sig_q05_global_6tissues"]),
                ),
            )
        )
        for tier in ("P1m", "P0m"):
            r = t2_p1.loc[(tissue, tier)]
            mat_cls = transport_class(float(r["q_bh_composite"]), float(r["q_bh_empirical"]), True)
            mat_both = mat_cls == "BOTH_PASS"
            verd, vnote = verdict_standard(orig_both, mat_both, float(o_p1["nb_p1_effect_lfc"]), float(r["nb_p1_effect_lfc"]))
            z_val = float(r["nb_p1_composite_z"])
            eff = float(r["nb_p1_effect_lfc"])
            se_nb = abs(eff / z_val) if z_val != 0 else np.nan
            mde_null = mde80_from_null_sd(float(r["nb_p1_null_sd"]))
            note_parts = [
                f"orig_class={orig_cls}",
                f"mature_class={mat_cls}",
                f"comp_z={z_val:+.3f}",
                f"draws={P1_ITERATIONS}",
                f"seed={ORIG_NULL_SEED}",
            ]
            if vnote:
                note_parts.append(vnote)
            out_rows.append(
                make_row(
                    COHORT, tissue, f"transport_{tier}_track2", n_mature, n_ex,
                    float(o_p1["nb_p1_composite_p"]), float(o_p1["q_bh_composite"]),
                    float(r["nb_p1_composite_p"]), float(r["q_bh_composite"]),
                    float(r["nb_p1_empirical_p"]), float(r["q_bh_empirical"]),
                    verd, ";".join(note_parts),
                    **_extra(
                        arm="track2_nbglm", residual_df=int(r["residual_df_approx"]),
                        original_effect=float(o_p1["nb_p1_effect_lfc"]),
                        mature_effect=eff, mature_se=se_nb,
                        null_sd=float(r["nb_p1_null_sd"]), mde80=mde_null, mde80_method="2.4865_times_null_sd",
                        original_empirical_p=float(o_p1["nb_p1_empirical_p"]),
                        original_empirical_q=float(o_p1["q_bh_empirical"]),
                        orig_transport_class=orig_cls, mature_transport_class=mat_cls,
                        n_up_detected=int(r["nb_p1_n_up"]), n_down_detected=int(r["nb_p1_n_down"]),
                    ),
                )
            )

    return order_cohort_df(out_rows)


def build_bh_check() -> pd.DataFrame:
    """Independent BH recomputation from the written detail CSVs.

    For each family: recompute BH over the FULL family with statsmodels
    multipletests and compare against the q values stored in the files
    (max diff must be <= 1e-12). Recorded rows are the per-tissue top gene
    for gene families (full families are 748,735 and 140,383 tests) and every
    tested unit for transport families, plus one ALL row per family carrying
    the full-family max diff.
    """
    rows: list[dict[str, Any]] = []

    def _full_max_diff(pvals: np.ndarray, q_stored: np.ndarray) -> tuple[float, np.ndarray]:
        valid = np.isfinite(pvals)
        q_sm = np.full(len(pvals), np.nan)
        if valid.any():
            q_sm[valid] = multipletests(pvals[valid], method="fdr_bh")[1]
        diffs = np.abs(q_stored - q_sm)
        max_diff = float(np.nanmax(diffs)) if np.isfinite(diffs).any() else 0.0
        if max_diff > 1e-12:
            raise AssertionError(f"BH check failed: max diff {max_diff:.3e} > 1e-12")
        return max_diff, q_sm

    def _record(family: str, family_size: int, units: list[str], pvals: np.ndarray, q_stored: np.ndarray, q_sm: np.ndarray, max_diff: float) -> None:
        for unit, p, qs, qm in zip(units, pvals, q_stored, q_sm):
            rows.append(
                {
                    "cohort": COHORT,
                    "family": family,
                    "unit": unit,
                    "family_size": int(family_size),
                    "raw_p": float(p),
                    "q_bh_explicit": float(qs),
                    "q_bh_statsmodels": float(qm),
                    "abs_diff": float(abs(qs - qm)),
                    "family_max_abs_diff": max_diff,
                }
            )
        rows.append(
            {
                "cohort": COHORT,
                "family": family,
                "unit": "ALL",
                "family_size": int(family_size),
                "raw_p": float("nan"),
                "q_bh_explicit": float("nan"),
                "q_bh_statsmodels": float("nan"),
                "abs_diff": float("nan"),
                "family_max_abs_diff": max_diff,
            }
        )

    t1_genes = pd.read_csv(T1_GENES_CSV)
    for cls, fam in (
        ("primary_tissue", "track1_47_primary_tissues_mature_genes"),
        ("isolated_cell_line", "track1_2_isolated_cell_lines_mature_genes"),
    ):
        sub = t1_genes[t1_genes["cohort_class"] == cls].reset_index(drop=True)
        max_diff, q_sm = _full_max_diff(sub["pvalue"].to_numpy(), sub["q_bh_family"].to_numpy())
        top_pos = sub.groupby("group")["pvalue"].idxmin()
        top = sub.loc[top_pos]
        _record(
            fam, len(sub),
            [f"{g}::{gn}" for g, gn in zip(top["group"], top["gene"])],
            top["pvalue"].to_numpy(), top["q_bh_family"].to_numpy(), q_sm[top_pos.to_numpy()],
            max_diff,
        )

    t1_p1 = pd.read_csv(T1_P1_CSV)
    t1_p1 = t1_p1[t1_p1["status"] == "tested"]
    for tier in ("P1m", "P0m"):
        for cls, kind, pcol, qcol in (
            ("primary_tissue", "composite", "composite_pvalue", "q_bh_composite"),
            ("primary_tissue", "empirical", "empirical_null_pvalue", "q_bh_empirical"),
            ("isolated_cell_line", "composite", "composite_pvalue", "q_bh_composite"),
            ("isolated_cell_line", "empirical", "empirical_null_pvalue", "q_bh_empirical"),
        ):
            sub = t1_p1[(t1_p1["tier"] == tier) & (t1_p1["cohort_class"] == cls)].reset_index(drop=True)
            fam = f"track1_mature_{tier}_{kind}_{cls}"
            max_diff, q_sm = _full_max_diff(sub[pcol].to_numpy(), sub[qcol].to_numpy())
            _record(fam, len(sub), sub["group"].astype(str).tolist(), sub[pcol].to_numpy(), sub[qcol].to_numpy(), q_sm, max_diff)

    t2_genes = pd.read_csv(T2_GENES_CSV)
    sub2 = t2_genes.reset_index(drop=True)
    max_diff2, q_sm2 = _full_max_diff(sub2["pvalue"].to_numpy(), sub2["q_bh_track2_global"].to_numpy())
    top_pos2 = sub2.groupby("tissue")["pvalue"].idxmin()
    top2 = sub2.loc[top_pos2]
    _record(
        "track2_6_tissues_mature_genes", len(sub2),
        [f"{t}::{s}" for t, s in zip(top2["tissue"], top2["symbol"])],
        top2["pvalue"].to_numpy(), top2["q_bh_track2_global"].to_numpy(), q_sm2[top_pos2.to_numpy()],
        max_diff2,
    )

    t2_p1 = pd.read_csv(T2_P1_CSV)
    for tier in ("P1m", "P0m"):
        for kind, pcol, qcol in (
            ("composite", "nb_p1_composite_p", "q_bh_composite"),
            ("empirical", "nb_p1_empirical_p", "q_bh_empirical"),
        ):
            sub = t2_p1[t2_p1["tier"] == tier].reset_index(drop=True)
            fam = f"track2_6_tissues_mature_{tier}_{kind}"
            max_diff, q_sm = _full_max_diff(sub[pcol].to_numpy(), sub[qcol].to_numpy())
            _record(fam, len(sub), sub["tissue"].astype(str).tolist(), sub[pcol].to_numpy(), sub[qcol].to_numpy(), q_sm, max_diff)

    return pd.DataFrame(rows)


def build_report(main_df: pd.DataFrame, file_audit: dict[str, Any], gate: dict[str, Any], t1_audit: pd.DataFrame) -> str:
    """Build the report from the written main CSV read back from disk."""
    m = pd.read_csv(MAIN_CSV)
    # pandas read_csv treats the literal verdict string "NULL" as a missing
    # value by default; restore it so verdict counts are complete.
    m["verdict"] = m["verdict"].fillna("NULL")
    assert len(m) == len(main_df) and set(m.columns) == set(main_df.columns)

    def _bothpass(df: pd.DataFrame, fam: str) -> pd.DataFrame:
        return df[(df["endpoint_family"] == fam) & (df["mature_transport_class"] == "BOTH_PASS")]

    def _verdicts(df: pd.DataFrame, fam: str) -> dict[str, int]:
        vc = df[df["endpoint_family"] == fam]["verdict"].value_counts().to_dict()
        return {k: int(v) for k, v in vc.items()}

    t1_genes = m[m["endpoint_family"] == "genes_track1"]
    t2_genes = m[m["endpoint_family"] == "genes_track2"]
    n_t1_tested = int((t1_genes["verdict"] != "UNTESTABLE").sum())
    n_t2_tested = int((t2_genes["verdict"] != "UNTESTABLE").sum())

    lines: list[str] = []
    lines.append("# Item 5 Batch C: GTEx v8 mature (30-plus bins) transport rerun")
    lines.append("")
    lines.append("## 1. Scope, input verification, and frozen rules")
    lines.append("")
    lines.append(
        "GTEx v8 open-access gene-read counts only, no dbGaP application, no exact "
        "ages. File verified before any analysis: "
        f"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz downloaded from "
        f"{GCT_URL} ; size {file_audit['size_bytes']} bytes (expected {EXPECTED_GCT_BYTES}); "
        f"MD5 base64 {file_audit['md5_base64']} (expected {EXPECTED_GCT_MD5_B64}); "
        f"md5_verified = {str(file_audit['md5_verified']).lower()}. "
        f"SHA-256 {file_audit['sha256']}. Count gate on the full uncompressed stream: "
        f"{gate['n_genes']} gene rows x {gate['n_samples']} RNASEQ sample columns, "
        f"{gate['nonfinite']} nonfinite, {gate['negative']} negative, "
        f"{gate['noninteger_lines']} noninteger lines, {gate['total_nonzero_entries']} nonzero entries."
    )
    lines.append("")
    lines.append(
        "Age restriction: GTEx reports only ordinal decade bins 20-29, 30-39, 40-49, "
        "50-59, 60-69, 70-79. Subject counts by bin verified from "
        "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt: 20-29: 84, 30-39: 78, "
        "40-49: 153, 50-59: 315, 60-69: 317, 70-79: 33, total 980 subjects. The rerun "
        "uses bins 30-39 through 70-79 only, with age_bin_order kept at its original "
        "values 2.0 to 6.0 so each slope stays per decade bin. No exact age was "
        "imputed. Donors in the 20-29 bin are excluded before count aggregation, "
        "detection filtering, and model fitting in both tracks."
    )
    lines.append("")
    lines.append(
        "Mature signatures from Item 1 (commit d3dc14c), loaded and asserted by "
        "load_verified_mature_programs: P1m 1,914 mouse genes (1,119 up, 795 down), "
        "1,724 human-mapped (1,013 up, 711 down); P0m 2,598 mouse genes (1,474 up, "
        "1,124 down), 2,049 human-mapped (1,157 up, 892 down). Seed 20261105. "
        "2000-draw expression-decile-matched null as primary. Donor is the "
        "inferential unit, never the cell. Bulk RNA-seq only: gene-level age tests "
        "and signature tests only, no fraction/CLR models and no cross-tissue "
        "concordance, exactly as the original GTEx report states."
    )
    lines.append("")
    lines.append("Multiple-testing families (declared before testing):")
    lines.append("1. T1B-primary: all mature gene OLS p-values across the 47 primary tissues.")
    lines.append("2. T1B-cellline: all mature gene OLS p-values across the 2 isolated cell lines.")
    lines.append("3. T1C-primary-P1m and T1C-primary-P0m: 47 composite p-values and 47 empirical p-values each.")
    lines.append("4. T1C-cellline-P1m and T1C-cellline-P0m: 2 composite and 2 empirical p-values each.")
    lines.append("5. T2B-global: all mature NB-GLM gene p-values across the 6 focal tissues.")
    lines.append("6. T2C-P1m and T2C-P0m: 6 composite Wald p-values and 6 empirical p-values each.")
    lines.append("BH by bh_checked (explicit BH cross-checked against statsmodels, assert diff at most 1e-12); the independent recomputation is in results/followup-mature-organ-transports-c-gtex-bh-check.csv.")
    lines.append("")
    lines.append("Verdict rule: verdict_standard from the Batch A script, applied to the original versus mature both-pass (transport) or any-family-hit (genes) status. SURVIVED means the original signal is still significant at mature q below 0.05; DIED means the original signal is gone; NEW means a mature-only signal; NULL means neither side is significant; UNTESTABLE means the row could not be tested at 30-plus. Note on the task wording: the Batch A partial report claimed 48 of 68 Gao rows have a blank verdict, but the Batch A file results/followup-mature-organ-transports-a-gao.csv as read here contains zero blank verdicts (48 NULL, 20 DIED), so no blank verdicts were invented here and the Batch A verdict_standard rule was used unchanged.")
    lines.append("")

    # Track 1 transport headline
    lines.append("## 2. Track 1 transport (OLS, P1m and P0m)")
    lines.append("")
    for tier in ("P1m", "P0m"):
        fam = f"transport_{tier}_track1"
        sub = m[m["endpoint_family"] == fam]
        prim = sub[sub["cell_type"].isin(t1_audit[t1_audit["cohort_class"] == "primary_tissue"]["SMTSD"])]
        cell = sub[sub["cell_type"].isin(t1_audit[t1_audit["cohort_class"] == "isolated_cell_line"]["SMTSD"])]
        bp_prim = _bothpass(prim, fam)
        bp_cell = _bothpass(cell, fam)
        lines.append(
            f"Under {tier}, {len(bp_prim)} of {len(prim)} primary tissues pass both "
            f"parametric and empirical q below 0.05, and {len(bp_cell)} of {len(cell)} "
            f"isolated cell lines pass both."
        )
        for _, r in bp_prim.iterrows():
            lines.append(
                f"  Both-pass {tier} primary: {r['cell_type']}: effect {float(r['mature_effect']):+.4f} per decade bin, "
                f"q_comp {float(r['mature_q']):.4g}, q_emp {float(r['mature_empirical_q']):.4g}, verdict {r['verdict']} "
                f"(original class {r['orig_transport_class']})."
            )
        for _, r in bp_cell.iterrows():
            lines.append(
                f"  Both-pass {tier} cell line: {r['cell_type']}: effect {float(r['mature_effect']):+.4f} per decade bin, "
                f"q_comp {float(r['mature_q']):.4g}, q_emp {float(r['mature_empirical_q']):.4g}, verdict {r['verdict']} "
                f"(original class {r['orig_transport_class']})."
            )
        lines.append(f"  Verdicts {tier} Track 1 transport (49 rows): {_verdicts(m, fam)}.")
        lines.append("")
    lines.append(
        "Original Track 1 dual-pass tissues (3 positive, 3 negative) and their mature verdicts:"
    )
    orig_dual = [
        ("Nerve - Tibial", "positive"), ("Adrenal Gland", "positive"), ("Thyroid", "positive"),
        ("Adipose - Visceral (Omentum)", "negative"), ("Colon - Sigmoid", "negative"), ("Minor Salivary Gland", "negative"),
    ]
    for tissue, direction in orig_dual:
        for tier in ("P1m", "P0m"):
            r = m[(m["cell_type"] == tissue) & (m["endpoint_family"] == f"transport_{tier}_track1")].iloc[0]
            lines.append(
                f"  {tissue} ({direction} originally), {tier}: mature effect {float(r['mature_effect']):+.4f}, "
                f"q_comp {float(r['mature_q']):.4g}, q_emp {float(r['mature_empirical_q']):.4g}, "
                f"class {r['mature_transport_class']}, verdict {r['verdict']}."
            )
    lines.append("")

    # Track 1 genes headline
    t1_tests = int(t1_genes["n_genes_tested"].sum())
    t1_sig = int(t1_genes["n_sig_q05"].sum())
    lines.append("## 3. Track 1 gene-level OLS")
    lines.append("")
    lines.append(
        f"Across {n_t1_tested} tested tissues, {t1_tests} mature gene OLS tests, "
        f"{t1_sig} pass the mature family q below 0.05 "
        f"(primary family and cell-line family separately). Original totals were "
        f"748,735 primary tests with 86,495 family-significant and 26,916 cell-line "
        f"tests with 4,026 family-significant. Verdicts across the 49 genes_track1 rows: "
        f"{_verdicts(m, 'genes_track1')}."
    )
    died = t1_genes[t1_genes["verdict"] == "DIED"]["cell_type"].tolist()
    new = t1_genes[t1_genes["verdict"] == "NEW"]["cell_type"].tolist()
    lines.append(f"  Tissues whose gene-level family signal DIED at 30-plus: {len(died)}: {'; '.join(died) if died else 'none'}.")
    lines.append(f"  Tissues with a NEW gene-level family signal at 30-plus: {len(new)}: {'; '.join(new) if new else 'none'}.")
    lines.append("")

    # Track 2
    lines.append("## 4. Track 2 NB-GLM (6 focal tissues)")
    lines.append("")
    t2_tests = int(t2_genes["n_genes_tested"].sum())
    t2_sig = int(t2_genes["n_sig_q05"].sum())
    lines.append(
        f"Across {n_t2_tested} tested tissues, {t2_tests} mature NB-GLM gene tests, "
        f"{t2_sig} pass the mature global q below 0.05. Original totals were 140,383 "
        f"tests with 23,759 global-significant. Verdicts across the 6 genes_track2 rows: "
        f"{_verdicts(m, 'genes_track2')}."
    )
    for tier in ("P1m", "P0m"):
        fam = f"transport_{tier}_track2"
        sub = m[m["endpoint_family"] == fam]
        bp = _bothpass(sub, fam)
        lines.append(
            f"Under {tier}, {len(bp)} of {len(sub)} focal tissues pass both composite "
            f"and empirical q below 0.05. Verdicts {tier} Track 2 transport (6 rows): {_verdicts(m, fam)}."
        )
        for _, r in bp.iterrows():
            lines.append(
                f"  Both-pass {tier} Track 2: {r['cell_type']}: effect {float(r['mature_effect']):+.4f} LFC per decade bin, "
                f"q_comp {float(r['mature_q']):.4g}, q_emp {float(r['mature_empirical_q']):.4g}, verdict {r['verdict']} "
                f"(original class {r['orig_transport_class']})."
            )
        muscle = sub[sub["cell_type"] == "Muscle - Skeletal"].iloc[0]
        lines.append(
            f"  Muscle - Skeletal, {tier}: mature effect {float(muscle['mature_effect']):+.4f} LFC per decade bin, "
            f"q_comp {float(muscle['mature_q']):.4g}, q_emp {float(muscle['mature_empirical_q']):.4g}, "
            f"class {muscle['mature_transport_class']}, verdict {muscle['verdict']} "
            f"(original P1 was a negative both-pass: effect -0.0560, q_comp 5.85e-12, q_emp 0.00300)."
        )
    lines.append("")

    lines.append("## 5. What was not done")
    lines.append("")
    lines.append(
        "No cross-tissue concordance: the original GTEx analysis ran none for bulk "
        "tissue, and none was run here. No Track 2 sensitivity arms (categorical-age "
        "LRT, oldest-bin exclusion, extended quality/death covariates, signed "
        "ischemic time, read-length restriction): the item asked for the primary "
        "per-gene and transport reruns at 30-plus, and the sensitivities were out of "
        "scope. No exact ages were imputed and no protected GTEx data were accessed. "
        "Tissues whose mature covariate-complete n falls below 50 would be UNTESTABLE; "
        "the per-tissue audit is in the Track 1 genes summary CSV."
    )
    lines.append("")
    lines.append(
        "Reading the empirical null: P1m (1,724 human-mapped genes) and P0m (2,049) "
        "are about 12 to 15 times larger than the original 139-gene P1_Clean, so the "
        "expression-matched random-gene null is much narrower and a small bulk shift "
        "can pass it where the same shift failed against the P1 null. The NEW "
        "both-pass verdicts above are statements about the mature signature under "
        "the frozen protocol, not proof of cell-intrinsic aging: the original "
        "report's central caveat still applies, because bulk expression changes "
        "conflate within-cell aging with shifts in cell-type proportions, and the "
        "Track 2 Whole Blood and Testis both-pass effects are small (+0.0119 to "
        "+0.0122 in blood and -0.0144 to -0.0132 in testis, LFC per decade bin)."
    )
    lines.append("")
    lines.append("## 6. Files")
    lines.append("")
    lines.append("  Main table (Batch A schema, 165 rows): results/followup-mature-organ-transports-c-gtex.csv")
    lines.append("  BH check (independent statsmodels recomputation): results/followup-mature-organ-transports-c-gtex-bh-check.csv")
    lines.append("  Runtime (per-step seconds, peak RSS, donor counts): results/followup-mature-organ-transports-c-gtex-runtime.csv")
    lines.append("  Track 1 gene-level OLS: results/followup-mature-organ-transports-c-gtex-track1-genes.csv")
    lines.append("  Track 1 per-tissue summary: results/followup-mature-organ-transports-c-gtex-track1-genes-summary.csv")
    lines.append("  Track 1 P1m/P0m transport: results/followup-mature-organ-transports-c-gtex-track1-p1.csv")
    lines.append("  Track 2 NB-GLM gene-level: results/followup-mature-organ-transports-c-gtex-track2-nbglm-genes.csv")
    lines.append("  Track 2 per-tissue summary: results/followup-mature-organ-transports-c-gtex-track2-nbglm-summary.csv")
    lines.append("  Track 2 P1m/P0m transport: results/followup-mature-organ-transports-c-gtex-track2-p1.csv")
    lines.append("  Script: scripts/run_followup_mature_organ_transports_c_gtex.py")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gct", default=GCT_PATH, type=Path)
    parser.add_argument("--work-dir", default=WORK_DIR, type=Path)
    parser.add_argument("--n-cpus-track2", default=N_CPUS_TRACK2, type=int)
    args = parser.parse_args()
    t0_all = time.monotonic()
    s0_all = utc_now()

    programs = load_verified_mature_programs()
    orig = load_original_tables()

    # 1. Verify GCT file (size + MD5, same as the original script).
    t0 = time.monotonic()
    s0 = utc_now()
    file_audit = verify_gct_file(args.gct)
    assert int(args.gct.stat().st_size) == EXPECTED_GCT_BYTES, (
        f"GCT size {args.gct.stat().st_size} != expected {EXPECTED_GCT_BYTES}"
    )
    file_audit["path"] = Path(file_audit["path"]).name
    file_audit["url"] = GCT_URL
    file_audit["sha256"] = sha256_file(args.gct)
    print(f"[GTEx C] Verified GCT: {file_audit}", flush=True)
    append_runtime("download_verify", s0, t0, "ok", f"size={file_audit['size_bytes']} md5={file_audit['md5_base64']} sha256={file_audit['sha256'][:16]}...")

    # 2. Metadata + subject bin-count audit.
    t0 = time.monotonic()
    s0 = utc_now()
    merged_meta, subj = load_gtex_sample_subject_metadata()
    bin_counts = subj["AGE"].astype(str).str.strip().value_counts().to_dict()
    assert bin_counts == EXPECTED_SUBJECT_BIN_COUNTS, f"Subject bin counts {bin_counts} != audit {EXPECTED_SUBJECT_BIN_COUNTS}"
    assert len(subj) == 980, f"Subject count {len(subj)} != 980"
    groups_49 = sorted(orig["t1_audit"][orig["t1_audit"]["cohort_class"].isin(["primary_tissue", "isolated_cell_line"])].index.tolist())
    assert len(groups_49) == 49, f"Expected 49 original groups, got {len(groups_49)}"
    print(f"[GTEx C] Metadata OK: 980 subjects, bins {bin_counts}, 49 groups", flush=True)
    append_runtime("metadata_audit", s0, t0, "ok", f"subjects=980 bins={bin_counts} groups49={len(groups_49)}")

    # 3. Stream GCT to memmap (reuses the original function and its gate).
    t0 = time.monotonic()
    s0 = utc_now()
    mmap_path, shape, genes_df, sample_ids, gate = stream_gct_to_memmap(
        args.gct, set(merged_meta.index.astype(str)), args.work_dir
    )
    assert gate == EXPECTED_GATE, f"GCT gate {gate} != expected {EXPECTED_GATE}"
    sample_id_to_col = {sid: idx for idx, sid in enumerate(sample_ids)}
    mm = np.memmap(mmap_path, dtype=np.int32, mode="r", shape=shape)
    print(f"[GTEx C] GCT stream OK: shape={shape}, gate={gate}", flush=True)
    append_runtime("stream_gct_memmap", s0, t0, "ok", f"shape={shape} nonzero={gate['total_nonzero_entries']}")

    # 4. Track 1.
    t0 = time.monotonic()
    s0 = utc_now()
    t1 = run_track1(mm, genes_df, sample_id_to_col, merged_meta, programs, orig, groups_49)
    t1_genes_out = t1["genes_all"][
        ["group", "gene", "n_donors", "predictor", "covariates", "omitted_covariates", "mean_log2cpm",
         "effect", "se", "ci95_low", "ci95_high", "t_stat", "pvalue", "residual_df", "r2",
         "q_bh_within_group", "cohort_class", "bh_family_name", "bh_family_size", "q_bh_family", "mde80"]
    ].sort_values(["group", "pvalue"]).reset_index(drop=True)
    write_clean_csv(t1_genes_out, T1_GENES_CSV)
    write_clean_csv(t1["summary"].sort_values("group").reset_index(drop=True), T1_GENES_SUMMARY_CSV)
    t1_p1_out = t1["p1_all"].sort_values(["tier", "group"]).reset_index(drop=True)
    write_clean_csv(t1_p1_out, T1_P1_CSV)
    append_runtime("track1_total", s0, t0, "ok", f"genes_tests={len(t1_genes_out)} p1_tested={len(t1['p1_tested'])} p1_rows={len(t1_p1_out)}")
    print(f"[GTEx C] Track 1 done: {len(t1_genes_out)} gene tests, {len(t1_p1_out)} tested P1 rows", flush=True)

    # 5. Track 2.
    t0 = time.monotonic()
    s0 = utc_now()
    t2 = run_track2(mm, genes_df, sample_id_to_col, merged_meta, programs, orig, n_cpus_track2=int(args.n_cpus_track2))
    t2_genes_out = t2["genes_all"][
        ["tissue", "ensembl_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
         "residual_df_approx", "symbol", "q_bh_within_tissue", "family_b_track2_size", "q_bh_track2_global"]
    ].sort_values(["tissue", "pvalue"]).reset_index(drop=True)
    write_clean_csv(t2_genes_out, T2_GENES_CSV)
    write_clean_csv(t2["summary"].sort_values("tissue").reset_index(drop=True), T2_SUMMARY_CSV)
    write_clean_csv(t2["p1"].sort_values(["tier", "tissue"]).reset_index(drop=True), T2_P1_CSV)
    append_runtime("track2_total", s0, t0, "ok", f"genes_tests={len(t2_genes_out)} p1_rows={len(t2['p1'])}")
    print(f"[GTEx C] Track 2 done: {len(t2_genes_out)} gene tests, {len(t2['p1'])} P1 rows", flush=True)
    del mm
    gc.collect()

    # 6. Main table (Batch A schema).
    t0 = time.monotonic()
    s0 = utc_now()
    main_df = assemble_main_rows(t1, t2, orig, groups_49)
    assert len(main_df) == 49 * 3 + 6 * 3, f"Main rows {len(main_df)} != 165"
    write_clean_csv(main_df, MAIN_CSV)
    append_runtime("assemble_main", s0, t0, "ok", f"rows={len(main_df)}")

    # 7. Independent BH check from the written files.
    t0 = time.monotonic()
    s0 = utc_now()
    bh_df = build_bh_check()
    write_clean_csv(bh_df, BH_CHECK_CSV)
    append_runtime("bh_check", s0, t0, "ok", f"families={int(bh_df['family'].nunique())} rows={len(bh_df)} max_diff={float(bh_df['family_max_abs_diff'].max()):.3g}")
    print(f"[GTEx C] BH check done: {int(bh_df['family'].nunique())} families, max diff {float(bh_df['family_max_abs_diff'].max()):.3g}", flush=True)

    # 8. Report, built from the main CSV read back from disk.
    t0 = time.monotonic()
    s0 = utc_now()
    report_text = build_report(main_df, file_audit, gate, t1["audit"])
    write_clean_markdown(report_text, REPORT_MD)
    append_runtime("report", s0, t0, "ok", f"report={REPORT_MD.name}")
    append_runtime("cohort_total", s0_all, t0_all, "ok", f"rows={len(main_df)}")
    print(f"[GTEx C DONE] rows={len(main_df)}, elapsed={time.monotonic() - t0_all:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
