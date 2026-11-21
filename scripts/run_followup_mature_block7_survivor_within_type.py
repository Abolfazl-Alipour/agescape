#!/usr/bin/env python3
"""Phase 2 Item 13: Mature-to-old (30-plus) Block 7 within-type death/stress scoring,
removal sensitivity, and narrowing localization.

Adapted from scripts/run_followup_11pt_block7_survivor_within_type.py with frozen protocol,
frozen 271-gene scoring list, and frozen verdict rules preserved verbatim, restricting all
cell eligibility, 90th-percentile thresholds, donor-cell eligibility, regressions, removal
comparisons, and Part B joins strictly to human mature donors aged >= 30 years.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Any

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_between_subject_gao_yazar import _prune_covariate_matrix
from src.followup_organs import (
    counts_to_log2cpm,
    run_family_b_genes_and_concordance,
    write_clean_csv,
    write_clean_markdown,
)

DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
BLOCK1_CACHE = DATA_ROOT / "block1_cache"
CACHE_ROOT = DATA_ROOT / "block7_mature_cache"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"

REPORT_PATH = REPORTS_DIR / "followup-mature-block7-survivor-tests.md"
DEATH_PATH = RESULTS_DIR / "followup-mature-block7-death-stress.csv"
REMOVAL_PATH = RESULTS_DIR / "followup-mature-block7-removal-sensitivity.csv"
NARROWING_PATH = RESULTS_DIR / "followup-mature-block7-narrowing-localization.csv"

ORIG_DEATH_PATH = RESULTS_DIR / "followup-11pt-block7-death-stress.csv"
ORIG_REMOVAL_PATH = RESULTS_DIR / "followup-11pt-block7-removal-sensitivity.csv"
ORIG_NARROWING_PATH = RESULTS_DIR / "followup-11pt-block7-narrowing-localization.csv"

MATURE_COMP_HUMAN_PATH = RESULTS_DIR / "followup-mature-composition-human.csv"
MATURE_FRACTION_PATHS = {
    "gao": RESULTS_DIR / "followup-mature-organ-transports-a-gao.csv",
    "yazar": RESULTS_DIR / "followup-mature-organ-transports-a-yazar.csv",
}
MATURE_SPREAD_PATHS = {
    "gao": RESULTS_DIR / "followup-mature-fulldepth-gao-spread.csv",
    "yazar": RESULTS_DIR / "followup-mature-fulldepth-yazar-spread.csv",
}

MIN_AGE_YEARS = 30.0
MASTER_SEED = 20261117
DATASET_CODE = {"gao": 1, "yazar": 2}
MIN_CELLS = 20
MIN_DONORS = 6
MIN_RESID_DF = {"gao": 5, "yazar": 10}
GENE_MIN_RESID_DF = {"gao": 2, "yazar": 10}
BATCH_REF = {"gao": "ENCODE v4 (Snyder)", "yazar": None}
MDE80_SMALL_SD = 0.50
MIN_AVAILABLE_MEM_GIB = 6.0

GENE_SOURCE = (
    "MSigDB HALLMARK_APOPTOSIS M5902 and HALLMARK_UNFOLDED_PROTEIN_RESPONSE M5922; "
    "PMID 26771021; DOI 10.1016/j.cels.2015.12.004"
)

APOPTOSIS = """ADD1 AIFM3 ANKH ANXA1 APP ATF3 AVPR1A BAX BCAP31 BCL10 BCL2L1 BCL2L10 BCL2L11 BCL2L2 BGN BID BIK BIRC3 BMF BMP2 BNIP3L BRCA1 BTG2 BTG3 CASP1 CASP2 CASP3 CASP4 CASP6 CASP7 CASP8 CASP9 CAV1 CCNA1 CCND1 CCND2 CD14 CD2 CD38 CD44 CD69 CDC25B CDK2 CDKN1A CDKN1B CFLAR CLU CREBBP CTH CTNNB1 CYLD DAP DAP3 DCN DDIT3 DFFA DIABLO DNAJA1 DNAJC3 DNM1L DPYD EBP EGR3 EMP1 ENO2 ERBB2 ERBB3 EREG ETF1 F2 F2R FAS FASLG FDXR FEZ1 GADD45A GADD45B GCH1 GNA15 GPX1 GPX3 GPX4 GSN GSR GSTM1 GUCY2D H1-0 HGF HMGB2 HMOX1 HSPB1 IER3 IFITM3 IFNB1 IFNGR1 IGF2R IGFBP6 IL18 IL1A IL1B IL6 IRF1 ISG20 JUN KRT18 LEF1 LGALS3 LMNA LUM MADD MCL1 MGMT MMP2 NEDD9 NEFH PAK1 PDCD4 PDGFRB PEA15 PLAT PLCB2 PLPPR4 PMAIP1 PPP2R5B PPP3R1 PPT1 PRF1 PSEN1 PSEN2 PTK2 RARA RELA RETSAT RHOB RHOT2 RNASEL ROCK1 SAT1 SATB1 SC5D SLC20A1 SMAD7 SOD1 SOD2 SPTAN1 SQSTM1 TAP1 TGFB2 TGFBR3 TIMP1 TIMP2 TIMP3 TNF TNFRSF12A TNFSF10 TOP2A TSPO TXNIP VDAC2 WEE1 XIAP""".split()
UPR = """ALDH18A1 ARFGAP1 ASNS ATF3 ATF4 ATF6 ATP6V0D1 BAG3 BANF1 CALR CCL2 CEBPB CEBPG CHAC1 CKS1B CNOT2 CNOT4 CNOT6 CXXC1 DCP1A DCP2 DCTN1 DDIT4 DDX10 DKC1 DNAJA4 DNAJB9 DNAJC3 EDC4 EDEM1 EEF2 EIF2AK3 EIF2S1 EIF4A1 EIF4A2 EIF4A3 EIF4E EIF4EBP1 EIF4G1 ERN1 ERO1A EXOC2 EXOSC1 EXOSC10 EXOSC2 EXOSC4 EXOSC5 EXOSC9 FKBP14 FUS GEMIN4 GOSR2 H2AX HERPUD1 HSP90B1 HSPA5 HSPA9 HYOU1 IARS1 IFIT1 IGFBP1 IMP3 KDELR3 KHSRP KIF5B LSM1 LSM4 MTHFD2 MTREX NABP1 NFYA NFYB NHP2 NOLC1 NOP14 NOP56 NPM1 PAIP1 PARN PDIA5 PDIA6 POP4 PREB PSAT1 RPS14 RRP9 SDAD1 SEC11A SEC31A SERP1 SHC1 SKIC3 SLC1A4 SLC30A5 SLC7A5 SPCS1 SPCS3 SRPRA SRPRB SSR1 STC2 TARS1 TATDN2 TSPYL2 TUBB2A VEGFA WFS1 WIPI1 XBP1 XPOT YIF1A YWHAZ ZBTB17""".split()
SOURCE_UNION = sorted(set(APOPTOSIS) | set(UPR))
MASKED_LIST_GENES = ["RPS14"]
SCORING_GENES = sorted(set(SOURCE_UNION) - set(MASKED_LIST_GENES))

assert len(APOPTOSIS) == 161
assert len(UPR) == 113
assert len(SOURCE_UNION) == 272
assert len(SCORING_GENES) == 271

GAO_STUDY_COLS = (
    "study_chaffin_2022",
    "study_koenig_2022",
    "study_penn",
    "study_simonson_2023",
    "study_kuppe_2022",
    "study_sim_2021",
)

REMOVAL_COLUMNS = [
    "dataset", "cell_type", "row_type", "status", "gene", "n_donors", "n_cells_unremoved",
    "n_cells_removed", "n_high_cells_removed", "score_threshold", "n_genes_unremoved",
    "n_genes_removed", "n_common_genes", "n_pass_unremoved", "n_pass_removed",
    "n_pass_to_fail_predicted", "n_pass_to_fail_not_attenuated", "n_fail_to_pass",
    "n_pass_to_pass", "n_fail_to_fail", "slope_unremoved", "se_unremoved",
    "ci95_low_unremoved", "ci95_high_unremoved", "pvalue_unremoved", "q_bh_unremoved",
    "mde80_unremoved", "slope_removed", "se_removed", "ci95_low_removed",
    "ci95_high_removed", "pvalue_removed", "q_bh_removed", "mde80_removed",
    "slope_change_removed_minus_unremoved", "absolute_slope_attenuated", "pass_fail_change",
    "predicted_pass_to_fail", "covariates_unremoved", "covariates_removed",
    "omitted_covariates_unremoved", "omitted_covariates_removed", "gene_list_source",
]


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def git_file_provenance(path: Path) -> str:
    try:
        rel = str(path.relative_to(REPO_ROOT))
    except ValueError:
        return "external_data_cache"
    try:
        out = subprocess.check_output(
            ["git", "log", "-n", "1", "--format=%h (%cs)", "--", rel],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if out:
            status = subprocess.check_output(
                ["git", "status", "--short", "--", rel],
                cwd=REPO_ROOT,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
            return f"committed {out}" + (" [modified in working tree]" if status else "")
        return "untracked working-tree output from parallel Phase 2 runner"
    except Exception:
        return "unknown"


def bh_manual(pvalues: np.ndarray) -> np.ndarray:
    """Independent Benjamini-Hochberg step-up adjusted p-values."""
    p = np.asarray(pvalues, dtype=np.float64)
    m = p.size
    if m == 0:
        return np.array([], dtype=np.float64)
    if not np.isfinite(p).all():
        raise ValueError("Non-finite p-value passed to bh_manual")
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adj = ranked * float(m) / np.arange(1, m + 1, dtype=np.float64)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty(m, dtype=np.float64)
    out[order] = adj
    return out


def compute_bh_checked(pvalues: np.ndarray) -> tuple[np.ndarray, float]:
    """Compute BH q-values with statsmodels and assert agreement with manual BH."""
    p = np.asarray(pvalues, dtype=np.float64)
    q_sm = multipletests(p, method="fdr_bh")[1]
    q_man = bh_manual(p)
    max_diff = float(np.max(np.abs(q_sm - q_man))) if p.size else 0.0
    if max_diff > 1e-12:
        raise AssertionError(f"BH recomputation mismatch: max_abs_diff={max_diff:.3e}")
    return q_sm, max_diff


def load_bh_audit() -> list[dict[str, Any]]:
    path = CACHE_ROOT / "bh_audit.json"
    if not path.exists():
        return []
    return json.loads(path.read_text(encoding="utf-8"))


def save_bh_audit(entries: list[dict[str, Any]], stage_key: str) -> None:
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    path = CACHE_ROOT / "bh_audit.json"
    existing = [e for e in load_bh_audit() if e.get("stage") != stage_key]
    existing.extend(entries)
    path.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")


def available_memory_gib() -> float:
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            return float(line.split()[1]) / (1024.0 ** 2)
    raise RuntimeError("MemAvailable is missing from /proc/meminfo")


def wait_for_memory(label: str, min_gib: float = MIN_AVAILABLE_MEM_GIB) -> None:
    while True:
        available = available_memory_gib()
        if available >= min_gib:
            log(f"Memory check for {label}: {available:.1f} GiB available")
            return
        log(f"Memory guard for {label}: only {available:.1f} GiB available (< {min_gib:.1f}); waiting 30 seconds")
        time.sleep(30)


def load_info(ds: str) -> dict[str, Any]:
    path = BLOCK1_CACHE / ds / "CACHE_DONE.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing required Block 1 cache metadata: {path}")
    info = json.loads(path.read_text(encoding="utf-8"))
    if not (BLOCK1_CACHE / ds / "donors_meta.csv").exists():
        raise FileNotFoundError(f"Missing donor metadata in {path.parent}")
    return info


def donor_metadata_mature(ds: str, info: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load Block 1 donor metadata and restrict strictly to age_years >= 30.0.

    Preserves each mature donor's original integer donor_idx from block1_cache so shard
    cell indexing and fixed-seed rarefaction SeedSequence([20261117, 7, dataset_code, donor_idx])
    match the frozen protocol exactly while excluding every under-30 donor before scoring.
    """
    all_donors = pd.read_csv(BLOCK1_CACHE / ds / "donors_meta.csv")
    all_donors["donor_id"] = all_donors["donor_id"].astype(str)
    if list(all_donors["donor_id"]) != list(info["donor_ids"]):
        raise AssertionError(f"Donor order mismatch in {ds} Block 1 cache")
    if all_donors["age_years"].isna().any():
        raise AssertionError(f"Missing age_years in {ds} Block 1 cache donor metadata")
    all_donors["donor_idx"] = np.arange(len(all_donors), dtype=np.int32)
    all_donors["study_or_pool"] = (
        all_donors["study"].astype(str)
        if ds == "gao"
        else all_donors["pool_number"].astype(int).astype(str)
    )
    mature = all_donors[all_donors["age_years"].astype(float) >= MIN_AGE_YEARS].copy().reset_index(drop=True)
    if mature.empty:
        raise RuntimeError(f"No mature donors (>= {MIN_AGE_YEARS}) in {ds}")
    assert float(mature["age_years"].min()) >= MIN_AGE_YEARS, f"Under-30 donor leaked in {ds}"
    audit = {
        "dataset": ds,
        "min_age_years": MIN_AGE_YEARS,
        "n_donors_original": int(len(all_donors)),
        "n_donors_under30_excluded": int(len(all_donors) - len(mature)),
        "n_donors_mature": int(len(mature)),
        "age_min_mature": float(mature["age_years"].min()),
        "age_max_mature": float(mature["age_years"].max()),
        "n_batches_mature": int(mature["study_or_pool"].nunique()),
    }
    return mature, audit


def fit_donor_ols(
    frame: pd.DataFrame,
    outcome: str,
    ds: str,
    adjusted: bool,
) -> dict[str, Any]:
    work = frame.dropna(subset=[outcome, "age_per_decade", "age_years"]).copy()
    if not work.empty and float(work["age_years"].min()) < MIN_AGE_YEARS:
        raise AssertionError(f"Under-30 donor leaked into fit_donor_ols for {ds}")
    base: dict[str, Any] = {
        "status": "EVALUATED",
        "n_donors": int(len(work)),
        "n_cells": int(work["n_cells"].sum()) if "n_cells" in work else 0,
        "min_donor_age": float(work["age_years"].min()) if len(work) else np.nan,
        "max_donor_age": float(work["age_years"].max()) if len(work) else np.nan,
        "residual_df": np.nan,
        "slope_per_decade": np.nan,
        "se": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "pvalue": np.nan,
        "outcome_sd": np.nan,
        "mde80": np.nan,
        "mde80_sd_per_decade": np.nan,
        "covariates_active": "none",
        "covariates_omitted": "none",
    }
    if len(work) < MIN_DONORS or work["age_per_decade"].nunique() < 2:
        base["status"] = "UNTESTABLE_DONOR_FLOOR"
        return base
    batch_col = "study_or_pool" if adjusted else None
    min_df = MIN_RESID_DF[ds] if adjusted else 5
    sub, active, omitted, status = _prune_covariate_matrix(
        work,
        nonbatch_covs=["sex_male"],
        batch_col=batch_col,
        min_batch_donors=2,
        min_residual_df=min_df,
        ref_level=BATCH_REF[ds] if adjusted else None,
        require_is_comp_collinear_check="age_per_decade",
        atomic_batch_block=True,
    )
    if status != "evaluated":
        base["status"] = f"UNTESTABLE_{status.upper()}"
        base["covariates_omitted"] = ";".join(omitted) if omitted else "none"
        return base
    assert float(sub["age_years"].min()) >= MIN_AGE_YEARS, f"Under-30 donor leaked after pruning in {ds}"
    y = sub[outcome].to_numpy(dtype=float)
    sd = float(np.std(y, ddof=1))
    base["n_donors"] = int(len(sub))
    base["n_cells"] = int(sub["n_cells"].sum())
    base["min_donor_age"] = float(sub["age_years"].min())
    base["max_donor_age"] = float(sub["age_years"].max())
    base["outcome_sd"] = sd
    if not np.isfinite(sd) or sd <= 0:
        base["status"] = "UNTESTABLE_ZERO_VAR"
        return base
    X = sm.add_constant(sub[["age_per_decade", *active]].astype(float), has_constant="add")
    fit = sm.OLS(y, X).fit()
    df = int(fit.df_resid)
    if df <= 0:
        base["status"] = "UNTESTABLE_INSUFFICIENT_DF"
        return base
    slope = float(fit.params["age_per_decade"])
    se = float(fit.bse["age_per_decade"])
    t975 = float(stats.t.ppf(0.975, df))
    t80 = float(stats.t.ppf(0.80, df))
    mde = (t975 + t80) * se
    base.update(
        residual_df=df,
        slope_per_decade=slope,
        se=se,
        ci95_low=slope - t975 * se,
        ci95_high=slope + t975 * se,
        pvalue=float(fit.pvalues["age_per_decade"]),
        mde80=mde,
        mde80_sd_per_decade=mde / sd,
        covariates_active=";".join(active) if active else "none",
        covariates_omitted=";".join(omitted) if omitted else "none",
    )
    return base


def run_scores(ds: str, force: bool = False) -> dict[str, Any]:
    """Compute per-cell full-depth and rarefied-1500 death/stress scores for mature donors only."""
    info = load_info(ds)
    donors, donor_audit = donor_metadata_mature(ds, info)
    mature_idx = np.sort(donors["donor_idx"].to_numpy(dtype=np.int32))
    out_dir = CACHE_ROOT / ds / "scores"
    out_dir.mkdir(parents=True, exist_ok=True)
    done_path = out_dir / "SCORES_DONE.json"

    symbols = list(info["retained_symbols"])
    symbol_to_idx = {gene: i for i, gene in enumerate(symbols)}
    score_idx = np.array([symbol_to_idx[g] for g in SCORING_GENES if g in symbol_to_idx], dtype=np.int32)
    absent = [g for g in SCORING_GENES if g not in symbol_to_idx]
    n_genes = len(symbols)
    n_shards = int(info["n_shards"])

    if done_path.exists() and not force:
        prev = json.loads(done_path.read_text(encoding="utf-8"))
        if prev.get("min_age_years") == MIN_AGE_YEARS and prev.get("n_donors_mature") == donor_audit["n_donors_mature"]:
            log(f"{ds}: reusing verified mature score cache ({prev['n_mature_cells_scored']} cells)")
            return prev

    total_cells_seen = 0
    total_cells_mature = 0
    total_cells_under30_dropped = 0
    min_umi_seen = None

    for shard in range(n_shards):
        out_path = out_dir / f"scores_shard{shard:02d}.npz"
        files = sorted((BLOCK1_CACHE / ds).glob(f"shard{shard:02d}_chunk*.npz"))
        if not files:
            continue
        mats: list[sp.csr_matrix] = []
        meta_parts: list[np.ndarray] = []
        for file in files:
            z = np.load(file)
            z_donor = z["donor"].astype(np.int32)
            total_cells_seen += len(z_donor)
            keep = np.isin(z_donor, mature_idx)
            n_keep = int(keep.sum())
            total_cells_under30_dropped += int(len(z_donor) - n_keep)
            if n_keep == 0:
                continue
            mat = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=(len(z_donor), n_genes))
            mats.append(mat[keep])
            meta_parts.append(
                np.column_stack([z_donor[keep], z["ct"][keep], z["umi"][keep], z["obs_idx"][keep]])
            )
        if not mats:
            np.savez_compressed(
                out_path,
                donor=np.array([], dtype=np.int32),
                ct=np.array([], dtype=np.int32),
                umi=np.array([], dtype=np.int64),
                obs_idx=np.array([], dtype=np.int32),
                list_counts=np.array([], dtype=np.int64),
                full_score=np.array([], dtype=np.float64),
                rare_score=np.array([], dtype=np.float64),
            )
            continue
        mat = sp.vstack(mats, format="csr")
        meta = np.vstack(meta_parts)
        order = np.argsort(meta[:, 3], kind="mergesort")
        mat = mat[order]
        meta = meta[order]
        del mats, meta_parts
        donor = meta[:, 0].astype(np.int32)
        ct = meta[:, 1].astype(np.int32)
        umi = meta[:, 2].astype(np.int64)
        obs_idx = meta[:, 3].astype(np.int32)

        assert np.isin(donor, mature_idx).all(), f"Under-30 donor leaked into {ds} shard {shard}"
        assert (umi >= 1500).all(), f"Cell with post-mask UMI < 1500 found in {ds} shard {shard}"
        shard_min_umi = int(umi.min())
        min_umi_seen = shard_min_umi if min_umi_seen is None else min(min_umi_seen, shard_min_umi)
        total_cells_mature += len(donor)

        list_counts = np.rint(np.asarray(mat[:, score_idx].sum(axis=1)).ravel()).astype(np.int64)
        full_score = list_counts.astype(np.float64) / (float(len(SCORING_GENES)) * umi.astype(np.float64))
        rare_score = np.empty(len(meta), dtype=np.float64)
        for d in np.unique(donor):
            rows = np.flatnonzero(donor == d)
            rng = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, 7, DATASET_CODE[ds], int(d)]))
            probability = 1500.0 / umi[rows].astype(np.float64)
            rare_counts = rng.binomial(list_counts[rows], probability)
            rare_score[rows] = rare_counts.astype(np.float64) / float(len(SCORING_GENES) * 1500)
        np.savez_compressed(
            out_path,
            donor=donor,
            ct=ct,
            umi=umi,
            obs_idx=obs_idx,
            list_counts=list_counts,
            full_score=full_score,
            rare_score=rare_score,
        )
        del mat, meta
        log(f"{ds} shard {shard + 1}/{n_shards}: scored {len(donor)} mature cells")

    audit = {
        **donor_audit,
        "source_union_genes": len(SOURCE_UNION),
        "masked_source_genes": MASKED_LIST_GENES,
        "scoring_genes": len(SCORING_GENES),
        "present_scoring_genes": int(len(score_idx)),
        "absent_scoring_genes": absent,
        "n_shards": n_shards,
        "n_cells_seen_in_block1_cache": int(total_cells_seen),
        "n_mature_cells_scored": int(total_cells_mature),
        "n_under30_cells_excluded": int(total_cells_under30_dropped),
        "min_post_mask_umi_scored": int(min_umi_seen) if min_umi_seen is not None else None,
        "master_seed": MASTER_SEED,
    }
    done_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    log(
        f"{ds}: mature score cache complete ({total_cells_mature} mature cells scored, "
        f"{total_cells_under30_dropped} under-30 cells excluded)"
    )
    return audit


def load_scores(ds: str, mature_idx: np.ndarray) -> dict[str, np.ndarray]:
    score_dir = CACHE_ROOT / ds / "scores"
    if not (score_dir / "SCORES_DONE.json").exists():
        raise FileNotFoundError(f"Missing completed mature score cache for {ds}")
    files = sorted(score_dir.glob("scores_shard*.npz"))
    pieces = {key: [] for key in ("donor", "ct", "umi", "obs_idx", "list_counts", "full_score", "rare_score")}
    for file in files:
        z = np.load(file)
        if len(z["donor"]) == 0:
            continue
        for key in pieces:
            pieces[key].append(z[key])
    out = {key: np.concatenate(value) for key, value in pieces.items()}
    assert np.isin(out["donor"], mature_idx).all(), f"Under-30 donor present in loaded mature scores for {ds}"
    assert (out["umi"] >= 1500).all(), f"UMI < 1500 present in loaded mature scores for {ds}"
    return out


def apply_test1_bh(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    out = frame.copy()
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    bh_entries: list[dict[str, Any]] = []
    eligible = (out["status"] == "EVALUATED") & np.isfinite(out["pvalue"].astype(float))
    for (ds, test), indices in out.loc[eligible].groupby(["dataset", "test"], sort=True).groups.items():
        idx = list(indices)
        pvals = out.loc[idx, "pvalue"].to_numpy(dtype=float)
        qvals, max_diff = compute_bh_checked(pvals)
        out.loc[idx, "family_size"] = len(idx)
        out.loc[idx, "q_bh"] = qvals
        bh_entries.append(
            {
                "stage": "test1",
                "dataset": str(ds),
                "family": str(test),
                "family_size": int(len(idx)),
                "n_q_lt_0_05": int((qvals < 0.05).sum()),
                "min_p": float(np.min(pvals)),
                "min_q": float(np.min(qvals)),
                "max_abs_diff_statsmodels_vs_manual": max_diff,
            }
        )
    return out, bh_entries


def run_test1() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for ds in ("gao", "yazar"):
        info = load_info(ds)
        donors, _ = donor_metadata_mature(ds, info)
        mature_idx = np.sort(donors["donor_idx"].to_numpy(dtype=np.int32))
        scores = load_scores(ds, mature_idx)
        (CACHE_ROOT / ds).mkdir(parents=True, exist_ok=True)
        for ci, cell_type in enumerate(info["cell_types"]):
            ct_mask = scores["ct"] == ci
            ct_donor = scores["donor"][ct_mask]
            ct_full = scores["full_score"][ct_mask]
            ct_rare = scores["rare_score"][ct_mask]
            all_donors, all_counts = np.unique(ct_donor, return_counts=True)
            count_map = dict(zip(all_donors.tolist(), all_counts.tolist()))
            eligible_donors = np.array(sorted(d for d, n in count_map.items() if n >= MIN_CELLS), dtype=np.int32)
            eligible_cell_mask = np.isin(ct_donor, eligible_donors)
            ed = ct_donor[eligible_cell_mask]
            ef = ct_full[eligible_cell_mask]
            er = ct_rare[eligible_cell_mask]
            base = {
                "dataset": ds,
                "cell_type": cell_type,
                "age_restriction": "age_ge_30",
                "n_donors_with_cells": int(len(all_donors)),
                "n_eligible_cells": int(len(ed)),
                "n_source_genes": len(SOURCE_UNION),
                "n_scoring_genes": len(SCORING_GENES),
                "masked_list_genes": ";".join(MASKED_LIST_GENES),
                "gene_list_source": GENE_SOURCE,
            }
            donor_rows: dict[str, pd.DataFrame] = {}
            for arm, values in (("full", ef), ("rare", er)):
                if len(ed) == 0:
                    donor_frame = pd.DataFrame()
                    threshold = np.nan
                    n_high = 0
                else:
                    threshold = float(np.quantile(values, 0.90, method="linear"))
                    high = values >= threshold
                    unique_d, inverse = np.unique(ed, return_inverse=True)
                    n_cells = np.bincount(inverse)
                    n_high_by_d = np.bincount(inverse, weights=high.astype(float)).astype(int)
                    dmeta = donors.set_index("donor_idx").loc[unique_d].reset_index()
                    assert float(dmeta["age_years"].min()) >= MIN_AGE_YEARS
                    donor_frame = dmeta.copy()
                    donor_frame["n_cells"] = n_cells
                    donor_frame["n_high"] = n_high_by_d
                    donor_frame["fraction_high"] = n_high_by_d / n_cells
                    donor_frame["score_threshold"] = threshold
                    n_high = int(high.sum())
                donor_rows[arm] = donor_frame
                for model, adjusted in (("primary", True), ("unadjusted", False)):
                    if arm == "rare" and model == "unadjusted":
                        continue
                    test = f"test1_{model}_{arm}" if arm == "full" else "test1_rarefied_1500"
                    if donor_frame.empty:
                        fit = {
                            "status": "UNTESTABLE_DONOR_FLOOR",
                            "n_donors": 0,
                            "n_cells": 0,
                            "min_donor_age": np.nan,
                            "max_donor_age": np.nan,
                        }
                    else:
                        fit = fit_donor_ols(donor_frame, "fraction_high", ds, adjusted)
                    rows.append(
                        {
                            **base,
                            "test": test,
                            "model": model,
                            "depth": "full_depth" if arm == "full" else "rarefied_1500",
                            "score_threshold": threshold,
                            "n_high_cells": n_high,
                            **fit,
                        }
                    )
            if not donor_rows["full"].empty:
                tmp = donor_rows["full"][["donor_idx", "age_years", "n_cells", "n_high"]].copy()
                tmp.to_csv(CACHE_ROOT / ds / f"donor_high_{ci:02d}.csv", index=False)

    result, bh_entries = apply_test1_bh(pd.DataFrame(rows))
    save_bh_audit(bh_entries, "test1")
    preferred = [
        "dataset", "cell_type", "age_restriction", "test", "model", "depth", "status",
        "slope_per_decade", "se", "ci95_low", "ci95_high", "pvalue", "q_bh", "family_size",
        "n_donors", "min_donor_age", "max_donor_age", "n_cells", "n_donors_with_cells",
        "n_eligible_cells", "n_high_cells", "score_threshold", "outcome_sd", "mde80",
        "mde80_sd_per_decade", "residual_df", "covariates_active", "covariates_omitted",
        "n_source_genes", "n_scoring_genes", "masked_list_genes", "gene_list_source",
    ]
    for col in preferred:
        if col not in result:
            result[col] = np.nan
    result = result[preferred].sort_values(["dataset", "test", "cell_type"]).reset_index(drop=True)
    write_clean_csv(result, DEATH_PATH)
    for ds in ("gao", "yazar"):
        sub = result[(result["dataset"] == ds) & (result["test"] == "test1_primary_full")]
        log(
            f"{ds}: mature test 1 evaluated={int((sub['status'] == 'EVALUATED').sum())}, "
            f"BH passing={int((sub['q_bh'] < 0.05).sum())}"
        )
    return result


def load_shard_for_removal(
    ds: str,
    shard: int,
    info: dict[str, Any],
    mature_idx: np.ndarray,
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray] | None:
    files = sorted((BLOCK1_CACHE / ds).glob(f"shard{shard:02d}_chunk*.npz"))
    if not files:
        return None
    mats: list[sp.csr_matrix] = []
    meta_parts: list[np.ndarray] = []
    n_genes = len(info["retained_symbols"])
    for file in files:
        z = np.load(file)
        z_donor = z["donor"].astype(np.int32)
        keep = np.isin(z_donor, mature_idx)
        if not keep.any():
            continue
        mat = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=(len(z_donor), n_genes))
        mats.append(mat[keep])
        meta_parts.append(np.column_stack([z_donor[keep], z["ct"][keep], z["obs_idx"][keep]]))
    if not mats:
        return None
    mat = sp.vstack(mats, format="csr")
    meta = np.vstack(meta_parts)
    order = np.argsort(meta[:, 2], kind="mergesort")
    mat = mat[order]
    meta = meta[order]
    donor = meta[:, 0].astype(np.int32)
    ct = meta[:, 1].astype(np.int32)
    obs_idx = meta[:, 2].astype(np.int32)
    score_path = CACHE_ROOT / ds / "scores" / f"scores_shard{shard:02d}.npz"
    if not score_path.exists():
        raise FileNotFoundError(f"Missing mature score shard: {score_path}")
    score = np.load(score_path)
    if not np.array_equal(obs_idx, score["obs_idx"]):
        raise AssertionError(f"Cell order mismatch between count and mature score shard: {ds} shard {shard}")
    if not np.array_equal(donor, score["donor"]) or not np.array_equal(ct, score["ct"]):
        raise AssertionError(f"Metadata mismatch between count and mature score shard: {ds} shard {shard}")
    assert np.isin(donor, mature_idx).all(), f"Under-30 donor leaked in removal shard {ds} {shard}"
    return mat, donor, ct, score["full_score"]


def prepare_removal_specs(
    ds: str,
    hits: pd.DataFrame,
    info: dict[str, Any],
    donors: pd.DataFrame,
) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    donor_lookup = donors.set_index("donor_idx")
    for _, hit in hits.iterrows():
        cell_type = str(hit["cell_type"])
        ci = info["cell_types"].index(cell_type)
        donor_high_path = CACHE_ROOT / ds / f"donor_high_{ci:02d}.csv"
        if not donor_high_path.exists():
            raise FileNotFoundError(f"Missing mature test 1 donor fractions: {donor_high_path}")
        high = pd.read_csv(donor_high_path)
        assert float(high["age_years"].min()) >= MIN_AGE_YEARS, f"Under-30 donor in {donor_high_path}"
        high["n_removed"] = high["n_cells"] - high["n_high"]
        eligible = high[(high["n_cells"] >= MIN_CELLS) & (high["n_removed"] >= MIN_CELLS)].copy()
        common = np.sort(eligible["donor_idx"].to_numpy(dtype=np.int32))
        if ds == "yazar" and len(common):
            pools = donor_lookup.loc[common, "pool_number"].astype(int)
            pool_counts = pools.value_counts()
            common = np.sort(pools[pools.map(pool_counts) >= 2].index.to_numpy(dtype=np.int32))
        elif ds == "gao" and len(common):
            studies = donor_lookup.loc[common, "study"].astype(str)
            study_counts = studies.value_counts()
            common = np.sort(studies[studies.map(study_counts) >= 2].index.to_numpy(dtype=np.int32))
        if len(common):
            assert float(donor_lookup.loc[common, "age_years"].min()) >= MIN_AGE_YEARS
        eligible = eligible.set_index("donor_idx").loc[common] if len(common) else eligible.iloc[0:0]
        base = {column: np.nan for column in REMOVAL_COLUMNS}
        base.update(
            dataset=ds,
            cell_type=cell_type,
            row_type="summary",
            status="EVALUATED" if len(common) >= MIN_DONORS else "UNTESTABLE_DONOR_FLOOR",
            gene="ALL_GENES",
            n_donors=int(len(common)),
            n_cells_unremoved=int(eligible["n_cells"].sum()) if len(common) else 0,
            n_cells_removed=int(eligible["n_removed"].sum()) if len(common) else 0,
            n_high_cells_removed=int(eligible["n_high"].sum()) if len(common) else 0,
            score_threshold=float(hit["score_threshold"]),
            n_genes_unremoved=0,
            n_genes_removed=0,
            n_common_genes=0,
            n_pass_unremoved=0,
            n_pass_removed=0,
            n_pass_to_fail_predicted=0,
            n_pass_to_fail_not_attenuated=0,
            n_fail_to_pass=0,
            n_pass_to_pass=0,
            n_fail_to_fail=0,
            gene_list_source=GENE_SOURCE,
        )
        specs.append(
            {
                "dataset": ds,
                "cell_type": cell_type,
                "ci": ci,
                "threshold": float(hit["score_threshold"]),
                "common": common,
                "summary": base,
                "unremoved_blocks": [],
                "removed_blocks": [],
            }
        )
    return specs


def aggregate_removal_blocks(
    ds: str,
    specs: list[dict[str, Any]],
    info: dict[str, Any],
    mature_idx: np.ndarray,
) -> None:
    active = [spec for spec in specs if spec["summary"]["status"] == "EVALUATED"]
    if not active:
        return
    n_donors = len(info["donor_ids"])
    n_shards = int(info["n_shards"])
    for shard in range(n_shards):
        loaded = load_shard_for_removal(ds, shard, info, mature_idx)
        if loaded is None:
            continue
        mat, donor, ct, full_score = loaded
        for spec in active:
            local = spec["common"][spec["common"] % n_shards == shard]
            if len(local) == 0:
                continue
            rows = np.flatnonzero((ct == spec["ci"]) & np.isin(donor, local))
            if len(rows) == 0:
                raise AssertionError(f"No mature cells for {ds} {spec['cell_type']} in expected shard {shard}")
            lookup = np.full(n_donors, -1, dtype=np.int32)
            lookup[local] = np.arange(len(local), dtype=np.int32)
            columns = lookup[donor[rows]]
            if (columns < 0).any():
                raise AssertionError("Donor mapping failure in mature removal aggregation")
            indicator = sp.csr_matrix(
                (np.ones(len(rows), dtype=np.float64), (np.arange(len(rows)), columns)),
                shape=(len(rows), len(local)),
            )
            ct_mat = mat[rows]
            unremoved = (indicator.T @ ct_mat).T.tocsr()
            high = full_score[rows] >= spec["threshold"]
            removed = (indicator[~high].T @ ct_mat[~high]).T.tocsr()
            spec["unremoved_blocks"].append((local, unremoved))
            spec["removed_blocks"].append((local, removed))
        del mat, donor, ct, full_score
        log(f"{ds}: aggregated mature removal counts for shard {shard + 1}/{n_shards}")


def combine_count_blocks(
    blocks: list[tuple[np.ndarray, sp.csr_matrix]],
    common: np.ndarray,
    n_genes: int,
) -> np.ndarray:
    if not blocks:
        raise ValueError("No pseudobulk blocks for an evaluable removal cell type")
    donor_order = np.concatenate([block[0] for block in blocks])
    matrix = sp.hstack([block[1] for block in blocks], format="csr")
    order = np.argsort(donor_order, kind="mergesort")
    donor_order = donor_order[order]
    matrix = matrix[:, order]
    if not np.array_equal(donor_order, common):
        raise AssertionError("Aggregated donor columns do not match the common mature donor set")
    if matrix.shape != (n_genes, len(common)):
        raise AssertionError(f"Unexpected pseudobulk shape {matrix.shape}")
    return matrix.toarray().astype(np.float64, copy=False)


def removal_meta(
    ds: str,
    donors: pd.DataFrame,
    common: np.ndarray,
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    meta = donors.set_index("donor_idx").loc[common].copy()
    assert float(meta["age_years"].min()) >= MIN_AGE_YEARS, f"Under-30 donor in removal_meta for {ds}"
    meta.index = meta["donor_id"].astype(str)
    if ds == "gao":
        return meta, ("sex_male", "log_lib", *GAO_STUDY_COLS)
    pools = meta["pool_number"].astype(int)
    counts = pools.value_counts()
    ordered = sorted(counts.index.tolist(), key=lambda p: (-int(counts[p]), int(p)))
    reference = ordered[0]
    all_pools = sorted(donors["pool_number"].dropna().astype(int).unique().tolist())
    pool_cols: list[str] = []
    for pool in all_pools:
        column = f"pool_{pool}"
        meta[column] = 0.0 if pool == reference else (pools == pool).astype(float)
        pool_cols.append(column)
    return meta, ("sex_male", "log_lib", *pool_cols)


def fit_removal_arm(
    counts: np.ndarray,
    symbols: list[str],
    donor_ids: list[str],
    meta: pd.DataFrame,
    covariates: tuple[str, ...],
    cell_type: str,
    ds: str,
    arm_label: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    assert float(meta["age_years"].min()) >= MIN_AGE_YEARS
    raw = pd.DataFrame(counts, index=symbols, columns=donor_ids)
    log2cpm, library = counts_to_log2cpm(raw)
    arm_meta = meta.copy()
    arm_meta["lib_size"] = library.loc[donor_ids].to_numpy(dtype=float)
    arm_meta["log_lib"] = np.log1p(arm_meta["lib_size"].to_numpy(dtype=float))
    genes, _, _ = run_family_b_genes_and_concordance(
        {cell_type: log2cpm},
        {cell_type: arm_meta},
        age_col="age_per_decade",
        candidate_covariates=covariates,
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=GENE_MIN_RESID_DF[ds],
        is_unbiased_single_cell=True,
        concordance_perms=0,
        seed=MASTER_SEED,
    )
    if genes.empty:
        return genes, {
            "stage": "test2",
            "dataset": ds,
            "family": f"{cell_type} ({arm_label})",
            "family_size": 0,
            "n_q_lt_0_05": 0,
            "min_p": None,
            "min_q": None,
            "max_abs_diff_statsmodels_vs_manual": 0.0,
        }
    pvals = genes["pvalue"].to_numpy(dtype=float)
    q_check, max_diff = compute_bh_checked(pvals)
    q_existing = genes["q_bh_within_group"].to_numpy(dtype=float)
    existing_diff = float(np.max(np.abs(q_check - q_existing)))
    if existing_diff > 1e-12:
        raise AssertionError(
            f"q_bh_within_group mismatch vs checked BH in {ds} {cell_type} {arm_label}: {existing_diff:.3e}"
        )
    audit = {
        "stage": "test2",
        "dataset": ds,
        "family": f"{cell_type} ({arm_label})",
        "family_size": int(len(genes)),
        "n_q_lt_0_05": int((q_existing < 0.05).sum()),
        "min_p": float(np.min(pvals)),
        "min_q": float(np.min(q_existing)),
        "max_abs_diff_statsmodels_vs_manual": max(max_diff, existing_diff),
    }
    return genes, audit


def rename_gene_arm(frame: pd.DataFrame, arm: str) -> pd.DataFrame:
    columns = {
        "gene": "gene",
        "effect": f"slope_{arm}",
        "se": f"se_{arm}",
        "ci95_low": f"ci95_low_{arm}",
        "ci95_high": f"ci95_high_{arm}",
        "pvalue": f"pvalue_{arm}",
        "q_bh_within_group": f"q_bh_{arm}",
        "covariates": f"covariates_{arm}",
        "omitted_covariates": f"omitted_covariates_{arm}",
        "residual_df": f"residual_df_{arm}",
    }
    if frame.empty:
        return pd.DataFrame(columns=list(columns.values()))
    out = frame[list(columns)].rename(columns=columns).copy()
    df = frame["residual_df"].to_numpy(dtype=float)
    critical = stats.t.ppf(0.975, df)
    t80 = stats.t.ppf(0.80, df)
    out[f"mde80_{arm}"] = (critical + t80) * out[f"se_{arm}"].to_numpy(dtype=float)
    return out


def compare_removal_arms(
    spec: dict[str, Any],
    unremoved: pd.DataFrame,
    removed: pd.DataFrame,
) -> pd.DataFrame:
    left = rename_gene_arm(unremoved, "unremoved")
    right = rename_gene_arm(removed, "removed")
    genes = left.merge(right, on="gene", how="outer", validate="one_to_one")
    both = genes["slope_unremoved"].notna() & genes["slope_removed"].notna()
    genes["row_type"] = "gene"
    genes["status"] = np.select(
        [both, genes["slope_unremoved"].notna(), genes["slope_removed"].notna()],
        ["EVALUATED_BOTH", "UNREMOVED_ONLY", "REMOVED_ONLY"],
        default="NOT_EXPRESSED_EITHER",
    )
    genes["slope_change_removed_minus_unremoved"] = genes["slope_removed"] - genes["slope_unremoved"]
    genes["absolute_slope_attenuated"] = np.where(
        both,
        genes["slope_removed"].abs() < genes["slope_unremoved"].abs(),
        np.nan,
    )
    unremoved_pass = genes["q_bh_unremoved"] < 0.05
    removed_pass = genes["q_bh_removed"] < 0.05
    pass_to_fail = both & unremoved_pass & ~removed_pass
    fail_to_pass = both & ~unremoved_pass & removed_pass
    predicted = pass_to_fail & (genes["absolute_slope_attenuated"] == True)  # noqa: E712
    genes["pass_fail_change"] = "NOT_COMPARABLE"
    genes.loc[both & ~unremoved_pass & ~removed_pass, "pass_fail_change"] = "FAIL_TO_FAIL"
    genes.loc[both & unremoved_pass & removed_pass, "pass_fail_change"] = "PASS_TO_PASS"
    genes.loc[pass_to_fail & ~predicted, "pass_fail_change"] = "PASS_TO_FAIL_NOT_ATTENUATED"
    genes.loc[predicted, "pass_fail_change"] = "PASS_TO_FAIL_PREDICTED"
    genes.loc[fail_to_pass, "pass_fail_change"] = "FAIL_TO_PASS"
    genes["predicted_pass_to_fail"] = np.where(both, predicted.astype(float), np.nan)
    summary = spec["summary"].copy()
    summary.update(
        status="EVALUATED" if int(both.sum()) > 0 else "UNTESTABLE_NO_COMMON_GENES",
        n_genes_unremoved=int(len(unremoved)),
        n_genes_removed=int(len(removed)),
        n_common_genes=int(both.sum()),
        n_pass_unremoved=int((both & unremoved_pass).sum()),
        n_pass_removed=int((both & removed_pass).sum()),
        n_pass_to_fail_predicted=int(predicted.sum()),
        n_pass_to_fail_not_attenuated=int((pass_to_fail & ~predicted).sum()),
        n_fail_to_pass=int(fail_to_pass.sum()),
        n_pass_to_pass=int((both & unremoved_pass & removed_pass).sum()),
        n_fail_to_fail=int((both & ~unremoved_pass & ~removed_pass).sum()),
    )
    for column in REMOVAL_COLUMNS:
        if column not in genes:
            genes[column] = np.nan
    shared_columns = {
        "dataset", "cell_type", "n_donors", "n_cells_unremoved", "n_cells_removed",
        "n_high_cells_removed", "score_threshold", "n_genes_unremoved", "n_genes_removed",
        "n_common_genes", "n_pass_unremoved", "n_pass_removed", "n_pass_to_fail_predicted",
        "n_pass_to_fail_not_attenuated", "n_fail_to_pass", "n_pass_to_pass", "n_fail_to_fail",
        "gene_list_source",
    }
    for column in shared_columns:
        genes[column] = summary[column]
    summary_row = pd.DataFrame([summary])
    out = pd.concat([summary_row, genes], ignore_index=True, sort=False)
    for column in REMOVAL_COLUMNS:
        if column not in out:
            out[column] = np.nan
    return out[REMOVAL_COLUMNS]


def no_eligible_removal_row() -> pd.DataFrame:
    row = {col: np.nan for col in REMOVAL_COLUMNS}
    row.update(
        dataset="NONE",
        cell_type="NONE_QUALIFIED_AT_30_PLUS",
        row_type="no_eligible_type",
        status="NO_TEST1_BH_HIT",
        gene="NONE",
        n_donors=0,
        n_cells_unremoved=0,
        n_cells_removed=0,
        n_high_cells_removed=0,
        n_genes_unremoved=0,
        n_genes_removed=0,
        n_common_genes=0,
        n_pass_unremoved=0,
        n_pass_removed=0,
        n_pass_to_fail_predicted=0,
        n_pass_to_fail_not_attenuated=0,
        n_fail_to_pass=0,
        n_pass_to_pass=0,
        n_fail_to_fail=0,
        gene_list_source=GENE_SOURCE,
    )
    return pd.DataFrame([row])[REMOVAL_COLUMNS]


def run_test2(death: pd.DataFrame | None = None) -> pd.DataFrame:
    if death is None:
        if not DEATH_PATH.exists():
            raise FileNotFoundError(f"Run test 1 first: {DEATH_PATH}")
        death = pd.read_csv(DEATH_PATH, float_precision="round_trip")
    pieces: list[pd.DataFrame] = []
    bh_entries: list[dict[str, Any]] = []
    for ds in ("gao", "yazar"):
        info = load_info(ds)
        donors, _ = donor_metadata_mature(ds, info)
        mature_idx = np.sort(donors["donor_idx"].to_numpy(dtype=np.int32))
        hits = (
            death[
                (death["dataset"] == ds)
                & (death["test"] == "test1_primary_full")
                & (death["status"] == "EVALUATED")
                & (death["q_bh"] < 0.05)
            ]
            .copy()
            .sort_values("cell_type")
        )
        log(f"{ds}: mature test 2 qualifying cell types={len(hits)}")
        if hits.empty:
            continue
        specs = prepare_removal_specs(ds, hits, info, donors)
        wait_for_memory(f"{ds} mature removal aggregation")
        aggregate_removal_blocks(ds, specs, info, mature_idx)
        symbols = list(info["retained_symbols"])
        for spec in specs:
            if spec["summary"]["status"] != "EVALUATED":
                pieces.append(pd.DataFrame([spec["summary"]]))
                continue
            common = spec["common"]
            donor_ids = donors.set_index("donor_idx").loc[common, "donor_id"].astype(str).tolist()
            unremoved_counts = combine_count_blocks(spec["unremoved_blocks"], common, len(symbols))
            removed_counts = combine_count_blocks(spec["removed_blocks"], common, len(symbols))
            meta, covariates = removal_meta(ds, donors, common)
            try:
                unremoved, bh_un = fit_removal_arm(
                    unremoved_counts, symbols, donor_ids, meta, covariates, spec["cell_type"], ds, "unremoved"
                )
                removed, bh_rm = fit_removal_arm(
                    removed_counts, symbols, donor_ids, meta, covariates, spec["cell_type"], ds, "removed"
                )
                bh_entries.extend([bh_un, bh_rm])
            except ValueError as exc:
                summary = spec["summary"].copy()
                summary["status"] = f"UNTESTABLE_MODEL_{exc}"
                pieces.append(pd.DataFrame([summary]))
                continue
            compared = compare_removal_arms(spec, unremoved, removed)
            pieces.append(compared)
            summary = compared.iloc[0]
            log(
                f"{ds} {spec['cell_type']}: n_donors={int(summary['n_donors'])}, "
                f"common genes={int(summary['n_common_genes'])}, "
                f"unremoved_pass={int(summary['n_pass_unremoved'])}, "
                f"removed_pass={int(summary['n_pass_removed'])}, "
                f"predicted pass-to-fail={int(summary['n_pass_to_fail_predicted'])}"
            )
            spec["unremoved_blocks"].clear()
            spec["removed_blocks"].clear()

    save_bh_audit(bh_entries, "test2")
    result = pd.concat(pieces, ignore_index=True, sort=False) if pieces else no_eligible_removal_row()
    for column in REMOVAL_COLUMNS:
        if column not in result:
            result[column] = np.nan
    result = result[REMOVAL_COLUMNS]
    if not result.empty and (result["row_type"] != "no_eligible_type").any():
        result["_row_order"] = result["row_type"].map({"summary": 0, "gene": 1}).fillna(2)
        result = result.sort_values(["dataset", "cell_type", "_row_order", "gene"]).drop(columns="_row_order")
    write_clean_csv(result.reset_index(drop=True), REMOVAL_PATH)
    return result.reset_index(drop=True)


def bool_column(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].astype(str).str.lower().isin(["true", "1"])


def load_mature_spread(ds: str) -> tuple[pd.DataFrame | None, str, dict[str, Any]]:
    path = MATURE_SPREAD_PATHS[ds]
    audit: dict[str, Any] = {
        "dataset": ds,
        "path": str(path.relative_to(REPO_ROOT)),
        "exists": path.exists(),
    }
    if not path.exists():
        return None, f"missing mature full-depth spread CSV ({path.relative_to(REPO_ROOT)})", audit
    df = pd.read_csv(path)
    req_cols = {
        "metric", "depth_mode", "analysis_role", "model_spec",
        "is_pooled", "is_stratified_global", "group",
        "slope_raw_per_decade", "se_raw", "pvalue", "q_bh_within_slice", "n_donors",
    }
    missing = sorted(req_cols - set(df.columns))
    if missing:
        return None, f"mature spread CSV missing required columns: {', '.join(missing)}", audit
    mask = (
        df["metric"].eq("shannon")
        & df["depth_mode"].eq("full_depth")
        & df["analysis_role"].eq("primary_full_depth")
        & df["model_spec"].eq("primary")
        & ~bool_column(df, "is_pooled")
        & ~bool_column(df, "is_stratified_global")
    )
    sub = df.loc[mask].copy()
    if sub.empty:
        return None, f"no matching mature primary full-depth Shannon rows in {path.relative_to(REPO_ROOT)}", audit
    if sub["group"].duplicated().any():
        raise ValueError(f"Duplicate mature spread cell-type keys in {path}")
    # Verify donor counts cannot exceed mature cohort sizes (109 Gao, 932 Yazar)
    max_mature = 109 if ds == "gao" else 932
    if int(sub["n_donors"].max()) > max_mature:
        return None, f"spread CSV n_donors max {int(sub['n_donors'].max())} exceeds mature donor ceiling {max_mature}", audit
    audit.update(
        n_rows_total=int(len(df)),
        n_rows_filtered=int(len(sub)),
        max_n_donors=int(sub["n_donors"].max()),
        sha256=sha256_file(path),
        provenance=git_file_provenance(path),
    )
    return sub.rename(columns={"group": "cell_type"}), "ok", audit


def load_mature_fractions(ds: str) -> tuple[pd.DataFrame | None, str, dict[str, Any]]:
    """Verify mature composition human CSV (c0a034a) and load mature cohort CLR fraction slopes.

    Never joins any legacy full-range fraction file.
    """
    audit: dict[str, Any] = {"dataset": ds, "checked_mature_comp_human": False}
    if MATURE_COMP_HUMAN_PATH.exists():
        comp = pd.read_csv(MATURE_COMP_HUMAN_PATH)
        audit["checked_mature_comp_human"] = True
        audit["mature_comp_human_rows"] = int(len(comp))
        audit["mature_comp_human_cohorts"] = sorted(comp["cohort"].astype(str).unique().tolist())
        # Verify whether Gao or Yazar exist in followup-mature-composition-human.csv
        target = "gao" if ds == "gao" else "yazar"
        hits = comp["cohort"].astype(str).str.lower().str.contains(target, na=False)
        if hits.any():
            raise RuntimeError(f"Unexpected {ds} rows in {MATURE_COMP_HUMAN_PATH}; inspect schema")
        audit["mature_comp_human_has_dataset"] = False
    else:
        audit["mature_comp_human_has_dataset"] = False

    path = MATURE_FRACTION_PATHS[ds]
    audit["path"] = str(path.relative_to(REPO_ROOT))
    audit["exists"] = path.exists()
    if not path.exists():
        return (
            None,
            (
                f"mature fraction slopes absent: {MATURE_COMP_HUMAN_PATH.relative_to(REPO_ROOT)} "
                f"contains only the 11 core brain/kidney/lung cohorts (no {ds} rows), and "
                f"{path.relative_to(REPO_ROOT)} is not yet present"
            ),
            audit,
        )
    df = pd.read_csv(path)
    expected_cohort = "Gao_Heart" if ds == "gao" else "Yazar_Blood"
    expected_arm = "pooled_7study_primary" if ds == "gao" else "pool_fe"
    expected_n = 109 if ds == "gao" else 932
    expected_excl = 8 if ds == "gao" else 49
    sub = df[
        df["cohort"].astype(str).eq(expected_cohort)
        & df["endpoint_family"].astype(str).eq("fractions_clr")
        & df["arm"].astype(str).eq(expected_arm)
    ].copy()
    if sub.empty:
        return (
            None,
            f"no matching mature fractions_clr ({expected_arm}) rows in {path.relative_to(REPO_ROOT)}",
            audit,
        )
    if sub["cell_type"].duplicated().any():
        raise ValueError(f"Duplicate mature fraction cell-type keys in {path}")
    if not (sub["n_donors"].astype(int) == expected_n).all():
        return (
            None,
            f"unexpected mature donor count in {path.relative_to(REPO_ROOT)} (expected {expected_n})",
            audit,
        )
    if not (sub["n_donors_30plus_excluded"].astype(int) == expected_excl).all():
        return (
            None,
            f"unexpected excluded under-30 donor count in {path.relative_to(REPO_ROOT)} (expected {expected_excl})",
            audit,
        )
    out = pd.DataFrame(
        {
            "cell_type": sub["cell_type"].astype(str),
            "effect_clr": sub["mature_effect"].astype(float),
            "se": sub["mature_se"].astype(float),
            "pvalue": sub["mature_p"].astype(float),
            "q_bh_family_a": sub["mature_q"].astype(float),
            "n_donors": sub["n_donors"].astype(int),
            "n_donors_30plus_excluded": sub["n_donors_30plus_excluded"].astype(int),
            "covariates": sub["covariates"].astype(str),
            "arm": sub["arm"].astype(str),
        }
    )
    audit.update(
        cohort=expected_cohort,
        arm=expected_arm,
        n_rows_filtered=int(len(out)),
        n_donors=expected_n,
        n_donors_30plus_excluded=expected_excl,
        sha256=sha256_file(path),
        provenance=git_file_provenance(path),
    )
    return out, "ok", audit


def run_narrowing() -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    part_b_audit: dict[str, Any] = {
        "mature_composition_human_path": str(MATURE_COMP_HUMAN_PATH.relative_to(REPO_ROOT)),
        "mature_composition_human_exists": MATURE_COMP_HUMAN_PATH.exists(),
        "mature_composition_human_sha256": sha256_file(MATURE_COMP_HUMAN_PATH) if MATURE_COMP_HUMAN_PATH.exists() else None,
        "mature_composition_human_provenance": git_file_provenance(MATURE_COMP_HUMAN_PATH) if MATURE_COMP_HUMAN_PATH.exists() else None,
        "datasets": {},
    }
    if MATURE_COMP_HUMAN_PATH.exists():
        comp = pd.read_csv(MATURE_COMP_HUMAN_PATH)
        part_b_audit["mature_composition_human_rows"] = int(len(comp))
        part_b_audit["mature_composition_human_cohorts"] = sorted(comp["cohort"].astype(str).unique().tolist())

    for ds in ("gao", "yazar"):
        spread, spread_msg, spread_audit = load_mature_spread(ds)
        fraction, frac_msg, frac_audit = load_mature_fractions(ds)
        part_b_audit["datasets"][ds] = {
            "spread": spread_audit,
            "fraction": frac_audit,
        }
        spread_src = str(MATURE_SPREAD_PATHS[ds].relative_to(REPO_ROOT)) if spread is not None else "ABSENT"
        frac_src = str(MATURE_FRACTION_PATHS[ds].relative_to(REPO_ROOT)) if fraction is not None else "ABSENT"

        if spread is None or fraction is None:
            reasons = [m for m in (spread_msg, frac_msg) if m != "ok"]
            if spread is None and fraction is None:
                status = "UNTESTABLE_MISSING_BOTH_MATURE_INPUTS"
            elif spread is None:
                status = "UNTESTABLE_MISSING_MATURE_SPREAD"
            else:
                status = "UNTESTABLE_MISSING_MATURE_FRACTIONS"
            rows.append(
                {
                    "dataset": ds,
                    "row_type": "cross_type_model",
                    "cell_type": "ALL_MATCHED_CELL_TYPES",
                    "status": status,
                    "untestable_reason": "; ".join(reasons),
                    "n_cell_types": 0,
                    "spread_source": spread_src,
                    "fraction_source": frac_src,
                    "weighting": "inverse_variance_1_over_spread_se_squared",
                    "family_size": 1,
                    "q_bh": np.nan,
                }
            )
            log(f"{ds}: Part B {status} ({'; '.join(reasons)})")
            continue

        joined = spread.merge(fraction, on="cell_type", how="inner", suffixes=("_spread", "_fraction"))
        for _, item in joined.sort_values("cell_type").iterrows():
            fraction_slope = float(item["effect_clr"])
            rows.append(
                {
                    "dataset": ds,
                    "row_type": "cell_type",
                    "cell_type": item["cell_type"],
                    "status": "MATCHED",
                    "untestable_reason": "",
                    "spread_slope_per_decade": float(item["slope_raw_per_decade"]),
                    "spread_se": float(item["se_raw"]),
                    "spread_pvalue": float(item["pvalue_spread"]),
                    "spread_q_bh": float(item["q_bh_within_slice"]) if pd.notna(item["q_bh_within_slice"]) else np.nan,
                    "spread_n_donors": int(item["n_donors_spread"]),
                    "fraction_slope_clr_per_decade": fraction_slope,
                    "fraction_se": float(item["se"]),
                    "fraction_pvalue": float(item["pvalue_fraction"]),
                    "fraction_q_bh": float(item["q_bh_family_a"]),
                    "fraction_n_donors": int(item["n_donors_fraction"]),
                    "fraction_direction": "losing" if fraction_slope < 0 else "gaining",
                    "spread_source": spread_src,
                    "fraction_source": frac_src,
                }
            )
        model = joined[np.isfinite(joined["slope_raw_per_decade"]) & np.isfinite(joined["effect_clr"])].copy()
        model = model[np.isfinite(model["se_raw"]) & (model["se_raw"] > 0)]
        losing = model[model["effect_clr"] < 0]
        gaining = model[model["effect_clr"] >= 0]
        base: dict[str, Any] = {
            "dataset": ds,
            "row_type": "cross_type_model",
            "cell_type": "ALL_MATCHED_CELL_TYPES",
            "status": "EVALUATED",
            "untestable_reason": "",
            "n_cell_types": int(len(model)),
            "n_fraction_losing": int(len(losing)),
            "n_fraction_gaining": int(len(gaining)),
            "mean_spread_slope_losing": float(losing["slope_raw_per_decade"].mean()) if len(losing) else np.nan,
            "mean_spread_slope_gaining": float(gaining["slope_raw_per_decade"].mean()) if len(gaining) else np.nan,
            "spread_source": spread_src,
            "fraction_source": frac_src,
            "weighting": "inverse_variance_1_over_spread_se_squared",
            "family_size": 1,
            "q_bh": np.nan,
        }
        if len(model) < 3 or model["effect_clr"].nunique() < 2:
            base["status"] = "UNTESTABLE_CELL_TYPE_FLOOR"
            base["untestable_reason"] = f"only {len(model)} matched cell types"
        else:
            y = model["slope_raw_per_decade"].to_numpy(dtype=float)
            x = model["effect_clr"].to_numpy(dtype=float)
            X = sm.add_constant(pd.DataFrame({"fraction_slope": x}), has_constant="add")
            weights = 1.0 / np.square(model["se_raw"].to_numpy(dtype=float))
            fit = sm.WLS(y, X, weights=weights).fit()
            df = int(fit.df_resid)
            if df <= 0:
                base["status"] = "UNTESTABLE_INSUFFICIENT_DF"
                base["untestable_reason"] = "zero residual df"
            else:
                slope = float(fit.params["fraction_slope"])
                se = float(fit.bse["fraction_slope"])
                t975 = float(stats.t.ppf(0.975, df))
                t80 = float(stats.t.ppf(0.80, df))
                mde = (t975 + t80) * se
                spread_sd = float(np.std(y, ddof=1))
                base.update(
                    association_slope=slope,
                    association_se=se,
                    ci95_low=slope - t975 * se,
                    ci95_high=slope + t975 * se,
                    pvalue=float(fit.pvalues["fraction_slope"]),
                    residual_df=df,
                    mde80=mde,
                    spread_slope_sd=spread_sd,
                    mde80_spread_sd_units=mde / spread_sd if spread_sd > 0 else np.nan,
                    weight_sum=float(weights.sum()),
                )
        rows.append(base)
        log(f"{ds}: mature narrowing join={len(joined)}, weighted model cell types={len(model)}")

    result = pd.DataFrame(rows)
    preferred = [
        "dataset", "row_type", "cell_type", "status", "untestable_reason",
        "association_slope", "association_se", "ci95_low", "ci95_high", "pvalue",
        "q_bh", "family_size", "mde80", "mde80_spread_sd_units", "residual_df",
        "n_cell_types", "n_fraction_losing", "n_fraction_gaining",
        "mean_spread_slope_losing", "mean_spread_slope_gaining", "spread_slope_sd",
        "weight_sum", "weighting", "spread_slope_per_decade", "spread_se",
        "spread_pvalue", "spread_q_bh", "spread_n_donors",
        "fraction_slope_clr_per_decade", "fraction_se", "fraction_pvalue",
        "fraction_q_bh", "fraction_n_donors", "fraction_direction",
        "spread_source", "fraction_source",
    ]
    for column in preferred:
        if column not in result:
            result[column] = np.nan
    result = result[preferred]
    result["_order"] = result["row_type"].map({"cell_type": 0, "cross_type_model": 1})
    result = result.sort_values(["dataset", "_order", "cell_type"]).drop(columns="_order").reset_index(drop=True)
    write_clean_csv(result, NARROWING_PATH)
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    (CACHE_ROOT / "part_b_audit.json").write_text(json.dumps(part_b_audit, indent=2) + "\n", encoding="utf-8")
    return result, part_b_audit


def part_a_verdicts(death: pd.DataFrame, removal: pd.DataFrame) -> dict[str, Any]:
    verdicts: dict[str, str] = {}
    details: dict[str, str] = {}
    for ds in ("gao", "yazar"):
        primary = death[(death["dataset"] == ds) & (death["test"] == "test1_primary_full")]
        evaluated = primary[primary["status"] == "EVALUATED"]
        hits = evaluated[evaluated["q_bh"] < 0.05]
        if len(hits):
            predicted = 0
            evaluated_removals = 0
            if not removal.empty:
                summaries = removal[(removal["dataset"] == ds) & (removal["row_type"] == "summary")]
                evaluated_removals = int((summaries["status"] == "EVALUATED").sum())
                predicted = int(summaries["n_pass_to_fail_predicted"].fillna(0).sum())
            if predicted >= 1:
                verdicts[ds] = "PROVES"
                details[ds] = (
                    f"{len(hits)} primary test 1 BH hit(s) and {predicted} predicted pass-to-fail "
                    f"gene change(s) with absolute-slope attenuation across {evaluated_removals} evaluated type(s)"
                )
            else:
                verdicts[ds] = "INCONCLUSIVE"
                details[ds] = (
                    f"{len(hits)} primary test 1 BH hit(s) but 0 predicted pass-to-fail gene changes "
                    f"with absolute-slope attenuation (evaluated removal types: {evaluated_removals})"
                )
        elif len(evaluated) and (evaluated["mde80_sd_per_decade"] <= MDE80_SMALL_SD).all():
            max_mde = float(evaluated["mde80_sd_per_decade"].max())
            verdicts[ds] = "DOES NOT PROVE"
            details[ds] = (
                f"0 primary test 1 BH hits across {len(evaluated)} evaluated types and every primary null "
                f"has MDE80 <= 0.50 SD/decade (max MDE80 = {max_mde:.3f} SD/decade)"
            )
        else:
            over = evaluated[evaluated["mde80_sd_per_decade"] > MDE80_SMALL_SD]
            max_mde = float(evaluated["mde80_sd_per_decade"].max()) if len(evaluated) else float("nan")
            over_names = ", ".join(over["cell_type"].astype(str).tolist()) if len(over) else "none"
            verdicts[ds] = "INCONCLUSIVE"
            details[ds] = (
                f"0 primary test 1 BH hits across {len(evaluated)} evaluated types, but {len(over)} primary "
                f"null(s) exceed MDE80 0.50 SD/decade (max MDE80 = {max_mde:.3f} SD/decade; over threshold: {over_names})"
            )
    if "PROVES" in (verdicts["gao"], verdicts["yazar"]):
        verdicts["overall"] = "PROVES"
    elif verdicts["gao"] == "DOES NOT PROVE" and verdicts["yazar"] == "DOES NOT PROVE":
        verdicts["overall"] = "DOES NOT PROVE"
    else:
        verdicts["overall"] = "INCONCLUSIVE"
    return {"verdicts": verdicts, "details": details}


def part_b_verdict(narrowing: pd.DataFrame) -> tuple[str, str]:
    models = narrowing[narrowing["row_type"] == "cross_type_model"]
    evaluated = models[models["status"] == "EVALUATED"]
    untestable = models[models["status"] != "EVALUATED"]
    if len(evaluated) == 0:
        reasons = "; ".join(f"{r['dataset']}: {r['status']}" for _, r in untestable.iterrows())
        return "UNTESTABLE", f"neither dataset had both verified mature inputs ({reasons})"
    if len(evaluated) < 2:
        ev = evaluated.iloc[0]
        un = untestable.iloc[0]
        return (
            "INCONCLUSIVE",
            (
                f"{ev['dataset']} evaluated (slope {float(ev['association_slope']):+.5f}, "
                f"p = {float(ev['pvalue']):.4f}), while {un['dataset']} is {un['status']} "
                f"({un['untestable_reason']}); both datasets are required to establish DOES NOT PROVE "
                f"or verify non-opposition for PROVES"
            ),
        )
    positive = (evaluated["association_slope"] > 0) & (evaluated["pvalue"] < 0.05)
    opposite = (evaluated["association_slope"] < 0) & (evaluated["pvalue"] < 0.05)
    if positive.any() and not opposite.any():
        return "PROVES", "positive significant cross-type association in at least one dataset and not significantly opposite in the other"
    if positive.any() and opposite.any():
        return "INCONCLUSIVE", "significantly positive in one dataset and significantly opposite in the other"
    null_or_opposite = (evaluated["pvalue"] >= 0.05) | (evaluated["association_slope"] < 0)
    if null_or_opposite.all():
        return "DOES NOT PROVE", "association is null or opposite in both mature datasets"
    return "INCONCLUSIVE", "mixed or non-evaluable Part B pattern"


def fmt(value: Any, digits: int = 4, signed: bool = False) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):+.{digits}f}" if signed else f"{float(value):.{digits}f}"


def fmt_p(value: Any) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.3e}"


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def generate_report() -> dict[str, Any]:
    death = pd.read_csv(DEATH_PATH)
    removal = pd.read_csv(REMOVAL_PATH)
    narrowing = pd.read_csv(NARROWING_PATH)
    orig_death = pd.read_csv(ORIG_DEATH_PATH)
    orig_removal = pd.read_csv(ORIG_REMOVAL_PATH)
    orig_narrowing = pd.read_csv(ORIG_NARROWING_PATH)

    score_audits = {
        ds: json.loads((CACHE_ROOT / ds / "scores" / "SCORES_DONE.json").read_text(encoding="utf-8"))
        for ds in ("gao", "yazar")
    }
    bh_audit = load_bh_audit()
    part_b_audit = json.loads((CACHE_ROOT / "part_b_audit.json").read_text(encoding="utf-8"))

    part_a_res = part_a_verdicts(death, removal)
    verdict_a = part_a_res["verdicts"]
    detail_a = part_a_res["details"]
    verdict_b, detail_b = part_b_verdict(narrowing)

    apoptosis = ", ".join(APOPTOSIS)
    upr = ", ".join(UPR)
    union = ", ".join(SCORING_GENES)

    lines: list[str] = [
        "# Mature-to-Old (30-plus) Block 7: Survivor-Bias Within-Type Tests",
        "",
        "## Status and scope",
        "",
        "Phase 2 Item 13 of the Agescape mature-to-old rebuild. This report reruns 11pt Block 7 strictly on human mature donors aged 30 years and older (`age_years >= 30`), reusing the frozen protocol, frozen 271-gene death/stress scoring list, and frozen verdict rules from `reports/followup-11pt-block7-survivor-within-type.md` verbatim. Mature signature P1m does not enter Block 7 because Block 7 scores the frozen Hallmark apoptosis plus unfolded protein response gene list.",
        "",
        "## Frozen protocol as applied at 30-plus",
        "",
        "### Frozen death/stress gene list and citation",
        "",
        "Source: Liberzon et al. 2015, Cell Systems, The Molecular Signatures Database Hallmark Gene Set Collection, DOI `10.1016/j.cels.2015.12.004`, PMID `26771021`. Official MSigDB sets: `HALLMARK_APOPTOSIS` (systematic ID `M5902`, 161 genes) and `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` (systematic ID `M5922`, 113 genes). Their union is 272 unique symbols (`ATF3` and `DNAJC3` occur in both). The shared human technical mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. `RPS14` is the sole excluded source-list gene, leaving 271 fixed scoring genes. Absent dataset symbols contribute zero counts and remain in the fixed denominator of 271 (Gao has 248 present and 23 absent; Yazar has 271 present and 0 absent).",
        "",
        f"Apoptosis source genes (161): {apoptosis}.",
        "",
        f"Unfolded protein response source genes (113): {upr}.",
        "",
        f"Final scoring symbols (271): {union}.",
        "",
        "### Part A specification (restricted to age 30-plus)",
        "",
        "- Per-cell program score is the mean within-cell proportion across the fixed 271 scoring genes: sum of post-mask raw counts divided by `271 * post-mask UMI`. Primary analyses use full-depth counts without downsampling. Cells with post-mask UMI below 1500 are excluded as an eligibility filter only.",
        "- All steps are restricted strictly to donors aged 30 years and older (`age_years >= 30.0`): cell eligibility, donor-cell eligibility (at least 20 eligible cells per mature donor per cell type), the per-cell-type 90th-percentile high-score threshold (pooled across eligible mature donor-cell-type pairs only, numpy linear quantile interpolation, high defined as `score >= threshold` including ties), and every regression. No full-range or young-derived threshold is reused.",
        "- Test 1 outcome is each mature donor's fraction of high-scoring cells among that donor's eligible cells in the cell type. Primary OLS fits `fraction_high ~ age_per_decade + sex_male + batch` (Gao 7-study fixed effects with `ENCODE v4 (Snyder)` reference; Yazar 75-pool fixed effects) requiring at least 6 mature donors, singleton-batch pruning (`min_batch_donors = 2`), rank-checked covariate pruning, and residual df at least 5 (Gao) or 10 (Yazar). Unadjusted sensitivity retains `age_per_decade + sex_male` without batch (never substituted for a failed batch-adjusted fit). Models without sufficient residual df or donor/batch floors are recorded as UNTESTABLE, never null.",
        "- Labeled 1500 sensitivity repeats test 1 only using exact marginal `Binomial(original_list_count, 1500 / original_UMI)` draws from `SeedSequence([20261117, 7, dataset_code, donor_idx])` with mature cells sorted by original h5ad row index, and recomputes the 90th-percentile threshold on mature rarefied scores only.",
        "- Test 1 BH families are each `(dataset, test)` slice across evaluable non-pooled cell types at 30-plus, separately for `test1_primary_full`, `test1_unadjusted_full`, and `test1_rarefied_1500`. Every BH family was independently recomputed using a manual step-up implementation and cross-checked against `statsmodels.stats.multitest.multipletests`.",
        "- Test 2 runs only for cell types passing primary test 1 BH `q < 0.05` at 30-plus. Unremoved and removed pseudobulk gene OLS models (`log2CPM ~ age_per_decade + sex_male + log_lib + batch`) use the same common mature donor set (at least 20 cells per donor in both arms, singleton batches pruned, at least 6 common mature donors), gene expression filter `CPM >= 1` in at least half of donors (minimum 3 donors), and within-arm gene BH (`q_bh_within_group < 0.05`). A predicted pass-to-fail change requires an unremoved pass becoming a removed fail with `|slope_removed| < |slope_unremoved|`. Slope change is `removed - unremoved`. Removing captured high-scoring cells is a lower bound because dying cells are under-captured.",
        "- MDE80 is `(t_0.975 + t_0.80) * SE` using residual df, reported on both the raw outcome scale and standardized by the donor-outcome SD (`mde80_sd_per_decade`). Small MDE80 is predefined as at most `0.50` outcome SD per decade.",
        "",
        "### Part B specification (mature slopes only)",
        "",
        "- Joins ONLY mature-recomputed per-cell-type batch-adjusted CLR fraction slopes per decade and mature-recomputed full-depth primary Shannon spread slopes per decade (`metric=shannon`, `depth_mode=full_depth`, `analysis_role=primary_full_depth`, `model_spec=primary`, non-pooled, non-stratified-global). Legacy full-range slopes are never joined.",
        "- Fraction slope negative means losing share; spread slope is Shannon change per decade (negative means narrowing). Frozen prediction: survivor bias predicts a positive cross-type association between fraction slope and spread slope, with more negative spread slopes in fraction-losing cell types.",
        "- Weighted least squares (`spread_slope ~ 1 + fraction_slope`) with inverse-variance weights `1 / SE_spread^2`, fitted separately per dataset, one test per dataset, no multiplicity correction across the two datasets, never pooled. Zero fraction slope counts as gaining.",
        "",
        "### Frozen verdict rules",
        "",
        "- Part A per dataset: `PROVES` only if test 1 passes BH correction (`q < 0.05`) and test 2 produces at least one predicted pass-to-fail gene change with absolute-slope attenuation. `DOES NOT PROVE` if no test 1 result passes BH and every evaluable primary null has MDE80 at most 0.50 SD per decade. `INCONCLUSIVE` if any primary null has MDE80 above 0.50 SD per decade, or test 1 is significant but removal produces zero predicted changes, or the required removal comparison is not evaluable. Overall Part A is `PROVES` if either dataset meets `PROVES`, `DOES NOT PROVE` only if both meet `DOES NOT PROVE`, and otherwise `INCONCLUSIVE`.",
        "- Part B: `PROVES` only if the cross-type association is positive and significant (`p < 0.05`) in at least one dataset and the other dataset is not significantly opposite. `DOES NOT PROVE` if the association is null or opposite in both datasets. `INCONCLUSIVE` if one is significantly positive and the other significantly opposite, or if a dataset is UNTESTABLE.",
        "- Neither Part A nor Part B alone proves survivor bias as a mechanism. The Nov 17 cross-type correlation is not rerun here.",
        "",
        "## Input provenance and mature cohort audit",
        "",
    ]

    prov_rows = [
        [
            "Gao counts cache (Block 1)",
            "`/home/reza/projects/agescape-revision-data/tmp_organs/block1_cache/gao`",
            "Reused Block 1 donor-sharded post-mask UMI >= 1500 cache",
            f"{score_audits['gao']['n_donors_mature']} of {score_audits['gao']['n_donors_original']} donors >= 30 "
            f"({score_audits['gao']['n_donors_under30_excluded']} excluded; ages {score_audits['gao']['age_min_mature']:.0f} to {score_audits['gao']['age_max_mature']:.0f}; "
            f"{score_audits['gao']['n_mature_cells_scored']} mature cells scored, {score_audits['gao']['n_under30_cells_excluded']} under-30 cells excluded)",
        ],
        [
            "Yazar counts cache (Block 1)",
            "`/home/reza/projects/agescape-revision-data/tmp_organs/block1_cache/yazar`",
            "Reused Block 1 donor-sharded post-mask UMI >= 1500 cache",
            f"{score_audits['yazar']['n_donors_mature']} of {score_audits['yazar']['n_donors_original']} donors >= 30 "
            f"({score_audits['yazar']['n_donors_under30_excluded']} excluded; ages {score_audits['yazar']['age_min_mature']:.0f} to {score_audits['yazar']['age_max_mature']:.0f}; "
            f"{score_audits['yazar']['n_mature_cells_scored']} mature cells scored, {score_audits['yazar']['n_under30_cells_excluded']} under-30 cells excluded)",
        ],
        [
            "Core mature human composition",
            "`results/followup-mature-composition-human.csv`",
            f"{part_b_audit.get('mature_composition_human_provenance', 'missing')}; sha256 `{str(part_b_audit.get('mature_composition_human_sha256'))[:12]}`",
            f"Checked first: {part_b_audit.get('mature_composition_human_rows', 0)} rows across 11 core brain/kidney/lung cohorts; 0 Gao and 0 Yazar rows (not joined)",
        ],
    ]
    for ds in ("gao", "yazar"):
        ds_label = "Gao heart" if ds == "gao" else "Yazar blood"
        sp_a = part_b_audit["datasets"][ds]["spread"]
        fr_a = part_b_audit["datasets"][ds]["fraction"]
        prov_rows.append(
            [
                f"{ds_label} mature spread",
                f"`{sp_a['path']}`",
                f"{sp_a.get('provenance', 'ABSENT')}; sha256 `{str(sp_a.get('sha256', 'NA'))[:12]}`" if sp_a["exists"] else "ABSENT at runtime",
                f"Verified mature primary full-depth Shannon rows: {sp_a.get('n_rows_filtered', 0)}, max N donors = {sp_a.get('max_n_donors', 0)}" if sp_a["exists"] else "Missing; dataset marked UNTESTABLE",
            ]
        )
        prov_rows.append(
            [
                f"{ds_label} mature CLR fractions",
                f"`{fr_a['path']}`",
                f"{fr_a.get('provenance', 'ABSENT')}; sha256 `{str(fr_a.get('sha256', 'NA'))[:12]}`" if fr_a["exists"] else "ABSENT at runtime",
                (
                    f"Verified mature `fractions_clr` (`arm={fr_a.get('arm')}`): {fr_a.get('n_rows_filtered', 0)} cell types, "
                    f"N = {fr_a.get('n_donors')} mature donors ({fr_a.get('n_donors_30plus_excluded')} under-30 excluded)"
                ) if fr_a["exists"] else "Missing; dataset marked UNTESTABLE",
            ]
        )
    lines.append(markdown_table(["Input", "Path", "Provenance / SHA-256 prefix", "Verification"], prov_rows))
    lines.extend(["", "## Results", ""])

    # Part A Test 1
    lines.extend(["### Part A test 1: high death/stress-score cell fractions at age 30-plus", ""])
    orig_primary_map = {
        (str(r["dataset"]), str(r["cell_type"])): r
        for _, r in orig_death[orig_death["test"] == "test1_primary_full"].iterrows()
    }
    for ds in ("gao", "yazar"):
        primary = death[(death["dataset"] == ds) & (death["test"] == "test1_primary_full")].sort_values("cell_type")
        evaluated = primary[primary["status"] == "EVALUATED"]
        untestable = primary[primary["status"] != "EVALUATED"]
        hits = evaluated[evaluated["q_bh"] < 0.05]
        over_mde = evaluated[evaluated["mde80_sd_per_decade"] > MDE80_SMALL_SD]
        lines.extend(
            [
                f"#### {ds.capitalize()} (primary full depth, batch-adjusted)",
                "",
                (
                    f"Evaluated {len(evaluated)} of {len(primary)} cell types ({len(untestable)} UNTESTABLE). "
                    f"BH family size {int(evaluated['family_size'].max()) if len(evaluated) else 0}; "
                    f"{len(hits)} passed BH q < 0.05. "
                    f"Evaluated primary nulls with MDE80 > 0.50 SD/decade: {len(over_mde)} "
                    f"(max MDE80 = {float(evaluated['mde80_sd_per_decade'].max()):.3f} SD/decade)."
                ),
                "",
                markdown_table(
                    ["Cell type", "Status", "N orig", "N mature", "Resid df", "Slope per decade", "95 percent CI", "p", "BH q", "MDE80 SD"],
                    [
                        [
                            row["cell_type"],
                            row["status"],
                            int(orig_primary_map[(ds, str(row["cell_type"]))]["n_donors"])
                            if (ds, str(row["cell_type"])) in orig_primary_map
                            else "NA",
                            int(row["n_donors"]),
                            int(row["residual_df"]) if pd.notna(row["residual_df"]) else "NA",
                            fmt(row["slope_per_decade"], 6, True),
                            f"{fmt(row['ci95_low'], 6)} to {fmt(row['ci95_high'], 6)}",
                            fmt_p(row["pvalue"]),
                            fmt_p(row["q_bh"]),
                            fmt(row["mde80_sd_per_decade"], 3),
                        ]
                        for _, row in primary.iterrows()
                    ],
                ),
                "",
            ]
        )
        if not untestable.empty:
            lines.append(f"UNTESTABLE primary cell types in {ds.capitalize()}:")
            lines.append("")
            lines.append(
                markdown_table(
                    ["Cell type", "Status / Reason", "N mature donors (>= 20 cells)", "N mature donors with >= 1 cell", "Covariates omitted"],
                    [
                        [
                            row["cell_type"],
                            row["status"],
                            int(row["n_donors"]),
                            int(row["n_donors_with_cells"]),
                            str(row["covariates_omitted"]) if pd.notna(row["covariates_omitted"]) else "none",
                        ]
                        for _, row in untestable.iterrows()
                    ],
                )
            )
            lines.append("")

        sensitivity = death[(death["dataset"] == ds) & (death["test"] != "test1_primary_full")]
        for test_name, group in sensitivity.groupby("test", sort=True):
            ev_s = group[group["status"] == "EVALUATED"]
            hits_s = ev_s[ev_s["q_bh"] < 0.05]
            hit_list = ", ".join(
                f"{r['cell_type']} (slope {float(r['slope_per_decade']):+.6f}, p={float(r['pvalue']):.3e}, q={float(r['q_bh']):.3e})"
                for _, r in hits_s.sort_values("q_bh").iterrows()
            )
            lines.append(
                f"- `{ds} {test_name}`: {len(ev_s)} evaluated of {len(group)}, {len(hits_s)} BH hit(s)"
                + (f": {hit_list}." if hit_list else ".")
            )
        lines.append("")

    # Independent BH recomputation check
    lines.extend(
        [
            "### Independent BH recomputation verification",
            "",
            "Every BH family in Part A test 1 and test 2 was independently recomputed from raw p-values using a standalone step-up implementation (`p_sorted * m / rank`, reverse cumulative minimum, clipped to `[0, 1]`) and checked against `statsmodels.stats.multitest.multipletests(..., method='fdr_bh')`.",
            "",
            markdown_table(
                ["Stage", "Dataset", "BH family", "Family size", "Passing q < 0.05", "Min raw p", "Min BH q", "Max abs diff"],
                [
                    [
                        e["stage"],
                        e["dataset"],
                        e["family"],
                        int(e["family_size"]),
                        int(e["n_q_lt_0_05"]),
                        fmt_p(e["min_p"]),
                        fmt_p(e["min_q"]),
                        f"{float(e['max_abs_diff_statsmodels_vs_manual']):.3e}",
                    ]
                    for e in sorted(bh_audit, key=lambda x: (x["stage"], x["dataset"], x["family"]))
                ],
            ),
            "",
        ]
    )

    # Part A Test 2
    lines.extend(["### Part A test 2: removal sensitivity at age 30-plus", ""])
    summaries = removal[removal["row_type"] == "summary"] if not removal.empty else pd.DataFrame()
    if summaries.empty:
        lines.append(
            "No cell type passed primary test 1 BH correction (`q < 0.05`) at age 30-plus, so test 2 was not run. "
            "`results/followup-mature-block7-removal-sensitivity.csv` contains the full schema header plus one explicit "
            "`no_eligible_type` record row (`status=NO_TEST1_BH_HIT`)."
        )
    else:
        lines.append(
            markdown_table(
                [
                    "Dataset", "Cell type", "Status", "Common mature donors",
                    "Unremoved cells", "Removed cells", "High cells removed",
                    "Common genes", "Unremoved pass", "Removed pass",
                    "Predicted pass to fail", "Other pass to fail", "Fail to pass",
                    "Pass to pass", "Fail to fail",
                ],
                [
                    [
                        row["dataset"],
                        row["cell_type"],
                        row["status"],
                        int(row["n_donors"]),
                        int(row["n_cells_unremoved"]),
                        int(row["n_cells_removed"]),
                        int(row["n_high_cells_removed"]),
                        int(row["n_common_genes"]),
                        int(row["n_pass_unremoved"]),
                        int(row["n_pass_removed"]),
                        int(row["n_pass_to_fail_predicted"]),
                        int(row["n_pass_to_fail_not_attenuated"]),
                        int(row["n_fail_to_pass"]),
                        int(row["n_pass_to_pass"]),
                        int(row["n_fail_to_fail"]),
                    ]
                    for _, row in summaries.iterrows()
                ],
            )
        )
        lines.append("")
        # Also note median absolute slope shift among unremoved-passing genes
        for _, srow in summaries.iterrows():
            ds_s = str(srow["dataset"])
            ct_s = str(srow["cell_type"])
            g_sub = removal[
                (removal["dataset"] == ds_s)
                & (removal["cell_type"] == ct_s)
                & (removal["row_type"] == "gene")
                & (removal["status"] == "EVALUATED_BOTH")
            ]
            if not g_sub.empty:
                pred_g = g_sub[g_sub["pass_fail_change"] == "PASS_TO_FAIL_PREDICTED"].copy()
                if not pred_g.empty:
                    pred_g["abs_atten"] = pred_g["slope_unremoved"].abs() - pred_g["slope_removed"].abs()
                    top_genes = ", ".join(pred_g.sort_values("abs_atten", ascending=False)["gene"].head(5).tolist())
                    med_atten = float(pred_g["abs_atten"].median())
                    lines.append(
                        f"- {ds_s} `{ct_s}`: {len(pred_g)} predicted pass-to-fail genes with absolute-slope attenuation "
                        f"(median |slope| reduction = {med_atten:.4f} log2CPM/decade; top attenuated genes: {top_genes})."
                    )
    lines.extend(
        [
            "",
            "Pass/fail change counts use genes evaluated in both arms (`EVALUATED_BOTH`) on the exact same common mature donor set. Slope change is `removed - unremoved`. Caveat: removing captured high-scoring cells provides a lower bound on stress/dying-cell contribution because dying cells are under-captured in droplet scRNA-seq/snRNA-seq.",
            "",
        ]
    )

    # Part B
    lines.extend(["### Part B: narrowing localization (mature slopes only)", ""])
    models_b = narrowing[narrowing["row_type"] == "cross_type_model"].sort_values("dataset")
    lines.append(
        markdown_table(
            [
                "Dataset", "Status", "N types", "Association slope", "95 percent CI",
                "p", "MDE80", "Losing N", "Gaining N", "Mean spread losing", "Mean spread gaining", "Reason / Notes",
            ],
            [
                [
                    row["dataset"],
                    row["status"],
                    int(row["n_cell_types"]) if pd.notna(row["n_cell_types"]) else 0,
                    fmt(row["association_slope"], 5, True),
                    f"{fmt(row['ci95_low'], 5)} to {fmt(row['ci95_high'], 5)}",
                    fmt_p(row["pvalue"]),
                    fmt(row["mde80"], 5),
                    int(row["n_fraction_losing"]) if pd.notna(row["n_fraction_losing"]) else "NA",
                    int(row["n_fraction_gaining"]) if pd.notna(row["n_fraction_gaining"]) else "NA",
                    fmt(row["mean_spread_slope_losing"], 5, True),
                    fmt(row["mean_spread_slope_gaining"], 5, True),
                    str(row["untestable_reason"]) if pd.notna(row["untestable_reason"]) and str(row["untestable_reason"]) else "evaluated on mature slopes",
                ]
                for _, row in models_b.iterrows()
            ],
        )
    )
    lines.extend(
        [
            "",
            "The association slope is the change in mature Shannon spread slope per unit mature CLR fraction slope per decade (WLS with weights `1 / SE_spread^2`). Fraction slope negative means losing share; zero fraction slope counts as gaining. Frozen prediction: positive association. One test per dataset, no multiplicity correction across datasets, never pooled.",
            "",
        ]
    )
    matched_b = narrowing[narrowing["row_type"] == "cell_type"].sort_values(["dataset", "cell_type"])
    if not matched_b.empty:
        lines.extend(["#### Part B matched per-cell-type mature slopes", ""])
        lines.append(
            markdown_table(
                [
                    "Dataset", "Cell type", "Spread N", "Mature spread slope", "Spread SE",
                    "Fraction N", "Mature CLR fraction slope", "Fraction SE", "Direction",
                ],
                [
                    [
                        r["dataset"],
                        r["cell_type"],
                        int(r["spread_n_donors"]),
                        fmt(r["spread_slope_per_decade"], 5, True),
                        fmt(r["spread_se"], 5),
                        int(r["fraction_n_donors"]),
                        fmt(r["fraction_slope_clr_per_decade"], 5, True),
                        fmt(r["fraction_se"], 5),
                        r["fraction_direction"],
                    ]
                    for _, r in matched_b.iterrows()
                ],
            )
        )
        lines.append("")

    # Verdicts
    lines.extend(
        [
            "## Verdicts under the frozen rules",
            "",
            f"- Part A Gao: **{verdict_a['gao']}** ({detail_a['gao']}).",
            f"- Part A Yazar: **{verdict_a['yazar']}** ({detail_a['yazar']}).",
            f"- Part A overall: **{verdict_a['overall']}**.",
            f"- Part B: **{verdict_b}** ({detail_b}).",
            "",
            "Neither Part A nor Part B alone proves survivor bias as a mechanism. These tests evaluate within-type high-stress cell dependence and cross-type narrowing localization only. The Nov 17 cross-type correlation was not rerun here.",
            "",
        ]
    )

    # Comparison vs original Block 7
    lines.extend(["## Comparison vs original full-range Block 7", ""])
    orig_verdict_a = part_a_verdicts(orig_death, orig_removal)["verdicts"]
    orig_verdict_b, _ = part_b_verdict(orig_narrowing)
    orig_hits_y = orig_death[
        (orig_death["dataset"] == "yazar")
        & (orig_death["test"] == "test1_primary_full")
        & (orig_death["status"] == "EVALUATED")
        & (orig_death["q_bh"] < 0.05)
    ].sort_values("cell_type")
    mat_y_primary = death[
        (death["dataset"] == "yazar") & (death["test"] == "test1_primary_full")
    ].set_index("cell_type")
    orig_rem_sum = orig_removal[orig_removal["row_type"] == "summary"].set_index(["dataset", "cell_type"])
    mat_rem_sum = (
        summaries.set_index(["dataset", "cell_type"])
        if not summaries.empty
        else pd.DataFrame()
    )

    comp_verdict_rows = [
        ["Part A Gao", orig_verdict_a["gao"], verdict_a["gao"], detail_a["gao"]],
        ["Part A Yazar", orig_verdict_a["yazar"], verdict_a["yazar"], detail_a["yazar"]],
        ["Part A overall", orig_verdict_a["overall"], verdict_a["overall"], "PROVES if either dataset PROVES"],
        ["Part B", orig_verdict_b, verdict_b, detail_b],
    ]
    lines.append(markdown_table(["Component", "Original full-range", "Mature 30-plus", "Mature summary"], comp_verdict_rows))
    lines.append("")

    # Per-hit survival table
    lines.append("### Original test 1 BH hits and any new mature test 1 BH hits")
    lines.append("")
    mat_hits_all = death[
        (death["test"] == "test1_primary_full")
        & (death["status"] == "EVALUATED")
        & (death["q_bh"] < 0.05)
    ]
    hit_keys = sorted(
        set((str(r["dataset"]), str(r["cell_type"])) for _, r in orig_hits_y.iterrows())
        | set((str(r["dataset"]), str(r["cell_type"])) for _, r in mat_hits_all.iterrows())
    )
    hit_comp_rows: list[list[Any]] = []
    for ds_h, ct_h in hit_keys:
        orow = orig_primary_map[(ds_h, ct_h)]
        mrow = death[
            (death["dataset"] == ds_h)
            & (death["test"] == "test1_primary_full")
            & (death["cell_type"] == ct_h)
        ].iloc[0]
        orig_pass = bool(pd.notna(orow["q_bh"]) and float(orow["q_bh"]) < 0.05)
        mat_pass = bool(pd.notna(mrow["q_bh"]) and float(mrow["q_bh"]) < 0.05)
        if orig_pass and mat_pass:
            fate = "SURVIVED"
        elif orig_pass and not mat_pass:
            fate = "DIED"
        elif not orig_pass and mat_pass:
            fate = "NEW_AT_30_PLUS"
        else:
            fate = "NULL_BOTH"
        orig_pred = (
            int(orig_rem_sum.loc[(ds_h, ct_h), "n_pass_to_fail_predicted"])
            if (ds_h, ct_h) in orig_rem_sum.index
            else "NA"
        )
        mat_pred = (
            int(mat_rem_sum.loc[(ds_h, ct_h), "n_pass_to_fail_predicted"])
            if not mat_rem_sum.empty and (ds_h, ct_h) in mat_rem_sum.index
            else "NA"
        )
        hit_comp_rows.append(
            [
                ds_h,
                ct_h,
                fate,
                int(orow["n_donors"]),
                int(mrow["n_donors"]),
                fmt(orow["slope_per_decade"], 6, True),
                fmt_p(orow["q_bh"]),
                fmt(mrow["slope_per_decade"], 6, True),
                fmt_p(mrow["q_bh"]),
                orig_pred,
                mat_pred,
            ]
        )
    lines.append(
        markdown_table(
            [
                "Dataset", "Cell type", "Status vs original", "N orig", "N mature",
                "Orig slope", "Orig BH q", "Mature slope", "Mature BH q",
                "Orig predicted pass-to-fail", "Mature predicted pass-to-fail",
            ],
            hit_comp_rows,
        )
    )
    lines.append("")

    # Part B comparison table
    orig_b_map = {
        str(r["dataset"]): r
        for _, r in orig_narrowing[orig_narrowing["row_type"] == "cross_type_model"].iterrows()
    }
    part_b_comp_rows = []
    for _, mr in models_b.iterrows():
        ds_b = str(mr["dataset"])
        or_b = orig_b_map[ds_b]
        part_b_comp_rows.append(
            [
                ds_b,
                int(or_b["n_cell_types"]),
                fmt(or_b["association_slope"], 5, True),
                fmt_p(or_b["pvalue"]),
                mr["status"],
                int(mr["n_cell_types"]) if pd.notna(mr["n_cell_types"]) else 0,
                fmt(mr["association_slope"], 5, True),
                fmt_p(mr["pvalue"]),
            ]
        )
    lines.extend(["### Part B comparison vs original full-range Block 7", ""])
    lines.append(
        markdown_table(
            ["Dataset", "Orig N types", "Orig slope", "Orig p", "Mature status", "Mature N types", "Mature slope", "Mature p"],
            part_b_comp_rows,
        )
    )
    lines.extend(
        [
            "",
            "## Caveats and power notes",
            "",
            "- Removing captured high-scoring cells in Part A test 2 is a lower bound because dying or severely stressed cells are under-captured in single-cell and single-nucleus RNA-seq.",
            "- Any primary test 1 null with MDE80 exceeding 0.50 outcome SD per decade is treated as underpowered (`INCONCLUSIVE`), not as a powered null.",
            "- Neither Part A nor Part B alone establishes survivor bias as a causal mechanism.",
            "- All randomness uses fixed master seed `20261117` and deterministic per-donor derived seeds (`SeedSequence([20261117, 7, dataset_code, donor_idx])`). Donor is the regression unit in every Part A model, and cell type is the regression unit in Part B.",
            "",
        ]
    )

    write_clean_markdown("\n".join(lines), REPORT_PATH)
    log(f"Report written to {REPORT_PATH}: Part A {verdict_a['overall']} (Gao {verdict_a['gao']}, Yazar {verdict_a['yazar']}); Part B {verdict_b}")
    return {
        "part_a_gao": verdict_a["gao"],
        "part_a_yazar": verdict_a["yazar"],
        "part_a_overall": verdict_a["overall"],
        "part_b": verdict_b,
    }


def validate_outputs() -> None:
    for path in (DEATH_PATH, REMOVAL_PATH, NARROWING_PATH, REPORT_PATH, Path(__file__)):
        if not path.exists():
            raise FileNotFoundError(path)
        text = path.read_text(encoding="utf-8")
        if re.search(r"[\u2013\u2014]", text):
            raise ValueError(f"Forbidden em/en dash in {path}")

    death = pd.read_csv(DEATH_PATH)
    evaluated = death[death["status"] == "EVALUATED"]
    if not np.isfinite(evaluated[["slope_per_decade", "se", "pvalue", "q_bh", "mde80"]].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite evaluated test 1 result")
    if not evaluated["pvalue"].between(0, 1).all() or not evaluated["q_bh"].between(0, 1).all():
        raise ValueError("Invalid test 1 p or q value")
    if (evaluated["n_donors"] < MIN_DONORS).any():
        raise ValueError("Evaluated test 1 row below donor floor")
    if (evaluated["min_donor_age"] < MIN_AGE_YEARS).any():
        raise ValueError("Under-30 donor detected in mature test 1 output")

    # Independent BH re-check directly from saved CSV
    for (ds, test), grp in evaluated.groupby(["dataset", "test"]):
        q_check, diff = compute_bh_checked(grp["pvalue"].to_numpy(dtype=float))
        csv_diff = float(np.max(np.abs(q_check - grp["q_bh"].to_numpy(dtype=float))))
        if max(diff, csv_diff) > 1e-12:
            raise ValueError(f"Saved CSV BH check failed for {ds} {test}: diff={max(diff, csv_diff):.3e}")

    narrowing = pd.read_csv(NARROWING_PATH)
    models = narrowing[narrowing["row_type"] == "cross_type_model"]
    if len(models) != 2:
        raise ValueError(f"Expected 2 cross_type_model rows in narrowing CSV, found {len(models)}")
    for _, row in narrowing.iterrows():
        for col in ("spread_source", "fraction_source"):
            val = str(row[col])
            if any(
                legacy in val
                for legacy in (
                    "followup-fulldepth-",
                    "followup-organ-17-heart-gao-fractions",
                    "followup-hcata-pool-sensitivity-fractions",
                    "reviewer-phase1-composition",
                )
            ):
                raise ValueError(f"Legacy full-range slope source detected in Part B: {val}")

    removal = pd.read_csv(REMOVAL_PATH)
    if not removal.empty:
        common = removal[removal["status"] == "EVALUATED_BOTH"]
        if len(common) and not np.isfinite(
            common[["slope_unremoved", "slope_removed", "slope_change_removed_minus_unremoved"]].to_numpy(dtype=float)
        ).all():
            raise ValueError("Non-finite paired removal result")

    report_text = REPORT_PATH.read_text(encoding="utf-8")
    if "Final scoring symbols (271):" not in report_text or "## Verdicts under the frozen rules" not in report_text:
        raise ValueError("Report missing frozen list or verdicts section")
    log("Output validation passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=["all", "part_a", "scores", "test1", "test2", "narrowing", "report", "validate"],
        default="all",
    )
    parser.add_argument("--dataset", choices=["both", "gao", "yazar"], default="both")
    parser.add_argument("--force-rescore", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    datasets = ("gao", "yazar") if args.dataset == "both" else (args.dataset,)
    if args.stage == "all":
        wait_for_memory("Mature Block 7 score cache")
        for ds in datasets:
            run_scores(ds, force=args.force_rescore)
        death = run_test1()
        run_test2(death)
        run_narrowing()
        generate_report()
        validate_outputs()
        return
    if args.stage == "part_a":
        wait_for_memory("Mature Block 7 Part A")
        for ds in datasets:
            run_scores(ds, force=args.force_rescore)
        death = run_test1()
        run_test2(death)
        return
    if args.stage == "scores":
        wait_for_memory("Mature Block 7 score cache")
        for ds in datasets:
            run_scores(ds, force=args.force_rescore)
    elif args.stage == "test1":
        run_test1()
    elif args.stage == "test2":
        run_test2()
    elif args.stage == "narrowing":
        run_narrowing()
    elif args.stage == "report":
        generate_report()
    elif args.stage == "validate":
        validate_outputs()


if __name__ == "__main__":
    main()
