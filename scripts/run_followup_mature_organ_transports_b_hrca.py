#!/usr/bin/env python3
"""Phase 2 Audit Item 5 Batch B (HRCA): mature (30-plus) follow-up retina transports.

Cohort: Human Retina Cell Atlas (HRCA), Li et al. Nat Genet 2026.
Nine suites, mirroring scripts/run_followup_organ_10_eye_hrca.py exactly:
  snRNA (rerun at 30-plus): Sn_Adult_MacularOnly, Sn_Adult_PeripheralOnly,
    Sn_Adult_Unenriched, Sn_Adult_ExactAge, Sn_ChenAncestry_Adult,
    Sn_AllAges_Unenriched (mature rerun; its 30-plus donor set is identical
    to Sn_Adult_Unenriched by construction, so the two mature runs are an
    empirical cross-check of each other).
  scRNA (no age change; all donors already 42-plus): Sc15_StudyAdj,
    Sc9_External_StudyAdj, Sc5_External_NoDisease.
Every suite is rerun with mature P1m and P0m (Item 1, commit d3dc14c).

Reuses the frozen Batch A protocol: pseudobulk log2(CPM+1), technical-gene
removal, CPM>=1 in >=50% mature donors, original covariate models,
2000-draw expression-decile-matched null with seed 20261105, 1000-perm
concordance with seed 20261105. Donors aged >= 30.0 are restricted BEFORE
pseudobulk aggregation and detection filtering (sn suites). Donor is the
inferential unit, never the cell.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time
from typing import Any

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_mature_organ_transports_a import (  # noqa: E402
    ALPHA,
    CONCORDANCE_PERMS,
    MIN_AGE,
    ORIG_NULL_SEED,
    P1_ITERATIONS,
    bh_checked,
    load_verified_mature_programs,
    make_row,
    mde80_from_null_sd,
    mde80_from_se,
    order_cohort_df,
    peak_rss_mb,
    run_family_b_mature,
    run_family_c_mature,
    transport_class,
    utc_now,
    verdict_standard,
)
from scripts.run_followup_organ_10_eye_hrca import (  # noqa: E402
    CLASSES,
    SC_PATH,
    SN_PATH,
    _aggregate,
    _donor_meta,
)
from src.followup_organs import (  # noqa: E402
    assert_no_forbidden_dashes,
    finalize_pseudobulk_group,
    run_family_a_clr,
    sanitize_dashes,
    write_clean_csv,
)

RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
TMP_WORK_DIR = DATA_ROOT / "tmp_organs" / "hrca_mature_batch_b_tmp"

PREFIX = "followup-mature-organ-transports-b-hrca"
COHORT_CSV = RESULTS_DIR / f"{PREFIX}.csv"
BH_CHECK_CSV = RESULTS_DIR / f"{PREFIX}-bh-check.csv"
RUNTIME_CSV = RESULTS_DIR / f"{PREFIX}-runtime.csv"

MODEL_USED = "Muse Spark (direct)"
MIN_AVAIL_RAM_GIB = 4.0
RAM_WAIT_TIMEOUT_SEC = 900

EXPECTED_COLS = [
    "cohort", "cell_type", "endpoint_family", "n_donors", "n_donors_30plus_excluded",
    "original_p", "original_q", "mature_p", "mature_q", "mature_empirical_p",
    "mature_empirical_q", "verdict", "notes", "arm", "residual_df",
    "original_effect", "mature_effect", "mature_se", "null_sd", "mde80",
    "mde80_method", "covariates", "omitted_covariates", "n_genes_tested",
    "n_sig_q05", "original_empirical_p", "original_empirical_q",
    "orig_transport_class", "mature_transport_class", "n_up_detected",
    "n_down_detected",
]

ORIGINAL_RESULT_FILES = [
    "results/followup-mature-p1m-genes.csv",
    "results/followup-mature-p1m-p0m-genes.csv",
    "results/followup-organ-10-eye-hrca-donors.csv",
    "results/followup-organ-10-eye-hrca-fractions.csv",
    "results/followup-organ-10-eye-hrca-gene-summary.csv",
    "results/followup-organ-10-eye-hrca-top-genes.csv",
    "results/followup-organ-10-eye-hrca-concordance.csv",
    "results/followup-organ-10-eye-hrca-p1-signature.csv",
]

SUITES = [
    {"name": "Sn_Adult_MacularOnly", "assay": "sn", "donors_key": "sn_adult_30",
     "regions": ["macular"], "covs": ("STUDY", "sex_male", "log_lib"),
     "min_donors": 15, "min_df": 5, "studies_ref": "Chen_ancestry"},
    {"name": "Sn_Adult_PeripheralOnly", "assay": "sn", "donors_key": "sn_adult_30",
     "regions": ["peripheral"], "covs": ("STUDY", "sex_male", "log_lib"),
     "min_donors": 15, "min_df": 5, "studies_ref": "Chen_ancestry"},
    {"name": "Sn_Adult_Unenriched", "assay": "sn", "donors_key": "sn_adult_30",
     "regions": ["macular", "peripheral"],
     "covs": ("STUDY", "sex_male", "frac_peripheral", "log_lib"),
     "min_donors": 20, "min_df": 5, "studies_ref": "Chen_ancestry"},
    {"name": "Sn_Adult_ExactAge", "assay": "sn", "donors_key": "sn_adult_exact_30",
     "regions": ["macular", "peripheral"],
     "covs": ("STUDY", "sex_male", "frac_peripheral", "log_lib"),
     "min_donors": 20, "min_df": 5, "studies_ref": "Chen_ancestry"},
    {"name": "Sn_ChenAncestry_Adult", "assay": "sn", "donors_key": "sn_anc_adult_30",
     "regions": ["macular", "peripheral"],
     "covs": ("sex_male", "frac_peripheral", "log_lib"),
     "min_donors": 15, "min_df": 5, "studies_ref": "Chen_ancestry"},
    {"name": "Sn_AllAges_Unenriched", "assay": "sn", "donors_key": "sn_all_30",
     "regions": ["macular", "peripheral"],
     "covs": ("STUDY", "sex_male", "frac_peripheral", "log_lib"),
     "min_donors": 20, "min_df": 5, "studies_ref": "Chen_ancestry"},
    {"name": "Sc15_StudyAdj", "assay": "sc", "donors_key": "sc_all",
     "regions": ["macular", "peripheral", "other"],
     "covs": ("STUDY", "disease_any", "sex_male", "log_lib"),
     "min_donors": 6, "min_df": 3, "studies_ref": "Sanes_GSE148077"},
    {"name": "Sc9_External_StudyAdj", "assay": "sc", "donors_key": "sc_ext",
     "regions": ["macular", "peripheral", "other"],
     "covs": ("STUDY", "disease_any", "sex_male", "log_lib"),
     "min_donors": 6, "min_df": 3, "studies_ref": "Hafler_GSE137537"},
    {"name": "Sc5_External_NoDisease", "assay": "sc", "donors_key": "sc_ext_nodis",
     "regions": ["macular", "peripheral", "other"],
     "covs": ("STUDY", "sex_male", "log_lib"),
     "min_donors": 5, "min_df": 2, "studies_ref": "Hafler_GSE137537"},
]
SUITE_BY_NAME = {s["name"]: s for s in SUITES}


def available_ram_gib() -> float:
    with open("/proc/meminfo", encoding="ascii") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024.0 ** 2)
    raise RuntimeError("Could not read MemAvailable from /proc/meminfo")


def wait_for_memory(stage: str, min_gib: float = MIN_AVAIL_RAM_GIB,
                    timeout_sec: int = RAM_WAIT_TIMEOUT_SEC) -> float:
    deadline = time.monotonic() + timeout_sec
    while True:
        avail = available_ram_gib()
        if avail >= min_gib:
            print(f"[RAM OK] {stage}: {avail:.2f} GiB available (>= {min_gib:.1f} GiB)", flush=True)
            return avail
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise MemoryError(
                f"{stage}: available RAM {avail:.2f} GiB stayed below {min_gib:.1f} GiB for {timeout_sec}s")
        sleep_s = min(30.0, remaining)
        print(f"[RAM WAIT] {stage}: {avail:.2f} GiB available < {min_gib:.1f} GiB; sleeping {sleep_s:.0f}s",
              flush=True)
        time.sleep(sleep_s)


def append_runtime(cohort: str, step: str, start_utc: str, t0: float, status: str, notes: str = "") -> None:
    row = pd.DataFrame([{
        "cohort": cohort,
        "step": step,
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "duration_seconds": round(time.monotonic() - t0, 3),
        "peak_rss_mb": round(peak_rss_mb(), 2),
        "model_used": MODEL_USED,
        "status": status,
        "notes": sanitize_dashes(notes),
    }])
    if RUNTIME_CSV.exists():
        prev = pd.read_csv(RUNTIME_CSV)
        prev = prev[~((prev["cohort"] == cohort) & (prev["step"] == step))]
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


def record_sha256_manifest() -> dict[str, str]:
    manifest: dict[str, str] = {}
    for rel in ORIGINAL_RESULT_FILES:
        p = REPO_ROOT / rel
        if not p.exists():
            raise FileNotFoundError(f"Expected original input missing: {p}")
        manifest[rel] = sha256_file(p)
    TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
    snap_path = TMP_WORK_DIR / "original_sha256_before.json"
    if not snap_path.exists():
        snap_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def verify_sha256_manifest() -> None:
    before_path = TMP_WORK_DIR / "original_sha256_before.json"
    before_map = json.loads(before_path.read_text(encoding="utf-8"))
    after_map: dict[str, str] = {}
    for rel, h_before in before_map.items():
        h_after = sha256_file(REPO_ROOT / rel)
        after_map[rel] = h_after
        if h_after != h_before:
            raise AssertionError(f"Original file modified! {rel}: before={h_before}, after={h_after}")
    (TMP_WORK_DIR / "original_sha256_verified.json").write_text(
        json.dumps({"verified_utc": utc_now(), "all_unchanged": True, "files": after_map}, indent=2) + "\n",
        encoding="utf-8")
    print(f"[SHA256] Verified {len(after_map)} original input files unchanged.", flush=True)


def load_originals() -> dict[str, Any]:
    orig_frac = pd.read_csv(RESULTS_DIR / "followup-organ-10-eye-hrca-fractions.csv")
    orig_gsum = pd.read_csv(RESULTS_DIR / "followup-organ-10-eye-hrca-gene-summary.csv")
    orig_top = pd.read_csv(RESULTS_DIR / "followup-organ-10-eye-hrca-top-genes.csv")
    orig_conc = pd.read_csv(RESULTS_DIR / "followup-organ-10-eye-hrca-concordance.csv")
    orig_p1 = pd.read_csv(RESULTS_DIR / "followup-organ-10-eye-hrca-p1-signature.csv")
    orig_donors = pd.read_csv(RESULTS_DIR / "followup-organ-10-eye-hrca-donors.csv")
    return {
        "frac": {(str(r.analysis), str(r.cell_type)): r for r in orig_frac.itertuples(index=False)},
        "frac_by_suite": {a: sorted(g["cell_type"].astype(str).tolist())
                          for a, g in orig_frac.groupby("analysis")},
        "gsum": {(str(r.analysis), str(r.group)): r for r in orig_gsum.itertuples(index=False)},
        "gsum_by_suite": {a: sorted(g["group"].astype(str).tolist())
                          for a, g in orig_gsum.groupby("analysis")},
        "top": {(str(r.analysis), str(r.group), str(r.gene)): r
                for r in orig_top.itertuples(index=False)},
        "conc": {str(r.analysis): r for r in orig_conc.itertuples(index=False)},
        "p1": {(str(r.analysis), str(r.group)): r for r in orig_p1.itertuples(index=False)},
        "p1_by_suite": {a: sorted(g["group"].astype(str).tolist())
                        for a, g in orig_p1.groupby("analysis")},
        "suite_n": orig_donors.groupby("analysis")["donor_id"].nunique().to_dict(),
    }

def _s(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x)


def run_suite(suite: dict[str, Any], agg: dict[str, Any], programs: dict[str, dict[str, list[str]]],
              orig: dict[str, Any], donor_lists: dict[str, list[str]]) -> dict[str, Any]:
    name = str(suite["name"])
    cohort = "HRCA_" + name
    t0_all = time.monotonic()
    s0_all = utc_now()
    print(f"\n================ SUITE: {name} ================", flush=True)
    wait_for_memory(f"{name} start")

    keys = agg["keys"]
    donors = donor_lists[str(suite["donors_key"])]
    regions = list(suite["regions"])
    k = keys[keys["donor_id"].isin(donors) & keys["regiongroup"].isin(regions)].copy()
    meta = _donor_meta(k, donors, str(suite["studies_ref"]))
    meta["age_years"] = meta["age"]
    assert len(meta) >= 4, f"{name}: only {len(meta)} donors in meta"
    assert float(meta["age_years"].min()) >= MIN_AGE, f"{name}: under-30 donor leaked into meta"
    orig_n_suite = int(orig["suite_n"][name])
    n_mature_suite = int(len(meta))
    n_excluded_suite = orig_n_suite - n_mature_suite
    print(f"[{name}] donors orig={orig_n_suite} mature={n_mature_suite} excluded={n_excluded_suite} "
          f"age {float(meta['age_years'].min()):.0f} to {float(meta['age_years'].max()):.0f}",
          flush=True)

    study_covs = tuple(c for c in meta.columns if c.startswith("study_"))
    cov_list: list[str] = []
    for c in suite["covs"]:
        if c == "STUDY":
            cov_list.extend(study_covs)
        else:
            cov_list.append(c)
    covariates = tuple(cov_list)
    cov_a = tuple(c for c in covariates if c != "log_lib")
    min_donors = int(suite["min_donors"])
    min_df = int(suite["min_df"])

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []

    # ---- Family A: CLR fractions ----
    t0_a = time.monotonic()
    s0_a = utc_now()
    cc = k.pivot_table(index="donor_id", columns="majorclass", values="n_cells", aggfunc="sum", fill_value=0)
    cc = cc.reindex(index=meta.index, columns=[c for c in CLASSES if c in cc.columns], fill_value=0)
    fam_a = run_family_a_clr(
        cc, meta, age_col="age_per_decade", candidate_covariates=cov_a,
        min_donors_present=min(min_donors, len(meta)), min_total_cells=20, min_residual_df=min_df)
    if fam_a.empty:
        q_a_exp = np.array([], dtype=np.float64)
        q_a_sm = np.array([], dtype=np.float64)
        bh_diff_a = 0.0
    else:
        q_a_exp, q_a_sm, bh_diff_a = bh_checked(fam_a["pvalue"])
        fam_a["q_bh_family_a"] = q_a_exp
    append_runtime(cohort, "family_a_fractions", s0_a, t0_a, "ok", f"n_tested={len(fam_a)}")

    for idx_r, r in fam_a.iterrows():
        ct = str(r["cell_type"])
        orow = orig["frac"][(name, ct)]
        orig_sig = bool(float(orow.q_bh_family_a) < ALPHA)
        mat_sig = bool(float(r["q_bh_family_a"]) < ALPHA)
        verd, vnote = verdict_standard(orig_sig, mat_sig, float(orow.effect_clr), float(r["effect_clr"]))
        mde = mde80_from_se(float(r["se"]), int(r["residual_df"]))
        out_rows.append(make_row(
            cohort, ct, "fractions_clr", int(r["n_donors"]), int(orow.n_donors) - int(r["n_donors"]),
            float(orow.pvalue), float(orow.q_bh_family_a),
            float(r["pvalue"]), float(r["q_bh_family_a"]), None, None, verd,
            f"n_present={int(r['n_donors_present'])}" + (f";{vnote}" if vnote else ""),
            arm="primary", residual_df=int(r["residual_df"]),
            original_effect=float(orow.effect_clr), mature_effect=float(r["effect_clr"]),
            mature_se=float(r["se"]), null_sd=np.nan, mde80=mde, mde80_method="t_quantile_se",
            covariates=str(r["covariates"]), omitted_covariates=str(r["omitted_covariates"])))
        bh_rows.append({
            "cohort": cohort, "family": f"{cohort}_fractions", "unit": ct,
            "family_size": int(len(fam_a)), "raw_p": float(r["pvalue"]),
            "q_bh_explicit": float(q_a_exp[idx_r]), "q_bh_statsmodels": float(q_a_sm[idx_r]),
            "abs_diff": float(abs(q_a_exp[idx_r] - q_a_sm[idx_r])), "family_max_abs_diff": bh_diff_a})
    tested_frac_cts = set(fam_a["cell_type"].astype(str).tolist()) if not fam_a.empty else set()
    for ct in orig["frac_by_suite"][name]:
        if ct in tested_frac_cts:
            continue
        orow = orig["frac"][(name, ct)]
        n_present = int((cc[ct] > 0).sum()) if ct in cc.columns else 0
        out_rows.append(make_row(
            cohort, ct, "fractions_clr", n_mature_suite, n_excluded_suite,
            float(orow.pvalue), float(orow.q_bh_family_a), None, None, None, None,
            "UNTESTABLE",
            f"untested_at_30plus;n_present={n_present};min_donors_present={min(min_donors, len(meta))}",
            arm="primary", residual_df=np.nan, original_effect=float(orow.effect_clr),
            mature_effect=np.nan, mature_se=np.nan, null_sd=np.nan, mde80=np.nan,
            mde80_method="", covariates="", omitted_covariates=""))

    # ---- Finalize pseudobulks strictly on mature donors ----
    genes = agg["genes"]
    expr_by: dict[str, pd.DataFrame] = {}
    meta_by: dict[str, pd.DataFrame] = {}
    for ct in cc.columns:
        eligible = [d for d in meta.index if int(cc.loc[d, ct]) >= 10]
        if len(eligible) < min_donors:
            continue
        cols: dict[str, np.ndarray] = {}
        for d in eligible:
            idx = k.index[(k["donor_id"] == d) & (k["majorclass"] == ct)]
            cols[d] = agg["sums"][idx].sum(axis=0)
        raw_df = pd.DataFrame(cols, index=genes)
        e, m, _ = finalize_pseudobulk_group(raw_df, meta.loc[eligible], cc.loc[eligible, ct])
        expr_by[ct], meta_by[ct] = e, m
    print(f"[{name}] pseudobulk groups: {sorted(expr_by.keys())}", flush=True)

    # ---- Family B: genes + concordance ----
    wait_for_memory(f"{name} Family B")
    t0_b = time.monotonic()
    s0_b = utc_now()
    genes_df, sum_df, conc_df = run_family_b_mature(
        expr_by, meta_by, age_col="age_per_decade", candidate_covariates=covariates,
        min_donors=min_donors, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=min_df,
        run_concordance=True, concordance_perms=CONCORDANCE_PERMS, seed=ORIG_NULL_SEED)
    if genes_df.empty:
        q_g_exp = np.array([], dtype=np.float64)
        q_g_sm = np.array([], dtype=np.float64)
        bh_diff_g = 0.0
        genes_idx: dict[tuple[str, str], pd.Series] = {}
        genes_pos: dict[tuple[str, str], int] = {}
    else:
        q_g_exp, q_g_sm, bh_diff_g = bh_checked(genes_df["pvalue"])
        stored_q = genes_df["q_bh_global"].to_numpy(dtype=np.float64)
        assert float(np.nanmax(np.abs(q_g_exp - stored_q))) <= 1e-12, f"{name}: genes q recompute mismatch"
        genes_idx = {(str(g), str(gn)): row for (g, gn), row in
                     zip(zip(genes_df["group"], genes_df["gene"]),
                         genes_df.itertuples(index=False))}
        genes_pos = {(str(g), str(gn)): i for i, (g, gn) in
                     enumerate(zip(genes_df["group"], genes_df["gene"]))}
    n_b_tests = int(len(genes_df))
    n_b_sig = int((q_g_exp < ALPHA).sum()) if n_b_tests else 0
    append_runtime(cohort, "family_b_genes_concordance", s0_b, t0_b, "ok",
                   f"tests={n_b_tests}, sig_global={n_b_sig}")

    sum_by_group: dict[str, Any] = {}
    if not sum_df.empty:
        sum_by_group = {str(r.group): r for r in sum_df.itertuples(index=False)}
    for grp in sorted(set(orig["gsum_by_suite"][name]) | set(sum_by_group.keys())):
        orow = orig["gsum"].get((name, grp))
        srec = sum_by_group.get(grp)
        if srec is None:
            out_rows.append(make_row(
                cohort, grp, "genes_celltype", 0,
                int(orow.n_donors) if orow is not None else orig_n_suite,
                float(orow.min_pvalue) if orow is not None else np.nan,
                float(orow.min_q_bh_global) if orow is not None else np.nan,
                None, None, None, None, "UNTESTABLE", "no_donors_with_min_cells_at_30plus",
                arm="primary", residual_df=np.nan, original_effect=np.nan, mature_effect=np.nan,
                mature_se=np.nan, null_sd=np.nan, mde80=np.nan, mde80_method="",
                n_genes_tested=0, n_sig_q05=0, covariates="", omitted_covariates=""))
            continue
        status = str(srec.status)
        if status != "tested":
            out_rows.append(make_row(
                cohort, grp, "genes_celltype", int(srec.n_donors),
                int(orow.n_donors) - int(srec.n_donors) if orow is not None else n_excluded_suite,
                float(orow.min_pvalue) if orow is not None else np.nan,
                float(orow.min_q_bh_global) if orow is not None else np.nan,
                None, None, None, None, "UNTESTABLE", f"status={status}",
                arm="primary", residual_df=np.nan, original_effect=np.nan, mature_effect=np.nan,
                mature_se=np.nan, null_sd=np.nan, mde80=np.nan, mde80_method="",
                n_genes_tested=0, n_sig_q05=0,
                covariates=_s(getattr(srec, "covariates", "")),
                omitted_covariates=_s(getattr(srec, "omitted_covariates", ""))))
            continue
        orig_nsig = int(orow.n_sig_global_q05) if orow is not None else 0
        mat_nsig = int(srec.n_sig_global_q05)
        orig_minq = float(orow.min_q_bh_global) if orow is not None else np.nan
        verd, _ = verdict_standard(orig_nsig > 0, mat_nsig > 0)
        out_rows.append(make_row(
            cohort, grp, "genes_celltype", int(srec.n_donors),
            int(orow.n_donors) - int(srec.n_donors) if orow is not None else n_excluded_suite,
            float(orow.min_pvalue) if orow is not None else np.nan, orig_minq,
            float(srec.min_pvalue), float(srec.min_q_bh_global), None, None, verd,
            f"n_genes={int(srec.n_genes_tested)};orig_sig_global={orig_nsig};"
            f"mature_sig_global={mat_nsig};orig_sig_within={int(orow.n_sig_within_q05) if orow is not None else 0};"
            f"mature_sig_within={int(srec.n_sig_within_q05)};top={srec.top_gene}",
            arm="primary", residual_df=int(srec.residual_df), original_effect=np.nan,
            mature_effect=float(srec.top_gene_effect), mature_se=float(srec.median_se),
            null_sd=np.nan, mde80=float(srec.median_mde80), mde80_method="median_t_quantile_se",
            n_genes_tested=int(srec.n_genes_tested), n_sig_q05=mat_nsig,
            covariates=str(srec.covariates), omitted_covariates=str(srec.omitted_covariates)))
        g_sub = genes_df[genes_df["group"] == grp].sort_values("pvalue").iloc[0]
        pos = genes_pos[(grp, str(g_sub["gene"]))]
        bh_rows.append({
            "cohort": cohort, "family": f"{cohort}_genes_global", "unit": f"{grp}::{g_sub['gene']}",
            "family_size": n_b_tests, "raw_p": float(g_sub["pvalue"]),
            "q_bh_explicit": float(q_g_exp[pos]), "q_bh_statsmodels": float(q_g_sm[pos]),
            "abs_diff": float(abs(q_g_exp[pos] - q_g_sm[pos])), "family_max_abs_diff": bh_diff_g})

    # ---- Individual gene hit tracking ----
    orig_sig_pairs = {(str(g), str(gn)) for (a, g, gn), r in orig["top"].items()
                      if a == name and float(r.q_bh_global) < ALPHA}
    mat_sig_pairs = {(str(g), str(gn)) for (g, gn) in genes_pos
                     if float(genes_df.iloc[genes_pos[(g, gn)]]["q_bh_global"]) < ALPHA} if n_b_tests else set()
    for grp, gene in sorted(orig_sig_pairs | mat_sig_pairs):
        orig_rec = orig["top"].get((name, grp, gene))
        orig_known = orig_rec is not None
        orig_sig = bool(orig_known and float(orig_rec.q_bh_global) < ALPHA)
        if (grp, gene) in genes_idx:
            mrec = genes_idx[(grp, gene)]
            mat_p = float(mrec.pvalue)
            mat_q = float(mrec.q_bh_global)
            mat_q_within = float(mrec.q_bh_within_group)
            mat_eff = float(mrec.effect)
            mat_se = float(mrec.se)
            mat_df = int(mrec.residual_df)
            mat_n = int(mrec.n_donors)
            mat_mde = float(mrec.mde80)
            testable = True
        else:
            mat_p = np.nan
            mat_q = np.nan
            mat_q_within = np.nan
            mat_eff = np.nan
            mat_se = np.nan
            mat_df = np.nan
            srec_g = sum_by_group.get(grp)
            mat_n = int(srec_g.n_donors) if srec_g is not None else 0
            mat_mde = np.nan
            testable = False
        if orig_known:
            orig_p = float(orig_rec.pvalue)
            orig_q = float(orig_rec.q_bh_global)
            orig_eff = float(orig_rec.effect)
            n_ex = int(orig_rec.n_donors) - mat_n
        else:
            orig_p = np.nan
            orig_q = np.nan
            orig_eff = np.nan
            orow_g = orig["gsum"].get((name, grp))
            n_ex = int(orow_g.n_donors) - mat_n if orow_g is not None else n_excluded_suite
        if not testable:
            verd, vnote = "UNTESTABLE", "gene_failed_mature_50pct_detection_filter"
        elif not orig_known:
            verd, vnote = "", "orig_status_unknown_outside_top50"
        else:
            verd, vnote = verdict_standard(orig_sig, bool(mat_q < ALPHA), orig_eff, mat_eff)
        note_str = (f"orig_sig_global={orig_sig};mature_q_global={mat_q:.6g}" if np.isfinite(mat_q)
                    else f"orig_sig_global={orig_sig};unfiltered_out")
        if np.isfinite(mat_q_within):
            note_str += f";mature_q_within={mat_q_within:.6g}"
        if vnote:
            note_str += f";{vnote}"
        out_rows.append(make_row(
            cohort, f"{grp}::{gene}", "genes_hit", mat_n, n_ex, orig_p, orig_q,
            mat_p, mat_q, None, None, verd, note_str,
            arm="primary", residual_df=mat_df, original_effect=orig_eff, mature_effect=mat_eff,
            mature_se=mat_se, null_sd=np.nan, mde80=mat_mde, mde80_method="t_quantile_se"))
        if testable:
            pos = genes_pos[(grp, gene)]
            bh_rows.append({
                "cohort": cohort, "family": f"{cohort}_genes_global",
                "unit": f"{grp}::{gene}", "family_size": n_b_tests, "raw_p": mat_p,
                "q_bh_explicit": float(q_g_exp[pos]), "q_bh_statsmodels": float(q_g_sm[pos]),
                "abs_diff": float(abs(q_g_exp[pos] - q_g_sm[pos])),
                "family_max_abs_diff": bh_diff_g})

    # ---- Concordance row ----
    orow_c = orig["conc"][name]
    orig_c_p = float(orow_c.perm_p_two_sided)
    if not conc_df.empty and str(conc_df.iloc[0].get("status", "")) == "tested":
        crow = conc_df.iloc[0]
        mat_c_p = float(crow["perm_p_two_sided"])
        verd_c, _ = verdict_standard(orig_c_p < ALPHA, mat_c_p < ALPHA,
                                     float(orow_c.mean_pairwise_pearson_r),
                                     float(crow["mean_pairwise_pearson_r"]))
        out_rows.append(make_row(
            cohort, f"All_{int(crow['n_cell_types'])}_CellTypes", "concordance",
            int(crow["n_donors_union"]), orig_n_suite - int(crow["n_donors_union"]),
            orig_c_p, orig_c_p, mat_c_p, mat_c_p, mat_c_p, mat_c_p, verd_c,
            f"n_pairs={int(crow['n_pairs'])};n_common_genes={int(crow['n_common_genes'])};"
            f"median_r={float(crow['median_pairwise_pearson_r']):.4f};"
            f"sign_agree={float(crow['mean_sign_agreement']):.4f};"
            f"perms={int(crow['n_valid_perms'])};seed={ORIG_NULL_SEED}",
            arm="primary", residual_df=np.nan,
            original_effect=float(orow_c.mean_pairwise_pearson_r),
            mature_effect=float(crow["mean_pairwise_pearson_r"]), mature_se=np.nan,
            null_sd=float(crow["perm_null_sd"]), mde80=float(crow["mde80_null"]),
            mde80_method="2.4865_times_null_sd"))
    else:
        reason = "concordance_untestable"
        if not conc_df.empty:
            reason = str(conc_df.iloc[0].get("status", reason))
        out_rows.append(make_row(
            cohort, "All_CellTypes", "concordance", n_mature_suite, n_excluded_suite,
            orig_c_p, orig_c_p, None, None, None, None, "UNTESTABLE", f"reason={reason}",
            arm="primary", residual_df=np.nan,
            original_effect=float(orow_c.mean_pairwise_pearson_r), mature_effect=np.nan,
            mature_se=np.nan, null_sd=np.nan, mde80=np.nan, mde80_method=""))

    # ---- Family C: transport P1m and P0m ----
    wait_for_memory(f"{name} Family C")
    t0_c = time.monotonic()
    s0_c = utc_now()
    fam_c_p1m, skip_p1m = run_family_c_mature(
        expr_by, meta_by, program=programs["P1m"], candidate_covariates=covariates,
        min_donors=min_donors, min_residual_df=min_df,
        iterations=P1_ITERATIONS, seed=ORIG_NULL_SEED)
    fam_c_p0m, skip_p0m = run_family_c_mature(
        expr_by, meta_by, program=programs["P0m"], candidate_covariates=covariates,
        min_donors=min_donors, min_residual_df=min_df,
        iterations=P1_ITERATIONS, seed=ORIG_NULL_SEED)
    append_runtime(cohort, "family_c_transport", s0_c, t0_c, "ok",
                   f"P1m_tested={len(fam_c_p1m)}, P0m_tested={len(fam_c_p0m)}, "
                   f"seed={ORIG_NULL_SEED}, draws={P1_ITERATIONS}")

    for tier, df_c, skips in [("P1m", fam_c_p1m, skip_p1m), ("P0m", fam_c_p0m, skip_p0m)]:
        fam_key = f"transport_{tier}"
        if df_c.empty:
            q_c_exp = np.array([], dtype=np.float64)
            q_c_sm = np.array([], dtype=np.float64)
            diff_c = 0.0
            q_e_exp = np.array([], dtype=np.float64)
            q_e_sm = np.array([], dtype=np.float64)
            diff_e = 0.0
            df_c_idx: dict[str, Any] = {}
            pos_map: dict[str, int] = {}
        else:
            q_c_exp, q_c_sm, diff_c = bh_checked(df_c["composite_pvalue"])
            q_e_exp, q_e_sm, diff_e = bh_checked(df_c["empirical_null_pvalue"])
            assert float(np.nanmax(np.abs(q_c_exp - df_c["q_bh_composite"].to_numpy(dtype=np.float64)))) <= 1e-12
            assert float(np.nanmax(np.abs(q_e_exp - df_c["q_bh_empirical"].to_numpy(dtype=np.float64)))) <= 1e-12
            df_c_idx = {str(r.group): r for r in df_c.itertuples(index=False)}
            pos_map = {str(g): i for i, g in enumerate(df_c["group"])}
        for grp in sorted(set(orig["p1_by_suite"][name]) | set(df_c_idx.keys())):
            orow_p1 = orig["p1"].get((name, grp))
            orig_n = int(orow_p1.n_donors) if orow_p1 is not None else orig_n_suite
            if tier == "P1m" and orow_p1 is not None:
                op_c = float(orow_p1.composite_pvalue)
                oq_c = float(orow_p1.q_bh_composite)
                op_e = float(orow_p1.empirical_null_pvalue)
                oq_e = float(orow_p1.q_bh_empirical)
                oeff = float(orow_p1.composite_effect)
                orig_cls = transport_class(oq_c, oq_e, True)
                orig_both: bool | None = orig_cls == "BOTH_PASS"
            else:
                op_c = np.nan
                oq_c = np.nan
                op_e = np.nan
                oq_e = np.nan
                oeff = np.nan
                orig_cls = ""
                orig_both = None
            if grp not in df_c_idx:
                reason = skips.get(grp, "ineligible_at_30plus")
                out_rows.append(make_row(
                    cohort, grp, fam_key, 0, orig_n,
                    (op_c if tier == "P1m" else None), (oq_c if tier == "P1m" else None),
                    None, None, None, None, "UNTESTABLE",
                    f"reason={reason}" + (f";orig_class={orig_cls}" if orig_cls else
                                           ";no_original_P0_baseline"),
                    arm="primary", original_effect=(oeff if tier == "P1m" else np.nan),
                    original_empirical_p=(op_e if tier == "P1m" else np.nan),
                    original_empirical_q=(oq_e if tier == "P1m" else np.nan),
                    orig_transport_class=orig_cls, mature_transport_class="UNTESTABLE"))
                continue
            r = df_c_idx[grp]
            pos = pos_map[grp]
            mat_cls = transport_class(float(r.q_bh_composite), float(r.q_bh_empirical), True)
            mat_both = mat_cls == "BOTH_PASS"
            if tier == "P1m":
                verd, vnote = verdict_standard(orig_both, mat_both, oeff, float(r.composite_effect))
            else:
                verd, vnote = "", "no_original_P0_baseline"
            mde_null = mde80_from_null_sd(float(r.null_sd_effect))
            mde_se = mde80_from_se(float(r.composite_se), int(r.residual_df))
            note_parts = [
                f"orig_class={orig_cls}" if orig_cls else "orig_class=none",
                f"mature_class={mat_cls}",
                f"up_eff={float(r.up_effect):+.4f}(p={float(r.up_pvalue):.3g})",
                f"dn_eff={float(r.down_effect):+.4f}(p={float(r.down_pvalue):.3g})",
                f"mde80_se={mde_se:.4f}",
                f"draws={P1_ITERATIONS}",
                f"seed={ORIG_NULL_SEED}",
                f"min_decile_margin={int(r.min_decile_margin)}",
            ]
            if vnote:
                note_parts.append(vnote)
            out_rows.append(make_row(
                cohort, grp, fam_key, int(r.n_donors), orig_n - int(r.n_donors),
                (op_c if tier == "P1m" else None), (oq_c if tier == "P1m" else None),
                float(r.composite_pvalue), float(r.q_bh_composite),
                float(r.empirical_null_pvalue), float(r.q_bh_empirical), verd, ";".join(note_parts),
                arm="primary", residual_df=int(r.residual_df),
                original_effect=(oeff if tier == "P1m" else np.nan),
                mature_effect=float(r.composite_effect), mature_se=float(r.composite_se),
                null_sd=float(r.null_sd_effect), mde80=mde_null,
                mde80_method="2.4865_times_null_sd",
                original_empirical_p=(op_e if tier == "P1m" else np.nan),
                original_empirical_q=(oq_e if tier == "P1m" else np.nan),
                orig_transport_class=orig_cls, mature_transport_class=mat_cls,
                n_up_detected=int(r.n_up_detected), n_down_detected=int(r.n_down_detected),
                covariates=str(r.covariates), omitted_covariates=str(r.omitted_covariates)))
            bh_rows.append({
                "cohort": cohort, "family": f"{cohort}_{fam_key}_composite", "unit": grp,
                "family_size": int(len(df_c)), "raw_p": float(r.composite_pvalue),
                "q_bh_explicit": float(q_c_exp[pos]), "q_bh_statsmodels": float(q_c_sm[pos]),
                "abs_diff": float(abs(q_c_exp[pos] - q_c_sm[pos])), "family_max_abs_diff": diff_c})
            bh_rows.append({
                "cohort": cohort, "family": f"{cohort}_{fam_key}_empirical", "unit": grp,
                "family_size": int(len(df_c)), "raw_p": float(r.empirical_null_pvalue),
                "q_bh_explicit": float(q_e_exp[pos]), "q_bh_statsmodels": float(q_e_sm[pos]),
                "abs_diff": float(abs(q_e_exp[pos] - q_e_sm[pos])), "family_max_abs_diff": diff_e})

    # ---- Diagnostics ----
    def _both_pass(df_c: pd.DataFrame) -> int:
        if df_c.empty:
            return 0
        return int(((df_c["q_bh_composite"] < ALPHA) & (df_c["q_bh_empirical"] < ALPHA)).sum())

    diag = {
        "suite": name,
        "cohort": cohort,
        "assay": str(suite["assay"]),
        "donors_original": orig_n_suite,
        "donors_mature": n_mature_suite,
        "donors_excluded": n_excluded_suite,
        "age_min": float(meta["age_years"].min()),
        "age_max": float(meta["age_years"].max()),
        "n_male": int(meta["sex_male"].sum()),
        "n_female": int((1 - meta["sex_male"]).sum()),
        "study_counts": {str(k): int(v) for k, v in meta["study"].value_counts().items()},
        "n_rows": int(len(out_rows)),
        "fractions_tested": int(len(fam_a)),
        "fractions_sig": int((fam_a["q_bh_family_a"] < ALPHA).sum()) if not fam_a.empty else 0,
        "fractions_min_q": float(fam_a["q_bh_family_a"].min()) if not fam_a.empty else None,
        "genes_tests": n_b_tests,
        "genes_sig_global": n_b_sig,
        "genes_groups_tested": int(len([1 for s in sum_by_group.values()
                                        if str(s.status) == "tested"])),
        "concordance_tested": bool(not conc_df.empty
                                   and str(conc_df.iloc[0].get("status", "")) == "tested"),
        "concordance_mean_r": (float(conc_df.iloc[0]["mean_pairwise_pearson_r"])
                               if not conc_df.empty
                               and str(conc_df.iloc[0].get("status", "")) == "tested" else None),
        "concordance_p_two_sided": (float(conc_df.iloc[0]["perm_p_two_sided"])
                                    if not conc_df.empty
                                    and str(conc_df.iloc[0].get("status", "")) == "tested" else None),
        "p1m_both_pass": _both_pass(fam_c_p1m),
        "p1m_tested": int(len(fam_c_p1m)),
        "p0m_both_pass": _both_pass(fam_c_p0m),
        "p0m_tested": int(len(fam_c_p0m)),
        "orig_p1_both_pass": int(sum(
            1 for (a, _g), r in orig["p1"].items() if a == name
            and float(r.q_bh_composite) < ALPHA and float(r.q_bh_empirical) < ALPHA)),
        "orig_fractions_sig": int(sum(
            1 for (a, _c), r in orig["frac"].items() if a == name
            and float(r.q_bh_family_a) < ALPHA)),
        "orig_genes_sig_groups": int(sum(
            1 for (a, _g), r in orig["gsum"].items() if a == name
            and int(r.n_sig_global_q05) > 0)),
        "orig_genes_sig_total": int(sum(
            int(r.n_sig_global_q05) for (a, _g), r in orig["gsum"].items() if a == name)),
    }
    append_runtime(cohort, "suite_total", s0_all, t0_all, "ok",
                   f"rows={len(out_rows)}, P1m_bothpass={diag['p1m_both_pass']}/{diag['p1m_tested']}, "
                   f"P0m_bothpass={diag['p0m_both_pass']}/{diag['p0m_tested']}")
    print(f"[{name} DONE] rows={len(out_rows)}, elapsed={time.monotonic() - t0_all:.1f}s, "
          f"P1m both-pass={diag['p1m_both_pass']}/{diag['p1m_tested']}, "
          f"P0m both-pass={diag['p0m_both_pass']}/{diag['p0m_tested']}", flush=True)
    payload = {"rows": out_rows, "bh_rows": bh_rows, "diag": diag}
    TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
    tmp_p = TMP_WORK_DIR / f"checkpoint_{name}.pkl.tmp"
    with tmp_p.open("wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    tmp_p.replace(TMP_WORK_DIR / f"checkpoint_{name}.pkl")
    return payload


def build_donor_lists(sn_keys: pd.DataFrame, sc_keys: pd.DataFrame) -> dict[str, list[str]]:
    dages = sn_keys.groupby("donor_id")["age"].first()
    dtop = sn_keys.groupby("donor_id")["donor_age_raw"].first().str.startswith(">89")
    sn_adult = sorted(dages.index[dages >= 18])
    sn_all = sorted(dages.index)
    sn_adult_exact = sorted(dages.index[(dages >= 18) & ~dtop])
    sn_anc_donors = set(sn_keys.loc[sn_keys["study"] == "Chen_ancestry", "donor_id"].unique())
    sn_anc_adult = sorted(d for d in sn_adult if d in sn_anc_donors)
    sc_all = sorted(sc_keys["donor_id"].unique())
    sc_ext = sorted(sc_keys.loc[sc_keys["study"] != "Sanes_GSE148077", "donor_id"].unique())
    sc_ext_nodis = sorted(
        sc_keys.loc[(sc_keys["study"] != "Sanes_GSE148077") & (sc_keys["disease"] == "none"),
                    "donor_id"].unique())
    return {
        "sn_adult_30": [d for d in sn_adult if float(dages[d]) >= MIN_AGE],
        "sn_all_30": [d for d in sn_all if float(dages[d]) >= MIN_AGE],
        "sn_adult_exact_30": [d for d in sn_adult_exact if float(dages[d]) >= MIN_AGE],
        "sn_anc_adult_30": [d for d in sn_anc_adult if float(dages[d]) >= MIN_AGE],
        "sc_all": sc_all,
        "sc_ext": sc_ext,
        "sc_ext_nodis": sc_ext_nodis,
    }


def load_checkpoints(names: list[str]) -> list[dict[str, Any]]:
    payloads = []
    for name in names:
        p = TMP_WORK_DIR / f"checkpoint_{name}.pkl"
        if not p.exists():
            raise FileNotFoundError(f"Checkpoint missing for suite {name}: {p}")
        with p.open("rb") as fh:
            payloads.append(pickle.load(fh))
    return payloads


def assemble_and_write(payloads: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    for payload in payloads:
        rows.extend(payload["rows"])
        bh_rows.extend(payload["bh_rows"])
    df = order_cohort_df(rows)
    df = df.reindex(columns=EXPECTED_COLS)
    assert list(df.columns) == EXPECTED_COLS, f"Column order mismatch: {list(df.columns)}"
    write_clean_csv(df, COHORT_CSV)
    bh_df = pd.DataFrame(bh_rows).drop_duplicates(
        subset=["cohort", "family", "unit"], keep="first").reset_index(drop=True)
    max_diff = float(bh_df["abs_diff"].max()) if not bh_df.empty else 0.0
    max_fam_diff = float(bh_df["family_max_abs_diff"].max()) if not bh_df.empty else 0.0
    assert max_diff <= 1e-12 and max_fam_diff <= 1e-12, (
        f"BH check failed: max_diff={max_diff}, max_fam_diff={max_fam_diff}")
    write_clean_csv(bh_df, BH_CHECK_CSV)
    assert_no_forbidden_dashes(COHORT_CSV.read_text(encoding="utf-8"), label=str(COHORT_CSV))
    assert_no_forbidden_dashes(BH_CHECK_CSV.read_text(encoding="utf-8"), label=str(BH_CHECK_CSV))
    diag = {p["diag"]["suite"]: p["diag"] for p in payloads}
    (TMP_WORK_DIR / "diagnostics.json").write_text(
        json.dumps(diag, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"[WRITE] cohort rows={len(df)}, bh rows={len(bh_df)}, max_bh_diff={max_diff:.3e}", flush=True)
    return df, bh_df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suites", default="all",
                    help="Comma-separated suite names to (re)compute, or 'all'.")
    ap.add_argument("--merge-checkpoints", action="store_true",
                    help="Only assemble outputs from existing checkpoint pickles.")
    args = ap.parse_args()

    t0_all = time.monotonic()
    s0_all = utc_now()
    record_sha256_manifest()

    if args.suites.strip().lower() == "all":
        wanted = [s["name"] for s in SUITES]
    else:
        wanted = [x.strip() for x in args.suites.split(",") if x.strip()]
        unknown = sorted(set(wanted) - set(SUITE_BY_NAME))
        if unknown:
            raise ValueError(f"Unknown suites: {unknown}")

    if not args.merge_checkpoints:
        programs = load_verified_mature_programs()
        wait_for_memory("HRCA load aggregates")
        sn = _aggregate(SN_PATH, "sn")
        sc = _aggregate(SC_PATH, "sc")
        sck = sc["keys"]
        sck = sck[~((sck["study"] == "Sanes_GSE148077") & (sck["regiongroup"] != "macular"))].copy()
        sc = {**sc, "keys": sck}
        donor_lists = build_donor_lists(sn["keys"], sck)
        print(f"[DONORS] sn_adult_30={len(donor_lists['sn_adult_30'])} "
              f"sn_all_30={len(donor_lists['sn_all_30'])} "
              f"sn_adult_exact_30={len(donor_lists['sn_adult_exact_30'])} "
              f"sn_anc_adult_30={len(donor_lists['sn_anc_adult_30'])} "
              f"sc_all={len(donor_lists['sc_all'])} sc_ext={len(donor_lists['sc_ext'])} "
              f"sc_ext_nodis={len(donor_lists['sc_ext_nodis'])}", flush=True)
        orig = load_originals()
        for suite_name in wanted:
            suite = SUITE_BY_NAME[suite_name]
            agg = sn if suite["assay"] == "sn" else sc
            run_suite(suite, agg, programs, orig, donor_lists)
            gc.collect()
        verify_sha256_manifest()

    payloads = load_checkpoints([s["name"] for s in SUITES])
    df, bh_df = assemble_and_write(payloads)
    append_runtime("HRCA_All", "cohort_total", s0_all, t0_all, "ok",
                   f"rows={len(df)}, bh_rows={len(bh_df)}, suites={len(payloads)}")
    print(f"[HRCA DONE] rows={len(df)}, elapsed={time.monotonic() - t0_all:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
