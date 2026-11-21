#!/usr/bin/env python3
"""Phase 2 Audit Item 5 Batch C: mature (30-plus) follow-up organ transports.

Cohorts in order:
  1. TS2 blood NK plus remaining tissues (Blood non-NK, Bladder, Fat, Lymph_Node, Spleen)
  2. Yan retina (organ 07, Fovea major class primary suite)
  3. RPE bulk (organ 08, GSE159435, primary batch and sex adjusted arm)
  4. Ovary postmenopausal (GSE329417, primary depth plex clinical arm)
  5. PBMC (GSE213516 Zhu and GSE214546 Thomson, mature 30 to 39 vs 60 plus)
  6. Lake kidney (scCv3 and snCv3, primary adult primary age sex arm)
  7. Perez muscle (continuous 30-plus whole-nucleus gene model only; transport blocked)

Frozen protocol (same as Batch A): mature signatures P1m/P0m from Item 1
(commit d3dc14c) loaded via load_verified_mature_programs; seed 20261105;
2000-draw expression-decile-matched null as primary; donor is the inferential
unit, never the cell; donors aged at least 30 before pseudobulk aggregation
and detection filtering, except PBMC where the mature design is the binary
30 to 39 vs 60 plus contrast prescribed by the audit, and Lake where the
primary-adult donor set is restricted to age at least 30; every
multiple-testing family is stated explicitly before BH correction, with BH
cross-checked against statsmodels (assert diff <= 1e-12).

Matrix rules follow each cohort's original construction: TS2, Yan and RPE
use log2CPM pseudobulks with the original CPM >= 1 in >= 50 percent donor
filter via run_family_c_mature; ovary, PBMC, Lake and TS2 NK use the
original mean-of-log per-cell matrices with presence-based program
detection and no CPM filter, via run_transport_meanlog in this script.
Null pools are all genes in the tested matrix in every case.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import pickle
import resource
import sys
import time
from typing import Any, Callable

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_mature_organ_transports_a import (  # noqa: E402
    bh_checked,
    fast_random_signature_scores,
    load_verified_mature_programs,
    mde80_from_null_sd,
    mde80_from_se,
    make_row,
    order_cohort_df,
    run_family_c_mature,
    transport_class,
    verdict_standard,
)
from src.followup_muscle import fit_ols, matched_null_effect_pvalue  # noqa: E402
from src.followup_organs import (  # noqa: E402
    assert_no_forbidden_dashes,
    sanitize_dashes,
    select_design_columns,
    write_clean_csv,
)

RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
TMP_WORK_DIR = DATA_ROOT / "tmp_organs" / "hcata_mature_batch_c_tmp"

MIN_AGE = 30.0
ORIG_NULL_SEED = 20261105
P1_ITERATIONS = 2000
ALPHA = 0.05
MODEL_USED = "Muse Spark (direct, no CLI)"
MIN_AVAIL_RAM_GIB = 6.0
RAM_WAIT_TIMEOUT_SEC = 900

PREFIX = "followup-mature-organ-transports-c"
RUNTIME_CSV = RESULTS_DIR / f"{PREFIX}-runtime.csv"
BH_CHECK_CSV = RESULTS_DIR / f"{PREFIX}-bh-check.csv"

REQUIRED_COHORT_COLS = [
    "cohort",
    "cell_type",
    "endpoint_family",
    "n_donors",
    "n_donors_30plus_excluded",
    "original_p",
    "original_q",
    "mature_p",
    "mature_q",
    "mature_empirical_p",
    "mature_empirical_q",
    "verdict",
    "notes",
]

# Canonical extra-column order, identical to the Batch A Gao CSV header.
EXTRA_COLS = [
    "arm",
    "residual_df",
    "original_effect",
    "mature_effect",
    "mature_se",
    "null_sd",
    "mde80",
    "mde80_method",
    "covariates",
    "omitted_covariates",
    "n_genes_tested",
    "n_sig_q05",
    "original_empirical_p",
    "original_empirical_q",
    "orig_transport_class",
    "mature_transport_class",
    "n_up_detected",
    "n_down_detected",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def peak_rss_mb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0


def available_ram_gib() -> float:
    with open("/proc/meminfo", encoding="ascii") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024.0 ** 2)
    raise RuntimeError("Could not read MemAvailable from /proc/meminfo")


def wait_for_memory(stage: str, min_gib: float = MIN_AVAIL_RAM_GIB, timeout_sec: int = RAM_WAIT_TIMEOUT_SEC) -> float:
    deadline = time.monotonic() + timeout_sec
    while True:
        avail = available_ram_gib()
        if avail >= min_gib:
            print(f"[RAM OK] {stage}: {avail:.2f} GiB available (>= {min_gib:.1f} GiB)", flush=True)
            return avail
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise MemoryError(f"{stage}: available RAM {avail:.2f} GiB stayed below {min_gib:.1f} GiB for {timeout_sec}s")
        sleep_s = min(30.0, remaining)
        print(f"[RAM WAIT] {stage}: {avail:.2f} GiB available < {min_gib:.1f} GiB; sleeping {sleep_s:.0f}s", flush=True)
        time.sleep(sleep_s)


def append_runtime(cohort: str, step: str, start_utc: str, t0: float, status: str, notes: str = "") -> None:
    row = pd.DataFrame(
        [
            {
                "cohort": cohort,
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
        prev = prev[~((prev["cohort"] == cohort) & (prev["step"] == step))]
        out = pd.concat([prev, row], ignore_index=True)
    else:
        out = row
    write_clean_csv(out, RUNTIME_CSV)


def append_bh_rows(new_rows: list[dict[str, Any]]) -> None:
    if not new_rows:
        return
    new_df = pd.DataFrame(new_rows)
    if BH_CHECK_CSV.exists():
        prev = pd.read_csv(BH_CHECK_CSV)
        out = pd.concat([prev, new_df], ignore_index=True)
        out = out.drop_duplicates(subset=["cohort", "family", "unit"], keep="last").reset_index(drop=True)
    else:
        out = new_df
    write_clean_csv(out, BH_CHECK_CSV)


def save_cohort_payload(cohort_key: str, payload: dict[str, Any]) -> None:
    TMP_WORK_DIR.mkdir(parents=True, exist_ok=True)
    with (TMP_WORK_DIR / f"payload_{cohort_key}.pkl").open("wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_cohort_payload(cohort_key: str) -> dict[str, Any]:
    with (TMP_WORK_DIR / f"payload_{cohort_key}.pkl").open("rb") as fh:
        return pickle.load(fh)


DesignBuilder = Callable[[pd.DataFrame], tuple[np.ndarray, list[str], list[str]]]


def run_transport_meanlog(
    groups: dict[str, dict[str, Any]],
    program: dict[str, list[str]],
    design_builder: DesignBuilder,
    min_donors: int = 4,
    min_residual_df: int = 2,
    iterations: int = P1_ITERATIONS,
    seed: int = ORIG_NULL_SEED,
    verify_first_parity: bool = True,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Transport runner for original mean-of-log matrices (ovary, PBMC, Lake, TS2 NK).

    groups maps a group name to a dict with:
      matrix: float ndarray, genes x donors, log scale (mean of per-cell log values)
      genes: ndarray of str gene symbols, one per matrix row
      meta: DataFrame indexed by donor id with age_years and design columns
      raw_counts: optional int ndarray, genes x donors, for the presence gate
    Program genes are detected by presence in the gene list (plus total raw
    count above zero when raw_counts is given, the ovary original rule).
    The null pool is every gene in the matrix. No CPM expression filter is
    applied, matching the original mean-log analyses. BH is NOT applied here;
    the caller applies bh_checked per declared family.
    """
    up_set = set(program["up"])
    down_set = set(program["down"])
    rows: list[dict[str, Any]] = []
    skipped: dict[str, str] = {}
    verified_once = False

    for group, spec in groups.items():
        matrix = np.asarray(spec["matrix"], dtype=np.float64)
        genes = np.asarray(spec["genes"]).astype(str)
        meta = spec["meta"]
        common_donors = [d for d in meta.index if d in spec.get("donors", list(meta.index))]
        if len(common_donors) < min_donors:
            skipped[group] = f"n_donors={len(common_donors)} < min_donors={min_donors}"
            continue
        m_meta = meta.loc[common_donors]
        assert float(m_meta["age_years"].min()) >= MIN_AGE, f"Under-30 donor leaked into {group}"
        donor_pos = {d: i for i, d in enumerate(spec["donors"])}
        cols = [donor_pos[d] for d in common_donors]
        mat = matrix[:, cols]
        raw = spec.get("raw_counts")
        if raw is not None:
            raw_sub = np.asarray(raw)[:, cols]
            detected = raw_sub.sum(axis=1) > 0
        else:
            detected = np.ones(mat.shape[0], dtype=bool)
        detectable = np.asarray(spec.get("detectable", np.ones(mat.shape[0], dtype=bool)), dtype=bool)
        up_idx = np.array([i for i, g in enumerate(genes) if g in up_set and detected[i] and detectable[i]], dtype=int)
        down_idx = np.array([i for i, g in enumerate(genes) if g in down_set and detected[i] and detectable[i]], dtype=int)
        if len(up_idx) == 0 or len(down_idx) == 0:
            skipped[group] = f"detected up={len(up_idx)}, down={len(down_idx)}"
            continue
        try:
            design, names, omitted = design_builder(m_meta)
        except ValueError as exc:
            skipped[group] = f"design_error: {exc}"
            continue
        if design.shape[0] - design.shape[1] < min_residual_df:
            skipped[group] = f"residual_df={design.shape[0] - design.shape[1]} < {min_residual_df}"
            continue

        up_score = mat[up_idx, :].mean(axis=0)
        down_score = mat[down_idx, :].mean(axis=0)
        comp_score = up_score - down_score
        ols = fit_ols(np.column_stack([up_score, down_score, comp_score]), design, effect_column=1)
        try:
            null_scores, min_decile_margin = fast_random_signature_scores(
                mat,
                up_idx,
                down_idx,
                iterations=iterations,
                seed=seed,
                verify_parity=(verify_first_parity and not verified_once),
            )
            verified_once = True
        except ValueError as exc:
            skipped[group] = f"matched_null_pool_error: {exc}"
            continue
        null_ols = fit_ols(null_scores.T, design, effect_column=1)
        null_effects = null_ols["effect"]
        p_emp = matched_null_effect_pvalue(float(ols["effect"][2]), null_effects)
        rows.append(
            {
                "group": str(group),
                "n_donors": int(len(common_donors)),
                "age_min": float(m_meta["age_years"].min()),
                "age_max": float(m_meta["age_years"].max()),
                "n_expressed_pool": int(mat.shape[0]),
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
            }
        )

    if not rows:
        return pd.DataFrame(), skipped
    return pd.DataFrame(rows).sort_values("empirical_null_pvalue", ascending=True).reset_index(drop=True), skipped


def apply_bh_families(
    df: pd.DataFrame,
    family_keys: dict[str, str],
    cohort_label: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Apply bh_checked per declared family and build BH-check rows.

    family_keys maps a group name to its family key (for example a cohort or
    assay name). Composite and empirical p values get separate BH families
    with the same key plus a _composite or _empirical suffix.
    """
    out = df.copy()
    out["q_bh_composite"] = np.nan
    out["q_bh_empirical"] = np.nan
    bh_rows: list[dict[str, Any]] = []
    fams = sorted(set(family_keys.values()))
    for fam in fams:
        members = [g for g in df["group"] if family_keys.get(g) == fam]
        sub = df[df["group"].isin(members)]
        for p_col, q_col, suffix in [
            ("composite_pvalue", "q_bh_composite", "composite"),
            ("empirical_null_pvalue", "q_bh_empirical", "empirical"),
        ]:
            q_exp, q_sm, diff = bh_checked(sub[p_col])
            pos = {g: i for i, g in enumerate(sub["group"])}
            for _, r in sub.iterrows():
                i = pos[r["group"]]
                out.loc[out["group"] == r["group"], q_col] = q_exp[i]
                bh_rows.append(
                    {
                        "cohort": cohort_label,
                        "family": f"{fam}_{suffix}_{len(sub)}",
                        "unit": r["group"],
                        "family_size": int(len(sub)),
                        "raw_p": float(r[p_col]),
                        "q_bh_explicit": float(q_exp[i]),
                        "q_bh_statsmodels": float(q_sm[i]),
                        "abs_diff": float(abs(q_exp[i] - q_sm[i])),
                        "family_max_abs_diff": float(diff),
                    }
                )
    return out, bh_rows


def build_transport_rows(
    cohort_label: str,
    fam_key: str,
    arm: str,
    mature_df: pd.DataFrame,
    skips: dict[str, str],
    orig_ref: pd.DataFrame,
    effect_scale_note: str = "per_decade",
) -> list[dict[str, Any]]:
    """Assemble Batch A schema rows for one transport tier from mature results.

    orig_ref is indexed by group with columns orig_p, orig_q, orig_emp_p,
    orig_emp_q, orig_effect, orig_n, orig_excluded (donors excluded count).
    A group absent from mature_df becomes an UNTESTABLE row.
    """
    rows: list[dict[str, Any]] = []
    df_idx = mature_df.set_index("group") if not mature_df.empty else pd.DataFrame().set_index(pd.Index([]))
    for grp in sorted(orig_ref.index.tolist()):
        o = orig_ref.loc[grp]
        orig_n = int(o["orig_n"])
        orig_cls = transport_class(float(o["orig_q"]), float(o["orig_emp_q"]), True)
        orig_both = orig_cls == "BOTH_PASS"
        oeff = float(o["orig_effect"])
        if grp not in df_idx.index:
            reason = skips.get(grp, "ineligible_at_30plus")
            rows.append(
                make_row(
                    cohort_label,
                    grp,
                    fam_key,
                    0,
                    orig_n,
                    float(o["orig_p"]),
                    float(o["orig_q"]),
                    None,
                    None,
                    None,
                    None,
                    "UNTESTABLE",
                    f"reason={reason};orig_class={orig_cls};{effect_scale_note}",
                    arm=arm,
                    residual_df=None,
                    original_effect=oeff,
                    mature_effect=None,
                    mature_se=None,
                    null_sd=None,
                    mde80=None,
                    mde80_method=None,
                    covariates=None,
                    omitted_covariates=None,
                    n_genes_tested=None,
                    n_sig_q05=None,
                    original_empirical_p=float(o["orig_emp_p"]),
                    original_empirical_q=float(o["orig_emp_q"]),
                    orig_transport_class=orig_cls,
                    mature_transport_class="UNTESTABLE",
                    n_up_detected=None,
                    n_down_detected=None,
                )
            )
            continue
        r = df_idx.loc[grp]
        mat_cls = transport_class(float(r["q_bh_composite"]), float(r["q_bh_empirical"]), True)
        mat_both = mat_cls == "BOTH_PASS"
        verd, vnote = verdict_standard(orig_both, mat_both, oeff, float(r["composite_effect"]))
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
            effect_scale_note,
        ]
        if vnote:
            note_parts.append(vnote)
        rows.append(
            make_row(
                cohort_label,
                grp,
                fam_key,
                int(r["n_donors"]),
                int(o["orig_excluded"]),
                float(o["orig_p"]),
                float(o["orig_q"]),
                float(r["composite_pvalue"]),
                float(r["q_bh_composite"]),
                float(r["empirical_null_pvalue"]),
                float(r["q_bh_empirical"]),
                verd,
                ";".join(note_parts),
                arm=arm,
                residual_df=int(r["residual_df"]),
                original_effect=oeff,
                mature_effect=float(r["composite_effect"]),
                mature_se=float(r["composite_se"]),
                null_sd=float(r["null_sd_effect"]),
                mde80=mde_null,
                mde80_method="2.4865_times_null_sd",
                covariates=str(r["covariates"]),
                omitted_covariates=str(r["omitted_covariates"]),
                n_genes_tested=None,
                n_sig_q05=None,
                original_empirical_p=float(o["orig_emp_p"]),
                original_empirical_q=float(o["orig_emp_q"]),
                orig_transport_class=orig_cls,
                mature_transport_class=mat_cls,
                n_up_detected=int(r["n_up_detected"]),
                n_down_detected=int(r["n_down_detected"]),
            )
        )
    return rows


def write_cohort_csv(cohort_key: str, rows: list[dict[str, Any]]) -> Path:
    df = order_cohort_df(rows)
    # Enforce the exact Batch A extra-column order; order_cohort_df already
    # puts REQUIRED first, then extras in first-appearance order, which is
    # EXTRA_COLS when every make_row call passes extras in that order.
    assert list(df.columns) == REQUIRED_COHORT_COLS + EXTRA_COLS, f"Schema mismatch: {list(df.columns)}"
    path = RESULTS_DIR / f"{PREFIX}-{cohort_key}.csv"
    write_clean_csv(df, path)
    return path


# ---------------------------------------------------------------------------
# Shared helper for run_family_c_mature outputs (TS2, Yan, RPE)
# ---------------------------------------------------------------------------

def bh_rows_from_family_c(
    df_c: pd.DataFrame,
    cohort_label: str,
    family_prefix: str,
) -> list[dict[str, Any]]:
    """Recompute BH with bh_checked and emit per-unit BH-check rows (Batch A pattern)."""
    bh_rows: list[dict[str, Any]] = []
    if df_c.empty:
        return bh_rows
    q_c_exp, q_c_sm, diff_c = bh_checked(df_c["composite_pvalue"])
    q_e_exp, q_e_sm, diff_e = bh_checked(df_c["empirical_null_pvalue"])
    pos_map = {str(g): i for i, g in enumerate(df_c["group"])}
    for grp in df_c["group"]:
        pos = pos_map[str(grp)]
        r = df_c.iloc[pos]
        bh_rows.append(
            {
                "cohort": cohort_label,
                "family": f"{family_prefix}_composite_{len(df_c)}",
                "unit": str(grp),
                "family_size": int(len(df_c)),
                "raw_p": float(r["composite_pvalue"]),
                "q_bh_explicit": float(q_c_exp[pos]),
                "q_bh_statsmodels": float(q_c_sm[pos]),
                "abs_diff": float(abs(q_c_exp[pos] - q_c_sm[pos])),
                "family_max_abs_diff": float(diff_c),
            }
        )
        bh_rows.append(
            {
                "cohort": cohort_label,
                "family": f"{family_prefix}_empirical_{len(df_c)}",
                "unit": str(grp),
                "family_size": int(len(df_c)),
                "raw_p": float(r["empirical_null_pvalue"]),
                "q_bh_explicit": float(q_e_exp[pos]),
                "q_bh_statsmodels": float(q_e_sm[pos]),
                "abs_diff": float(abs(q_e_exp[pos] - q_e_sm[pos])),
                "family_max_abs_diff": float(diff_e),
            }
        )
    return bh_rows


def count_both_pass(df: pd.DataFrame) -> dict[str, int]:
    if df.empty:
        return {"both": 0, "emp_only": 0, "comp_only": 0, "n": 0}
    both = int(((df["q_bh_composite"] < ALPHA) & (df["q_bh_empirical"] < ALPHA)).sum())
    emp_only = int(((df["q_bh_composite"] >= ALPHA) & (df["q_bh_empirical"] < ALPHA)).sum())
    comp_only = int(((df["q_bh_composite"] < ALPHA) & (df["q_bh_empirical"] >= ALPHA)).sum())
    return {"both": both, "emp_only": emp_only, "comp_only": comp_only, "n": int(len(df))}


# ---------------------------------------------------------------------------
# Cohort 1: TS2 (blood NK + remaining tissues)
# ---------------------------------------------------------------------------

def run_cohort_ts2(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    from scripts.run_followup_ts2_remaining import (
        CACHE_ROOT as TS2_CACHE_ROOT,
        MIN_CELLS_PER_DONOR as TS2_MIN_CELLS,
        MIN_DONORS as TS2_MIN_DONORS,
        load_curated_donor_table,
    )
    from src.followup_organs import finalize_pseudobulk_group

    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("ts2")
    tissues = ["Blood", "Bladder", "Fat", "Lymph_Node", "Spleen"]
    cache_paths = {}
    for tissue in tissues:
        t_dir = TS2_CACHE_ROOT / f"{tissue}_cache"
        for name in ("obs_audit.json", "cell_counts_10x3v3.csv", "raw_pb_by_group.pkl"):
            p = t_dir / name
            if not p.exists():
                raise FileNotFoundError(f"TS2 {tissue} input missing: {p}")
        cache_paths[tissue] = t_dir
    nk_cache = DATA_ROOT / "ts2" / "blood_nk_pseudobulk.npz"
    if not nk_cache.exists():
        raise FileNotFoundError(f"TS2 NK input missing: {nk_cache}")

    curated_donors = load_curated_donor_table()
    orig_rem = pd.read_csv(RESULTS_DIR / "followup-ts2-remaining-p1-signature.csv")
    orig_nk = pd.read_csv(RESULTS_DIR / "followup-ts2-blood-nk-p1-model.csv")
    assert len(orig_rem) == 41, f"TS2 remaining original rows {len(orig_rem)} != 41"
    assert len(orig_nk) == 1, f"TS2 NK original rows {len(orig_nk)} != 1"

    # Build log2CPM expression + meta per tissue, exactly as the original
    # run_single_tissue, then pool all 41 groups into one Family C call.
    expr_all: dict[str, pd.DataFrame] = {}
    meta_all: dict[str, pd.DataFrame] = {}
    group_to_tissue: dict[str, str] = {}
    n_donors_orig: dict[str, int] = {}
    for tissue in tissues:
        t_dir = cache_paths[tissue]
        audit = json.loads((t_dir / "obs_audit.json").read_text(encoding="utf-8"))
        cell_counts = pd.read_csv(t_dir / "cell_counts_10x3v3.csv", index_col=0)
        with open(t_dir / "raw_pb_by_group.pkl", "rb") as pf:
            raw_pb_by_group = pickle.load(pf)
        donor_rows = []
        for rec in audit["donors_10x3v3"]:
            did = rec["donor_id"]
            c_row = curated_donors.loc[did]
            age_h5ad = int(str(rec["development_stage"]).split("-")[0])
            age_curated = int(c_row["age_years"])
            assert age_h5ad == age_curated, f"{tissue} {did}: H5AD age {age_h5ad} != curated {age_curated}"
            sex_h5ad = str(rec["sex"])
            assert ("M" if sex_h5ad == "male" else "F") == str(c_row["sex"])
            donor_rows.append(
                {
                    "tissue": tissue,
                    "donor_id": did,
                    "age_years": float(age_curated),
                    "age_per_decade": float(age_curated) / 10.0,
                    "sex": sex_h5ad,
                    "sex_male": int(sex_h5ad == "male"),
                    "n_cells_10x3v3": int(rec["n_cells_10x3v3"]),
                }
            )
        donor_meta = pd.DataFrame(donor_rows).set_index("donor_id")
        assert float(donor_meta["age_years"].min()) >= MIN_AGE, f"TS2 {tissue} has under-30 donors"
        all_grp_name = "All_NonNK_Cells" if tissue == "Blood" else "All_Cells"
        if tissue == "Blood" and "natural killer cell" in cell_counts.columns:
            non_nk_counts = cell_counts.drop(columns=["natural killer cell"]).sum(axis=1)
        else:
            non_nk_counts = cell_counts.sum(axis=1)
        eligible_groups: list[str] = []
        for ct in sorted(cell_counts.columns):
            if tissue == "Blood" and ct == "natural killer cell":
                continue
            if int((cell_counts[ct] >= TS2_MIN_CELLS).sum()) >= TS2_MIN_DONORS:
                eligible_groups.append(ct)
        if int((non_nk_counts >= TS2_MIN_CELLS).sum()) >= TS2_MIN_DONORS:
            eligible_groups.append(all_grp_name)
        for grp in eligible_groups:
            c_ser = non_nk_counts if grp == all_grp_name else cell_counts[grp]
            elig_donors = [d for d in donor_meta.index if int(c_ser.loc[d]) >= TS2_MIN_CELLS]
            raw_sub = raw_pb_by_group[grp][elig_donors]
            log2cpm, grp_meta, _stats = finalize_pseudobulk_group(
                raw_sub, donor_meta.loc[elig_donors], c_ser.loc[elig_donors]
            )
            key = f"{tissue}::{grp}"
            expr_all[key] = log2cpm
            meta_all[key] = grp_meta
            group_to_tissue[key] = tissue
            n_donors_orig[key] = len(elig_donors)
    assert len(expr_all) == 41, f"TS2 remaining groups {len(expr_all)} != 41"

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    tier_results: dict[str, Any] = {}
    for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
        t0 = time.monotonic()
        s0 = utc_now()
        fam_c, skips = run_family_c_mature(
            expr_all,
            meta_all,
            programs[prog_key],
            age_col="age_per_decade",
            candidate_covariates=("sex_male", "log_lib"),
            min_donors=TS2_MIN_DONORS,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=2,
            iterations=P1_ITERATIONS,
            seed=ORIG_NULL_SEED,
        )
        append_runtime("TS2", f"family_c_{tier}_remaining_41", s0, t0, "ok",
                       f"tested={len(fam_c)}, skipped={len(skips)}")
        tier_results[tier] = {"df": fam_c, "skips": skips}
        bh_rows.extend(bh_rows_from_family_c(fam_c, "TS2", f"ts2_{fam_key}"))
        # Original reference for the remaining 41 groups, keyed by Tissue::group.
        orig_ref_rows = []
        for _, o in orig_rem.iterrows():
            key = f"{o['tissue']}::{o['group']}"
            assert int(o["n_donors"]) == n_donors_orig[key], f"TS2 {key}: original n {o['n_donors']} != rebuilt {n_donors_orig[key]}"
            orig_ref_rows.append(
                {
                    "group": key,
                    "orig_p": float(o["composite_pvalue"]),
                    "orig_q": float(o["q_bh_ts2_all41_composite"]),
                    "orig_emp_p": float(o["empirical_null_pvalue"]),
                    "orig_emp_q": float(o["q_bh_ts2_all41_empirical"]),
                    "orig_effect": float(o["composite_effect"]),
                    "orig_n": int(o["n_donors"]),
                    "orig_excluded": 0,
                }
            )
        orig_ref = pd.DataFrame(orig_ref_rows).set_index("group")
        # Build rows per tissue with tissue-specific cohort labels and plain
        # group names as cell_type; the 41-group BH family spans tissues.
        for tissue in tissues:
            tissue_keys = [k for k in orig_ref.index if group_to_tissue.get(k) == tissue]
            sub_ref = orig_ref.loc[tissue_keys].copy()
            sub_ref.index = [k.split("::", 1)[1] for k in tissue_keys]
            sub_df = fam_c[fam_c["group"].isin(tissue_keys)].copy() if not fam_c.empty else fam_c
            if not sub_df.empty:
                sub_df = sub_df.copy()
                sub_df["group"] = sub_df["group"].str.split("::", n=1).str[1]
            tissue_skips = {k.split("::", 1)[1]: v for k, v in skips.items() if group_to_tissue.get(k) == tissue}
            out_rows.extend(
                build_transport_rows(
                    f"TS2_{tissue}", fam_key, "primary_age_sex_loglib",
                    sub_df, tissue_skips, sub_ref,
                )
            )

    # NK single group via the mean-log runner (matrix is mean of log1p CPM).
    nk = np.load(nk_cache, allow_pickle=False)
    nk_donors = pd.read_csv(RESULTS_DIR / "followup-ts2-blood-nk-donor-summary.csv")
    nk_donors = nk_donors.set_index("donor_id")
    nk_matrix = np.asarray(nk["pseudobulk"], dtype=np.float64).T  # genes x donors
    nk_genes = np.asarray(nk["gene_names"]).astype(str)
    nk_order = [str(d) for d in nk["donors"]]
    nk_meta = pd.DataFrame(
        {
            "age_years": nk_donors.loc[nk_order, "age_years"].to_numpy(float),
            "age_per_decade": nk_donors.loc[nk_order, "age_years"].to_numpy(float) / 10.0,
            "male": nk_donors.loc[nk_order, "male"].to_numpy(float),
            "assay_5prime_present": nk_donors.loc[nk_order, "assay_5prime_present"].to_numpy(float),
            "mean_log_umi": nk_donors.loc[nk_order, "mean_log_umi"].to_numpy(float),
        },
        index=nk_order,
    )
    assert float(nk_meta["age_years"].min()) >= MIN_AGE

    def nk_design(m: pd.DataFrame) -> tuple[np.ndarray, list[str], list[str]]:
        design = np.column_stack([
            np.ones(len(m)),
            m["age_per_decade"].to_numpy(float),
            m["male"].to_numpy(float),
            m["assay_5prime_present"].to_numpy(float),
            m["mean_log_umi"].to_numpy(float),
        ])
        assert np.linalg.matrix_rank(design) == design.shape[1]
        return design, ["intercept", "age_per_decade", "male", "assay_5prime_present", "mean_log_umi"], []

    nk_orig = orig_nk.iloc[0]
    nk_ref = pd.DataFrame(
        [
            {
                "group": "natural killer cell",
                "orig_p": float(nk_orig["pvalue"]),
                "orig_q": float(nk_orig["pvalue"]),
                "orig_emp_p": float(nk_orig["empirical_pvalue"]),
                "orig_emp_q": float(nk_orig["q_bh_single_endpoint"]),
                "orig_effect": float(nk_orig["effect_age_10y"]),
                "orig_n": int(nk_orig["n_donors"]),
                "orig_excluded": 0,
            }
        ]
    ).set_index("group")
    for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
        t0 = time.monotonic()
        s0 = utc_now()
        nk_df, nk_skips = run_transport_meanlog(
            {
                "natural killer cell": {
                    "matrix": nk_matrix,
                    "genes": nk_genes,
                    "meta": nk_meta,
                    "donors": nk_order,
                    "raw_counts": None,
                }
            },
            programs[prog_key],
            nk_design,
            min_donors=4,
            min_residual_df=2,
        )
        append_runtime("TS2", f"family_c_{tier}_nk", s0, t0, "ok",
                       f"tested={len(nk_df)}, skipped={len(nk_skips)}")
        if not nk_df.empty:
            nk_df, nk_bh = apply_bh_families(nk_df, {"natural killer cell": f"ts2_nk_{tier}"}, "TS2")
            bh_rows.extend(nk_bh)
        tier_results[f"{tier}_nk"] = {"df": nk_df, "skips": nk_skips}
        out_rows.extend(
            build_transport_rows(
                "TS2_Blood_NK", fam_key, "age_sex_assay5_meanlogumi",
                nk_df, nk_skips, nk_ref,
            )
        )

    path = write_cohort_csv("ts2", out_rows)
    append_bh_rows(bh_rows)
    append_runtime("TS2", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}")
    summary = {
        "cohort": "TS2",
        "rows": len(out_rows),
        "remaining_P1m": count_both_pass(tier_results["P1m"]["df"]),
        "remaining_P0m": count_both_pass(tier_results["P0m"]["df"]),
        "nk_P1m": count_both_pass(tier_results["P1m_nk"]["df"]),
        "nk_P0m": count_both_pass(tier_results["P0m_nk"]["df"]),
        "remaining_skips_P1m": tier_results["P1m"]["skips"],
        "remaining_skips_P0m": tier_results["P0m"]["skips"],
        "nk_skips_P1m": tier_results["P1m_nk"]["skips"],
        "nk_skips_P0m": tier_results["P0m_nk"]["skips"],
    }
    save_cohort_payload("ts2", {"summary": summary, "bh_rows": bh_rows})
    return summary


# ---------------------------------------------------------------------------
# Cohort 2: Yan retina
# ---------------------------------------------------------------------------

def run_cohort_yan(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    from scripts.run_followup_organ_07_eye_yan import (
        _build_donor_meta,
        _load_sanes_from_hrca,
    )
    from src.followup_organs import aggregate_pseudobulk_sparse

    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("yan")
    hrca_path = DATA_ROOT / "tmp_organs" / "shared_hrca_yan" / "hrca.h5ad"
    if not hrca_path.exists():
        raise FileNotFoundError(f"Yan input missing: {hrca_path}")
    orig = pd.read_csv(RESULTS_DIR / "followup-organ-07-eye-yan-gse148077-p1-signature.csv")
    prim = orig[orig["analysis"] == "Fovea_MajorClass_Primary"].copy()
    assert len(prim) == 7, f"Yan primary groups {len(prim)} != 7"

    counts_mat, genes, obs = _load_sanes_from_hrca(hrca_path)
    donor_labels = obs["donor_id"].to_numpy(dtype=object)
    major_labels = obs["majorclass_str"].to_numpy(dtype=object)
    fovea_mask = obs["is_unbiased_fovea"].to_numpy()
    meta_fovea = _build_donor_meta(obs, fovea_mask)
    meta_fovea["age_years"] = meta_fovea["age"].astype(float)
    assert float(meta_fovea["age_years"].min()) >= MIN_AGE, "Yan has under-30 donors"
    _cell_counts, expr_by_group, meta_by_group, _stats = aggregate_pseudobulk_sparse(
        counts_mat,
        genes,
        donor_labels=donor_labels,
        group_labels=major_labels,
        donor_meta=meta_fovea,
        cell_mask=fovea_mask,
        min_cells_per_donor_group=5,
        min_donors_per_group=4,
        allow_fractional_expected=False,
    )
    del counts_mat
    gc.collect()
    assert set(expr_by_group.keys()) == set(prim["group"].tolist()), (
        f"Yan groups {sorted(expr_by_group.keys())} != original {sorted(prim['group'].tolist())}"
    )

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    tier_results: dict[str, Any] = {}
    orig_ref = pd.DataFrame(
        [
            {
                "group": r["group"],
                "orig_p": float(r["composite_pvalue"]),
                "orig_q": float(r["q_bh_composite"]),
                "orig_emp_p": float(r["empirical_null_pvalue"]),
                "orig_emp_q": float(r["q_bh_empirical"]),
                "orig_effect": float(r["composite_effect"]),
                "orig_n": int(r["n_donors"]),
                "orig_excluded": 0,
            }
            for _, r in prim.iterrows()
        ]
    ).set_index("group")
    for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
        t0 = time.monotonic()
        s0 = utc_now()
        fam_c, skips = run_family_c_mature(
            expr_by_group,
            meta_by_group,
            programs[prog_key],
            age_col="age_per_decade",
            candidate_covariates=("log_lib",),
            min_donors=4,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=2,
            iterations=P1_ITERATIONS,
            seed=ORIG_NULL_SEED,
        )
        append_runtime("Yan_Retina", f"family_c_{tier}", s0, t0, "ok",
                       f"tested={len(fam_c)}, skipped={len(skips)}")
        tier_results[tier] = {"df": fam_c, "skips": skips}
        bh_rows.extend(bh_rows_from_family_c(fam_c, "Yan_Retina", f"yan_{fam_key}"))
        out_rows.extend(
            build_transport_rows("Yan_Retina", fam_key, "Fovea_MajorClass_Primary", fam_c, skips, orig_ref)
        )

    path = write_cohort_csv("yan", out_rows)
    append_bh_rows(bh_rows)
    append_runtime("Yan_Retina", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}")
    summary = {
        "cohort": "Yan_Retina",
        "rows": len(out_rows),
        "P1m": count_both_pass(tier_results["P1m"]["df"]),
        "P0m": count_both_pass(tier_results["P0m"]["df"]),
        "skips_P1m": tier_results["P1m"]["skips"],
        "skips_P0m": tier_results["P0m"]["skips"],
    }
    save_cohort_payload("yan", {"summary": summary, "bh_rows": bh_rows})
    return summary


# ---------------------------------------------------------------------------
# Cohort 3: RPE bulk
# ---------------------------------------------------------------------------

def run_cohort_rpe(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    from scripts.run_followup_organ_08_eye_rpe_bulk import (
        _build_ensembl_map,
        _fetch_geo_donor_table,
        _load_tar_counts,
    )
    from src.followup_organs import counts_to_log2cpm, prepare_counts_matrix

    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("rpe")
    tar_path = DATA_ROOT / "tmp_organs" / "organ08_rpe" / "GSE159435_RAW.tar"
    if not tar_path.exists():
        raise FileNotFoundError(f"RPE input missing: {tar_path}")
    orig = pd.read_csv(RESULTS_DIR / "followup-organ-08-eye-rpe-gse159435-p1-signature.csv")
    prim = orig[orig["analysis"] == "Primary_BatchSexAdj"]
    assert len(prim) == 1, f"RPE primary rows {len(prim)} != 1"
    prim = prim.iloc[0]

    try:
        donor_meta = _fetch_geo_donor_table()
        donor_source = "live GEO fetch"
    except Exception as exc:  # noqa: BLE001 - network fallback to frozen table
        print(f"[RPE] Live GEO fetch failed ({exc}); using frozen donor table", flush=True)
        donor_meta = pd.read_csv(RESULTS_DIR / "followup-organ-08-eye-rpe-gse159435-donors.csv")
        donor_meta = donor_meta.set_index("donor_id", drop=False)
        donor_meta["age_per_decade"] = donor_meta["age"] / 10.0
        donor_meta["sex_male"] = (donor_meta["sex"].astype(str).str.upper() == "M").astype(float)
        donor_meta["batch_s2"] = (donor_meta["batch"].astype(str) == "s2").astype(float)
        donor_source = "frozen donors CSV"
    donor_meta["age_years"] = donor_meta["age"].astype(float)
    assert len(donor_meta) == 13, f"RPE donors {len(donor_meta)} != 13"
    assert float(donor_meta["age_years"].min()) >= MIN_AGE, "RPE has under-30 donors"

    sym_by_ens = _build_ensembl_map()
    raw_counts = _load_tar_counts(tar_path, donor_meta, sym_by_ens)
    clean_counts, _count_stats = prepare_counts_matrix(raw_counts, allow_fractional_expected=False)
    log2cpm_all, lib_all = counts_to_log2cpm(clean_counts)
    donor_meta["lib_size"] = lib_all.loc[donor_meta.index]
    donor_meta["log_lib"] = np.log1p(donor_meta["lib_size"])

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    tier_results: dict[str, Any] = {}
    orig_ref = pd.DataFrame(
        [
            {
                "group": "RPE_Bulk",
                "orig_p": float(prim["composite_pvalue"]),
                "orig_q": float(prim["q_bh_composite"]),
                "orig_emp_p": float(prim["empirical_null_pvalue"]),
                "orig_emp_q": float(prim["q_bh_empirical"]),
                "orig_effect": float(prim["composite_effect"]),
                "orig_n": int(prim["n_donors"]),
                "orig_excluded": 0,
            }
        ]
    ).set_index("group")
    for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
        t0 = time.monotonic()
        s0 = utc_now()
        fam_c, skips = run_family_c_mature(
            {"RPE_Bulk": log2cpm_all},
            {"RPE_Bulk": donor_meta},
            programs[prog_key],
            age_col="age_per_decade",
            candidate_covariates=("sex_male", "batch_s2", "log_lib"),
            min_donors=13,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=3,
            iterations=P1_ITERATIONS,
            seed=ORIG_NULL_SEED,
        )
        append_runtime("RPE_Bulk", f"family_c_{tier}", s0, t0, "ok",
                       f"tested={len(fam_c)}, skipped={len(skips)}, donor_source={donor_source}")
        tier_results[tier] = {"df": fam_c, "skips": skips}
        bh_rows.extend(bh_rows_from_family_c(fam_c, "RPE_Bulk", f"rpe_{fam_key}"))
        out_rows.extend(
            build_transport_rows("RPE_Bulk", fam_key, "Primary_BatchSexAdj", fam_c, skips, orig_ref)
        )

    path = write_cohort_csv("rpe", out_rows)
    append_bh_rows(bh_rows)
    append_runtime("RPE_Bulk", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}, donor_source={donor_source}")
    summary = {
        "cohort": "RPE_Bulk",
        "rows": len(out_rows),
        "donor_source": donor_source,
        "P1m": count_both_pass(tier_results["P1m"]["df"]),
        "P0m": count_both_pass(tier_results["P0m"]["df"]),
        "skips_P1m": tier_results["P1m"]["skips"],
        "skips_P0m": tier_results["P0m"]["skips"],
    }
    save_cohort_payload("rpe", {"summary": summary, "bh_rows": bh_rows})
    return summary


# ---------------------------------------------------------------------------
# Cohort 4: Ovary postmenopausal
# ---------------------------------------------------------------------------

def run_cohort_ovary(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    from src.followup_ovary import build_postmenopausal_design

    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("ovary")
    cache_path = DATA_ROOT / "ovary" / "GSE329417" / "interim" / "postmenopausal_donor_aggregates.npz"
    if not cache_path.exists():
        raise FileNotFoundError(f"Ovary input missing: {cache_path}")
    orig = pd.read_csv(RESULTS_DIR / "followup-ovary-postmenopausal-p1-models.csv")
    prim = orig[orig["model"] == "primary_depth_plex_clinical"].copy()
    assert len(prim) == 9, f"Ovary primary cell types {len(prim)} != 9"
    cell_types = sorted(prim["cell_type"].tolist())

    with np.load(cache_path, allow_pickle=False) as data:
        cache = {key: data[key] for key in data.files}
    gene_names = cache["gene_names"].astype(str)
    donors_all = cache["donor"].astype(str)
    assert len(set(donors_all)) == 28, f"Ovary donors {len(set(donors_all))} != 28"

    def unit_frame(indices: np.ndarray) -> pd.DataFrame:
        frame = pd.DataFrame(
            {
                "donor": cache["donor"][indices].astype(str),
                "cell_type": cache["cell_type"][indices].astype(str),
                "age_years": cache["age_years"][indices],
                "BMI": cache["BMI"][indices],
                "parity_live_births": cache["parity_live_births"][indices],
                "non_ovarian_abnormal": cache["non_ovarian_abnormal"][indices],
                "n_cells": cache["n_cells"][indices],
                "mean_log_umi": cache["mean_log_umi"][indices],
            }
        )
        for plex in range(1, 5):
            frame[f"plex_fraction_{plex}"] = (
                cache["unit_plex_counts"][indices, plex - 1] / frame["n_cells"].to_numpy(float)
            )
        return frame

    groups: dict[str, dict[str, Any]] = {}
    for ct in cell_types:
        indices = np.flatnonzero(
            (cache["cell_type"].astype(str) == ct) & (cache["n_cells"] >= 20)
        )
        indices = indices[np.argsort(cache["donor"][indices].astype(str))]
        frame = unit_frame(indices)
        assert float(frame["age_years"].min()) >= MIN_AGE, f"Ovary {ct} has under-30 donors"
        groups[ct] = {
            "matrix": (cache["log_sums"][indices] / cache["n_cells"][indices, None]).astype(float).T,
            "genes": gene_names,
            "meta": frame.set_index("donor"),
            "donors": frame["donor"].tolist(),
            "raw_counts": cache["raw_counts"][indices].T,
        }

    def ovary_design(m: pd.DataFrame) -> tuple[np.ndarray, list[str], list[str]]:
        design, names = build_postmenopausal_design(m.reset_index(), "primary_depth_plex_clinical")
        return design, names, []

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    tier_results: dict[str, Any] = {}
    orig_ref = pd.DataFrame(
        [
            {
                "group": r["cell_type"],
                "orig_p": float(r["pvalue"]),
                "orig_q": float(r["q_bh"]),
                "orig_emp_p": float(r["random_pvalue"]),
                "orig_emp_q": float(r["random_q_bh"]),
                "orig_effect": float(r["effect"]),
                "orig_n": int(r["n_donors"]),
                "orig_excluded": 0,
            }
            for _, r in prim.iterrows()
        ]
    ).set_index("group")
    for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
        t0 = time.monotonic()
        s0 = utc_now()
        df_c, skips = run_transport_meanlog(
            groups,
            programs[prog_key],
            ovary_design,
            min_donors=20,
            min_residual_df=2,
        )
        append_runtime("Ovary_Postmenopausal", f"family_c_{tier}", s0, t0, "ok",
                       f"tested={len(df_c)}, skipped={len(skips)}")
        if not df_c.empty:
            df_c, fam_bh = apply_bh_families(df_c, {g: f"ovary_{tier}" for g in df_c["group"]}, "Ovary_Postmenopausal")
            bh_rows.extend(fam_bh)
        tier_results[tier] = {"df": df_c, "skips": skips}
        out_rows.extend(
            build_transport_rows(
                "Ovary_Postmenopausal", fam_key, "primary_depth_plex_clinical", df_c, skips, orig_ref
            )
        )

    path = write_cohort_csv("ovary", out_rows)
    append_bh_rows(bh_rows)
    append_runtime("Ovary_Postmenopausal", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}")
    summary = {
        "cohort": "Ovary_Postmenopausal",
        "rows": len(out_rows),
        "P1m": count_both_pass(tier_results["P1m"]["df"]),
        "P0m": count_both_pass(tier_results["P0m"]["df"]),
        "skips_P1m": tier_results["P1m"]["skips"],
        "skips_P0m": tier_results["P0m"]["skips"],
    }
    save_cohort_payload("ovary", {"summary": summary, "bh_rows": bh_rows})
    return summary


# ---------------------------------------------------------------------------
# Cohort 5: PBMC (Zhu GSE213516 and Thomson GSE214546)
# ---------------------------------------------------------------------------

def run_cohort_pbmc(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("pbmc")
    cohorts = ["GSE213516", "GSE214546"]
    for cohort in cohorts:
        p = DATA_ROOT / "pbmc" / "interim" / cohort / "donor-aggregates.npz"
        if not p.exists():
            raise FileNotFoundError(f"PBMC input missing: {p}")
    orig = pd.read_csv(RESULTS_DIR / "followup-pbmc-p1-models.csv")
    prim = orig[(orig["design"] == "strict_group") & (orig["model"] == "primary")].copy()
    assert len(prim) == 36, f"PBMC strict primary rows {len(prim)} != 36"

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    tier_summaries: dict[str, Any] = {}
    for cohort in cohorts:
        cache = dict(np.load(DATA_ROOT / "pbmc" / "interim" / cohort / "donor-aggregates.npz", allow_pickle=True))
        symbols = cache["symbol"].astype(str)
        unique_symbols = ~np.asarray(pd.Index(symbols).duplicated(keep=False))
        donor_ids = cache["donor"].astype(str)
        ages = cache["age"].astype(float)
        sexes = cache["sex"].astype(str)
        donor_age = {}
        donor_sex = {}
        for d, a, s in zip(donor_ids, ages, sexes):
            donor_age[str(d)] = float(a)
            donor_sex[str(d)] = str(s)
        mature_donors = [d for d in sorted(donor_age) if 30 <= donor_age[d] <= 39 or donor_age[d] >= 60]
        assert mature_donors, f"PBMC {cohort}: no mature donors"

        orig_cohort = prim[prim["cohort"] == cohort]
        cell_types = sorted(orig_cohort["cell_type"].tolist())
        groups: dict[str, dict[str, Any]] = {}
        excluded: dict[str, int] = {}
        for ct in cell_types:
            mask = (cache["cell_type"].astype(str) == ct) & (cache["n_cells"] >= 20)
            indices = np.flatnonzero(mask)
            frame_donors = [str(d) for d in cache["donor"][indices]]
            keep = [i for i, d in zip(indices, frame_donors) if d in mature_donors]
            # Deterministic order: sort kept donor rows by donor id.
            kept_pairs = sorted((str(cache["donor"][i]), i) for i in keep)
            order = np.asarray([i for _, i in kept_pairs], dtype=int)
            donors = [str(cache["donor"][i]) for i in order]
            old = np.asarray([1.0 if donor_age[d] >= 60 else 0.0 for d in donors])
            sex_male = np.asarray([1.0 if donor_sex[d].lower() == "male" else 0.0 for d in donors])
            meta = pd.DataFrame(
                {
                    "age_years": [donor_age[d] for d in donors],
                    "old": old,
                    "sex_male": sex_male,
                    "mean_log_umi": cache["mean_log_umi"][order],
                },
                index=donors,
            )
            # Donors in the original strict design but absent from the mature design.
            orig_row = orig_cohort[orig_cohort["cell_type"] == ct].iloc[0]
            orig_strict_donors = set(
                d for d in donor_age
                if (donor_age[d] <= 40 or donor_age[d] > 60)
                and d in [str(x) for x in cache["donor"][indices]]
            )
            excluded[ct] = len(orig_strict_donors - set(donors))
            groups[ct] = {
                "matrix": (cache["log_sums"][order] / cache["n_cells"][order, None]).astype(float).T,
                "genes": symbols,
                "meta": meta,
                "donors": donors,
                "raw_counts": None,
                "detectable": unique_symbols,
            }

        def pbmc_design(m: pd.DataFrame) -> tuple[np.ndarray, list[str], list[str]]:
            cols = [np.ones(len(m)), m["old"].to_numpy(float)]
            names = ["intercept", "old"]
            omitted: list[str] = []
            if m["sex_male"].nunique() == 2:
                if len(m) - (len(cols) + 1) < 2:
                    omitted.append("sex_male(df_limit)")
                else:
                    cols.append(m["sex_male"].to_numpy(float))
                    names.append("sex_male")
            else:
                omitted.append("sex_male(constant_or_na)")
            design = np.column_stack(cols)
            if np.linalg.matrix_rank(design) != design.shape[1]:
                raise ValueError("PBMC mature design is not full rank")
            return design, names, omitted

        # Restrict program detection to unique symbols, as in the original.
        orig_ref = pd.DataFrame(
            [
                {
                    "group": r["cell_type"],
                    "orig_p": float(r["pvalue"]),
                    "orig_q": float(r["q_bh"]),
                    "orig_emp_p": float(r["random_pvalue"]),
                    "orig_emp_q": float(r["random_q_bh"]),
                    "orig_effect": float(r["effect"]),
                    "orig_n": int(r["n_donors"]),
                    "orig_excluded": int(excluded[r["cell_type"]]),
                }
                for _, r in orig_cohort.iterrows()
            ]
        ).set_index("group")
        for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
            t0 = time.monotonic()
            s0 = utc_now()
            # run_transport_meanlog detects program genes by presence; apply
            # run_transport_meanlog detects program genes by presence; the
            # detectable mask restricts detection to unique symbols, as in
            # the original. The null pool keeps every matrix row.
            df_c, skips = run_transport_meanlog(
                groups,
                programs[prog_key],
                pbmc_design,
                min_donors=6,
                min_residual_df=2,
            )
            append_runtime(f"PBMC_{cohort}", f"family_c_{tier}", s0, t0, "ok",
                           f"tested={len(df_c)}, skipped={len(skips)}")
            if not df_c.empty:
                df_c, fam_bh = apply_bh_families(df_c, {g: f"pbmc_{cohort}_{tier}" for g in df_c["group"]}, f"PBMC_{cohort}")
                bh_rows.extend(fam_bh)
            tier_summaries[f"{cohort}_{tier}"] = {"df": df_c, "skips": skips}
            out_rows.extend(
                build_transport_rows(
                    f"PBMC_{cohort}", fam_key, "mature_30_39_vs_60_plus_primary",
                    df_c, skips, orig_ref, effect_scale_note="old_minus_young_30_39_vs_60plus",
                )
            )

    path = write_cohort_csv("pbmc", out_rows)
    append_bh_rows(bh_rows)
    append_runtime("PBMC", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}")
    summary = {
        "cohort": "PBMC",
        "rows": len(out_rows),
        "GSE213516_P1m": count_both_pass(tier_summaries["GSE213516_P1m"]["df"]),
        "GSE213516_P0m": count_both_pass(tier_summaries["GSE213516_P0m"]["df"]),
        "GSE214546_P1m": count_both_pass(tier_summaries["GSE214546_P1m"]["df"]),
        "GSE214546_P0m": count_both_pass(tier_summaries["GSE214546_P0m"]["df"]),
        "skips": {k: v["skips"] for k, v in tier_summaries.items()},
    }
    save_cohort_payload("pbmc", {"summary": summary, "bh_rows": bh_rows})
    return summary


# ---------------------------------------------------------------------------
# Cohort 6: Lake kidney (scCv3 and snCv3)
# ---------------------------------------------------------------------------

def run_cohort_lake(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("lake")
    assays = ["scCv3", "snCv3"]
    for assay in assays:
        p = DATA_ROOT / "lake_kidney" / "interim" / f"{assay}_donor_aggregates.npz"
        if not p.exists():
            raise FileNotFoundError(f"Lake input missing: {p}")
    orig = pd.read_csv(RESULTS_DIR / "followup-lake-kidney-p1-models.csv")
    prim = orig[(orig["donor_set"] == "primary_adult") & (orig["model"] == "primary_age_sex")].copy()
    assert len(prim) == 22, f"Lake primary rows {len(prim)} != 22"

    out_rows: list[dict[str, Any]] = []
    bh_rows: list[dict[str, Any]] = []
    tier_summaries: dict[str, Any] = {}
    for assay in assays:
        cache = dict(np.load(DATA_ROOT / "lake_kidney" / "interim" / f"{assay}_donor_aggregates.npz", allow_pickle=True))
        genes = cache["genes"].astype(str)
        orig_assay = prim[prim["assay_release"] == assay]
        cell_types = sorted(orig_assay["cell_type"].tolist())
        groups: dict[str, dict[str, Any]] = {}
        excluded: dict[str, int] = {}
        for ct in cell_types:
            indices = np.flatnonzero(cache["cell_type"].astype(str) == ct)
            keep = indices[
                (cache["n_cells"][indices] >= 20)
                & (cache["primary"][indices])
                & (cache["age"][indices] >= MIN_AGE)
            ]
            order = keep[np.argsort(cache["donor_id"][keep].astype(str))]
            donors = cache["donor_id"][order].astype(str).tolist()
            ages = cache["age"][order].astype(float)
            sex_male = (pd.Series(cache["sex"][order].astype(str)).str.lower().to_numpy() == "male").astype(float)
            meta = pd.DataFrame(
                {"age_years": ages, "sex_male": sex_male, "mean_log_umi": cache["mean_log_umi"][order]},
                index=donors,
            )
            orig_row = orig_assay[orig_assay["cell_type"] == ct].iloc[0]
            orig_n = int(orig_row["n_donors"])
            excluded[ct] = orig_n - len(donors)
            assert excluded[ct] >= 0, f"Lake {assay} {ct}: mature n exceeds original n"
            groups[ct] = {
                "matrix": (cache["log_sums"][order] / cache["n_cells"][order, None]).astype(float).T,
                "genes": genes,
                "meta": meta,
                "donors": donors,
                "raw_counts": None,
            }

        def lake_design(m: pd.DataFrame) -> tuple[np.ndarray, list[str], list[str]]:
            cols = [np.ones(len(m)), m["age_years"].to_numpy(float)]
            names = ["intercept", "age_years"]
            omitted: list[str] = []
            if m["sex_male"].nunique() == 2:
                if len(m) - (len(cols) + 1) < 2:
                    omitted.append("sex_male(df_limit)")
                else:
                    cols.append(m["sex_male"].to_numpy(float))
                    names.append("sex_male")
            else:
                omitted.append("sex_male(constant_or_na)")
            design = np.column_stack(cols)
            if np.linalg.matrix_rank(design) != design.shape[1]:
                raise ValueError("Lake mature design is not full rank")
            return design, names, omitted

        orig_ref = pd.DataFrame(
            [
                {
                    "group": r["cell_type"],
                    "orig_p": float(r["pvalue"]),
                    "orig_q": float(r["q_bh"]),
                    "orig_emp_p": float(r["random_pvalue"]),
                    "orig_emp_q": float(r["random_q_bh"]),
                    "orig_effect": float(r["slope_per_year"]),
                    "orig_n": int(r["n_donors"]),
                    "orig_excluded": int(excluded[r["cell_type"]]),
                }
                for _, r in orig_assay.iterrows()
            ]
        ).set_index("group")
        for tier, prog_key, fam_key in [("P1m", "P1m", "transport_P1m"), ("P0m", "P0m", "transport_P0m")]:
            t0 = time.monotonic()
            s0 = utc_now()
            df_c, skips = run_transport_meanlog(
                groups,
                programs[prog_key],
                lake_design,
                min_donors=6,
                min_residual_df=2,
            )
            append_runtime(f"Lake_Kidney_{assay}", f"family_c_{tier}", s0, t0, "ok",
                           f"tested={len(df_c)}, skipped={len(skips)}")
            if not df_c.empty:
                df_c, fam_bh = apply_bh_families(df_c, {g: f"lake_{assay}_{tier}" for g in df_c["group"]}, f"Lake_Kidney_{assay}")
                bh_rows.extend(fam_bh)
            tier_summaries[f"{assay}_{tier}"] = {"df": df_c, "skips": skips}
            out_rows.extend(
                build_transport_rows(
                    f"Lake_Kidney_{assay}", fam_key, "primary_adult_primary_age_sex",
                    df_c, skips, orig_ref, effect_scale_note="per_year",
                )
            )

    path = write_cohort_csv("lake", out_rows)
    append_bh_rows(bh_rows)
    append_runtime("Lake_Kidney", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}")
    summary = {
        "cohort": "Lake_Kidney",
        "rows": len(out_rows),
        "scCv3_P1m": count_both_pass(tier_summaries["scCv3_P1m"]["df"]),
        "scCv3_P0m": count_both_pass(tier_summaries["scCv3_P0m"]["df"]),
        "snCv3_P1m": count_both_pass(tier_summaries["snCv3_P1m"]["df"]),
        "snCv3_P0m": count_both_pass(tier_summaries["snCv3_P0m"]["df"]),
        "skips": {k: v["skips"] for k, v in tier_summaries.items()},
    }
    save_cohort_payload("lake", {"summary": summary, "bh_rows": bh_rows})
    return summary


# ---------------------------------------------------------------------------
# Cohort 7: Perez muscle (continuous 30-plus whole-nucleus gene model only)
# ---------------------------------------------------------------------------

def run_cohort_muscle(programs: dict[str, dict[str, list[str]]]) -> dict[str, Any]:
    from scipy import stats as scipy_stats
    from scripts.run_followup_perez import donor_path, read_donor

    t_start = time.monotonic()
    s_start = utc_now()
    wait_for_memory("muscle")
    manifest = json.loads((REPO_ROOT / "data" / "manifest" / "GSE167186.json").read_text())
    assert len(manifest["donors"]) == 17, f"Perez manifest donors {len(manifest['donors'])} != 17"
    for donor in manifest["donors"]:
        if not donor_path(donor).is_file():
            raise FileNotFoundError(f"Perez muscle input missing: {donor_path(donor)}")

    orig_genes = pd.read_csv(RESULTS_DIR / "followup-muscle-perez-gene-models.csv")
    orig_exact = orig_genes[orig_genes["model"] == "exact_age_10y_sensitivity"].copy()
    assert not orig_exact.empty, "Perez original exact-age gene rows missing"
    orig_q_col = "q_bh" if "q_bh" in orig_exact.columns else "q"
    orig_nsig = int((orig_exact[orig_q_col] < ALPHA).sum())
    orig_n_tested = int(len(orig_exact))
    orig_min = orig_exact.loc[orig_exact["pvalue"].idxmin()]

    # Read donors in gsm-sorted order so count rows align with the donor table.
    summaries, count_rows = [], []
    expected_gene_ids = expected_gene_names = None
    for donor in sorted(manifest["donors"], key=lambda d: d["gsm"]):
        summary, gene_sums, gene_ids, gene_names = read_donor(donor)
        if expected_gene_ids is None:
            expected_gene_ids, expected_gene_names = gene_ids, gene_names
        else:
            assert np.array_equal(gene_ids, expected_gene_ids), "Perez feature order differs"
            assert np.array_equal(gene_names, expected_gene_names), "Perez feature names differ"
        summaries.append(summary)
        count_rows.append(gene_sums)
    donors = pd.DataFrame(summaries).reset_index(drop=True)
    assert len(donors) == 17, f"Perez donors {len(donors)} != 17"
    mature_mask = donors["age_years"].to_numpy() >= MIN_AGE
    assert int((~mature_mask).sum()) == 6, "Perez young donors != 6"
    mature = donors[mature_mask].reset_index(drop=True)
    assert len(mature) == 11, f"Perez mature donors {len(mature)} != 11"
    assert float(mature["age_years"].min()) == 70.0 and float(mature["age_years"].max()) == 90.0

    counts = np.vstack([count_rows[i] for i in np.flatnonzero(mature_mask)]).astype(float)
    totals = counts.sum(axis=1)
    expression = np.log1p(counts * (1e6 / totals[:, None]))
    gene_filter = (counts.sum(axis=0) >= 10) & ((counts > 0).sum(axis=0) >= 2)
    design = np.column_stack([
        np.ones(len(mature)),
        mature["age_years"].to_numpy(float) / 10,
        mature["mean_log_umi"].to_numpy(float),
    ])
    assert np.linalg.matrix_rank(design) == 3, "Perez mature continuous design is not full rank"
    assert len(mature) - 3 >= 2, "Perez mature design lacks residual df"
    fitted = fit_ols(expression[:, gene_filter], design)
    testable = np.isfinite(fitted["pvalue"]) & np.isfinite(fitted["se"]) & (fitted["se"] > 0)
    gene_indices = np.flatnonzero(gene_filter)[testable]
    pvals = fitted["pvalue"][testable]
    q_glob, _q_sm, bh_diff = bh_checked(pvals)
    assert bh_diff <= 1e-12, f"Perez muscle BH diff {bh_diff}"
    n_sig = int((q_glob < ALPHA).sum())
    effects = fitted["effect"][testable]
    ses = fitted["se"][testable]
    df_val = int(fitted["residual_df"][testable][0])
    mde_factor = float(scipy_stats.t.ppf(0.975, df_val) + scipy_stats.t.ppf(0.80, df_val))
    mde_vec = mde_factor * ses
    top_pos = int(np.argmin(pvals))
    top_gene = str(expected_gene_names[gene_indices[top_pos]])
    verd, _vnote = verdict_standard(orig_nsig > 0, n_sig > 0)

    out_rows = [
        make_row(
            "Perez_Muscle",
            "whole_nucleus",
            "genes",
            int(len(mature)),
            6,
            float(orig_min["pvalue"]),
            float(orig_min[orig_q_col]),
            float(pvals[top_pos]),
            float(q_glob[top_pos]),
            None,
            None,
            verd,
            f"n_genes={int(len(pvals))};orig_sig={orig_nsig};mature_sig={n_sig};"
            f"top={top_gene}(eff={float(effects[top_pos]):+.4f},p={float(pvals[top_pos]):.3g},q={float(q_glob[top_pos]):.3g});"
            f"continuous_30plus_70_to_90;composition_confounded",
            arm="exact_age_10y_continuous_30plus",
            residual_df=df_val,
            original_effect=None,
            mature_effect=float(effects[top_pos]),
            mature_se=float(ses[top_pos]),
            null_sd=None,
            mde80=float(np.median(mde_vec)),
            mde80_method="median_t_quantile_se",
            covariates="age_per_10y;mean_log_umi",
            omitted_covariates="none",
            n_genes_tested=int(len(pvals)),
            n_sig_q05=n_sig,
            original_empirical_p=None,
            original_empirical_q=None,
            orig_transport_class=None,
            mature_transport_class=None,
            n_up_detected=None,
            n_down_detected=None,
        )
    ]
    # Transport stays blocked: no released per-nucleus annotations.
    for fam_key in ("transport_P1m", "transport_P0m"):
        out_rows.append(
            make_row(
                "Perez_Muscle",
                "whole_nucleus",
                fam_key,
                0,
                17,
                None,
                None,
                None,
                None,
                None,
                None,
                "UNTESTABLE",
                "reason=blocked_no_released_annotations;no whole-nucleus transport claim",
                arm="exact_age_10y_continuous_30plus",
                residual_df=None,
                original_effect=None,
                mature_effect=None,
                mature_se=None,
                null_sd=None,
                mde80=None,
                mde80_method=None,
                covariates=None,
                omitted_covariates=None,
                n_genes_tested=None,
                n_sig_q05=None,
                original_empirical_p=None,
                original_empirical_q=None,
                orig_transport_class=None,
                mature_transport_class="UNTESTABLE",
                n_up_detected=None,
                n_down_detected=None,
            )
        )

    path = write_cohort_csv("muscle", out_rows)
    append_runtime("Perez_Muscle", "cohort_total", s_start, t_start, "ok",
                   f"rows={len(out_rows)}, csv={path.name}, genes_tested={len(pvals)}, n_sig={n_sig}")
    summary = {
        "cohort": "Perez_Muscle",
        "rows": len(out_rows),
        "n_donors_mature": int(len(mature)),
        "n_genes_tested": int(len(pvals)),
        "n_sig_q05": n_sig,
        "orig_n_genes_tested": orig_n_tested,
        "orig_n_sig_q05": orig_nsig,
        "top_gene": top_gene,
        "top_effect": float(effects[top_pos]),
        "top_p": float(pvals[top_pos]),
        "top_q": float(q_glob[top_pos]),
        "verdict": verd,
        "bh_diff": float(bh_diff),
    }
    save_cohort_payload("muscle", {"summary": summary, "bh_rows": []})
    return summary


# ---------------------------------------------------------------------------
# Finalize + main
# ---------------------------------------------------------------------------

COHORT_KEYS = ["ts2", "yan", "rpe", "ovary", "pbmc", "lake", "muscle"]


FAMILY_BASE = {
    "Ovary_Postmenopausal": "ovary",
    "PBMC_GSE213516": "pbmc_GSE213516",
    "PBMC_GSE214546": "pbmc_GSE214546",
    "Lake_Kidney_scCv3": "lake_scCv3",
    "Lake_Kidney_snCv3": "lake_snCv3",
    "Yan_Retina": "yan_transport",
    "RPE_Bulk": "rpe_transport",
}


def rebuild_bh_check_from_csvs() -> pd.DataFrame:
    """Rebuild the BH-check CSV from the cohort CSVs (the verified artifacts).

    Some cohort runners wrote tier-less family keys (ovary, pbmc, lake,
    ts2_nk), so the incremental append deduped P1m rows away when P0m rows
    with the same (cohort, family, unit) arrived. Recompute BH fresh from
    the CSV raw p values per declared family and assert the recomputed q
    equals the CSV mature q for every row. Family names are tier-distinct.
    """
    repaired: list[dict[str, Any]] = []
    n_tested = 0
    for key in COHORT_KEYS:
        df = pd.read_csv(RESULTS_DIR / f"{PREFIX}-{key}.csv")
        trans = df[df["endpoint_family"].isin(("transport_P1m", "transport_P0m"))]
        trans = trans[trans["mature_q"].notna()].copy()
        n_tested += len(trans)

        def family_id(cohort_label: str) -> str:
            if cohort_label == "TS2_Blood_NK":
                return "ts2_nk"
            if str(cohort_label).startswith("TS2_"):
                return "ts2_transport"
            return str(cohort_label)

        trans["family_id"] = trans["cohort"].astype(str).map(family_id)
        for (fam_id, fam_key), sub in trans.groupby(
            ["family_id", "endpoint_family"], sort=False
        ):
            tier = "P1m" if fam_key == "transport_P1m" else "P0m"
            if fam_id == "ts2_nk":
                base = "ts2_nk"
                units = sub["cell_type"].astype(str).tolist()
                cohort_out = "TS2"
            elif fam_id == "ts2_transport":
                base = "ts2_transport"
                units = [
                    f"{str(c)[len('TS2_'):]}::{ct}"
                    for c, ct in zip(sub["cohort"].astype(str), sub["cell_type"].astype(str))
                ]
                cohort_out = "TS2"
            else:
                base = FAMILY_BASE[fam_id]
                units = sub["cell_type"].astype(str).tolist()
                cohort_out = fam_id
            for p_col, q_col, suffix in [
                ("mature_p", "mature_q", "composite"),
                ("mature_empirical_p", "mature_empirical_q", "empirical"),
            ]:
                q_exp, q_sm, diff = bh_checked(sub[p_col].to_numpy(float))
                assert np.allclose(q_exp, sub[q_col].to_numpy(float), atol=1e-12), (
                    f"BH recomputation disagrees with CSV q for {key} "
                    f"{fam_id} {fam_key} {suffix}"
                )
                assert diff <= 1e-12, f"BH explicit-vs-statsmodels diff {diff} too large"
                for i, (_, r) in enumerate(sub.iterrows()):
                    repaired.append(
                        {
                            "cohort": cohort_out,
                            "family": f"{base}_{tier}_{suffix}_{len(sub)}",
                            "unit": units[i],
                            "family_size": int(len(sub)),
                            "raw_p": float(r[p_col]),
                            "q_bh_explicit": float(q_exp[i]),
                            "q_bh_statsmodels": float(q_sm[i]),
                            "abs_diff": float(abs(q_exp[i] - q_sm[i])),
                            "family_max_abs_diff": float(diff),
                        }
                    )
    out = pd.DataFrame(repaired)
    assert len(out) == 2 * n_tested, (
        f"Repaired BH rows {len(out)} != 2 x tested transport rows {2 * n_tested}"
    )
    assert not out.duplicated(subset=["cohort", "family", "unit"]).any(), (
        "Repaired BH-check CSV still has duplicate (cohort, family, unit) keys"
    )
    write_clean_csv(out, BH_CHECK_CSV)
    return out


def run_finalize() -> None:
    t0 = time.monotonic()
    s0 = utc_now()
    print("\n================ FINALIZING BATCH C OUTPUTS ================", flush=True)
    summaries: dict[str, Any] = {}
    bh_all: list[dict[str, Any]] = []
    for key in COHORT_KEYS:
        payload = load_cohort_payload(key)
        summaries[key] = payload["summary"]
        for r in payload.get("bh_rows", []):
            r = dict(r)
            r["_cohort_key"] = key
            bh_all.append(r)

    # Rebuild the BH-check CSV from the cohort CSVs (tier-disambiguated);
    # this repairs the tier-less family keys written by ovary/pbmc/lake/nk runs.
    rebuilt = rebuild_bh_check_from_csvs()
    print(f"[FINALIZE] rebuilt BH-check CSV: {len(rebuilt)} rows", flush=True)

    # Verify every cohort CSV: exact Batch A schema, no forbidden dashes.
    # Verify every cohort CSV: exact Batch A schema, no forbidden dashes.
    for key in COHORT_KEYS:
        p = RESULTS_DIR / f"{PREFIX}-{key}.csv"
        if not p.exists():
            raise FileNotFoundError(f"Missing cohort CSV: {p}")
        df = pd.read_csv(p)
        assert list(df.columns) == REQUIRED_COHORT_COLS + EXTRA_COLS, f"Schema mismatch in {p}"
        assert_no_forbidden_dashes(p.read_text(encoding="utf-8"), label=str(p))

    # Verify BH-check rows already written incrementally: max diff and
    # cross-check of q values against the cohort CSVs.
    bh_df = pd.read_csv(BH_CHECK_CSV)
    assert not bh_df.empty, "BH-check CSV is empty"
    max_diff = float(bh_df["abs_diff"].max())
    max_fam_diff = float(bh_df["family_max_abs_diff"].max())
    assert max_diff <= 1e-12 and max_fam_diff <= 1e-12, (
        f"BH check failed: max_diff={max_diff}, max_fam_diff={max_fam_diff}"
    )
    # Cross-check: every cohort CSV mature q must match a BH-check row with
    # the same raw p value (composite and empirical separately).
    for key in COHORT_KEYS:
        df = pd.read_csv(RESULTS_DIR / f"{PREFIX}-{key}.csv")
        for _, r in df.iterrows():
            if r["endpoint_family"] not in ("transport_P1m", "transport_P0m"):
                continue
            if pd.isna(r["mature_q"]):
                continue
            for p_col, q_col, suffix in [
                ("mature_p", "mature_q", "composite"),
                ("mature_empirical_p", "mature_empirical_q", "empirical"),
            ]:
                raw_p = float(r[p_col])
                q_val = float(r[q_col])
                hits = bh_df[
                    (bh_df["family"].str.contains(f"_{suffix}_"))
                    & (np.abs(bh_df["raw_p"] - raw_p) < 1e-12)
                    & (np.abs(bh_df["q_bh_explicit"] - q_val) < 1e-12)
                ]
                assert len(hits) >= 1, (
                    f"BH cross-check failed for {key} {r['cell_type']} {r['endpoint_family']} "
                    f"{suffix}: p={raw_p}, q={q_val}"
                )
    (TMP_WORK_DIR / "summaries.json").write_text(
        json.dumps(summaries, indent=2, default=str) + "\n", encoding="utf-8"
    )
    append_runtime("Batch_C_All", "finalize_and_verify", s0, t0, "ok",
                   f"cohorts={len(summaries)}, bh_rows={len(bh_df)}, max_bh_diff={max_diff:.3e}")
    print(
        f"[FINALIZE OK] cohorts={len(summaries)}, bh_check={len(bh_df)} rows "
        f"(max_diff={max_diff:.3e}, max_fam_diff={max_fam_diff:.3e})",
        flush=True,
    )
    for key in COHORT_KEYS:
        print(f"  {key}: {json.dumps(summaries[key], default=str)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cohort",
        required=True,
        choices=["ts2", "yan", "rpe", "ovary", "pbmc", "lake", "muscle", "finalize", "all"],
    )
    args = parser.parse_args()

    if args.cohort == "finalize":
        run_finalize()
        return 0

    programs = load_verified_mature_programs()
    runners = {
        "ts2": run_cohort_ts2,
        "yan": run_cohort_yan,
        "rpe": run_cohort_rpe,
        "ovary": run_cohort_ovary,
        "pbmc": run_cohort_pbmc,
        "lake": run_cohort_lake,
        "muscle": run_cohort_muscle,
    }
    if args.cohort == "all":
        for key in COHORT_KEYS:
            summary = runners[key](programs)
            print(f"[{key}] done: {json.dumps(summary, default=str)}", flush=True)
            gc.collect()
        run_finalize()
        return 0
    summary = runners[args.cohort](programs)
    print(f"[{args.cohort}] done: {json.dumps(summary, default=str)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
