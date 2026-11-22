#!/usr/bin/env python3
"""Mature M10 PsychAD pipeline at 30-plus with P1m/P0m.

This is the M10-pipeline rerun the Phase 1 audit prescribes (audit line 188)
that Phase 2 Item 2 did NOT produce. Item 2 (commit 3f5b57a) already ran the
core PsychAD transport with P1m/P0m at 30-plus: the fully adjusted per-cell-type
slopes, the 2000-draw matched null, and the strict P0m rescue comparison.
Those results are reused, never recomputed as new claims.

What this script adds, mirroring the original run_m10_pipeline.py outputs:
  1. Per-cell-type unadjusted slopes and raw-count sensitivity slopes
     (m10-celltype-regressions equivalent).
  2. Stratified fits at 30-plus: age >= 40, HBCC-only, MSSM-only
     (m10-confounder-audit equivalent). The adult >= 20 and developmental
     < 20 strata are not applicable: every fitted donor is aged >= 30.
  3. Whole-cohort pooled summary at 30-plus (m10-cohort-summary equivalent).
  4. A 200-draw matched null as HISTORICAL SENSITIVITY ONLY, same seed and
     procedure as Item 2's 2000-draw null, compared against it. The 2000-draw
     null stays primary; the 200-draw run exists only to show the Monte Carlo
     floor effect (empirical p floor 1/201 vs 1/2001).

The fully adjusted slope in this script is asserted to match the committed
Item 2 transport CSV exactly (rtol 1e-10, atol 1e-12), the same tolerance
Item 2's own null script used. The controlled DUA input stays on beast and
is never copied to git or uploaded anywhere.
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_m10_pipeline import ols as m10_ols
from scripts.run_mature_p1m_transport import (
    MIN_AGE,
    RDS,
    RDS_MD5,
    TIERS,
    TRANSPORT_CSV,
    load_programs,
    md5sum,
    psychad_design,
)
from scripts.run_reviewer_phase1_signature_null import (
    random_signature_scores as m10_random_scores,
)
from src.followup_muscle import matched_null_effect_pvalue

RESULTS = REPO / "results"
FIGURES = REPO / "figures"
REGRESSIONS_CSV = RESULTS / "followup-mature-psychad-m10-regressions.csv"
AUDIT_CSV = RESULTS / "followup-mature-psychad-m10-confounder-audit.csv"
SUMMARY_CSV = RESULTS / "followup-mature-psychad-m10-cohort-summary.csv"
NULL200_CSV = RESULTS / "followup-mature-psychad-m10-null200.csv"
RUNTIME_CSV = RESULTS / "followup-mature-psychad-m10-runtime.csv"
NULL2000_CSV = RESULTS / "followup-mature-p1m-transport-null.csv"
N_DRAWS_200 = 200
SEED = 20261105
MODEL = "direct-runner-muse-subagent"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def ols_fit(y: np.ndarray, columns: list[np.ndarray]):
    """OLS with intercept + age as the second column; None if ineligible."""
    X = np.column_stack(columns)
    if X.shape[0] <= X.shape[1] or np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    try:
        return m10_ols(y, X)
    except np.linalg.LinAlgError:
        return None


def values(fit, prefix: str) -> dict:
    return {
        f"slope_{prefix}": fit["slope"] if fit else np.nan,
        f"se_{prefix}": fit["se"] if fit else np.nan,
        f"p_{prefix}": fit["p"] if fit else np.nan,
    }


def main() -> int:
    warnings.filterwarnings("ignore")
    t0 = time.time()
    stages: list[dict] = []

    digest = md5sum(RDS)
    assert digest == RDS_MD5, f"PsychAD MD5 mismatch: {digest}"
    print(f"PsychAD M10 mature: MD5 verified {digest}", flush=True)

    programs = load_programs()
    transport = pd.read_csv(TRANSPORT_CSV)
    psy = transport[transport["cohort"] == "DLPFC_Lifespan_Atlas"].copy()
    expected_units = sorted(psy[psy["gene_set"] == "P1m"]["unit"].tolist())
    assert len(expected_units) == 27, f"expected 27 PsychAD units, got {len(expected_units)}"
    expected_slope = {
        (row["unit"], row["gene_set"]): row["composite_slope"]
        for _, row in psy.iterrows()
    }
    null2000 = pd.read_csv(NULL2000_CSV)
    null2000 = null2000[null2000["cohort"] == "DLPFC_Lifespan_Atlas"].copy()
    expected_null = {
        (row["unit"], row["gene_set"]): (row["empirical_pvalue"], row["empirical_q"])
        for _, row in null2000.iterrows()
    }
    assert len(expected_null) == 54, f"expected 54 PsychAD null rows, got {len(expected_null)}"

    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    donor_ids = [str(donor) for donor in cd.rownames]
    metadata = {
        name: np.asarray(cd.listData[name]).ravel()
        for name in ("Age", "Sex", "Brain_bank", "PMI", "TechPC1", "TechPC2", "TechPC3")
    }
    age_all = metadata["Age"].astype(float)
    mature = age_all >= MIN_AGE
    assert int(mature.sum()) == 201, f"expected 201 mature donors, got {int(mature.sum())}"
    cell_types = sorted(obj.assays.data.listData.keys())
    assert cell_types == expected_units, "RDS cell types do not match the Item 2 PsychAD family"

    regression_rows: list[dict] = []
    audit_rows: list[dict] = []
    null_rows: list[dict] = []
    donor_scores: dict[str, dict[str, list[float]]] = {tier: {} for tier in TIERS}

    for tier in TIERS:
        program = programs[tier]
        rng = np.random.default_rng(SEED)
        t_tier = time.time()
        for cell_type in cell_types:
            source = obj.assays.data.listData[cell_type]
            genes = np.array([str(gene) for gene in source.Dimnames[0]])
            assert [str(donor) for donor in source.Dimnames[1]] == donor_ids
            matrix = sp.csc_matrix(
                (np.asarray(source.x), np.asarray(source.i), np.asarray(source.p)),
                shape=tuple(source.Dim),
            )
            library = np.asarray(matrix.sum(axis=0)).ravel()
            keep = mature & (library > 0)
            n_donors = int(keep.sum())
            n_zero = int((mature & (library == 0)).sum())
            mat = matrix[:, keep]
            normalized = mat @ sp.diags(1e6 / library[keep])
            normalized = normalized.tocsc()
            normalized.data = np.log1p(normalized.data)
            raw = mat.tocsc()
            raw.data = np.log1p(raw.data.astype(np.float64))
            lookup = {gene: i for i, gene in enumerate(genes)}
            up = np.array([lookup[g] for g in program["up"] if g in lookup], dtype=int)
            down = np.array([lookup[g] for g in program["down"] if g in lookup], dtype=int)
            score = (
                np.asarray(normalized[up, :].mean(axis=0)).ravel()
                - np.asarray(normalized[down, :].mean(axis=0)).ravel()
            )
            raw_score = (
                np.asarray(raw[up, :].mean(axis=0)).ravel()
                - np.asarray(raw[down, :].mean(axis=0)).ravel()
            )
            coverage = (len(up) + len(down)) / max(len(set(program["up"]) | set(program["down"])), 1)

            a = age_all[keep]
            male = (metadata["Sex"][keep] == "Male").astype(float)
            hbcc = (metadata["Brain_bank"][keep] == "HBCC").astype(float)
            pmi = metadata["PMI"][keep].astype(float)
            pc1 = metadata["TechPC1"][keep].astype(float)
            pc2 = metadata["TechPC2"][keep].astype(float)
            pc3 = metadata["TechPC3"][keep].astype(float)
            ones = np.ones(len(a))
            unadjusted = ols_fit(score, [ones, a])
            covar = ols_fit(score, [ones, a, male, hbcc, pmi, pc1, pc2, pc3])
            raw_fit = ols_fit(raw_score, [ones, a])

            # The fully adjusted slope must reproduce the committed Item 2
            # transport slope exactly; same tolerance Item 2's null used.
            assert covar is not None, f"{tier} {cell_type}: fully adjusted fit ineligible"
            assert np.isclose(
                covar["slope"], expected_slope[(cell_type, tier)], rtol=1e-10, atol=1e-12
            ), f"{tier} {cell_type}: slope {covar['slope']} != transport {expected_slope[(cell_type, tier)]}"

            def stratum(mask, include_bank=True):
                if int(mask.sum()) < 10:
                    return None, int(mask.sum())
                cols = [np.ones(int(mask.sum())), a[mask], male[mask]]
                if include_bank:
                    cols.append(hbcc[mask])
                cols.extend([pmi[mask], pc1[mask], pc2[mask], pc3[mask]])
                return ols_fit(score[mask], cols), int(mask.sum())

            age40_fit, n_age40 = stratum(a >= 40)
            hbcc_fit, n_hbcc = stratum(hbcc == 1, include_bank=False)
            mssm_fit, n_mssm = stratum(hbcc == 0, include_bank=False)

            for donor, value in zip(np.asarray(donor_ids)[keep], score):
                donor_scores[tier].setdefault(str(donor), []).append(float(value))

            regression_rows.append({
                "cell_type": cell_type,
                "n_donors": n_donors,
                "n_zero_donors": n_zero,
                "coverage": round(float(coverage), 4),
                "n_up_genes": len(up),
                "n_down_genes": len(down),
                "slope": unadjusted["slope"],
                "se": unadjusted["se"],
                "p": unadjusted["p"],
                "ci_lo": unadjusted["ci_lo"],
                "ci_hi": unadjusted["ci_hi"],
                "rho": unadjusted["rho"],
                "rho_p": unadjusted["rho_p"],
                "slope_raw_counts": raw_fit["slope"],
                "p_raw_counts": raw_fit["p"],
                "slope_covar_adj": covar["slope"],
                "p_covar_adj": covar["p"],
                "ci_lo_covar_adj": covar["ci_lo"],
                "ci_hi_covar_adj": covar["ci_hi"],
                "n_donors_30plus_cohort": int(mature.sum()),
                "age_min": float(a.min()),
                "age_max": float(a.max()),
                "program_tier": tier,
                "normalization": "log1p_CPM_pseudobulk",
                "donor_filter": "age>=30",
            })
            audit_rows.append({
                "cell_type": cell_type,
                "n_all": n_donors, **values(unadjusted, "all"),
                "n_covar_adj": n_donors, **values(covar, "covar_adj"),
                "n_age_ge40": n_age40, **values(age40_fit, "age_ge40"),
                "n_hbcc": n_hbcc, **values(hbcc_fit, "hbcc"),
                "n_mssm": n_mssm, **values(mssm_fit, "mssm"),
                "program_tier": tier,
                "normalization": "log1p_CPM_pseudobulk",
                "donor_filter": "age>=30",
            })

            # 200-draw matched null, historical sensitivity only. Same seed,
            # same per-tier stream order, same procedure as Item 2's null.
            available = np.ones(len(genes), dtype=bool)
            available[np.concatenate([up, down])] = False
            random_scores = m10_random_scores(
                {"matrix": normalized, "genes": genes}, up, down, N_DRAWS_200, rng, available
            )
            design = psychad_design(metadata, keep)
            random_slopes = np.array([m10_ols(s, design)["slope"] for s in random_scores])
            emp_p = matched_null_effect_pvalue(covar["slope"], random_slopes)
            p2000, q2000 = expected_null[(cell_type, tier)]
            null_rows.append({
                "cohort": "DLPFC_Lifespan_Atlas",
                "unit": cell_type,
                "gene_set": tier,
                "source": "M10",
                "model": "fully_adjusted",
                "status": "TESTED",
                "empirical_pvalue_200": emp_p,
                "n_donors": n_donors,
                "n_up_detected": len(up),
                "n_down_detected": len(down),
                "n_genes_detected": len(up) + len(down),
                "null_iterations": N_DRAWS_200,
                "random_seed": SEED,
                "random_mean_slope": float(random_slopes.mean()),
                "random_sd": float(random_slopes.std(ddof=1)),
                "empirical_pvalue_2000": p2000,
                "empirical_q_2000": q2000,
                "multiplicity_family": "mature_PsychAD_30plus_matched_null",
            })
            print(f"  {tier} {cell_type}: n={n_donors} covar_p={covar['p']:.4g} null200_p={emp_p:.4g}", flush=True)
        stages.append({"stage": f"tier_{tier}", "wall_minutes": (time.time() - t_tier) / 60.0})

    regressions = pd.DataFrame(regression_rows)
    for pcol, qcol in (("p", "q_bh"), ("p_covar_adj", "q_bh_covar_adj"), ("p_raw_counts", "q_bh_raw")):
        regressions[qcol] = np.nan
        for tier in TIERS:
            mask = regressions["program_tier"] == tier
            regressions.loc[mask, qcol] = multipletests(regressions.loc[mask, pcol], method="fdr_bh")[1]
    # Cross-reference: Item 2's composite q over the 170-unit mature family.
    q170 = psy.set_index(["unit", "gene_set"])["composite_q"].to_dict()
    regressions["q_composite_mature170"] = [
        q170[(row["cell_type"], row["program_tier"])] for _, row in regressions.iterrows()
    ]
    regressions.to_csv(REGRESSIONS_CSV, index=False)

    audit = pd.DataFrame(audit_rows)
    for pcol in ("p_all", "p_covar_adj", "p_age_ge40", "p_hbcc", "p_mssm"):
        qcol = pcol.replace("p_", "q_")
        audit[qcol] = np.nan
        for tier in TIERS:
            mask = audit["program_tier"] == tier
            valid = mask & audit[pcol].notna()
            audit.loc[valid, qcol] = multipletests(audit.loc[valid, pcol], method="fdr_bh")[1]
    audit.to_csv(AUDIT_CSV, index=False)

    summary_rows = []
    for tier in TIERS:
        donors = sorted(donor_scores[tier].keys())
        whole = np.array([np.mean(donor_scores[tier][d]) for d in donors])
        idx = np.array([donor_ids.index(d) for d in donors])
        a = age_all[idx]
        male = (metadata["Sex"][idx] == "Male").astype(float)
        hbcc = (metadata["Brain_bank"][idx] == "HBCC").astype(float)
        unadjusted = m10_ols(whole, np.column_stack([np.ones(len(a)), a]))
        covar = m10_ols(whole, np.column_stack([
            np.ones(len(a)), a, male, hbcc,
            metadata["PMI"][idx].astype(float),
            metadata["TechPC1"][idx].astype(float),
            metadata["TechPC2"][idx].astype(float),
            metadata["TechPC3"][idx].astype(float),
        ]))
        summary_rows.append({
            "cohort": "DLPFC_Lifespan_Atlas",
            "n_donors": len(donors),
            "age_min": float(a.min()),
            "age_max": float(a.max()),
            "age_sd": float(a.std(ddof=1)),
            "n_cell_types": len(cell_types),
            "slope": unadjusted["slope"],
            "se": unadjusted["se"],
            "p": unadjusted["p"],
            "ci_lo": unadjusted["ci_lo"],
            "ci_hi": unadjusted["ci_hi"],
            "rho": unadjusted["rho"],
            "rho_p": unadjusted["rho_p"],
            "slope_covar_adj": covar["slope"],
            "p_covar_adj": covar["p"],
            "ci_lo_covar_adj": covar["ci_lo"],
            "ci_hi_covar_adj": covar["ci_hi"],
            "program_tier": tier,
            "normalization": "log1p_CPM_pseudobulk",
            "donor_filter": "age>=30",
        })
    pd.DataFrame(summary_rows).to_csv(SUMMARY_CSV, index=False)

    nulls = pd.DataFrame(null_rows)
    nulls["empirical_q_200"] = np.nan
    for tier in TIERS:
        mask = nulls["gene_set"] == tier
        nulls.loc[mask, "empirical_q_200"] = multipletests(
            nulls.loc[mask, "empirical_pvalue_200"], method="fdr_bh"
        )[1]
    nulls["empirical_pass_200"] = nulls["empirical_q_200"] < 0.05
    nulls["empirical_pass_2000"] = nulls["empirical_q_2000"] < 0.05

    def change(row):
        if row["empirical_pass_200"] and row["empirical_pass_2000"]:
            return "SAME_PASS"
        if not row["empirical_pass_200"] and not row["empirical_pass_2000"]:
            return "SAME_FAIL"
        if row["empirical_pass_2000"]:
            return "PASS_2000_ONLY"
        return "PASS_200_ONLY"

    nulls["verdict_change"] = nulls.apply(change, axis=1)
    nulls.to_csv(NULL200_CSV, index=False)

    # Forest plots of the fully adjusted slopes, mirroring the original M10 figure.
    try:
        import matplotlib.pyplot as plt

        for tier in TIERS:
            sub = regressions[regressions["program_tier"] == tier].sort_values("slope_covar_adj")
            fig, ax = plt.subplots(figsize=(6.5, 7.5))
            ax.errorbar(
                sub["slope_covar_adj"], sub["cell_type"],
                xerr=[sub["slope_covar_adj"] - sub["ci_lo_covar_adj"],
                      sub["ci_hi_covar_adj"] - sub["slope_covar_adj"]],
                fmt="o", color="#0369a1", ecolor="#9ec3dd", ms=4,
            )
            ax.axvline(0, color="#444", lw=0.8, linestyle="--")
            ax.set_xlabel(f"{tier} score slope per year, donors 30-plus (95% CI)")
            ax.set_title("DLPFC lifespan association by cell type, mature only")
            fig.tight_layout()
            fig.savefig(FIGURES / f"followup-mature-psychad-m10-forest-{tier.lower()}.png", dpi=150)
            plt.close(fig)
    except Exception as exc:  # plots are supplementary; never fail the run on them
        print(f"forest plot failed (non-fatal): {exc}", flush=True)

    stages.append({"stage": "total", "wall_minutes": (time.time() - t0) / 60.0})
    runtime = pd.DataFrame([{
        "stage": s["stage"],
        "wall_minutes": round(s["wall_minutes"], 3),
        "model": MODEL,
        "rds_md5": digest,
        "null_seed": SEED,
        "null_iterations_sensitivity": N_DRAWS_200,
        "null_iterations_primary": 2000,
        "started_utc": utc_now(),
    } for s in stages])
    runtime.to_csv(RUNTIME_CSV, index=False)
    print(f"DONE in {(time.time() - t0) / 60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
