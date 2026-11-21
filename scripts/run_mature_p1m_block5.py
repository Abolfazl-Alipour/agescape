#!/usr/bin/env python3
"""Run the core 2000-draw P1m/P0m matched nulls and strict Block 5 comparison."""

from __future__ import annotations

import gc
import os
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_m10_pipeline import ols as m10_ols
from scripts.run_reviewer_phase1_signature_null import random_signature_scores as m10_random_scores
from scripts.run_rev2_transport_null import fit_age_sex, random_signature_scores as mssm_random_scores
from scripts.run_mature_p1m_transport import (
    MIN_AGE,
    MODEL,
    MSSM_H5AD,
    RDS,
    RDS_MD5,
    TIERS,
    TRANSPORT_CSV,
    donor_metadata,
    load_programs,
    md5sum,
    psychad_design,
    read_csr_block,
    read_obs_columns,
)
from src.followup_muscle import matched_null_effect_pvalue

RESULTS = REPO / "results"
NULL_CSV = RESULTS / "followup-mature-p1m-transport-null.csv"
COMPARISON_CSV = RESULTS / "followup-mature-p1m-block5-p0m-comparison.csv"
SUMMARY_CSV = RESULTS / "followup-mature-p1m-block5-summary.csv"
RUNTIME_CSV = RESULTS / "followup-mature-p1m-block5-runtime.csv"
ITERATIONS = 2000
SEED = 20261105


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_mssm_pseudobulk(cell_types: list[str], scores: pd.DataFrame):
    pairs = []
    for cell_type in cell_types:
        sub = scores[scores["cell_type"] == cell_type]
        donors = sorted(sub["donor_id"].astype(str))
        pairs.extend((cell_type, donor) for donor in donors)
    cache = REPO / "scratch" / "mature-p1m-block5" / "mssm_pseudobulk.npz"
    if cache.exists():
        with np.load(cache, allow_pickle=False) as cached:
            cached_pairs = [tuple(row) for row in cached["pairs"].astype(str)]
            assert cached_pairs == pairs
            print(f"MSSM: loading {len(pairs)} cached mature pseudobulks", flush=True)
            return cached["genes"].astype(str), pairs, cached["pseudobulk"]
    print(f"MSSM: building {len(pairs)} mature donor-cell-type pseudobulks", flush=True)
    pair_counts = np.zeros(len(pairs), dtype=np.int64)
    with h5py.File(MSSM_H5AD, "r") as handle:
        columns = read_obs_columns(handle)
        donor_meta = donor_metadata(columns)
        donor_arr = columns["donor_id"].astype(str)
        ct_arr = columns["cell_type"].astype(str)
        primary = columns["is_primary_data"].astype(bool)
        donor_codes, donor_levels = pd.factorize(donor_arr, sort=False)
        ct_codes, ct_levels = pd.factorize(ct_arr, sort=False)
        donor_lookup = {donor: i for i, donor in enumerate(donor_levels.astype(str))}
        ct_lookup = {cell_type: i for i, cell_type in enumerate(ct_levels.astype(str))}
        lookup = np.full((len(donor_levels), len(ct_levels)), -1, dtype=np.int64)
        for i, (cell_type, donor) in enumerate(pairs):
            assert donor in donor_meta.index and donor_meta.loc[donor, "age"] >= MIN_AGE
            lookup[donor_lookup[donor], ct_lookup[cell_type]] = i
        pair_ids = lookup[donor_codes, ct_codes]
        pair_ids[~primary] = -1
        genes = np.asarray(__import__("anndata").io.read_elem(handle["var"]["feature_name"])).astype(str)
        pair_sums = np.zeros((len(pairs), len(genes)), dtype=np.float64)
        chunk = 5000
        for start in range(0, len(donor_arr), chunk):
            stop = min(start + chunk, len(donor_arr))
            local_pairs = pair_ids[start:stop]
            keep = local_pairs >= 0
            if not keep.any():
                continue
            matrix = read_csr_block(handle["X"], start, stop, len(genes)).astype(np.float64)
            library_totals = np.asarray(matrix.sum(axis=1)).ravel()
            scale = 1e6 / np.where(library_totals > 0, library_totals, 1.0)
            matrix.data *= np.repeat(scale, np.diff(matrix.indptr))
            np.log1p(matrix.data, out=matrix.data)
            matrix = matrix[keep]
            local_pairs = local_pairs[keep]
            unique, inverse, cell_counts = np.unique(local_pairs, return_inverse=True, return_counts=True)
            indicator = sp.csr_matrix(
                (np.ones(len(local_pairs)), (inverse, np.arange(len(local_pairs)))),
                shape=(len(unique), len(local_pairs)),
            )
            pair_sums[unique, :] += (indicator @ matrix).toarray()
            np.add.at(pair_counts, unique, cell_counts)
            del matrix, indicator
            if start % 200000 == 0:
                print(f"  MSSM cells {start} to {stop} of {len(donor_arr)}", flush=True)
    expected = scores.set_index(["cell_type", "donor_id"])["n_cells"]
    for i, pair in enumerate(pairs):
        assert pair_counts[i] == int(expected.loc[pair]), f"MSSM {pair}: cell count mismatch"
    pseudobulk = pair_sums / pair_counts[:, None]
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, genes=genes, pairs=np.asarray(pairs, dtype=str), pseudobulk=pseudobulk)
    return genes, pairs, pseudobulk


def run_mssm_nulls(programs: dict[str, dict[str, list[str]]], transport: pd.DataFrame) -> pd.DataFrame:
    scores = pd.read_csv(RESULTS / "followup-mature-p1m-m9-donor-scores.csv")
    scores = scores[scores["cohort"] == "MSSM_DLPFC"].copy()
    scores["donor_id"] = scores["donor_id"].astype(str)
    cell_types = sorted(transport[
        (transport["cohort"] == "MSSM_DLPFC") & (transport["gene_set"] == "P1m")
    ]["unit"])
    genes, pairs, pseudobulk = build_mssm_pseudobulk(cell_types, scores)
    lookup = {gene: i for i, gene in enumerate(genes)}
    rows = []
    for tier in TIERS:
        program = programs[tier]
        rng = np.random.default_rng(SEED)
        tier_rows = []
        for cell_type in cell_types:
            indices = np.array([i for i, pair in enumerate(pairs) if pair[0] == cell_type], dtype=int)
            donors = [pairs[i][1] for i in indices]
            matrix = pseudobulk[indices, :].T
            meta = scores[scores["cell_type"] == cell_type].set_index("donor_id").loc[donors]
            up = np.array([lookup[g] for g in program["up"] if g in lookup], dtype=int)
            down = np.array([lookup[g] for g in program["down"] if g in lookup], dtype=int)
            actual_scores = matrix[up, :].mean(axis=0) - matrix[down, :].mean(axis=0)
            max_diff = float(np.max(np.abs(actual_scores - meta[f"score_{tier}"].to_numpy(float))))
            assert max_diff < 1e-7, f"MSSM {tier} {cell_type}: score mismatch {max_diff}"
            fit = fit_age_sex(actual_scores, meta["age"].to_numpy(float), meta["sex"].astype(str).to_numpy())
            expected = transport[
                (transport["cohort"] == "MSSM_DLPFC") & (transport["unit"] == cell_type)
                & (transport["gene_set"] == tier)
            ]["composite_slope"].iloc[0]
            assert np.isclose(fit["slope"], expected, rtol=1e-7, atol=1e-9)
            available = np.ones(len(genes), dtype=bool)
            available[np.concatenate([up, down])] = False
            random_scores = mssm_random_scores(
                {"matrix": matrix, "genes": genes}, up, down, ITERATIONS, rng, available
            )
            random_slopes = np.array([
                fit_age_sex(score, meta["age"].to_numpy(float), meta["sex"].astype(str).to_numpy())["slope"]
                for score in random_scores
            ])
            tier_rows.append({
                "cohort": "MSSM_DLPFC",
                "unit": cell_type,
                "gene_set": tier,
                "source": "M9",
                "model": "age_plus_sex",
                "status": "TESTED",
                "empirical_pvalue": matched_null_effect_pvalue(fit["slope"], random_slopes),
                "n_donors": len(donors),
                "n_up_detected": len(up),
                "n_down_detected": len(down),
                "n_genes_detected": len(up) + len(down),
                "null_iterations": ITERATIONS,
                "random_seed": SEED,
                "random_mean_slope": float(random_slopes.mean()),
                "random_sd": float(random_slopes.std(ddof=1)),
                "max_abs_score_diff_vs_transport": max_diff,
                "multiplicity_family": "mature_MSSM_30plus_matched_null",
            })
            print(f"  MSSM {tier} {cell_type}: p={tier_rows[-1]['empirical_pvalue']:.12g}", flush=True)
        q = multipletests([r["empirical_pvalue"] for r in tier_rows], method="fdr_bh")[1]
        for row, qvalue in zip(tier_rows, q):
            row["empirical_q"] = qvalue
            row["empirical_pass"] = bool(qvalue < 0.05)
        rows.extend(tier_rows)
    del pseudobulk
    gc.collect()
    return pd.DataFrame(rows)


def run_m10_nulls(programs: dict[str, dict[str, list[str]]], transport: pd.DataFrame) -> pd.DataFrame:
    digest = md5sum(RDS)
    assert digest == RDS_MD5, f"PsychAD MD5 mismatch: {digest}"
    print(f"PsychAD null: MD5 verified {digest}", flush=True)
    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    metadata = {
        name: np.asarray(cd.listData[name]).ravel()
        for name in ("Age", "Sex", "Brain_bank", "PMI", "TechPC1", "TechPC2", "TechPC3")
    }
    mature = metadata["Age"].astype(float) >= MIN_AGE
    cell_types = sorted(transport[
        (transport["cohort"] == "DLPFC_Lifespan_Atlas") & (transport["gene_set"] == "P1m")
    ]["unit"])
    rows = []
    for tier in TIERS:
        program = programs[tier]
        rng = np.random.default_rng(SEED)
        tier_rows = []
        for cell_type in cell_types:
            source = obj.assays.data.listData[cell_type]
            genes = np.array([str(gene) for gene in source.Dimnames[0]])
            matrix = sp.csc_matrix(
                (np.asarray(source.x), np.asarray(source.i), np.asarray(source.p)), shape=tuple(source.Dim)
            )
            library = np.asarray(matrix.sum(axis=0)).ravel()
            keep = mature & (library > 0)
            matrix = matrix[:, keep]
            normalized = matrix @ sp.diags(1e6 / library[keep])
            normalized = normalized.tocsc()
            normalized.data = np.log1p(normalized.data)
            lookup = {gene: i for i, gene in enumerate(genes)}
            up = np.array([lookup[g] for g in program["up"] if g in lookup], dtype=int)
            down = np.array([lookup[g] for g in program["down"] if g in lookup], dtype=int)
            actual_scores = np.asarray(normalized[up, :].mean(axis=0)).ravel() - np.asarray(normalized[down, :].mean(axis=0)).ravel()
            design = psychad_design(metadata, keep)
            fit = m10_ols(actual_scores, design)
            expected = transport[
                (transport["cohort"] == "DLPFC_Lifespan_Atlas") & (transport["unit"] == cell_type)
                & (transport["gene_set"] == tier)
            ]["composite_slope"].iloc[0]
            assert np.isclose(fit["slope"], expected, rtol=1e-10, atol=1e-12)
            available = np.ones(len(genes), dtype=bool)
            available[np.concatenate([up, down])] = False
            random_scores = m10_random_scores(
                {"matrix": normalized, "genes": genes}, up, down, ITERATIONS, rng, available
            )
            random_slopes = np.array([m10_ols(score, design)["slope"] for score in random_scores])
            tier_rows.append({
                "cohort": "DLPFC_Lifespan_Atlas",
                "unit": cell_type,
                "gene_set": tier,
                "source": "M10",
                "model": "fully_adjusted",
                "status": "TESTED",
                "empirical_pvalue": matched_null_effect_pvalue(fit["slope"], random_slopes),
                "n_donors": int(keep.sum()),
                "n_up_detected": len(up),
                "n_down_detected": len(down),
                "n_genes_detected": len(up) + len(down),
                "null_iterations": ITERATIONS,
                "random_seed": SEED,
                "random_mean_slope": float(random_slopes.mean()),
                "random_sd": float(random_slopes.std(ddof=1)),
                "max_abs_score_diff_vs_transport": 0.0,
                "multiplicity_family": "mature_PsychAD_30plus_matched_null",
            })
            print(f"  PsychAD {tier} {cell_type}: p={tier_rows[-1]['empirical_pvalue']:.12g}", flush=True)
        q = multipletests([r["empirical_pvalue"] for r in tier_rows], method="fdr_bh")[1]
        for row, qvalue in zip(tier_rows, q):
            row["empirical_q"] = qvalue
            row["empirical_pass"] = bool(qvalue < 0.05)
        rows.extend(tier_rows)
    return pd.DataFrame(rows)


def unit_status(composite_pass: bool, empirical_pass: bool) -> str:
    if composite_pass and empirical_pass:
        return "BOTH_PASS"
    if empirical_pass:
        return "EMPIRICAL_ONLY"
    if composite_pass:
        return "COMPOSITE_ONLY"
    return "FAIL_BOTH"


def build_comparison(transport: pd.DataFrame, nulls: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    merged = transport.merge(
        nulls[["cohort", "unit", "gene_set", "empirical_pvalue", "empirical_q", "empirical_pass"]],
        on=["cohort", "unit", "gene_set"], how="inner", validate="one_to_one"
    )
    merged["both_pass"] = merged["composite_pass"] & merged["empirical_pass"]
    merged["status"] = [unit_status(c, e) for c, e in zip(merged["composite_pass"], merged["empirical_pass"])]
    key = ["cohort", "unit"]
    p0 = merged[merged["gene_set"] == "P0m"].set_index(key)
    p1 = merged[merged["gene_set"] == "P1m"].set_index(key)
    assert p0.index.equals(p1.index)
    rows = []
    for index in p0.index:
        a = p0.loc[index]
        b = p1.loc[index]
        rescue = bool(a["both_pass"] and not b["composite_pass"] and not b["empirical_pass"])
        rows.append({
            "cohort": index[0],
            "unit": index[1],
            "n_donors": int(b["n_donors"]),
            "p0m_composite_pvalue": a["composite_pvalue"],
            "p0m_composite_q": a["composite_q"],
            "p0m_empirical_pvalue": a["empirical_pvalue"],
            "p0m_empirical_q": a["empirical_q"],
            "p0m_status": a["status"],
            "p0m_both_pass": a["both_pass"],
            "p1m_composite_pvalue": b["composite_pvalue"],
            "p1m_composite_q": b["composite_q"],
            "p1m_empirical_pvalue": b["empirical_pvalue"],
            "p1m_empirical_q": b["empirical_q"],
            "p1m_status": b["status"],
            "p1m_both_pass": b["both_pass"],
            "p1m_fails_both_strict": bool(not b["composite_pass"] and not b["empirical_pass"]),
            "strict_p0m_rescue": rescue,
        })
    comparison = pd.DataFrame(rows).sort_values(key).reset_index(drop=True)
    n_rescue = int(comparison["strict_p0m_rescue"].sum())
    verdict = "PROVES" if n_rescue else "DOES NOT PROVE"
    summary_rows = []
    for tier in ("P0m", "P1m"):
        sub = merged[merged["gene_set"] == tier]
        summary_rows.append({
            "gene_set": tier,
            "n_units": len(sub),
            "both_pass": int(sub["both_pass"].sum()),
            "empirical_only": int((sub["status"] == "EMPIRICAL_ONLY").sum()),
            "composite_only": int((sub["status"] == "COMPOSITE_ONLY").sum()),
            "fail_both": int((sub["status"] == "FAIL_BOTH").sum()),
            "strict_p0m_rescue_units": n_rescue,
            "strict_rescue_verdict": verdict,
        })
    return comparison, pd.DataFrame(summary_rows), verdict


def main() -> int:
    start = time.monotonic()
    start_utc = utc_now()
    print(f"Mature Block 5 started {start_utc}; model {MODEL}", flush=True)
    programs = load_programs()
    transport = pd.read_csv(TRANSPORT_CSV)
    nulls = pd.concat([
        run_mssm_nulls(programs, transport),
        run_m10_nulls(programs, transport),
    ], ignore_index=True)
    assert (nulls["null_iterations"] == ITERATIONS).all()
    assert (nulls["random_seed"] == SEED).all()
    assert nulls["n_donors"].min() > 0
    nulls.to_csv(NULL_CSV, index=False)
    comparison, summary, verdict = build_comparison(transport, nulls)
    comparison.to_csv(COMPARISON_CSV, index=False)
    summary.to_csv(SUMMARY_CSV, index=False)
    end_utc = utc_now()
    duration = time.monotonic() - start
    pd.DataFrame([{
        "model": MODEL,
        "start_utc": start_utc,
        "end_utc": end_utc,
        "duration_seconds": duration,
        "null_iterations": ITERATIONS,
        "random_seed": SEED,
        "n_null_rows": len(nulls),
        "n_comparison_units": len(comparison),
        "strict_rescue_verdict": verdict,
    }]).to_csv(RUNTIME_CSV, index=False)
    print(summary.to_string(index=False), flush=True)
    print(f"Strict P0m rescue verdict: {verdict}", flush=True)
    print(f"Mature Block 5 ended {end_utc}; duration {duration:.3f} seconds", flush=True)
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    raise SystemExit(main())
