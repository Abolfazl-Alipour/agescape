#!/usr/bin/env python3
"""Item 7 mature MSSM and HBCC interneuron replication."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import run_followup_11pt_block9_mature_interneuron as block9
from scripts import run_rev2_de_interneuron as original
from scripts.run_rev3_composition_fly_worm_zebrafish import fit_f_dist_smyth
from src.followup_organs import write_clean_csv

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
MSSM_PATH = block9.MSSM_PATH
HBCC_PATH = Path(
    "/home/reza/projects/agescape-revision-data/human-h5ads/"
    "5c97eeeb-7e52-44b3-b010-b832b1f5424c.h5ad"
)
DEFAULT_SCRATCH = block9.DEFAULT_SCRATCH
MATURE_AGE = block9.MATURE_AGE
SHARED_GENES = block9.SHARED_GENES
SEED = original.SEED
INTERNEURON_TYPES = original.INTERNEURON_TYPES
SHORT = original.SHORT
SUBTYPES = INTERNEURON_TYPES[:3]
N_PERMS = 1000
EXPECTED_DONORS = {"VIP": 193, "pvalb": 185, "sst": 174, "GABA": 153}
EXPECTED_TESTED = {"VIP": 1932, "pvalb": 1873, "sst": 1906, "GABA": 1883}
EXPECTED_HBCC_CELLS = {"VIP": 74520, "pvalb": 67895, "sst": 71641, "GABA": 31968}
GENE_COLUMNS = [
    "cohort", "cell_type", "short_type", "model", "primary", "status", "gene",
    "slope", "se", "ci95_low", "ci95_high", "pvalue", "q_bh", "n_donors",
    "n_genes_tested", "mean_cp10k", "detected_frac", "age_min", "age_median",
    "age_max", "design",
]
GENE_OUTPUT = RESULTS / "followup-mature-interneuron.csv"
REPLICATION_OUTPUT = RESULTS / "followup-mature-interneuron-replication.csv"
GENES_OUTPUT = RESULTS / "followup-mature-interneuron-genes.csv"
MODERATED_GENES_OUTPUT = RESULTS / "followup-mature-interneuron-genes-moderated.csv"
MODERATED_DE_OUTPUT = RESULTS / "followup-mature-interneuron-de-moderated.csv"
DISCOVERY_OUTPUT = RESULTS / "followup-mature-interneuron-discovery-ols.csv"
POWER_OUTPUT = RESULTS / "followup-mature-interneuron-power.csv"
SIX_GENE_OUTPUT = RESULTS / "followup-mature-interneuron-six-gene-status.csv"
GATE_OUTPUT = RESULTS / "followup-mature-interneuron-moderation-gate.csv"
BH_OUTPUT = RESULTS / "followup-mature-interneuron-bh-audit.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def wait_for_memory(minimum: float = 20.0, timeout_seconds: int = 3600) -> float:
    start = time.monotonic()
    while True:
        available = block9.available_memory_gib()
        if available >= minimum:
            return available
        if time.monotonic() - start >= timeout_seconds:
            raise MemoryError(
                f"Only {available:.1f} GiB available after waiting 60 minutes; "
                f"at least {minimum:.0f} GiB is required."
            )
        original.log(f"Waiting 60 seconds for at least {minimum:.0f} GiB available memory")
        time.sleep(60)


block9.require_memory = wait_for_memory


def validate_hbcc() -> dict:
    import anndata as ad

    wait_for_memory()
    if not HBCC_PATH.exists():
        raise FileNotFoundError(f"Required HBCC UUID object is absent: {HBCC_PATH}")
    backed = ad.read_h5ad(HBCC_PATH, backed="r")
    try:
        assert backed.shape == (1486324, 34176), backed.shape
        obs = backed.obs
        required = {"cell_type", "disease", "is_primary_data", "n_counts", "development_stage", "donor_id"}
        assert required.issubset(obs.columns), sorted(required - set(obs.columns))
        cell_counts = {}
        for cell_type in INTERNEURON_TYPES:
            count = int(original._cat_eq(obs["cell_type"], cell_type).sum())
            cell_counts[SHORT[cell_type]] = count
            assert count == EXPECTED_HBCC_CELLS[SHORT[cell_type]], (cell_type, count)
        normal_count = int(original._cat_eq(obs["disease"], "normal").sum())
        assert normal_count > 0
        primary = original._cat_eq(obs["is_primary_data"], "true")
        assert bool(primary.all())
        age = original._age_from_stage(obs["development_stage"])
        assert float(age.min()) == 2.0 and float(age.max()) == 85.0
        donor_count = int(obs["donor_id"].nunique())
        assert donor_count == 300
        validation = {
            "path": str(HBCC_PATH),
            "size_bytes": HBCC_PATH.stat().st_size,
            "shape": list(backed.shape),
            "obs_required_columns": sorted(required),
            "cell_type_counts": cell_counts,
            "disease_normal_count": normal_count,
            "is_primary_data_all_true": True,
            "n_counts_present": True,
            "raw_sum_present": "raw_sum" in obs.columns,
            "development_stage_year_min": float(age.min()),
            "development_stage_year_max": float(age.max()),
            "month_stages_parse_as_missing_age": bool(age.isna().any()),
            "donor_count": donor_count,
            "status": "VALID",
        }
    finally:
        backed.file.close()
    return validation


def _assert_close(left, right, label: str) -> None:
    if not np.allclose(np.asarray(left, float), np.asarray(right, float), atol=1e-9, rtol=0):
        difference = float(np.nanmax(np.abs(np.asarray(left, float) - np.asarray(right, float))))
        raise AssertionError(f"{label} differs from cache; maximum absolute difference {difference}")


def validate_mssm_unit(cell_type: str, scratch: Path, block_csv: pd.DataFrame) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    short = SHORT[cell_type]
    counts, meta, genes = original.load_counts_npz(scratch, "mssm_mature30", cell_type)
    assert len(meta) == EXPECTED_DONORS[short]
    assert meta["age"].ge(MATURE_AGE).all()
    assert not meta["donor"].duplicated().any()
    de = original.load_de_npz(scratch, "mssm_mature30", cell_type)
    assert len(de) == EXPECTED_TESTED[short]
    unit_csv = block_csv[block_csv["cell_type"].eq(cell_type)].copy()
    wls_csv = unit_csv[unit_csv["model"].eq("rev2_voom_wls")]
    ols_csv = unit_csv[unit_csv["model"].eq("age_plus_sex_depth_ols")]
    assert list(wls_csv["gene"]) == list(de["gene"])
    assert list(ols_csv["gene"]) == list(de["gene"])
    _assert_close(wls_csv["slope"], de["slope"], f"{short} MSSM WLS slope")
    _assert_close(wls_csv["pvalue"], de["pvalue"], f"{short} MSSM WLS pvalue")
    _assert_close(ols_csv["slope"], de["slope_ols"], f"{short} MSSM OLS slope")
    _assert_close(ols_csv["pvalue"], de["pvalue_ols"], f"{short} MSSM OLS pvalue")
    summary = {
        "cohort": "MSSM_DLPFC",
        "cell_type": cell_type,
        "short_type": short,
        "status": "FIT",
        "cache_validation": "PASS",
        "n_donors": int(len(meta)),
        "n_genes_tested": int(len(de)),
        "n_fdr_wls": int(de["q_bh"].lt(0.05).sum()),
        "n_fdr_ols_same_genes": int(de["q_bh_ols"].lt(0.05).sum()),
        "n_fdr_depth_unadjusted_ols": int(unit_csv.loc[
            unit_csv["model"].eq("age_plus_sex_ols_depth_unadjusted"), "q_bh"
        ].lt(0.05).sum()),
        "age_min": float(meta["age"].min()),
        "age_median": float(meta["age"].median()),
        "age_max": float(meta["age"].max()),
        "age_sd": float(meta["age"].std()),
        "design": str(wls_csv["design"].iloc[0]),
        "median_n_cells": float(meta["n_cells"].median()),
    }
    assert summary["n_fdr_wls"] == {"VIP": 148, "pvalb": 133, "sst": 366, "GABA": 39}[short]
    assert summary["n_fdr_ols_same_genes"] == {"VIP": 145, "pvalb": 147, "sst": 204, "GABA": 32}[short]
    assert summary["n_fdr_depth_unadjusted_ols"] == {"VIP": 385, "pvalb": 353, "sst": 359, "GABA": 153}[short]
    return unit_csv.reindex(columns=GENE_COLUMNS), summary, de


def refit_mssm_unit(cell_type: str, scratch: Path) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    wait_for_memory()
    counts, meta, genes, _ = block9.fetch_mature_unit_h5ad(
        MSSM_PATH, cell_type, scratch, "mssm_mature30"
    )
    if counts is None or len(meta) < original.MIN_DONORS:
        raise RuntimeError(f"MSSM {cell_type} has fewer than 10 mature donors after refit")
    summary, merged = original.run_unit_de(
        counts, meta, genes, "MSSM_DLPFC", cell_type, scratch, "mssm_mature30"
    )
    sensitivity = original.ols_arm(counts, meta, genes, columns=("age", "sex_male"))
    frame = block9.model_frames(merged, sensitivity)
    frame["cohort"] = "MSSM_DLPFC"
    frame["cell_type"] = cell_type
    frame["short_type"] = SHORT[cell_type]
    frame["status"] = "FIT"
    frame["age_min"] = float(meta["age"].min())
    frame["age_median"] = float(meta["age"].median())
    frame["age_max"] = float(meta["age"].max())
    frame["n_genes_tested"] = int(summary["n_genes_tested"])
    frame["design"] = summary["design"]
    summary.update({
        "short_type": SHORT[cell_type], "status": "FIT", "cache_validation": "REFIT_AFTER_FAILURE",
        "age_median": float(meta["age"].median()),
        "n_fdr_depth_unadjusted_ols": int(sensitivity["q_bh_ols"].lt(0.05).sum()),
    })
    return frame.reindex(columns=GENE_COLUMNS), summary, original.load_de_npz(scratch, "mssm_mature30", cell_type)


def fit_hbcc_unit(cell_type: str, scratch: Path) -> tuple[pd.DataFrame, dict, pd.DataFrame] | None:
    available = wait_for_memory()
    counts, meta, genes, _ = block9.fetch_mature_unit_h5ad(
        HBCC_PATH, cell_type, scratch, "hbcc_mature30"
    )
    if counts is None or len(meta) < original.MIN_DONORS:
        return None
    assert meta["age"].ge(MATURE_AGE).all() and not meta["donor"].duplicated().any()
    summary, merged = original.run_unit_de(
        counts, meta, genes, "HBCC_DLPFC", cell_type, scratch, "hbcc_mature30"
    )
    sensitivity = original.ols_arm(counts, meta, genes, columns=("age", "sex_male"))
    frame = block9.model_frames(merged, sensitivity)
    frame["cohort"] = "HBCC_DLPFC"
    frame["cell_type"] = cell_type
    frame["short_type"] = SHORT[cell_type]
    frame["status"] = "FIT"
    frame["age_min"] = float(meta["age"].min())
    frame["age_median"] = float(meta["age"].median())
    frame["age_max"] = float(meta["age"].max())
    frame["n_genes_tested"] = int(summary["n_genes_tested"])
    frame["design"] = summary["design"]
    summary.update({
        "short_type": SHORT[cell_type], "status": "FIT", "status_reason": "none",
        "age_median": float(meta["age"].median()),
        "n_fdr_depth_unadjusted_ols": int(sensitivity["q_bh_ols"].lt(0.05).sum()),
        "min_available_gib_before_unit": available,
    })
    original.log(
        f"HBCC mature {SHORT[cell_type]}: {len(meta)} donors, WLS FDR {summary['n_fdr_wls']}, "
        f"OLS FDR {summary['n_fdr_ols_same_genes']}"
    )
    de = original.load_de_npz(scratch, "hbcc_mature30", cell_type)
    del counts, meta, genes, merged, sensitivity
    gc.collect()
    return frame.reindex(columns=GENE_COLUMNS), summary, de


def voom_wls_counts_with_s2(counts, metadata, gene_names, columns=("age", "sex_male", "mean_log_umi")):
    keep, gene_stats = original.filter_genes(counts, gene_names)
    C = counts[:, keep].astype(float)
    genes = gene_names[keep]
    lib = counts.sum(axis=1).astype(float)
    Y = np.log1p(C / lib[:, None] * 1e4)
    Y = Y.T
    n_cells = metadata["n_cells"].to_numpy(float)
    cell_factor = n_cells / np.median(n_cells)
    X, colnames = original.design_matrix(metadata, list(columns))
    age_idx = colnames.index("age")
    D, k = X.shape

    def ols_pass(Ym):
        beta, *_ = np.linalg.lstsq(X, Ym.T, rcond=None)
        resid = Ym.T - X @ beta
        dof = D - k
        s2 = np.sum(resid ** 2, axis=0) / max(dof, 1)
        return s2

    mean_logcpm = Y.mean(axis=1)
    s2_ols = ols_pass(Y)
    v_g, _, _ = original._fit_trend(mean_logcpm, s2_ols)
    W = cell_factor[None, :] / v_g[:, None]
    med = np.median(W, axis=1, keepdims=True)
    W = np.clip(W, 0.25 * med, 4.0 * med)
    beta, se, s2_w, dof = original._wls_batched(Y, X, W)
    resid_w = Y - beta @ X.T
    s2_wls = np.sum(resid_w ** 2, axis=1) / max(dof, 1)
    v_g2, _, _ = original._fit_trend(mean_logcpm, s2_wls)
    W2 = cell_factor[None, :] / v_g2[:, None]
    med2 = np.median(W2, axis=1, keepdims=True)
    W2 = np.clip(W2, 0.25 * med2, 4.0 * med2)
    beta, se, s2_w, dof = original._wls_batched(Y, X, W2)

    slope = beta[:, age_idx]
    se_age = se[:, age_idx]
    tstat = np.divide(slope, se_age, out=np.zeros_like(slope), where=se_age > 0)
    pvalue = 2 * stats.t.sf(np.abs(tstat), dof)
    qvalue = multipletests(pvalue, method="fdr_bh")[1]
    tcrit = stats.t.ppf(0.975, dof)
    out = pd.DataFrame({
        "gene": genes,
        "n_donors": D,
        "slope": slope,
        "se": se_age,
        "t": tstat,
        "pvalue": pvalue,
        "q_bh": qvalue,
        "ci95_low": slope - tcrit * se_age,
        "ci95_high": slope + tcrit * se_age,
        "mean_cp10k": gene_stats.loc[gene_stats["kept"], "mean_cp10k"].to_numpy(),
        "detected_frac": gene_stats.loc[gene_stats["kept"], "detected_frac"].to_numpy(),
        "median_n_cells": float(np.median(n_cells)),
        "mean_n_cells": float(np.mean(n_cells)),
    })
    out.attrs["n_genes_total"] = int(len(gene_names))
    out.attrs["n_genes_tested"] = int(keep.sum())
    out.attrs["n_genes_technical_dropped"] = int(original.technical_gene_mask(gene_names).sum())
    out.attrs["design"] = " + ".join(colnames)
    out.attrs["dof"] = int(dof)
    out.attrs["gene_stats"] = gene_stats
    return out, s2_w, dof


def recover_moderation_gate() -> tuple[pd.DataFrame, bool]:
    unmoderated = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv")
    moderated = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes-moderated.csv")
    priors = pd.read_csv(RESULTS / "reviewer-rev2-human-de-moderated.csv")
    rows = []
    for _, prior in priors.iterrows():
        prefix = "mssm" if prior["cohort"] == "MSSM_DLPFC" else "hbcc"
        left = unmoderated.loc[unmoderated["cell_type"].eq(prior["cell_type"]), ["gene", f"se_{prefix}"]]
        right = moderated.loc[moderated["cell_type"].eq(prior["cell_type"]), ["gene", f"se_{prefix}_mod"]]
        joined = left.merge(right, on="gene", validate="one_to_one")
        ratio = (joined[f"se_{prefix}_mod"] / joined[f"se_{prefix}"]).to_numpy(float) ** 2
        denominator = ratio * (float(prior["resid_dof"]) + float(prior["d0_prior"])) - float(prior["resid_dof"])
        s2 = float(prior["d0_prior"]) * float(prior["s0_sq_prior"]) / denominator
        assert np.isfinite(s2).all() and (s2 > 0).all()
        recovered_s0, recovered_d0 = fit_f_dist_smyth(s2, float(prior["resid_dof"]))
        s0_relative = abs(recovered_s0 / float(prior["s0_sq_prior"]) - 1)
        d0_relative = abs(recovered_d0 / float(prior["d0_prior"]) - 1)
        rows.append({
            "cohort": prior["cohort"], "cell_type": prior["cell_type"],
            "short_type": SHORT[prior["cell_type"]], "n_recovered_genes": len(s2),
            "frozen_s0_sq": float(prior["s0_sq_prior"]), "recovered_s0_sq": recovered_s0,
            "s0_sq_relative_difference": s0_relative,
            "frozen_d0": float(prior["d0_prior"]), "recovered_d0": recovered_d0,
            "d0_relative_difference": d0_relative,
            "pass_5_percent": bool(s0_relative <= 0.05 and d0_relative <= 0.05),
            "recovery_universe": "frozen selected-gene rows available for this unit",
        })
    gate = pd.DataFrame(rows)
    return gate, bool(gate["pass_5_percent"].all())


def fit_s2_and_validate(cell_type: str, tag: str, de: pd.DataFrame, scratch: Path) -> tuple[np.ndarray, int, dict]:
    wait_for_memory()
    counts, meta, genes = original.load_counts_npz(scratch, tag, cell_type)
    fitted, s2, dof = voom_wls_counts_with_s2(counts, meta, genes)
    assert list(fitted["gene"]) == list(de["gene"])
    for column in ("slope", "se", "pvalue", "q_bh"):
        _assert_close(fitted[column], de[column], f"{tag} {SHORT[cell_type]} {column}")
    validation = {
        "cell_type": cell_type, "tag": tag, "n_genes": len(de), "resid_dof": int(dof),
        "max_abs_slope_difference": float(np.max(np.abs(fitted["slope"] - de["slope"]))),
        "max_abs_se_difference": float(np.max(np.abs(fitted["se"] - de["se"]))),
        "max_abs_pvalue_difference": float(np.max(np.abs(fitted["pvalue"] - de["pvalue"]))),
        "max_abs_q_difference": float(np.max(np.abs(fitted["q_bh"] - de["q_bh"]))),
    }
    del counts, meta, genes, fitted
    gc.collect()
    return s2, dof, validation


def moderate_fit(de: pd.DataFrame, s2: np.ndarray, dof: int) -> tuple[pd.DataFrame, dict]:
    s0_sq, d0 = fit_f_dist_smyth(s2, dof)
    posterior = (dof * s2 + d0 * s0_sq) / (dof + d0)
    se_mod = de["se"].to_numpy(float) * np.sqrt(posterior / s2)
    df_post = dof + d0
    p_mod = 2 * stats.t.sf(np.abs(de["slope"].to_numpy(float) / se_mod), df_post)
    q_mod = multipletests(p_mod, method="fdr_bh")[1]
    fitted = de.copy()
    fitted["se_mod"] = se_mod
    fitted["pvalue_mod"] = p_mod
    fitted["q_bh_mod"] = q_mod
    fitted["s2_w"] = s2
    fitted["resid_dof"] = dof
    fitted["d0_prior"] = d0
    fitted["s0_sq_prior"] = s0_sq
    fitted["df_post"] = df_post
    prior = {
        "resid_dof": int(dof), "d0_prior": float(d0), "s0_sq_prior": float(s0_sq),
        "df_post": float(df_post), "n_de_unmod": int(de["q_bh"].lt(0.05).sum()),
        "n_de_mod": int((q_mod < 0.05).sum()),
        "n_de_unmod_fc025": int(((de["q_bh"] < 0.05) & (de["slope"].abs() / np.log(2) >= 0.25)).sum()),
        "n_de_mod_fc025": int(((q_mod < 0.05) & (de["slope"].abs() / np.log(2) >= 0.25)).sum()),
    }
    return fitted, prior


def build_replication(scratch: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries, gene_rows = [], []
    for cell_type in INTERNEURON_TYPES:
        m = original.load_de_npz(scratch, "mssm_mature30", cell_type)
        h = original.load_de_npz(scratch, "hbcc_mature30", cell_type)
        merged = m.merge(h, on="gene", suffixes=("_mssm", "_hbcc"))
        merged["same_direction"] = merged["slope_mssm"] * merged["slope_hbcc"] > 0
        merged["mssm_fdr"] = merged["q_bh_mssm"] < 0.05
        merged["hbcc_fdr"] = merged["q_bh_hbcc"] < 0.05
        merged["both_fdr"] = merged["mssm_fdr"] & merged["hbcc_fdr"] & merged["same_direction"]
        cm, mm, gm = original.load_counts_npz(scratch, "mssm_mature30", cell_type)
        ch, mh, gh = original.load_counts_npz(scratch, "hbcc_mature30", cell_type)
        gw = merged
        pearson = float(stats.pearsonr(gw["slope_mssm"], gw["slope_hbcc"]).statistic)
        spearman = float(stats.spearmanr(gw["slope_mssm"], gw["slope_hbcc"]).statistic)
        same_dir = int(gw["same_direction"].sum())
        binom_p = float(stats.binomtest(same_dir, len(gw), 0.5, alternative="greater").pvalue)
        sel = gw[gw["mssm_fdr"]]
        sel_pearson = float(stats.pearsonr(sel["slope_mssm"], sel["slope_hbcc"]).statistic) if len(sel) > 2 else np.nan
        sel_spearman = float(stats.spearmanr(sel["slope_mssm"], sel["slope_hbcc"]).statistic) if len(sel) > 2 else np.nan
        sel_same = int(sel["same_direction"].sum())
        sel_nominal = int((sel["same_direction"] & (sel["pvalue_hbcc"] < 0.05)).sum())
        both = int(gw["both_fdr"].sum())
        n_mssm_fdr = int(gw["mssm_fdr"].sum())
        n_hbcc_fdr = int(gw["hbcc_fdr"].sum())
        hyper_p = float(stats.hypergeom.sf(both - 1, len(gw), n_mssm_fdr, n_hbcc_fdr))
        hyper_exp = n_mssm_fdr * n_hbcc_fdr / len(gw)
        dec = pd.qcut(gw["mean_cp10k_mssm"], 10, labels=False, duplicates="drop")
        gw = gw.assign(decile=dec)
        exp_matched = 0.0
        for _, group in gw.groupby("decile", observed=True):
            exp_matched += float(group["mssm_fdr"].sum()) * float(group["hbcc_fdr"].mean())
        gidx_m = {gene: index for index, gene in enumerate(gm)}
        midx = np.array([gidx_m[gene] for gene in merged["gene"]])
        blocks = original.coexpression_blocks(cm[:, midx], merged["gene"].to_numpy(), seed=SEED)
        rng = np.random.default_rng(SEED)
        mssm_fdr_idx = np.flatnonzero(gw["mssm_fdr"].to_numpy())
        hbcc_fdr = gw["hbcc_fdr"].to_numpy()
        hbcc_sign = np.sign(gw["slope_hbcc"].to_numpy())
        mssm_sign = np.sign(gw["slope_mssm"].to_numpy())
        null_overlaps = []
        unique_blocks = np.unique(blocks)
        for _ in range(N_PERMS):
            perm_hbcc_fdr = hbcc_fdr.copy()
            perm_hbcc_sign = hbcc_sign.copy()
            for block in unique_blocks:
                idx = np.flatnonzero(blocks == block)
                perm = rng.permutation(idx)
                perm_hbcc_fdr[idx] = hbcc_fdr[perm]
                perm_hbcc_sign[idx] = hbcc_sign[perm]
            overlap = int(np.sum(
                perm_hbcc_fdr[mssm_fdr_idx]
                & (perm_hbcc_sign[mssm_fdr_idx] == mssm_sign[mssm_fdr_idx])
            ))
            null_overlaps.append(overlap)
        null_overlaps = np.asarray(null_overlaps)
        block_null_mean = float(null_overlaps.mean())
        block_p = float((1 + int(np.sum(null_overlaps >= both))) / (N_PERMS + 1))
        summary = {
            "cell_type": cell_type, "short_type": SHORT[cell_type], "n_tested_both": len(gw),
            "gw_pearson_r": pearson, "gw_spearman_rho": spearman,
            "gw_same_direction": same_dir, "gw_same_direction_frac": same_dir / len(gw),
            "gw_same_direction_binomial_p": binom_p, "mssm_fdr_genes": n_mssm_fdr,
            "hbcc_fdr_genes": n_hbcc_fdr, "sel_pearson_r": sel_pearson,
            "sel_spearman_rho": sel_spearman, "sel_same_direction": sel_same,
            "sel_same_direction_frac": sel_same / max(len(sel), 1),
            "sel_hbcc_nominal_same_direction": sel_nominal, "both_fdr_same_direction": both,
            "hypergeometric_expected": hyper_exp, "hypergeometric_p": hyper_p,
            "matched_background_expected": exp_matched, "block_perm_null_mean": block_null_mean,
            "block_perm_p": block_p, "block_perm_n": N_PERMS,
        }
        summaries.append(summary)
        original.emit_json("mature_replication_unit", summary)
        cm_sum = cm.sum(axis=0)
        ch_sum = ch.sum(axis=0)
        gidx_h = {gene: index for index, gene in enumerate(gh)}
        for _, row in merged[merged["mssm_fdr"]].iterrows():
            gi_m = gidx_m[row["gene"]]
            gi_h = gidx_h[row["gene"]]
            tier = (
                "both_FDR" if bool(row["both_fdr"])
                else "mssm_FDR_hbcc_nominal" if row["same_direction"] and row["pvalue_hbcc"] < 0.05
                else "mssm_FDR_only"
            )
            gene_rows.append({
                "cell_type": cell_type, "short_type": SHORT[cell_type], "gene": row["gene"], "tier": tier,
                "slope_mssm": row["slope_mssm"], "se_mssm": row["se_mssm"],
                "ci95_low_mssm": row["ci95_low_mssm"], "ci95_high_mssm": row["ci95_high_mssm"],
                "pvalue_mssm": row["pvalue_mssm"], "q_bh_mssm": row["q_bh_mssm"],
                "slope_hbcc": row["slope_hbcc"], "se_hbcc": row["se_hbcc"],
                "ci95_low_hbcc": row["ci95_low_hbcc"], "ci95_high_hbcc": row["ci95_high_hbcc"],
                "pvalue_hbcc": row["pvalue_hbcc"], "q_bh_hbcc": row["q_bh_hbcc"],
                "same_direction": bool(row["same_direction"]), "mssm_total_cells": int(cm_sum[gi_m]),
                "mssm_n_donors": len(mm), "mssm_median_cells_per_donor": float(mm["n_cells"].median()),
                "hbcc_total_cells": int(ch_sum[gi_h]), "hbcc_n_donors": len(mh),
                "hbcc_median_cells_per_donor": float(mh["n_cells"].median()),
                "mean_cp10k_mssm": row["mean_cp10k_mssm"],
                "detected_frac_mssm": row["detected_frac_mssm"],
            })
    return pd.DataFrame(summaries), pd.DataFrame(gene_rows)


def build_moderated_outputs(moderated: dict, scratch: Path, gate_pass: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_rows, de_rows = [], []
    for cell_type in INTERNEURON_TYPES:
        m = moderated[("MSSM_DLPFC", cell_type)]
        h = moderated[("HBCC_DLPFC", cell_type)]
        merged = m.merge(h, on="gene", suffixes=("_mssm", "_hbcc"))
        same = merged["slope_mssm"] * merged["slope_hbcc"] > 0
        mssm_fdr = merged["q_bh_mod_mssm"] < 0.05
        hbcc_fdr = merged["q_bh_mod_hbcc"] < 0.05
        both = same & mssm_fdr & hbcc_fdr
        cm, mm, _ = original.load_counts_npz(scratch, "mssm_mature30", cell_type)
        ch, mh, _ = original.load_counts_npz(scratch, "hbcc_mature30", cell_type)
        for cohort, frame, meta in (("MSSM_DLPFC", m, mm), ("HBCC_DLPFC", h, mh)):
            prior = {key: frame[key].iloc[0] for key in ("resid_dof", "d0_prior", "s0_sq_prior", "df_post")}
            de_rows.append({
                "cohort": cohort, "cell_type": cell_type, "n_donors": len(meta),
                "n_genes_tested": len(frame), "design": "intercept + age + sex_male + mean_log_umi",
                **prior,
                "n_de_unmod": int(frame["q_bh"].lt(0.05).sum()),
                "n_de_mod": int(frame["q_bh_mod"].lt(0.05).sum()),
                "n_de_unmod_fc025": int(((frame["q_bh"] < 0.05) & (frame["slope"].abs() / np.log(2) >= 0.25)).sum()),
                "n_de_mod_fc025": int(((frame["q_bh_mod"] < 0.05) & (frame["slope"].abs() / np.log(2) >= 0.25)).sum()),
            })
        for index in np.flatnonzero(mssm_fdr.to_numpy()):
            row = merged.iloc[index]
            tier = (
                "both_FDR" if bool(both.iloc[index])
                else "mssm_FDR_hbcc_nominal" if bool(same.iloc[index]) and row["pvalue_mod_hbcc"] < 0.05
                else "mssm_FDR_only"
            )
            gene_rows.append({
                "cell_type": cell_type, "short_type": SHORT[cell_type], "gene": row["gene"], "tier": tier,
                "slope_mssm": row["slope_mssm"], "se_mssm_mod": row["se_mod_mssm"],
                "pvalue_mssm_mod": row["pvalue_mod_mssm"], "q_bh_mssm_mod": row["q_bh_mod_mssm"],
                "slope_hbcc": row["slope_hbcc"], "se_hbcc_mod": row["se_mod_hbcc"],
                "pvalue_hbcc_mod": row["pvalue_mod_hbcc"], "q_bh_hbcc_mod": row["q_bh_mod_hbcc"],
                "same_direction": bool(same.iloc[index]), "mean_cp10k_mssm": row["mean_cp10k_mssm"],
                "mean_cp10k_hbcc": row["mean_cp10k_hbcc"],
            })
    genes = pd.DataFrame(gene_rows)
    de = pd.DataFrame(de_rows)
    genes.attrs["moderation_gate_pass"] = gate_pass
    return genes, de, build_six_gene_status(moderated, gate_pass)


def build_six_gene_status(moderated: dict, gate_pass: bool) -> pd.DataFrame:
    rows = []
    for gene in SHARED_GENES:
        for cell_type in SUBTYPES:
            for cohort in ("MSSM_DLPFC", "HBCC_DLPFC"):
                frame = moderated[(cohort, cell_type)]
                selected = frame[frame["gene"].eq(gene)]
                if selected.empty:
                    unmoderated_status = "NOT_TESTED"
                    moderated_status = "NOT_TESTED"
                    unmoderated_q = np.nan
                    moderated_q = np.nan
                else:
                    unmoderated_q = float(selected["q_bh"].iloc[0])
                    moderated_q = float(selected["q_bh_mod"].iloc[0])
                    unmoderated_status = "PASS" if unmoderated_q < 0.05 else "NOT_FDR"
                    moderated_status = "PASS" if moderated_q < 0.05 else "NOT_FDR"
                rows.append({
                    "gene": gene, "cell_type": cell_type, "short_type": SHORT[cell_type], "cohort": cohort,
                    "unmoderated_status": unmoderated_status, "unmoderated_q_bh": unmoderated_q,
                    "moderated_status": moderated_status, "moderated_q_bh": moderated_q,
                    "moderation_verification": "VERIFIED" if gate_pass else "UNVERIFIED_GATE_FAILED",
                })
    return pd.DataFrame(rows)


def bh_audit(complete_genes: pd.DataFrame, moderated: dict) -> pd.DataFrame:
    rows = []
    for (cohort, cell_type, model), frame in complete_genes.groupby(["cohort", "cell_type", "model"], sort=False):
        recomputed = multipletests(frame["pvalue"].to_numpy(float), method="fdr_bh")[1]
        difference = float(np.max(np.abs(recomputed - frame["q_bh"].to_numpy(float))))
        rows.append({"cohort": cohort, "cell_type": cell_type, "model": model, "n_genes": len(frame), "max_abs_q_difference": difference})
    for (cohort, cell_type), frame in moderated.items():
        recomputed = multipletests(frame["pvalue_mod"].to_numpy(float), method="fdr_bh")[1]
        difference = float(np.max(np.abs(recomputed - frame["q_bh_mod"].to_numpy(float))))
        rows.append({"cohort": cohort, "cell_type": cell_type, "model": "rev2_voom_wls_moderated", "n_genes": len(frame), "max_abs_q_difference": difference})
    audit = pd.DataFrame(rows)
    assert audit["max_abs_q_difference"].max() < 1e-9, audit.to_string()
    return audit


def fisher_mde80(n: int) -> float:
    return float(np.tanh((stats.norm.ppf(0.975) + stats.norm.ppf(0.80)) / np.sqrt(n - 3)))


def build_power(complete: pd.DataFrame, replication: pd.DataFrame) -> pd.DataFrame:
    frozen_se = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-all-celltypes.csv")
    frozen_de = pd.read_csv(RESULTS / "reviewer-rev2-human-de-moderated.csv")
    frozen_rep = pd.read_csv(RESULTS / "reviewer-rev2-replication-summary.csv")
    rows = []
    for cell_type in INTERNEURON_TYPES:
        for cohort in ("MSSM_DLPFC", "HBCC_DLPFC"):
            for age_range in ("full_range", "mature30"):
                if age_range == "full_range":
                    se_row = frozen_se[frozen_se["cohort"].eq(cohort) & frozen_se["cell_type"].eq(cell_type)].iloc[0]
                    de_row = frozen_de[frozen_de["cohort"].eq(cohort) & frozen_de["cell_type"].eq(cell_type)].iloc[0]
                    median_se = float(se_row["median_se"])
                    dof = int(de_row["resid_dof"])
                    n_donors = int(de_row["n_donors"])
                    n_shared = int(frozen_rep.loc[frozen_rep["cell_type"].eq(cell_type), "n_tested_both"].iloc[0])
                    source = "reviewer-rev2-interneuron-all-celltypes.csv median_se"
                else:
                    frame = complete[
                        complete["cohort"].eq(cohort) & complete["cell_type"].eq(cell_type)
                        & complete["model"].eq("rev2_voom_wls")
                    ]
                    median_se = float(frame["se"].median())
                    n_donors = int(frame["n_donors"].iloc[0])
                    dof = n_donors - 4
                    n_shared = int(replication.loc[replication["cell_type"].eq(cell_type), "n_tested_both"].iloc[0])
                    source = "followup-mature-interneuron.csv WLS SE"
                mde = (stats.t.ppf(0.975, dof) + stats.t.ppf(0.80, dof)) * median_se
                rows.append({
                    "cohort": cohort, "cell_type": cell_type, "short_type": SHORT[cell_type],
                    "age_range": age_range, "n_donors": n_donors, "resid_dof": dof,
                    "n_tested_genes": int(frozen_de.loc[frozen_de["cohort"].eq(cohort) & frozen_de["cell_type"].eq(cell_type), "n_genes_tested"].iloc[0]) if age_range == "full_range" else len(frame),
                    "median_wls_se": median_se, "mde80_slope_per_year": float(mde),
                    "mde80_slope_per_10_years": float(10 * mde), "n_shared_genes": n_shared,
                    "mde80_pearson_r": fisher_mde80(n_shared), "se_source": source,
                })
    return pd.DataFrame(rows)


def build_discovery(complete: pd.DataFrame) -> pd.DataFrame:
    output = RESULTS / "followup-mature-interneuron-discovery-ols.csv"
    rows = []
    for cell_type in INTERNEURON_TYPES:
        frame = complete[complete["cell_type"].eq(cell_type) & complete["model"].eq("age_plus_sex_depth_ols")]
        m = frame[frame["cohort"].eq("MSSM_DLPFC")]
        h = frame[frame["cohort"].eq("HBCC_DLPFC")]
        joined = m.merge(h, on="gene", suffixes=("_mssm", "_hbcc"), validate="one_to_one")
        rows.append(pd.DataFrame({
            "cell_type": cell_type, "model": "age_plus_sex_depth", "gene": joined["gene"],
            "slope_log_cp10k_per_year_mssm": joined["slope_mssm"],
            "pvalue_mssm": joined["pvalue_mssm"], "q_bh_mssm": joined["q_bh_mssm"],
            "slope_log_cp10k_per_year_hbcc": joined["slope_hbcc"],
            "pvalue_hbcc": joined["pvalue_hbcc"], "q_bh_hbcc": joined["q_bh_hbcc"],
        }))
    discovery = pd.concat(rows, ignore_index=True)
    write_clean_csv(discovery, output)
    return discovery


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    args = parser.parse_args()
    started = time.time()
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    assert MSSM_PATH.exists()
    hbcc_validation = validate_hbcc()
    original.log(json.dumps({"hbcc_validation": hbcc_validation}, default=float))
    block_csv = pd.read_csv(RESULTS / "followup-11pt-block9-mature-interneuron.csv")
    complete_frames, summaries, de_fits = [], [], {}
    mssm_cache_rows = []
    for cell_type in INTERNEURON_TYPES:
        try:
            frame, summary, de = validate_mssm_unit(cell_type, args.scratch_dir, block_csv)
        except (AssertionError, FileNotFoundError, KeyError, ValueError) as exc:
            original.log(f"MSSM cache validation failed for {cell_type}: {exc}; refitting")
            mssm_cache_rows.append({"cell_type": cell_type, "status": "REFIT", "reason": str(exc)})
            frame, summary, de = refit_mssm_unit(cell_type, args.scratch_dir)
        else:
            mssm_cache_rows.append({"cell_type": cell_type, "status": "PASS", "reason": "none"})
        complete_frames.append(frame)
        summaries.append(summary)
        de_fits[("MSSM_DLPFC", cell_type)] = de
    for cell_type in INTERNEURON_TYPES:
        result = fit_hbcc_unit(cell_type, args.scratch_dir)
        if result is None:
            summaries.append({
                "cohort": "HBCC_DLPFC", "cell_type": cell_type, "short_type": SHORT[cell_type],
                "status": "UNTESTABLE", "status_reason": "fewer_than_10_mature_donors", "n_donors": 0,
            })
            raise RuntimeError(f"HBCC {cell_type} is UNTESTABLE; downstream replication cannot be completed")
        frame, summary, de = result
        complete_frames.append(frame)
        summaries.append(summary)
        de_fits[("HBCC_DLPFC", cell_type)] = de
    complete = pd.concat(complete_frames, ignore_index=True).reindex(columns=GENE_COLUMNS)
    write_clean_csv(complete, GENE_OUTPUT)
    gate, gate_pass = recover_moderation_gate()
    write_clean_csv(gate, GATE_OUTPUT)
    original.log(f"Moderation recovery gate pass: {gate_pass}")
    moderated = {}
    s2_validations = []
    for (cohort, cell_type), de in de_fits.items():
        tag = "mssm_mature30" if cohort == "MSSM_DLPFC" else "hbcc_mature30"
        s2, dof, validation = fit_s2_and_validate(cell_type, tag, de, args.scratch_dir)
        s2_validations.append(validation)
        fitted, _ = moderate_fit(de, s2, dof)
        moderated[(cohort, cell_type)] = fitted
    moderated_genes, moderated_de, six_gene = build_moderated_outputs(moderated, args.scratch_dir, gate_pass)
    write_clean_csv(moderated_genes, MODERATED_GENES_OUTPUT)
    write_clean_csv(moderated_de, MODERATED_DE_OUTPUT)
    write_clean_csv(six_gene, SIX_GENE_OUTPUT)
    replication, genes = build_replication(args.scratch_dir)
    write_clean_csv(replication, REPLICATION_OUTPUT)
    write_clean_csv(genes, GENES_OUTPUT)
    stored_complete = pd.read_csv(GENE_OUTPUT)
    audit = bh_audit(stored_complete, moderated)
    write_clean_csv(audit, BH_OUTPUT)
    discovery = build_discovery(stored_complete)
    power = build_power(stored_complete, replication)
    write_clean_csv(power, POWER_OUTPUT)
    outputs = [
        GENE_OUTPUT, REPLICATION_OUTPUT, GENES_OUTPUT, MODERATED_GENES_OUTPUT,
        MODERATED_DE_OUTPUT, DISCOVERY_OUTPUT, POWER_OUTPUT, SIX_GENE_OUTPUT, GATE_OUTPUT, BH_OUTPUT,
    ]
    payload = {
        "analysis": "Item 7 mature interneuron replication",
        "seed": SEED, "block_permutations": N_PERMS, "mature_age_minimum": MATURE_AGE,
        "hbcc_validation": hbcc_validation, "mssm_cache_validation": mssm_cache_rows,
        "units": summaries, "moderation_gate_pass": gate_pass,
        "moderation_gate": gate.to_dict("records"),
        "moderated_results_verified": gate_pass,
        "s2_function_validation": s2_validations,
        "bh_max_abs_difference": float(audit["max_abs_q_difference"].max()),
        "mature_discovery_sha256": sha256(DISCOVERY_OUTPUT),
        "script_sha256": sha256(Path(__file__)),
        "runtime_seconds": time.time() - started,
        "outputs": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size} for path in outputs},
        "discovery_rows": len(discovery),
    }
    (args.scratch_dir / "run_summary_mature.json").write_text(json.dumps(payload, indent=2) + "\n")
    original.log(f"Wrote {len(outputs)} mature result files and run_summary_mature.json")


if __name__ == "__main__":
    main()
