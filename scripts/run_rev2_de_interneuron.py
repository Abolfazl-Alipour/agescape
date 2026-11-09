#!/usr/bin/env python3
"""Rev2 fixes 8-9: pseudobulk DE upgrade (voom-style precision-weighted
regression) and human interneuron replication redo.

Reanalysis worker for the Agescape Rev2 revision. CPU only; data already in
the repo or re-pulled from CELLxGENE Census 2025-11-08 (same source and
version as the original pipelines). No new experiments, no new data
collections. Writes new scripts/result files only; never touches .tex or
report/manuscript-aging-mechanisms.md; never commits or pushes.

Model (documented in full in results/reviewer-rev2-de-methods.md):
  * Donor pseudobulk = summed raw counts per donor x cell type
    (selection identical to the original pipelines: primary data, disease
    normal, cell raw_sum >= 1500, >= 20 cells per donor-cell-type, age
    20-120 y from development_stage).
  * Low-expression filter (explicit, new): drop technical genes
    (MT-/RP[LS]/MRP[LS] prefixes, XIST/TSIX/DDX3Y/EIF2S3Y/UTY/KDM5D) and
    keep genes with total pseudobulk counts >= 10 across donors AND
    detected (count > 0) in >= 10% of donors AND mean CP10k >= 1.
  * Response: log1p(CP10k) of the summed counts (same scale as the original
    t-regressions, so slopes are directly comparable).
  * Precision model (limma-voom style, implemented in numpy because rpy2 /
    limma are unavailable here): per-gene residual variance from an OLS
    pass is smoothed against mean logCPM (60 quantile bins, median
    sqrt-variance, moving-average smooth, linear interpolation) to give a
    fitted variance v_g; precision weights are
    w_gd = (n_cells_d / median n_cells) / v_g, clipped so the within-gene
    weight ratio stays in [0.25, 4]. The trend is refit once from WLS
    residuals (one voom iteration). Weighted least squares is then fit by
    batched normal equations. This is a precision-weighted linear model on
    log counts, NOT a negative-binomial GLM; the mouse arm additionally
    runs a statsmodels negative-binomial GLM on signature genes as a
    cross-check, and the original PyDESeq2 NB-GLM tables are audited
    (PyDESeq2 itself is not installed in this environment).
  * Design (human): logCPM ~ 1 + age + sex_male + mean_log_umi
    (mean_log_umi = donor mean log library size, the original depth
    covariate). Sex is dropped if constant; depth is dropped only if the
    design is rank-deficient. BH q within each cohort x cell type over the
    filtered tested genes.
  * Design (mouse): logCPM ~ 1 + age_group18m + sex, on mouse
    pseudobulk (TMS droplet, Lung + Kidney, 3m vs 18m).
  * Leave-one-donor-out (human) / leave-one-mouse-out (mouse) sensitivity
    for headline genes: refit WLS dropping each donor/mouse once; report
    slope range, max |delta slope| / SE, and how many drops flip the sign
    or lose nominal significance.

Outputs (results/):
  reviewer-rev2-interneuron-genes.csv       replicated gene list
  reviewer-rev2-interneuron-figure-data.csv Figure 4 data, all tested genes
  reviewer-rev2-interneuron-all-celltypes.csv per-cell-type DEG counts
  reviewer-rev2-de-methods.md               (written separately, from these)
Scratch (resume/audit): --scratch-dir, default ~/workspace/rev2-de-scratch,
  one .npz per cohort x cell type with raw counts, metadata, and DE tables.

Every stage prints a JSON summary to stdout (sandbox writes flap EPERM, so
stdout is the durable record) and saves files with retries.

Usage:
  python scripts/run_rev2_de_interneuron.py --stage probe
  python scripts/run_rev2_de_interneuron.py --stage human-de
  python scripts/run_rev2_de_interneuron.py --stage replication
  python scripts/run_rev2_de_interneuron.py --stage lodo
  python scripts/run_rev2_de_interneuron.py --stage all-celltypes
  python scripts/run_rev2_de_interneuron.py --stage mouse
  python scripts/run_rev2_de_interneuron.py --stage psychad
  python scripts/run_rev2_de_interneuron.py --stage emit --artifact genes
  python scripts/run_rev2_de_interneuron.py --stage emit --artifact figure-data --chunk 0 --n-chunks 8
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

RESULTS = REPO / "results"
CENSUS_VERSION = "2025-11-08"
MSSM_ID = "37a17b78-4864-4a42-b67b-31c00962795a"
HBCC_ID = "5c97eeeb-7e52-44b3-b010-b832b1f5424c"
AGE_RE = re.compile(r"(\d+)-year-old")
MIN_CELLS_PER_DONOR_CT = 20
MIN_DONORS = 10
SEED = 20261108

INTERNEURON_TYPES = [
    "VIP GABAergic cortical interneuron",
    "pvalb GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron",
    "GABAergic neuron",
]
SHORT = {
    "VIP GABAergic cortical interneuron": "VIP",
    "pvalb GABAergic cortical interneuron": "pvalb",
    "sst GABAergic cortical interneuron": "sst",
    "GABAergic neuron": "GABA",
}

TECHNICAL_RE = re.compile(r"^(MT-|RP[LS]|MRP[LS])")


def technical_gene_mask(gene_names):
    genes = pd.Index(gene_names).astype(str).str.upper()
    return np.array([
        bool(TECHNICAL_RE.match(g))
        or g in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for g in genes
    ])


def log(*args):
    print(*args, flush=True)


def tiledb_config_from_env():
    """TileDB S3 config for the sandbox egress proxy.

    Adapted from the sibling Rev2 worker's probe
    (~/workspace/rev2-mssm-null-scratch/probe_soma_proxy.py): TileDB's S3
    client does not read the proxy env vars, so they are passed explicitly
    as vfs.s3.* settings (host/port/scheme/user/password from HTTPS_PROXY,
    CA bundle from CURL_CA_BUNDLE). Credentials are never hardcoded or
    logged.
    """
    import os
    from urllib.parse import urlparse

    cfg = {"vfs.s3.region": "us-west-2"}
    proxy = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
    if proxy:
        u = urlparse(proxy)
        cfg["vfs.s3.proxy_host"] = u.hostname
        cfg["vfs.s3.proxy_port"] = str(u.port or 3128)
        cfg["vfs.s3.proxy_scheme"] = u.scheme or "http"
        if u.username:
            cfg["vfs.s3.proxy_username"] = u.username
        if u.password:
            cfg["vfs.s3.proxy_password"] = u.password
    ca = os.environ.get("CURL_CA_BUNDLE") or os.environ.get("AWS_CA_BUNDLE")
    if ca and os.path.exists(ca):
        cfg["vfs.s3.ca_file"] = ca
    # Memory: TileDB defaults to many parallel S3 readers; on this 8 GB VM
    # with sibling jobs running, a whole-cell-type pull OOM-killed twice
    # (exit 137). Cap reader/concurrency so S3 buffers stay small.
    # Raised from 2 to 8/4 after the 2-stream setting proved too slow:
    # a single VIP pull had not finished in 20 min while the sibling's
    # 16-stream h5ad download held most of the ~4 MB/s S3 pipe. Hellinger
    # and transport-null have since finished, so the headroom exists now;
    # watch RSS on the first unit and drop back to 2 if it nears 2 GB.
    cfg["sm.io_concurrency_level"] = "4"
    cfg["sm.compute_concurrency_level"] = "4"
    cfg["vfs.s3.max_parallel_ops"] = "8"
    # Timeouts: the egress proxy stalls S3 reads; TileDB defaults
    # (3 s connect/request timeout) abort fragment reads with curlCode 28
    # ("Operation too slow"). Raise both and retry more.
    cfg["vfs.s3.connect_timeout_ms"] = "120000"
    cfg["vfs.s3.request_timeout_ms"] = "600000"
    cfg["vfs.s3.max_retries"] = "10"
    cfg["vfs.s3.retry_delay_ms"] = "2000"
    cfg["vfs.s3.connect_max_tries"] = "5"
    return cfg


def open_census():
    """Open CELLxGENE Census 2025-11-08 through the egress proxy."""
    import cellxgene_census
    return cellxgene_census.open_soma(
        census_version=CENSUS_VERSION, tiledb_config=tiledb_config_from_env())


def open_census_retry(max_attempts: int = 24, sleep_s: float = 300.0):
    """Open Census, retrying through transient S3/proxy outages.

    The egress proxy's S3 path intermittently stalls (curlCode 28 low-speed
    timeouts, or 60 s with 0 bytes on a plain curl). When it is down, retry
    every 5 min for up to 2 h rather than failing the whole stage.
    """
    import time
    last = None
    for attempt in range(1, max_attempts + 1):
        try:
            census = open_census()
            log(f"Census open OK on attempt {attempt}/{max_attempts}")
            return census
        except Exception as e:
            last = e
            log(f"Census open attempt {attempt}/{max_attempts} failed: "
                f"{type(e).__name__}: {str(e)[:200]}")
            if attempt < max_attempts:
                log(f"  sleeping {sleep_s:.0f}s before retry")
                time.sleep(sleep_s)
    raise last


def save_npz_retry(path: Path, tries: int = 8, **arrays) -> bool:
    """np.savez_compressed with retries: sandbox file writes flap EPERM."""
    last = None
    for attempt in range(1, tries + 1):
        try:
            np.savez_compressed(path, **arrays)
            return True
        except OSError as exc:
            last = exc
            log(f"  save {path.name} attempt {attempt} failed ({exc}); retrying")
            time.sleep(10 * attempt)
    log(f"  save {path.name} FAILED after {tries} attempts: {last}")
    return False


def save_csv_retry(df: pd.DataFrame, path: Path, tries: int = 8) -> bool:
    last = None
    for attempt in range(1, tries + 1):
        try:
            df.to_csv(path, index=False)
            return True
        except OSError as exc:
            last = exc
            log(f"  save {path.name} attempt {attempt} failed ({exc}); retrying")
            time.sleep(10 * attempt)
    log(f"  save {path.name} FAILED after {tries} attempts: {last}")
    return False


def emit_json(tag: str, payload: dict) -> None:
    print(f"@@JSON {tag} " + json.dumps(payload, default=float), flush=True)

# ---------------------------------------------------------------------------
# Pseudobulk + filtering + voom-style WLS
# ---------------------------------------------------------------------------

def pseudobulk_counts(adata):
    """Summed raw counts per donor. Selection mirrors the original pipelines.

    Expects adata.obs with donor_id, development_stage, sex and adata.X raw
    counts. Returns (counts donors x genes int64, metadata DataFrame,
    gene names). Donors with < 20 cells are dropped; age restricted to
    20-120 y.
    """
    adata.obs["donor_id"] = adata.obs["donor_id"].astype(str)
    adata.obs["age"] = adata.obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
    adata = adata[adata.obs["age"].between(20, 120)].copy()
    adata.var["symbol"] = adata.var["feature_name"].astype(str)
    adata.var_names = adata.var["symbol"].values
    adata.var_names_make_unique()
    matrix = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    count_rows, meta_rows = [], []
    for donor in adata.obs["donor_id"].unique():
        mask = (adata.obs["donor_id"] == donor).to_numpy()
        if mask.sum() < MIN_CELLS_PER_DONOR_CT:
            continue
        dm = matrix[mask]
        count_rows.append(np.asarray(dm.sum(axis=0)).ravel())
        lib = np.asarray(dm.sum(axis=1)).ravel()
        meta_rows.append({
            "donor": donor,
            "age": float(adata.obs.loc[mask, "age"].iloc[0]),
            "sex_male": float(adata.obs.loc[mask, "sex"].astype(str).str.lower().iloc[0] == "male"),
            "n_cells": int(mask.sum()),
            "mean_log_umi": float(np.log(lib).mean()),
            "library": float(lib.sum()),
        })
    counts = np.vstack(count_rows).astype(np.int64)
    metadata = pd.DataFrame(meta_rows)
    return counts, metadata, adata.var_names.to_numpy()


def filter_genes(counts, gene_names):
    """Explicit low-expression filter (see module docstring).

    Returns keep mask and a per-gene stats frame.
    """
    lib = counts.sum(axis=1, keepdims=True).astype(float)
    cp10k = counts / lib * 1e4
    total = counts.sum(axis=0)
    detected_frac = (counts > 0).mean(axis=0)
    mean_cp10k = cp10k.mean(axis=0)
    tech = technical_gene_mask(gene_names)
    keep = (~tech) & (total >= 10) & (detected_frac >= 0.10) & (mean_cp10k >= 1.0)
    gene_stats = pd.DataFrame({
        "gene": gene_names,
        "technical": tech,
        "total_counts": total,
        "detected_frac": detected_frac,
        "mean_cp10k": mean_cp10k,
        "kept": keep,
    })
    return keep, gene_stats


def design_matrix(metadata, columns):
    """Design matrix with intercept; drops constant columns (sex if needed)."""
    cols = [c for c in columns if metadata[c].nunique() > 1] if len(metadata) > 1 else list(columns)
    X = np.column_stack([np.ones(len(metadata))] + [metadata[c].astype(float).to_numpy() for c in cols])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        # drop the last non-age column and retry (mirrors original fallback)
        for drop in reversed(range(len(cols))):
            if cols[drop] == "age":
                continue
            cols2 = cols[:drop] + cols[drop + 1:]
            X2 = np.column_stack([np.ones(len(metadata))] + [metadata[c].astype(float).to_numpy() for c in cols2])
            if np.linalg.matrix_rank(X2) == X2.shape[1]:
                return X2, ["intercept"] + cols2
        raise ValueError("rank-deficient design even after dropping columns")
    return X, ["intercept"] + cols


def _fit_trend(mean_logcpm, resid_var, n_bins=60):
    """Voom-style mean-variance trend: binned median sqrt-variance, smoothed."""
    order = np.argsort(mean_logcpm)
    m_sorted = mean_logcpm[order]
    v_sorted = np.sqrt(np.maximum(resid_var[order], 1e-12))
    bins = np.array_split(np.arange(len(order)), n_bins)
    centers, vals = [], []
    for b in bins:
        if len(b) == 0:
            continue
        centers.append(np.median(m_sorted[b]))
        vals.append(np.median(v_sorted[b]))
    centers = np.array(centers)
    vals = np.array(vals)
    # centered moving average (window 5) for smoothing
    k = 5
    pad = np.concatenate([[vals[0]] * (k // 2), vals, [vals[-1]] * (k // 2)])
    smooth = np.convolve(pad, np.ones(k) / k, mode="valid")
    fitted_sd = np.interp(mean_logcpm, centers, smooth)
    return fitted_sd ** 2, centers, smooth


def _wls_batched(Y, X, W):
    """Batched WLS for many genes sharing one design.

    Y: genes x donors logCPM; X: donors x k; W: genes x donors weights.
    Returns beta (genes x k), se (genes x k), resid_var (genes,).
    """
    G, D = Y.shape
    k = X.shape[1]
    # A_g = X' W_g X ; b_g = X' W_g y_g
    A = np.einsum("gd,dk,dl->gkl", W, X, X)
    b = np.einsum("gd,dk,gd->gk", W, X, Y)
    Ainv = np.linalg.inv(A)
    beta = np.einsum("gkl,gl->gk", Ainv, b)
    resid = Y - beta @ X.T
    rss_w = np.einsum("gd,gd->g", W, resid ** 2)
    dof = D - k
    s2 = rss_w / max(dof, 1)
    se = np.sqrt(np.maximum(s2[:, None] * np.diagonal(Ainv, axis1=1, axis2=2), 0.0))
    return beta, se, s2, dof


def voom_wls_counts(counts, metadata, gene_names, columns=("age", "sex_male", "mean_log_umi")):
    """Full voom-WLS fit on raw pseudobulk counts (donors x genes)."""
    keep, gene_stats = filter_genes(counts, gene_names)
    C = counts[:, keep].astype(float)
    genes = gene_names[keep]
    lib = counts.sum(axis=1).astype(float)  # total library (all genes), as in the original
    Y = np.log1p(C / lib[:, None] * 1e4)  # genes x donors after transpose
    Y = Y.T  # genes x donors
    n_cells = metadata["n_cells"].to_numpy(float)
    cell_factor = n_cells / np.median(n_cells)
    X, colnames = design_matrix(metadata, list(columns))
    age_idx = colnames.index("age")
    D, k = X.shape

    def ols_pass(Ym):
        beta, *_ = np.linalg.lstsq(X, Ym.T, rcond=None)
        resid = Ym.T - X @ beta
        dof = D - k
        s2 = np.sum(resid ** 2, axis=0) / max(dof, 1)
        return s2

    mean_logcpm = Y.mean(axis=1)
    # pass 1: OLS trend
    s2_ols = ols_pass(Y)
    v_g, _, _ = _fit_trend(mean_logcpm, s2_ols)
    W = cell_factor[None, :] / v_g[:, None]
    # clip within-gene weight ratios to [0.25, 4]
    med = np.median(W, axis=1, keepdims=True)
    W = np.clip(W, 0.25 * med, 4.0 * med)
    beta, se, s2_w, dof = _wls_batched(Y, X, W)
    # pass 2: refit trend from WLS residuals, refit once (one voom iteration)
    resid_w = Y - beta @ X.T
    s2_wls = np.sum(resid_w ** 2, axis=1) / max(dof, 1)
    v_g2, _, _ = _fit_trend(mean_logcpm, s2_wls)
    W2 = cell_factor[None, :] / v_g2[:, None]
    med2 = np.median(W2, axis=1, keepdims=True)
    W2 = np.clip(W2, 0.25 * med2, 4.0 * med2)
    beta, se, s2_w, dof = _wls_batched(Y, X, W2)

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
    out.attrs["n_genes_technical_dropped"] = int(technical_gene_mask(gene_names).sum())
    out.attrs["design"] = " + ".join(colnames)
    out.attrs["dof"] = int(dof)
    out.attrs["gene_stats"] = gene_stats
    return out


def ols_arm(counts, metadata, gene_names, columns=("age", "sex_male", "mean_log_umi")):
    """Original-method OLS arm on the same filtered genes (head-to-head)."""
    keep, gene_stats = filter_genes(counts, gene_names)
    C = counts[:, keep].astype(float)
    genes = gene_names[keep]
    lib = counts.sum(axis=1).astype(float)  # total library (all genes), as in the original
    Y = np.log1p(C / lib[:, None] * 1e4).T  # genes x donors
    X, colnames = design_matrix(metadata, list(columns))
    age_idx = colnames.index("age")
    D, k = X.shape
    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    resid = Y.T - X @ beta
    dof = D - k
    s2 = np.sum(resid ** 2, axis=0) / max(dof, 1)
    cov_age = np.linalg.pinv(X.T @ X)[age_idx, age_idx]
    slope = beta[age_idx]
    se_age = np.sqrt(s2 * cov_age)
    tstat = np.divide(slope, se_age, out=np.zeros_like(slope), where=se_age > 0)
    pvalue = 2 * stats.t.sf(np.abs(tstat), dof)
    qvalue = multipletests(pvalue, method="fdr_bh")[1]
    tcrit = stats.t.ppf(0.975, dof)
    return pd.DataFrame({
        "gene": genes,
        "slope_ols": slope,
        "se_ols": se_age,
        "pvalue_ols": pvalue,
        "q_bh_ols": qvalue,
        "ci95_low_ols": slope - tcrit * se_age,
        "ci95_high_ols": slope + tcrit * se_age,
    })


def lodo_sensitivity(counts, metadata, gene_names, gene_list, columns=("age", "sex_male", "mean_log_umi")):
    """Leave-one-donor-out WLS for headline genes.

    Returns per-gene: full slope/se, min/max LODO slope, max |delta slope|/SE,
    n drops where the sign flips, n drops where nominal p > 0.05.
    """
    keep, _ = filter_genes(counts, gene_names)
    gidx = {g: i for i, g in enumerate(gene_names)}
    sel = [gidx[g] for g in gene_list if g in gidx and keep[gidx[g]]]
    genes = [gene_names[i] for i in sel]
    if not genes:
        return pd.DataFrame(columns=["gene", "slope_full", "se_full", "lodo_min", "lodo_max",
                                     "max_abs_delta_over_se", "n_sign_flip", "n_lose_nominal"])
    C = counts[:, sel].astype(float)
    lib = counts.sum(axis=1).astype(float)  # total library (all genes)
    # NB: library sizes from ALL genes (same as main fit)
    Y = np.log1p(C / lib[:, None] * 1e4).T  # genes x donors
    n_cells = metadata["n_cells"].to_numpy(float)
    cell_factor = n_cells / np.median(n_cells)
    X, colnames = design_matrix(metadata, list(columns))
    age_idx = colnames.index("age")
    D, k = X.shape
    mean_logcpm = Y.mean(axis=1)
    beta_ols, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    resid = Y.T - X @ beta_ols
    s2 = np.sum(resid ** 2, axis=0) / max(D - k, 1)
    v_g, _, _ = _fit_trend(mean_logcpm, s2)
    W = cell_factor[None, :] / v_g[:, None]
    med = np.median(W, axis=1, keepdims=True)
    W = np.clip(W, 0.25 * med, 4.0 * med)
    beta, se, s2_w, dof = _wls_batched(Y, X, W)
    slope_full = beta[:, age_idx]
    se_full = se[:, age_idx]
    tcrit = stats.t.ppf(0.975, dof)
    rows = []
    for d in range(D):
        keep_d = np.ones(D, dtype=bool)
        keep_d[d] = False
        b_d, se_d, _, _ = _wls_batched(Y[:, keep_d], X[keep_d], W[:, keep_d])
        rows.append((b_d[:, age_idx], se_d[:, age_idx]))
    slopes_d = np.array([r[0] for r in rows]).T  # genes x donors
    ses_d = np.array([r[1] for r in rows]).T
    dof_d = D - 1 - k
    t_d = np.divide(slopes_d, ses_d, out=np.zeros_like(slopes_d), where=ses_d > 0)
    p_d = 2 * stats.t.sf(np.abs(t_d), dof_d)
    out_rows = []
    for i, g in enumerate(genes):
        sd = slopes_d[i]
        out_rows.append({
            "gene": g,
            "slope_full": slope_full[i],
            "se_full": se_full[i],
            "lodo_min": float(sd.min()),
            "lodo_max": float(sd.max()),
            "max_abs_delta_over_se": float(np.max(np.abs(sd - slope_full[i]) / se_full[i])),
            "n_sign_flip": int(np.sum(np.sign(sd) != np.sign(slope_full[i]))),
            "n_lose_nominal": int(np.sum(p_d[i] > 0.05)),
            "n_donors": D,
        })
    return pd.DataFrame(out_rows)

# ---------------------------------------------------------------------------
# Census fetching + human DE stages
# ---------------------------------------------------------------------------

def fetch_unit(census, dataset_id, cell_type, scratch: Path, tag: str, donor_batch: int = 100000):
    """Fetch one cohort x cell type from Census; cache raw counts as .npz.

    Selection mirrors the original pipelines exactly (primary data, disease
    normal, cell raw_sum >= 1500, this cell type, age 20-120 y, >= 20 cells
    per donor). One get_anndata call per unit (donor_batch effectively
    unlimited): 5-donor batches were tried first and abandoned, each call
    re-scans S3 fragments and the first 5-donor batch had not finished in
    7 min, while the original pipelines' whole-unit pulls did complete.
    Peak RSS for a whole-unit pull is ~2 GB; the earlier OOM kill happened
    only while sibling jobs held ~4.7 GB. Counts are cached per unit, so a
    restart resumes at the next unit.
    Returns (counts, metadata, gene_names, cache_path) or
    (None, None, None, None) if the unit has too few donors.
    """
    safe = f"{tag}_{cell_type}".replace(" ", "_").replace("/", "_")
    cache = scratch / f"counts_{safe}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        meta = pd.DataFrame({k: z[k] for k in ("donor", "age", "sex_male", "n_cells", "mean_log_umi", "library")})
        return z["counts"], meta, np.array([str(g) for g in z["genes"]]), cache
    import cellxgene_census
    obs = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=(f"dataset_id == '{dataset_id}' and is_primary_data == True "
                      f"and disease == 'normal' and raw_sum >= 1500 "
                      f"and cell_type == '{cell_type}'"),
        column_names=["donor_id", "development_stage", "sex", "raw_sum"],
    )
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["age"] = obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
    obs = obs[obs["age"].between(20, 120)].copy()
    per_donor = obs.groupby("donor_id", observed=True).agg(
        n_cells=("donor_id", "size"),
        age=("age", "first"),
        sex=("sex", lambda s: str(s.iloc[0]))).reset_index()
    donors = sorted(per_donor.loc[per_donor["n_cells"] >= MIN_CELLS_PER_DONOR_CT, "donor_id"])
    if len(donors) < MIN_DONORS:
        log(f"  [{tag} | {cell_type}] too few donors ({len(donors)} < {MIN_DONORS}); skipped")
        return None, None, None, None
    donor_meta = per_donor.set_index("donor_id").loc[donors]
    counts = None
    genes = None
    meta_rows = []
    import gc
    import time
    for start in range(0, len(donors), donor_batch):
        batch = donors[start:start + donor_batch]
        donor_list = ",".join(f"'{d}'" for d in batch)
        filt = (f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                f"and raw_sum >= 1500 and donor_id in [{donor_list}] and cell_type == '{cell_type}'")
        adata = None
        for attempt in range(1, 4):
            try:
                adata = cellxgene_census.get_anndata(
                    census, "Homo sapiens",
                    obs_value_filter=filt,
                    obs_column_names=["donor_id", "development_stage", "sex"],
                    var_column_names=["feature_name"],
                )
                break
            except Exception as e:
                log(f"  [{tag} | {cell_type}] batch {start + 1}-{min(start + donor_batch, len(donors))} "
                    f"attempt {attempt}/3 failed: {type(e).__name__}: {e}")
                if attempt == 3:
                    raise
                time.sleep(10 * attempt)
        bc, bm, bg = pseudobulk_counts(adata)
        del adata
        import gc
        gc.collect()
        if genes is None:
            genes = bg
            counts = np.zeros((len(donors), len(genes)), dtype=np.int64)
        else:
            assert list(bg) == list(genes), f"var changed between batches for {cell_type}"
        pos = {d: i for i, d in enumerate(donors)}
        for j, d in enumerate(bm["donor"]):
            counts[pos[d]] = bc[j]
        meta_rows.append(bm)
        log(f"  [{tag} | {cell_type}] donors {start + 1}-{min(start + donor_batch, len(donors))}"
            f" of {len(donors)}", flush=True)
    meta = pd.concat(meta_rows, ignore_index=True)
    meta = meta.set_index("donor").loc[donors].reset_index()
    log(f"  [{tag} | {cell_type}] {len(meta)} donors, {len(genes)} genes, "
        f"{int(counts.sum())} total counts")
    ok = save_npz_retry(cache, counts=counts, genes=np.array(genes),
                        donor=meta["donor"].to_numpy(),
                        age=meta["age"].to_numpy(float),
                        sex_male=meta["sex_male"].to_numpy(float),
                        n_cells=meta["n_cells"].to_numpy(np.int64),
                        mean_log_umi=meta["mean_log_umi"].to_numpy(float),
                        library=meta["library"].to_numpy(float))
    if not ok:
        log("  WARNING: counts cache save failed; continuing in-memory")
    return counts, meta, genes, cache


def fetch_unit_h5ad(h5ad_path: Path, cell_type: str, scratch: Path, tag: str):
    """Same unit fetch as fetch_unit, but from a local Census h5ad (backed).

    The SOMA X-read path through the egress proxy stalls (a single VIP pull
    unfinished after 20 min, twice; the sibling Rev2 worker hit the same
    wall and switched to h5ad downloads). The h5ad is the same Census
    2025-11-08 dataset, same cells/genes; only the access path differs.
    Selection mirrors fetch_unit exactly (disease normal, raw_sum >= 1500,
    this cell type, age 20-120 y, >= 20 cells per donor). Columns that the
    h5ad obs lacks are handled by fallback: is_primary_data is skipped (the
    h5ad IS the primary data); raw_sum, if absent, is computed from X after
    the slice is loaded. Cache filename is identical to fetch_unit's, so a
    unit completed by either path is reused by the other.
    Returns (counts, metadata, gene_names, cache_path) or
    (None, None, None, None) if the unit has too few donors.
    """
    import anndata as ad
    safe = f"{tag}_{cell_type}".replace(" ", "_").replace("/", "_")
    cache = scratch / f"counts_{safe}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        meta = pd.DataFrame({k: z[k] for k in ("donor", "age", "sex_male", "n_cells", "mean_log_umi", "library")})
        return z["counts"], meta, np.array([str(g) for g in z["genes"]]), cache
    b = ad.read_h5ad(h5ad_path, backed="r")
    obs = b.obs
    cols = set(obs.columns)
    mask = pd.Series(True, index=obs.index)
    mask &= obs["cell_type"].astype(str) == cell_type
    if "disease" in cols:
        mask &= obs["disease"].astype(str) == "normal"
    if "is_primary_data" in cols:
        mask &= obs["is_primary_data"].astype(str).str.lower() == "true"
    if "raw_sum" in cols:
        mask &= pd.to_numeric(obs["raw_sum"], errors="coerce") >= 1500
    age = obs["development_stage"].astype(str).str.extract(AGE_RE).iloc[:, 0].astype(float)
    mask &= age.between(20, 120).to_numpy()
    sub_obs = obs[mask]
    per_donor = sub_obs.groupby("donor_id", observed=True).size()
    donors = sorted(per_donor[per_donor >= MIN_CELLS_PER_DONOR_CT].index.astype(str))
    if len(donors) < MIN_DONORS:
        log(f"  [{tag} | {cell_type}] too few donors ({len(donors)} < {MIN_DONORS}); skipped")
        return None, None, None, None
    keep = mask & obs["donor_id"].astype(str).isin(donors)
    sub = b[keep].to_memory()
    del b
    if "raw_sum" not in cols:
        X = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
        rs = np.asarray(X.sum(axis=1)).ravel()
        sub = sub[rs >= 1500].copy()
        log(f"  [{tag} | {cell_type}] raw_sum computed from X; "
            f"{int((rs >= 1500).sum())}/{len(rs)} cells kept")
    counts, meta, genes = pseudobulk_counts(sub)
    del sub
    log(f"  [{tag} | {cell_type}] {len(meta)} donors, {len(genes)} genes, "
        f"{int(counts.sum())} total counts")
    if len(meta) < MIN_DONORS:
        log(f"  [{tag} | {cell_type}] too few donors after pseudobulk; skipped")
        return None, None, None, None
    ok = save_npz_retry(cache, counts=counts, genes=np.array(genes),
                        donor=meta["donor"].to_numpy(),
                        age=meta["age"].to_numpy(float),
                        sex_male=meta["sex_male"].to_numpy(float),
                        n_cells=meta["n_cells"].to_numpy(np.int64),
                        mean_log_umi=meta["mean_log_umi"].to_numpy(float),
                        library=meta["library"].to_numpy(float))
    if not ok:
        log("  WARNING: counts cache save failed; continuing in-memory")
    return counts, meta, genes, cache


def run_unit_de(counts, metadata, gene_names, cohort, cell_type, scratch: Path, tag: str):
    """WLS + OLS DE for one unit; saves result .npz; returns summary dict."""
    wls = voom_wls_counts(counts, metadata, gene_names)
    ols = ols_arm(counts, metadata, gene_names)
    merged = wls.merge(ols, on="gene", how="left")
    merged["cohort"] = cohort
    merged["cell_type"] = cell_type
    merged["model"] = "rev2_voom_wls"
    safe = f"{tag}_{cell_type}".replace(" ", "_").replace("/", "_")
    out = scratch / f"de_{safe}.npz"
    save_npz_retry(out,
                   gene=np.array(merged["gene"]),
                   slope=merged["slope"].to_numpy(float),
                   se=merged["se"].to_numpy(float),
                   pvalue=merged["pvalue"].to_numpy(float),
                   q_bh=merged["q_bh"].to_numpy(float),
                   ci95_low=merged["ci95_low"].to_numpy(float),
                   ci95_high=merged["ci95_high"].to_numpy(float),
                   slope_ols=merged["slope_ols"].to_numpy(float),
                   se_ols=merged["se_ols"].to_numpy(float),
                   pvalue_ols=merged["pvalue_ols"].to_numpy(float),
                   q_bh_ols=merged["q_bh_ols"].to_numpy(float),
                   mean_cp10k=merged["mean_cp10k"].to_numpy(float),
                   detected_frac=merged["detected_frac"].to_numpy(float))
    n_fdr = int((merged["q_bh"] < 0.05).sum())
    n_fdr_ols = int((merged["q_bh_ols"] < 0.05).sum())
    rho = float(stats.spearmanr(merged["slope"], merged["slope_ols"]).statistic)
    summary = {
        "cohort": cohort, "cell_type": cell_type,
        "n_donors": int(len(metadata)),
        "n_genes_total": int(wls.attrs["n_genes_total"]),
        "n_genes_tested": int(wls.attrs["n_genes_tested"]),
        "n_genes_technical_dropped": int(wls.attrs["n_genes_technical_dropped"]),
        "design": str(wls.attrs["design"]),
        "n_fdr_wls": n_fdr,
        "n_fdr_ols_same_genes": n_fdr_ols,
        "spearman_slope_wls_vs_ols": rho,
        "median_n_cells": float(np.median(metadata["n_cells"])),
        "age_sd": float(metadata["age"].std()),
        "age_min": float(metadata["age"].min()),
        "age_max": float(metadata["age"].max()),
    }
    return summary, merged


def stage_lodo(scratch: Path):
    """Leave-one-donor-out sensitivity for the replicated (both-FDR) genes.

    Reads results/reviewer-rev2-interneuron-genes.csv (tier == both_FDR) and
    the cached counts .npz files; refits WLS dropping each donor once.
    """
    genes_df = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv")
    repl = genes_df[genes_df["tier"] == "both_FDR"]
    rows = []
    for (cohort, tag) in (("MSSM_DLPFC", "mssm"), ("HBCC_DLPFC", "hbcc")):
        for cell_type in INTERNEURON_TYPES:
            gl = repl[repl["cell_type"] == cell_type]["gene"].tolist()
            if not gl:
                continue
            counts, meta, gene_names = load_counts_npz(scratch, tag, cell_type)
            lodo = lodo_sensitivity(counts, meta, gene_names, gl)
            lodo["cohort"] = cohort
            lodo["cell_type"] = cell_type
            rows.append(lodo)
            log(f"[LODO {cohort} | {cell_type}] {len(lodo)} genes; "
                f"max |delta|/SE {lodo['max_abs_delta_over_se'].max():.2f}; "
                f"sign flips {(lodo['n_sign_flip'] > 0).sum()}; "
                f"lose-nominal {(lodo['n_lose_nominal'] > 0).sum()}")
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(df):
        save_csv_retry(df, scratch / "lodo_summary.csv")
        save_csv_retry(df, RESULTS / "reviewer-rev2-lodo-summary.csv")
    emit_json("lodo", {"genes": len(df),
                       "max_delta_over_se": float(df["max_abs_delta_over_se"].max()) if len(df) else None,
                       "n_sign_flip": int((df["n_sign_flip"] > 0).sum()) if len(df) else 0,
                       "n_lose_nominal": int((df["n_lose_nominal"] > 0).sum()) if len(df) else 0})
    return df


def stage_human_de(scratch: Path):
    """MSSM + HBCC interneuron DE (4 types x 2 cohorts)."""
    import cellxgene_census
    census = open_census_retry()
    summaries = []
    for cohort, dataset_id, tag in (("MSSM_DLPFC", MSSM_ID, "mssm"), ("HBCC_DLPFC", HBCC_ID, "hbcc")):
        for cell_type in INTERNEURON_TYPES:
            log(f"[{cohort} | {cell_type}] fetching")
            counts, meta, genes, _ = fetch_unit(census, dataset_id, cell_type, scratch, tag)
            if counts is None:
                continue
            log(f"[{cohort} | {cell_type}] fitting voom-WLS + OLS")
            summary, _ = run_unit_de(counts, meta, genes, cohort, cell_type, scratch, tag)
            summaries.append(summary)
            log(f"  -> FDR WLS {summary['n_fdr_wls']} / {summary['n_genes_tested']} tested; "
                f"OLS same genes {summary['n_fdr_ols_same_genes']}; "
                f"slope Spearman WLS-vs-OLS {summary['spearman_slope_wls_vs_ols']:.4f}")
            emit_json("human_de_unit", summary)
    df = pd.DataFrame(summaries)
    save_csv_retry(df, scratch / "human_de_summary.csv")
    save_csv_retry(df, RESULTS / "reviewer-rev2-human-de-summary.csv")
    emit_json("human_de", {"units": summaries})
    return df


def stage_human_de_h5ad(scratch: Path):
    """MSSM + HBCC interneuron DE from local Census h5ad files (same data
    as the SOMA path; see fetch_unit_h5ad). MSSM comes from the sibling
    Rev2 worker's download (~/workspace/rev2-mssm-null-scratch/mssm.h5ad);
    HBCC from this worker's download (scratch/hbcc.h5ad). A cohort whose
    file is not yet assembled is skipped with a plain log line, not a
    failure, so the stage can be re-run as downloads land.
    """
    h5ad_paths = {
        "MSSM_DLPFC": Path.home() / "workspace" / "rev2-mssm-null-scratch" / "mssm.h5ad",
        "HBCC_DLPFC": scratch / "hbcc.h5ad",
    }
    summaries = []
    for cohort, tag in (("MSSM_DLPFC", "mssm"), ("HBCC_DLPFC", "hbcc")):
        path = h5ad_paths[cohort]
        if not path.exists():
            log(f"[{cohort}] h5ad not present yet ({path}); skipped")
            continue
        for cell_type in INTERNEURON_TYPES:
            log(f"[{cohort} | {cell_type}] fetching from h5ad")
            counts, meta, genes, _ = fetch_unit_h5ad(path, cell_type, scratch, tag)
            if counts is None:
                continue
            log(f"[{cohort} | {cell_type}] fitting voom-WLS + OLS")
            summary, _ = run_unit_de(counts, meta, genes, cohort, cell_type, scratch, tag)
            summaries.append(summary)
            log(f"  -> FDR WLS {summary['n_fdr_wls']} / {summary['n_genes_tested']} tested; "
                f"OLS same genes {summary['n_fdr_ols_same_genes']}; "
                f"slope Spearman WLS-vs-OLS {summary['spearman_slope_wls_vs_ols']:.4f}")
            emit_json("human_de_unit", summary)
    df = pd.DataFrame(summaries)
    save_csv_retry(df, scratch / "human_de_summary.csv")
    save_csv_retry(df, RESULTS / "reviewer-rev2-human-de-summary.csv")
    emit_json("human_de", {"units": summaries})
    return df


def enumerate_cell_types(census, dataset_id):
    """All cell types in a Census dataset with donor counts (obs only)."""
    import cellxgene_census
    obs = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=f"dataset_id == '{dataset_id}'",
        column_names=["donor_id", "development_stage", "cell_type", "raw_sum",
                      "sex", "disease", "is_primary_data"],
    )
    obs["age"] = obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
    obs = obs[(obs["is_primary_data"] == True) & (obs["disease"] == "normal")
              & obs["age"].between(20, 120) & (obs["raw_sum"] >= 1500)].copy()
    counts = (obs.groupby(["cell_type", "donor_id"], observed=True).size()
              .rename("n_cells").reset_index())
    per_type = (counts[counts["n_cells"] >= MIN_CELLS_PER_DONOR_CT]
                .groupby("cell_type", observed=True)
                .agg(n_donors=("donor_id", "nunique"),
                     median_cells=("n_cells", "median"),
                     total_cells=("n_cells", "sum"))
                .reset_index())
    ages = obs.groupby("cell_type", observed=True)["age"].agg(["std", "min", "max"]).reset_index()
    per_type = per_type.merge(ages, on="cell_type", how="left")
    return per_type.sort_values("n_donors", ascending=False), obs


def enumerate_cell_types_h5ad(h5ad_path: Path):
    """All cell types in a local Census h5ad with donor counts (obs only).

    Same filters as enumerate_cell_types (primary, disease normal, age
    20-120 y, cell raw_sum >= 1500). raw_sum, if absent from obs, is
    computed from X in chunks (as in stage_mouse).
    """
    import anndata as ad
    b = ad.read_h5ad(h5ad_path, backed="r")
    obs = b.obs
    cols = set(obs.columns)
    mask = pd.Series(True, index=obs.index)
    if "is_primary_data" in cols:
        mask &= obs["is_primary_data"].astype(str).str.lower() == "true"
    if "disease" in cols:
        mask &= obs["disease"].astype(str) == "normal"
    age = obs["development_stage"].astype(str).str.extract(AGE_RE).iloc[:, 0].astype(float)
    mask &= age.between(20, 120).to_numpy()
    if "raw_sum" in cols:
        mask &= pd.to_numeric(obs["raw_sum"], errors="coerce") >= 1500
        obs = obs[mask].copy()
        obs["age"] = age[mask].to_numpy()
    else:
        obs = obs[mask].copy()
        obs["age"] = age[mask].to_numpy()
        # raw_sum not in obs: compute per-cell row sums from backed X.
        idx = np.where(mask.to_numpy())[0]
        rs = np.zeros(len(idx), dtype=np.float64)
        for i in range(0, len(idx), 10000):
            chunk = idx[i:i + 10000]
            Xi = b.X[chunk]
            Xi = Xi.tocsr() if sp.issparse(Xi) else sp.csr_matrix(Xi)
            rs[i:i + len(chunk)] = np.asarray(Xi.sum(axis=1)).ravel()
        obs = obs[rs >= 1500].copy()
    counts = (obs.groupby(["cell_type", "donor_id"], observed=True).size()
              .rename("n_cells").reset_index())
    per_type = (counts[counts["n_cells"] >= MIN_CELLS_PER_DONOR_CT]
                .groupby("cell_type", observed=True)
                .agg(n_donors=("donor_id", "nunique"),
                     median_cells=("n_cells", "median"),
                     total_cells=("n_cells", "sum"))
                .reset_index())
    ages = obs.groupby("cell_type", observed=True)["age"].agg(["std", "min", "max"]).reset_index()
    per_type = per_type.merge(ages, on="cell_type", how="left")
    return per_type.sort_values("n_donors", ascending=False), obs


def stage_all_celltypes_h5ad(scratch: Path, n_perms: int = 5):
    """Run the interneuron model in EVERY MSSM/HBCC cell type, from local
    h5ad files (same data as the SOMA path; see fetch_unit_h5ad). A cohort
    whose file is not yet assembled is skipped, not failed.
    """
    h5ad_paths = {
        "MSSM_DLPFC": Path.home() / "workspace" / "rev2-mssm-null-scratch" / "mssm.h5ad",
        "HBCC_DLPFC": scratch / "hbcc.h5ad",
    }
    rows = []
    for cohort, tag in (("MSSM_DLPFC", "mssm"), ("HBCC_DLPFC", "hbcc")):
        path = h5ad_paths[cohort]
        if not path.exists():
            log(f"[{cohort}] h5ad not present yet ({path}); skipped")
            continue
        per_type, _ = enumerate_cell_types_h5ad(path)
        log(f"[{cohort}] {len(per_type)} cell types with >= {MIN_CELLS_PER_DONOR_CT} cells "
            f"in >= 1 donor; {int((per_type['n_donors'] >= MIN_DONORS).sum())} with "
            f">= {MIN_DONORS} donors")
        for _, r in per_type[per_type["n_donors"] >= MIN_DONORS].iterrows():
            cell_type = str(r["cell_type"])
            counts, meta, genes, _ = fetch_unit_h5ad(path, cell_type, scratch, tag)
            if counts is None:
                continue
            wls = voom_wls_counts(counts, meta, genes)
            n_fdr = int((wls["q_bh"] < 0.05).sum())
            rng = np.random.default_rng(SEED)
            null_counts = []
            for _ in range(n_perms):
                meta_p = meta.copy()
                meta_p["age"] = rng.permutation(meta["age"].to_numpy())
                wls_p = voom_wls_counts(counts, meta_p, genes)
                null_counts.append(int((wls_p["q_bh"] < 0.05).sum()))
            rows.append({
                "cohort": cohort, "cell_type": cell_type,
                "short_type": SHORT.get(cell_type, cell_type),
                "n_donors": int(len(meta)),
                "n_genes_tested": int(wls.attrs["n_genes_tested"]),
                "median_n_cells": float(np.median(meta["n_cells"])),
                "age_sd": float(meta["age"].std()),
                "n_fdr_wls": n_fdr,
                "deg_rate_per_1000": 1000.0 * n_fdr / max(int(wls.attrs["n_genes_tested"]), 1),
                "median_se": float(wls["se"].median()),
                "null_fdr_mean": float(np.mean(null_counts)),
                "null_fdr_values": ";".join(str(int(x)) for x in null_counts),
                "obs_over_expected": float(n_fdr / max(np.mean(null_counts), 1e-9)),
                "is_interneuron": cell_type in INTERNEURON_TYPES,
            })
            log(f"  [{cohort} | {cell_type}] donors {len(meta)} tested "
                f"{int(wls.attrs['n_genes_tested'])} FDR {n_fdr} null-mean "
                f"{np.mean(null_counts):.1f} obs/exp {n_fdr / max(np.mean(null_counts), 1e-9):.2f}")
            emit_json("all_celltypes_unit", rows[-1])
    df = pd.DataFrame(rows)
    save_csv_retry(df, scratch / "all_celltypes.csv")
    save_csv_retry(df, RESULTS / "reviewer-rev2-interneuron-all-celltypes.csv")
    emit_json("all_celltypes", {"rows": len(df)})
    return df


def stage_all_celltypes(scratch: Path, n_perms: int = 5):
    """Run the interneuron model in EVERY MSSM/HBCC cell type.

    Power-adjusted comparison: per cell type, DEG count, DEG rate, median
    age-slope SE (precision), and a permutation null (age labels shuffled
    across donors, same WLS refit) giving expected DEG count under no age
    effect. 'Targeted' is judged on observed/expected DEG ratio and DEG
    rate, not on raw DEG counts.
    """
    import cellxgene_census
    census = open_census_retry()
    rows = []
    for cohort, dataset_id, tag in (("MSSM_DLPFC", MSSM_ID, "mssm"), ("HBCC_DLPFC", HBCC_ID, "hbcc")):
        per_type, _ = enumerate_cell_types(census, dataset_id)
        log(f"[{cohort}] {len(per_type)} cell types with >= {MIN_CELLS_PER_DONOR_CT} cells "
            f"in >= 1 donor; {int((per_type['n_donors'] >= MIN_DONORS).sum())} with "
            f">= {MIN_DONORS} donors")
        for _, r in per_type[per_type["n_donors"] >= MIN_DONORS].iterrows():
            cell_type = str(r["cell_type"])
            counts, meta, genes, _ = fetch_unit(census, dataset_id, cell_type, scratch, tag)
            if counts is None:
                continue
            wls = voom_wls_counts(counts, meta, genes)
            n_fdr = int((wls["q_bh"] < 0.05).sum())
            # permutation null: shuffle age labels, refit, count FDR hits
            rng = np.random.default_rng(SEED)
            null_counts = []
            for _ in range(n_perms):
                meta_p = meta.copy()
                meta_p["age"] = rng.permutation(meta["age"].to_numpy())
                wls_p = voom_wls_counts(counts, meta_p, genes)
                null_counts.append(int((wls_p["q_bh"] < 0.05).sum()))
            rows.append({
                "cohort": cohort, "cell_type": cell_type,
                "short_type": SHORT.get(cell_type, cell_type),
                "n_donors": int(len(meta)),
                "n_genes_tested": int(wls.attrs["n_genes_tested"]),
                "median_n_cells": float(np.median(meta["n_cells"])),
                "age_sd": float(meta["age"].std()),
                "n_fdr_wls": n_fdr,
                "deg_rate_per_1000": 1000.0 * n_fdr / max(int(wls.attrs["n_genes_tested"]), 1),
                "median_se": float(wls["se"].median()),
                "null_fdr_mean": float(np.mean(null_counts)),
                "null_fdr_values": ";".join(str(int(x)) for x in null_counts),
                "obs_over_expected": float(n_fdr / max(np.mean(null_counts), 1e-9)),
                "is_interneuron": cell_type in INTERNEURON_TYPES,
            })
            log(f"  [{cohort} | {cell_type}] donors {len(meta)} tested "
                f"{int(wls.attrs['n_genes_tested'])} FDR {n_fdr} null-mean "
                f"{np.mean(null_counts):.1f} obs/exp {n_fdr / max(np.mean(null_counts), 1e-9):.2f}")
            emit_json("all_celltypes_unit", rows[-1])
    df = pd.DataFrame(rows)
    save_csv_retry(df, scratch / "all_celltypes.csv")
    save_csv_retry(df, RESULTS / "reviewer-rev2-interneuron-all-celltypes.csv")
    emit_json("all_celltypes", {"rows": len(df)})
    return df

# ---------------------------------------------------------------------------
# Replication statistics (genome-wide first, selected second)
# ---------------------------------------------------------------------------

def load_de_npz(scratch: Path, tag: str, cell_type: str):
    safe = f"{tag}_{cell_type}".replace(" ", "_").replace("/", "_")
    z = np.load(scratch / f"de_{safe}.npz", allow_pickle=False)
    return pd.DataFrame({k: z[k] for k in
                         ("gene", "slope", "se", "pvalue", "q_bh", "ci95_low", "ci95_high",
                          "slope_ols", "se_ols", "pvalue_ols", "q_bh_ols",
                          "mean_cp10k", "detected_frac")})


def load_counts_npz(scratch: Path, tag: str, cell_type: str):
    safe = f"{tag}_{cell_type}".replace(" ", "_").replace("/", "_")
    z = np.load(scratch / f"counts_{safe}.npz", allow_pickle=False)
    meta = pd.DataFrame({k: z[k] for k in ("donor", "age", "sex_male", "n_cells", "mean_log_umi", "library")})
    return z["counts"], meta, np.array([str(g) for g in z["genes"]])


def coexpression_blocks(counts, gene_names, n_blocks=30, seed=SEED):
    """Coexpression blocks from top-PC loadings of logCPM pseudobulk.

    Genes are assigned to the PC (of the top 30 PCs of the donor x gene
    logCPM matrix) on which they load most strongly in absolute value.
    Genes sharing a dominant PC are coexpressed by construction; permuting
    gene labels within blocks preserves that coarse coexpression structure.
    """
    lib = counts.sum(axis=1, keepdims=True).astype(float)
    Y = np.log1p(counts / lib * 1e4)  # donors x genes
    Yc = Y - Y.mean(axis=0, keepdims=True)
    # SVD of donors x genes via covariance on donors (D << G)
    U, S, _ = np.linalg.svd(Yc, full_matrices=False)
    n_pc = min(n_blocks, U.shape[1])
    # gene loadings: V = Yc' U / S  (genes x PCs)
    loadings = (Yc.T @ U[:, :n_pc]) / np.maximum(S[:n_pc], 1e-12)[None, :]
    block = np.argmax(np.abs(loadings), axis=1)
    return block


def stage_replication(scratch: Path, n_block_perms: int = 1000):
    """Genome-wide replication stats + matched background + block null."""
    summaries, gene_rows, fig_rows = [], [], []
    for cell_type in INTERNEURON_TYPES:
        m = load_de_npz(scratch, "mssm", cell_type)
        h = load_de_npz(scratch, "hbcc", cell_type)
        merged = m.merge(h, on="gene", suffixes=("_mssm", "_hbcc"))
        merged["same_direction"] = merged["slope_mssm"] * merged["slope_hbcc"] > 0
        merged["mssm_fdr"] = merged["q_bh_mssm"] < 0.05
        merged["hbcc_fdr"] = merged["q_bh_hbcc"] < 0.05
        merged["both_fdr"] = merged["mssm_fdr"] & merged["hbcc_fdr"] & merged["same_direction"]
        merged["cell_type"] = cell_type
        merged["short_type"] = SHORT.get(cell_type, cell_type)
        # counts for the gene list (each cohort keeps its own gene order)
        cm, mm, gm = load_counts_npz(scratch, "mssm", cell_type)
        ch, mh, gh = load_counts_npz(scratch, "hbcc", cell_type)
        # -- genome-wide (PRIMARY): all tested genes in both cohorts --
        gw = merged
        pearson = float(stats.pearsonr(gw["slope_mssm"], gw["slope_hbcc"]).statistic)
        spearman = float(stats.spearmanr(gw["slope_mssm"], gw["slope_hbcc"]).statistic)
        same_dir = int(gw["same_direction"].sum())
        binom_p = float(stats.binomtest(same_dir, len(gw), 0.5, alternative="greater").pvalue)
        # -- selected-gene (SECONDARY): MSSM FDR only --
        sel = gw[gw["mssm_fdr"]]
        sel_pearson = float(stats.pearsonr(sel["slope_mssm"], sel["slope_hbcc"]).statistic) if len(sel) > 2 else float("nan")
        sel_spearman = float(stats.spearmanr(sel["slope_mssm"], sel["slope_hbcc"]).statistic) if len(sel) > 2 else float("nan")
        sel_same = int(sel["same_direction"].sum())
        sel_nominal = int((sel["same_direction"] & (sel["pvalue_hbcc"] < 0.05)).sum())
        both = int(gw["both_fdr"].sum())
        n_mssm_fdr = int(gw["mssm_fdr"].sum())
        n_hbcc_fdr = int(gw["hbcc_fdr"].sum())
        # hypergeometric (reference only; genes are NOT independent)
        from scipy.stats import hypergeom
        hyper_p = float(hypergeom.sf(both - 1, len(gw), n_mssm_fdr, n_hbcc_fdr))
        hyper_exp = n_mssm_fdr * n_hbcc_fdr / len(gw)
        # detected-gene-matched background: MSSM mean-CP10k decile matched
        # HBCC FDR rate -> expected both-FDR count among MSSM FDR genes.
        dec = pd.qcut(gw["mean_cp10k_mssm"], 10, labels=False, duplicates="drop")
        gw = gw.assign(decile=dec)
        exp_matched = 0.0
        for _, grp in gw.groupby("decile", observed=True):
            rate = float(grp["hbcc_fdr"].mean())
            exp_matched += float(grp["mssm_fdr"].sum()) * rate
        # coexpression-block permutation null for the both-FDR overlap.
        # Blocks are computed on the merged (tested) genes in merged order.
        block = None
        block_p, block_null_mean = float("nan"), float("nan")
        try:
            gidx_m = {g: i for i, g in enumerate(gm)}
            midx = np.array([gidx_m[g] for g in merged["gene"]])
            block = coexpression_blocks(cm[:, midx], merged["gene"].to_numpy())
            rng = np.random.default_rng(SEED)
            mssm_fdr_idx = np.flatnonzero(gw["mssm_fdr"].to_numpy())
            hbcc_fdr = gw["hbcc_fdr"].to_numpy()
            hbcc_sign = np.sign(gw["slope_hbcc"].to_numpy())
            mssm_sign = np.sign(gw["slope_mssm"].to_numpy())
            null_overlaps = []
            blocks = np.unique(block)
            for _ in range(n_block_perms):
                perm_hbcc_fdr = hbcc_fdr.copy()
                perm_hbcc_sign = hbcc_sign.copy()
                for b in blocks:
                    idx = np.flatnonzero(block == b)
                    perm = rng.permutation(idx)
                    perm_hbcc_fdr[idx] = hbcc_fdr[perm]
                    perm_hbcc_sign[idx] = hbcc_sign[perm]
                ov = int(np.sum(perm_hbcc_fdr[mssm_fdr_idx]
                                & (perm_hbcc_sign[mssm_fdr_idx] == mssm_sign[mssm_fdr_idx])))
                null_overlaps.append(ov)
            null_overlaps = np.array(null_overlaps)
            block_null_mean = float(null_overlaps.mean())
            block_p = float((1 + int(np.sum(null_overlaps >= both))) / (n_block_perms + 1))
        except Exception as exc:  # noqa: BLE001
            log(f"  [{cell_type}] block permutation failed: {exc}")
        summary = {
            "cell_type": cell_type, "short_type": SHORT.get(cell_type, cell_type),
            "n_tested_both": int(len(gw)),
            "gw_pearson_r": pearson, "gw_spearman_rho": spearman,
            "gw_same_direction": same_dir,
            "gw_same_direction_frac": same_dir / max(len(gw), 1),
            "gw_same_direction_binomial_p": binom_p,
            "mssm_fdr_genes": n_mssm_fdr,
            "hbcc_fdr_genes": n_hbcc_fdr,
            "sel_pearson_r": sel_pearson, "sel_spearman_rho": sel_spearman,
            "sel_same_direction": sel_same,
            "sel_same_direction_frac": sel_same / max(len(sel), 1),
            "sel_hbcc_nominal_same_direction": sel_nominal,
            "both_fdr_same_direction": both,
            "hypergeometric_expected": hyper_exp,
            "hypergeometric_p": hyper_p,
            "matched_background_expected": exp_matched,
            "block_perm_null_mean": block_null_mean,
            "block_perm_p": block_p,
            "block_perm_n": n_block_perms,
        }
        summaries.append(summary)
        emit_json("replication_unit", summary)
        # gene list rows (both-FDR primary + MSSM-FDR/HBCC-nominal secondary)
        cm_sum = cm.sum(axis=0)
        ch_sum = ch.sum(axis=0)
        gidx_m = {g: i for i, g in enumerate(gm)}
        gidx_h = {g: i for i, g in enumerate(gh)}
        keep_rows = merged[merged["mssm_fdr"]].copy()
        for _, r in keep_rows.iterrows():
            gi_m = gidx_m[r["gene"]]
            gi_h = gidx_h[r["gene"]]
            gene_rows.append({
                "cell_type": cell_type, "short_type": SHORT.get(cell_type, cell_type),
                "gene": r["gene"],
                "tier": "both_FDR" if bool(r["both_fdr"]) else "mssm_FDR_hbcc_nominal" if (r["same_direction"] and r["pvalue_hbcc"] < 0.05) else "mssm_FDR_only",
                "slope_mssm": r["slope_mssm"], "se_mssm": r["se_mssm"],
                "ci95_low_mssm": r["ci95_low_mssm"], "ci95_high_mssm": r["ci95_high_mssm"],
                "pvalue_mssm": r["pvalue_mssm"], "q_bh_mssm": r["q_bh_mssm"],
                "slope_hbcc": r["slope_hbcc"], "se_hbcc": r["se_hbcc"],
                "ci95_low_hbcc": r["ci95_low_hbcc"], "ci95_high_hbcc": r["ci95_high_hbcc"],
                "pvalue_hbcc": r["pvalue_hbcc"], "q_bh_hbcc": r["q_bh_hbcc"],
                "same_direction": bool(r["same_direction"]),
                "mssm_total_cells": int(cm_sum[gi_m]),
                "mssm_n_donors": int(len(mm)),
                "mssm_median_cells_per_donor": float(np.median(mm["n_cells"])),
                "hbcc_total_cells": int(ch_sum[gi_h]),
                "hbcc_n_donors": int(len(mh)),
                "hbcc_median_cells_per_donor": float(np.median(mh["n_cells"])),
                "mean_cp10k_mssm": r["mean_cp10k_mssm"],
                "detected_frac_mssm": r["detected_frac_mssm"],
            })
        # figure data rows (ALL tested genes, both cohorts)
        for _, r in merged.iterrows():
            fig_rows.append({
                "cell_type": cell_type, "short_type": SHORT.get(cell_type, cell_type),
                "gene": r["gene"],
                "slope_mssm": r["slope_mssm"], "se_mssm": r["se_mssm"],
                "ci95_low_mssm": r["ci95_low_mssm"], "ci95_high_mssm": r["ci95_high_mssm"],
                "pvalue_mssm": r["pvalue_mssm"], "q_bh_mssm": r["q_bh_mssm"],
                "slope_hbcc": r["slope_hbcc"], "se_hbcc": r["se_hbcc"],
                "ci95_low_hbcc": r["ci95_low_hbcc"], "ci95_high_hbcc": r["ci95_high_hbcc"],
                "pvalue_hbcc": r["pvalue_hbcc"], "q_bh_hbcc": r["q_bh_hbcc"],
                "mssm_fdr": bool(r["mssm_fdr"]), "hbcc_fdr": bool(r["hbcc_fdr"]),
                "both_fdr_same_direction": bool(r["both_fdr"]),
                "same_direction": bool(r["same_direction"]),
                "mean_cp10k_mssm": r["mean_cp10k_mssm"],
                "mean_cp10k_hbcc": r["mean_cp10k_hbcc"],
            })
    rep_df = pd.DataFrame(summaries)
    genes_df = pd.DataFrame(gene_rows)
    fig_df = pd.DataFrame(fig_rows)
    save_csv_retry(rep_df, scratch / "replication_summary.csv")
    save_csv_retry(rep_df, RESULTS / "reviewer-rev2-replication-summary.csv")
    save_csv_retry(genes_df, RESULTS / "reviewer-rev2-interneuron-genes.csv")
    save_csv_retry(fig_df, scratch / "figure_data.csv")
    emit_json("replication", {
        "units": summaries,
        "n_gene_rows": len(genes_df),
        "n_both_fdr": int((genes_df["tier"] == "both_FDR").sum()) if len(genes_df) else 0,
        "n_fig_rows": len(fig_df),
    })
    return rep_df, genes_df, fig_df

# ---------------------------------------------------------------------------
# Mouse stage (TMS droplet, Lung + Kidney, 3m vs 18m)
# ---------------------------------------------------------------------------

def stage_mouse(scratch: Path):
    """Re-derive the mouse signature DE from the raw TMS droplet object.

    QC mirrors M1 (per tissue: 500-5000 genes/cell, >= 1000 counts/cell).
    Scrublet doublet removal is NOT rerun: M1's per-cell doublet calls are
    not saved in the repo, and rerunning Scrublet would change the cell set
    in ways that cannot be validated against M1; the difference is reported
    as eligible-unit overlap vs M2. Eligibility mirrors M2: >= 3 mice and
    >= 50 cells per age group per tissue x cell type.
    Primary model: voom-WLS on mouse pseudobulk (logCPM ~ 1 + age18m + sex),
    same precision machinery as the human arm. Cross-checks: (a) statsmodels
    negative-binomial GLM on the P1 signature genes only; (b) audit of the
    saved PyDESeq2 tables (m2_de_*.csv) with threshold sweeps.
    """
    import anndata as ad
    raw = REPO / "data" / "raw" / "tms_droplet_raw.h5ad"
    if not raw.exists():
        emit_json("mouse", {"status": "blocked", "reason": "data/raw/tms_droplet_raw.h5ad missing"})
        log("mouse: raw TMS file missing; stage blocked")
        return None
    log("mouse: opening TMS droplet (backed)")
    b = ad.read_h5ad(raw, backed="r")
    log(f"mouse: raw object {b.n_obs} cells x {b.n_vars} genes")
    obs = b.obs[["tissue", "cell_ontology_class", "age", "sex", "mouse.id"]].copy()
    # Backed X is a _CSRDataset (no .sum); slice it in chunks. Do NOT use
    # obs n_genes: it is NaN for all 45,602 3m cells in this file, which
    # silently drops the entire 3m arm (that bug produced 0 eligible units
    # on the first run). Compute both per-cell quantities from X.
    n_counts = np.zeros(len(obs), dtype=np.float64)
    n_genes = np.zeros(len(obs), dtype=np.int64)
    for i in range(0, len(obs), 10000):
        Xi = b.X[i:i + 10000]
        Xi = Xi.tocsr() if sp.issparse(Xi) else sp.csr_matrix(Xi)
        n_counts[i:i + 10000] = np.asarray(Xi.sum(axis=1)).ravel()
        n_genes[i:i + 10000] = np.asarray((Xi > 0).sum(axis=1)).ravel()
    obs["n_counts"] = n_counts
    obs["n_genes"] = n_genes
    keep_cells = (obs["age"].isin(["3m", "18m"])
                  & obs["tissue"].isin(["Lung", "Kidney"])
                  & obs["cell_ontology_class"].notna()
                  & obs["n_genes"].between(500, 5000)
                  & (obs["n_counts"] >= 1000))
    log(f"mouse: cells passing age/tissue/QC filters: {int(keep_cells.sum())} / {len(obs)}")
    sub = b[keep_cells].to_memory()
    del b
    sub.obs.reset_index(drop=True, inplace=True)  # labels == positions from here on
    sub.obs["age_group"] = sub.obs["age"].astype(str)
    sub.obs["mouse"] = sub.obs["mouse.id"].astype(str)
    sub.obs["sex"] = sub.obs["sex"].astype(str)
    sub.obs["tissue"] = sub.obs["tissue"].astype(str)
    sub.obs["ctype"] = sub.obs["cell_ontology_class"].astype(str)
    genes = np.array([str(g) for g in sub.var_names])
    X = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    # eligibility per unit
    units = []
    for (tissue, ctype), idx in sub.obs.groupby(["tissue", "ctype"], observed=True).groups.items():
        d = sub.obs.loc[idx]
        m3 = d[d["age_group"] == "3m"]
        m18 = d[d["age_group"] == "18m"]
        ok = (m3["mouse"].nunique() >= 3 and m18["mouse"].nunique() >= 3
              and len(m3) >= 50 and len(m18) >= 50)
        units.append({"unit": f"{tissue}::{ctype}", "tissue": tissue, "ctype": ctype,
                      "n_cells_3m": len(m3), "n_cells_18m": len(m18),
                      "n_mice_3m": int(m3["mouse"].nunique()),
                      "n_mice_18m": int(m18["mouse"].nunique()), "eligible": bool(ok),
                      "idx": np.array(sorted(idx))})
    elig = [u for u in units if u["eligible"]]
    log(f"mouse: {len(units)} units evaluated, {len(elig)} eligible")
    # compare eligibility vs M2
    m2elig = pd.read_csv(RESULTS / "m2-unit-eligibility.csv")
    m2set = set(m2elig[m2elig["eligible"]]["unit"])
    myset = {u["unit"] for u in elig}
    log(f"mouse: eligible-unit overlap vs M2: {len(myset & m2set)} shared, "
        f"{len(myset - m2set)} new, {len(m2set - myset)} lost")
    emit_json("mouse_eligibility", {
        "n_units": len(units), "n_eligible": len(elig),
        "overlap_with_m2": len(myset & m2set),
        "new_units": sorted(myset - m2set), "lost_units": sorted(m2set - myset),
    })
    # per-unit DE
    import statsmodels.api as sm
    rows = []
    for u in elig:
        idx = u["idx"]
        o = sub.obs.iloc[idx]
        Xu = X[idx]
        # pseudobulk per mouse
        mice = sorted(o["mouse"].unique())
        crow, mrow = [], []
        for mouse in mice:
            mm = (o["mouse"] == mouse).to_numpy()
            crow.append(np.asarray(Xu[mm].sum(axis=0)).ravel())
            mrow.append({
                "mouse": mouse,
                "age18m": float(o.loc[mm, "age_group"].iloc[0] == "18m"),
                "sex_male": float(o.loc[mm, "sex"].astype(str).str.lower().iloc[0] == "male"),
                "n_cells": int(mm.sum()),
                "mean_log_umi": float(np.log(np.asarray(Xu[mm].sum(axis=1)).ravel()).mean()),
            })
        counts = np.vstack(crow).astype(np.int64)
        meta = pd.DataFrame(mrow)
        # gene filter: same rule as human (total>=10, detected>=10%, mean CP10k>=1, non-technical)
        keep, gs = filter_genes(counts, genes)
        C = counts[:, keep].astype(float)
        gk = genes[keep]
        lib = counts.sum(axis=1).astype(float)  # total library (all genes)
        Y = np.log1p(C / lib[:, None] * 1e4).T
        Xd, cols = design_matrix(meta.rename(columns={"age18m": "age"}),
                                 ["age", "sex_male", "mean_log_umi"])
        # NOTE: design_matrix looks for column 'age'; age18m renamed above.
        age_idx = cols.index("age")
        D, k = Xd.shape
        n_cells = meta["n_cells"].to_numpy(float)
        cell_factor = n_cells / np.median(n_cells)
        mean_logcpm = Y.mean(axis=1)
        beta_ols, *_ = np.linalg.lstsq(Xd, Y.T, rcond=None)
        resid = Y.T - Xd @ beta_ols
        s2 = np.sum(resid ** 2, axis=0) / max(D - k, 1)
        v_g, _, _ = _fit_trend(mean_logcpm, s2)
        W = cell_factor[None, :] / v_g[:, None]
        med = np.median(W, axis=1, keepdims=True)
        W = np.clip(W, 0.25 * med, 4.0 * med)
        beta, se, _, dof = _wls_batched(Y, Xd, W)
        slope = beta[:, age_idx]
        se_age = se[:, age_idx]
        tstat = np.divide(slope, se_age, out=np.zeros_like(slope), where=se_age > 0)
        pvalue = 2 * stats.t.sf(np.abs(tstat), dof)
        qvalue = multipletests(pvalue, method="fdr_bh")[1]
        # LODO-mouse for P1 genes in this unit
        hier = pd.read_csv(RESULTS / "m5-program-hierarchy.csv")
        p1 = set(hier[hier["tier"] == "P1_Clean"]["gene"])
        p1_here = [g for g in gk if g in p1]
        lodo = lodo_sensitivity(counts, meta.rename(columns={"age18m": "age"}), genes,
                                p1_here, columns=("age", "sex_male", "mean_log_umi")) \
            if p1_here else pd.DataFrame()
        # NB GLM cross-check on P1 genes (statsmodels, log link, offset log lib)
        nb_rows = []
        lib_all = counts.sum(axis=1).astype(float)
        nb_cols = ["age18m"] + (["sex_male"] if meta["sex_male"].nunique() > 1 else [])
        for g in p1_here:
            gi = int(np.flatnonzero(gk == g)[0])
            y = counts[:, keep][:, gi].astype(float)
            try:
                Xn = sm.add_constant(meta[nb_cols].astype(float), has_constant="add")
                # NB GLM with fixed dispersion alpha=1 (cross-check only;
                # dispersion is not estimated per gene here).
                fit = sm.GLM(y, Xn, family=sm.families.NegativeBinomial(),
                             offset=np.log(lib_all)).fit(maxiter=100)
                nb_rows.append({"gene": g, "nb_logfc": float(fit.params["age18m"]),
                                "nb_se": float(fit.bse["age18m"]),
                                "nb_p": float(fit.pvalues["age18m"])})
            except Exception as exc:  # noqa: BLE001
                nb_rows.append({"gene": g, "nb_logfc": float("nan"),
                                "nb_se": float("nan"), "nb_p": float("nan")})
        nb = pd.DataFrame(nb_rows)
        wls_p1 = pd.DataFrame({"gene": gk, "wls_log2fc": slope / np.log(2),
                               "wls_q": qvalue})
        cmp = nb.merge(wls_p1, on="gene", how="left")
        agree = float(np.mean(np.sign(cmp["nb_logfc"]) == np.sign(cmp["wls_log2fc"]))) if len(cmp) else float("nan")
        n_de = int(((qvalue < 0.05) & (np.abs(slope / np.log(2)) >= 0.25)).sum())
        rows.append({
            "unit": u["unit"], "tissue": u["tissue"], "cell_type": u["ctype"],
            "n_mice": D, "n_cells_3m": u["n_cells_3m"], "n_cells_18m": u["n_cells_18m"],
            "n_genes_tested": int(keep.sum()),
            "n_de_wls_fdr05_log2fc025": n_de,
            "n_p1_genes_in_unit": len(p1_here),
            "nb_wls_sign_agreement_p1": agree,
            "lodo_max_delta_over_se_p1": float(lodo["max_abs_delta_over_se"].max()) if len(lodo) else None,
            "lodo_p1_any_sign_flip": int((lodo["n_sign_flip"] > 0).sum()) if len(lodo) else 0,
        })
        # save per-unit WLS table
        unit_df = pd.DataFrame({"gene": gk, "log2fc_wls": slope / np.log(2),
                                "se_log2fc": se_age / np.log(2), "pvalue": pvalue,
                                "q_bh": qvalue})
        save_csv_retry(unit_df, scratch / f"mouse_de_{u['unit']}.csv".replace("::", "__").replace(" ", "_"))
        log(f"  [{u['unit']}] mice {D} tested {int(keep.sum())} DE {n_de} "
            f"P1-in-unit {len(p1_here)} NB/WLS sign agreement {agree:.2f}")
        emit_json("mouse_unit", rows[-1])
    df = pd.DataFrame(rows)
    save_csv_retry(df, scratch / "mouse_summary.csv")
    save_csv_retry(df, RESULTS / "reviewer-rev2-mouse-de-summary.csv")
    # audit saved PyDESeq2 tables: threshold sweeps for signature stability
    audit_rows = []
    hier = pd.read_csv(RESULTS / "m5-program-hierarchy.csv")
    p0 = set(hier[hier["tier"] == "P0_Raw257"]["gene"])
    p1 = set(hier[hier["tier"] == "P1_Clean"]["gene"])
    for fdr_t, fc_t, bm_t in ((0.05, 0.25, 0.0), (0.01, 0.5, 10.0), (0.05, 0.25, 10.0),
                             (0.05, 0.5, 0.0), (0.10, 0.25, 0.0)):
        got = set()
        for fp in sorted(RESULTS.glob("m2_de_*.csv")):
            d = pd.read_csv(fp)
            s = d[(d["FDR"] < fdr_t) & (d["log2FC"].abs() >= fc_t) & (d["baseMean"] >= bm_t)]
            got |= set(s["index"])
        audit_rows.append({
            "fdr_threshold": fdr_t, "log2fc_threshold": fc_t, "basemean_threshold": bm_t,
            "n_genes": len(got),
            "p0_overlap": len(got & p0), "p0_jaccard": len(got & p0) / max(len(got | p0), 1),
            "p1_overlap": len(got & p1), "p1_jaccard": len(got & p1) / max(len(got | p1), 1),
        })
    audit = pd.DataFrame(audit_rows)
    save_csv_retry(audit, RESULTS / "reviewer-rev2-mouse-deseq2-threshold-audit.csv")
    emit_json("mouse", {"units": rows, "deseq2_threshold_audit": audit_rows})
    return df


# ---------------------------------------------------------------------------
# PsychAD adjudication
# ---------------------------------------------------------------------------

def stage_psychad(scratch: Path):
    """Test the replicated interneuron genes in PsychAD as adjudication.

    PsychAD's published pseudobulk RDS (syn62064518) is not in the repo and
    needs Synapse authentication, so the first attempt is CELLxGENE Census:
    look for a PsychAD dataset, pull its IN_VIP / IN_PVALB / IN_SST types,
    and run the same voom-WLS model on the replicated genes (genome-wide
    slopes are computed too, for the correlation). If PsychAD is not in
    Census, report blocked with the exact reason.
    """
    import cellxgene_census
    census = open_census_retry()
    datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    hits = datasets[datasets["dataset_title"].astype(str).str.contains("PsychAD", case=False, na=False)]
    log(f"psychad: {len(hits)} PsychAD datasets in Census")
    emit_json("psychad_search", {"n_hits": int(len(hits)),
                                 "titles": hits["dataset_title"].astype(str).tolist() if len(hits) else []})
    if len(hits) == 0:
        emit_json("psychad", {"status": "blocked",
                              "reason": "No PsychAD dataset in CELLxGENE Census 2025-11-08; "
                                        "the PsychAD pseudobulk RDS (syn62064518) is not in the "
                                        "repo and requires Synapse authentication, which this "
                                        "worker does not have."})
        return None
    # If a PsychAD dataset exists, map its interneuron labels and run WLS.
    ds_id = str(hits.iloc[0]["dataset_id"])
    genes_df = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv")
    repl = genes_df[genes_df["tier"] == "both_FDR"]
    label_map = {"IN_VIP": "VIP", "IN_PVALB": "pvalb", "IN_SST": "sst",
                 "VIP GABAergic cortical interneuron": "VIP",
                 "pvalb GABAergic cortical interneuron": "pvalb",
                 "sst GABAergic cortical interneuron": "sst"}
    rows = []
    for label, short in label_map.items():
        try:
            counts, meta, genes, _ = fetch_unit(census, ds_id, label, scratch, "psychad")
        except Exception as exc:  # noqa: BLE001
            log(f"psychad: label {label} fetch failed: {exc}")
            continue
        if counts is None:
            continue
        wls = voom_wls_counts(counts, meta, genes)
        repl_sub = repl[repl["short_type"] == short][["gene", "slope_mssm"]]
        joined = wls[["gene", "slope"]].merge(repl_sub, on="gene", how="inner")
        sub = wls[wls["gene"].isin(set(repl_sub["gene"]))]
        for _, r in sub.iterrows():
            rows.append({"cell_type": label, "short_type": short, "gene": r["gene"],
                         "slope_psychad": r["slope"], "se_psychad": r["se"],
                         "ci95_low_psychad": r["ci95_low"], "ci95_high_psychad": r["ci95_high"],
                         "pvalue_psychad": r["pvalue"], "q_bh_psychad": r["q_bh"],
                         "n_donors": int(len(meta))})
        gw_r = float(stats.pearsonr(joined["slope"], joined["slope_mssm"]).statistic) \
            if len(joined) > 2 else float("nan")
        emit_json("psychad_unit", {"cell_type": label, "n_donors": int(len(meta)),
                                   "n_replicated_tested": int(len(sub)),
                                   "pearson_vs_mssm": gw_r})
    df = pd.DataFrame(rows)
    if len(df):
        save_csv_retry(df, RESULTS / "reviewer-rev2-psychad-adjudication.csv")
    emit_json("psychad", {"status": "done", "rows": len(df)})
    return df


# ---------------------------------------------------------------------------
# Emit (stdout chunking for large artifacts; the caller reassembles)
# ---------------------------------------------------------------------------

def stage_emit(scratch: Path, artifact: str, chunk: int, n_chunks: int, fmt: str):
    if artifact == "figure-data":
        df = pd.read_csv(scratch / "figure_data.csv")
        df = df.sort_values(["short_type", "gene"]).reset_index(drop=True)
        cols = ["short_type", "gene", "slope_mssm", "se_mssm", "q_bh_mssm",
                "slope_hbcc", "se_hbcc", "q_bh_hbcc", "both_fdr_same_direction",
                "mean_cp10k_mssm", "mean_cp10k_hbcc"]
        df = df[cols]
        for c in ("slope_mssm", "se_mssm", "slope_hbcc", "se_hbcc",
                  "mean_cp10k_mssm", "mean_cp10k_hbcc"):
            df[c] = df[c].round(6)
        for c in ("q_bh_mssm", "q_bh_hbcc"):
            df[c] = df[c].apply(lambda v: f"{v:.3g}")
        text = df.to_csv(index=False)
    elif artifact == "genes":
        text = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv").to_csv(index=False)
    elif artifact == "all-celltypes":
        text = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-all-celltypes.csv").to_csv(index=False)
    else:
        raise ValueError(f"unknown artifact {artifact}")
    if fmt == "b64gz":
        import base64, gzip
        blob = base64.b64encode(gzip.compress(text.encode(), compresslevel=9)).decode()
        lines = [blob[i:i + 100] for i in range(0, len(blob), 100)]
        size = (len(lines) + n_chunks - 1) // n_chunks
        total = n_chunks
        print(f"@@CHUNK {artifact} {chunk}/{total} fmt=b64gz nlines={len(lines)}")
        for line in lines[chunk * size:(chunk + 1) * size]:
            print(line)
        print(f"@@ENDCHUNK {artifact} {chunk}/{total}")
    else:
        lines = text.splitlines(keepends=True)
        size = (len(lines) + n_chunks - 1) // n_chunks
        total = n_chunks
        print(f"@@CHUNK {artifact} {chunk}/{total} fmt=raw nlines={len(lines)}")
        for line in lines[chunk * size:(chunk + 1) * size]:
            print(line, end="")
        print(f"@@ENDCHUNK {artifact} {chunk}/{total}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True,
                        choices=["probe", "human-de", "human-de-h5ad", "replication", "lodo", "all-celltypes",
                                 "all-celltypes-h5ad", "mouse", "psychad", "emit"])
    parser.add_argument("--scratch-dir", type=Path,
                        default=Path.home() / "workspace" / "rev2-de-scratch")
    parser.add_argument("--artifact", default="figure-data")
    parser.add_argument("--chunk", type=int, default=0)
    parser.add_argument("--n-chunks", type=int, default=1)
    parser.add_argument("--format", default="raw", choices=["raw", "b64gz"])
    args = parser.parse_args()
    scratch: Path = args.scratch_dir
    scratch.mkdir(parents=True, exist_ok=True)
    if args.stage == "probe":
        try:
            import cellxgene_census  # noqa: F401
            print("PROBE cellxgene_census import OK")
        except ImportError as exc:
            print(f"PROBE cellxgene_census import FAILED: {exc}")
        try:
            import anndata, scanpy  # noqa: F401
            print("PROBE anndata/scanpy import OK")
        except ImportError as exc:
            print(f"PROBE anndata/scanpy import FAILED: {exc}")
        print(f"PROBE TMS raw exists: {(REPO / 'data' / 'raw' / 'tms_droplet_raw.h5ad').exists()}")
        return
    if args.stage == "human-de":
        stage_human_de(scratch)
    elif args.stage == "human-de-h5ad":
        stage_human_de_h5ad(scratch)
    elif args.stage == "replication":
        stage_replication(scratch)
    elif args.stage == "lodo":
        stage_lodo(scratch)
    elif args.stage == "all-celltypes":
        stage_all_celltypes(scratch)
    elif args.stage == "all-celltypes-h5ad":
        stage_all_celltypes_h5ad(scratch)
    elif args.stage == "mouse":
        stage_mouse(scratch)
    elif args.stage == "psychad":
        stage_psychad(scratch)
    elif args.stage == "emit":
        stage_emit(scratch, args.artifact, args.chunk, args.n_chunks, args.format)


if __name__ == "__main__":
    main()
