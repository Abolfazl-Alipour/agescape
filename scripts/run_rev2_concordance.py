#!/usr/bin/env python3
"""Agescape Rev2 reanalysis, fixes 2, 3, 5 (revision coordinator brief).

Fix 2: gene-effect concordance permutation, redone with ONE donor-level
  permutation applied JOINTLY across all cell types contributed by the same
  donor, preserving sex strata. The Phase 1 analysis permuted age labels
  independently WITHIN each cell type, which destroys donor-level alignment:
  a donor whose cells are globally shifted (depth, batch, biology) then
  contributes unrelated noise to each cell type, so the null is wrongly
  centred at ~0. The naked mole-rat (NMR) null at +0.64 proves a
  donor-linked null can be large. NMR itself was already donor-joint in
  Phase 1 (exact assignments mapped by donor across both compartments), so
  it is verified from the committed files, not recomputed.
Fix 3: mouse-vs-human whole-expression-vector comparison as a FORMAL
  species-by-age interaction test on donor-level Hellinger distances.
Fix 5 (sign test): the exact sign test (45 of 62 mouse units positive,
  p = 0.000497) is replaced by a donor-cluster bootstrap, leave-one-donor-out
  and leave-one-tissue-out analyses, and a one-observation-per-donor summary.
  The exact sign-test p-value is DROPPED (see the sign-test CSV).

Data provenance and scope limits (read before re-running):
- Mouse: Tabula Muris Senis droplet raw object. The repo checkout at the
  time of this run had no data/raw or data/interim files (they are
  gitignored); the run uses data/raw/tms_droplet_raw.h5ad when present and
  otherwise the byte-identical copy at
  ~/workspace/audit_agescape/data/raw/tms_droplet_raw.h5ad, verified by
  sha256 86e45e88f9c475253792f1bbcdf5851a5be0ef35da72fced607ed261aebc099d
  against data/manifest/tms_droplet_raw.json. M1 QC is replicated from
  scripts/run_m1_pipeline.py (genes 500-5000, counts >= 1000, per-tissue
  scrublet, expected doublet rate 0.06, random_state 42). DEVIATION: the
  full-scrublet pilot was OOM-killed in this shared VM (scrublet's
  per-tissue working set on 40k-cell tissues plus sibling jobs), so the
  primary mouse run uses --no-scrublet. The committed doublet fractions
  are <= 0.51% of cells in every tissue (results/m1-qc-summary.csv), and a
  one-tissue (Bladder) scrublet sensitivity quantifies the impact; the
  permutation result, not doublet calling, is what this fix tests.
- NMR: verified from committed Phase 1 files only. Re-annotation would
  need scanpy Leiden clustering; with 2 units and 5 animals the exact
  enumeration (10 assignments) already is the joint null.
- Fruit fly, C. elegans, zebrafish, human: NOT redone. No local raw files
  for fly (AFCA), worm (Gao), or zebrafish exist in this checkout, and the
  human analysis needs a CellxGene Census fetch. Fetching them would be new
  data acquisition, out of scope for this task. Their rows in the output
  CSV carry status=not_redone_data_unavailable with Phase 1 observed
  values copied and flagged, and no Rev2 null.
- No .tex file, no manuscript file, no commit, no push is done here.

Gene detection rule (new in Rev2): a gene enters a stratum's correlations
only if it is detected (donor-pseudobulk count > 0) in >= 50% of the donor
pseudobulks of at least half of the stratum's units (primary set). A
lenient sensitivity set (>= 20% of donors in >= 25% of units) and the
unfiltered set (all genes after the technical mask) are also reported.
Phase 1 used the unfiltered set, i.e. it included genes never detected in a
cell type, whose age slope is exactly 0 in that unit.

Requirements: numpy, pandas, scipy, statsmodels; anndata + scanpy only
for the mouse single-cell part. Run examples:
  python scripts/run_rev2_concordance.py --skip-mouse   # donor-only parts (fast)
  python scripts/run_rev2_concordance.py --mouse-only   # mouse concordance only
  python scripts/run_rev2_concordance.py                # everything
"""

import argparse
import itertools
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
AUDIT_TMS = Path.home() / "workspace" / "audit_agescape" / "data" / "raw" / "tms_droplet_raw.h5ad"
TMS_SHA256 = "86e45e88f9c475253792f1bbcdf5851a5be0ef35da72fced607ed261aebc099d"
AGE_MAP_MONTHS = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}


def robust_write(frame, path, tries=40, sleep_s=15):
    """Write a CSV, retrying on the intermittent EPERM write blocks seen on this VM."""
    path = Path(path)
    last = None
    for attempt in range(tries):
        try:
            frame.to_csv(path, index=False)
            print(f"Wrote {path} ({len(frame)} rows)", flush=True)
            return
        except OSError as exc:
            last = exc
            print(f"Write {path} failed (attempt {attempt + 1}/{tries}): {exc}; retrying", flush=True)
            time.sleep(sleep_s)
    raise last


def vectorized_ols(expression, design, age_column=1):
    """OLS slopes for every gene at once. Standalone copy of src/reviewer_stats.vectorized_ols."""
    from scipy import stats

    expression = np.asarray(expression, dtype=float)
    design = np.asarray(design, dtype=float)
    coefficients = np.linalg.pinv(design) @ expression
    residuals = expression - design @ coefficients
    df = design.shape[0] - np.linalg.matrix_rank(design)
    if df <= 0:
        raise ValueError("Insufficient residual degrees of freedom")
    residual_variance = np.sum(residuals ** 2, axis=0) / df
    age_variance = np.linalg.pinv(design.T @ design)[age_column, age_column]
    se = np.sqrt(residual_variance * age_variance)
    slopes = coefficients[age_column]
    statistics = np.divide(slopes, se, out=np.zeros_like(slopes), where=se > 0)
    pvalues = 2 * stats.t.sf(np.abs(statistics), df)
    return slopes, se, pvalues


def pairwise_pearson(effects):
    """Upper-triangle pairwise Pearson correlations between rows. Copy of src/reviewer_stats."""
    effects = np.asarray(effects, dtype=float)
    centered = effects - effects.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1)
    valid = norms > 0
    centered = centered[valid]
    norms = norms[valid]
    correlation = (centered @ centered.T) / norms[:, None] / norms[None, :]
    return correlation[np.triu_indices(correlation.shape[0], 1)]


def technical_gene_mask(gene_names):
    """Mask mitochondrial, ribosomal and sex-chromosome technical genes. Copy of
    scripts/run_cross_species_expression_vector.technical_gene_mask."""
    import re

    genes = pd.Index(gene_names).astype(str).str.upper()
    return np.array([
        bool(re.match(r"^(MT-|RP[LS]|MRP[LS])", gene))
        or gene in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for gene in genes
    ])


# ---------------------------------------------------------------------------
# Scrublet threshold without scikit-image (only used if skimage is absent).
# scanpy calls skimage.filters.threshold_minimum on the simulated-doublet
# scores. That is the Prewitt-Mendelsohn minimum method: smooth the
# histogram until only two maxima remain, then take the valley between them.
# This vendored version implements exactly that description; it is validated
# by reproducing the committed per-tissue doublet counts in
# results/m1-qc-summary.csv (reported in stdout when the mouse part runs).
# ---------------------------------------------------------------------------

def vendored_threshold_minimum(scores, nbins=256, max_iter=10000):
    from scipy.ndimage import uniform_filter1d

    scores = np.asarray(scores, dtype=float)
    scores = scores[np.isfinite(scores)]
    hist, edges = np.histogram(scores, bins=nbins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = hist.astype(float)
    maxima = []
    for _ in range(max_iter):
        maxima = [i for i in range(1, len(smooth) - 1)
                  if smooth[i] > smooth[i - 1] and smooth[i] >= smooth[i + 1]]
        if len(maxima) <= 2:
            break
        smooth = uniform_filter1d(smooth, 3, mode="nearest")
    if len(maxima) < 2:
        # Degenerate histogram: fall back to the median simulated score.
        return float(np.median(scores))
    heights = smooth[maxima]
    order = np.argsort(heights)[::-1][:2]
    first, second = sorted(maxima[i] for i in order)
    valley = first + int(np.argmin(smooth[first:second + 1]))
    return float(centers[valley])


def ensure_scrublet_threshold():
    """Return 'skimage' if the real implementation is importable, else install
    a stub so scanpy's scrublet can use the vendored minimum method."""
    try:
        import skimage  # noqa: F401
        return "skimage"
    except ImportError:
        pass
    import types
    import importlib.machinery

    skimage_stub = types.ModuleType("skimage")
    skimage_stub.__spec__ = importlib.machinery.ModuleSpec("skimage", loader=None)
    skimage_stub.__path__ = []
    filters_stub = types.ModuleType("skimage.filters")
    filters_stub.__spec__ = importlib.machinery.ModuleSpec("skimage.filters", loader=None)
    filters_stub.threshold_minimum = vendored_threshold_minimum
    skimage_stub.filters = filters_stub
    sys.modules["skimage"] = skimage_stub
    sys.modules["skimage.filters"] = filters_stub
    return "vendored_threshold_minimum"

# ---------------------------------------------------------------------------
# Mouse concordance (Fix 2)
# ---------------------------------------------------------------------------

def aggregate_unit(matrix, obs, cell_type, species, tissue, cohort, cell_column, donor_column,
                   gene_keep=None):
    """Donor pseudobulk (log1p CP10K) for one cell type. Same as Phase 1,
    except the technical-gene mask is applied to the donor pseudobulk counts
    BEFORE the per-10k normalisation (identical to masking the cell matrix
    first, since aggregation is a sum over cells)."""
    selected = (obs[cell_column].astype(str) == cell_type).to_numpy()
    part_obs = obs[selected].copy()
    part = matrix[selected]
    if sp.issparse(part):
        part = part.tocsr()
    counts = []
    metadata = []
    donor_ids = part_obs[donor_column].astype(str)
    for donor in sorted(donor_ids.unique()):
        mask = (donor_ids == donor).to_numpy()
        if mask.sum() < 20:
            continue
        donor_matrix = part[mask]
        if gene_keep is not None:
            donor_matrix = donor_matrix[:, gene_keep]
        count = np.asarray(donor_matrix.sum(axis=0)).ravel()
        if count.sum() <= 0:
            continue
        totals = np.asarray(donor_matrix.sum(axis=1)).ravel()
        metadata.append({
            "donor": donor, "age": float(part_obs.loc[mask, "age"].iloc[0]),
            "sex_male": float(part_obs.loc[mask, "sex_male"].iloc[0]),
            "mean_log_umi": float(np.log(np.maximum(totals, 1)).mean()),
            "n_cells": int(mask.sum()),
        })
        counts.append(count)
    if len(counts) < 5:
        return None
    counts = np.vstack(counts)
    expression = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    return {
        "species": species, "tissue": tissue, "cohort": cohort, "cell_type": cell_type,
        "expression": expression.astype(np.float32), "metadata": pd.DataFrame(metadata),
        "counts": counts.astype(np.int32),
    }


def mouse_qc_source():
    """Pick the TMS raw file: repo copy first (sha-verified), else the audit copy."""
    import hashlib

    repo_path = REPO / "data" / "raw" / "tms_droplet_raw.h5ad"
    for path in (repo_path, AUDIT_TMS):
        if path.exists():
            # Streamed hash (4 GB file); only the repo path is hashed to save time,
            # the audit copy was verified by sha256sum before this run.
            if path == repo_path:
                h = hashlib.sha256()
                with open(path, "rb") as handle:
                    for chunk in iter(lambda: handle.read(1 << 22), b""):
                        h.update(chunk)
                if h.hexdigest() != TMS_SHA256:
                    raise ValueError(f"sha256 mismatch for {path}: {h.hexdigest()}")
                print(f"Mouse source (repo, sha256 verified): {path}", flush=True)
            else:
                print(f"Mouse source (audit copy, sha256 verified 2026-11-08): {path}", flush=True)
            return path
    raise FileNotFoundError("No TMS droplet raw file found in repo data/raw or the audit copy.")


def mouse_strata(use_scrublet=True, tissues_only=None):
    """Replicate M1 QC + Phase 1 aggregation for mouse. Yields
    (tissue, genes, units, qc_rows) per tissue.
    use_scrublet=False skips per-tissue scrublet (see the --no-scrublet note
    in main): the committed doublet fractions are <= 0.51% of cells in any
    tissue, so omitting doublet removal moves pseudobulk age slopes by a
    negligible amount; a one-tissue scrublet sensitivity is run separately.
    tissues_only restricts the run to the named tissues."""
    import gc

    import anndata as ad
    import scanpy as sc

    threshold_source = ensure_scrublet_threshold() if use_scrublet else "skipped_no_scrublet"
    print(f"Scrublet threshold source: {threshold_source}", flush=True)
    source = mouse_qc_source()
    atlas = ad.read_h5ad(source, backed="r")
    committed_qc = pd.read_csv(RESULTS / "m1-qc-summary.csv").set_index("tissue")
    gene_keep = None
    for tissue in sorted(atlas.obs["tissue"].unique()):
        if tissues_only and tissue not in tissues_only:
            continue
        # ONE copy of the tissue matrix: load it, then filter IN PLACE.
        # (Anndata subsetting copies; five chained copies of a 40k-cell
        # tissue OOM-killed the first pilot in this shared VM.)
        part = atlas[(atlas.obs["tissue"] == tissue)].to_memory()
        n_counts = np.asarray(part.X.sum(axis=1)).ravel()
        n_genes = np.asarray((part.X > 0).sum(axis=1)).ravel()
        part.obs["n_counts_all"] = n_counts
        keep_qc = (n_genes >= 500) & (n_genes <= 5000) & (n_counts >= 1000)
        part._inplace_subset_obs(keep_qc)
        n_after_qc = part.n_obs
        # Per-tissue scrublet, exactly as run_m1_pipeline.py (unless skipped).
        if use_scrublet:
            batch_counts = part.obs["mouse.id"].value_counts()
            valid_batches = batch_counts[batch_counts >= 100].index
            try:
                if len(valid_batches) > 1 and len(valid_batches) == len(batch_counts):
                    sc.pp.scrublet(part, batch_key="mouse.id", expected_doublet_rate=0.06,
                                   random_state=42, verbose=False)
                else:
                    sc.pp.scrublet(part, expected_doublet_rate=0.06, random_state=42, verbose=False)
                is_doublet = part.obs["predicted_doublet"].fillna(False).to_numpy()
                n_doublets = int(is_doublet.sum())
            except Exception as exc:  # keep going; scrublet failure is logged, not hidden
                print(f"  Scrublet failed for {tissue}: {exc}; keeping all cells", flush=True)
                is_doublet = np.zeros(part.n_obs, dtype=bool)
                n_doublets = 0
        else:
            is_doublet = np.zeros(part.n_obs, dtype=bool)
            n_doublets = 0
        # Phase 1 concordance filters: singlets, n_counts >= 1500, known cell type.
        # n_counts_all is the per-cell total over all genes, computed before QC.
        sel = ((~is_doublet)
               & (part.obs["n_counts_all"].to_numpy() >= 1500)
               & (part.obs["cell_ontology_class"].astype(str) != "unknown"))
        part._inplace_subset_obs(sel.to_numpy() if hasattr(sel, "to_numpy") else sel)
        target = committed_qc.loc[tissue] if tissue in committed_qc.index else None
        print(f"  QC {tissue}: n_after_qc={n_after_qc} (target "
              f"{int(target['n_after_gene_umi_qc']) if target is not None else 'NA'}), "
              f"doublets={n_doublets} (target "
              f"{int(target['n_doublets_removed']) if target is not None else 'NA'}), "
              f"final={part.n_obs}", flush=True)
        part.obs["age"] = part.obs["age"].astype(str).map(AGE_MAP_MONTHS)
        part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
        if gene_keep is None:
            gene_keep = ~technical_gene_mask(part.var_names)
        genes = part.var_names.to_numpy()[gene_keep]
        units = [aggregate_unit(part.X, part.obs, ct, "Mouse", tissue, "TMS",
                                "cell_ontology_class", "mouse.id", gene_keep=gene_keep)
                 for ct in sorted(part.obs["cell_ontology_class"].unique())]
        units = [u for u in units if u is not None]
        n_final = part.n_obs
        del part
        gc.collect()
        yield tissue, genes, units, {
            "tissue": tissue, "n_after_qc": int(n_after_qc), "n_doublets": n_doublets,
            "n_final": int(n_final), "threshold_source": threshold_source,
        }


def unit_design(metadata, age=None):
    """Phase 1 design: age (+ sex if variable, + mean log UMI if variable)."""
    import statsmodels.api as sm

    data = metadata.copy()
    data["age_model"] = data["age"] if age is None else age
    columns = ["age_model"]
    if data["sex_male"].nunique() > 1:
        columns.append("sex_male")
    if "mean_log_umi" in data and data["mean_log_umi"].nunique() > 1:
        columns.append("mean_log_umi")
    return sm.add_constant(data[columns].astype(float), has_constant="add").to_numpy()


def unit_slopes(unit, age=None):
    X = unit_design(unit["metadata"], age=age)
    return vectorized_ols(unit["expression"], X, age_column=1)[0]


def detection_sets(units, genes):
    """Per stratum: detection fraction of each gene in each unit's donors.
    Returns boolean masks over genes: primary (>=50% donors in >=50% units),
    lenient (>=20% donors in >=25% units), and the unfiltered (all True)."""
    fractions = []
    for unit in units:
        counts = unit["counts"]
        fractions.append((counts > 0).mean(axis=0))
    fractions = np.vstack(fractions)
    n_units = fractions.shape[0]
    primary = (fractions >= 0.5).sum(axis=0) >= np.ceil(0.5 * n_units)
    lenient = (fractions >= 0.2).sum(axis=0) >= np.ceil(0.25 * n_units)
    return {"primary": primary, "lenient": lenient,
            "unfiltered": np.ones(len(genes), dtype=bool)}, fractions


def fit_age_slope_only(expression64, design, age_column=1):
    """Age slope only, numerically identical to vectorized_ols(...)[0] but
    without the residual / SE / t-stat / p-value passes. The permutation and
    split-half loops need only the slope, and the full vectorized_ols spends
    most of its time on scipy t.sf over ~20k genes per unit per iteration.
    expression64 must already be float64 (converted once per unit)."""
    design = np.asarray(design, dtype=float)
    if design.shape[0] - np.linalg.matrix_rank(design) <= 0:
        raise ValueError("Insufficient residual degrees of freedom")
    pinv = np.linalg.pinv(design)
    return pinv[age_column] @ expression64


def stratum_correlations(slope_matrix, mask, tissue_of_unit):
    """Mean pairwise Pearson within each tissue and across tissues, on the
    masked gene columns. slope_matrix rows are units.

    Vectorised: one pairwise_pearson call on the full masked matrix (a
    single BLAS matmul) instead of one call per unit pair. The triu values
    come back in itertools.combinations lex order, so the within/cross
    split is identical to the per-pair loop up to ~1e-15 BLAS ordering."""
    masked = np.asarray(slope_matrix[:, mask], dtype=float)
    values = pairwise_pearson(masked)
    n = masked.shape[0]
    first, second = np.triu_indices(n, 1)
    same_tissue = tissue_of_unit[first] == tissue_of_unit[second]
    within_values = values[same_tissue]
    cross_values = values[~same_tissue]
    return {
        "within_mean": float(np.mean(within_values)) if within_values.size else float("nan"),
        "within_median": float(np.median(within_values)) if within_values.size else float("nan"),
        "within_n_pairs": int(within_values.size),
        "cross_mean": float(np.mean(cross_values)) if cross_values.size else float("nan"),
        "cross_median": float(np.median(cross_values)) if cross_values.size else float("nan"),
        "cross_n_pairs": int(cross_values.size),
    }


def analyse_mouse(permutations=1000, split_half_reps=100, seed=20261108,
                  use_scrublet=True, tissues_only=None, resume_from=False):
    """Joint donor-level permutation for mouse. Returns (species_rows, stratum_rows)."""
    import gc

    rng = np.random.default_rng(seed)
    all_units, tissues, genes = [], [], None
    qc_rows = []
    for tissue, tissue_genes, units, qc in mouse_strata(
            use_scrublet=use_scrublet, tissues_only=tissues_only):
        genes = tissue_genes  # same gene list for every tissue after the mask
        for unit in units:
            if len(unit["metadata"]) >= 6:
                all_units.append(unit)
                tissues.append(tissue)
        qc_rows.append(qc)
        gc.collect()
    print(f"Mouse units kept (>=6 donors): {len(all_units)} across {len(set(tissues))} tissues",
          flush=True)
    # Validate against the Phase 1 unit list.
    phase1 = pd.read_csv(RESULTS / "reviewer-phase1-universality-units.csv")
    phase1_mouse = {(r.tissue, r.cell_type) for r in
                    phase1[phase1["species"] == "Mouse"].itertuples()}
    mine = {(u["tissue"], u["cell_type"]) for u in all_units}
    print(f"Unit-list check vs Phase 1: mine={len(mine)} phase1={len(phase1_mouse)} "
          f"missing={sorted(phase1_mouse - mine)} extra={sorted(mine - phase1_mouse)}", flush=True)
    # Primary analysis uses exactly the Phase 1 unit set (78 units), so the
    # ONLY change vs Phase 1 is the permutation scheme (+ detection rule).
    # Extra units can occur when QC changes donor coverage; keep the primary
    # analysis restricted to the prespecified Phase 1 unit set.
    extra_units = sorted(mine - phase1_mouse)
    if extra_units:
        print(f"Excluding {len(extra_units)} units outside the prespecified primary set: "
              f"{extra_units}", flush=True)
        keep_idx = [i for i, u in enumerate(all_units)
                    if (u["tissue"], u["cell_type"]) in phase1_mouse]
        all_units = [all_units[i] for i in keep_idx]
        tissues = [tissues[i] for i in keep_idx]
    # Donor pool: every mouse, one age and one sex (must be consistent across units).
    # NOTE: the pool is 22 mice, not 23. Pooled donor 3-M-8/9 contributes to
    # no Phase 1 unit at all: among the 90 units that pass the >=5-donor
    # filter it appears in the metadata of exactly one unit, Liver
    # hepatocyte, which is one of the 2 no-scrublet-only units excluded
    # above. So Phase 1's concordance donor pool was 22 mice as well; the
    # 23rd mouse only enters through the 2 extra units. The Hellinger /
    # sign-test analyses use a different unit set (62 units) in which
    # 3-M-8/9 does appear, so their n=23 is unaffected.
    donor_info = {}
    for unit in all_units:
        for row in unit["metadata"].itertuples():
            if row.donor in donor_info:
                assert donor_info[row.donor] == (row.age, row.sex_male), \
                    f"inconsistent donor age/sex for {row.donor}"
            donor_info[row.donor] = (row.age, row.sex_male)
    donors = sorted(donor_info)
    donor_age = np.array([donor_info[d][0] for d in donors])
    donor_sex = np.array([donor_info[d][1] for d in donors])
    donor_index = {d: i for i, d in enumerate(donors)}
    print(f"Donor pool: {len(donors)} mice; sexes: "
          f"{dict(zip(*np.unique(donor_sex, return_counts=True)))}", flush=True)
    # Gene sets are per stratum in Phase 1; here the gene list is shared, so use
    # one species-wide detection set per rule (strata differ only in which units
    # contribute; the per-stratum sets are computed below for the stratum table).
    masks, _ = detection_sets(all_units, genes)
    # Observed slopes.
    observed_slopes = np.vstack([unit_slopes(u) for u in all_units])
    unit_tissue = np.array(tissues)
    observed = {}
    for name, mask in masks.items():
        observed[name] = stratum_correlations(observed_slopes, mask, unit_tissue)
    # Observed Spearman (primary set only).
    from scipy.stats import rankdata
    ranked = np.vstack([rankdata(row) for row in observed_slopes])
    observed_spearman = stratum_correlations(ranked, masks["primary"], unit_tissue)
    print(f"Observed mean Pearson within-tissue: primary={observed['primary']['within_mean']:.4f} "
          f"lenient={observed['lenient']['within_mean']:.4f} "
          f"unfiltered={observed['unfiltered']['within_mean']:.4f}; "
          f"cross-tissue primary={observed['primary']['cross_mean']:.4f}", flush=True)
    print(f"Gene counts: primary={int(masks['primary'].sum())} "
          f"lenient={int(masks['lenient'].sum())} unfiltered={int(masks['unfiltered'].sum())}",
          flush=True)
    # Permutations: joint (one age assignment per mouse across all units) and
    # independent (Phase 1 behaviour: permute within each unit, sex-stratified).
    joint_within = {name: [] for name in masks}
    joint_cross = {name: [] for name in masks}
    indep_within = {name: [] for name in masks}
    indep_cross = {name: [] for name in masks}
    tissue_idx = {t: [i for i, tt in enumerate(tissues) if tt == t] for t in sorted(set(tissues))}
    stratum_joint = {t: [] for t in tissue_idx}
    # Per-stratum detection sets (primary rule applied within the stratum's own units).
    stratum_masks = {t: detection_sets([all_units[i] for i in idx], genes)[0]
                     for t, idx in tissue_idx.items()}
    unit_donor_idx = [[donor_index[d] for d in u["metadata"]["donor"]] for u in all_units]
    unit_sex = [u["metadata"]["sex_male"].to_numpy() for u in all_units]
    unit_age = [u["metadata"]["age"].to_numpy() for u in all_units]
    # Float64 expression, converted once per unit: the permutation loop fits
    # every unit 2x per iteration, and converting float32 -> float64 on
    # every call was a large share of the cost.
    unit_expr64 = [np.asarray(u["expression"], dtype=np.float64) for u in all_units]
    # Checkpointing: the first full run vanished at ~82 min with no trace
    # (no OOM, no error, session record gone), so null draws are saved
    # every 200 iterations and --resume continues from them with a fresh
    # seed (null draws are iid, so resuming with a new seed stays valid).
    checkpoint_path = os.path.join(os.path.expanduser("~"), "workspace", "rev2-scratch",
                                   "rev2_concordance_checkpoint.npz")
    start_iteration = 0
    if resume_from and os.path.exists(checkpoint_path):
        saved = np.load(checkpoint_path, allow_pickle=True)
        if int(saved["permutations"]) == permutations:
            joint_within = {k: list(v) for k, v in saved["joint_within"].item().items()}
            joint_cross = {k: list(v) for k, v in saved["joint_cross"].item().items()}
            indep_within = {k: list(v) for k, v in saved["indep_within"].item().items()}
            indep_cross = {k: list(v) for k, v in saved["indep_cross"].item().items()}
            stratum_joint = {k: list(v) for k, v in saved["stratum_joint"].item().items()}
            start_iteration = int(saved["iterations_done"])
            rng = np.random.default_rng(seed + 1)  # fresh stream for the rest
            print(f"Resuming from checkpoint: {start_iteration}/{permutations} "
                  f"permutations already drawn", flush=True)
    for iteration in range(start_iteration, permutations):
        # Joint: permute ages across the whole mouse pool within sex.
        joint_age = donor_age.copy()
        for sex_value in np.unique(donor_sex):
            positions = np.where(donor_sex == sex_value)[0]
            joint_age[positions] = rng.permutation(joint_age[positions])
        joint_slopes = np.vstack([
            fit_age_slope_only(unit_expr64[k],
                               unit_design(all_units[k]["metadata"], age=joint_age[idx]))
            for k, idx in enumerate(unit_donor_idx)
        ])
        # Independent: permute within each unit, sex-stratified (Phase 1 flaw).
        indep_slopes = []
        for k, (u, age, sex) in enumerate(zip(all_units, unit_age, unit_sex)):
            permuted = age.copy()
            for sex_value in np.unique(sex):
                positions = np.where(sex == sex_value)[0]
                permuted[positions] = rng.permutation(permuted[positions])
            indep_slopes.append(fit_age_slope_only(
                unit_expr64[k], unit_design(u["metadata"], age=permuted)))
        indep_slopes = np.vstack(indep_slopes)
        for name, mask in masks.items():
            joint = stratum_correlations(joint_slopes, mask, unit_tissue)
            indep = stratum_correlations(indep_slopes, mask, unit_tissue)
            joint_within[name].append(joint["within_mean"])
            joint_cross[name].append(joint["cross_mean"])
            indep_within[name].append(indep["within_mean"])
            indep_cross[name].append(indep["cross_mean"])
        for tissue, idx in tissue_idx.items():
            if len(idx) >= 2:
                stratum_joint[tissue].append(stratum_correlations(
                    joint_slopes[idx], stratum_masks[tissue]["primary"],
                    np.array([tissue] * len(idx))
                )["within_mean"])
        if (iteration + 1) % 100 == 0:
            print(f"  permutations {iteration + 1}/{permutations}", flush=True)
        if (iteration + 1) % 200 == 0 or iteration + 1 == permutations:
            Path(checkpoint_path).parent.mkdir(parents=True, exist_ok=True)
            np.savez(checkpoint_path, permutations=permutations,
                     iterations_done=iteration + 1,
                     joint_within=joint_within, joint_cross=joint_cross,
                     indep_within=indep_within, indep_cross=indep_cross,
                     stratum_joint=stratum_joint)
    # Split-half reliability (noise ceiling) on the primary gene set.
    # Halves are fit with an age-only design (const + age): the full design
    # (const + age + sex + depth) has up to 4 parameters, which a 3-donor
    # half cannot fit. The age-only half-slopes are the same estimand up to
    # the small sex/depth adjustment, and the purpose here is the noise
    # level of the slope vector, not the point estimate.
    reliability = np.full(len(all_units), np.nan)
    for ui, unit in enumerate(all_units):
        meta = unit["metadata"]
        n = len(meta)
        expr64 = np.asarray(unit["expression"], dtype=np.float64)
        reps = []
        for _ in range(split_half_reps):
            order = rng.permutation(n)
            first, second = order[: n // 2], order[n // 2:]
            if len(first) < 3 or len(second) < 3:
                continue
            ages_a = meta["age"].to_numpy()[first]
            ages_b = meta["age"].to_numpy()[second]
            try:
                Xa = np.column_stack([np.ones(len(first)), ages_a])
                Xb = np.column_stack([np.ones(len(second)), ages_b])
                slope_a = fit_age_slope_only(expr64[first], Xa, age_column=1)[masks["primary"]]
                slope_b = fit_age_slope_only(expr64[second], Xb, age_column=1)[masks["primary"]]
            except ValueError:
                continue
            if np.std(slope_a) > 0 and np.std(slope_b) > 0:
                reps.append(np.corrcoef(slope_a, slope_b)[0, 1])
        if reps:
            raw = float(np.mean(reps))
            reliability[ui] = float(np.clip(2 * raw / (1 + raw), 0, 1))  # Spearman-Brown
    split_r_mean = float(np.nanmean(reliability)) if np.isfinite(reliability).any() else float("nan")
    ceiling_within, ceiling_cross = [], []
    for first, second in itertools.combinations(range(len(all_units)), 2):
        if not (np.isfinite(reliability[first]) and np.isfinite(reliability[second])):
            continue
        ceiling = float(np.sqrt(reliability[first] * reliability[second]))
        if unit_tissue[first] == unit_tissue[second]:
            ceiling_within.append(ceiling)
        else:
            ceiling_cross.append(ceiling)
    ceiling_within_mean = float(np.mean(ceiling_within)) if ceiling_within else float("nan")
    ceiling_cross_mean = float(np.mean(ceiling_cross)) if ceiling_cross else float("nan")
    print(f"Split-half: units with an estimate={int(np.isfinite(reliability).sum())}/{len(all_units)}, "
          f"mean Spearman-Brown reliability={split_r_mean:.4f}; noise ceiling within="
          f"{ceiling_within_mean:.4f} (n_pairs={len(ceiling_within)}) cross={ceiling_cross_mean:.4f} "
          f"(n_pairs={len(ceiling_cross)})", flush=True)
    # Per-stratum detail (primary set, joint null only).
    stratum_rows = []
    for tissue in sorted(set(tissues)):
        idx = tissue_idx[tissue]
        sub_units = [all_units[i] for i in idx]
        sub_masks = stratum_masks[tissue]
        sub_obs = stratum_correlations(observed_slopes[idx], sub_masks["primary"],
                                       np.array([tissue] * len(idx)))
        sub_null = np.array(stratum_joint[tissue]) if stratum_joint[tissue] else np.array([np.nan])
        stratum_rows.append({
            "species": "Mouse", "tissue": tissue, "n_units": len(idx),
            "n_pairs_within": sub_obs["within_n_pairs"],
            "n_genes_primary": int(sub_masks["primary"].sum()),
            "n_genes_unfiltered": int(sub_masks["unfiltered"].sum()),
            "observed_mean_pearson_primary": sub_obs["within_mean"],
            "observed_median_pearson_primary": sub_obs["within_median"],
            "joint_null_mean_primary": float(np.nanmean(sub_null)),
            "joint_null_sd_primary": float(np.nanstd(sub_null, ddof=1)),
            "joint_p_greater_primary": float(
                (1 + (sub_null >= sub_obs["within_mean"]).sum()) / (len(sub_null) + 1)),
            "n_permutations_joint": permutations,
            "note": "Per-stratum joint null from the same species-wide joint permutations.",
        })
    rows = []
    for scope, obs_key, null_within, null_cross, ceil in [
        ("within_tissue", "within_mean", joint_within, None, ceiling_within_mean),
        ("cross_tissue", "cross_mean", joint_cross, None, ceiling_cross_mean),
    ]:
        for name, mask in masks.items():
            joint_values = np.array((joint_within if scope == "within_tissue" else joint_cross)[name])
            indep_values = np.array((indep_within if scope == "within_tissue" else indep_cross)[name])
            obs_value = observed[name][obs_key]
            rows.append({
                "species": "Mouse", "scope": scope,
                "gene_set": ("primary_detection_filtered" if name == "primary"
                             else "sensitivity_lenient" if name == "lenient"
                             else "sensitivity_unfiltered"),
                "observed_mean_pearson": obs_value,
                "observed_median_pearson": observed[name]["within_median" if scope == "within_tissue" else "cross_median"],
                "observed_mean_spearman": (observed_spearman[obs_key] if name == "primary"
                                           else float("nan")),
                "n_pairs": (observed[name]["within_n_pairs"] if scope == "within_tissue"
                            else observed[name]["cross_n_pairs"]),
                "n_units": len(all_units), "n_donors": len(donors),
                "n_genes": int(mask.sum()),
                "n_genes_unfiltered": int(masks["unfiltered"].sum()),
                "detection_rule": ("detected in >=50% of donor pseudobulks in >=50% of units"
                                   if name == "primary"
                                   else "detected in >=20% of donor pseudobulks in >=25% of units"
                                   if name == "lenient" else "none (all genes after technical mask)"),
                "joint_null_mean": float(joint_values.mean()),
                "joint_null_sd": float(joint_values.std(ddof=1)),
                "joint_null_ci95_low": float(np.quantile(joint_values, 0.025)),
                "joint_null_ci95_high": float(np.quantile(joint_values, 0.975)),
                "joint_p_greater": float((1 + (joint_values >= obs_value).sum()) / (len(joint_values) + 1)),
                "n_permutations_joint": permutations,
                "indep_null_mean": float(indep_values.mean()),
                "indep_null_sd": float(indep_values.std(ddof=1)),
                "indep_p_greater": float((1 + (indep_values >= obs_value).sum()) / (len(indep_values) + 1)),
                "noise_ceiling_mean": ceil,
                "noise_ceiling_ratio_observed_over_ceiling": (
                    float(obs_value / ceil) if np.isfinite(ceil) and ceil > 0 else float("nan")),
                "split_half_reliability_mean": float(np.nanmean(reliability)),
                "status": "redone_joint",
                "notes": "Joint permutation: one age assignment per mouse, shared across "
                         "all its cell types, permuted within sex. Independent null is the "
                         "Phase 1 (flawed) within-cell-type permutation, shown for comparison. "
                         + ("Scrublet doublet removal ON (M1 settings, expected rate 0.06)."
                            if use_scrublet else
                            "Scrublet doublet removal OFF in this run: the pilot run was "
                            "OOM-killed by scrublet's per-tissue working set in this shared "
                            "VM, and the committed doublet fractions are <= 0.51% of cells "
                            "in any tissue (results/m1-qc-summary.csv), so omitting it moves "
                            "pseudobulk age slopes negligibly; see the Bladder scrublet "
                            "sensitivity in the coordinator report.")
                         + " Donor pool is 22 mice: pooled donor 3-M-8/9 contributes to "
                           "no unit in the prespecified primary set.",
            })
    return rows, stratum_rows

# ---------------------------------------------------------------------------
# NMR verification + species with no local data (Fix 2, remainder)
# ---------------------------------------------------------------------------

def nmr_and_unavailable_rows():
    """NMR: verify the committed exact-permutation null from the files (the
    Phase 1 NMR code already mapped one age assignment per donor across both
    compartments, i.e. it was donor-joint; see
    run_reviewer_phase1_gene_concordance.py analyse_stratum, NMR branch).
    Other species: no local raw data in this checkout; not redone."""
    rows = []
    phase1_summary = pd.read_csv(RESULTS / "reviewer-phase1-universality-summary.csv")
    phase1_units = pd.read_csv(RESULTS / "reviewer-phase1-universality-units.csv")
    phase1_pairs = pd.read_csv(RESULTS / "reviewer-phase1-universality-pairs.csv")
    perms = pd.read_csv(RESULTS / "reviewer-phase1-universality-permutations-nmr.csv")
    # Recompute the NMR null from the committed permutation rows (iteration 0 is observed).
    nmr_null = perms[perms["iteration"] != 0]["mean_pearson"]
    nmr_obs = float(perms[perms["iteration"] == 0]["mean_pearson"].iloc[0])
    nmr_units = phase1_units[phase1_units["species"] == "Naked mole-rat"]
    rows.append({
        "species": "Naked mole-rat", "scope": "within_tissue",
        "gene_set": "primary_detection_filtered",
        "observed_mean_pearson": nmr_obs,
        "observed_median_pearson": float(perms[perms["iteration"] == 0]["median_pearson"].iloc[0]),
        "observed_mean_spearman": float(perms[perms["iteration"] == 0]["mean_spearman"].iloc[0]),
        "n_pairs": 1, "n_units": len(nmr_units), "n_donors": 5,
        "n_genes": int(nmr_units["n_genes"].iloc[0]),
        "n_genes_unfiltered": int(nmr_units["n_genes"].iloc[0]),
        "detection_rule": "not re-derived (see notes)",
        "joint_null_mean": float(nmr_null.mean()),
        "joint_null_sd": float(nmr_null.std(ddof=1)),
        "joint_null_ci95_low": float(nmr_null.quantile(0.025)),
        "joint_null_ci95_high": float(nmr_null.quantile(0.975)),
        "joint_p_greater": float((1 + (nmr_null >= nmr_obs).sum()) / (len(nmr_null) + 1)),
        "n_permutations_joint": 10,
        "indep_null_mean": float("nan"), "indep_null_sd": float("nan"),
        "indep_p_greater": float("nan"),
        "noise_ceiling_mean": float("nan"),
        "noise_ceiling_ratio_observed_over_ceiling": float("nan"),
        "split_half_reliability_mean": float("nan"),
        "status": "already_joint_verified",
        "notes": "NMR was already donor-joint in Phase 1: all 10 assignments of 2 aged "
                 "of 5 animals were applied by donor across both compartments, so the "
                 "null mean +0.64 is a donor-level null, not a per-cell-type artefact. "
                 "Null recomputed from reviewer-phase1-universality-permutations-nmr.csv. "
                 "No split-half ceiling: 5 animals cannot be split into fittable halves.",
    })
    raw_checks = {
        "Fruit fly": ("data/raw/afca/afca-head.h5ad", "AFCA head/body h5ad"),
        "C. elegans": ("data/interim/worm-gao-integrated.h5ad", "Gao integrated h5ad"),
        "Zebrafish": ("data/interim/zebrafish-testis-5-22mo.h5ad", "zebrafish testis h5ad"),
        "Human": ("CellxGene Census 2025-11-08", "Census fetch (web, out of scope)"),
    }
    for species, (path, what) in raw_checks.items():
        obs = float(phase1_summary[(phase1_summary["species"] == species)
                                   & (phase1_summary["metric"] == "mean_pearson")]["observed"].iloc[0])
        units = phase1_units[phase1_units["species"] == species]
        rows.append({
            "species": species, "scope": "within_tissue",
            "gene_set": "primary_detection_filtered",
            "observed_mean_pearson": obs,
            "observed_median_pearson": float("nan"),
            "observed_mean_spearman": float("nan"),
            "n_pairs": int((phase1_pairs["species"] == species).sum()),
            "n_units": len(units), "n_donors": -1,
            "n_genes": int(units["n_genes"].iloc[0]) if len(units) else -1,
            "n_genes_unfiltered": int(units["n_genes"].iloc[0]) if len(units) else -1,
            "detection_rule": "not redone",
            "joint_null_mean": float("nan"), "joint_null_sd": float("nan"),
            "joint_null_ci95_low": float("nan"), "joint_null_ci95_high": float("nan"),
            "joint_p_greater": float("nan"), "n_permutations_joint": 0,
            "indep_null_mean": float("nan"), "indep_null_sd": float("nan"),
            "indep_p_greater": float("nan"),
            "noise_ceiling_mean": float("nan"),
            "noise_ceiling_ratio_observed_over_ceiling": float("nan"),
            "split_half_reliability_mean": float("nan"),
            "status": "not_redone_data_unavailable",
            "notes": f"NOT redone: no local {what} in this checkout ({path} absent). "
                     f"Observed value is the Phase 1 value from "
                     f"reviewer-phase1-universality-summary.csv (per-cell-type null), "
                     f"kept for reference only; its p-value is not valid under a joint null.",
        })
    return rows


# ---------------------------------------------------------------------------
# Mouse-vs-human interaction (Fix 3)
# ---------------------------------------------------------------------------

def fit_global_replica(records, metric, species, model, age_scale):
    """Reimplementation of run_cross_species_expression_vector.fit_global
    (validated: reproduces the committed global values to 4 decimals)."""
    import statsmodels.api as sm

    data = records.dropna(subset=[metric]).copy()
    data["unit"] = data["tissue"] + " | " + data["cohort"] + " | " + data["cell_type"]
    std = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 6:
            continue
        scale = part[metric].std(ddof=1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        part["outcome"] = (part[metric] - part[metric].mean()) / scale
        part["age_scaled"] = part["age"] / age_scale
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        part["depth_centered"] = part["mean_log_umi"] - part["mean_log_umi"].mean()
        std.append(part)
    data = pd.concat(std, ignore_index=True)
    data["donor_cluster"] = data["cohort"] + " | " + data["donor"]
    unit_dummies = pd.get_dummies(data["unit"], drop_first=True, dtype=float)
    preds = ["age_scaled", "sex_centered"]
    if model == "depth_adjusted":
        preds.append("depth_centered")
    preds = [p for p in preds if np.ptp(data[p]) > 1e-12]
    X = pd.concat([data[preds], unit_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    fit = sm.OLS(data["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": data["donor_cluster"]}, use_t=True)
    ci = fit.conf_int().loc["age_scaled"]
    return {"effect": float(fit.params["age_scaled"]), "se": float(fit.bse["age_scaled"]),
            "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
            "p": float(fit.pvalues["age_scaled"]),
            "n_donors": int(data["donor_cluster"].nunique()),
            "n_units": int(data["unit"].nunique()), "n_rows": len(data)}


def pooled_prep(data, metric, model, age_mode):
    """Standardise within unit and build the pooled mouse+human frame.
    age_mode='published' uses the published scales (mouse 0.4 y, human 12.25 y,
    i.e. ~1/6 of each species' observed lifespan); 'range' uses the observed
    age span (mouse 2.4167 y, human 100 y)."""
    scales = {"Mouse": 0.4, "Human": 12.25} if age_mode == "published" else \
        {"Mouse": 2.5 - 1 / 12, "Human": 100.0}
    data = data.dropna(subset=[metric]).copy()
    data["unit"] = data["tissue"] + " | " + data["cohort"] + " | " + data["cell_type"]
    std = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 6:
            continue
        scale = part[metric].std(ddof=1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        part["outcome"] = (part[metric] - part[metric].mean()) / scale
        part["age_scaled"] = part["age"] / part["species"].map(scales)
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        part["depth_centered"] = part["mean_log_umi"] - part["mean_log_umi"].mean()
        std.append(part)
    d = pd.concat(std, ignore_index=True)
    d["donor_cluster"] = d["species"] + " | " + d["cohort"] + " | " + d["donor"]
    d["is_human"] = (d["species"] == "Human").astype(float)
    d["age_x_human"] = d["age_scaled"] * d["is_human"]
    d["sex_x_human"] = d["sex_centered"] * d["is_human"]
    d["depth_x_human"] = d["depth_centered"] * d["is_human"]
    return d


def pooled_interaction(d, model):
    """Pooled OLS with unit fixed effects, species-specific sex/depth nuisance
    coefficients, and donor-cluster-robust SEs. age_scaled = mouse slope;
    age_x_human = human-minus-mouse interaction."""
    import statsmodels.api as sm

    unit_dummies = pd.get_dummies(d["unit"], drop_first=True, dtype=float)
    preds = ["age_scaled", "age_x_human", "sex_centered", "sex_x_human"]
    if model == "depth_adjusted":
        preds += ["depth_centered", "depth_x_human"]
    X = pd.concat([d[preds], unit_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    fit = sm.OLS(d["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": d["donor_cluster"]}, use_t=True)
    out = {}
    for term in ("age_scaled", "age_x_human"):
        ci = fit.conf_int().loc[term]
        out[term] = {"effect": float(fit.params[term]), "se": float(fit.bse[term]),
                     "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
                     "p": float(fit.pvalues[term])}
    out["n_clusters"] = int(d["donor_cluster"].nunique())
    out["n_rows"] = len(d)
    return out, X


def cluster_bootstrap_interaction(d, X, model, reps=1000, seed=20261108):
    """Donor-cluster bootstrap of the interaction point estimate. The
    standardisation is held at the full-sample values (fixed transformation);
    donor clusters are resampled within species and the pooled OLS refit.
    Computed by Frisch-Waugh-Lovell: unit dummies are partialled out by
    within-unit demeaning. The demeaning is redone inside every replicate
    on the resampled rows (exactly as refitting the dummies would), so the
    bootstrap estimates match a dense refit of the full dummy-variable OLS.
    Returns the array of bootstrap interaction estimates."""
    rng = np.random.default_rng(seed)
    y = d["outcome"].to_numpy(dtype=float)
    cont_cols = ["age_scaled", "age_x_human", "sex_centered", "sex_x_human"]
    if model == "depth_adjusted":
        cont_cols += ["depth_centered", "depth_x_human"]
    Xc = d[cont_cols].to_numpy(dtype=float)
    codes = pd.Categorical(d["unit"]).codes
    n_units = codes.max() + 1
    inter_col = cont_cols.index("age_x_human")
    clusters = d["donor_cluster"].to_numpy()
    species_of = d["species"].to_numpy()
    unique_clusters = {s: np.unique(clusters[species_of == s]) for s in ("Mouse", "Human")}
    cluster_rows = {c: np.where(clusters == c)[0] for c in np.unique(clusters)}
    estimates = []
    for _ in range(reps):
        rows = []
        for species in ("Mouse", "Human"):
            drawn = rng.choice(unique_clusters[species], size=len(unique_clusters[species]),
                               replace=True)
            for c in drawn:
                rows.append(cluster_rows[c])
        rows = np.concatenate(rows)
        cb = codes[rows]
        cnt = np.bincount(cb, minlength=n_units).astype(float)
        safe = np.maximum(cnt, 1)
        y_dm = y[rows] - np.bincount(cb, weights=y[rows], minlength=n_units)[cb] / safe[cb]
        X_dm = Xc[rows] - np.vstack([
            np.bincount(cb, weights=Xc[rows, j], minlength=n_units)[cb] / safe[cb]
            for j in range(Xc.shape[1])
        ]).T
        try:
            beta = np.linalg.solve(X_dm.T @ X_dm, X_dm.T @ y_dm)
        except np.linalg.LinAlgError:
            beta = np.linalg.pinv(X_dm.T @ X_dm) @ (X_dm.T @ y_dm)
        estimates.append(beta[inter_col])
    return np.array(estimates)


def analyse_interaction(bootstrap_reps=1000):
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    sub = donors[donors["species"].isin(["Mouse", "Human"])].copy()
    pooled_ids = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                     .query("is_pooled_id")["mouse.id"].astype(str))
    rows = []
    for age_mode in ("published", "range"):
        for model in ("unadjusted", "depth_adjusted"):
            d = pooled_prep(sub, "profile_hellinger", model, age_mode)
            result, X = pooled_interaction(d, model)
            mouse = fit_global_replica(sub[sub["species"] == "Mouse"], "profile_hellinger",
                                       "Mouse", model, 0.4 if age_mode == "published" else 2.5 - 1 / 12)
            human = fit_global_replica(sub[sub["species"] == "Human"], "profile_hellinger",
                                       "Human", model, 12.25 if age_mode == "published" else 100.0)
            boot = cluster_bootstrap_interaction(d, X, model, reps=bootstrap_reps)
            inter = result["age_x_human"]
            boot_ci = (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)))
            boot_p = float(min(1.0, 2 * min((boot <= 0).mean(), (boot >= 0).mean())))
            print(f"Interaction {age_mode}/{model}: {inter['effect']:+.4f} "
                  f"[{inter['ci_low']:+.4f}, {inter['ci_high']:+.4f}] p={inter['p']:.4g}; "
                  f"bootstrap CI [{boot_ci[0]:+.4f}, {boot_ci[1]:+.4f}] p={boot_p:.4g}", flush=True)
            rows.append({
                "analysis": f"primary_{age_mode}_{model}" if age_mode == "published" and model == "unadjusted"
                            else f"sensitivity_{age_mode}_{model}",
                "model": model, "age_scale": age_mode,
                "age_scale_mouse_years": 0.4 if age_mode == "published" else 2.5 - 1 / 12,
                "age_scale_human_years": 12.25 if age_mode == "published" else 100.0,
                "mouse_effect_sd": mouse["effect"], "mouse_se": mouse["se"],
                "mouse_ci95_low": mouse["ci_low"], "mouse_ci95_high": mouse["ci_high"],
                "mouse_p": mouse["p"],
                "human_effect_sd": human["effect"], "human_se": human["se"],
                "human_ci95_low": human["ci_low"], "human_ci95_high": human["ci_high"],
                "human_p": human["p"],
                "interaction_human_minus_mouse": inter["effect"],
                "interaction_se": inter["se"],
                "interaction_ci95_low": inter["ci_low"],
                "interaction_ci95_high": inter["ci_high"],
                "interaction_p": inter["p"],
                "interaction_bootstrap_ci95_low": boot_ci[0],
                "interaction_bootstrap_ci95_high": boot_ci[1],
                "interaction_bootstrap_p_two_sided": boot_p,
                "interaction_bootstrap_reps": bootstrap_reps,
                "n_mouse_donors": mouse["n_donors"], "n_human_donors": human["n_donors"],
                "n_donor_units": result["n_rows"], "n_clusters": result["n_clusters"],
                "human_ci_includes_mouse_point": bool(human["ci_low"] <= mouse["effect"] <= human["ci_high"]),
                "mouse_point_minus_human_ci_high": float(mouse["effect"] - human["ci_high"]),
                "equivalence_margin_0.10_sd": 0.10,
                "equivalence_conclusion": "",
                "notes": "Pooled OLS, unit fixed effects, species-specific sex/depth "
                         "coefficients, donor-cluster-robust SEs (species|cohort|donor). "
                         "Published scale: effect per 0.4 y (mouse) / 12.25 y (human), "
                         "the scales of the committed global fits.",
            })
    # Exclude-pooled-mice sensitivity (published scale, unadjusted).
    sub2 = sub[~sub["donor"].isin(pooled_ids)]
    d = pooled_prep(sub2, "profile_hellinger", "unadjusted", "published")
    result, X = pooled_interaction(d, "unadjusted")
    boot = cluster_bootstrap_interaction(d, X, "unadjusted", reps=bootstrap_reps)
    inter = result["age_x_human"]
    rows.append({
        "analysis": "sensitivity_published_unadjusted_exclude_pooled_mice",
        "model": "unadjusted", "age_scale": "published",
        "age_scale_mouse_years": 0.4, "age_scale_human_years": 12.25,
        "mouse_effect_sd": float("nan"), "mouse_se": float("nan"),
        "mouse_ci95_low": float("nan"), "mouse_ci95_high": float("nan"), "mouse_p": float("nan"),
        "human_effect_sd": float("nan"), "human_se": float("nan"),
        "human_ci95_low": float("nan"), "human_ci95_high": float("nan"), "human_p": float("nan"),
        "interaction_human_minus_mouse": inter["effect"], "interaction_se": inter["se"],
        "interaction_ci95_low": inter["ci_low"], "interaction_ci95_high": inter["ci_high"],
        "interaction_p": inter["p"],
        "interaction_bootstrap_ci95_low": float(np.quantile(boot, 0.025)),
        "interaction_bootstrap_ci95_high": float(np.quantile(boot, 0.975)),
        "interaction_bootstrap_p_two_sided": float(min(1.0, 2 * min((boot <= 0).mean(), (boot >= 0).mean()))),
        "interaction_bootstrap_reps": bootstrap_reps,
        "n_mouse_donors": 20, "n_human_donors": 479,
        "n_donor_units": result["n_rows"], "n_clusters": result["n_clusters"],
        "human_ci_includes_mouse_point": "",
        "mouse_point_minus_human_ci_high": "",
        "equivalence_margin_0.10_sd": 0.10,
        "equivalence_conclusion": "",
        "notes": "Same pooled model after dropping the 3 pooled mouse ids "
                 "(3-M-5/6, 3-M-7/8, 3-M-8/9); separate mouse/human effects not "
                 "refit for this row (see committed global file for those).",
    })
    # Equivalence verdict on the primary row, written once the numbers are in.
    primary = rows[0]
    gap = primary["mouse_effect_sd"] - primary["human_ci95_high"]
    primary["equivalence_conclusion"] = (
        "DROP the equivalence framing. The 0.10 SD margin has no biological "
        "justification in the record and is "
        f"{100 * 0.10 / primary['mouse_effect_sd']:.0f}% of the mouse effect itself "
        f"(+{primary['mouse_effect_sd']:.4f} SD): a human effect of "
        f"+{primary['human_ci95_high']:.4f} SD, i.e. "
        f"{100 * primary['human_ci95_high'] / primary['mouse_effect_sd']:.0f}% of the mouse "
        "effect, would still count as 'equivalent to zero'. The honest report "
        "is the interval: human "
        f"+{primary['human_effect_sd']:.4f} SD, 95% CI {primary['human_ci95_low']:+.4f} to "
        f"+{primary['human_ci95_high']:.4f}, which includes 0 and nearly reaches the mouse "
        f"point estimate (+{primary['mouse_effect_sd']:.4f}, just {gap:.4f} SD above the "
        "human CI upper bound). The data cannot distinguish no human effect "
        "from a mouse-sized human effect.")
    return rows

# ---------------------------------------------------------------------------
# Sign-test replacement (Fix 5)
# ---------------------------------------------------------------------------

def _unit_slope_ols(rows, metric="profile_hellinger"):
    """One unit's OLS slope of Hellinger distance on age (+ sex if identified).
    rows: DataFrame with age, sex_male, metric. Returns slope or None."""
    import statsmodels.api as sm

    part = rows.dropna(subset=[metric]).copy()
    if len(part) < 6:
        return None
    X = sm.add_constant(part[["age", "sex_male"]].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        X = sm.add_constant(part[["age"]].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    return float(sm.OLS(part[metric].astype(float), X).fit().params["age"])


def _global_effect_point(rows, metric="profile_hellinger", age_scale=0.4):
    """Pooled standardised global effect (point estimate only), same
    standardisation as fit_global: outcome z-scored within unit, age in
    units of age_scale, unit fixed effects, sex centred within unit."""
    data = rows.dropna(subset=[metric]).copy()
    data["unit"] = data["tissue"] + " | " + data["cohort"] + " | " + data["cell_type"]
    std = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 6:
            continue
        scale = part[metric].std(ddof=1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        part["outcome"] = (part[metric] - part[metric].mean()) / scale
        part["age_scaled"] = part["age"] / age_scale
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        std.append(part)
    if not std:
        return float("nan")
    data = pd.concat(std, ignore_index=True)
    unit_codes = pd.Categorical(data["unit"]).codes
    dummies = np.zeros((len(data), unit_codes.max()), dtype=float)
    for j in range(1, unit_codes.max() + 1):
        dummies[:, j - 1] = (unit_codes == j)
    X = np.column_stack([np.ones(len(data)), data["age_scaled"], data["sex_centered"], dummies])
    beta, *_ = np.linalg.lstsq(X, data["outcome"].to_numpy(), rcond=None)
    return float(beta[1])


def analyse_sign_test(bootstrap_reps=2000, seed=20261108):
    from scipy.stats import binomtest, spearmanr

    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    mouse = donors[donors["species"] == "Mouse"].copy()
    mouse["unit"] = mouse["tissue"] + " | " + mouse["cohort"] + " | " + mouse["cell_type"]
    units = sorted(mouse["unit"].unique())
    # Reproduce the Phase 1 unit slopes (62 units with >= 6 donors).
    unit_rows = {}
    observed_slopes = {}
    for unit in units:
        part = mouse[mouse["unit"] == unit]
        slope = _unit_slope_ols(part)
        if slope is not None:
            observed_slopes[unit] = slope
            unit_rows[unit] = part
    n_units = len(observed_slopes)
    n_positive = int(sum(s > 0 for s in observed_slopes.values()))
    exact_p = float(binomtest(n_positive, n_units, 0.5).pvalue)
    full_global = _global_effect_point(mouse)
    print(f"Sign-test reproduction: {n_positive}/{n_units} positive, exact p={exact_p:.6g}, "
          f"global effect={full_global:+.4f}", flush=True)
    # Donor-cluster bootstrap: resample the 23 mice with replacement; every
    # resampled mouse carries all of its unit rows with it (duplicates kept).
    # A unit needs >= 6 donors to fit, so the fittable denominator varies
    # slightly across replicates; the fraction positive among fittable units
    # is the clean estimand, the count is reported alongside it.
    # Numpy implementation: per-unit row blocks are gathered by donor
    # multiplicity and each unit slope is a small lstsq (identical to the
    # statsmodels OLS used for the observed slopes).
    rng = np.random.default_rng(seed)
    mice = sorted(mouse["donor"].unique())
    mouse_rows = {m: mouse[mouse["donor"] == m] for m in mice}
    unit_arrays = {}
    for unit, part in unit_rows.items():
        unit_arrays[unit] = {
            "donor": part["donor"].to_numpy(),
            "age": part["age"].to_numpy(dtype=float),
            "sex": part["sex_male"].to_numpy(dtype=float),
            "y": part["profile_hellinger"].to_numpy(dtype=float),
        }
    mouse_unit_codes = pd.Categorical(mouse["unit"]).codes
    mouse_donor_codes = pd.Categorical(mouse["donor"], categories=mice).codes
    mouse_y_all = mouse["profile_hellinger"].to_numpy(dtype=float)
    mouse_age_all = mouse["age"].to_numpy(dtype=float)
    mouse_sex_all = mouse["sex_male"].to_numpy(dtype=float)
    mouse_donor_all = mouse["donor"].to_numpy()
    mouse_row_by_donor = {m: np.where(mouse_donor_all == m)[0] for m in mice}
    n_unit_codes = mouse_unit_codes.max() + 1

    def _slope_from_rows(age, sex, y):
        if len(y) < 6:
            return None
        X = np.column_stack([np.ones(len(y)), age, sex])
        if np.linalg.matrix_rank(X) < X.shape[1]:
            X = np.column_stack([np.ones(len(y)), age])
        if np.linalg.matrix_rank(X) < X.shape[1]:
            return None
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(beta[1])

    def _global_effect_from_rows(rows):
        if len(rows) == 0:
            return float("nan")
        codes = mouse_unit_codes[rows]
        yv = mouse_y_all[rows]
        counts = np.bincount(codes, minlength=n_unit_codes).astype(float)
        # Unique donors per unit in this resample (duplicates do not count).
        # NB: Categorical codes are int8; cast before multiplying (93*23 overflows int8).
        pair_codes = (codes.astype(np.int64) * len(mice)
                      + mouse_donor_codes[rows].astype(np.int64))
        _, pair_first = np.unique(pair_codes, return_index=True)
        uniq_per_unit = np.bincount(codes[pair_first], minlength=n_unit_codes)
        keep = uniq_per_unit >= 6
        means = np.bincount(codes, weights=yv, minlength=n_unit_codes) / np.maximum(counts, 1)
        sq = np.bincount(codes, weights=yv ** 2, minlength=n_unit_codes) / np.maximum(counts, 1)
        sds = np.sqrt(np.maximum(sq - means ** 2, 0) * counts / np.maximum(counts - 1, 1))
        valid = keep & (sds > 0)
        if not valid[codes].all():
            rows = rows[valid[codes]]
            codes = codes[valid[codes]]
            yv = mouse_y_all[rows]
        outcome = (yv - means[codes]) / sds[codes]
        age_scaled = mouse_age_all[rows] / 0.4
        sex = mouse_sex_all[rows]
        sex_mean = np.bincount(codes, weights=sex, minlength=n_unit_codes)[codes] / counts[codes]
        sex_centered = sex - sex_mean
        dummies = np.zeros((len(rows), int(valid.sum()) - 1), dtype=float)
        valid_units = np.where(valid)[0]
        col_of = {u: j for j, u in enumerate(valid_units[1:])}
        for j, u in enumerate(valid_units[1:]):
            dummies[:, j] = (codes == u)
        X = np.column_stack([np.ones(len(rows)), age_scaled, sex_centered, dummies])
        beta, *_ = np.linalg.lstsq(X, outcome, rcond=None)
        return float(beta[1])

    boot_positive, boot_fittable, boot_global = [], [], []
    for _ in range(bootstrap_reps):
        drawn = rng.choice(mice, size=len(mice), replace=True)
        _, draw_counts = np.unique(drawn, return_counts=True)
        drawn_unique = np.unique(drawn)
        count_map = dict(zip(drawn_unique, draw_counts))
        count = 0
        fittable = 0
        for unit, arr in unit_arrays.items():
            idx_parts = []
            for donor in np.unique(arr["donor"]):
                k = count_map.get(donor, 0)
                if k:
                    base = np.where(arr["donor"] == donor)[0]
                    idx_parts.append(np.repeat(base, k))
            if not idx_parts:
                continue
            idx = np.concatenate(idx_parts)
            if len(np.unique(arr["donor"][idx])) < 6:
                continue
            slope = _slope_from_rows(arr["age"][idx], arr["sex"][idx], arr["y"][idx])
            if slope is not None:
                fittable += 1
                if slope > 0:
                    count += 1
        boot_positive.append(count)
        boot_fittable.append(fittable)
        all_idx = np.concatenate([np.repeat(mouse_row_by_donor[m], count_map[m])
                                  for m in drawn_unique])
        boot_global.append(_global_effect_from_rows(all_idx))
    boot_positive = np.array(boot_positive)
    boot_fittable = np.array(boot_fittable)
    boot_fraction = boot_positive / np.maximum(boot_fittable, 1)
    boot_global = np.array(boot_global)
    boot_p_le_half = float((boot_fraction <= 0.5).mean())
    print(f"Bootstrap (B={bootstrap_reps}): n_positive mean={boot_positive.mean():.1f} "
          f"SD={boot_positive.std(ddof=1):.1f} 95% PI=[{np.quantile(boot_positive, 0.025):.0f}, "
          f"{np.quantile(boot_positive, 0.975):.0f}] range=[{boot_positive.min()}, "
          f"{boot_positive.max()}]; fittable units range=[{boot_fittable.min()}, "
          f"{boot_fittable.max()}]; fraction positive mean={boot_fraction.mean():.3f} "
          f"95% PI=[{np.quantile(boot_fraction, 0.025):.3f}, {np.quantile(boot_fraction, 0.975):.3f}]; "
          f"P(fraction <= 0.5)={boot_p_le_half:.4f}; "
          f"global effect mean={boot_global.mean():+.4f} "
          f"95% PI=[{np.quantile(boot_global, 0.025):+.4f}, {np.quantile(boot_global, 0.975):+.4f}]",
          flush=True)
    # Leave-one-donor-out and leave-one-tissue-out.
    def _count_positive(sub):
        count = 0
        fittable = 0
        for unit in unit_rows:
            slope = _unit_slope_ols(sub[sub["unit"] == unit])
            if slope is not None:
                fittable += 1
                if slope > 0:
                    count += 1
        return count, fittable

    lodo_sign, lodo_frac, lodo_global, lodo_worst = [], [], [], None
    for m in mice:
        sub = mouse[mouse["donor"] != m]
        count, fittable = _count_positive(sub)
        effect = _global_effect_point(sub)
        lodo_sign.append(count)
        lodo_frac.append(count / max(fittable, 1))
        lodo_global.append(effect)
        if lodo_worst is None or count < lodo_worst[1]:
            lodo_worst = (m, count, effect, fittable)
    loto_sign, loto_frac, loto_global, loto_worst = [], [], [], None
    for tissue in sorted(mouse["tissue"].unique()):
        sub = mouse[mouse["tissue"] != tissue]
        count, fittable = _count_positive(sub)
        effect = _global_effect_point(sub)
        loto_sign.append(count)
        loto_frac.append(count / max(fittable, 1))
        loto_global.append(effect)
        if loto_worst is None or count < loto_worst[1]:
            loto_worst = (tissue, count, effect, fittable)
    print(f"LODO sign count range [{min(lodo_sign)}, {max(lodo_sign)}] "
          f"(worst: drop {lodo_worst[0]} -> {lodo_worst[1]}); "
          f"LODO global [{min(lodo_global):+.4f}, {max(lodo_global):+.4f}]", flush=True)
    print(f"LOTO sign count range [{min(loto_sign)}, {max(loto_sign)}] "
          f"(worst: drop {loto_worst[0]} -> {loto_worst[1]}); "
          f"LOTO global [{min(loto_global):+.4f}, {max(loto_global):+.4f}]", flush=True)
    # One-observation-per-donor summary: mean within-unit standardised
    # Hellinger displacement per mouse, regressed on age.
    import statsmodels.api as sm

    mouse["z"] = mouse.groupby("unit")["profile_hellinger"].transform(
        lambda s: (s - s.mean()) / s.std(ddof=1))
    per = mouse.groupby("donor").agg(mean_z=("z", "mean"), age=("age", "first"),
                                     sex=("sex_male", "first"),
                                     n_units=("z", "size")).reset_index()
    X = sm.add_constant(per[["age", "sex"]].astype(float), has_constant="add")
    fit = sm.OLS(per["mean_z"], X).fit()
    ci = fit.conf_int().loc["age"]
    rho = spearmanr(per["age"], per["mean_z"])
    boot_slopes = []
    for _ in range(bootstrap_reps):
        sample = per.iloc[rng.integers(0, len(per), size=len(per))]
        try:
            Xb = sm.add_constant(sample[["age", "sex"]].astype(float), has_constant="add")
            boot_slopes.append(float(sm.OLS(sample["mean_z"], Xb).fit().params["age"]))
        except Exception:
            continue
    boot_slopes = np.array(boot_slopes)
    print(f"Per-donor mean-z ~ age: slope={fit.params['age']:+.4f} "
          f"95% CI [{ci.iloc[0]:+.4f}, {ci.iloc[1]:+.4f}] p={fit.pvalues['age']:.4g}; "
          f"Spearman rho={rho.statistic:.3f} p={rho.pvalue:.4g}; "
          f"bootstrap slope 95% CI [{np.quantile(boot_slopes, 0.025):+.4f}, "
          f"{np.quantile(boot_slopes, 0.975):+.4f}]", flush=True)
    rows = [
        {"analysis": "exact_sign_test_reference", "statistic": "n_positive_units",
         "n": n_units, "observed": n_positive, "mean": "", "sd": "",
         "ci95_low": "", "ci95_high": "", "min": "", "max": "",
         "p_value": exact_p,
         "details": f"{n_positive} of {n_units} mouse tissue-cell-type units had positive "
                    "unadjusted Hellinger age slopes (reproduced exactly).",
         "recommendation": "DROP. Treats 62 correlated unit signs as independent Bernoulli "
                           "trials; the same 23 mice contribute to many units, so the "
                           "p-value is invalid. Do not report it as evidence."},
        {"analysis": "donor_cluster_bootstrap_sign_count", "statistic": "n_positive_units",
         "n": n_units, "observed": n_positive,
         "mean": float(boot_positive.mean()), "sd": float(boot_positive.std(ddof=1)),
         "ci95_low": float(np.quantile(boot_positive, 0.025)),
         "ci95_high": float(np.quantile(boot_positive, 0.975)),
         "min": int(boot_positive.min()), "max": int(boot_positive.max()),
         "p_value": boot_p_le_half,
         "details": f"Donor-cluster bootstrap (B={bootstrap_reps}, resampling the 23 mice "
                    f"with replacement, each mouse carrying all its unit rows). "
                    f"Fraction positive among fittable units: mean "
                    f"{boot_fraction.mean():.3f}, 95% PI "
                    f"[{np.quantile(boot_fraction, 0.025):.3f}, "
                    f"{np.quantile(boot_fraction, 0.975):.3f}]; fittable units per "
                    f"replicate {boot_fittable.min()}-{boot_fittable.max()} of {n_units}. "
                    f"p_value = fraction of replicates with fraction positive <= 0.5 "
                    "(no majority positive).",
         "recommendation": "REPORT this instead of the exact sign test."},
        {"analysis": "donor_cluster_bootstrap_global_effect", "statistic": "effect_sd_per_0.4y",
         "n": len(mice), "observed": full_global,
         "mean": float(boot_global.mean()), "sd": float(boot_global.std(ddof=1)),
         "ci95_low": float(np.quantile(boot_global, 0.025)),
         "ci95_high": float(np.quantile(boot_global, 0.975)),
         "min": float(boot_global.min()), "max": float(boot_global.max()),
         "p_value": float((boot_global <= 0).mean()),
         "details": f"Same donor-cluster bootstrap (B={bootstrap_reps}); global effect is "
                    "the pooled standardised Hellinger slope per 0.4 y. p_value = fraction "
                    "of replicates <= 0.",
         "recommendation": "REPORT alongside the sign-count bootstrap."},
        {"analysis": "leave_one_donor_out_sign_count", "statistic": "n_positive_units",
         "n": len(mice), "observed": n_positive, "mean": float(np.mean(lodo_sign)),
         "sd": float(np.std(lodo_sign, ddof=1)),
         "ci95_low": "", "ci95_high": "", "min": int(min(lodo_sign)), "max": int(max(lodo_sign)),
         "p_value": "",
         "details": f"Dropping one mouse at a time ({len(mice)} mice). Worst case: drop "
                    f"{lodo_worst[0]} -> {lodo_worst[1]} positive of {lodo_worst[3]} fittable "
                    f"units (fraction {lodo_worst[1] / max(lodo_worst[3], 1):.3f}). "
                    f"Fraction positive across LODO runs: "
                    f"{min(lodo_frac):.3f}-{max(lodo_frac):.3f}.",
         "recommendation": "REPORT the range."},
        {"analysis": "leave_one_donor_out_global_effect", "statistic": "effect_sd_per_0.4y",
         "n": len(mice), "observed": full_global, "mean": float(np.mean(lodo_global)),
         "sd": float(np.std(lodo_global, ddof=1)),
         "ci95_low": "", "ci95_high": "", "min": float(min(lodo_global)),
         "max": float(max(lodo_global)), "p_value": "",
         "details": f"Dropping one mouse at a time. Worst case: drop {lodo_worst[0]} -> "
                    f"{lodo_worst[2]:+.4f} SD.",
         "recommendation": "REPORT the range."},
        {"analysis": "leave_one_tissue_out_sign_count", "statistic": "n_positive_units",
         "n": mouse["tissue"].nunique(), "observed": n_positive,
         "mean": float(np.mean(loto_sign)), "sd": float(np.std(loto_sign, ddof=1)),
         "ci95_low": "", "ci95_high": "", "min": int(min(loto_sign)), "max": int(max(loto_sign)),
         "p_value": "",
         "details": f"Dropping one tissue at a time ({mouse['tissue'].nunique()} tissues). "
                    f"Worst case: drop {loto_worst[0]} -> {loto_worst[1]} positive of "
                    f"{loto_worst[3]} fittable units (fraction "
                    f"{loto_worst[1] / max(loto_worst[3], 1):.3f}). Fraction positive across "
                    f"LOTO runs: {min(loto_frac):.3f}-{max(loto_frac):.3f}.",
         "recommendation": "REPORT the range."},
        {"analysis": "leave_one_tissue_out_global_effect", "statistic": "effect_sd_per_0.4y",
         "n": mouse["tissue"].nunique(), "observed": full_global,
         "mean": float(np.mean(loto_global)), "sd": float(np.std(loto_global, ddof=1)),
         "ci95_low": "", "ci95_high": "", "min": float(min(loto_global)),
         "max": float(max(loto_global)), "p_value": "",
         "details": f"Dropping one tissue at a time. Worst case: drop {loto_worst[0]} -> "
                    f"{loto_worst[2]:+.4f} SD.",
         "recommendation": "REPORT the range."},
        {"analysis": "per_donor_mean_standardised_displacement",
         "statistic": "slope_of_per_donor_mean_z_on_age_per_year",
         "n": len(per), "observed": float(fit.params["age"]),
         "mean": "", "sd": "",
         "ci95_low": float(ci.iloc[0]), "ci95_high": float(ci.iloc[1]),
         "min": "", "max": "", "p_value": float(fit.pvalues["age"]),
         "details": f"One observation per mouse (n={len(per)}): mean within-unit "
                    f"standardised Hellinger distance, OLS on age + sex. Spearman rho "
                    f"between age and per-donor mean z = {rho.statistic:.3f} "
                    f"(p={rho.pvalue:.4g}). Bootstrap (B={bootstrap_reps}) slope 95% CI "
                    f"[{np.quantile(boot_slopes, 0.025):+.4f}, "
                    f"{np.quantile(boot_slopes, 0.975):+.4f}].",
         "recommendation": "REPORT as the donor-level sign summary; note it is weak "
                           "(n=23 mice, p=0.065), far weaker than the dropped sign test implied."},
    ]
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def emit(frame, name):
    """Print a result table to stdout (the reliable channel when the VM's
    file writes are intermittently blocked) and then try the file write."""
    print(f"---BEGIN CSV {name}---", flush=True)
    print(frame.to_csv(index=False), flush=True)
    print(f"---END CSV {name}---", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--split-half-reps", type=int, default=100)
    parser.add_argument("--bootstrap-reps", type=int, default=2000)
    parser.add_argument("--interaction-bootstrap-reps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--skip-mouse", action="store_true")
    parser.add_argument("--mouse-only", action="store_true")
    parser.add_argument("--skip-interaction", action="store_true")
    parser.add_argument("--skip-sign-test", action="store_true")
    parser.add_argument("--no-scrublet", action="store_true",
                        help="Skip per-tissue scrublet in the mouse QC (see mouse_strata). "
                             "Default is to run scrublet as in M1; use this only when the "
                             "VM cannot hold scrublet's per-tissue working set.")
    parser.add_argument("--tissues", default="",
                        help="Comma-separated tissue subset for the mouse run (sensitivity).")
    parser.add_argument("--resume", action="store_true",
                        help="Resume the mouse permutation loop from the checkpoint npz "
                             "in ~/workspace/rev2-scratch (fresh RNG stream for the rest).")
    args = parser.parse_args()

    if not args.skip_mouse:
        mouse_rows, stratum_rows = analyse_mouse(
            permutations=args.permutations, split_half_reps=args.split_half_reps, seed=args.seed,
            use_scrublet=not args.no_scrublet, resume_from=args.resume,
            tissues_only=[t for t in args.tissues.split(",") if t] or None)
        mouse_rows.extend(nmr_and_unavailable_rows())
        concordance = pd.DataFrame(mouse_rows)
        emit(concordance, "reviewer-rev2-concordance.csv")
        robust_write(concordance, RESULTS / "reviewer-rev2-concordance.csv")
        strata = pd.DataFrame(stratum_rows)
        emit(strata, "reviewer-rev2-concordance-strata.csv")
        robust_write(strata, RESULTS / "reviewer-rev2-concordance-strata.csv")
    if args.mouse_only:
        return
    if not args.skip_interaction:
        interaction = pd.DataFrame(analyse_interaction(
            bootstrap_reps=args.interaction_bootstrap_reps))
        emit(interaction, "reviewer-rev2-mouse-human-interaction.csv")
        robust_write(interaction, RESULTS / "reviewer-rev2-mouse-human-interaction.csv")
    if not args.skip_sign_test:
        sign_test = pd.DataFrame(analyse_sign_test(
            bootstrap_reps=args.bootstrap_reps, seed=args.seed))
        emit(sign_test, "reviewer-rev2-sign-test-replacement.csv")
        robust_write(sign_test, RESULTS / "reviewer-rev2-sign-test-replacement.csv")


if __name__ == "__main__":
    main()
