#!/usr/bin/env python3
"""Rev2 fix 4: Hellinger whole-expression-vector robustness package.

Reanalysis of the cross-species whole-expression-vector result
(results/h3-cross-species-vector-global.csv: mouse profile-Hellinger global
effect +0.1028 SD per 10% of maximum lifespan unadjusted, +0.0679 with a
linear mean-log-UMI depth covariate; human +0.0365 unadjusted, +0.0231
depth-adjusted).

What this script does (CPU only, no new data, no new experiments):
  Part A (donor level, from results/h3-cross-species-vector-donors.csv):
    - reproduces the published baseline global fits (sanity check);
    - nonlinear depth adjustment: 1/sqrt(UMI), quadratic log-UMI, and a
      3-df cubic spline on log-UMI, compared with the linear adjustment;
    - per-organ (mouse) and per-cohort / per-tissue (human) decomposition;
    - young-only age slope (mouse 1m vs 3m; human 20-40y) as a
      donor-level developmental-drift check.
  Part B (cell level, mouse Tabula Muris Senis):
    - rebuilds the M1 QC-cleaned object from the Figshare droplet raw file
      (data/raw/tms_droplet_raw.h5ad, sha256 per data/manifest/
      tms_droplet_raw.json) using only the M1 Step-3 filters
      (500-5000 genes, >=1000 UMIs, per-tissue scrublet, seed 42);
    - fixed-UMI downsampled vector analysis at 500 / 1500 / 3000 UMIs;
    - distance-to-median dispersion and median distance-to-centroid;
    - donor-by-cell-type pseudobulk Hellinger distances;
    - young-vs-young negative controls (split-half young reference;
      1m-only reference vs 3m donors);
    - alternative young references, including a mature-adult 18m reference.
  Part C (cell level, naked mole-rat skin, GSE152886, 28 MB raw file):
    - same downsampled / pseudobulk / dispersion endpoints where the
      5-sample design allows them.

Human cell-level re-fetch (CELLxGENE census) is NOT attempted here: it
requires re-downloading the 7 census datasets donor-by-donor (the original
run's scale) and is out of scope for this worker; human results below are
donor-level only, using the committed donor file. Fly (AFCA, ~1 GB),
worm (Gao 8.2 GB, Roux 4.1 GB) and zebrafish (528 MB RAW tar) cell-level
re-fetches are likewise not attempted; see the final report for reasons.

Outputs (new files only):
  results/reviewer-rev2-hellinger.csv          one row per analysis variant
  results/reviewer-rev2-hellinger-donors.csv   donor-unit rows behind Part B/C
  (Part A uses the committed results/h3-cross-species-vector-donors.csv.)

Usage:
  python scripts/run_rev2_hellinger.py --donor-only
  python scripts/run_rev2_hellinger.py --cell-level   (needs TMS raw + scanpy)
  python scripts/run_rev2_hellinger.py --nmr          (needs GSE152886 raw)
  python scripts/run_rev2_hellinger.py                (all available parts)

No .tex file and no manuscript file is touched. Nothing is committed.
"""

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

RESULTS = REPO / "results"
MOUSE_AGES = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
RNG_SEED = 42


# --------------------------------------------------------------------------
# Shared global-fit machinery (generalised from
# scripts/run_cross_species_expression_vector.py::fit_global)
# --------------------------------------------------------------------------

def _depth_features(data, adjustment):
    """Return (column frame, column names) for a depth-adjustment variant.

    Depth is centred within unit (tissue | cohort | cell_type). 1/sqrt(UMI)
    uses the donor-unit arithmetic mean UMI when available, else the
    geometric mean exp(mean_log_umi).
    """
    if adjustment == "unadjusted":
        return pd.DataFrame(index=data.index), []
    depth_c = data["mean_log_umi"] - data.groupby("unit")["mean_log_umi"].transform("mean")
    if "mean_umi" in data.columns and data["mean_umi"].notna().all():
        inv_sqrt = 1.0 / np.sqrt(data["mean_umi"])
    else:
        inv_sqrt = np.exp(-0.5 * data["mean_log_umi"])
    inv_sqrt = pd.Series(np.asarray(inv_sqrt, dtype=float), index=data.index)
    inv_sqrt_c = inv_sqrt - inv_sqrt.groupby(data["unit"].values).transform("mean")
    cols = {}
    if adjustment == "linear_log_umi":
        cols["depth_c"] = depth_c
    elif adjustment == "inv_sqrt_umi":
        cols["inv_sqrt_c"] = inv_sqrt_c
    elif adjustment == "linear_plus_inv_sqrt":
        cols["depth_c"] = depth_c
        cols["inv_sqrt_c"] = inv_sqrt_c
    elif adjustment == "quadratic_log_umi":
        cols["depth_c"] = depth_c
        cols["depth_c2"] = depth_c ** 2
    elif adjustment == "spline_log_umi":
        try:
            from patsy import dmatrix
            basis = dmatrix(
                "cr(depth_c, df=4) - 1",
                {"depth_c": depth_c.to_numpy()},
                return_type="dataframe",
            )
            basis.index = data.index
            # B-spline bases sum to a constant (partition of unity), which is
            # collinear with the intercept: drop one column, leaving the
            # 3 free nonlinear directions of the 4-df spline.
            basis = basis.iloc[:, 1:]
            for name in basis.columns:
                cols[f"spline_{name}"] = basis[name]
        except Exception:
            cols["depth_c"] = depth_c
            cols["depth_c2"] = depth_c ** 2
            cols["depth_c3"] = depth_c ** 3
    else:
        raise ValueError(f"unknown adjustment {adjustment}")
    return pd.DataFrame(cols, index=data.index), list(cols)


def fit_global2(records, metric, species, tissue_label, cohort_label, adjustment,
                age_scale, variant, depth_label="", reference="young_1_3m",
                donor_set="all_donor_ids", endpoint=None, extra=None,
                source="donor_file"):
    """Standardised global fit: outcome SD within unit, unit intercepts,
    age scaled by age_scale, cluster-robust SEs by donor (cohort | donor).

    extra: optional dict {name: column in records} added as additional
    predictors, centred within unit (e.g. a platform indicator)."""
    data = records.dropna(subset=[metric]).copy()
    if data.empty:
        return None
    data["unit"] = data["tissue"].astype(str) + " | " + data["cohort"].astype(str) \
        + " | " + data["cell_type"].astype(str)
    standardized = []
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
        standardized.append(part)
    if not standardized:
        return None
    data = pd.concat(standardized, ignore_index=True)
    data["donor_cluster"] = data["cohort"].astype(str) + " | " + data["donor"].astype(str)
    if data["donor_cluster"].nunique() < 6:
        return None
    depth_frame, depth_cols = _depth_features(data, adjustment)
    keep_depth = [c for c in depth_cols if np.ptp(depth_frame[c].to_numpy()) > 1e-12]
    extra_frame = pd.DataFrame(index=data.index)
    if extra:
        for name, column in extra.items():
            values = pd.Series(np.asarray(data[column], dtype=float), index=data.index)
            centred = values - values.groupby(data["unit"].values).transform("mean")
            if np.ptp(centred.to_numpy()) > 1e-12:
                extra_frame[name] = centred
    unit_dummies = pd.get_dummies(data["unit"], drop_first=True, dtype=float)
    base = data[["age_scaled", "sex_centered"]]
    base = base[[c for c in base.columns if np.ptp(base[c].to_numpy()) > 1e-12]]
    X = pd.concat([base, depth_frame[keep_depth],
                   extra_frame, unit_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    fit = sm.OLS(data["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": data["donor_cluster"]}, use_t=True
    )
    ci = fit.conf_int().loc["age_scaled"]
    return {
        "endpoint": endpoint or metric,
        "species": species, "tissue": tissue_label, "cohort": cohort_label,
        "variant": variant, "depth": depth_label, "reference": reference,
        "adjustment": adjustment, "donor_set": donor_set,
        "age_scale_years": age_scale,
        "n_donors": int(data["donor_cluster"].nunique()),
        "n_units": int(data["unit"].nunique()),
        "n_donor_units": int(len(data)),
        "effect_sd": float(fit.params["age_scaled"]),
        "se": float(fit.bse["age_scaled"]),
        "ci95_low": float(ci.iloc[0]), "ci95_high": float(ci.iloc[1]),
        "pvalue": float(fit.pvalues["age_scaled"]),
        "source": source,
        "notes": "",
    }


def young_only_slope(records, metric, species, age_unit_label, age_scale_months=None):
    """Age slope restricted to young donors only (developmental-drift check).

    Mouse: 1m vs 3m donors (ages in years in `age`; slope reported per month).
    Human: donors aged 20-40y (slope reported per year). Outcome standardised
    within unit; unit intercepts; cluster SEs by donor.
    """
    data = records.dropna(subset=[metric]).copy()
    data["unit"] = data["tissue"].astype(str) + " | " + data["cohort"].astype(str) \
        + " | " + data["cell_type"].astype(str)
    rows = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 4 or part["age"].nunique() < 2:
            continue
        scale = part[metric].std(ddof=1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        part["outcome"] = (part[metric] - part[metric].mean()) / scale
        rows.append(part)
    if not rows:
        return None
    data = pd.concat(rows, ignore_index=True)
    data["donor_cluster"] = data["cohort"].astype(str) + " | " + data["donor"].astype(str)
    if data["donor_cluster"].nunique() < 6 or data["age"].nunique() < 2:
        return None
    data["sex_centered"] = data["sex_male"] - data["sex_male"].mean()
    data["depth_centered"] = data["mean_log_umi"] - data.groupby("unit")["mean_log_umi"].transform("mean")
    unit_dummies = pd.get_dummies(data["unit"], drop_first=True, dtype=float)
    X = pd.concat([data[["age", "sex_centered", "depth_centered"]], unit_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    fit = sm.OLS(data["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": data["donor_cluster"]}, use_t=True
    )
    ci = fit.conf_int().loc["age"]
    return {
        "endpoint": metric + "_young_only_slope",
        "species": species, "tissue": "All tested tissues", "cohort": "All cohorts",
        "variant": "young_only", "depth": "full", "reference": "young_1_3m" if species == "Mouse" else "young_le40",
        "adjustment": "linear_log_umi", "donor_set": "all_donor_ids",
        "age_scale_years": age_scale_months,
        "n_donors": int(data["donor_cluster"].nunique()),
        "n_units": int(data["unit"].nunique()),
        "n_donor_units": int(len(data)),
        "effect_sd": float(fit.params["age"]),
        "se": float(fit.bse["age"]),
        "ci95_low": float(ci.iloc[0]), "ci95_high": float(ci.iloc[1]),
        "pvalue": float(fit.pvalues["age"]),
        "source": "donor_file",
        "notes": f"Age slope within young donors only ({age_unit_label}); "
                 "a nonzero slope here would indicate drift inside the reference window.",
    }


# --------------------------------------------------------------------------
# Part A: donor-level analyses on the committed donor file
# --------------------------------------------------------------------------

def part_a_donor_level(out_rows):
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    donors = donors[donors["species"] != "Naked mole-rat"].copy()
    adjustments = ["unadjusted", "linear_log_umi", "inv_sqrt_umi",
                   "linear_plus_inv_sqrt", "quadratic_log_umi", "spline_log_umi"]
    for species, part in donors.groupby("species"):
        age_scale = 0.4 if species == "Mouse" else 12.25
        # --- global, all tissues/cohorts, every depth adjustment ---
        for adj in adjustments:
            row = fit_global2(part, "profile_hellinger", species,
                              "All tested tissues", "All cohorts", adj,
                              age_scale, variant="global", depth_label="full",
                              endpoint="profile_hellinger_centroid_mean")
            if row:
                out_rows.append(row)
        # --- per-organ (mouse) / per-tissue (human) decomposition ---
        for tissue, tpart in part.groupby("tissue"):
            for adj in ("unadjusted", "linear_log_umi", "inv_sqrt_umi", "spline_log_umi"):
                row = fit_global2(tpart, "profile_hellinger", species, tissue,
                                  "All cohorts", adj, age_scale,
                                  variant="per_tissue", depth_label="full",
                                  endpoint="profile_hellinger_centroid_mean")
                if row:
                    out_rows.append(row)
        # --- per-cohort (human study/platform view; each cohort is one
        #     CELLxGENE dataset, so cohort FE = dataset/platform FE) ---
        for cohort, cpart in part.groupby("cohort"):
            for adj in ("unadjusted", "linear_log_umi"):
                row = fit_global2(cpart, "profile_hellinger", species,
                                  "All tested tissues", cohort, adj, age_scale,
                                  variant="per_cohort", depth_label="full",
                                  endpoint="profile_hellinger_centroid_mean")
                if row:
                    out_rows.append(row)
        # --- young-only slope (developmental drift inside reference window) ---
        if species == "Mouse":
            young = part[part["age"] <= 0.25].copy()
            row = young_only_slope(young, "profile_hellinger", species,
                                   "per month (1m vs 3m)", age_scale_months=1 / 12)
            if row:
                # age is stored in years, so the fitted slope is per year;
                # rescale to per month for the 1m-vs-3m contrast.
                for key in ("effect_sd", "se", "ci95_low", "ci95_high"):
                    row[key] = row[key] / 12.0
                out_rows.append(row)
        else:
            young = part[part["age"] <= 40].copy()
            row = young_only_slope(young, "profile_hellinger", species,
                                   "per year (20-40y)", age_scale_months=1.0)
            if row:
                out_rows.append(row)
    part_a_extended(out_rows)
    return out_rows


def part_a_extended(out_rows):
    """Donor-level nonlinear-depth fits for worm / fly / zebrafish, whose
    donor files are already in the repo (same schema as the cross-species
    donor file). No cell-level re-fetch is attempted for these species.
    Young references (from the scripts that built the donor files):
    worm day 1, fly 30 days, zebrafish 5 months."""
    specs = [
        ("results/extended-worm-vector-donors.csv", "C. elegans", 14.0,
         "young_day1", "Gao_GSE229022"),
        ("results/extended-afca-vector-donors.csv", "Fruit fly", 40.0,
         "young_30d", "AFCA"),
        ("results/extended-zebrafish-vector-donors.csv", "Zebrafish", 22.0,
         "young_5mo", "GSE275361"),
    ]
    adjustments = ["unadjusted", "linear_log_umi", "inv_sqrt_umi",
                   "linear_plus_inv_sqrt", "quadratic_log_umi", "spline_log_umi"]
    for path, species, age_scale, reference, cohort in specs:
        donors = pd.read_csv(REPO / path)
        donors["cohort"] = cohort
        for adj in adjustments:
            row = fit_global2(donors, "profile_hellinger", species,
                              "All tested tissues", "All cohorts", adj,
                              age_scale, variant="global", depth_label="full",
                              reference=reference,
                              endpoint="profile_hellinger_centroid_mean")
            if row:
                out_rows.append(row)
        for tissue, tpart in donors.groupby("tissue"):
            row = fit_global2(tpart, "profile_hellinger", species, tissue,
                              "All cohorts", "linear_log_umi", age_scale,
                              variant="per_tissue", depth_label="full",
                              reference=reference,
                              endpoint="profile_hellinger_centroid_mean")
            if row:
                out_rows.append(row)
        if species == "C. elegans" and "platform_nextseq" in donors.columns:
            row = fit_global2(donors, "profile_hellinger", species,
                              "All tested tissues", "All cohorts",
                              "linear_log_umi", age_scale, variant="global",
                              depth_label="full", reference=reference,
                              endpoint="profile_hellinger_centroid_mean",
                              extra={"platform_nextseq": "platform_nextseq"})
            if row:
                row["adjustment"] = "linear_log_umi+platform_nextseq"
                row["notes"] = ("Adds the NextSeq platform indicator "
                                "(centred within unit) to the linear depth model.")
                out_rows.append(row)
    return out_rows


# --------------------------------------------------------------------------
# Part B: mouse cell-level (Tabula Muris Senis)
# --------------------------------------------------------------------------

def hellinger_distance(matrix, reference_proportions):
    """Hellinger distance from each count row to a reference probability vector."""
    matrix = sp.csr_matrix(matrix, dtype=float)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError("Hellinger distance requires positive row sums")
    reference_proportions = np.asarray(reference_proportions, dtype=float)
    if reference_proportions.sum() <= 0:
        raise ValueError("Reference profile has no counts")
    reference_proportions = reference_proportions / reference_proportions.sum()
    root_matrix = matrix.copy()
    root_matrix.data = np.sqrt(root_matrix.data)
    affinity = np.asarray(root_matrix @ np.sqrt(reference_proportions)).ravel() / np.sqrt(totals)
    return np.sqrt(np.maximum(0.0, 1.0 - affinity))


def technical_gene_mask(gene_names):
    import re
    genes = pd.Index(gene_names).astype(str).str.upper()
    return np.array([
        bool(re.match(r"^(MT-|RP[LS]|MRP[LS])", gene))
        or gene in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for gene in genes
    ])


def build_reference(matrix, obs, cell_column, donor_column):
    stats = {}
    for cell_type in obs[cell_column].unique():
        ct_mask = (obs[cell_column] == cell_type).to_numpy()
        ct_matrix = matrix[ct_mask]
        ct_obs = obs.loc[ct_mask]
        per_donor = {}
        for donor in ct_obs[donor_column].unique():
            donor_mask = (ct_obs[donor_column] == donor).to_numpy()
            donor_matrix = ct_matrix[donor_mask]
            per_donor[str(donor)] = {
                "counts": np.asarray(donor_matrix.sum(axis=0)).ravel(),
                "n_cells": int(donor_matrix.shape[0]),
            }
        stats[str(cell_type)] = {
            "counts": np.asarray(ct_matrix.sum(axis=0)).ravel(),
            "n_cells": int(ct_matrix.shape[0]),
            "donors": per_donor,
        }
    return stats


def reference_without_donor(stats, donor):
    result = {"counts": stats["counts"].copy(), "n_cells": stats["n_cells"]}
    if donor in stats["donors"]:
        result["counts"] -= stats["donors"][donor]["counts"]
        result["n_cells"] -= stats["donors"][donor]["n_cells"]
    return result


def qc_tissue_path(tissue):
    return REPO / "data" / "interim" / f"m1-qc-cleaned-{tissue}.h5ad"


def rebuild_qc(force=False):
    """Rebuild the M1 QC-cleaned object from the Figshare raw file.

    Only M1 Step 3 is reproduced: per-tissue 500-5000 genes, >=1000 UMIs,
    per-tissue scrublet (expected doublet rate 0.06, seed 42). No PCA/UMAP/
    clustering is needed for the vector analysis. Verified against
    results/m1-qc-summary.csv (per-tissue final cell counts).

    Each cleaned tissue is written to its own file as it finishes
    (data/interim/m1-qc-cleaned-<tissue>.h5ad) so peak memory stays low
    and a later failure does not lose earlier tissues. A merged
    data/interim/m1-qc-cleaned.h5ad is assembled at the end when memory
    allows; the vector analysis reads the per-tissue files.
    """
    import gc
    import scanpy as sc
    raw_path = REPO / "data" / "raw" / "tms_droplet_raw.h5ad"
    out_path = REPO / "data" / "interim" / "m1-qc-cleaned.h5ad"
    if out_path.exists() and not force:
        print(f"QC-cleaned object already present: {out_path}", flush=True)
        return out_path
    if not raw_path.exists():
        raise FileNotFoundError(
            f"{raw_path} missing. Re-fetch the Figshare file per "
            "data/manifest/tms_droplet_raw.json and verify its sha256.")
    adata = sc.read_h5ad(raw_path, backed="r")
    print(f"Opened raw TMS droplet object (backed): {adata.shape}", flush=True)
    tissues = sorted(adata.obs["tissue"].astype(str).unique())
    decisions = []
    for tissue in tissues:
        tissue_path = qc_tissue_path(tissue)
        if tissue_path.exists() and not force:
            print(f"  {tissue}: cleaned file already present, skipping", flush=True)
            decisions.append({"tissue": tissue, "n_final_cells": int(sc.read_h5ad(tissue_path, backed="r").n_obs)})
            continue
        # per-tissue load keeps peak memory low (full object is ~250k cells)
        t_adata = adata[adata.obs["tissue"] == tissue].to_memory()
        for col in ("tissue", "age", "sex", "mouse.id", "cell_ontology_class"):
            t_adata.obs[col] = t_adata.obs[col].astype(str).str.strip()
        n_init = t_adata.n_obs
        t_adata.obs["n_counts"] = np.asarray(t_adata.X.sum(axis=1)).ravel()
        t_adata.obs["n_genes_by_counts"] = np.asarray((t_adata.X > 0).sum(axis=1)).ravel()
        keep = (t_adata.obs["n_genes_by_counts"] >= 500) & \
               (t_adata.obs["n_genes_by_counts"] <= 5000) & \
               (t_adata.obs["n_counts"] >= 1000)
        t_filtered = t_adata[keep].copy()
        del t_adata
        gc.collect()
        n_after_qc = t_filtered.n_obs
        batch_counts = t_filtered.obs["mouse.id"].value_counts()
        valid_batches = batch_counts[batch_counts >= 100].index
        try:
            if len(valid_batches) > 1 and len(valid_batches) == len(batch_counts):
                sc.pp.scrublet(t_filtered, batch_key="mouse.id",
                               expected_doublet_rate=0.06, random_state=42, verbose=False)
            else:
                sc.pp.scrublet(t_filtered, expected_doublet_rate=0.06,
                               random_state=42, verbose=False)
            is_doublet = t_filtered.obs["predicted_doublet"].fillna(False).to_numpy()
        except Exception as exc:  # noqa: BLE001 - keep the M1 fallback behaviour
            print(f"  [Scrublet notice for {tissue}]: {exc}", flush=True)
            is_doublet = np.zeros(n_after_qc, dtype=bool)
        t_cleaned = t_filtered[~is_doublet].copy()
        del t_filtered
        gc.collect()
        # Scrublet stores per-batch results under uns['scrublet']['batches']
        # keyed by mouse.id; one Lung mouse.id is '3-M-5/6', whose slash
        # is illegal as an HDF5 group key in current anndata/h5py. The
        # doublet calls themselves live in obs['predicted_doublet'] (kept);
        # the uns block is not needed downstream, so drop it before writing.
        t_cleaned.uns.pop("scrublet", None)
        tmp = tissue_path.with_suffix(".tmp.h5ad")
        t_cleaned.write(tmp)
        tmp.replace(tissue_path)
        decisions.append({"tissue": tissue, "n_initial": n_init,
                          "n_after_gene_umi_qc": n_after_qc,
                          "n_doublets_removed": int(is_doublet.sum()),
                          "n_final_cells": int(t_cleaned.n_obs)})
        print(f"  {tissue}: {n_init} -> {n_after_qc} (QC) -> {t_cleaned.n_obs} (singlets)", flush=True)
        del t_cleaned
        gc.collect()
    expected = pd.read_csv(RESULTS / "m1-qc-summary.csv").set_index("tissue")["n_final_cells"]
    for row in decisions:
        exp = expected.get(row["tissue"])
        if exp is None or "n_final_cells" not in row:
            continue
        flag = "OK" if exp == row["n_final_cells"] else f"MISMATCH (expected {exp})"
        print(f"  verify {row['tissue']}: {row['n_final_cells']} {flag}", flush=True)
    # No merged file is assembled here: concatenating all 16 tissues
    # materialises ~2.8 GB at once and was OOM-killed twice on this
    # shared VM. The per-tissue files are the working input for the
    # vector analysis (load_mouse_tissue); a merged m1-qc-cleaned.h5ad
    # can be assembled on a bigger machine with `ad.concat` over them.
    return out_path


def mouse_tissue_names():
    """Tissue list for the mouse vector analysis, from per-tissue QC files
    when present, else from the merged object."""
    import glob
    files = sorted(glob.glob(str(REPO / "data" / "interim" / "m1-qc-cleaned-*.h5ad")))
    names = [Path(f).name[len("m1-qc-cleaned-"):-len(".h5ad")] for f in files]
    return [name for name in names if name != ""]


def load_mouse_tissue(tissue):
    """Load one QC-cleaned tissue (singlets only) into memory."""
    import anndata as ad
    tissue_path = qc_tissue_path(tissue)
    if tissue_path.exists():
        part = ad.read_h5ad(tissue_path)
    else:
        merged = REPO / "data" / "interim" / "m1-qc-cleaned.h5ad"
        import scanpy as sc
        adata = sc.read_h5ad(merged, backed="r")
        mask = ((adata.obs["tissue"] == tissue)
                & ~adata.obs["predicted_doublet"]
                & (adata.obs["cell_ontology_class"].astype(str) != "unknown"))
        part = adata[mask].to_memory()
    part.obs["age"] = part.obs["age"].astype(str).map(MOUSE_AGES) / 12.0
    part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
    part.obs["mouse.id"] = part.obs["mouse.id"].astype(str)
    part.X = part.X.tocsr()
    return part


def mouse_reference_mask(obs_age_years, reference):
    if reference == "young_1_3m":
        return obs_age_years <= 0.25
    if reference == "1m_only":
        return np.isclose(obs_age_years, 1.0 / 12.0)
    if reference == "3m_only":
        return np.isclose(obs_age_years, 0.25)
    if reference == "18m":
        return np.isclose(obs_age_years, 1.5)
    raise ValueError(f"unknown reference {reference}")


def mouse_tissue_endpoints(X, obs, var_names, tissue_name, depth, reference):
    """Donor-unit endpoints for one tissue.

    X: raw-count CSR matrix (cells x genes); obs: DataFrame with age
    (years, float), sex_male, mouse.id, cell_ontology_class. depth None =
    full depth with the original >=1500 UMI inclusion; otherwise fixed-UMI
    downsampling (same sc.pp.downsample_counts, seed 42, as the Phase 2
    depth sweep). Takes plain arrays, not AnnData, so the caller can free
    the tissue object before computing.
    """
    import scanpy as sc
    X = sp.csr_matrix(X)
    totals_actual = np.asarray(X.sum(axis=1)).ravel()
    obs = obs.copy()
    obs["actual_n_counts"] = totals_actual
    obs["mouse.id"] = obs["mouse.id"].astype(str)
    if depth is None:
        keep_cells = totals_actual >= 1500
        matrix = X[keep_cells][:, ~technical_gene_mask(var_names)].tocsr()
        obs = obs.loc[keep_cells]
    else:
        keep_cells = totals_actual >= depth
        X = X[keep_cells].tocsr()
        import anndata as ad
        selected = ad.AnnData(X=X)
        sc.pp.downsample_counts(selected, counts_per_cell=depth, replace=False,
                                random_state=RNG_SEED)
        matrix = selected.X.tocsr()
        assert np.all(np.asarray(matrix.sum(axis=1)).ravel() == depth)
        matrix = matrix[:, ~technical_gene_mask(var_names)].tocsr()
        obs = obs.loc[keep_cells]
    ref_mask = mouse_reference_mask(obs["age"].to_numpy(dtype=float), reference)
    ref_obs = obs.loc[ref_mask]
    reference_stats = build_reference(matrix[ref_mask], ref_obs,
                                      "cell_ontology_class", "mouse.id")
    rows = []
    for cell_type in obs["cell_ontology_class"].unique():
        cell_type = str(cell_type)
        if cell_type == "unknown" or cell_type not in reference_stats:
            continue
        ct_mask = (obs["cell_ontology_class"].astype(str) == cell_type).to_numpy()
        ct_obs = obs.loc[ct_mask]
        ct_matrix = matrix[ct_mask]
        # actual (pre-downsampling) UMI totals for this cell-type block
        ct_actual = ct_obs["actual_n_counts"].to_numpy(dtype=float)
        for donor in ct_obs["mouse.id"].unique():
            donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
            if donor_mask.sum() < 20:
                continue
            donor_obs = ct_obs.loc[donor_mask]
            donor_matrix = ct_matrix[donor_mask]
            ref = reference_without_donor(reference_stats[cell_type], str(donor))
            if ref["n_cells"] < 50 or ref["counts"].sum() <= 0:
                continue
            ref_props = ref["counts"] / ref["counts"].sum()
            dists = hellinger_distance(donor_matrix, ref_props)
            # pseudobulk: one distance from the donor's summed counts
            pseudo_counts = np.asarray(donor_matrix.sum(axis=0)).ravel()
            pseudo = float(hellinger_distance(
                sp.csr_matrix(pseudo_counts.reshape(1, -1)), ref_props)[0])
            # dispersion to donor's own centroid (mean profile)
            own_props = pseudo_counts / pseudo_counts.sum()
            disp_centroid = float(hellinger_distance(donor_matrix, own_props).mean())
            # dispersion to donor's own median profile (gene-wise median cell)
            dense = donor_matrix.toarray() if donor_matrix.shape[0] <= 400 else \
                donor_matrix[np.random.RandomState(RNG_SEED).choice(
                    donor_matrix.shape[0], size=400, replace=False)].toarray()
            median_profile = np.median(dense, axis=0)
            if median_profile.sum() > 0:
                disp_median = float(hellinger_distance(
                    donor_matrix, median_profile / median_profile.sum()).mean())
            else:
                disp_median = np.nan
            donor_totals = ct_actual[donor_mask]
            rows.append({
                "species": "Mouse", "tissue": tissue_name,
                "cohort": "TMS", "cell_type": cell_type, "donor": str(donor),
                "age": float(donor_obs["age"].iloc[0]),
                "sex_male": float(donor_obs["sex_male"].iloc[0]),
                "n_cells": int(donor_mask.sum()),
                "mean_umi": float(donor_totals.mean()),
                "mean_log_umi": float(np.log(donor_totals).mean()),
                "profile_hellinger_centroid_mean": float(dists.mean()),
                "profile_hellinger_centroid_median": float(np.median(dists)),
                "dispersion_to_median": disp_median,
                "dispersion_to_centroid": disp_centroid,
                "pseudobulk_hellinger": pseudo,
            })
    return rows


def load_tissue_arrays(tissue):
    """Load one cleaned tissue and return (X, obs, var_names), freeing
    the AnnData immediately so at most one copy of the tissue matrix is
    alive during the endpoint computation (OOM guard on this VM)."""
    part = load_mouse_tissue(tissue)
    part = part[part.obs["cell_ontology_class"].astype(str) != "unknown"].copy()
    X = part.X.tocsr()
    obs = part.obs.copy()
    var_names = part.var_names.copy()
    del part
    gc.collect()
    return X, obs, var_names


def run_mouse_cell_level(out_rows, donor_tables, arms=None):
    qc_path = rebuild_qc()
    print(f"QC source: {qc_path}", flush=True)
    tissues = mouse_tissue_names()
    print(f"Mouse tissues: {tissues}", flush=True)
    endpoints = ["profile_hellinger_centroid_mean", "profile_hellinger_centroid_median",
                 "dispersion_to_median", "dispersion_to_centroid", "pseudobulk_hellinger"]
    adjustments = ["unadjusted", "linear_log_umi", "inv_sqrt_umi", "spline_log_umi"]
    arms = set(arms) if arms else {"full", "500", "1500", "3000", "refs", "controls"}
    # --- fixed-UMI depths x young reference (checkpointed per arm) ---
    for depth_label, depth in (("full", None), ("500", 500), ("1500", 1500), ("3000", 3000)):
        if depth_label not in arms:
            continue
        all_rows = []
        for tissue in tissues:
            X, obs, var_names = load_tissue_arrays(tissue)
            rows = mouse_tissue_endpoints(X, obs, var_names, tissue, depth, "young_1_3m")
            all_rows.extend(rows)
            print(f"Mouse depth {depth_label} {tissue}: {len(rows)} donor units", flush=True)
            del X, obs, var_names
            gc.collect()
        donors = pd.DataFrame(all_rows)
        donor_tables[f"mouse_depth_{depth_label}"] = donors
        for endpoint in endpoints:
            for adj in adjustments:
                row = fit_global2(donors, endpoint, "Mouse", "All tested tissues",
                                  "All cohorts", adj, 0.4, variant="global",
                                  depth_label=depth_label, reference="young_1_3m",
                                  endpoint=endpoint, source="cell_level")
                if row:
                    out_rows.append(row)
        # per-organ decomposition at this depth (linear depth adjustment)
        for tissue, tpart in donors.groupby("tissue"):
            row = fit_global2(tpart, "profile_hellinger_centroid_mean", "Mouse",
                              tissue, "All cohorts", "linear_log_umi", 0.4,
                              variant="per_tissue", depth_label=depth_label,
                              reference="young_1_3m",
                              endpoint="profile_hellinger_centroid_mean", source="cell_level")
            if row:
                out_rows.append(row)
        save_results(out_rows, donor_tables)
    # --- alternative references at full depth ---
    if "refs" in arms:
        for reference in ("1m_only", "3m_only", "18m"):
            all_rows = []
            for tissue in tissues:
                X, obs, var_names = load_tissue_arrays(tissue)
                rows = mouse_tissue_endpoints(X, obs, var_names, tissue, None, reference)
                all_rows.extend(rows)
                print(f"Mouse ref {reference} {tissue}: {len(rows)} donor units", flush=True)
                del X, obs, var_names
                gc.collect()
            donors = pd.DataFrame(all_rows)
            donor_tables[f"mouse_ref_{reference}"] = donors
            for endpoint in ("profile_hellinger_centroid_mean", "pseudobulk_hellinger"):
                for adj in ("unadjusted", "linear_log_umi", "spline_log_umi"):
                    row = fit_global2(donors, endpoint, "Mouse", "All tested tissues",
                                      "All cohorts", adj, 0.4, variant="global",
                                      depth_label="full", reference=reference,
                                      endpoint=endpoint, source="cell_level")
                    if row:
                        row["notes"] = ("Reference is " +
                                        {"1m_only": "1-month mice only",
                                         "3m_only": "3-month mice only",
                                         "18m": "mature-adult 18-month mice only"}[reference] +
                                        " (cross-fitted by donor).")
                        out_rows.append(row)
            save_results(out_rows, donor_tables)
    # --- young-vs-young negative controls (split-half + 1m-reference) ---
    if "controls" in arms:
        young_controls(out_rows, donor_tables)
        save_results(out_rows, donor_tables)
    return out_rows


def _young_donor_split(donor_ids_by_age):
    """Deterministic stratified split of young donors into halves A/B."""
    group_a, group_b = [], []
    for age in sorted(donor_ids_by_age):
        donors = sorted(donor_ids_by_age[age])
        for i, donor in enumerate(donors):
            (group_a if i % 2 == 0 else group_b).append(donor)
    return group_a, group_b


def young_controls(out_rows, donor_tables):
    """Split-half young reference and 1m-only reference controls.

    Statistic: mean standardised (within tissue-cell-type) distance of the
    held-out young group minus that of the reference group (LODO). Null
    distribution by permuting donor labels within age strata (split-half)
    or across the 9 young donors (1m vs 3m). 95% CI by donor bootstrap.
    """
    rng = np.random.RandomState(RNG_SEED)
    tissues = mouse_tissue_names()
    # collect young donor ids per age group across tissues
    young_by_age = {}
    for tissue in tissues:
        part = load_mouse_tissue(tissue)
        part = part[part.obs["cell_ontology_class"].astype(str) != "unknown"]
        young = part[part.obs["age"] <= 0.25]
        for age, sub in young.obs.groupby("age"):
            young_by_age.setdefault(float(age), set()).update(
                sub["mouse.id"].astype(str).unique().tolist())
        del part
    young_by_age = {age: sorted(donors) for age, donors in young_by_age.items()}
    print(f"Young donors by age (years): {young_by_age}", flush=True)

    def distances_for_reference(ref_donor_set, test_donor_set, ref_age_mask):
        """Mean per-cell Hellinger per donor-unit; reference donors are
        leave-one-donor-out, test donors use the full reference."""
        rows = []
        for tissue in tissues:
            X, obs, var_names = load_tissue_arrays(tissue)
            totals = np.asarray(X.sum(axis=1)).ravel()
            keep = totals >= 1500
            X = X[keep].tocsr()
            obs = obs.loc[keep]
            keep_genes = ~technical_gene_mask(var_names)
            matrix = X[:, keep_genes].tocsr()
            del X
            gc.collect()
            ref_mask = ref_age_mask(obs["age"].to_numpy(dtype=float)) & \
                obs["mouse.id"].astype(str).isin(ref_donor_set).to_numpy()
            if ref_mask.sum() == 0:
                continue
            ref_stats = build_reference(matrix[ref_mask], obs.loc[ref_mask],
                                        "cell_ontology_class", "mouse.id")
            wanted = set(ref_donor_set) | set(test_donor_set)
            for cell_type in obs["cell_ontology_class"].unique():
                cell_type = str(cell_type)
                if cell_type == "unknown" or cell_type not in ref_stats:
                    continue
                ct_mask = (obs["cell_ontology_class"].astype(str) == cell_type).to_numpy()
                ct_obs = obs.loc[ct_mask]
                ct_matrix = matrix[ct_mask]
                for donor in ct_obs["mouse.id"].unique():
                    donor = str(donor)
                    if donor not in wanted:
                        continue
                    donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
                    if donor_mask.sum() < 20:
                        continue
                    ref = reference_without_donor(ref_stats[cell_type], donor) \
                        if donor in ref_donor_set else ref_stats[cell_type]
                    if ref["n_cells"] < 50 or ref["counts"].sum() <= 0:
                        continue
                    dists = hellinger_distance(ct_matrix[donor_mask],
                                               ref["counts"] / ref["counts"].sum())
                    rows.append({"tissue": tissue, "cell_type": cell_type,
                                 "donor": donor, "group": "ref" if donor in ref_donor_set else "test",
                                 "age": float(ct_obs.loc[donor_mask, "age"].iloc[0]),
                                 "distance": float(dists.mean())})
        return pd.DataFrame(rows)

    def control_statistic(frame):
        """Standardised test-minus-ref mean difference, pooled over units."""
        if frame.empty:
            return np.nan
        diffs = []
        for _, part in frame.groupby(["tissue", "cell_type"]):
            if part["group"].nunique() < 2:
                continue
            scale = part["distance"].std(ddof=1)
            if not np.isfinite(scale) or scale <= 0:
                continue
            ref_mean = part.loc[part["group"] == "ref", "distance"].mean()
            test_mean = part.loc[part["group"] == "test", "distance"].mean()
            diffs.append((test_mean - ref_mean) / scale)
        return float(np.mean(diffs)) if diffs else np.nan

    def run_control(name, ref_donors, test_donors, ref_age_mask, permute_labels):
        frame = distances_for_reference(set(ref_donors), set(test_donors), ref_age_mask)
        donor_tables[f"mouse_young_control_{name}"] = frame
        observed = control_statistic(frame)
        # permutation null: re-split the same donor pool
        pool = sorted(set(ref_donors) | set(test_donors))
        perm_stats = []
        for _ in range(2000):
            shuffled = rng.permutation(pool)
            fake_ref = set(shuffled[:len(ref_donors)])
            fake_frame = frame.copy()
            fake_frame["group"] = np.where(
                fake_frame["donor"].isin(fake_ref), "ref", "test")
            perm_stats.append(control_statistic(fake_frame))
        perm_stats = np.array([s for s in perm_stats if np.isfinite(s)])
        pvalue = float(np.mean(np.abs(perm_stats) >= abs(observed) - 1e-12)) \
            if len(perm_stats) else np.nan
        # bootstrap CI over donor-units within each group (row resampling)
        boot = []
        ref_rows = frame[frame["group"] == "ref"]
        test_rows = frame[frame["group"] == "test"]
        for _ in range(1000):
            sampled = pd.concat([
                ref_rows.sample(n=len(ref_rows), replace=True,
                                random_state=int(rng.randint(0, 2 ** 31 - 1))),
                test_rows.sample(n=len(test_rows), replace=True,
                                 random_state=int(rng.randint(0, 2 ** 31 - 1))),
            ], ignore_index=True)
            boot.append(control_statistic(sampled))
        boot = np.array([s for s in boot if np.isfinite(s)])
        out_rows.append({
            "endpoint": "young_vs_young_control", "species": "Mouse",
            "tissue": "All tested tissues", "cohort": "TMS",
            "variant": name, "depth": "full", "reference": "split_young",
            "adjustment": "none", "donor_set": "all_donor_ids",
            "age_scale_years": np.nan,
            "n_donors": int(frame["donor"].nunique()),
            "n_units": int(frame.groupby(["tissue", "cell_type"]).ngroups),
            "n_donor_units": int(len(frame)),
            "effect_sd": observed,
            "se": float(np.std(boot, ddof=1)) if len(boot) > 1 else np.nan,
            "ci95_low": float(np.percentile(boot, 2.5)) if len(boot) else np.nan,
            "ci95_high": float(np.percentile(boot, 97.5)) if len(boot) else np.nan,
            "pvalue": pvalue,
            "source": "cell_level",
            "notes": ("Permutation p over 2000 donor-label reshuffles (the "
                      "inference); CI is a row bootstrap over donor-units "
                      "and ignores donor clustering, so it is anti-conservative. "
                      "A value near 0 means no spurious young-vs-young displacement."),
        })
        print(f"Young control {name}: stat={observed:.4f} p={pvalue:.4f} "
              f"n_donor_units={len(frame)}", flush=True)

    # (a) split-half young reference
    group_a, group_b = _young_donor_split(young_by_age)
    run_control("split_half_young_ref", group_a, group_b,
                lambda age: age <= 0.25, permute_labels=True)
    # (b) 1m-only reference: 3m donors (test) vs 1m donors (ref, LODO)
    donors_1m = young_by_age.get(1.0 / 12.0, [])
    donors_3m = young_by_age.get(0.25, [])
    run_control("ref_1m_test_3m", donors_1m, donors_3m,
                lambda age: np.isclose(age, 1.0 / 12.0), permute_labels=True)


# --------------------------------------------------------------------------
# Part C: naked mole-rat skin (GSE152886)
# --------------------------------------------------------------------------

def nmr_annotated_data_light():
    """Memory-light re-implementation of nmr_annotated_data() from
    scripts/run_cross_species_expression_vector.py (same QC thresholds,
    HVG/PCA/neighbors/Leiden parameters, marker compartments, and
    technical-gene removal), with two deviations forced by this VM:
    (1) the raw TSV is loaded in 1000-barcode column chunks into CSR
    instead of one dense pandas DataFrame (the dense load peaks above
    4 GB and is OOM-killed here; the values are identical, int64);
    (2) the Ensembl-to-symbol TSV is read from data/reference/ (the
    source the interim copy is made from) because the interim copy
    data/interim/nmr-ensg2symbol.tsv is absent in this checkout."""
    import scanpy as sc
    import sys as _sys
    _sys.path.insert(0, str(REPO))
    from scripts.run_cross_species_expression_vector import (
        NMR_AGES, NMR_GROUPS, technical_gene_mask)
    raw_path = REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz"
    # barcode (column) positions, dropping sample 2V before loading.
    # NOTE: the file's header tokens carry literal single quotes, and
    # integer-position usecols misaligns against index_col=0 here (it
    # silently returns a trailing all-NaN column), so select by name
    # with the quotes stripped; the returned columns keep the quoted
    # names, exactly as in the original full read.
    header = pd.read_csv(raw_path, sep="\t", index_col=0, nrows=0)
    barcodes = [str(c) for c in header.columns]
    keep_pos = [i for i, b in enumerate(barcodes) if b.split("_")[0] != "2V"]
    gene_ids = None
    blocks = []
    for start in range(0, len(keep_pos), 1000):
        pos = keep_pos[start:start + 1000]
        cols = [barcodes[p].strip("'") for p in pos]
        chunk = pd.read_csv(raw_path, sep="\t", index_col=0, usecols=cols)
        # A usecols read returns the header tokens with their literal
        # single quotes kept, while a plain read strips them; normalise
        # to the plain-read (unquoted) names so sample parsing matches
        # the original loader exactly.
        chunk.columns = [str(c).strip("'") for c in chunk.columns]
        if gene_ids is None:
            gene_ids = chunk.index.astype(str)
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.int64)).T)
        del chunk
        gc.collect()
    X = sp.vstack(blocks, format="csr")
    del blocks
    gc.collect()
    obs = pd.DataFrame({"barcode": [barcodes[i] for i in keep_pos]})
    obs["sample"] = obs["barcode"].str.split("_").str[0]
    obs["age"] = obs["sample"].map(NMR_AGES)
    obs["group"] = obs["sample"].map(NMR_GROUPS)
    adata = sc.AnnData(X=X, obs=obs.set_index("barcode"))
    adata.var_names = gene_ids
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata = adata[(adata.obs["n_counts"] >= 2000) & (adata.obs["n_genes"] >= 700)].copy()
    raw_counts = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    markers = {
        "immune": ["ENSHGLG00000017278", "ENSHGLG00000015756", "ENSHGLG00000005722"],
        "basal": ["ENSHGLG00000000726", "ENSHGLG00000016277", "ENSHGLG00000017682"],
        "corneum": ["ENSHGLG00000010577", "ENSHGLG00000001162"],
    }
    marker_genes = {gene for genes in markers.values() for gene in genes}
    retain = adata.var["highly_variable"] | adata.var_names.isin(marker_genes)
    adata = adata[:, retain].copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, random_state=42)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
    sc.tl.leiden(adata, resolution=0.4, random_state=42)
    annotations = []
    for cluster in sorted(adata.obs["leiden"].unique()):
        cells = adata[adata.obs["leiden"] == cluster]
        scores = {}
        for name, genes in markers.items():
            present = [gene for gene in genes if gene in cells.var_names]
            scores[name] = float(np.asarray(cells[:, present].X.mean()).ravel()[0]) if present else 0.0
        if scores["immune"] > 0.1 and scores["immune"] >= max(scores["basal"], scores["corneum"]):
            compartment = "immune"
        elif scores["corneum"] > scores["basal"]:
            compartment = "corneum"
        elif scores["basal"] > 0:
            compartment = "basal"
        else:
            compartment = "intermediate"
        annotations.append({"cluster": cluster, "compartment": compartment})
    annotation_map = pd.DataFrame(annotations).set_index("cluster")["compartment"]
    adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(annotation_map)
    adata.X = raw_counts[:, retain.to_numpy()]
    symbol_path = REPO / "data" / "interim" / "nmr-ensg2symbol.tsv"
    if not symbol_path.exists():
        symbol_path = REPO / "data" / "reference" / "nmr-ensg2symbol.tsv"
    symbols = pd.read_csv(
        symbol_path, sep="\t", skiprows=1,
        names=["ensembl_id", "symbol"],
    ).set_index("ensembl_id")["symbol"]
    gene_symbols = adata.var_names.map(symbols).fillna("").astype(str)
    keep = ~technical_gene_mask(gene_symbols)
    return adata[:, keep].copy()


def run_nmr(out_rows, donor_tables):
    import sys as _sys
    _sys.path.insert(0, str(REPO))
    from scripts.run_cross_species_expression_vector import NMR_AGES, NMR_GROUPS
    import anndata as ad
    raw_path = REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz"
    if not raw_path.exists():
        print(f"NMR raw file missing ({raw_path}); skipping Part C.", flush=True)
        return out_rows
    adata = nmr_annotated_data_light()
    adata.obs["sex_male"] = 0.0
    adata.obs["age"] = adata.obs["sample"].map(NMR_AGES)
    adata.obs["group"] = adata.obs["sample"].map(NMR_GROUPS)
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    print(f"NMR: {adata.n_obs} cells, median UMI {np.median(totals):.0f}, "
          f"fraction >=1500 UMIs: {(totals >= 1500).mean():.2f}", flush=True)
    import scanpy as sc
    for depth in (None, 500):
        depth_label = "full" if depth is None else str(depth)
        if depth is None:
            selected = adata.copy()
            matrix = selected.X.tocsr()
        else:
            selected = adata[totals >= depth].copy()
            sc.pp.downsample_counts(selected, counts_per_cell=depth, replace=False,
                                    random_state=RNG_SEED)
            matrix = selected.X.tocsr()
        obs = selected.obs.copy()
        # nmr_annotated_data_light already removed technical genes (MT/RP/sex genes
        # on gene symbols) before returning, so no second filter is applied.
        matrix = matrix
        young_mask = (obs["group"] == "young").to_numpy()
        ref_stats = build_reference(matrix[young_mask], obs.loc[young_mask],
                                    "cell_type", "sample")
        rows = []
        for cell_type in obs["cell_type"].unique():
            ct_mask = (obs["cell_type"] == cell_type).to_numpy()
            ct_obs = obs.loc[ct_mask]
            ct_matrix = matrix[ct_mask]
            for donor in ct_obs["sample"].unique():
                donor_mask = (ct_obs["sample"] == donor).to_numpy()
                if donor_mask.sum() < 20:
                    continue
                ref = reference_without_donor(ref_stats[str(cell_type)], str(donor))
                if ref["n_cells"] < 50 or ref["counts"].sum() <= 0:
                    continue
                dists = hellinger_distance(ct_matrix[donor_mask],
                                           ref["counts"] / ref["counts"].sum())
                pseudo = float(hellinger_distance(
                    sp.csr_matrix(np.asarray(ct_matrix[donor_mask].sum(axis=0)).ravel().reshape(1, -1)),
                    ref["counts"] / ref["counts"].sum())[0])
                rows.append({"species": "Naked mole-rat", "tissue": "Skin",
                             "cohort": "GSE152886", "cell_type": str(cell_type),
                             "donor": str(donor),
                             "age": float(ct_obs.loc[donor_mask, "age"].iloc[0]),
                             "sex_male": 0.0, "n_cells": int(donor_mask.sum()),
                             "mean_umi": float(np.asarray(
                                 ct_matrix[donor_mask].sum(axis=1)).ravel().mean()),
                             "mean_log_umi": float(np.log(np.asarray(
                                 ct_matrix[donor_mask].sum(axis=1)).ravel()).mean()),
                             "profile_hellinger_centroid_mean": float(dists.mean()),
                             "pseudobulk_hellinger": pseudo})
        donors = pd.DataFrame(rows)
        donor_tables[f"nmr_depth_{depth_label}"] = donors
        # exact permutation test: 2 aged of 5 samples (as in nmr_tests)
        import itertools
        for endpoint in ("profile_hellinger_centroid_mean", "pseudobulk_hellinger"):
            scales, effects = {}, []
            for cell_type, part in donors.groupby("cell_type"):
                if len(part) != 5:
                    continue
                aged = part["age"] > 5
                pooled_sd = np.sqrt((part.loc[aged, endpoint].var(ddof=1)
                                     + part.loc[~aged, endpoint].var(ddof=1)) / 2)
                scales[cell_type] = pooled_sd
                effects.append((part.loc[aged, endpoint].mean()
                                - part.loc[~aged, endpoint].mean()) / pooled_sd)
            if not effects:
                continue
            observed = float(np.mean(effects))
            sample_ids = sorted(donors["donor"].unique())
            permuted = []
            for indices in itertools.combinations(range(5), 2):
                aged_samples = {sample_ids[i] for i in indices}
                effs = []
                for cell_type in scales:
                    ct = donors[donors["cell_type"] == cell_type]
                    aged = ct["donor"].isin(aged_samples)
                    effs.append((ct.loc[aged, endpoint].mean()
                                 - ct.loc[~aged, endpoint].mean()) / scales[cell_type])
                permuted.append(np.mean(effs))
            pvalue = float(np.mean(np.abs(permuted) >= abs(observed) - 1e-12))
            out_rows.append({
                "endpoint": endpoint, "species": "Naked mole-rat",
                "tissue": "Skin", "cohort": "GSE152886",
                "variant": "global", "depth": depth_label,
                "reference": "young", "adjustment": "none",
                "donor_set": "all_donor_ids", "age_scale_years": 3.1,
                "n_donors": 5, "n_units": len(scales),
                "n_donor_units": int(len(donors)),
                "effect_sd": observed * 3.1 / 10.2, "se": np.nan,
                "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": pvalue,
                "source": "cell_level",
                "notes": "Exact permutation over 2-aged-of-5 sample splits; "
                         "effect rescaled as in the original nmr_tests.",
            })
    return out_rows


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--donor-only", action="store_true")
    parser.add_argument("--cell-level", action="store_true")
    parser.add_argument("--nmr", action="store_true")
    parser.add_argument("--qc-only", action="store_true")
    parser.add_argument("--rebuild-qc", action="store_true")
    parser.add_argument("--arms", default=None,
                        help="Comma-separated subset of Part B arms to run: "
                             "full,500,1500,3000,refs,controls (default: all).")
    args = parser.parse_args()
    if args.qc_only:
        rebuild_qc(force=args.rebuild_qc)
        return
    run_all = not (args.donor_only or args.cell_level or args.nmr)

    out_rows = []
    donor_tables = {}

    print("=== Part A: donor-level ===", flush=True)
    part_a_donor_level(out_rows)
    save_results(out_rows, donor_tables)

    if run_all or args.cell_level:
        print("=== Part B: mouse cell-level ===", flush=True)
        if args.rebuild_qc:
            rebuild_qc(force=True)
        arms = args.arms.split(",") if args.arms else None
        run_mouse_cell_level(out_rows, donor_tables, arms=arms)
        save_results(out_rows, donor_tables)

    if run_all or args.nmr:
        print("=== Part C: NMR cell-level ===", flush=True)
        run_nmr(out_rows, donor_tables)
        save_results(out_rows, donor_tables)

    print(results_summary(out_rows), flush=True)


def results_summary(out_rows):
    return pd.DataFrame(out_rows).to_string(index=False) if out_rows else "(no rows)"


def save_results(out_rows, donor_tables):
    """Checkpoint the output CSVs (also called at the end of each part
    and each Part B arm, so a later failure does not lose earlier work).

    Upserts by the row key, because arms run as separate processes:
    each process re-runs Part A and its selected arms, so a plain
    overwrite would delete other arms' rows."""
    cols = ["endpoint", "species", "tissue", "cohort", "variant", "depth",
            "reference", "adjustment", "donor_set", "age_scale_years",
            "n_donors", "n_units", "n_donor_units", "effect_sd", "se",
            "ci95_low", "ci95_high", "pvalue", "source", "notes"]
    results = pd.DataFrame(out_rows)
    for col in cols:
        if col not in results.columns:
            results[col] = np.nan
    results = results[cols]
    key_cols = ["endpoint", "species", "tissue", "cohort", "variant",
                "depth", "reference", "adjustment", "donor_set", "source"]
    out_path = RESULTS / "reviewer-rev2-hellinger.csv"
    if out_path.exists():
        old = pd.read_csv(out_path, dtype=str)
        if "source" not in old.columns:
            # one-time migration for CSVs written before source existed:
            # downsampled-depth rows came from the cell-level recompute.
            old["source"] = np.where(
                old["depth"].astype(str).isin(["500", "1500", "3000"]),
                "cell_level", "donor_file")
        old = old[~old.set_index(key_cols).index.isin(
            results.set_index(key_cols).index)]
        combined = pd.concat([old, results.astype(str)], ignore_index=True)
    else:
        combined = results.astype(str)
    combined.to_csv(out_path, index=False)
    print(f"Checkpoint {out_path}: {len(combined)} rows", flush=True)
    if donor_tables:
        frames = []
        for name, frame in donor_tables.items():
            frame = frame.copy()
            frame["table"] = name
            frames.append(frame)
        new = pd.concat(frames, ignore_index=True, sort=False)
        donors_path = RESULTS / "reviewer-rev2-hellinger-donors.csv"
        if donors_path.exists():
            old_d = pd.read_csv(donors_path)
            old_d = old_d[~old_d["table"].isin(new["table"].unique())]
            new = pd.concat([old_d, new], ignore_index=True, sort=False)
        new.to_csv(donors_path, index=False)
        print(f"Checkpoint {donors_path}: {len(new)} rows", flush=True)


if __name__ == "__main__":
    main()
