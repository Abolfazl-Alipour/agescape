#!/usr/bin/env python3
"""Rev2 fix 10 (h5ad route): matched random-signature null in MSSM.

Same analysis as scripts/run_rev2_transport_null.py, but reads expression
from the Census source h5ad for the MSSM dataset instead of via SOMA.
This route exists because the SOMA path proved unusable on this VM:
`cellxgene_census.get_anndata` balloons past 2.7 GB RSS for a 193-cell
slice (OOM-killed three times) and S3 reads through the egress proxy
stall below TileDB's low-speed timeout, while a plain curl download
streams fine.

Data: <scratch>/mssm.h5ad (Census 2025-11-08 source h5ad for dataset
37a17b78-4864-4a42-b67b-31c00962795a, 36,092,670,231 bytes).
Cell selection mirrors M9 exactly: primary data, donors with at least one
normal cell, >= 20 cells per donor-cell-type, all genes in var.
Donor pseudobulk = mean over cells of log1p(CPM) per gene; library per
cell is the X row sum of raw counts (== obs raw_sum; verified on a sample
before the null runs). The null itself (200 expression-decile-matched
random signatures, seed 20261105, M9 age-plus-sex OLS, centred two-sided
empirical p, BH over the 21 MSSM cell types) is imported verbatim from
run_rev2_transport_null so the two routes cannot diverge.

Writes results/reviewer-rev2-transport-null-mssm.csv and caches per-cell-type
pseudobulk .npz files in the same scratch dir / format as the SOMA route,
so the two routes share intermediates and the null phase.

Usage:
    python scripts/run_rev2_transport_null_h5ad.py [--h5ad PATH]
        [--scratch-dir DIR] [--rebuild-pseudobulk]
"""

from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.run_rev2_transport_null import (  # noqa: E402
    COHORT,
    ITERATIONS,
    MIN_CELLS_PER_DONOR_CT,
    RESULTS,
    SEED,
    fit_age_sex,
    load_program_tiers,
    random_signature_scores,
    save_npz_retry,
)

DEFAULT_SCRATCH = Path.home() / "workspace" / "rev2-mssm-null-scratch"
DEFAULT_H5AD = DEFAULT_SCRATCH / "mssm.h5ad"
CELL_CHUNK = 2000


def build_pseudobulk_h5ad(adata, obs_index, cell_type, donors, scratch: Path) -> Path:
    """Donor pseudobulk (mean logCPM per gene) for one MSSM cell type.

    obs_index: DataFrame with one row per needed cell: columns
    row (h5ad row index), donor_id, raw_sum. Already restricted to this
    cell type and to donors with >= 20 cells. Reads X in row chunks.
    """
    out = scratch / f"mssm_pseudobulk_{cell_type}.npz".replace(" ", "_").replace("/", "_")
    if out.exists():
        return out
    genes = np.array([str(g) for g in adata.var["feature_name"]]) if "feature_name" in adata.var.columns \
        else np.array([str(g) for g in adata.var.index])
    n_genes = len(genes)
    sums = np.zeros((len(donors), n_genes), dtype=np.float64)
    counts = np.zeros(len(donors), dtype=np.int64)
    library = np.zeros(len(donors), dtype=np.float64)
    donor_pos = {d: i for i, d in enumerate(donors)}
    sub = obs_index.sort_values("row")
    rows = sub["row"].to_numpy()
    dpos = np.array([donor_pos[d] for d in sub["donor_id"]])
    lib = sub["raw_sum"].to_numpy(float)
    for start in range(0, len(rows), CELL_CHUNK):
        chunk_rows = rows[start:start + CELL_CHUNK]
        X = sp.csr_matrix(adata.X[chunk_rows]).astype(np.float32)
        # Library from the X row sums (raw counts); must match obs raw_sum.
        row_lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        scale = (1e6 / np.where(row_lib > 0, row_lib, 1.0)).astype(np.float32)
        np.multiply(X.data, np.repeat(scale, np.diff(X.indptr)), out=X.data)
        np.log1p(X.data, out=X.data)
        block_lib = lib[start:start + CELL_CHUNK]
        if np.isfinite(block_lib).all():
            assert np.allclose(row_lib, block_lib, rtol=1e-6, atol=1e-6), \
                f"{cell_type}: X row sums != obs raw_sum"
        ind = sp.csr_matrix(
            (np.ones(len(chunk_rows), dtype=np.float32),
             (np.arange(len(chunk_rows)), dpos[start:start + CELL_CHUNK])),
            shape=(len(chunk_rows), len(donors)),
        )
        # ind.T @ X is already per-donor (rows align with `sums` rows).
        block = (ind.T @ X).toarray().astype(np.float64)
        sums += block
        np.add.at(counts, dpos[start:start + CELL_CHUNK], 1)
        np.add.at(library, dpos[start:start + CELL_CHUNK], row_lib)
        print(f"  [{cell_type}] rows {start + 1}-{min(start + CELL_CHUNK, len(rows))}"
              f" of {len(rows)}", flush=True)
        del X, block, ind
        gc.collect()
    assert (counts >= MIN_CELLS_PER_DONOR_CT).all(), f"short donor in {cell_type}"
    pseudobulk = (sums / counts[:, None]).astype(np.float32)
    meta = obs_index.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"))
    meta = meta.loc[donors]
    save_npz_retry(
        out, genes=np.asarray(genes, dtype=str), donors=np.asarray(donors, dtype=str),
        pseudobulk=pseudobulk, n_cells=counts, library=library,
        age=meta["age"].to_numpy(float),
        sex=meta["sex"].astype(str).to_numpy(dtype=str),
    )
    print(f"  [{cell_type}] saved {out.name}: {len(donors)} donors x {n_genes} genes", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    parser.add_argument("--rebuild-pseudobulk", action="store_true")
    args = parser.parse_args()
    scratch: Path = args.scratch_dir
    scratch.mkdir(parents=True, exist_ok=True)
    if args.rebuild_pseudobulk:
        for f in scratch.glob("mssm_pseudobulk_*.npz"):
            f.unlink()
    if not args.h5ad.exists():
        raise SystemExit(f"h5ad not found: {args.h5ad} (run download_h5ad.sh first)")

    import anndata

    program = load_program_tiers(REPO, "human")["P1_Clean"]
    print(f"P1_Clean: {len(program['up'])} up, {len(program['down'])} down human genes")
    m9_regs = pd.read_csv(RESULTS / "m9-regressions.csv")
    mssm_regs = m9_regs[m9_regs["cohort"] == COHORT].copy()
    cell_types = sorted(mssm_regs["cell_type"])
    print(f"MSSM primary cell types (M9 family): {len(cell_types)}")
    m9_scores = pd.read_csv(RESULTS / "m9-donor-scores.csv")
    m9_scores = m9_scores[m9_scores["cohort"] == COHORT]

    print(f"Opening {args.h5ad} (backed) ...", flush=True)
    adata = anndata.read_h5ad(str(args.h5ad), backed="r")
    print(f"h5ad: {adata.n_obs} cells x {adata.n_vars} genes; "
          f"obs columns: {list(adata.obs.columns)}", flush=True)
    obs = adata.obs
    # Slim immediately: keep only needed columns (one load; obs_keys() lists
    # columns without loading the 4.1M-row frame a second time).
    keep = ["donor_id", "development_stage", "disease", "cell_type", "sex", "is_primary_data"]
    has_raw_sum = "raw_sum" in adata.obs_keys()
    if has_raw_sum:
        keep.append("raw_sum")
    missing = [c for c in keep if c not in obs.columns]
    if missing:
        raise SystemExit(f"h5ad obs missing columns: {missing}")
    obs = obs[keep].copy()
    obs["_row"] = np.arange(len(obs))  # positional row in the full h5ad
    # Keep categorical dtypes: astype(str) on 4.1M rows materialises ~1 GB
    # of Python strings per column and was the OOM spike that killed the
    # sibling MSSM fetch twice (exit 137). Categorical == / groupby give
    # identical selection; age is extracted on the ~62 distinct stage
    # categories and mapped back via codes.
    if isinstance(obs["development_stage"].dtype, pd.CategoricalDtype):
        _cats = (obs["development_stage"].cat.categories.astype(str)
                 .str.extract(r"(\d+)-year-old").iloc[:, 0].astype(float))
        _codes = obs["development_stage"].cat.codes.to_numpy()
        _age = np.full(len(obs), np.nan)
        _ok = _codes >= 0
        _age[_ok] = _cats.to_numpy()[_codes[_ok]]
        obs["age"] = _age
    else:
        obs["age"] = obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
    _primary = obs["is_primary_data"]
    _primary_ok = (_primary.to_numpy() if _primary.dtype == bool else
                   _primary.astype(str).isin(["True", "true", "1"]).to_numpy())
    obs = obs[_primary_ok & (obs["disease"] == "normal").to_numpy()
              & obs["age"].notna().to_numpy()].copy()
    obs["row"] = obs["_row"].to_numpy()
    # Library per cell from X row sums (computed per chunk); obs raw_sum, if
    # present, is used only as a cross-check inside the builder.
    if not has_raw_sum:
        obs["raw_sum"] = np.nan
    donor_meta = (obs.groupby("donor_id", observed=True)
                  .agg(age=("age", "first"), sex=("sex", lambda s: str(s.iloc[0])))
                  .reset_index())
    print(f"Normal control donors with numeric age: {len(donor_meta)}", flush=True)
    counts = (obs[obs["cell_type"] != "unknown"]
              .groupby(["donor_id", "cell_type"], observed=True).size()
              .rename("n_cells").reset_index())
    obs_small = obs[["row", "donor_id", "cell_type", "age", "sex", "raw_sum"]].copy()
    del obs
    gc.collect()

    rng = np.random.default_rng(SEED)
    rows = []
    for cell_type in cell_types:
        ct_counts = counts[counts["cell_type"] == cell_type]
        donors = sorted(ct_counts.loc[ct_counts["n_cells"] >= MIN_CELLS_PER_DONOR_CT, "donor_id"])
        n_pub = int(mssm_regs.loc[mssm_regs["cell_type"] == cell_type, "n_donors"].iloc[0])
        assert len(donors) == n_pub, f"{cell_type}: {len(donors)} donors, M9 reports {n_pub}"
        print(f"[{cell_type}] {len(donors)} donors", flush=True)
        sub = obs_small[(obs_small["cell_type"] == cell_type)
                        & (obs_small["donor_id"].isin(donors))].copy()
        # raw_sum is NaN if the h5ad lacks the column; the builder then
        # skips the cross-check and uses X row sums as the library.
        path = build_pseudobulk_h5ad(adata, sub, cell_type, donors, scratch)
        z = np.load(path, allow_pickle=False)
        genes = np.array([str(g) for g in z["genes"]])
        matrix = z["pseudobulk"].astype(np.float64).T  # genes x donors
        age = z["age"].astype(float)
        sex = np.array([str(s) for s in z["sex"]])
        library = z["library"].astype(float)
        assert list(z["donors"]) == donors
        lookup = {gene: i for i, gene in enumerate(genes)}
        p1_up = np.array([lookup[g] for g in program["up"] if g in lookup], dtype=int)
        p1_down = np.array([lookup[g] for g in program["down"] if g in lookup], dtype=int)

        actual_scores = (matrix[p1_up, :].mean(axis=0) - matrix[p1_down, :].mean(axis=0))
        published = m9_scores[m9_scores["cell_type"] == cell_type].set_index("donor_id").loc[donors]
        max_diff = float(np.max(np.abs(actual_scores - published["score"].to_numpy(float))))
        obs_fit = fit_age_sex(actual_scores, age, sex)
        pub_row = mssm_regs[mssm_regs["cell_type"] == cell_type].iloc[0]
        print(f"  [{cell_type}] P1 validation: max |recomputed - M9 score| = {max_diff:.3g}; "
              f"recomputed slope {obs_fit['slope']:.6g} vs M9 {pub_row['slope_per_year']:.6g}; "
              f"p {obs_fit['p']:.3g} vs M9 {pub_row['p']:.3g}", flush=True)

        available = np.ones(len(genes), dtype=bool)
        available[np.concatenate([p1_up, p1_down])] = False
        random_scores = random_signature_scores(
            {"matrix": matrix, "genes": genes}, p1_up, p1_down, ITERATIONS, rng, available
        )
        random_slopes = np.array([fit_age_sex(s, age, sex)["slope"] for s in random_scores])
        observed = obs_fit["slope"]
        null_mean = float(random_slopes.mean())
        emp_p = (1 + int(np.sum(np.abs(random_slopes - null_mean) >= abs(observed - null_mean)))) / (ITERATIONS + 1)
        depth = np.log(library)
        actual_depth_r = stats.pearsonr(actual_scores, depth).statistic
        random_depth_rs = [stats.pearsonr(s, depth).statistic for s in random_scores]
        save_npz_retry(scratch / f"mssm_null_slopes_{cell_type}.npz".replace(" ", "_").replace("/", "_"),
                       random_slopes=random_slopes)
        rows.append({
            "cohort": COHORT, "cell_type": cell_type, "n_donors": len(donors),
            "n_genes_universe": len(genes), "n_cells_total": int(z["n_cells"].sum()),
            "n_random_signatures": ITERATIONS, "random_seed": SEED,
            "observed_slope": observed, "observed_slope_m9": float(pub_row["slope_per_year"]),
            "observed_p_m9": float(pub_row["p"]),
            "q_global_bh_252": float(pub_row["q_global_bh"]),
            # m9-regressions.csv carries no survives_global_bh_05 column;
            # derive it from the global-BH q (q < 0.05), same definition.
            "survives_global_bh_05": bool(float(pub_row["q_global_bh"]) < 0.05),
            "p1_up_detected": len(p1_up), "p1_down_detected": len(p1_down),
            "p1_coverage": (len(p1_up) + len(p1_down)) / (len(program["up"]) + len(program["down"])),
            "max_abs_score_diff_vs_m9": max_diff,
            "random_mean_slope": null_mean,
            "random_sd": float(random_slopes.std(ddof=1)),
            "random_ci95_low": float(np.quantile(random_slopes, 0.025)),
            "random_ci95_high": float(np.quantile(random_slopes, 0.975)),
            "random_q025": float(np.quantile(random_slopes, 0.025)),
            "random_q25": float(np.quantile(random_slopes, 0.25)),
            "random_median": float(np.quantile(random_slopes, 0.5)),
            "random_q75": float(np.quantile(random_slopes, 0.75)),
            "random_q975": float(np.quantile(random_slopes, 0.975)),
            "empirical_p_two_sided": float(emp_p),
            "mean_log_library": float(depth.mean()),
            "actual_score_depth_r": float(actual_depth_r),
            "random_score_depth_r_mean": float(np.mean(random_depth_rs)),
        })
        print(f"  [{cell_type}] null mean {null_mean:.6g} sd {random_slopes.std(ddof=1):.6g} "
              f"empirical p {emp_p:.4f}", flush=True)

    result = pd.DataFrame(rows)
    result["q_bh_null_family"] = multipletests(result["empirical_p_two_sided"], method="fdr_bh")[1]
    result["survives_null_bh_05"] = result["q_bh_null_family"] < 0.05
    result["survivor"] = result["survives_global_bh_05"] & result["survives_null_bh_05"]
    cols = ["cohort", "cell_type", "n_donors", "n_genes_universe", "n_cells_total",
            "n_random_signatures", "random_seed", "observed_slope", "observed_slope_m9",
            "observed_p_m9", "q_global_bh_252", "survives_global_bh_05",
            "p1_up_detected", "p1_down_detected", "p1_coverage",
            "max_abs_score_diff_vs_m9", "random_mean_slope", "random_sd",
            "random_ci95_low", "random_ci95_high", "random_q025", "random_q25",
            "random_median", "random_q75", "random_q975",
            "empirical_p_two_sided", "q_bh_null_family", "survives_null_bh_05",
            "survivor", "mean_log_library", "actual_score_depth_r",
            "random_score_depth_r_mean"]
    result = result[cols].sort_values("empirical_p_two_sided").reset_index(drop=True)
    out = RESULTS / "reviewer-rev2-transport-null-mssm.csv"
    import time as _time
    for _attempt in range(1, 9):
        try:
            result.to_csv(out, index=False)
            break
        except OSError as _exc:
            print(f"CSV save attempt {_attempt} failed ({_exc}); retrying", flush=True)
            _time.sleep(10 * _attempt)
    else:
        result.to_csv(out, index=False)
    print(f"\nWrote {out} ({len(result)} MSSM cell types)")
    print(result[result["survivor"]].to_string(index=False) if result["survivor"].any()
          else "No cell type survives both the global BH and the matched random-signature null.")
    print(result[["cell_type", "observed_slope", "random_mean_slope", "random_sd",
                  "empirical_p_two_sided", "q_bh_null_family", "survivor"]].to_string(index=False))


if __name__ == "__main__":
    main()
