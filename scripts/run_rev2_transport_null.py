#!/usr/bin/env python3
"""Rev2 fix 10: matched random-signature null in MSSM (M9 cohort).

Replicates, for the CELLxGENE Census MSSM_DLPFC cohort, the matched
random-signature null that Phase 1 ran in PsychAD
(`scripts/run_reviewer_phase1_signature_null.py`, seed 20261105):

* P1_Clean signature: 76 human up genes + 42 down genes (118 total;
  `src.aging_programs.load_program_tiers`, human mapping).
* Per cell type, donor pseudobulk = mean over cells of log1p(CPM) per gene.
  By linearity of means this equals the M9 donor score
  (mean over cells of per-cell up-mean minus down-mean); the recomputed
  P1 slopes are validated against `results/m9-regressions.csv` before the
  null runs.
* Expression matching: deciles of mean logCPM across donors
  (`pd.qcut(means, 10, duplicates="drop")`). For each direction (up, down
  separately) the random set draws, per decile, as many genes as P1 has in
  that decile, sampled without replacement from genes in that decile
  excluding ALL detected P1 genes. 200 draws, `np.random.default_rng(20261105)`,
  one RNG stream across cell types in sorted order, exactly as in Phase 1.
* Random score = mean random-up minus mean random-down per donor; each
  random score is regressed with the cohort's own primary model, the M9
  age-plus-sex OLS (`score ~ 1 + age + male`, male only if both sexes
  present), the same OLS implementation as `scripts/run_m9_pipeline.py`.
* Empirical two-sided p, centred on the null mean:
  (1 + #{|rand - mean| >= |obs - mean|}) / (200 + 1).
* BH q over the MSSM null family (the 21 MSSM cell types with M9 primary
  tests), mirroring Phase 1, which applied BH over its own 27 PsychAD cell
  types. A cell type "survives" only if it survived the original 252-test
  global BH (q < 0.05) AND survives this null BH (q < 0.05).

Why Census MSSM and not the M10 MSSM arm: the two MSSM survivors are M9
`MSSM_DLPFC` rows (214 donors, CELLxGENE Census 2025-11-08, dataset
37a17b78-4864-4a42-b67b-31c00962795a). The M10 atlas MSSM arm is a different
donor set (112 donors; 103 exact identifiers overlap the 214, see
`results/m9-m10-donor-overlap-audit.csv`) with different cell-type labels,
so running the null there would not test the reported survivors. This script
therefore re-reads the Census MSSM dataset with the same cell selection as
M9 (primary data, donors with at least one normal cell, >= 20 cells per
donor-cell-type) and all genes in the dataset var.

Network note: TileDB's bundled libcurl does not honour the sandbox
HTTPS_PROXY env vars, so a direct `open_soma()` fails at the S3 listing
with an SSL "wrong version number" error. `tiledb_config_from_env()`
therefore builds a TileDB S3 proxy config from the HTTPS_PROXY env var at
runtime (host/port/scheme/user/password parsed, never hardcoded) plus the
egress CA bundle from CURL_CA_BUNDLE, and passes it to `open_soma()`.
Verified 2026-11-09: SOMA open + MSSM obs read succeed with this config.

Ambiguities in the Phase 1 construction and the choices made here:
1. Phase 1 did not state whether the MSSM analogue should use the M9
   age-plus-sex model or the M10 fully adjusted model. Choice: the M9
   model, because the MSSM survivors are M9 tests and the null must test
   the statistic that produced them. (The M10 covariates PMI/TechPCs do
   not exist in Census obs.)
2. Phase 1's BH family was its own 27 PsychAD cell types. Choice: BH over
   the 21 MSSM cell types with M9 primary tests (the MSSM analogue).
   Raw empirical p-values are reported so any other family can be derived.
3. Phase 1 computed deciles from the mean of the log-normalised pseudobulk
   matrix across donors. Choice: same, with donors weighted equally
   (mean over donors of donor pseudobulk), not weighted by cell count.
4. Phase 1's pool excluded detected P1 genes only. Choice: same; the pool
   is every gene in the Census MSSM var (no detection filter), minus the
   detected P1 genes. All 118 P1 genes are present in the Census var
   (M9 coverage 1.0).

Usage (CPU only; downloads the MSSM slice of CELLxGENE Census 2025-11-08):
    python scripts/run_rev2_transport_null.py [--scratch-dir DIR]
        [--rebuild-pseudobulk] [--skip-download]
Writes results/reviewer-rev2-transport-null.csv. Pseudobulk intermediates
are cached as .npz per cell type under the scratch dir (default
~/workspace/rev2-mssm-null-scratch) so a rerun resumes without re-download.
No .tex or manuscript file is touched; nothing is committed or pushed.
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

import cellxgene_census

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers  # noqa: E402

RESULTS = REPO / "results"
CENSUS_VERSION = "2025-11-08"
MSSM_DATASET_ID = "37a17b78-4864-4a42-b67b-31c00962795a"
COHORT = "MSSM_DLPFC"
ITERATIONS = 200
SEED = 20261105
MIN_CELLS_PER_DONOR_CT = 20
DONOR_BATCH = 5  # small on purpose: the VM is shared with other Rev2 workers
DEFAULT_SCRATCH = Path.home() / "workspace" / "rev2-mssm-null-scratch"


def ols(y, X):
    """OLS slope on X[:, 1]; identical to run_m9_pipeline.ols."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    s2 = (resid @ resid) / max(dof, 1)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(cov), 1e-300))
    tstat = beta[1] / se[1]
    p = 2 * stats.t.sf(abs(tstat), max(dof, 1))
    tcrit = stats.t.ppf(0.975, max(dof, 1))
    return {
        "slope": beta[1], "se": se[1], "p": p,
        "ci_lo": beta[1] - tcrit * se[1], "ci_hi": beta[1] + tcrit * se[1],
        "rho": stats.spearmanr(X[:, 1], y).statistic,
        "rho_p": stats.spearmanr(X[:, 1], y).pvalue,
    }


def fit_age_sex(score, age, sex):
    """M9 cell-type model: score ~ 1 + age (+ male if both sexes present)."""
    male = (pd.Series(sex).astype(str).str.lower() == "male").astype(float).values
    columns = [np.ones(len(age)), np.asarray(age, float)]
    if len(np.unique(male)) > 1:
        columns.append(male)
    return ols(np.asarray(score, float), np.column_stack(columns))


def random_signature_scores(expression, p1_up, p1_down, iterations, rng, available):
    """Verbatim logic of Phase 1 (run_reviewer_phase1_signature_null.py).

    expression: dict with genes x donors logCPM pseudobulk matrix.
    Deciles of mean expression across donors; per direction and decile,
    draw as many genes as P1 has there, without replacement, from that
    decile's pool excluding all P1 genes.
    """
    genes = np.asarray(expression["genes"])
    matrix = expression["matrix"]
    means = np.asarray(matrix.mean(axis=1)).ravel()
    deciles = pd.qcut(means, 10, labels=False, duplicates="drop")
    deciles = pd.Series(deciles).fillna(-1).to_numpy()
    rows = []
    for _ in range(iterations):
        selected_by_direction = []
        for program_genes in (p1_up, p1_down):
            selected = []
            for decile, count in pd.Series(deciles[program_genes]).value_counts().items():
                pool = np.flatnonzero((deciles == decile) & available)
                selected.extend(rng.choice(pool, size=int(count), replace=False))
            selected_by_direction.append(selected)
        up_score = np.asarray(matrix[selected_by_direction[0], :].mean(axis=0)).ravel()
        down_score = np.asarray(matrix[selected_by_direction[1], :].mean(axis=0)).ravel()
        rows.append(up_score - down_score)
    return np.vstack(rows)


def tiledb_config_from_env():
    """TileDB S3 config routing through the sandbox egress proxy.

    TileDB's bundled libcurl ignores HTTPS_PROXY, so without this the S3
    listing fails with an SSL 'wrong version number' error. Credentials
    come from the environment at runtime and are never written to disk.
    """
    import os
    from urllib.parse import urlparse

    cfg = {"vfs.s3.region": "us-west-2",
           # The egress proxy throttles S3 heavily (HEAD latency ~20 s);
           # keep TileDB's S3 parallelism low and give requests room.
           "vfs.s3.max_parallel_ops": "2",
           "vfs.s3.connect_timeout_ms": "60000",
           "vfs.s3.request_timeout_ms": "600000",
           "vfs.s3.connect_max_tries": "5"}
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
    return cfg


def save_npz_retry(path: Path, tries: int = 8, **arrays) -> None:
    """np.savez_compressed with retries: sandbox file writes flap EPERM."""
    import time
    last = None
    for attempt in range(1, tries + 1):
        try:
            np.savez_compressed(path, **arrays)
            return
        except OSError as exc:  # noqa: BLE001 - sandbox write flap, retry
            last = exc
            print(f"  save {path.name} attempt {attempt} failed ({exc}); retrying", flush=True)
            time.sleep(10 * attempt)
    raise last


def build_pseudobulk(census, donor_meta, cell_type, donors, scratch: Path) -> Path:
    """Donor pseudobulk (mean logCPM per gene) for one MSSM cell type.

    donors: donor_ids with >= 20 cells in this cell type, sorted.
    Cell selection mirrors M9 exactly: primary data, donors from the
    normal-cell donor list, this cell type; all genes in the var.
    """
    out = scratch / f"mssm_pseudobulk_{cell_type}.npz".replace(" ", "_").replace("/", "_")
    if out.exists():
        return out
    n_genes = None
    gene_names = None
    sums = None
    counts = np.zeros(len(donors), dtype=np.int64)
    library = np.zeros(len(donors), dtype=np.float64)
    donor_pos = {d: i for i, d in enumerate(donors)}
    for start in range(0, len(donors), DONOR_BATCH):
        batch = donors[start:start + DONOR_BATCH]
        donor_list = ",".join(f"'{d}'" for d in batch)
        filt = (f"dataset_id == '{MSSM_DATASET_ID}' and is_primary_data == True "
                f"and donor_id in [{donor_list}] and cell_type == '{cell_type}'")
        adata = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=filt,
            obs_column_names=["donor_id", "cell_type", "raw_sum"],
            var_column_names=["feature_name"],
        )
        names = list(adata.var["feature_name"])
        if gene_names is None:
            gene_names = names
            n_genes = len(names)
            sums = np.zeros((len(donors), n_genes), dtype=np.float64)
        else:
            assert names == gene_names, f"var changed between batches for {cell_type}"
        raw_sum = np.asarray(adata.obs["raw_sum"], float)
        batch_donor_idx = np.array([donor_pos[str(d)] for d in adata.obs["donor_id"]])
        # float32 + in-place per-cell CPM scaling: the full-gene slice is
        # the memory spike, so avoid float64 copies and diag-matmul copies.
        X = sp.csr_matrix(adata.X).astype(np.float32)
        scale = (1e6 / np.where(raw_sum > 0, raw_sum, 1.0)).astype(np.float32)
        np.multiply(X.data, np.repeat(scale, np.diff(X.indptr)), out=X.data)
        np.log1p(X.data, out=X.data)
        ind = sp.csr_matrix(
            (np.ones(len(batch_donor_idx), dtype=np.float32),
             (np.arange(len(batch_donor_idx)), batch_donor_idx)),
            shape=(len(batch_donor_idx), len(batch)),
        )
        block = (ind.T @ X).toarray().astype(np.float64)  # batch donors x genes
        for j, d in enumerate(batch):
            sums[donor_pos[d]] += block[j]
        np.add.at(counts, batch_donor_idx, 1)
        np.add.at(library, batch_donor_idx, raw_sum)
        print(f"  [{cell_type}] donors {start + 1}-{min(start + DONOR_BATCH, len(donors))}"
              f" of {len(donors)}: {adata.n_obs} cells", flush=True)
        del adata, X, block, ind
        gc.collect()
    assert (counts >= MIN_CELLS_PER_DONOR_CT).all(), f"short donor in {cell_type}"
    pseudobulk = (sums / counts[:, None]).astype(np.float32)
    meta = donor_meta.set_index("donor_id").loc[donors]
    save_npz_retry(
        out, genes=np.array(gene_names), donors=np.array(donors),
        pseudobulk=pseudobulk, n_cells=counts, library=library,
        age=meta["age"].to_numpy(float),
        sex=meta["sex"].astype(str).to_numpy(),
    )
    print(f"  [{cell_type}] saved {out.name}: {len(donors)} donors x {n_genes} genes", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    parser.add_argument("--rebuild-pseudobulk", action="store_true",
                        help="delete cached pseudobulk .npz files and re-download")
    args = parser.parse_args()
    scratch: Path = args.scratch_dir
    scratch.mkdir(parents=True, exist_ok=True)
    if args.rebuild_pseudobulk:
        for f in scratch.glob("mssm_pseudobulk_*.npz"):
            f.unlink()

    program = load_program_tiers(REPO, "human")["P1_Clean"]
    print(f"P1_Clean: {len(program['up'])} up, {len(program['down'])} down human genes")

    m9_regs = pd.read_csv(RESULTS / "m9-regressions.csv")
    mssm_regs = m9_regs[m9_regs["cohort"] == COHORT].copy()
    cell_types = sorted(mssm_regs["cell_type"])
    print(f"MSSM primary cell types (M9 family): {len(cell_types)}")
    m9_scores = pd.read_csv(RESULTS / "m9-donor-scores.csv")
    m9_scores = m9_scores[m9_scores["cohort"] == COHORT]

    print(f"Opening CELLxGENE Census {CENSUS_VERSION}")
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION,
                                         tiledb_config=tiledb_config_from_env())
    # Filter server-side to primary normal cells: the unfiltered MSSM obs
    # is 4.1M rows and its Arrow-to-pandas conversion spikes memory on a
    # VM shared with other Rev2 workers (an earlier run was OOM-killed
    # while loading it). The filtered frame is 842,830 rows.
    obs = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=(f"dataset_id == '{MSSM_DATASET_ID}' and is_primary_data == True "
                      f"and disease == 'normal'"),
        column_names=["donor_id", "development_stage", "cell_type", "sex"],
    )
    obs["donor_id"] = obs["donor_id"].astype(str)
    import re
    obs["age"] = obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
    obs = obs[obs["age"].notna()].copy()
    donor_meta = (obs.groupby(["donor_id"], observed=True)
                  .agg(age=("age", "first"), sex=("sex", lambda s: str(s.iloc[0])))
                  .reset_index())
    print(f"Normal control donors with numeric age: {len(donor_meta)}")

    counts = (obs[obs["cell_type"].astype(str) != "unknown"]
              .groupby(["donor_id", "cell_type"], observed=True).size()
              .rename("n_cells").reset_index())
    # The unfiltered obs is 4.1M rows; drop it now to stay lean on a VM
    # shared with other Rev2 workers (an earlier run was OOM-killed).
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
        path = build_pseudobulk(census, donor_meta, cell_type, donors, scratch)
        z = np.load(path, allow_pickle=False)
        genes = np.array([str(g) for g in z["genes"]])
        matrix = z["pseudobulk"].astype(np.float64)  # donors x genes
        matrix = matrix.T  # genes x donors, as in Phase 1
        age = z["age"].astype(float)
        sex = np.array([str(s) for s in z["sex"]])
        library = z["library"].astype(float)
        assert list(z["donors"]) == donors
        lookup = {gene: i for i, gene in enumerate(genes)}
        p1_up = np.array([lookup[g] for g in program["up"] if g in lookup], dtype=int)
        p1_down = np.array([lookup[g] for g in program["down"] if g in lookup], dtype=int)

        # Observed P1 donor scores from pseudobulk; validate against M9.
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
            "survives_global_bh_05": bool(pub_row["survives_global_bh_05"]),
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
    out = RESULTS / "reviewer-rev2-transport-null.csv"
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
