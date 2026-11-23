#!/usr/bin/env python3
"""Held-out validation of the cleaned mature-to-old signature (P1m).

Question: does P1m (1,914 genes, 18m vs 30m, kidney/lung/marrow) capture real
age signal, or is it noise fitted to 8 mice? The cross-species transport test
mostly fails, and without a positive control that failure is uninterpretable.

Four variants, all with the frozen derivation rules (pseudobulk counts per
mouse per cell type, design sex plus age group, PyDESeq2, contrast 30m vs 18m,
eligibility at least 3 mice and at least 50 cells on both sides, FDR below
0.05 and absolute log2 fold change at least 0.25, same clean-gene filtering
and floor-10 rule) and the frozen transport rules (score = mean log1p-CPM of
up genes minus mean log1p-CPM of down genes per donor per cell type;
composite test = OLS of score on age group plus sex; empirical test =
2,000-draw matched-random null, random sets matched on up/down counts and
expression-decile distribution with probe genes excluded; BH within family):

1. CROSS-ORGAN 3-FOLD: derive a fresh signature from two organs, test it in
   the third organ it never saw.
2. INTERMEDIATE AGES: score existing P1m in the 21m (2 mice) and 24m
   (4 mice) animals, which never entered the derivation.
3. DONOR SPLIT (descriptive only): derive on 3 mice per side, check the
   direction of the held-out mouse's score. No p-value is possible.
4. FACS ARM: NOT RUN (no FACS files locally; and TMS FACS has no 30m mice,
   so the 18m-vs-30m test is impossible in FACS even after download).

Scoring and null code replicate the frozen implementations verbatim in
method: run_m9_pipeline.ols, the matched_null_effect_pvalue centering rule
from src/followup_muscle.py, and the decile-matching rule of
random_signature_scores from scripts/run_rev2_transport_null.py. The random
draws are generated with a vectorized permutation (argsort of uniforms),
which is the same sampling distribution as choice-without-replacement, with
a different random stream. Seed is fixed and recorded.

Model used for this implementation and report: Muse Spark (subagent).
"""

from __future__ import annotations

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import json
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import scripts.run_mature_p1m_derivation as deriv
from scripts.run_m9_pipeline import ols as m9_ols

RESULTS = REPO / "results"
SCRATCH = REPO / "scratch" / "heldout"
SCRATCH.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(exist_ok=True)

SEED = 20261122
ITERATIONS = 2000
MIN_CELLS_PER_DONOR_CT = 20
TISSUES = ["Kidney", "Lung", "Marrow"]
ALL_AGES = ["18m", "21m", "24m", "30m"]
AGE_MONTHS = {"18m": 18, "21m": 21, "24m": 24, "30m": 30}
MODEL = "Muse Spark (subagent)"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def matched_null_effect_pvalue(observed: float, null_effects: np.ndarray) -> float:
    """Verbatim rule from src/followup_muscle.py (two-sided, centered)."""
    centered = null_effects - null_effects.mean()
    return float(
        (1 + np.sum(np.abs(centered) >= abs(observed - null_effects.mean())))
        / (len(null_effects) + 1)
    )


def build_expression_matrix(tissue: str):
    """Mean log1p-CPM per (mouse, cell type) for one tissue, all four ages.

    Per-cell counts are scaled to counts per million, log1p-transformed, then
    averaged over the cells of each mouse-by-cell-type pair. This is the same
    object the Block 5 matched null operates on (genes x donors).
    Returns (genes, pairs DataFrame, matrix genes x pairs).
    """
    cache = SCRATCH / f"expr-{tissue}.npz"
    if cache.exists():
        with np.load(cache, allow_pickle=False) as z:
            genes = z["genes"].astype(str)
            pairs = pd.DataFrame(
                {k: z[k] for k in ("mouse_id", "age", "sex", "cell_type", "n_cells")}
            )
            pairs["n_cells"] = pairs["n_cells"].astype(int)
            return genes, pairs, z["matrix"]
    path = deriv.DATA_DIR / f"m1-qc-cleaned-{tissue}.h5ad"
    adata = ad.read_h5ad(path)
    obs = adata.obs
    keep = obs["age"].isin(ALL_AGES) & obs["cell_ontology_class"].notna()
    adata = adata[keep].copy()
    X = adata.X.tocsr().astype(np.float64)
    lib = np.asarray(X.sum(axis=1)).ravel()
    scale = 1e6 / np.where(lib > 0, lib, 1.0)
    X.data *= np.repeat(scale, np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    genes = np.asarray(adata.var_names).astype(str)
    pair_keys = pd.DataFrame(
        {
            "mouse_id": adata.obs["mouse.id"].astype(str).to_numpy(),
            "age": adata.obs["age"].astype(str).to_numpy(),
            "sex": adata.obs["sex"].astype(str).to_numpy(),
            "cell_type": adata.obs["cell_ontology_class"].astype(str).to_numpy(),
        }
    )
    codes, uniques = pd.factorize(
        list(zip(pair_keys["mouse_id"], pair_keys["age"], pair_keys["sex"], pair_keys["cell_type"])),
        sort=True,
    )
    n_pairs = len(uniques)
    indicator = sp.csr_matrix(
        (np.ones(len(codes)), (codes, np.arange(len(codes)))),
        shape=(n_pairs, len(codes)),
    )
    pair_sums = (indicator @ X).toarray().T  # genes x pairs
    pair_counts = np.bincount(codes, minlength=n_pairs).astype(float)
    matrix = pair_sums / pair_counts[None, :]
    pairs = pd.DataFrame(
        [dict(zip(("mouse_id", "age", "sex", "cell_type"), u)) for u in uniques]
    )
    pairs["n_cells"] = pair_counts.astype(int)
    np.savez_compressed(
        cache,
        genes=np.asarray(genes, dtype=str),
        matrix=matrix,
        mouse_id=np.asarray(pairs["mouse_id"], dtype=str),
        age=np.asarray(pairs["age"], dtype=str),
        sex=np.asarray(pairs["sex"], dtype=str),
        cell_type=np.asarray(pairs["cell_type"], dtype=str),
        n_cells=pairs["n_cells"].to_numpy(),
    )
    print(f"  {tissue}: {len(genes)} genes, {n_pairs} mouse-cell-type pairs", flush=True)
    return genes, pairs, matrix


def score_from_matrix(matrix, genes, up_genes, down_genes):
    """Up-minus-down mean logCPM score per column (donor)."""
    lookup = {g: i for i, g in enumerate(genes)}
    up_idx = np.array([lookup[g] for g in up_genes if g in lookup], dtype=int)
    down_idx = np.array([lookup[g] for g in down_genes if g in lookup], dtype=int)
    if len(up_idx) == 0 or len(down_idx) == 0:
        return None, len(up_idx), len(down_idx)
    scores = matrix[up_idx, :].mean(axis=0) - matrix[down_idx, :].mean(axis=0)
    return np.asarray(scores).ravel(), len(up_idx), len(down_idx)


def fit_unit(scores, is_old, male):
    """OLS score ~ 1 + is_old (+ male when it varies). Returns dict."""
    cols = [np.ones(len(scores)), is_old.astype(float)]
    sex_included = len(np.unique(male)) > 1
    if sex_included:
        cols.append(male.astype(float))
    X = np.column_stack(cols)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        X = X[:, :2]
        sex_included = False
    fit = m9_ols(scores, X)
    return {
        "slope": float(fit["slope"]),
        "se": float(fit["se"]),
        "pvalue": float(fit["p"]),
        "sex_included": bool(sex_included),
        "n_donors": len(scores),
    }


def random_scores_vectorized(matrix, genes, up_idx, down_idx, iterations, rng):
    """2000 matched-random signature scores, same distribution as
    random_signature_scores in run_rev2_transport_null.py: deciles of mean
    expression across donors; per direction and decile, draw as many genes as
    the program has there, without replacement, from that decile's pool
    excluding all probe genes. Draws use argsort of uniforms (a random
    permutation), equivalent to choice without replacement.
    Returns array (iterations, n_donors)."""
    means = np.asarray(matrix.mean(axis=1)).ravel()
    deciles = pd.qcut(means, 10, labels=False, duplicates="drop")
    deciles = pd.Series(deciles).fillna(-1).to_numpy()
    available = np.ones(len(genes), dtype=bool)
    available[np.concatenate([up_idx, down_idx])] = False
    n_donors = matrix.shape[1]
    out_up = np.zeros((iterations, n_donors))
    out_down = np.zeros((iterations, n_donors))
    for target, up in ((out_up, True), (out_down, False)):
        prog = up_idx if up else down_idx
        total = len(prog)
        counts = pd.Series(deciles[prog]).value_counts()
        for decile, count in counts.items():
            pool = np.flatnonzero((deciles == decile) & available)
            k = int(count)
            if len(pool) < k:
                raise ValueError(f"decile {decile} pool {len(pool)} < needed {k}")
            u = rng.random((iterations, len(pool)))
            order = np.argsort(u, axis=1)[:, :k]
            sel = pool[order]  # (iterations, k)
            # sum over selected genes per draw, in chunks to bound memory
            acc = np.zeros((iterations, n_donors))
            for start in range(0, iterations, 250):
                stop = min(start + 250, iterations)
                acc[start:stop] = matrix[sel[start:stop]].sum(axis=1)
            target += acc / total
    return out_up - out_down


def empirical_pvalue(matrix, genes, up_genes, down_genes, is_old, male, observed_slope, rng):
    """Empirical p for one unit from the matched-random null."""
    lookup = {g: i for i, g in enumerate(genes)}
    up_idx = np.array([lookup[g] for g in up_genes if g in lookup], dtype=int)
    down_idx = np.array([lookup[g] for g in down_genes if g in lookup], dtype=int)
    rand_scores = random_scores_vectorized(matrix, genes, up_idx, down_idx, ITERATIONS, rng)
    rand_slopes = np.array(
        [fit_unit(rand_scores[i], is_old, male)["slope"] for i in range(ITERATIONS)]
    )
    p = matched_null_effect_pvalue(observed_slope, rand_slopes)
    return p, float(rand_slopes.mean()), float(rand_slopes.std(ddof=1))


def derive_signature(derive_tissues, allow_mice=None, tag=""):
    """Derive a P1m-style signature from the given tissues (and, if given,
    only the given mice), using the frozen derivation functions and rules.
    Returns (gene DataFrame with P1m columns, ranking DataFrame)."""
    chrom_df = pd.read_csv(deriv.REF_DIR / "mouse-gene-chromosomes.csv").set_index(
        "mouse_gene_symbol"
    )
    biotypes = chrom_df["gene_biotype"].to_dict()
    chromosomes = chrom_df["chromosome_name"].to_dict()
    orthologs = pd.read_csv(deriv.REF_DIR / "orthologs-mouse-human.csv")
    ortholog_map = dict(zip(orthologs["mouse_gene_symbol"], orthologs["human_gene_symbol"]))

    eligibility, _, _ = get_eligibility()
    if allow_mice is not None:
        # Recompute eligibility on the allowed mice from the expression pairs.
        rows = []
        for tissue in derive_tissues:
            _, pairs, _ = EXPR[tissue]
            sub = pairs[pairs["mouse_id"].isin(allow_mice)]
            for ctype, df in sub.groupby("cell_type"):
                row = {"tissue": tissue, "cell_ontology_class": ctype,
                       "unit": f"{tissue}::{ctype}"}
                ok = True
                for age in deriv.PRIMARY_AGES:
                    adf = df[df["age"] == age]
                    row[f"n_mice_{age}"] = adf["mouse_id"].nunique()
                    row[f"n_cells_{age}"] = int(adf["n_cells"].sum())
                    row[f"mice_ge10_{age}"] = int((adf["n_cells"] >= deriv.FLOOR10_CELLS).sum())
                    if row[f"n_mice_{age}"] < deriv.MIN_MICE or row[f"n_cells_{age}"] < deriv.MIN_CELLS:
                        ok = False
                row["eligible"] = "yes" if ok else "no"
                row["p1m_floor10_eligible"] = bool(
                    ok and all(row[f"mice_ge10_{age}"] >= deriv.MIN_MICE for age in deriv.PRIMARY_AGES)
                )
                rows.append(row)
        eligibility = pd.DataFrame(rows)
    eligible = eligibility[
        (eligibility["tissue"].isin(derive_tissues)) & (eligibility["eligible"] == "yes")
    ].copy()
    print(f"  {tag}: {len(eligible)} eligible units in {derive_tissues}", flush=True)

    p1_entries = []
    ranking_rows = []
    for tissue, tissue_units in eligible.groupby("tissue", sort=False):
        source_file = tissue_units["source_file"].iloc[0] if "source_file" in tissue_units else f"m1-qc-cleaned-{tissue}.h5ad"
        adata = ad.read_h5ad(deriv.DATA_DIR / source_file, backed="r")
        for row in tissue_units.itertuples(index=False):
            unit = row.unit
            if allow_mice is not None:
                counts, metadata, n_cells = pseudobulk_unit_filtered(
                    adata, row.tissue, row.cell_ontology_class, allow_mice
                )
            else:
                counts, metadata, n_cells = deriv.pseudobulk_unit(
                    adata, row.tissue, row.cell_ontology_class
                )
            result, diagnostics, _ = deriv.run_de(counts, metadata)
            ranking_rows.append(
                {
                    "unit": unit,
                    "tissue": row.tissue,
                    "cell_ontology_class": row.cell_ontology_class,
                    "n_mice_total": len(metadata),
                    "n_cells_total": n_cells,
                    "n_de_uniform_before_clean": int(result["P1m_significant_before_filters"].sum()),
                    "p1m_floor10_eligible": bool(row.p1m_floor10_eligible),
                    **diagnostics,
                }
            )
            if row.p1m_floor10_eligible:
                sig = result[result["P1m_significant_before_filters"]]
                for gene, gene_row in sig.iterrows():
                    if deriv.is_p1_clean(gene, unit, biotypes):
                        p1_entries.append(
                            {
                                "unit": unit,
                                "mouse_gene_symbol": gene,
                                "log2FC": gene_row["log2FC"],
                                "FDR_uniform": gene_row["FDR_uniform"],
                                "pvalue": gene_row["pvalue"],
                            }
                        )
            print(f"    {unit}: uniform-BH DE before clean="
                  f"{int(result['P1m_significant_before_filters'].sum())}", flush=True)
        adata.file.close()
    genes = deriv.aggregate_genes(p1_entries, "FDR_uniform", ortholog_map, biotypes, chromosomes)
    ranking = pd.DataFrame(ranking_rows)
    print(f"  {tag}: signature {len(genes)} genes "
          f"({int((genes['direction'] == 'UP').sum())} up, "
          f"{int((genes['direction'] == 'DOWN').sum())} down)", flush=True)
    return genes, ranking


def pseudobulk_unit_filtered(adata, tissue, ctype, allow_mice):
    """Copy of deriv.pseudobulk_unit with a mouse allowlist applied first."""
    mask = (
        (adata.obs["tissue"] == tissue)
        & (adata.obs["cell_ontology_class"] == ctype)
        & (adata.obs["age"].isin(deriv.PRIMARY_AGES))
        & (adata.obs["mouse.id"].isin(allow_mice))
    )
    sub = adata[mask].to_memory()
    mouse_ids = []
    metadata_rows = []
    bulk_rows = []
    for mouse_id, mouse_df in sub.obs.groupby("mouse.id", observed=True):
        cell_idx = np.where((sub.obs["mouse.id"] == mouse_id).to_numpy())[0]
        counts = np.asarray(sub.X[cell_idx].sum(axis=0)).reshape(-1)
        bulk_rows.append(counts)
        mouse_ids.append(str(mouse_id))
        metadata_rows.append(
            {
                "mouse.id": str(mouse_id),
                "age_group": str(mouse_df["age"].iloc[0]),
                "sex": str(mouse_df["sex"].iloc[0]),
            }
        )
    counts_df = pd.DataFrame(bulk_rows, index=mouse_ids, columns=sub.var_names)
    metadata = pd.DataFrame(metadata_rows, index=mouse_ids)
    metadata["age_group"] = pd.Categorical(metadata["age_group"], categories=deriv.PRIMARY_AGES)
    metadata["sex"] = pd.Categorical(metadata["sex"], categories=["female", "male"])
    return counts_df, metadata, len(sub)


def test_signature_in_tissue(genes_df, tissue, rng, tag):
    """Run both transport tests for a signature in every eligible unit of one
    tissue. Returns (transport DataFrame, scores DataFrame)."""
    up_genes = genes_df.loc[genes_df["direction"] == "UP", "mouse_gene_symbol"].tolist()
    down_genes = genes_df.loc[genes_df["direction"] == "DOWN", "mouse_gene_symbol"].tolist()
    genes, pairs, matrix = EXPR[tissue]
    eligibility, _, _ = get_eligibility()
    units = eligibility[
        (eligibility["tissue"] == tissue) & (eligibility["eligible"] == "yes")
    ]["cell_ontology_class"].tolist()
    rows = []
    score_rows = []
    for ctype in units:
        unit = f"{tissue}::{ctype}"
        sub_pairs = pairs[
            (pairs["cell_type"] == ctype)
            & (pairs["age"].isin(deriv.PRIMARY_AGES))
            & (pairs["n_cells"] >= MIN_CELLS_PER_DONOR_CT)
        ].copy()
        n_18 = int(((sub_pairs["age"] == "18m")).sum())
        n_30 = int(((sub_pairs["age"] == "30m")).sum())
        base = {"unit": unit, "tissue": tissue, "cell_type": ctype,
                "n_mice_18m": n_18, "n_mice_30m": n_30}
        if n_18 < deriv.MIN_MICE or n_30 < deriv.MIN_MICE:
            rows.append({**base, "status": "UNTESTABLE",
                         "reason": "fewer than 3 mice per side after 20-cell floor"})
            continue
        idx = sub_pairs.index.to_numpy()
        sub_matrix = matrix[:, idx]
        scores, n_up, n_down = score_from_matrix(sub_matrix, genes, up_genes, down_genes)
        if scores is None:
            rows.append({**base, "status": "UNTESTABLE",
                         "reason": "no detected up or down genes in tissue"})
            continue
        is_old = (sub_pairs["age"].to_numpy() == "30m")
        male = (sub_pairs["sex"].to_numpy() == "male").astype(float)
        fit = fit_unit(scores, is_old, male)
        emp_p, rand_mean, rand_sd = empirical_pvalue(
            sub_matrix, genes, up_genes, down_genes, is_old, male, fit["slope"], rng
        )
        rows.append({**base, "status": "TESTED",
                     "n_up_detected": n_up, "n_down_detected": n_down,
                     "composite_slope": fit["slope"], "composite_se": fit["se"],
                     "composite_pvalue": fit["pvalue"], "sex_included": fit["sex_included"],
                     "empirical_pvalue": emp_p, "random_mean_slope": rand_mean,
                     "random_sd_slope": rand_sd, "null_iterations": ITERATIONS})
        for (mouse_id, age, sex, n_cells), score in zip(
            sub_pairs[["mouse_id", "age", "sex", "n_cells"]].itertuples(index=False), scores
        ):
            score_rows.append({"unit": unit, "tissue": tissue, "cell_type": ctype,
                               "mouse_id": mouse_id, "age": age, "sex": sex,
                               "n_cells": int(n_cells), "score": float(score)})
        print(f"    {unit}: slope={fit['slope']:.4f} p={fit['pvalue']:.3g} "
              f"emp_p={emp_p:.3g}", flush=True)
    transport = pd.DataFrame(rows)
    tested = transport[transport["status"] == "TESTED"].copy()
    if len(tested):
        tested["composite_q"] = multipletests(tested["composite_pvalue"], method="fdr_bh")[1]
        tested["empirical_q"] = multipletests(tested["empirical_pvalue"], method="fdr_bh")[1]
        tested["composite_pass"] = tested["composite_q"] < 0.05
        tested["empirical_pass"] = tested["empirical_q"] < 0.05
        tested["both_pass"] = tested["composite_pass"] & tested["empirical_pass"]
        tested["both_pass_predicted"] = tested["both_pass"] & (tested["composite_slope"] > 0)
        transport = transport.set_index("unit")
        transport.loc[tested["unit"], ["composite_q", "empirical_q", "composite_pass",
                                       "empirical_pass", "both_pass", "both_pass_predicted"]] = \
            tested[["composite_q", "empirical_q", "composite_pass", "empirical_pass",
                    "both_pass", "both_pass_predicted"]].to_numpy()
        transport = transport.reset_index()
    return transport, pd.DataFrame(score_rows)


EXPR = {}
_ELIGIBILITY = None


def get_eligibility():
    global _ELIGIBILITY
    if _ELIGIBILITY is None:
        _ELIGIBILITY = deriv.discover_units()
    return _ELIGIBILITY


def run_variant1():
    print("VARIANT 1: cross-organ 3-fold", flush=True)
    rng = np.random.default_rng(SEED)
    all_transport = []
    for fold, held_out in enumerate(TISSUES, start=1):
        derive_tissues = [t for t in TISSUES if t != held_out]
        tag = f"fold{fold}-leave-{held_out}-out"
        print(f"  {tag}: deriving from {derive_tissues}", flush=True)
        genes, ranking = derive_signature(derive_tissues, tag=tag)
        genes.to_csv(RESULTS / f"followup-mature-heldout-{tag}-genes.csv", index=False)
        ranking.to_csv(RESULTS / f"followup-mature-heldout-{tag}-unit-ranking.csv", index=False)
        transport, scores = test_signature_in_tissue(genes, held_out, rng, tag)
        transport.insert(0, "fold", fold)
        transport.insert(1, "held_out_tissue", held_out)
        transport.insert(2, "derive_tissues", ";".join(derive_tissues))
        scores.insert(0, "fold", fold)
        scores.insert(1, "held_out_tissue", held_out)
        transport.to_csv(RESULTS / f"followup-mature-heldout-{tag}-transport.csv", index=False)
        scores.to_csv(RESULTS / f"followup-mature-heldout-{tag}-scores.csv", index=False)
        tested = transport[transport["status"] == "TESTED"]
        print(f"  {tag}: {int(tested['both_pass_predicted'].sum())} of {len(tested)} "
              f"pass both tests in the predicted direction; "
              f"{int(tested['both_pass'].sum())} pass both in either direction", flush=True)
        all_transport.append(transport)
    combined = pd.concat(all_transport, ignore_index=True)
    combined.to_csv(RESULTS / "followup-mature-heldout-cross-organ-transport.csv", index=False)
    return combined


def run_variant2():
    print("VARIANT 2: intermediate ages (existing P1m)", flush=True)
    p1m = pd.read_csv(RESULTS / "followup-mature-p1m-genes.csv")
    assert len(p1m) == 1914 and int((p1m["direction"] == "UP").sum()) == 1119
    assert int((p1m["direction"] == "DOWN").sum()) == 795
    up_genes = p1m.loc[p1m["direction"] == "UP", "mouse_gene_symbol"].tolist()
    down_genes = p1m.loc[p1m["direction"] == "DOWN", "mouse_gene_symbol"].tolist()
    eligibility, _, _ = get_eligibility()
    rows = []
    for tissue in TISSUES:
        genes, pairs, matrix = EXPR[tissue]
        units = eligibility[
            (eligibility["tissue"] == tissue) & (eligibility["eligible"] == "yes")
        ]["cell_ontology_class"].tolist()
        for ctype in units:
            sub = pairs[(pairs["cell_type"] == ctype) & (pairs["n_cells"] >= MIN_CELLS_PER_DONOR_CT)]
            if len(sub) == 0:
                continue
            scores, n_up, n_down = score_from_matrix(
                matrix[:, sub.index.to_numpy()], genes, up_genes, down_genes
            )
            if scores is None:
                continue
            for (mouse_id, age, sex, n_cells), score in zip(
                sub[["mouse_id", "age", "sex", "n_cells"]].itertuples(index=False), scores
            ):
                rows.append({"unit": f"{tissue}::{ctype}", "tissue": tissue,
                             "cell_type": ctype, "mouse_id": mouse_id, "age": age,
                             "age_months": AGE_MONTHS[age], "sex": sex,
                             "n_cells": int(n_cells), "score": float(score),
                             "n_up_detected": n_up, "n_down_detected": n_down})
    per_unit = pd.DataFrame(rows)
    per_unit.to_csv(RESULTS / "followup-mature-heldout-intermediate-ages-per-unit.csv", index=False)
    per_mouse = per_unit.groupby(["mouse_id", "age", "age_months", "sex"], as_index=False).agg(
        mean_score=("score", "mean"), n_units=("score", "size")
    )
    per_mouse.to_csv(RESULTS / "followup-mature-heldout-intermediate-ages-per-mouse.csv", index=False)
    print(per_mouse.sort_values(["age_months", "mouse_id"]).to_string(index=False), flush=True)
    age_means = per_mouse.groupby("age_months")["mean_score"].agg(["mean", "count"])
    print(age_means.to_string(), flush=True)
    from scipy.stats import spearmanr
    rho = spearmanr(per_mouse["age_months"], per_mouse["mean_score"])
    print(f"  Spearman score vs age-months: rho={rho.statistic:.3f} p={rho.pvalue:.3g}", flush=True)
    return per_unit, per_mouse


def run_variant3():
    print("VARIANT 3: donor split (descriptive only)", flush=True)
    # Hold out the highest-ID mouse per side: 18-M-53 and 30-M-5. The 18m
    # derive set keeps both sexes (2 female + 1 male); the 30m derive set
    # keeps 3 males. Eligibility is recomputed on the 3v3 subset.
    held_out = ["18-M-53", "30-M-5"]
    eligibility, source_mice, _ = get_eligibility()
    all_mice = sorted(source_mice["mouse.id"].astype(str).unique())
    allow = [m for m in all_mice if m not in held_out]
    print(f"  held out: {held_out}; derive mice: {allow}", flush=True)
    genes, ranking = derive_signature(TISSUES, allow_mice=allow, tag="donor-split")
    genes.to_csv(RESULTS / "followup-mature-heldout-donor-split-genes.csv", index=False)
    ranking.to_csv(RESULTS / "followup-mature-heldout-donor-split-unit-ranking.csv", index=False)
    up_genes = genes.loc[genes["direction"] == "UP", "mouse_gene_symbol"].tolist()
    down_genes = genes.loc[genes["direction"] == "DOWN", "mouse_gene_symbol"].tolist()
    rows = []
    for tissue in TISSUES:
        genes_arr, pairs, matrix = EXPR[tissue]
        units = ranking[ranking["tissue"] == tissue]["cell_ontology_class"].tolist()
        for ctype in units:
            sub = pairs[(pairs["cell_type"] == ctype)
                        & (pairs["age"].isin(deriv.PRIMARY_AGES))
                        & (pairs["n_cells"] >= MIN_CELLS_PER_DONOR_CT)]
            if len(sub) == 0:
                continue
            scores, n_up, n_down = score_from_matrix(
                matrix[:, sub.index.to_numpy()], genes_arr, up_genes, down_genes
            )
            if scores is None:
                continue
            frame = sub.copy()
            frame["score"] = scores
            for held in held_out:
                hrow = frame[frame["mouse_id"] == held]
                if len(hrow) == 0:
                    continue
                held_age = hrow["age"].iloc[0]
                held_score = float(hrow["score"].iloc[0])
                derive = frame[frame["mouse_id"] != held]
                mean_18 = float(derive.loc[derive["age"] == "18m", "score"].mean())
                mean_30 = float(derive.loc[derive["age"] == "30m", "score"].mean())
                # Direction check: a held-out 30m mouse should score above
                # the derive 18m mean; a held-out 18m mouse below the
                # derive 30m mean.
                if held_age == "30m":
                    match = held_score > mean_18
                else:
                    match = held_score < mean_30
                rows.append({"unit": f"{tissue}::{ctype}", "tissue": tissue,
                             "cell_type": ctype, "held_out_mouse": held,
                             "held_out_age": held_age, "held_out_score": held_score,
                             "derive_mean_18m": mean_18, "derive_mean_30m": mean_30,
                             "derive_n_18m": int((derive["age"] == "18m").sum()),
                             "derive_n_30m": int((derive["age"] == "30m").sum()),
                             "direction_match": bool(match)})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "followup-mature-heldout-donor-split.csv", index=False)
    print(f"  donor split: {int(out['direction_match'].sum())} of {len(out)} "
          f"held-out unit checks match the predicted direction", flush=True)
    return out


def main():
    print(f"Held-out validation started {utc_now()}", flush=True)
    for tissue in TISSUES:
        print(f"Building expression matrix: {tissue}", flush=True)
        EXPR[tissue] = build_expression_matrix(tissue)
    # Mouse census sanity check (all four ages).
    for tissue in TISSUES:
        _, pairs, _ = EXPR[tissue]
        census = pairs.groupby(["age", "sex"])["mouse_id"].nunique()
        print(f"  {tissue} mice by age/sex:\n{census.to_string()}", flush=True)
    t0 = time.time()
    transport = run_variant1()
    per_unit, per_mouse = run_variant2()
    split = run_variant3()
    tested = transport[transport["status"] == "TESTED"]
    summary = {
        "run_started_utc": utc_now(),
        "model": MODEL,
        "seed": SEED,
        "null_iterations": ITERATIONS,
        "cross_organ_units_tested": int(len(tested)),
        "cross_organ_both_pass_predicted": int(tested["both_pass_predicted"].sum()),
        "cross_organ_both_pass_either_direction": int(tested["both_pass"].sum()),
        "cross_organ_composite_pass": int(tested["composite_pass"].sum()),
        "cross_organ_empirical_pass": int(tested["empirical_pass"].sum()),
        "intermediate_age_mice": int(per_mouse["mouse_id"].nunique()),
        "donor_split_checks": int(len(split)),
        "donor_split_direction_matches": int(split["direction_match"].sum()),
        "facs_arm": "NOT RUN: no FACS files locally; TMS FACS has no 30m mice",
        "runtime_s": time.time() - t0,
    }
    with open(RESULTS / "followup-mature-heldout-summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(pd.Series(summary).to_string(), flush=True)
    print(f"Held-out validation finished {utc_now()}", flush=True)


if __name__ == "__main__":
    main()
