#!/usr/bin/env python3
"""Run P0_Raw257 and P1_Clean through the existing human transport tests."""

from __future__ import annotations

import gc
import os
import sys
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp
from anndata.io import read_elem
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_m10_pipeline import RDS, ols as m10_ols
from scripts.run_reviewer_phase1_signature_null import random_signature_scores as m10_random_scores
from scripts.run_rev2_transport_null import (
    fit_age_sex,
    random_signature_scores as mssm_random_scores,
)
from scripts.run_followup_hcata_pool_sensitivity import add_yazar_pool_dummies, load_cache
from scripts.run_followup_organ_17_heart_gao import (
    H5AD_PATH,
    STUDY_DUMMY_MAP,
    build_donor_audit,
    build_donor_meta_table,
    stream_pseudobulk_counts,
)
from src.aging_programs import load_program_tiers
from src.followup_muscle import matched_null_effect_pvalue
from src.followup_organs import (
    finalize_pseudobulk_group,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

RESULTS = REPO / "results"
REPORT = REPO / "reports" / "followup-11pt-block5-p0-transport.md"
TRANSPORT_CSV = RESULTS / "followup-11pt-block5-p0-transport.csv"
COVERAGE_CSV = RESULTS / "followup-11pt-block5-p0-coverage.csv"
YAZAR_CACHE = Path("/home/reza/projects/agescape-revision-data/tmp_organs/hcata/yazar_cache")
MSSM_H5AD = Path("/home/reza/projects/agescape-revision-data/mssm.h5ad")
ITERATIONS = 2000
SEED = 20261105
TIERS = ("P0_Raw257", "P1_Clean")


def verify_and_load_programs() -> tuple[dict[str, dict[str, list[str]]], pd.DataFrame]:
    hierarchy = pd.read_csv(RESULTS / "m5-program-hierarchy.csv")
    gene_level = pd.read_csv(RESULTS / "m5-gene-level-programs.csv")
    catalog = pd.read_csv(RESULTS / "m5-mouse-program-catalog.csv")
    assert len(gene_level) == int(catalog["n_de"].sum()) == 312
    assert ((gene_level["FDR"] < 0.05) & (gene_level["log2FC"].abs() >= 0.25)).all()
    assert gene_level["gene"].nunique() == 257
    p0 = hierarchy[hierarchy["tier"] == "P0_Raw257"]
    p1 = hierarchy[hierarchy["tier"] == "P1_Clean"]
    assert p0["gene"].nunique() == 257
    assert p1["gene"].nunique() == 139
    assert (p1["direction"] == "UP").sum() == 89
    assert (p1["direction"] == "DOWN").sum() == 50
    programs = load_program_tiers(REPO, "human")
    assert (len(programs["P0_Raw257"]["up"]), len(programs["P0_Raw257"]["down"])) == (120, 68)
    assert (len(programs["P1_Clean"]["up"]), len(programs["P1_Clean"]["down"])) == (76, 42)
    return programs, hierarchy


def build_coverage(hierarchy: pd.DataFrame) -> pd.DataFrame:
    p0 = hierarchy[hierarchy["tier"] == "P0_Raw257"].copy()
    orthologs = pd.read_csv(REPO / "data" / "reference" / "orthologs-mouse-human.csv")
    orthologs = orthologs[orthologs["orthology_type"] == "ortholog_one2one"]
    mapper = orthologs.drop_duplicates("mouse_gene_symbol").set_index("mouse_gene_symbol")["human_gene_symbol"]
    p0["mapped"] = p0["gene"].isin(mapper.index)
    p0["biotype_class"] = np.where(p0["biotype"] == "protein_coding", "protein_coding", "noncoding_or_unknown")
    rows = []
    for label, sub in [("all", p0), *list(p0.groupby("biotype_class", sort=True))]:
        mapped = sub[sub["mapped"]]
        rows.append({
            "biotype_class": label,
            "n_mouse_genes": int(sub["gene"].nunique()),
            "n_mapped": int(mapped["gene"].nunique()),
            "n_unmapped": int((~sub["mapped"]).sum()),
            "mapping_coverage": float(sub["mapped"].mean()),
            "n_mouse_up": int((sub["direction"] == "UP").sum()),
            "n_mouse_down": int((sub["direction"] == "DOWN").sum()),
            "n_human_up": int((mapped["direction"] == "UP").sum()),
            "n_human_down": int((mapped["direction"] == "DOWN").sum()),
            "unmapped_genes": ";".join(sorted(sub.loc[~sub["mapped"], "gene"])),
        })
    coverage = pd.DataFrame(rows)
    assert int(coverage.loc[coverage["biotype_class"] == "all", "n_mapped"].iloc[0]) == 188
    return coverage


def empty_row(cohort: str, unit: str, tier: str, source: str, model: str) -> dict[str, object]:
    programs = load_program_tiers(REPO, "human")
    program = programs[tier]
    return {
        "cohort": cohort,
        "unit": unit,
        "gene_set": tier,
        "arm": "primary",
        "source": source,
        "model": model,
        "composite_slope": np.nan,
        "composite_slope_per_year": np.nan,
        "age_scale": "per_year",
        "composite_pvalue": np.nan,
        "composite_q": np.nan,
        "empirical_pvalue": np.nan,
        "empirical_q": np.nan,
        "composite_pass": False,
        "empirical_pass": pd.NA,
        "both_pass": pd.NA,
        "n_donors": np.nan,
        "n_genes_mapped": len(set(program["up"]) | set(program["down"])),
        "n_up_mapped": len(program["up"]),
        "n_down_mapped": len(program["down"]),
        "n_genes_detected": np.nan,
        "n_up_detected": np.nan,
        "n_down_detected": np.nan,
        "null_iterations": np.nan,
        "random_seed": np.nan,
        "multiplicity_family": "",
    }


def build_m9_m10(programs: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    rows = []
    scores = pd.read_csv(RESULTS / "m9-donor-scores.csv")
    for tier in TIERS:
        m9 = pd.read_csv(RESULTS / ("m9-regressions-p0.csv" if tier == "P0_Raw257" else "m9-regressions.csv"))
        coverage_column = "coverage_p0" if tier == "P0_Raw257" else "coverage"
        coverage_by_cohort = scores.groupby("cohort")[coverage_column].first()
        for r in m9.itertuples(index=False):
            row = empty_row(r.cohort, r.cell_type, tier, "M9", "age_plus_sex")
            coverage = float(coverage_by_cohort[r.cohort])
            row.update({
                "composite_slope": r.slope_per_year,
                "composite_slope_per_year": r.slope_per_year,
                "composite_pvalue": r.p,
                "n_donors": r.n_donors,
                "n_genes_detected": int(round(coverage * row["n_genes_mapped"])),
                "n_up_detected": int(round(coverage * row["n_up_mapped"])),
                "n_down_detected": int(round(coverage * row["n_down_mapped"])),
                "multiplicity_family": "primary_M9_M10_cell_type_transport",
            })
            rows.append(row)
        m10 = pd.read_csv(RESULTS / ("m10-celltype-regressions-p0.csv" if tier == "P0_Raw257" else "m10-celltype-regressions.csv"))
        for r in m10.itertuples(index=False):
            row = empty_row("DLPFC_Lifespan_Atlas", r.cell_type, tier, "M10", "fully_adjusted")
            row.update({
                "composite_slope": r.slope_covar_adj,
                "composite_slope_per_year": r.slope_covar_adj,
                "composite_pvalue": r.p_covar_adj,
                "n_donors": r.n_donors,
                "n_genes_detected": r.n_up_genes + r.n_down_genes,
                "n_up_detected": r.n_up_genes,
                "n_down_detected": r.n_down_genes,
                "multiplicity_family": "primary_M9_M10_cell_type_transport",
            })
            rows.append(row)
    out = pd.DataFrame(rows)
    assert len(out) == 504
    for tier in TIERS:
        mask = out["gene_set"] == tier
        assert mask.sum() == 252
        out.loc[mask, "composite_q"] = multipletests(out.loc[mask, "composite_pvalue"], method="fdr_bh")[1]
    p1 = out[out["gene_set"] == "P1_Clean"].set_index(["cohort", "unit"])
    original = pd.read_csv(RESULTS / "human-transport-global-multiplicity.csv").set_index(["cohort", "cell_type"])
    original = original.rename_axis(["cohort", "unit"])
    joined = p1.join(original[["q_global_bh"]], how="left")
    assert np.allclose(joined["composite_q"], joined["q_global_bh"], rtol=1e-10, atol=1e-12)
    return out


def family_frame(df: pd.DataFrame, cohort: str, tier: str, arm: str, source: str) -> pd.DataFrame:
    rows = []
    for r in df.itertuples(index=False):
        row = empty_row(cohort, r.group, tier, source, "batch_adjusted" if arm == "primary" else "unadjusted")
        row.update({
            "arm": arm,
            "composite_slope": r.composite_effect,
            "composite_slope_per_year": r.composite_effect / 10.0,
            "age_scale": "per_decade",
            "composite_pvalue": r.composite_pvalue,
            "empirical_pvalue": r.empirical_null_pvalue,
            "n_donors": r.n_donors,
            "n_genes_detected": r.n_p1_total_detected,
            "n_up_detected": r.n_p1_up_detected,
            "n_down_detected": r.n_p1_down_detected,
            "null_iterations": r.null_iterations,
            "random_seed": SEED,
            "multiplicity_family": f"{cohort}_{arm}_{tier}",
        })
        rows.append(row)
    return pd.DataFrame(rows)


def run_gao(programs: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    print("Gao: reading metadata and streaming counts", flush=True)
    with h5py.File(H5AD_PATH, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
    audit = build_donor_audit(obs)
    donor_meta = build_donor_meta_table(audit)
    included = set(donor_meta.index)
    keep = obs["donor_id"].astype(str).isin(included).to_numpy()
    raw, counts, _ = stream_pseudobulk_counts(
        H5AD_PATH,
        keep,
        obs["donor_id"].astype(str).to_numpy(),
        obs["final_cell_type"].astype(str).to_numpy(),
        var.index.astype(str).to_numpy(),
    )
    expr = {}
    metas = {}
    for group, raw_group in raw.items():
        donors = list(raw_group.columns)
        cell_count = counts.loc[donors, "Total_Cells"] if group == "All_LV_Nuclei" else counts.loc[donors, group]
        expr[group], metas[group], _ = finalize_pseudobulk_group(raw_group, donor_meta.loc[donors], cell_count)
    covariates = ["sex_male", "log_lib", *STUDY_DUMMY_MAP.values()]
    p0_primary = run_family_c_p1_signature(
        expr, metas, p1_program=programs["P0_Raw257"], candidate_covariates=tuple(covariates),
        min_donors=4, min_residual_df=2, iterations=ITERATIONS, seed=SEED,
    )
    p0_unadjusted = run_family_c_p1_signature(
        expr, metas, p1_program=programs["P0_Raw257"], candidate_covariates=("sex_male", "log_lib"),
        min_donors=4, min_residual_df=2, iterations=ITERATIONS, seed=SEED,
    )
    p1_unadjusted = run_family_c_p1_signature(
        expr, metas, p1_program=programs["P1_Clean"], candidate_covariates=("sex_male", "log_lib"),
        min_donors=4, min_residual_df=2, iterations=ITERATIONS, seed=SEED,
    )
    p1_primary = pd.read_csv(RESULTS / "followup-organ-17-heart-gao-p1-signature.csv")
    parts = [
        family_frame(p0_primary, "Gao_Heart", "P0_Raw257", "primary", "Gao"),
        family_frame(p1_primary, "Gao_Heart", "P1_Clean", "primary", "Gao"),
        family_frame(p0_unadjusted, "Gao_Heart", "P0_Raw257", "unadjusted", "Gao"),
        family_frame(p1_unadjusted, "Gao_Heart", "P1_Clean", "unadjusted", "Gao"),
    ]
    out = pd.concat(parts, ignore_index=True)
    assert len(out) == 56
    del raw, expr, metas, obs, var
    gc.collect()
    return out


def yazar_pool_inputs():
    donor_meta, cell_counts, expr, metas = load_cache(YAZAR_CACHE, drop_count_cols=("peripheral blood mononuclear cell",))
    all_pools = sorted(donor_meta["pool_number"].dropna().astype(int).unique())
    all_pool_cols = [f"pool_{p}" for p in all_pools]
    fe_expr = {}
    fe_meta = {}
    for group, meta in metas.items():
        pool_counts = meta["pool_number"].value_counts()
        ge2 = pool_counts[pool_counts >= 2]
        keep = meta.index[meta["pool_number"].isin(ge2.index)].tolist()
        estimated_df = len(keep) - (4 + max(0, len(ge2) - 1))
        if len(ge2) < 2 or estimated_df < 10:
            continue
        sub = meta.loc[keep].copy()
        reference = int(ge2.index[0])
        for pool in all_pools:
            sub[f"pool_{pool}"] = (sub["pool_number"].astype(int) == pool).astype(float) if pool != reference else 0.0
        fe_expr[group] = expr[group][keep]
        fe_meta[group] = sub
    return donor_meta, cell_counts, expr, metas, fe_expr, fe_meta, all_pool_cols


def run_yazar(programs: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    print("Yazar: running pool-adjusted and unadjusted P0 transport", flush=True)
    _, _, expr, metas, fe_expr, fe_meta, pool_cols = yazar_pool_inputs()
    covariates = ("sex_male", "log_lib", *pool_cols)
    p0_primary = run_family_c_p1_signature(
        fe_expr, fe_meta, p1_program=programs["P0_Raw257"], candidate_covariates=covariates,
        min_donors=4, min_residual_df=10, iterations=ITERATIONS, seed=SEED,
    )
    p0_unadjusted = run_family_c_p1_signature(
        expr, metas, p1_program=programs["P0_Raw257"], candidate_covariates=("sex_male", "log_lib"),
        min_donors=4, min_residual_df=2, iterations=ITERATIONS, seed=SEED,
    )
    p1_pool = pd.read_csv(RESULTS / "followup-hcata-pool-sensitivity-p1-signature.csv")
    p1_primary = p1_pool[(p1_pool["study"] == "yazar") & (p1_pool["arm"] == "pool_fe")].copy()
    p1_unadjusted_all = pd.read_csv(RESULTS / "followup-hcata-yazar-p1-signature.csv")
    p1_unadjusted = p1_unadjusted_all[p1_unadjusted_all["arm"] == "primary"].copy()
    parts = [
        family_frame(p0_primary, "Yazar_Blood", "P0_Raw257", "primary", "Yazar"),
        family_frame(p1_primary, "Yazar_Blood", "P1_Clean", "primary", "Yazar"),
        family_frame(p0_unadjusted, "Yazar_Blood", "P0_Raw257", "unadjusted", "Yazar"),
        family_frame(p1_unadjusted, "Yazar_Blood", "P1_Clean", "unadjusted", "Yazar"),
    ]
    out = pd.concat(parts, ignore_index=True)
    assert len(out) == 74
    return out


def run_m10_nulls(programs: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    print("PsychAD: running 2000 matched nulls for P0 and P1", flush=True)
    assert RDS.exists(), f"Missing {RDS}"
    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    metadata = {
        name: np.asarray(cd.listData[name]).ravel()
        for name in ("Age", "Sex", "Brain_bank", "PMI", "TechPC1", "TechPC2", "TechPC3")
    }
    rows = []
    for tier in TIERS:
        program = programs[tier]
        actual = pd.read_csv(RESULTS / ("m10-celltype-regressions-p0.csv" if tier == "P0_Raw257" else "m10-celltype-regressions.csv")).set_index("cell_type")
        rng = np.random.default_rng(SEED)
        tier_rows = []
        for cell_type in sorted(obj.assays.data.listData.keys()):
            source = obj.assays.data.listData[cell_type]
            genes = np.array([str(gene) for gene in source.Dimnames[0]])
            mat = sp.csc_matrix(
                (np.asarray(source.x), np.asarray(source.i), np.asarray(source.p)),
                shape=tuple(source.Dim),
            )
            library = np.asarray(mat.sum(axis=0)).ravel()
            keep = library > 0
            mat = mat[:, keep]
            normalized = mat @ sp.diags(1e6 / library[keep])
            normalized = normalized.tocsc()
            normalized.data = np.log1p(normalized.data)
            lookup = {gene: i for i, gene in enumerate(genes)}
            up = np.array([lookup[g] for g in program["up"] if g in lookup], dtype=int)
            down = np.array([lookup[g] for g in program["down"] if g in lookup], dtype=int)
            available = np.ones(len(genes), dtype=bool)
            available[np.concatenate([up, down])] = False
            actual_scores = np.asarray(normalized[up, :].mean(axis=0)).ravel() - np.asarray(normalized[down, :].mean(axis=0)).ravel()
            age = metadata["Age"][keep].astype(float)
            design = np.column_stack([
                np.ones(len(age)), age, (metadata["Sex"][keep] == "Male").astype(float),
                (metadata["Brain_bank"][keep] == "HBCC").astype(float),
                metadata["PMI"][keep].astype(float), metadata["TechPC1"][keep].astype(float),
                metadata["TechPC2"][keep].astype(float), metadata["TechPC3"][keep].astype(float),
            ])
            observed_fit = m10_ols(actual_scores, design)
            expected = float(actual.loc[cell_type, "slope_covar_adj"])
            assert np.isclose(observed_fit["slope"], expected, rtol=1e-9, atol=1e-12)
            random_scores = m10_random_scores(
                {"matrix": normalized, "genes": genes}, up, down, ITERATIONS, rng, available
            )
            random_slopes = np.array([m10_ols(score, design)["slope"] for score in random_scores])
            tier_rows.append({
                "cohort": "DLPFC_Lifespan_Atlas",
                "unit": cell_type,
                "gene_set": tier,
                "empirical_pvalue": matched_null_effect_pvalue(observed_fit["slope"], random_slopes),
                "n_up_detected": len(up),
                "n_down_detected": len(down),
                "n_genes_detected": len(up) + len(down),
                "null_iterations": ITERATIONS,
                "random_seed": SEED,
            })
            print(f"  PsychAD {tier} {cell_type}: p={tier_rows[-1]['empirical_pvalue']:.6g}", flush=True)
        q = multipletests([r["empirical_pvalue"] for r in tier_rows], method="fdr_bh")[1]
        for row, qvalue in zip(tier_rows, q):
            row["empirical_q"] = qvalue
        rows.extend(tier_rows)
    out = pd.DataFrame(rows)
    assert len(out) == 54
    return out


def build_mssm_pseudobulk(cell_types: list[str], scores: pd.DataFrame):
    pairs = []
    for cell_type in cell_types:
        sub = scores[scores["cell_type"] == cell_type]
        donors = sorted(sub["donor_id"].astype(str))
        pairs.extend((cell_type, donor) for donor in donors)
    cache = Path("/home/reza/projects/agescape-revision-data/tmp_organs/block5_mssm_pseudobulk.npz")
    if cache.exists():
        print("MSSM: loading cached pseudobulk", flush=True)
        with np.load(cache, allow_pickle=False) as cached:
            cached_pairs = [tuple(row) for row in cached["pairs"].astype(str)]
            assert cached_pairs == pairs
            return cached["genes"].astype(str), pairs, cached["pseudobulk"]
    print("MSSM: streaming local h5ad once for all null units", flush=True)
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    pair_counts = np.zeros(len(pairs), dtype=np.int64)
    with h5py.File(MSSM_H5AD, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        genes = var["feature_name"].astype(str).to_numpy()
        donor_arr = obs["donor_id"].astype(str).to_numpy()
        ct_arr = obs["cell_type"].astype(str).to_numpy()
        primary = obs["is_primary_data"].astype(bool).to_numpy()
        pair_ids = np.full(len(obs), -1, dtype=np.int64)
        for i, (donor, cell_type) in enumerate(zip(donor_arr, ct_arr)):
            if primary[i]:
                pair_ids[i] = pair_index.get((cell_type, donor), -1)
        pair_sums = np.zeros((len(pairs), len(genes)), dtype=np.float64)
        group = handle["X"]
        indptr = group["indptr"][:]
        data_ds = group["data"]
        indices_ds = group["indices"]
        chunk = 2000
        for start in range(0, len(obs), chunk):
            stop = min(start + chunk, len(obs))
            local_pairs = pair_ids[start:stop]
            keep = local_pairs >= 0
            if not keep.any():
                continue
            lo, hi = int(indptr[start]), int(indptr[stop])
            matrix = sp.csr_matrix(
                (data_ds[lo:hi], indices_ds[lo:hi], indptr[start:stop + 1] - lo),
                shape=(stop - start, len(genes)),
            )[keep]
            local_pairs = local_pairs[keep]
            library_totals = np.asarray(matrix.sum(axis=1)).ravel()
            matrix = matrix.astype(np.float64)
            scale = 1e6 / np.where(library_totals > 0, library_totals, 1.0)
            matrix.data *= np.repeat(scale, np.diff(matrix.indptr))
            np.log1p(matrix.data, out=matrix.data)
            unique, inverse, cell_counts = np.unique(local_pairs, return_inverse=True, return_counts=True)
            indicator = sp.csr_matrix(
                (np.ones(len(local_pairs)), (inverse, np.arange(len(local_pairs)))),
                shape=(len(unique), len(local_pairs)),
            )
            pair_sums[unique, :] += (indicator @ matrix).toarray()
            np.add.at(pair_counts, unique, cell_counts)
            if start % 200000 == 0:
                print(f"  MSSM cells {start} to {stop} of {len(obs)}", flush=True)
    for i, (cell_type, donor) in enumerate(pairs):
        expected = int(scores[(scores["cell_type"] == cell_type) & (scores["donor_id"].astype(str) == donor)]["n_cells"].iloc[0])
        assert pair_counts[i] == expected, f"MSSM {cell_type} {donor}: {pair_counts[i]} cells, expected {expected}"
    pseudobulk = pair_sums / pair_counts[:, None]
    np.savez(cache, genes=np.asarray(genes, dtype=str), pairs=np.asarray(pairs, dtype=str), pseudobulk=pseudobulk)
    return genes, pairs, pseudobulk


def run_mssm_nulls(programs: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    m9_scores = pd.read_csv(RESULTS / "m9-donor-scores.csv")
    scores = m9_scores[m9_scores["cohort"] == "MSSM_DLPFC"].copy()
    scores["donor_id"] = scores["donor_id"].astype(str)
    cell_types = sorted(pd.read_csv(RESULTS / "m9-regressions.csv").query("cohort == 'MSSM_DLPFC'")["cell_type"])
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
            score_column = "score_p0" if tier == "P0_Raw257" else "score"
            max_diff = float(np.max(np.abs(actual_scores - meta[score_column].to_numpy(float))))
            fit = fit_age_sex(actual_scores, meta["age"].to_numpy(float), meta["sex"].astype(str).to_numpy())
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
                "empirical_pvalue": matched_null_effect_pvalue(fit["slope"], random_slopes),
                "n_up_detected": len(up),
                "n_down_detected": len(down),
                "n_genes_detected": len(up) + len(down),
                "null_iterations": ITERATIONS,
                "random_seed": SEED,
                "max_abs_score_diff_vs_m9": max_diff,
            })
            print(f"  MSSM {tier} {cell_type}: p={tier_rows[-1]['empirical_pvalue']:.6g}", flush=True)
        q = multipletests([r["empirical_pvalue"] for r in tier_rows], method="fdr_bh")[1]
        for row, qvalue in zip(tier_rows, q):
            row["empirical_q"] = qvalue
        rows.extend(tier_rows)
    out = pd.DataFrame(rows)
    assert len(out) == 42
    return out


def finalize_transport(base: pd.DataFrame, followups: pd.DataFrame, nulls: pd.DataFrame) -> pd.DataFrame:
    base = base.copy()
    null_columns = [
        "cohort", "unit", "gene_set", "empirical_pvalue", "empirical_q",
        "n_up_detected", "n_down_detected", "n_genes_detected", "null_iterations", "random_seed",
    ]
    updates = nulls[null_columns]
    base = base.merge(updates, on=["cohort", "unit", "gene_set"], how="left", suffixes=("", "_new"), validate="one_to_one")
    for column in null_columns[3:]:
        mask = base[f"{column}_new"].notna()
        base.loc[mask, column] = base.loc[mask, f"{column}_new"]
        base = base.drop(columns=[f"{column}_new"])
    out = pd.concat([base, followups], ignore_index=True)
    for (_, _, _), idx in out[out["source"].isin(["Gao", "Yazar"])].groupby(["cohort", "arm", "gene_set"]).groups.items():
        out.loc[idx, "composite_q"] = multipletests(out.loc[idx, "composite_pvalue"], method="fdr_bh")[1]
        out.loc[idx, "empirical_q"] = multipletests(out.loc[idx, "empirical_pvalue"], method="fdr_bh")[1]
    out["composite_pass"] = out["composite_q"] < 0.05
    out["empirical_pass"] = out["empirical_q"].lt(0.05).astype("boolean")
    out["both_pass"] = out["composite_pass"].astype("boolean") & out["empirical_pass"]
    out = out.sort_values(["cohort", "unit", "arm", "gene_set"]).reset_index(drop=True)
    primary = out[out["arm"] == "primary"]
    assert len(primary) == 566
    assert primary.groupby("gene_set").size().to_dict() == {"P0_Raw257": 283, "P1_Clean": 283}
    assert primary["empirical_pvalue"].notna().sum() == 158
    return out


def fmt(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.4g}"


def write_report(transport: pd.DataFrame, coverage: pd.DataFrame) -> str:
    primary = transport[transport["arm"] == "primary"].copy()
    key = ["cohort", "unit"]
    p0 = primary[primary["gene_set"] == "P0_Raw257"].set_index(key)
    p1 = primary[primary["gene_set"] == "P1_Clean"].set_index(key)
    assert p0.index.equals(p1.index)
    changed = p0["both_pass"].fillna(False) & (~p1["composite_pass"]) & (~p1["empirical_pass"].fillna(True))
    coverage_fraction = float(coverage.loc[coverage["biotype_class"] == "all", "mapping_coverage"].iloc[0])
    if changed.any():
        verdict = "PROVES"
        reason = "P0 passes both tests in at least one unit where P1 fails both: " + ", ".join([f"{a} {b}" for a, b in changed[changed].index]) + "."
    elif coverage_fraction < 0.5 and not p0["both_pass"].fillna(False).any():
        verdict = "INCONCLUSIVE"
        reason = "P0 mapping coverage is below 50 percent and P0 passes nowhere, so filtering cannot be separated from mapping loss."
    else:
        verdict = "DOES NOT PROVE"
        reason = "P0 does not pass both tests in any unit where P1 fails both."

    lines = [
        "P0 and P1 were evaluated in 283 primary units each, including 252 M9/M10 units, 14 Gao units, and 17 Yazar units. The primary output also contains 34 unadjusted Gao/Yazar sensitivity rows per gene set.",
        "",
        "### Mapping coverage",
        "",
        "| Biotype class | Mouse genes | Mapped | Unmapped | Coverage | Mapped up | Mapped down |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in coverage.itertuples(index=False):
        lines.append(f"| {r.biotype_class} | {r.n_mouse_genes} | {r.n_mapped} | {r.n_unmapped} | {100*r.mapping_coverage:.2f}% | {r.n_human_up} | {r.n_human_down} |")
    lines.extend([
        "",
        "All 30 P0 genes with noncoding or unknown biotypes lack a one-to-one human mapping. The mapped P0 set has 120 up and 68 down human genes, 188 unique genes total. P1 has 76 up and 42 down human genes, 118 total.",
        "",
        "### Global comparison",
        "",
        "Nominal hit means composite p < 0.05. Composite pass means BH q < 0.05 in the stated composite family. Random pass means empirical BH q < 0.05 among the 79 units per gene set with an existing P1 null family. Both pass requires both q-values below 0.05.",
        "",
        "| Gene set | Primary units | Nominal composite hits | Composite passes | Null evaluated | Random passes | Both pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for tier in TIERS:
        sub = primary[primary["gene_set"] == tier]
        lines.append(
            f"| {tier} | {len(sub)} | {int((sub['composite_pvalue'] < 0.05).sum())} | "
            f"{int(sub['composite_pass'].sum())} | {int(sub['empirical_pvalue'].notna().sum())} | "
            f"{int(sub['empirical_pass'].fillna(False).sum())} | {int(sub['both_pass'].fillna(False).sum())} |"
        )
    lines.extend([
        "",
        "### Cohort comparison",
        "",
        "| Cohort | Gene set | Units | Nominal composite | Composite q passes | Null evaluated | Random q passes | Both pass |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for (cohort, tier), sub in primary.groupby(["cohort", "gene_set"], sort=True):
        lines.append(
            f"| {cohort} | {tier} | {len(sub)} | {int((sub['composite_pvalue'] < 0.05).sum())} | "
            f"{int(sub['composite_pass'].sum())} | {int(sub['empirical_pvalue'].notna().sum())} | "
            f"{int(sub['empirical_pass'].fillna(False).sum())} | {int(sub['both_pass'].fillna(False).sum())} |"
        )
    interesting = primary["both_pass"].fillna(False) | primary["composite_pass"]
    interesting_idx = set(zip(primary.loc[interesting, "cohort"], primary.loc[interesting, "unit"]))
    interesting_idx.update(changed[changed].index)
    lines.extend([
        "",
        "### Side-by-side units with a corrected composite pass or a verdict change",
        "",
        "The required CSV contains the complete side-by-side table for every unit. Slopes below are per year for M9/M10 and per decade for Gao/Yazar, matching each source pipeline.",
        "",
        "| Cohort | Unit | P0 slope | P0 composite q | P0 empirical q | P0 both | P1 slope | P1 composite q | P1 empirical q | P1 both |",
        "| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |",
    ])
    for index in sorted(interesting_idx):
        a = p0.loc[index]
        b = p1.loc[index]
        lines.append(
            f"| {index[0]} | {index[1]} | {fmt(a['composite_slope'])} | {fmt(a['composite_q'])} | {fmt(a['empirical_q'])} | {a['both_pass']} | "
            f"{fmt(b['composite_slope'])} | {fmt(b['composite_q'])} | {fmt(b['empirical_q'])} | {b['both_pass']} |"
        )
    lines.extend([
        "",
        "### Unadjusted sensitivity",
        "",
        "| Cohort | Gene set | Units | Composite q passes | Random q passes | Both pass |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ])
    sensitivity = transport[transport["arm"] == "unadjusted"]
    for (cohort, tier), sub in sensitivity.groupby(["cohort", "gene_set"], sort=True):
        lines.append(
            f"| {cohort} | {tier} | {len(sub)} | {int(sub['composite_pass'].sum())} | "
            f"{int(sub['empirical_pass'].sum())} | {int(sub['both_pass'].sum())} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        f"**{verdict}**. {reason}",
        "",
        "This verdict concerns the specified P0/P1 filter comparison, not universal cross-species transport. Empirical results are confined to the cohorts and units that already had P1 matched-null families; no new empirical family was defined for the other M9 units.",
        "",
        "## Outputs and reproducibility",
        "",
        "- `results/followup-11pt-block5-p0-transport.csv`: complete P0/P1 primary and unadjusted sensitivity comparison.",
        "- `results/followup-11pt-block5-p0-coverage.csv`: mapping coverage overall and by biotype class.",
        "- `scripts/run_followup_11pt_block5_p0_transport.py`: verification, transport, matched nulls, and report generation.",
        "- Fixed seed: 20261105. Random iterations: 2000 per evaluated unit and gene set.",
    ])
    spec = REPORT.read_text(encoding="utf-8").split("## Results", 1)[0]
    write_clean_markdown(spec + "## Results\n\n" + "\n".join(lines) + "\n", REPORT)
    return verdict


def main() -> int:
    programs, hierarchy = verify_and_load_programs()
    coverage = build_coverage(hierarchy)
    write_clean_csv(coverage, COVERAGE_CSV)
    base = build_m9_m10(programs)
    gao = run_gao(programs)
    yazar = run_yazar(programs)
    followups = pd.concat([gao, yazar], ignore_index=True)
    nulls = pd.concat([run_mssm_nulls(programs), run_m10_nulls(programs)], ignore_index=True)
    transport = finalize_transport(base, followups, nulls)
    write_clean_csv(transport, TRANSPORT_CSV)
    verdict = write_report(transport, coverage)
    print(f"Block 5 verdict: {verdict}", flush=True)
    print(f"Wrote {TRANSPORT_CSV} ({len(transport)} rows)", flush=True)
    print(f"Wrote {COVERAGE_CSV} ({len(coverage)} rows)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
