#!/usr/bin/env python3
"""
Derive mature P0m and P1m from exactly 18m versus 30m TMS mice.

The 30m group is all male, verified directly from every per-tissue h5ad input.
The design is ~ sex + age_group only when sex varies and the design matrix is
full rank. Otherwise sex is dropped and the reason is recorded. The age
contrast is always ["age_group", "30m", "18m"], so positive log2FC means up
with age. Animals aged 21m or 24m are excluded before unit discovery,
pseudobulk aggregation, and model fitting.

P0m preserves the original M2 independent-filtering DE rule. P1m preserves the
original uniform-BH, floor-10, clean-gene, unit-pooling, and mean-log2FC
 direction rules. P2/P3 sex-dimorphic and immediate-early filters do not apply.

Model used for this implementation and report: kimi-code/k3.1-pro.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "42")

import re
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.set_seeds import set_seeds

DATA_DIR = REPO / "data" / "interim"
REF_DIR = REPO / "data" / "reference"
RESULTS = REPO / "results"
RESULTS.mkdir(exist_ok=True)

REF_LEVEL = "18m"
TARGET_LEVEL = "30m"
PRIMARY_AGES = [REF_LEVEL, TARGET_LEVEL]
EXCLUDED_AGES = ["21m", "24m"]
CONTRAST = ["age_group", TARGET_LEVEL, REF_LEVEL]
FDR_THRESHOLD = 0.05
FC_THRESHOLD = 0.25
MIN_MICE = 3
MIN_CELLS = 50
FLOOR10_CELLS = 10
POSITIVE_CONTROLS = ["Cdkn2a", "Spp1", "Apoe", "Igfbp3", "Lgals3"]
MODEL = "kimi-code/k3.1-pro"

set_seeds(42)
warnings.filterwarnings("ignore", category=FutureWarning)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def safe_unit_name(unit):
    return unit.replace("::", "__").replace(" ", "_").replace(",", "").replace("-", "_")


def ineligibility_reason(row):
    reasons = []
    for age in PRIMARY_AGES:
        n_mice = int(row[f"n_mice_{age}"])
        n_cells = int(row[f"n_cells_{age}"])
        if n_mice == 0 and n_cells == 0:
            reasons.append(f"no {age} mice or cells")
        else:
            if n_mice < MIN_MICE:
                reasons.append(f"{age} has {n_mice} mice, requires at least {MIN_MICE}")
            if n_cells < MIN_CELLS:
                reasons.append(f"{age} has {n_cells} cells, requires at least {MIN_CELLS}")
    return "; ".join(reasons) if reasons else "eligible"


def discover_units():
    rows = []
    mouse_rows = []
    paths = sorted(DATA_DIR.glob("m1-qc-cleaned-*.h5ad"))
    if not paths:
        raise FileNotFoundError("No per-tissue m1-qc-cleaned-*.h5ad inputs found")
    for path in paths:
        adata = ad.read_h5ad(path, backed="r")
        required = {"age", "sex", "mouse.id", "tissue", "cell_ontology_class"}
        missing = required.difference(adata.obs.columns)
        if missing:
            raise KeyError(f"{path.name} lacks required obs columns: {sorted(missing)}")
        primary_obs = adata.obs[
            adata.obs["age"].isin(PRIMARY_AGES)
            & adata.obs["cell_ontology_class"].notna()
        ].copy()
        if not set(primary_obs["age"].unique()).issubset(set(PRIMARY_AGES)):
            raise AssertionError(f"Age leakage in {path.name}")
        mouse_meta = primary_obs[["mouse.id", "age", "sex"]].drop_duplicates()
        counts_per_mouse = mouse_meta.groupby("mouse.id", observed=True)[["age", "sex"]].nunique()
        if (counts_per_mouse > 1).any().any():
            raise AssertionError(f"Inconsistent mouse age or sex in {path.name}")
        for record in mouse_meta.drop_duplicates("mouse.id").to_dict("records"):
            mouse_rows.append({**record, "source_file": path.name})
        for (tissue, ctype), df in primary_obs.groupby(
            ["tissue", "cell_ontology_class"], observed=True
        ):
            row = {
                "tissue": str(tissue),
                "cell_ontology_class": str(ctype),
                "unit": f"{tissue}::{ctype}",
                "source_file": path.name,
            }
            for age in PRIMARY_AGES:
                age_df = df[df["age"] == age]
                row[f"n_cells_{age}"] = len(age_df)
                row[f"n_mice_{age}"] = age_df["mouse.id"].nunique()
                row[f"n_{age}_female_mice"] = age_df.loc[
                    age_df["sex"] == "female", "mouse.id"
                ].nunique()
                row[f"n_{age}_male_mice"] = age_df.loc[
                    age_df["sex"] == "male", "mouse.id"
                ].nunique()
                cells_per_mouse = age_df.groupby("mouse.id", observed=True).size()
                row[f"mice_ge10_{age}"] = int((cells_per_mouse >= FLOOR10_CELLS).sum())
            row["eligible"] = "yes" if all(
                row[f"n_mice_{age}"] >= MIN_MICE and row[f"n_cells_{age}"] >= MIN_CELLS
                for age in PRIMARY_AGES
            ) else "no"
            row["ineligibility_reason"] = ineligibility_reason(row)
            row["p1m_floor10_eligible"] = bool(
                row["eligible"] == "yes"
                and all(row[f"mice_ge10_{age}"] >= MIN_MICE for age in PRIMARY_AGES)
            )
            rows.append(row)
        adata.file.close()
    eligibility = pd.DataFrame(rows)
    eligibility = eligibility.sort_values(
        ["eligible", "n_cells_30m", "unit"], ascending=[False, False, True]
    ).reset_index(drop=True)
    source_mice = pd.DataFrame(mouse_rows).drop_duplicates(["mouse.id", "age", "sex"])
    sex_30m = set(source_mice.loc[source_mice["age"] == "30m", "sex"].astype(str))
    if sex_30m != {"male"}:
        raise AssertionError(f"Expected all-male 30m source group, observed {sorted(sex_30m)}")
    return eligibility, source_mice, len(paths)


def choose_design(metadata):
    design = "~ sex + age_group"
    sex_status = "kept"
    sex_reason = "sex varies and the intercept, sex, and age columns are full rank"
    if metadata["sex"].nunique() < 2:
        design = "~ age_group"
        sex_status = "dropped"
        sex_reason = "sex has no variation"
    else:
        design_matrix = np.column_stack(
            [
                np.ones(len(metadata)),
                (metadata["sex"].astype(str) == "male").astype(int),
                (metadata["age_group"].astype(str) == TARGET_LEVEL).astype(int),
            ]
        )
        if np.linalg.matrix_rank(design_matrix) < design_matrix.shape[1]:
            design = "~ age_group"
            sex_status = "dropped"
            sex_reason = "sex is collinear with age group"
    return design, sex_status, sex_reason


def independent_uniform_bh(pvalues):
    adjusted = np.full(len(pvalues), np.nan, dtype=float)
    valid = pvalues.notna().to_numpy()
    if valid.any():
        adjusted[valid] = multipletests(
            pvalues.loc[valid].to_numpy(), method="fdr_bh"
        )[1]
    return adjusted


def run_de(counts, metadata):
    keep_genes = counts.columns[counts.sum(axis=0) >= 10]
    if len(keep_genes) == 0:
        raise ValueError("No genes pass the total-count filter of at least 10")
    sub_counts = counts[keep_genes].astype(int)
    design, sex_status, sex_reason = choose_design(metadata)
    dds = DeseqDataSet(
        counts=sub_counts,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        n_cpus=4,
        quiet=True,
    )
    dds.deseq2()

    stat_p0 = DeseqStats(
        dds,
        contrast=CONTRAST,
        independent_filter=True,
        n_cpus=4,
        quiet=True,
    )
    stat_p0.summary()
    result = stat_p0.results_df.copy().rename(
        columns={"log2FoldChange": "log2FC", "padj": "FDR"}
    )

    stat_p1 = DeseqStats(
        dds,
        contrast=CONTRAST,
        independent_filter=False,
        n_cpus=4,
        quiet=True,
    )
    stat_p1.summary()
    uniform_result = stat_p1.results_df.copy()
    if not result.index.equals(uniform_result.index):
        raise AssertionError("P0 and P1 PyDESeq2 result indexes differ")
    if not np.allclose(
        result["log2FC"].to_numpy(),
        uniform_result["log2FoldChange"].to_numpy(),
        equal_nan=True,
    ):
        raise AssertionError("P0 and P1 log2FC values differ")
    if not np.allclose(
        result["pvalue"].to_numpy(), uniform_result["pvalue"].to_numpy(), equal_nan=True
    ):
        raise AssertionError("P0 and P1 raw p-values differ")

    result["FDR_uniform_pydeseq2"] = uniform_result["padj"]
    result["FDR_uniform"] = independent_uniform_bh(result["pvalue"])
    reported = result["FDR_uniform_pydeseq2"]
    recomputed = result["FDR_uniform"]
    if not (reported.isna() == recomputed.isna()).all():
        raise AssertionError("PyDESeq2 and recomputed BH missing-value patterns differ")
    valid = reported.notna()
    max_discrepancy = float(
        np.max(np.abs(reported.loc[valid] - recomputed.loc[valid]))
    ) if valid.any() else 0.0
    result["BH_abs_discrepancy"] = np.abs(reported - recomputed)
    result["P0m_significant"] = (
        (result["FDR"] < FDR_THRESHOLD)
        & (result["log2FC"].abs() >= FC_THRESHOLD)
    )
    result["P1m_significant_before_filters"] = (
        (result["FDR_uniform"] < FDR_THRESHOLD)
        & (result["log2FC"].abs() >= FC_THRESHOLD)
    )
    if "lfcSE" in result.columns:
        result["log2FC_ci_low"] = result["log2FC"] - 1.96 * result["lfcSE"]
        result["log2FC_ci_high"] = result["log2FC"] + 1.96 * result["lfcSE"]
    result["design_actual"] = design
    result["sex_term_status"] = sex_status
    result["sex_term_reason"] = sex_reason
    result["age_contrast"] = "30m_vs_18m"
    result.index.name = "index"

    size_factors = np.asarray(dds.obs["size_factors"], dtype=float).reshape(-1)
    normalized = sub_counts.to_numpy(dtype=float) / size_factors[:, None]
    normalized_means = {}
    for gene in POSITIVE_CONTROLS:
        if gene not in sub_counts.columns:
            continue
        gene_idx = sub_counts.columns.get_loc(gene)
        means = {}
        for age in PRIMARY_AGES:
            mask = (metadata["age_group"].astype(str) == age).to_numpy()
            means[age] = float(normalized[mask, gene_idx].mean())
        normalized_means[gene] = means
    diagnostics = {
        "design_actual": design,
        "sex_term_status": sex_status,
        "sex_term_reason": sex_reason,
        "n_tested_genes": len(result),
        "n_bh_pvalues": int(valid.sum()),
        "bh_max_abs_discrepancy": max_discrepancy,
        "bh_pass": bool(max_discrepancy <= 1e-10),
    }
    return result, diagnostics, normalized_means


def pseudobulk_unit(adata, tissue, ctype):
    mask = (
        (adata.obs["tissue"] == tissue)
        & (adata.obs["cell_ontology_class"] == ctype)
        & (adata.obs["age"].isin(PRIMARY_AGES))
    )
    sub = adata[mask].to_memory()
    if not set(sub.obs["age"].unique()).issubset(set(PRIMARY_AGES)):
        raise AssertionError(f"Age leakage in {tissue}::{ctype}")
    if sub.obs["mouse.id"].isna().any() or sub.obs["sex"].isna().any():
        raise ValueError(f"Missing mouse or sex metadata in {tissue}::{ctype}")
    mouse_ids = []
    metadata_rows = []
    bulk_rows = []
    for mouse_id, mouse_df in sub.obs.groupby("mouse.id", observed=True):
        if mouse_df["age"].nunique() != 1 or mouse_df["sex"].nunique() != 1:
            raise AssertionError(f"Inconsistent mouse metadata in {tissue}::{ctype}")
        cell_idx = np.where((sub.obs["mouse.id"] == mouse_id).to_numpy())[0]
        counts = np.asarray(sub.X[cell_idx].sum(axis=0)).reshape(-1)
        if (counts < 0).any() or not np.allclose(counts, np.rint(counts)):
            raise ValueError(f"Non-count expression values in {tissue}::{ctype}")
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
    metadata["age_group"] = pd.Categorical(
        metadata["age_group"], categories=PRIMARY_AGES
    )
    metadata["sex"] = pd.Categorical(metadata["sex"], categories=["female", "male"])
    if set(metadata["age_group"].astype(str)) != set(PRIMARY_AGES):
        raise AssertionError(f"Missing age group in {tissue}::{ctype}")
    return counts_df, metadata, len(sub)


def is_p1_clean(gene, unit, biotypes):
    gene = str(gene)
    if re.match(r"^Gm\d+$", gene):
        return False
    if re.search(r"-ps\d*$", gene):
        return False
    if re.match(r"^[0-9A-Z]+Rik$", gene):
        return False
    if gene.startswith(("Rpl", "Rps", "mt-", "Mrp", "MT-")) or gene == "Rn45s":
        return False
    if biotypes.get(gene, "unknown") != "protein_coding":
        return False
    if (
        "Lung" in unit
        and any(
            immune in unit
            for immune in ["B cell", "CD4", "NK", "monocyte", "macrophage"]
        )
        and gene in ["Mgp", "Col3a1"]
    ):
        return False
    if "Kidney" in unit and gene in ["Ins1", "Ins2"]:
        return False
    return True


def aggregate_genes(entries, fdr_column, ortholog_map, biotypes, chromosomes):
    if not entries:
        return pd.DataFrame(
            columns=[
                "mouse_gene_symbol", "human_ortholog", "direction", "n_units_significant",
                "n_units_up", "n_units_down", "mean_log2FC", "min_FDR",
                "max_abs_log2FC", "units_list", "biotype", "chromosome",
            ]
        )
    frame = pd.DataFrame(entries)
    rows = []
    for gene, group in frame.groupby("mouse_gene_symbol", sort=True):
        units = sorted(group["unit"].unique())
        mean_lfc = float(group["log2FC"].mean())
        rows.append(
            {
                "mouse_gene_symbol": gene,
                "human_ortholog": ortholog_map.get(gene, pd.NA),
                "direction": "UP" if mean_lfc > 0 else "DOWN",
                "n_units_significant": len(units),
                "n_units_up": int((group["log2FC"] > 0).sum()),
                "n_units_down": int((group["log2FC"] < 0).sum()),
                "mean_log2FC": mean_lfc,
                "min_FDR": float(group[fdr_column].min()),
                "max_abs_log2FC": float(group["log2FC"].abs().max()),
                "units_list": ";".join(units),
                "biotype": biotypes.get(gene, "unknown"),
                "chromosome": chromosomes.get(gene, "unknown"),
            }
        )
    return pd.DataFrame(rows)


def direction_counts(frame):
    if frame.empty:
        return {"total": 0, "up": 0, "down": 0}
    return {
        "total": len(frame),
        "up": int((frame["direction"] == "UP").sum()),
        "down": int((frame["direction"] == "DOWN").sum()),
    }


def main():
    start_time = utc_now()
    start_clock = time.time()
    print(f"Mature P1m/P0m derivation started {start_time}", flush=True)

    eligibility, source_mice, n_input_files = discover_units()
    eligibility_path = RESULTS / "followup-mature-p1m-eligibility.csv"
    eligibility.to_csv(eligibility_path, index=False)
    eligible = eligibility[eligibility["eligible"] == "yes"].copy()
    print(
        f"Evaluated {len(eligibility)} units; {len(eligible)} are eligible; "
        f"{int(eligibility['p1m_floor10_eligible'].sum())} pass the P1 floor-10 rule",
        flush=True,
    )
    print("Source mouse sex counts by age:", flush=True)
    print(
        source_mice.drop_duplicates("mouse.id")
        .groupby(["age", "sex"], observed=True)["mouse.id"]
        .nunique()
        .to_string(),
        flush=True,
    )

    chrom_df = pd.read_csv(REF_DIR / "mouse-gene-chromosomes.csv").set_index(
        "mouse_gene_symbol"
    )
    biotypes = chrom_df["gene_biotype"].to_dict()
    chromosomes = chrom_df["chromosome_name"].to_dict()
    orthologs = pd.read_csv(REF_DIR / "orthologs-mouse-human.csv")
    if set(orthologs["orthology_type"]) != {"ortholog_one2one"}:
        raise ValueError("Ortholog table is not exclusively one-to-one")
    if orthologs["mouse_gene_symbol"].duplicated().any():
        raise ValueError("Ortholog table has duplicate mouse symbols")
    ortholog_map = dict(
        zip(orthologs["mouse_gene_symbol"], orthologs["human_gene_symbol"])
    )

    de_results = {}
    ranking_rows = []
    bh_rows = []
    used_mouse_rows = []
    orientation_rows = []
    p0_entries = []
    p1_entries = []

    for tissue, tissue_units in eligible.groupby("tissue", sort=False):
        source_file = tissue_units["source_file"].iloc[0]
        adata = ad.read_h5ad(DATA_DIR / source_file, backed="r")
        for row in tissue_units.itertuples(index=False):
            unit = row.unit
            unit_start = time.time()
            counts, metadata, n_cells = pseudobulk_unit(
                adata, row.tissue, row.cell_ontology_class
            )
            result, diagnostics, normalized_means = run_de(counts, metadata)
            result.insert(0, "unit", unit)
            result.insert(1, "tissue", row.tissue)
            result.insert(2, "cell_ontology_class", row.cell_ontology_class)
            de_results[unit] = result
            unit_path = RESULTS / f"followup-mature-p1m-de_{safe_unit_name(unit)}.csv"
            result.to_csv(unit_path)
            for mouse_id, mouse in metadata.iterrows():
                used_mouse_rows.append(
                    {
                        "mouse.id": mouse_id,
                        "age": str(mouse["age_group"]),
                        "sex": str(mouse["sex"]),
                        "unit": unit,
                    }
                )
            n_p0 = int(result["P0m_significant"].sum())
            n_p1_before = int(result["P1m_significant_before_filters"].sum())
            ranking_rows.append(
                {
                    "unit": unit,
                    "tissue": row.tissue,
                    "cell_ontology_class": row.cell_ontology_class,
                    "n_mice_total": len(metadata),
                    "n_cells_total": n_cells,
                    "n_de_p0m": n_p0,
                    "n_de_uniform_before_clean": n_p1_before,
                    "p1m_floor10_eligible": bool(row.p1m_floor10_eligible),
                    **diagnostics,
                    "runtime_s": time.time() - unit_start,
                }
            )
            bh_rows.append({"unit": unit, **diagnostics})
            for gene, means in normalized_means.items():
                if gene not in result.index:
                    continue
                lfc = float(result.loc[gene, "log2FC"])
                if not np.isfinite(lfc) or lfc == 0:
                    continue
                unadjusted_lfc = float(
                    np.log2((means[TARGET_LEVEL] + 1e-12) / (means[REF_LEVEL] + 1e-12))
                )
                orientation_rows.append(
                    {
                        "unit": unit,
                        "gene": gene,
                        "log2FC": lfc,
                        "mean_normalized_18m": means[REF_LEVEL],
                        "mean_normalized_30m": means[TARGET_LEVEL],
                        "unadjusted_normalized_log2FC": unadjusted_lfc,
                        "sign_matches_normalized_means": bool(np.sign(lfc) == np.sign(unadjusted_lfc)),
                    }
                )
            sig_p0 = result[result["P0m_significant"]]
            for gene, gene_row in sig_p0.iterrows():
                p0_entries.append(
                    {
                        "unit": unit,
                        "mouse_gene_symbol": gene,
                        "log2FC": gene_row["log2FC"],
                        "FDR": gene_row["FDR"],
                        "pvalue": gene_row["pvalue"],
                    }
                )
            if row.p1m_floor10_eligible:
                sig_p1 = result[result["P1m_significant_before_filters"]]
                for gene, gene_row in sig_p1.iterrows():
                    if is_p1_clean(gene, unit, biotypes):
                        p1_entries.append(
                            {
                                "unit": unit,
                                "mouse_gene_symbol": gene,
                                "log2FC": gene_row["log2FC"],
                                "FDR_uniform": gene_row["FDR_uniform"],
                                "pvalue": gene_row["pvalue"],
                            }
                        )
            print(
                f"{unit}: P0m DE={n_p0}, uniform-BH DE before clean={n_p1_before}, "
                f"design={diagnostics['design_actual']}, "
                f"BH max diff={diagnostics['bh_max_abs_discrepancy']:.3g}, "
                f"runtime={time.time() - unit_start:.1f}s",
                flush=True,
            )
        adata.file.close()

    if not de_results:
        raise RuntimeError("No eligible unit models were fit")
    ranking = pd.DataFrame(ranking_rows).sort_values(
        ["n_de_p0m", "unit"], ascending=[False, True]
    )
    ranking.to_csv(RESULTS / "followup-mature-p1m-unit-ranking.csv", index=False)
    bh_check = pd.DataFrame(bh_rows)
    bh_check.to_csv(RESULTS / "followup-mature-p1m-bh-check.csv", index=False)

    design_columns = ranking.set_index("unit")[
        ["design_actual", "sex_term_status", "sex_term_reason"]
    ]
    for column in design_columns:
        eligibility[column] = eligibility["unit"].map(design_columns[column])
    eligibility.to_csv(eligibility_path, index=False)

    p0m = aggregate_genes(p0_entries, "FDR", ortholog_map, biotypes, chromosomes)
    p1m = aggregate_genes(p1_entries, "FDR_uniform", ortholog_map, biotypes, chromosomes)
    p0m.to_csv(RESULTS / "followup-mature-p1m-p0m-genes.csv", index=False)
    p1m.to_csv(RESULTS / "followup-mature-p1m-genes.csv", index=False)

    top_unit = ranking.iloc[0]["unit"]
    top_result = de_results[top_unit]
    positive_rows = []
    for gene in POSITIVE_CONTROLS:
        if gene in top_result.index:
            gene_row = top_result.loc[gene]
            lfc = float(gene_row["log2FC"])
            positive_rows.append(
                {
                    "top_unit": top_unit,
                    "gene": gene,
                    "detected": True,
                    "log2FC": lfc,
                    "pvalue": gene_row["pvalue"],
                    "FDR": gene_row["FDR"],
                    "FDR_uniform": gene_row["FDR_uniform"],
                    "direction": "UP with age" if lfc > 0 else "DOWN with age" if lfc < 0 else "NO CHANGE",
                }
            )
        else:
            positive_rows.append(
                {
                    "top_unit": top_unit,
                    "gene": gene,
                    "detected": False,
                    "log2FC": np.nan,
                    "pvalue": np.nan,
                    "FDR": np.nan,
                    "FDR_uniform": np.nan,
                    "direction": "UNTESTABLE: gene absent after total-count filtering",
                }
            )
    pd.DataFrame(positive_rows).to_csv(
        RESULTS / "followup-mature-p1m-positive-controls.csv", index=False
    )

    orientation = pd.DataFrame(orientation_rows)
    orientation.to_csv(RESULTS / "followup-mature-p1m-orientation-check.csv", index=False)
    orientation_pass = orientation[orientation["sign_matches_normalized_means"]]
    if orientation_pass.empty:
        raise AssertionError("Failed log2FC orientation assertion")
    orientation_example = orientation_pass.iloc[0]
    assert CONTRAST == ["age_group", "30m", "18m"]
    print(
        f"Orientation passed: {orientation_example['unit']} {orientation_example['gene']} "
        f"log2FC={orientation_example['log2FC']:.4f}; positive means 30m versus 18m up",
        flush=True,
    )

    original = pd.read_csv(RESULTS / "m5-program-hierarchy.csv")
    original_p1 = original[original["tier"] == "P1_Clean"].set_index("gene")
    new_p1 = p1m.set_index("mouse_gene_symbol")
    original_genes = set(original_p1.index)
    new_genes = set(new_p1.index)
    common = original_genes & new_genes
    sign_flips = sorted(
        gene for gene in common
        if original_p1.loc[gene, "direction"] != new_p1.loc[gene, "direction"]
    )
    overlap_rows = []
    for gene in sorted(original_genes | new_genes):
        old_direction = original_p1.loc[gene, "direction"] if gene in original_genes else pd.NA
        new_direction = new_p1.loc[gene, "direction"] if gene in new_genes else pd.NA
        if gene not in original_genes:
            status = "gained_in_P1m"
        elif gene not in new_genes:
            status = "lost_from_P1"
        elif gene in sign_flips:
            status = "common_sign_flip"
        else:
            status = "common_same_direction"
        overlap_rows.append(
            {
                "mouse_gene_symbol": gene,
                "original_P1_direction": old_direction,
                "P1m_direction": new_direction,
                "status": status,
            }
        )
    pd.DataFrame(overlap_rows).to_csv(
        RESULTS / "followup-mature-p1m-p1-overlap.csv", index=False
    )

    used_mice = pd.DataFrame(used_mouse_rows).drop_duplicates("mouse.id")
    source_unique = source_mice.drop_duplicates("mouse.id")
    p0_counts = direction_counts(p0m)
    p1_counts = direction_counts(p1m)
    p0_human_counts = direction_counts(p0m.dropna(subset=["human_ortholog"]))
    p1_human_counts = direction_counts(p1m.dropna(subset=["human_ortholog"]))
    max_bh_discrepancy = float(bh_check["bh_max_abs_discrepancy"].max())
    bh_pass = bool(bh_check["bh_pass"].all())
    end_time = utc_now()
    duration = time.time() - start_clock
    summary = {
        "model": MODEL,
        "derivation_start_utc": start_time,
        "derivation_end_utc": end_time,
        "derivation_duration_s": duration,
        "reference_age": REF_LEVEL,
        "target_age": TARGET_LEVEL,
        "excluded_ages": ";".join(EXCLUDED_AGES),
        "n_input_files": n_input_files,
        "n_units_evaluated": len(eligibility),
        "n_eligible_units": len(eligible),
        "n_ineligible_units": int((eligibility["eligible"] == "no").sum()),
        "n_p1m_floor10_units": int(eligibility["p1m_floor10_eligible"].sum()),
        "p0m_total": p0_counts["total"],
        "p0m_up": p0_counts["up"],
        "p0m_down": p0_counts["down"],
        "p0m_human_mapped_total": p0_human_counts["total"],
        "p0m_human_mapped_up": p0_human_counts["up"],
        "p0m_human_mapped_down": p0_human_counts["down"],
        "p1m_total": p1_counts["total"],
        "p1m_up": p1_counts["up"],
        "p1m_down": p1_counts["down"],
        "p1m_human_mapped_total": p1_human_counts["total"],
        "p1m_human_mapped_up": p1_human_counts["up"],
        "p1m_human_mapped_down": p1_human_counts["down"],
        "source_18m_female_mice": int(((source_unique["age"] == "18m") & (source_unique["sex"] == "female")).sum()),
        "source_18m_male_mice": int(((source_unique["age"] == "18m") & (source_unique["sex"] == "male")).sum()),
        "source_30m_female_mice": int(((source_unique["age"] == "30m") & (source_unique["sex"] == "female")).sum()),
        "source_30m_male_mice": int(((source_unique["age"] == "30m") & (source_unique["sex"] == "male")).sum()),
        "used_18m_female_mice": int(((used_mice["age"] == "18m") & (used_mice["sex"] == "female")).sum()),
        "used_18m_male_mice": int(((used_mice["age"] == "18m") & (used_mice["sex"] == "male")).sum()),
        "used_30m_female_mice": int(((used_mice["age"] == "30m") & (used_mice["sex"] == "female")).sum()),
        "used_30m_male_mice": int(((used_mice["age"] == "30m") & (used_mice["sex"] == "male")).sum()),
        "n_units_sex_kept": int((ranking["sex_term_status"] == "kept").sum()),
        "n_units_sex_dropped": int((ranking["sex_term_status"] == "dropped").sum()),
        "bh_recomputation_pass": bh_pass,
        "bh_max_abs_discrepancy": max_bh_discrepancy,
        "p1_common": len(common),
        "p1_lost": len(original_genes - new_genes),
        "p1m_gained": len(new_genes - original_genes),
        "p1_sign_flips": len(sign_flips),
        "top_unit": top_unit,
        "orientation_assertion_pass": True,
        "orientation_example_unit": orientation_example["unit"],
        "orientation_example_gene": orientation_example["gene"],
        "orientation_example_log2FC": orientation_example["log2FC"],
        "all_30m_male": True,
    }
    pd.DataFrame([summary]).to_csv(
        RESULTS / "followup-mature-p1m-summary.csv", index=False
    )
    if not bh_pass:
        raise AssertionError(f"BH recomputation mismatch: max difference {max_bh_discrepancy}")
    print(f"Mature P1m/P0m derivation finished {end_time}", flush=True)
    print(pd.Series(summary).to_string(), flush=True)


if __name__ == "__main__":
    main()
