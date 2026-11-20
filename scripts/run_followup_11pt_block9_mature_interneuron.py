#!/usr/bin/env python3
"""Block 9 mature-only MSSM and HBCC interneuron refit."""

from __future__ import annotations

import argparse
import gc
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import run_rev2_de_interneuron as original
from src.followup_organs import write_clean_csv, write_clean_markdown

RESULTS = REPO / "results"
REPORT = REPO / "reports" / "followup-11pt-block9-mature-interneuron.md"
GENE_OUTPUT = RESULTS / "followup-11pt-block9-mature-interneuron.csv"
REPLICATION_OUTPUT = RESULTS / "followup-11pt-block9-mature-replication.csv"
MSSM_PATH = Path("/home/reza/workspace/rev2-mssm-null-scratch/mssm.h5ad")
HBCC_PATH = REPO / "scratch" / "hbcc.h5ad"
DEFAULT_SCRATCH = Path(
    "/home/reza/projects/agescape-revision-data/tmp_organs/followup-11pt-block9-mature-interneuron"
)
MATURE_AGE = 30.0
MIN_AVAILABLE_GIB = 20.0
SHARED_GENES = ["FBN1", "GABRB1", "IGF1R", "PHF3", "TRMT11", "ZBTB16"]
SUBTYPES = [
    "VIP GABAergic cortical interneuron",
    "pvalb GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron",
]
SEED = original.SEED


def available_memory_gib() -> float:
    human = subprocess.run(["free", "-h"], check=True, capture_output=True, text=True)
    print(human.stdout, flush=True)
    raw = subprocess.run(["free", "-b"], check=True, capture_output=True, text=True)
    lines = raw.stdout.splitlines()
    header = lines[0].split()
    mem = next(line.split() for line in lines if line.startswith("Mem:"))
    available = int(mem[header.index("available") + 1])
    return available / (1024 ** 3)


def require_memory() -> float:
    available = available_memory_gib()
    if available < MIN_AVAILABLE_GIB:
        raise MemoryError(
            f"Only {available:.1f} GiB available. Wait and retry when at least "
            f"{MIN_AVAILABLE_GIB:.0f} GiB is available."
        )
    return available


def mssm_source_cohort_summary() -> dict:
    import anndata as ad

    require_memory()
    backed = ad.read_h5ad(MSSM_PATH, backed="r")
    obs = backed.obs
    age = original._age_from_stage(obs["development_stage"])
    common = original._cat_eq(obs["disease"], "normal")
    common &= original._cat_eq(obs["is_primary_data"], "true")
    if "raw_sum" in obs.columns:
        common &= pd.to_numeric(obs["raw_sum"], errors="coerce").ge(1500)
    else:
        common &= pd.to_numeric(obs["n_counts"], errors="coerce").ge(1500)
    donors = pd.DataFrame({
        "donor": obs["donor_id"].astype(str).to_numpy(),
        "age": age.to_numpy(float),
    })
    donors = donors.loc[common.to_numpy() & donors["age"].between(20, 120)]
    backed.file.close()
    donors = donors.drop_duplicates("donor")
    mature = donors.loc[donors["age"].ge(MATURE_AGE)]
    return {
        "full_n_donors": int(len(donors)),
        "n_donors": int(len(mature)),
        "age_min": float(mature["age"].min()),
        "age_median": float(mature["age"].median()),
        "age_max": float(mature["age"].max()),
    }


def fetch_mature_unit_h5ad(h5ad_path: Path, cell_type: str, scratch: Path, tag: str):
    """Original fetch_unit_h5ad with the sole selection change age 30 to 120."""
    import anndata as ad

    safe = f"{tag}_{cell_type}".replace(" ", "_").replace("/", "_")
    cache = scratch / f"counts_{safe}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        meta = pd.DataFrame({
            k: z[k]
            for k in ("donor", "age", "sex_male", "n_cells", "mean_log_umi", "library")
        })
        return z["counts"], meta, np.array([str(g) for g in z["genes"]]), cache
    b = ad.read_h5ad(h5ad_path, backed="r")
    obs = b.obs
    cols = set(obs.columns)
    mask = pd.Series(True, index=obs.index)
    mask &= original._cat_eq(obs["cell_type"], cell_type)
    if "disease" in cols:
        mask &= original._cat_eq(obs["disease"], "normal")
    if "is_primary_data" in cols:
        mask &= original._cat_eq(obs["is_primary_data"], "true")
    if "raw_sum" in cols:
        mask &= pd.to_numeric(obs["raw_sum"], errors="coerce") >= 1500
    elif "n_counts" in cols:
        mask &= pd.to_numeric(obs["n_counts"], errors="coerce") >= 1500
    age = original._age_from_stage(obs["development_stage"])
    mask &= age.between(MATURE_AGE, 120).to_numpy()
    sub_obs = obs[mask]
    per_donor = sub_obs.groupby("donor_id", observed=True).size()
    donors = sorted(
        per_donor[per_donor >= original.MIN_CELLS_PER_DONOR_CT].index.astype(str)
    )
    if len(donors) < original.MIN_DONORS:
        original.log(
            f"  [{tag} | {cell_type}] too few donors "
            f"({len(donors)} < {original.MIN_DONORS}); skipped"
        )
        b.file.close()
        return None, None, None, None
    keep = mask & obs["donor_id"].isin(donors)
    n_per_donor = {str(k): int(v) for k, v in per_donor.items()}
    batches, current, current_n = [], [], 0
    for donor in donors:
        n_donor = n_per_donor[donor]
        if current and current_n + n_donor > 15000:
            batches.append(current)
            current, current_n = [], 0
        current.append(donor)
        current_n += n_donor
    if current:
        batches.append(current)
    batch_masks = [
        (keep & obs["donor_id"].isin(batch)).to_numpy() for batch in batches
    ]
    del obs, mask, sub_obs, per_donor, age, keep
    gc.collect()
    counts = None
    genes = None
    meta_rows = []
    for batch_index, batch in enumerate(batches):
        require_memory()
        sub = b[batch_masks[batch_index]].to_memory()
        if "raw_sum" not in cols and "n_counts" not in cols:
            matrix = sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
            raw_sum = np.asarray(matrix.sum(axis=1)).ravel()
            sub = sub[raw_sum >= 1500].copy()
            original.log(
                f"  [{tag} | {cell_type}] raw_sum computed from X; "
                f"{int((raw_sum >= 1500).sum())}/{len(raw_sum)} cells kept"
            )
        batch_counts, batch_meta, batch_genes = original.pseudobulk_counts(sub)
        del sub
        gc.collect()
        if genes is None:
            genes = batch_genes
            counts = np.zeros((len(donors), len(genes)), dtype=np.int64)
        else:
            assert list(batch_genes) == list(genes), f"var changed between batches for {cell_type}"
        positions = {donor: index for index, donor in enumerate(donors)}
        for row_index, donor in enumerate(batch_meta["donor"]):
            counts[positions[donor]] = batch_counts[row_index]
        meta_rows.append(batch_meta)
        original.log(
            f"  [{tag} | {cell_type}] batch {batch_index + 1}/{len(batches)}: "
            f"donors {sum(len(x) for x in batches[:batch_index + 1])}/{len(donors)}"
        )
    b.file.close()
    del b
    meta = pd.concat(meta_rows, ignore_index=True)
    meta = meta.set_index("donor").loc[donors].reset_index()
    original.log(
        f"  [{tag} | {cell_type}] {len(meta)} donors, {len(genes)} genes, "
        f"{int(counts.sum())} total counts"
    )
    if len(meta) < original.MIN_DONORS:
        original.log(f"  [{tag} | {cell_type}] too few donors after pseudobulk; skipped")
        return None, None, None, None
    ok = original.save_npz_retry(
        cache,
        counts=counts,
        genes=np.array(genes),
        donor=meta["donor"].to_numpy(),
        age=meta["age"].to_numpy(float),
        sex_male=meta["sex_male"].to_numpy(float),
        n_cells=meta["n_cells"].to_numpy(np.int64),
        mean_log_umi=meta["mean_log_umi"].to_numpy(float),
        library=meta["library"].to_numpy(float),
    )
    if not ok:
        original.log("  WARNING: counts cache save failed; continuing in-memory")
    return counts, meta, genes, cache


def model_frames(merged: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    wls = merged[[
        "gene", "slope", "se", "pvalue", "q_bh", "ci95_low", "ci95_high",
        "mean_cp10k", "detected_frac", "n_donors",
    ]].copy()
    wls["model"] = "rev2_voom_wls"
    wls["primary"] = True
    ols = pd.DataFrame({
        "gene": merged["gene"],
        "slope": merged["slope_ols"],
        "se": merged["se_ols"],
        "pvalue": merged["pvalue_ols"],
        "q_bh": merged["q_bh_ols"],
        "ci95_low": merged["ci95_low_ols"],
        "ci95_high": merged["ci95_high_ols"],
        "mean_cp10k": merged["mean_cp10k"],
        "detected_frac": merged["detected_frac"],
        "n_donors": merged["n_donors"],
        "model": "age_plus_sex_depth_ols",
        "primary": False,
    })
    unadjusted = sensitivity.rename(columns={
        "slope_ols": "slope",
        "se_ols": "se",
        "pvalue_ols": "pvalue",
        "q_bh_ols": "q_bh",
        "ci95_low_ols": "ci95_low",
        "ci95_high_ols": "ci95_high",
    }).copy()
    unadjusted["model"] = "age_plus_sex_ols_depth_unadjusted"
    unadjusted["primary"] = False
    unadjusted["n_donors"] = int(merged["n_donors"].iloc[0])
    unadjusted = unadjusted.merge(
        merged[["gene", "mean_cp10k", "detected_frac"]], on="gene", how="left",
        validate="one_to_one",
    )
    return pd.concat([wls, ols, unadjusted], ignore_index=True)


def run_mssm(scratch: Path) -> tuple[pd.DataFrame, list[dict]]:
    if not MSSM_PATH.exists():
        raise FileNotFoundError(f"Required MSSM file is missing: {MSSM_PATH}")
    gene_frames = []
    summaries = []
    all_donors = []
    for cell_type in original.INTERNEURON_TYPES:
        available = require_memory()
        original.log(f"[MSSM_DLPFC | {cell_type}] available memory {available:.1f} GiB")
        counts, meta, genes, _ = fetch_mature_unit_h5ad(
            MSSM_PATH, cell_type, scratch, "mssm_mature30"
        )
        if counts is None:
            summaries.append({
                "cohort": "MSSM_DLPFC", "cell_type": cell_type,
                "short_type": original.SHORT[cell_type], "status": "UNTESTABLE",
                "status_reason": "fewer_than_10_mature_donors", "n_donors": 0,
            })
            continue
        if len(meta) < original.MIN_DONORS:
            raise AssertionError(f"Donor floor violated for {cell_type}")
        if meta["donor"].duplicated().any() or not meta["age"].ge(MATURE_AGE).all():
            raise AssertionError(f"Mature donor validation failed for {cell_type}")
        if counts.shape != (len(meta), len(genes)) or not meta["n_cells"].ge(20).all():
            raise AssertionError(f"Pseudobulk validation failed for {cell_type}")
        original.log(f"[MSSM_DLPFC | {cell_type}] fitting original voom-WLS + OLS")
        summary, merged = original.run_unit_de(
            counts, meta, genes, "MSSM_DLPFC", cell_type, scratch, "mssm_mature30"
        )
        sensitivity = original.ols_arm(
            counts, meta, genes, columns=("age", "sex_male")
        )
        frame = model_frames(merged, sensitivity)
        frame["cohort"] = "MSSM_DLPFC"
        frame["cell_type"] = cell_type
        frame["short_type"] = original.SHORT[cell_type]
        frame["status"] = "FIT"
        frame["age_min"] = float(meta["age"].min())
        frame["age_median"] = float(meta["age"].median())
        frame["age_max"] = float(meta["age"].max())
        frame["n_genes_tested"] = int(summary["n_genes_tested"])
        frame["design"] = summary["design"]
        gene_frames.append(frame)
        summary.update({
            "short_type": original.SHORT[cell_type],
            "status": "FIT",
            "status_reason": "none",
            "age_median": float(meta["age"].median()),
            "n_fdr_depth_unadjusted_ols": int((sensitivity["q_bh_ols"] < 0.05).sum()),
            "min_available_gib_before_unit": available,
        })
        summaries.append(summary)
        all_donors.append(meta[["donor", "age"]])
        original.log(
            f"  -> mature WLS FDR {summary['n_fdr_wls']}; "
            f"depth-adjusted OLS FDR {summary['n_fdr_ols_same_genes']}; "
            f"depth-unadjusted OLS FDR {summary['n_fdr_depth_unadjusted_ols']}"
        )
        del counts, meta, genes, merged, sensitivity, frame
        gc.collect()
    columns = [
        "cohort", "cell_type", "short_type", "model", "primary", "status", "gene",
        "slope", "se", "ci95_low", "ci95_high", "pvalue", "q_bh", "n_donors",
        "n_genes_tested", "mean_cp10k", "detected_frac", "age_min", "age_median",
        "age_max", "design",
    ]
    genes_out = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame(columns=columns)
    genes_out = genes_out.reindex(columns=columns)
    write_clean_csv(genes_out, GENE_OUTPUT)
    if all_donors:
        cohort_donors = pd.concat(all_donors).drop_duplicates("donor")
        cohort_summary = {
            "n_donors": int(len(cohort_donors)),
            "age_min": float(cohort_donors["age"].min()),
            "age_median": float(cohort_donors["age"].median()),
            "age_max": float(cohort_donors["age"].max()),
        }
    else:
        cohort_summary = {"n_donors": 0}
    payload = {"units": summaries, "mssm_cohort": cohort_summary, "seed": SEED}
    (scratch / "run_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return genes_out, summaries


def validate_shared_genes() -> None:
    moderated = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes-moderated.csv")
    moderated_sets = {
        subtype: set(moderated.loc[
            moderated["cell_type"].eq(subtype) & moderated["tier"].eq("both_FDR"), "gene"
        ])
        for subtype in SUBTYPES
    }
    intersection = set.intersection(*moderated_sets.values())
    if intersection != set(SHARED_GENES):
        raise AssertionError(f"Moderated full-range intersection is not the frozen six: {intersection}")
    unmoderated = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv")
    unmoderated_sets = {
        subtype: set(unmoderated.loc[
            unmoderated["cell_type"].eq(subtype) & unmoderated["tier"].eq("both_FDR"), "gene"
        ])
        for subtype in SUBTYPES
    }
    if set.intersection(*unmoderated_sets.values()) != set(SHARED_GENES + ["DYNC1I1"]):
        raise AssertionError("The unmoderated seven-gene audit does not match the frozen target")


def write_replication() -> pd.DataFrame:
    source = pd.read_csv(RESULTS / "reviewer-phase1-interneuron-replication-statistics.csv")
    rows = []
    for _, reference in source.iterrows():
        full = reference.to_dict()
        full.update({
            "full_range": True,
            "status": "FULL_RANGE_REFERENCE",
            "status_reason": "read_from_reviewer_phase1_target",
        })
        rows.append(full)
        mature = {column: np.nan for column in source.columns}
        mature["cell_type"] = reference["cell_type"]
        mature["subset"] = reference["subset"]
        mature.update({
            "full_range": False,
            "status": "UNTESTABLE",
            "status_reason": f"HBCC_h5ad_absent:{HBCC_PATH}",
        })
        rows.append(mature)
    out = pd.DataFrame(rows, columns=list(source.columns) + [
        "full_range", "status", "status_reason",
    ])
    write_clean_csv(out, REPLICATION_OUTPUT)
    return out


def fmt(value, digits=4) -> str:
    if value is None or pd.isna(value):
        return "UNTESTABLE"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.{digits}f}"


def build_report(genes: pd.DataFrame, summaries: list[dict], scratch: Path) -> None:
    validate_shared_genes()
    target = pd.read_csv(RESULTS / "reviewer-phase1-interneuron-replication-statistics.csv")
    full_de = pd.read_csv(RESULTS / "reviewer-rev2-human-de-summary.csv")
    full_rev2 = pd.read_csv(RESULTS / "reviewer-rev2-replication-summary.csv")
    run_summary = json.loads((scratch / "run_summary.json").read_text(encoding="utf-8"))
    if "mssm_source_cohort" not in run_summary:
        run_summary["mssm_source_cohort"] = mssm_source_cohort_summary()
        (scratch / "run_summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
    source_cohort = run_summary["mssm_source_cohort"]
    summary_by_type = {row["cell_type"]: row for row in summaries}
    primary = genes[genes["model"].eq("rev2_voom_wls")]
    ols = genes[genes["model"].eq("age_plus_sex_depth_ols")]
    unadjusted = genes[genes["model"].eq("age_plus_sex_ols_depth_unadjusted")]

    donor_lines = []
    fdr_lines = []
    replication_lines = []
    power_lines = []
    for cell_type in original.INTERNEURON_TYPES:
        short = original.SHORT[cell_type]
        summary = summary_by_type.get(cell_type, {"status": "UNTESTABLE", "n_donors": 0})
        full_row = full_de[
            full_de["cohort"].eq("MSSM_DLPFC") & full_de["cell_type"].eq(cell_type)
        ].iloc[0]
        full_hbcc = full_de[
            full_de["cohort"].eq("HBCC_DLPFC") & full_de["cell_type"].eq(cell_type)
        ].iloc[0]
        donor_lines.append(
            f"| {short} | {int(full_row['n_donors'])} | {int(summary['n_donors'])} | "
            f"{fmt(summary.get('age_min'), 0)} to {fmt(summary.get('age_max'), 0)} | "
            f"{fmt(summary.get('age_median'), 1)} | {int(full_hbcc['n_donors'])} | UNTESTABLE |"
        )
        selected = target[
            target["cell_type"].eq(cell_type) & target["subset"].eq("mssm_fdr_genes")
        ].iloc[0]
        all_shared = target[
            target["cell_type"].eq(cell_type) & target["subset"].eq("all_shared_genes")
        ].iloc[0]
        mature_primary = primary[primary["cell_type"].eq(cell_type)]
        mature_ols = ols[ols["cell_type"].eq(cell_type)]
        fdr_lines.append(
            f"| {short} | {int(all_shared['n_mssm_fdr'])} | {len(mature_primary)} | "
            f"{int(mature_primary['q_bh'].lt(0.05).sum()) if len(mature_primary) else 'UNTESTABLE'} | "
            f"{int(mature_ols['q_bh'].lt(0.05).sum()) if len(mature_ols) else 'UNTESTABLE'} | "
            f"{int(all_shared['n_hbcc_fdr'])} | UNTESTABLE |"
        )
        replication_lines.append(
            f"| {short} | {fmt(selected['pearson_effect_r'])} | UNTESTABLE | "
            f"{int(selected['same_direction'])} | UNTESTABLE | "
            f"{int(all_shared['both_fdr_same_direction'])} | "
            f"{fmt(all_shared['expected_both_fdr_under_independence'], 3)} | UNTESTABLE |"
        )
        power_lines.append(
            f"- {short}: MSSM retained {int(summary['n_donors'])} donors, above the original "
            f"10-donor floor. HBCC has no mature refit because the input file is absent. "
            f"Cross-cohort correlation, overlap power, and MDE80 are therefore not estimable. "
            f"This is UNTESTABLE, not a null."
        )

    shared_lines = []
    for gene in SHARED_GENES:
        cells = [gene]
        for cell_type in SUBTYPES:
            part = primary[
                primary["cell_type"].eq(cell_type) & primary["gene"].eq(gene)
            ]
            if len(part) == 0:
                status = "NOT_TESTED"
            else:
                status = "PASS" if float(part["q_bh"].iloc[0]) < 0.05 else "NOT_FDR"
            cells.extend([status, "UNTESTABLE"])
        shared_lines.append("| " + " | ".join(cells) + " |")

    full_rev2_lines = []
    for cell_type in SUBTYPES:
        row = full_rev2[full_rev2["cell_type"].eq(cell_type)].iloc[0]
        full_rev2_lines.append(
            f"| {original.SHORT[cell_type]} | {fmt(row['gw_pearson_r'])} | "
            f"{int(row['both_fdr_same_direction'])} | {fmt(row['hypergeometric_expected'], 3)} | "
            f"{fmt(row['block_perm_p'])} |"
        )

    cohort = run_summary["mssm_cohort"]
    sensitivity_lines = []
    for cell_type in original.INTERNEURON_TYPES:
        part = unadjusted[unadjusted["cell_type"].eq(cell_type)]
        sensitivity_lines.append(
            f"| {original.SHORT[cell_type]} | "
            f"{int(part['q_bh'].lt(0.05).sum()) if len(part) else 'UNTESTABLE'} |"
        )

    text = f"""# Block 9 mature-only interneuron refit

## Verdict

INCONCLUSIVE. MSSM was refit for donors age 30 years and older, but every subtype is UNTESTABLE for mature-only MSSM-to-HBCC replication because `/home/reza/projects/agescape-revision/scratch/hbcc.h5ad` is absent. Under the frozen rule, any subtype UNTESTABLE in either cohort makes the block INCONCLUSIVE. No cross-cohort null is claimed. No gene among FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16 can be confirmed as a surviving mature-only cross-cohort shared gene. All six are UNTESTABLE in HBCC for VIP, pvalb, and sst. In the supplemental MSSM-only WLS fit, GABRB1, IGF1R, and TRMT11 pass FDR in all three subtypes; FBN1 and PHF3 pass in pvalb and sst but not VIP; ZBTB16 passes only in sst. These MSSM-only patterns are not substitutes for replication.

## Frozen specification

The only change to the original donor selection is restriction to donors age 30 years and older in both cohorts. The source of this cutoff is `scripts/run_mature_only_composition.py`: Human mature = 30 years and older, Reza's maturity definition, based on brain development completing around age 30. Mouse maturity is not used.

The original `run_unit_de` implementation is reused for donor-pseudobulk voom-WLS and OLS. QC remains primary data, disease normal, at least 1500 raw UMIs per cell, at least 20 eligible cells per donor and cell type, and at least 10 donors per unit. Pseudobulk is summed raw counts. The response is log1p counts per 10,000. The design is age in years plus sex plus donor mean log UMI. Low-expression filtering, technical-gene masking, one-iteration voom precision weighting, and donor unit inference are unchanged. There is no downsampling.

The technical mask excludes uppercase MT-, RPL, RPS, MRPL, and MRPS prefixes, plus XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D. The original human model has no study, pool, or cohort indicator within either separately fit cohort. The shared-rule unadjusted sensitivity therefore removes depth and fits the original OLS function with age plus sex. It is labeled depth-unadjusted and does not replace the depth-adjusted primary.

The required replication statistics are MSSM-to-HBCC effect-size correlations, same-direction counts among MSSM FDR genes, both-FDR same-direction counts against hypergeometric expectation, and the moderated three-subtype shared-gene intersection. These require both cohorts and are UNTESTABLE in the mature-only arm. The original full-range rows are copied from the specified targets, not recomputed.

The frozen verdict rule is: PROVES only if at least one subtype retains significant MSSM-to-HBCC replication, with the effect-size correlation passing the original replication criterion and both-FDR same-direction count above hypergeometric expectation. DOES NOT PROVE requires no replication with adequate donor counts and per-subtype power reasoning. INCONCLUSIVE applies if donor counts are too small or any subtype is UNTESTABLE in either cohort. The last condition applies here.

## Families fixed before correction

The original report states: Benjamini-Hochberg correction was applied within each cohort and cell type. The original implementation further defines this as BH q within each cohort by cell type over the filtered tested genes. The mature WLS and depth-adjusted OLS families retain this definition separately for each model. The depth-unadjusted OLS sensitivity is BH corrected separately within each cohort by subtype and is never pooled with the primary family. No cross-cohort mature family was created after seeing results.

The fixed seed inherited from the original analysis is {SEED}. The MSSM fits are deterministic. No replication permutation was run because HBCC is absent.

## Inputs and blockers

MSSM input: `/home/reza/workspace/rev2-mssm-null-scratch/mssm.h5ad`, present. Its four mature pseudobulk caches were built without copying the h5ad. Each backed h5ad batch was gated on at least 20 GiB available memory.

HBCC input: `/home/reza/projects/agescape-revision/scratch/hbcc.h5ad`, absent. The HBCC arm stopped at this path. No download or source substitution was attempted.

The MSSM source-cohort QC route retains {int(source_cohort['full_n_donors'])} donors at age 20 and older and {int(source_cohort['n_donors'])} donors at age 30 and older before cell-type eligibility. The mature source-cohort age range is {fmt(source_cohort.get('age_min'), 0)} to {fmt(source_cohort.get('age_max'), 0)}, median {fmt(source_cohort.get('age_median'), 1)}. The original comparison target describes 212 MSSM and 88 HBCC adult donors. Across the union of tested MSSM interneuron units, {int(cohort['n_donors'])} mature donors met the 20-cell floor in at least one unit, age {fmt(cohort.get('age_min'), 0)} to {fmt(cohort.get('age_max'), 0)}, median {fmt(cohort.get('age_median'), 1)}. The table below uses original Rev2 per-unit eligible counts for like-for-like donor comparison.

## Donor retention

| Subtype | MSSM full N | MSSM mature N | Mature age range | Mature median | HBCC full N | HBCC mature N |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(donor_lines)}

All MSSM units remain above the original 10-donor floor. No MSSM subtype is UNTESTABLE on donor count.

## Full-range and mature-only FDR counts

Full-range values are the specified original Phase 1 OLS comparison targets. Mature values are the required Rev2 WLS primary and its original depth-adjusted OLS arm. This distinction is retained rather than presenting different estimators as identical. The mature tested-gene column is the BH family size for each MSSM subtype and model.

| Subtype | MSSM full OLS FDR | MSSM mature tested genes | MSSM mature WLS FDR | MSSM mature OLS FDR | HBCC full OLS FDR | HBCC mature FDR |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(fdr_lines)}

Depth-unadjusted OLS sensitivity FDR counts, each using its own cohort-by-subtype family:

| Subtype | MSSM mature age plus sex OLS FDR |
|---|---:|
{chr(10).join(sensitivity_lines)}

## Replication comparison

Full-range Pearson r and same-direction counts below are among original MSSM FDR genes. Both-FDR counts and expectations use all original shared tested genes. Every mature entry is UNTESTABLE because HBCC was not available.

| Subtype | Full r | Mature r | Full same direction | Mature same direction | Full both FDR | Full expected | Mature both FDR |
|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(replication_lines)}

For additional context, the original full-range Rev2 genome-wide statistics were:

| Subtype | Genome-wide r | Unmoderated both FDR | Hypergeometric expected | Block permutation p |
|---|---:|---:|---:|---:|
{chr(10).join(full_rev2_lines)}

The original full-range moderated both-FDR counts were 77 VIP, 65 pvalb, and 239 sst. The original Phase 1 targets were 24 VIP, 20 pvalb, and 62 sst. These are different original analysis versions and are not merged.

## Six moderated shared genes

The full-range moderated intersection is exactly 6 genes: FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16. DYNC1I1 is not in this intersection. It appears only when the unmoderated file is used. The mature moderated intersection itself is UNTESTABLE. The table gives supplemental MSSM mature WLS FDR status beside the required HBCC status. NOT_FDR means tested but q >= 0.05. NOT_TESTED means absent after the unchanged expression filter.

| Gene | MSSM VIP | HBCC VIP | MSSM pvalb | HBCC pvalb | MSSM sst | HBCC sst |
|---|---|---|---|---|---|---|
{chr(10).join(shared_lines)}

## Power and interpretation

{chr(10).join(power_lines)}

The high MSSM retention does not rescue cross-cohort replication. With no HBCC mature model there are zero eligible HBCC donor inputs for every subtype, no shared mature effect vectors, and no meaningful correlation or overlap MDE80. Accordingly, these results neither prove survival from maturity nor prove loss of the program. The MSSM-only FDR counts describe one cohort and cannot satisfy the two-cohort verdict rule.

## Outputs and reproduction

- `results/followup-11pt-block9-mature-interneuron.csv`: mature per-gene MSSM WLS, depth-adjusted OLS, and depth-unadjusted OLS results.
- `results/followup-11pt-block9-mature-replication.csv`: original columns plus `full_range`, `status`, and `status_reason`. Full-range reference rows are paired with UNTESTABLE mature rows.
- `scripts/run_followup_11pt_block9_mature_interneuron.py`: original model reuse, age 30 donor filter, memory guard, validation, and report generation.

Run with `/home/reza/miniconda3/envs/agescape/bin/python scripts/run_followup_11pt_block9_mature_interneuron.py`. Resume uses the existing pseudobulk caches. No manuscript or `report/` file was changed.
"""
    write_clean_markdown(text, REPORT)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    if not HBCC_PATH.exists():
        original.log(f"[HBCC_DLPFC] required h5ad absent: {HBCC_PATH}; arm stopped")
    else:
        raise RuntimeError(
            f"Unexpected HBCC file at {HBCC_PATH}. This frozen run was written for the verified missing-file arm."
        )
    validate_shared_genes()
    if args.report_only:
        genes = pd.read_csv(GENE_OUTPUT)
        summaries = json.loads((args.scratch_dir / "run_summary.json").read_text(encoding="utf-8"))["units"]
    else:
        genes, summaries = run_mssm(args.scratch_dir)
    write_replication()
    build_report(genes, summaries, args.scratch_dir)
    original.log(f"Wrote {GENE_OUTPUT}")
    original.log(f"Wrote {REPLICATION_OUTPUT}")
    original.log(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
