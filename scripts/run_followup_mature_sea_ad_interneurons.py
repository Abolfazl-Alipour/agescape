#!/usr/bin/env python3
"""Item 7 mature-derived SEA-AD interneuron replication test."""

from __future__ import annotations

import argparse
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

from src.aging_programs import load_program_tiers
from src.followup_organs import write_clean_markdown
from scripts import run_followup_sea_ad_interneurons as sea
from scripts.run_followup_mature_interneuron_replication import wait_for_memory
from scripts.run_rev2_de_interneuron import (
    coexpression_blocks,
    filter_genes,
    lodo_sensitivity,
    ols_arm,
)

RESULTS = REPO / "results"
MANIFEST = REPO / "data" / "manifest"
REPORT = REPO / "reports" / "followup-mature-interneuron.md"
DISCOVERY_SOURCE = RESULTS / "followup-mature-interneuron-discovery-ols.csv"
TARGET_SOURCE = RESULTS / "followup-mature-interneuron-genes.csv"
SUBCLASS = sea.SUBCLASS
CELL_TYPE = sea.CELL_TYPE
PRIMARY_COLUMNS = sea.PRIMARY_COLUMNS
FAMILY_ORDER = (
    "mature_both_fdr",
    "mature_discovery_all",
    "mature_discovery_without_both_fdr",
    "mouse_p1_all_young_derived",
    "mouse_p1_mature_discovery_overlap_young_derived",
    "mouse_p1_outside_mature_discovery_young_derived",
)
DISCOVERY_MODEL_LABEL = "mature MSSM age_plus_sex_depth OLS (30-plus)"
OUTPUT_NAMES = (
    "followup-mature-interneuron-sea-ad-gene-effects.csv",
    "followup-mature-interneuron-sea-ad-genome-wide.csv",
    "followup-mature-interneuron-sea-ad-replication-summary.csv",
    "followup-mature-interneuron-sea-ad-sensitivities.csv",
    "followup-mature-interneuron-sea-ad-lodo-summary.csv",
    "followup-mature-interneuron-sea-ad-lodo-genes.csv",
)


def load_discovery(source: Path, short: str) -> pd.DataFrame:
    columns = [
        "cell_type", "model", "gene", "slope_log_cp10k_per_year_mssm",
        "pvalue_mssm", "q_bh_mssm", "slope_log_cp10k_per_year_hbcc",
        "pvalue_hbcc", "q_bh_hbcc",
    ]
    frame = pd.read_csv(source, usecols=columns)
    frame = frame[(frame["cell_type"] == CELL_TYPE[short]) & (frame["model"] == "age_plus_sex_depth")].copy()
    assert not frame["gene"].duplicated().any()
    assert len(frame) > 0
    return frame.drop(columns=["cell_type", "model"]).rename(columns={
        "slope_log_cp10k_per_year_mssm": "mssm_slope_original_ols",
        "pvalue_mssm": "mssm_p_original_ols",
        "q_bh_mssm": "mssm_q_original_ols",
        "slope_log_cp10k_per_year_hbcc": "hbcc_slope_original_ols",
        "pvalue_hbcc": "hbcc_p_original_ols",
        "q_bh_hbcc": "hbcc_q_original_ols",
    })


def load_targets(short: str):
    frame = pd.read_csv(TARGET_SOURCE)
    frame = frame[frame["cell_type"] == CELL_TYPE[short]].copy()
    assert len(frame) > 0
    assert not frame["gene"].duplicated().any()
    both = set(frame.loc[frame["tier"] == "both_FDR", "gene"])
    discovery = set(frame["gene"])
    assert len(both) > 0
    assert both <= discovery
    directions = dict(zip(frame["gene"], np.sign(frame["slope_mssm"])))
    p1 = load_program_tiers(REPO, target="human")["P1_Clean"]
    p1_direction = {**{gene: 1.0 for gene in p1["up"]}, **{gene: -1.0 for gene in p1["down"]}}
    assert len(p1_direction) == 118
    return frame, both, discovery, directions, p1_direction


def family_membership(genes: pd.Index, short: str):
    target, both, discovery, directions, p1_direction = load_targets(short)
    target = target.set_index("gene")
    membership = {
        "mature_both_fdr": genes.isin(both),
        "mature_discovery_all": genes.isin(discovery),
        "mature_discovery_without_both_fdr": genes.isin(discovery - both),
        "mouse_p1_all_young_derived": genes.isin(p1_direction),
        "mouse_p1_mature_discovery_overlap_young_derived": genes.isin(set(p1_direction) & discovery),
        "mouse_p1_outside_mature_discovery_young_derived": genes.isin(set(p1_direction) - discovery),
    }
    membership = {key: np.asarray(value, dtype=bool) for key, value in membership.items()}
    source_sign = {}
    source_size = {}
    for family in FAMILY_ORDER:
        if family.startswith("mouse_p1"):
            selected = set(p1_direction)
            if family == "mouse_p1_mature_discovery_overlap_young_derived":
                selected &= discovery
            elif family == "mouse_p1_outside_mature_discovery_young_derived":
                selected -= discovery
            source_sign[family] = np.array([p1_direction.get(gene, 0.0) for gene in genes])
            source_size[family] = len(selected)
        else:
            selected = discovery
            if family == "mature_both_fdr":
                selected = both
            elif family == "mature_discovery_without_both_fdr":
                selected = discovery - both
            source_sign[family] = np.array([directions.get(gene, 0.0) for gene in genes])
            source_size[family] = len(selected)
    return membership, source_sign, source_size, target


def validate_eligibility(data_root: Path) -> dict:
    eligibility = pd.read_csv(RESULTS / "followup-sea-ad-eligibility.csv")
    assert set(eligibility["subtype"]) == set(SUBCLASS)
    assert eligibility["metadata_eligible_donors"].eq(15).all()
    assert eligibility["eligible_donors_20_nuclei"].eq(15).all()
    assert eligibility["age_min"].eq(72).all() and eligibility["age_max"].eq(102).all()
    workbook = data_root / "metadata" / sea.WORKBOOK_NAME
    donors = sea.eligible_donors(workbook)
    expected_ids = set(donors["Donor ID"].astype(str))
    assert len(expected_ids) == 15
    cache_files = sorted((data_root / "interim" / "donor-pseudobulk").glob("*.npz"))
    assert len(cache_files) == 15
    assert {path.stem for path in cache_files} == expected_ids
    subtype_validation = {}
    for short in SUBCLASS:
        _, metadata, _ = sea.load_counts(data_root, short)
        assert len(metadata) == 15
        assert set(metadata["donor"]) == expected_ids
        assert metadata["age"].between(72, 102).all()
        assert metadata["n_cells"].ge(20).all()
        subtype_validation[short] = {
            "n_donors": len(metadata), "age_min": float(metadata["age"].min()),
            "age_max": float(metadata["age"].max()), "same_donor_ids": True,
            "minimum_nuclei": int(metadata["n_cells"].min()),
        }
    return {"eligibility_file": "results/followup-sea-ad-eligibility.csv", "n_cached_donors": 15, "subtypes": subtype_validation}


def analyze_subtype(data_root: Path, discovery_source: Path, short: str,
                    permutations: int, bootstraps: int):
    wait_for_memory()
    counts, metadata, genes = sea.load_counts(data_root, short)
    wls = sea.fit_wls(counts, metadata, genes)
    ols = ols_arm(counts, metadata, genes, columns=wls.attrs["active_columns"])
    frame = wls.merge(ols, on="gene", how="left")
    discovery = load_discovery(discovery_source, short)
    frame = frame.merge(discovery, on="gene", how="left")
    membership, source_sign, source_size, target = family_membership(frame["gene"], short)
    frame["in_mature_both_fdr"] = membership["mature_both_fdr"]
    frame["in_mature_discovery"] = membership["mature_discovery_all"]
    frame["in_mouse_p1_young_derived"] = membership["mouse_p1_all_young_derived"]
    for column in ("slope_mssm", "se_mssm", "ci95_low_mssm", "ci95_high_mssm",
                   "pvalue_mssm", "q_bh_mssm", "slope_hbcc", "q_bh_hbcc"):
        if column in target.columns:
            frame[f"mature_{column}"] = frame["gene"].map(target[column])
    frame["subtype"] = short
    frame["model"] = "primary_voom_wls"
    frame["n_donors"] = len(metadata)
    frame["resid_dof"] = wls.attrs["dof"]
    frame["design"] = " + ".join(["intercept"] + wls.attrs["active_columns"])
    frame["slope_per_10_years"] = 10 * frame["slope"]
    frame["ci95_low_per_10_years"] = 10 * frame["ci95_low"]
    frame["ci95_high_per_10_years"] = 10 * frame["ci95_high"]

    shared = frame["mssm_slope_original_ols"].notna().to_numpy()
    shared_frame = frame.loc[shared].reset_index(drop=True)
    shared_genes = shared_frame["gene"]
    shared_membership = {key: value[shared] for key, value in membership.items()}
    shared_sign = {key: value[shared] for key, value in source_sign.items()}
    n_shared, pearson, spearman = sea.correlation_rows(shared_frame, "mssm_slope_original_ols", "slope_ols")
    _, pearson_wls, spearman_wls = sea.correlation_rows(shared_frame, "mssm_slope_original_ols", "slope")
    shared_positions = np.flatnonzero(shared)
    keep, _ = filter_genes(counts, genes)
    blocks = coexpression_blocks(counts[:, keep][:, shared_positions], shared_genes.to_numpy())
    rng = np.random.default_rng(20261108)
    null_correlation = np.empty(permutations)
    null_same = {family: np.empty(permutations) for family in FAMILY_ORDER}
    null_fdr = {family: np.empty(permutations) for family in FAMILY_ORDER}
    seaad_sign = np.sign(shared_frame["slope"].to_numpy())
    seaad_fdr = shared_frame["q_bh"].to_numpy() < 0.05
    seaad_ols = shared_frame["slope_ols"].to_numpy()
    mssm_ols = shared_frame["mssm_slope_original_ols"].to_numpy()
    for permutation in range(permutations):
        perm_slope = sea.permute_within_blocks(seaad_ols, blocks, rng)
        perm_sign = sea.permute_within_blocks(seaad_sign, blocks, rng)
        perm_fdr = sea.permute_within_blocks(seaad_fdr, blocks, rng)
        null_correlation[permutation] = stats.pearsonr(mssm_ols, perm_slope).statistic
        for family in FAMILY_ORDER:
            selected = shared_membership[family]
            same = perm_sign[selected] * shared_sign[family][selected] > 0
            null_same[family][permutation] = same.sum()
            null_fdr[family][permutation] = (same & perm_fdr[selected]).sum()

    boot = np.random.default_rng(20261208)
    boot_correlation = np.empty(bootstraps)
    boot_spearman = np.empty(bootstraps)
    boot_correlation_wls = np.empty(bootstraps)
    boot_spearman_wls = np.empty(bootstraps)
    boot_same = {family: np.empty(bootstraps) for family in FAMILY_ORDER}
    observed_sign = np.sign(shared_frame["slope_ols"].to_numpy())
    seaad_wls = shared_frame["slope"].to_numpy()
    for index in range(bootstraps):
        idx = sea.bootstrap_indices(blocks, boot)
        boot_correlation[index] = stats.pearsonr(mssm_ols[idx], seaad_ols[idx]).statistic
        boot_spearman[index] = stats.spearmanr(mssm_ols[idx], seaad_ols[idx]).statistic
        boot_correlation_wls[index] = stats.pearsonr(mssm_ols[idx], seaad_wls[idx]).statistic
        boot_spearman_wls[index] = stats.spearmanr(mssm_ols[idx], seaad_wls[idx]).statistic
        for family in FAMILY_ORDER:
            selected = shared_membership[family]
            family_idx = np.flatnonzero(selected)
            target_blocks = blocks[family_idx]
            if len(family_idx):
                target_idx = sea.bootstrap_indices(target_blocks, boot)
                same = observed_sign[family_idx][target_idx] * shared_sign[family][family_idx][target_idx] > 0
                boot_same[family][index] = same.mean() if len(same) else np.nan
            else:
                boot_same[family][index] = np.nan
    genome = {
        "subtype": short, "n_donors": len(metadata), "n_shared_genes": n_shared,
        "pearson_ols": pearson,
        "pearson_ols_ci95_low": float(np.nanpercentile(boot_correlation, 2.5)),
        "pearson_ols_ci95_high": float(np.nanpercentile(boot_correlation, 97.5)),
        "spearman_ols": spearman,
        "spearman_ols_ci95_low": float(np.nanpercentile(boot_spearman, 2.5)),
        "spearman_ols_ci95_high": float(np.nanpercentile(boot_spearman, 97.5)),
        "pearson_wls_secondary": pearson_wls,
        "pearson_wls_ci95_low": float(np.nanpercentile(boot_correlation_wls, 2.5)),
        "pearson_wls_ci95_high": float(np.nanpercentile(boot_correlation_wls, 97.5)),
        "spearman_wls_secondary": spearman_wls,
        "spearman_wls_ci95_low": float(np.nanpercentile(boot_spearman_wls, 2.5)),
        "spearman_wls_ci95_high": float(np.nanpercentile(boot_spearman_wls, 97.5)),
        "block_perm_null_mean": float(null_correlation.mean()),
        "block_perm_observed_minus_null": pearson - float(null_correlation.mean()),
        "block_perm_p_greater": float((1 + (null_correlation >= pearson).sum()) / (permutations + 1)),
        "block_perm_n": permutations, "block_bootstrap_n": bootstraps,
        "n_coexpression_blocks": int(len(np.unique(blocks))), "discovery_model": DISCOVERY_MODEL_LABEL,
    }
    targeted = []
    for family in FAMILY_ORDER:
        row = sea.targeted_summary(
            shared_frame, family, shared_membership[family], shared_sign[family], source_size[family]
        )
        observed_same = row["same_direction"]
        observed_fdr = row["same_direction_seaad_fdr"]
        row.update({
            "subtype": short,
            "block_perm_same_null_mean": float(null_same[family].mean()),
            "block_perm_same_observed_minus_null": observed_same - float(null_same[family].mean()),
            "block_perm_same_p_greater": float((1 + (null_same[family] >= observed_same).sum()) / (permutations + 1)),
            "block_perm_fdr_null_mean": float(null_fdr[family].mean()),
            "block_perm_fdr_observed_minus_null": observed_fdr - float(null_fdr[family].mean()),
            "block_perm_fdr_p_greater": float((1 + (null_fdr[family] >= observed_fdr).sum()) / (permutations + 1)),
            "same_direction_fraction_ci95_low": float(np.nanpercentile(boot_same[family], 2.5)),
            "same_direction_fraction_ci95_high": float(np.nanpercentile(boot_same[family], 97.5)),
            "block_perm_n": permutations, "block_bootstrap_n": bootstraps,
            "target_derivation": "young-derived mouse 3 versus 18 months" if family.startswith("mouse_p1") else "mature human 30-plus",
        })
        targeted.append(row)

    sensitivity = []
    arms = [
        ("primary_wls", counts, metadata, wls, wls.attrs["active_columns"]),
        ("ols_same_design", counts, metadata, ols, wls.attrs["active_columns"]),
    ]
    no_path = [column for column in PRIMARY_COLUMNS if column != "adnc_low"]
    arms.append(("wls_without_adnc", counts, metadata, sea.fit_wls(counts, metadata, genes, no_path), no_path))
    single_counts, single_meta, single_genes = sea.load_counts(data_root, short, arm="singleome")
    single_columns = ("age", "sex_male", "mean_log_umi", "adnc_low")
    arms.append(("singleome_wls", single_counts, single_meta,
                 sea.fit_wls(single_counts, single_meta, single_genes, single_columns), single_columns))
    floor_counts, floor_meta, floor_genes = sea.load_counts(data_root, short, floor=50)
    arms.append(("floor_50_nuclei_wls", floor_counts, floor_meta,
                 sea.fit_wls(floor_counts, floor_meta, floor_genes), PRIMARY_COLUMNS))
    for name, arm_counts, arm_meta, fitted, active in arms:
        merged = fitted.merge(discovery[["gene", "mssm_slope_original_ols"]], on="gene", how="left")
        slope_column = "slope" if "slope" in fitted else "slope_ols"
        p_column = "pvalue" if "pvalue" in fitted else "pvalue_ols"
        q_column = "q_bh" if "q_bh" in fitted else "q_bh_ols"
        n, r, rho = sea.correlation_rows(merged, "mssm_slope_original_ols", slope_column)
        arm_membership, arm_sign, arm_size, _ = family_membership(merged["gene"], short)
        both_row = sea.targeted_summary(
            merged, "mature_both_fdr", arm_membership["mature_both_fdr"],
            arm_sign["mature_both_fdr"], arm_size["mature_both_fdr"],
            slope=slope_column, pvalue=p_column, qvalue=q_column,
        )
        sensitivity.append({
            "subtype": short, "arm": name, "status": "tested", "n_donors": len(arm_meta),
            "design": " + ".join(["intercept"] + list(active)), "n_shared_genes": n,
            "pearson_vs_mature_mssm_ols": r, "spearman_vs_mature_mssm_ols": rho,
            **{f"both_fdr_{key}": value for key, value in both_row.items() if key != "family"},
        })
    not_ad_counts, not_ad_meta, not_ad_genes = sea.load_counts(data_root, short, not_ad_only=True)
    try:
        not_ad_fit = sea.fit_wls(not_ad_counts, not_ad_meta, not_ad_genes)
        not_ad_frame = not_ad_fit.merge(discovery[["gene", "mssm_slope_original_ols"]], on="gene", how="left")
        n, r, rho = sea.correlation_rows(not_ad_frame, "mssm_slope_original_ols", "slope")
        sensitivity.append({
            "subtype": short, "arm": "not_ad_only_wls", "status": "tested", "n_donors": len(not_ad_meta),
            "n_shared_genes": n, "pearson_vs_mature_mssm_ols": r, "spearman_vs_mature_mssm_ols": rho,
        })
    except (ValueError, AssertionError) as exc:
        sensitivity.append({
            "subtype": short, "arm": "not_ad_only_wls", "status": "untestable",
            "n_donors": len(not_ad_meta), "reason": str(exc),
        })

    lodo_rows = []
    lodo_gene_frames = []
    both_genes = sorted(load_targets(short)[1])
    lodo_genes = lodo_sensitivity(counts, metadata, genes, both_genes, columns=wls.attrs["active_columns"])
    if len(lodo_genes):
        lodo_genes["subtype"] = short
        lodo_gene_frames.append(lodo_genes)
    for drop in range(len(metadata)):
        keep_donor = np.arange(len(metadata)) != drop
        try:
            dropped = sea.fit_wls(counts[keep_donor], metadata.loc[keep_donor].reset_index(drop=True), genes)
            merged = dropped.merge(discovery[["gene", "mssm_slope_original_ols"]], on="gene", how="left")
            n, r, _ = sea.correlation_rows(merged, "mssm_slope_original_ols", "slope")
            drop_membership, drop_sign, drop_size, _ = family_membership(merged["gene"], short)
            row = sea.targeted_summary(
                merged, "mature_both_fdr", drop_membership["mature_both_fdr"],
                drop_sign["mature_both_fdr"], drop_size["mature_both_fdr"],
            )
            lodo_rows.append({"drop_index": drop, "pearson": r, **row})
        except ValueError:
            continue
    lodo_summary = {
        "subtype": short, "n_successful_donor_drops": len(lodo_rows),
        "pearson_min": min((row["pearson"] for row in lodo_rows), default=np.nan),
        "pearson_max": max((row["pearson"] for row in lodo_rows), default=np.nan),
        "pearson_sign_flips": int(sum(np.sign(row["pearson"]) != np.sign(pearson_wls) for row in lodo_rows)),
        "both_fdr_same_direction_fraction_min": min((row["same_direction_fraction"] for row in lodo_rows), default=np.nan),
        "both_fdr_same_direction_fraction_max": max((row["same_direction_fraction"] for row in lodo_rows), default=np.nan),
        "both_fdr_genes_evaluated": len(lodo_genes),
        "both_fdr_genes_with_any_slope_sign_flip": int((lodo_genes["n_sign_flip"] > 0).sum()) if len(lodo_genes) else 0,
        "both_fdr_genes_losing_nominal_in_any_drop": int((lodo_genes["n_lose_nominal"] > 0).sum()) if len(lodo_genes) else 0,
    }
    sea.log(f"mature analysis {short}: {len(metadata)} donors, {len(frame)} genes, {n_shared} shared")
    return frame, genome, targeted, sensitivity, lodo_summary, lodo_gene_frames


def seaad_bh_audit(genes: pd.DataFrame, genome: pd.DataFrame, targeted: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subtype, frame in genes.groupby("subtype", sort=False):
        for model, p_column, q_column in (
            ("primary_voom_wls", "pvalue", "q_bh"),
            ("ols_same_design", "pvalue_ols", "q_bh_ols"),
        ):
            recomputed = multipletests(frame[p_column].to_numpy(float), method="fdr_bh")[1]
            rows.append({
                "family": model, "subtype": subtype, "n_tests": len(frame),
                "max_abs_q_difference": float(np.max(np.abs(recomputed - frame[q_column].to_numpy(float)))),
            })
    genome_q = multipletests(genome["block_perm_p_greater"], method="fdr_bh")[1]
    rows.append({
        "family": "genome_wide_block_permutation", "subtype": "three_subtype_BH", "n_tests": 3,
        "max_abs_q_difference": float(np.max(np.abs(genome_q - genome["block_perm_q_bh_across_3_subtypes"]))),
    })
    primary = targeted["family"].eq("mature_both_fdr")
    targeted_q = multipletests(targeted.loc[primary, "block_perm_same_p_greater"], method="fdr_bh")[1]
    rows.append({
        "family": "mature_both_fdr_block_permutation", "subtype": "three_subtype_BH", "n_tests": 3,
        "max_abs_q_difference": float(np.max(np.abs(
            targeted_q - targeted.loc[primary, "block_perm_same_q_bh_across_3_subtypes"]
        ))),
    })
    audit = pd.DataFrame(rows)
    assert audit["max_abs_q_difference"].max() < 1e-9, audit.to_string()
    return audit


def analyze(data_root: Path, discovery_source: Path, permutations: int, bootstraps: int):
    plan = sea.load_plan()
    eligibility_validation = validate_eligibility(data_root)
    assert sea.sha256(discovery_source) == json.loads(
        (Path("/home/reza/projects/agescape-revision-data/tmp_organs/followup-11pt-block9-mature-interneuron")
         / "run_summary_mature.json").read_text()
    )["mature_discovery_sha256"]
    frames, genomes, targeted, sensitivities, lodo, lodo_genes = [], [], [], [], [], []
    for short in SUBCLASS:
        outputs = analyze_subtype(data_root, discovery_source, short, permutations, bootstraps)
        frames.append(outputs[0])
        genomes.append(outputs[1])
        targeted.extend(outputs[2])
        sensitivities.extend(outputs[3])
        lodo.append(outputs[4])
        lodo_genes.extend(outputs[5])
    genes = pd.concat(frames, ignore_index=True)
    genome = pd.DataFrame(genomes)
    genome["block_perm_q_bh_across_3_subtypes"] = multipletests(genome["block_perm_p_greater"], method="fdr_bh")[1]
    targeted = pd.DataFrame(targeted)
    primary = targeted["family"] == "mature_both_fdr"
    targeted.loc[primary, "block_perm_same_q_bh_across_3_subtypes"] = multipletests(
        targeted.loc[primary, "block_perm_same_p_greater"], method="fdr_bh"
    )[1]
    audit = seaad_bh_audit(genes, genome, targeted)
    sea.atomic_csv(genes, RESULTS / OUTPUT_NAMES[0])
    sea.atomic_csv(genome, RESULTS / OUTPUT_NAMES[1])
    sea.atomic_csv(targeted, RESULTS / OUTPUT_NAMES[2])
    sea.atomic_csv(pd.DataFrame(sensitivities), RESULTS / OUTPUT_NAMES[3])
    sea.atomic_csv(pd.DataFrame(lodo), RESULTS / OUTPUT_NAMES[4])
    sea.atomic_csv(pd.concat(lodo_genes, ignore_index=True), RESULTS / OUTPUT_NAMES[5])
    audit_name = "followup-mature-interneuron-sea-ad-bh-audit.csv"
    sea.atomic_csv(audit, RESULTS / audit_name)
    output_manifest = {
        "analysis": "Item 7 mature-derived SEA-AD interneuron test",
        "frozen_plan_sha256": plan["frozen_plan_sha256"],
        "script_sha256": sea.sha256(Path(__file__)),
        "mature_discovery_source": str(discovery_source.relative_to(REPO)),
        "mature_discovery_source_sha256": sea.sha256(discovery_source),
        "mature_gene_file": str(TARGET_SOURCE.relative_to(REPO)),
        "mature_gene_file_sha256": sea.sha256(TARGET_SOURCE),
        "source_manifest_sha256": sea.sha256(MANIFEST / "sea-ad-mtg-2024-02-13.json"),
        "eligibility_validation": eligibility_validation,
        "permutations": permutations, "block_bootstraps": bootstraps,
        "permutation_seed": 20261108, "bootstrap_seed": 20261208,
        "discovery_model": DISCOVERY_MODEL_LABEL,
        "mouse_p1_derivation": "young-derived mouse 3 versus 18 months",
        "bh_max_abs_difference": float(audit["max_abs_q_difference"].max()),
        "outputs": {},
    }
    for name in OUTPUT_NAMES + (audit_name,):
        output_manifest["outputs"][name] = {"sha256": sea.sha256(RESULTS / name)}
    (MANIFEST / "sea_ad_mature_interneuron_analysis_outputs.json").write_text(
        json.dumps(output_manifest, indent=2) + "\n"
    )


def fmt(value, digits=4) -> str:
    if value is None or pd.isna(value):
        return "UNTESTABLE"
    return f"{float(value):.{digits}f}"


def scientific(value) -> str:
    return f"{float(value):.3e}"


def integer(value) -> str:
    if value is None or pd.isna(value):
        return "UNTESTABLE"
    return str(int(value))


def intersection(genes: pd.DataFrame) -> list[str]:
    sets = []
    for cell_type in list(CELL_TYPE.values()):
        sets.append(set(genes.loc[genes["cell_type"].eq(cell_type) & genes["tier"].eq("both_FDR"), "gene"]))
    return sorted(set.intersection(*sets))


def build_report() -> None:
    complete = pd.read_csv(RESULTS / "followup-mature-interneuron.csv")
    replication = pd.read_csv(RESULTS / "followup-mature-interneuron-replication.csv")
    full_replication = pd.read_csv(RESULTS / "reviewer-rev2-replication-summary.csv")
    full_de = pd.read_csv(RESULTS / "reviewer-rev2-human-de-summary.csv")
    genes = pd.read_csv(RESULTS / "followup-mature-interneuron-genes.csv")
    moderated_genes = pd.read_csv(RESULTS / "followup-mature-interneuron-genes-moderated.csv")
    six = pd.read_csv(RESULTS / "followup-mature-interneuron-six-gene-status.csv")
    power = pd.read_csv(RESULTS / "followup-mature-interneuron-power.csv")
    gate = pd.read_csv(RESULTS / "followup-mature-interneuron-moderation-gate.csv")
    bh = pd.read_csv(RESULTS / "followup-mature-interneuron-bh-audit.csv")
    genome = pd.read_csv(RESULTS / OUTPUT_NAMES[1])
    targeted = pd.read_csv(RESULTS / OUTPUT_NAMES[2])
    sensitivities = pd.read_csv(RESULTS / OUTPUT_NAMES[3])
    lodo = pd.read_csv(RESULTS / OUTPUT_NAMES[4])
    seaad_genes = pd.read_csv(RESULTS / OUTPUT_NAMES[0])
    seaad_bh = pd.read_csv(RESULTS / "followup-mature-interneuron-sea-ad-bh-audit.csv")
    gate_pass = bool(gate["pass_5_percent"].all())
    mature_unmoderated_intersection = intersection(genes)
    mature_moderated_intersection = intersection(moderated_genes)
    passing = replication[
        replication["gw_pearson_r"].gt(0) & replication["block_perm_p"].lt(0.05)
        & replication["both_fdr_same_direction"].gt(replication["hypergeometric_expected"])
    ]
    all_fit = complete.groupby(["cohort", "cell_type"]).size().size == 8
    if len(passing):
        verdict = "PROVES"
        verdict_reason = "Significant positive genome-wide replication and both-FDR overlap above hypergeometric expectation survives in " + ", ".join(passing["short_type"]) + "."
    elif all_fit:
        verdict = "DOES NOT PROVE"
        verdict_reason = "No subtype meets both frozen replication criteria despite all eight cohort-by-subtype fits having at least 10 mature donors. The per-subtype MDE80 values below bound detectable effects but do not establish equivalence."
    else:
        verdict = "INCONCLUSIVE"
        verdict_reason = "At least one required cohort-by-subtype fit is UNTESTABLE."

    donor_lines = []
    fdr_lines = []
    replication_lines = []
    for cell_type in list(CELL_TYPE.values()) + ["GABAergic neuron"]:
        short = complete.loc[complete["cell_type"].eq(cell_type), "short_type"].iloc[0]
        donor_cells = [short]
        for cohort in ("MSSM_DLPFC", "HBCC_DLPFC"):
            full = full_de[full_de["cohort"].eq(cohort) & full_de["cell_type"].eq(cell_type)].iloc[0]
            mature = complete[complete["cohort"].eq(cohort) & complete["cell_type"].eq(cell_type) & complete["model"].eq("rev2_voom_wls")]
            donor_cells.extend([
                integer(full["n_donors"]), integer(mature["n_donors"].iloc[0]),
                f"{fmt(mature['age_min'].iloc[0], 0)} to {fmt(mature['age_max'].iloc[0], 0)}",
                fmt(mature["age_median"].iloc[0], 1),
            ])
            fdr_lines.append(
                f"| {short} | {cohort} | {integer(full['n_fdr_wls'])} | {integer(full['n_fdr_ols_same_genes'])} | "
                f"{integer(mature['q_bh'].lt(0.05).sum())} | "
                f"{integer(complete[complete['cohort'].eq(cohort) & complete['cell_type'].eq(cell_type) & complete['model'].eq('age_plus_sex_depth_ols')]['q_bh'].lt(0.05).sum())} | "
                f"{integer(complete[complete['cohort'].eq(cohort) & complete['cell_type'].eq(cell_type) & complete['model'].eq('age_plus_sex_ols_depth_unadjusted')]['q_bh'].lt(0.05).sum())} |"
            )
        donor_lines.append("| " + " | ".join(donor_cells) + " |")
        for label, frame in (("Full range", full_replication), ("Mature 30-plus", replication)):
            row = frame[frame["cell_type"].eq(cell_type)].iloc[0]
            replication_lines.append(
                f"| {short} | {label} | {integer(row['n_tested_both'])} | {fmt(row['gw_pearson_r'])} | "
                f"{fmt(row['sel_pearson_r'])} | {integer(row['sel_same_direction'])}/{integer(row['mssm_fdr_genes'])} | "
                f"{integer(row['both_fdr_same_direction'])} | {fmt(row['hypergeometric_expected'], 3)} | "
                f"{fmt(row['matched_background_expected'], 3)} | {fmt(row['block_perm_p'])} |"
            )

    six_lines = []
    for gene in ["FBN1", "GABRB1", "IGF1R", "PHF3", "TRMT11", "ZBTB16"]:
        for cohort in ("MSSM_DLPFC", "HBCC_DLPFC"):
            for short in ("VIP", "pvalb", "sst"):
                row = six[six["gene"].eq(gene) & six["cohort"].eq(cohort) & six["short_type"].eq(short)].iloc[0]
                moderated_status = row["moderated_status"] if gate_pass else f"{row['moderated_status']} (UNVERIFIED)"
                six_lines.append(
                    f"| {gene} | {cohort} | {short} | {row['unmoderated_status']} | {moderated_status} |"
                )
    power_lines = []
    for _, row in power.iterrows():
        power_lines.append(
            f"| {row['short_type']} | {row['cohort']} | {row['age_range']} | {integer(row['n_donors'])} | "
            f"{fmt(row['mde80_slope_per_year'], 6)} | {fmt(row['mde80_slope_per_10_years'], 6)} | "
            f"{integer(row['n_shared_genes'])} | {fmt(row['mde80_pearson_r'])} |"
        )
    genome_lines = []
    for _, row in genome.iterrows():
        genome_lines.append(
            f"| {row['subtype']} | {integer(row['n_shared_genes'])} | "
            f"{fmt(row['pearson_ols'])} ({fmt(row['pearson_ols_ci95_low'])} to {fmt(row['pearson_ols_ci95_high'])}) | "
            f"{fmt(row['pearson_wls_secondary'])} ({fmt(row['pearson_wls_ci95_low'])} to {fmt(row['pearson_wls_ci95_high'])}) | "
            f"{fmt(row['block_perm_p_greater'])} | {fmt(row['block_perm_q_bh_across_3_subtypes'])} |"
        )
    target_lines = []
    for _, row in targeted.iterrows():
        q = fmt(row["block_perm_same_q_bh_across_3_subtypes"]) if row["family"] == "mature_both_fdr" else "secondary"
        target_lines.append(
            f"| {row['subtype']} | {row['family']} | {integer(row['n_tested'])}/{integer(row['n_source_genes'])} | "
            f"{integer(row['same_direction'])} ({fmt(100 * row['same_direction_fraction'], 1)}%) | "
            f"{integer(row['same_direction_nominal'])} | {integer(row['same_direction_seaad_fdr'])} | "
            f"{fmt(row['block_perm_same_p_greater'])} | {q} |"
        )
    sensitivity_lines = []
    for _, row in sensitivities.iterrows():
        same = "UNTESTABLE" if row["status"] != "tested" else f"{integer(row['both_fdr_same_direction'])}/{integer(row['both_fdr_n_tested'])}"
        sensitivity_lines.append(
            f"| {row['subtype']} | {row['arm']} | {row['status']} | {integer(row['n_donors'])} | "
            f"{fmt(row.get('pearson_vs_mature_mssm_ols'))} | {same} |"
        )
    lodo_lines = []
    for _, row in lodo.iterrows():
        lodo_lines.append(
            f"| {row['subtype']} | {integer(row['n_successful_donor_drops'])} | "
            f"{fmt(row['pearson_min'])} to {fmt(row['pearson_max'])} | {integer(row['pearson_sign_flips'])} | "
            f"{fmt(row['both_fdr_same_direction_fraction_min'], 3)} to {fmt(row['both_fdr_same_direction_fraction_max'], 3)} | "
            f"{integer(row['both_fdr_genes_with_any_slope_sign_flip'])}/{integer(row['both_fdr_genes_evaluated'])} |"
        )
    gate_lines = []
    for _, row in gate.iterrows():
        gate_lines.append(
            f"| {row['cohort']} | {row['short_type']} | {integer(row['n_recovered_genes'])} | "
            f"{fmt(row['frozen_s0_sq'], 6)} | {fmt(row['recovered_s0_sq'], 6)} | {fmt(100 * row['s0_sq_relative_difference'], 2)}% | "
            f"{fmt(row['frozen_d0'], 6)} | {fmt(row['recovered_d0'], 6)} | {fmt(100 * row['d0_relative_difference'], 2)}% |"
        )
    comparison_parts = []
    for short in ("VIP", "pvalb", "sst", "GABA"):
        full_row = full_replication[full_replication["short_type"].eq(short)].iloc[0]
        mature_row = replication[replication["short_type"].eq(short)].iloc[0]
        ratios = {}
        for cohort in ("MSSM_DLPFC", "HBCC_DLPFC"):
            unit = power[power["short_type"].eq(short) & power["cohort"].eq(cohort)].set_index("age_range")
            ratios[cohort] = float(unit.loc["mature30", "mde80_slope_per_10_years"] / unit.loc["full_range", "mde80_slope_per_10_years"])
        comparison_parts.append(
            f"{short}: genome-wide r {fmt(full_row['gw_pearson_r'])} to {fmt(mature_row['gw_pearson_r'])}; "
            f"both-FDR {integer(full_row['both_fdr_same_direction'])} to {integer(mature_row['both_fdr_same_direction'])}; "
            f"mature/full-range slope MDE80 ratio {fmt(ratios['MSSM_DLPFC'], 3)} in MSSM and {fmt(ratios['HBCC_DLPFC'], 3)} in HBCC"
        )
    comparison_text = ". ".join(comparison_parts) + "."
    bh_family_lines = []
    for model, frame in bh.groupby("model", sort=False):
        bh_family_lines.append(f"| MSSM/HBCC {model} | {scientific(frame['max_abs_q_difference'].max())} |")
    for family, frame in seaad_bh.groupby("family", sort=False):
        bh_family_lines.append(f"| SEA-AD {family} | {scientific(frame['max_abs_q_difference'].max())} |")
    seaad_fdr = seaad_genes.groupby("subtype")["q_bh"].apply(lambda values: int(values.lt(0.05).sum())).to_dict()
    unmoderated_text = ", ".join(mature_unmoderated_intersection) if mature_unmoderated_intersection else "empty"
    moderated_text = ", ".join(mature_moderated_intersection) if mature_moderated_intersection else "empty"
    moderation_note = "verified" if gate_pass else "UNVERIFIED because the frozen recovery gate failed"
    text = f"""# Item 7 mature interneuron replication and mature-derived SEA-AD test

## Verdict

{verdict}. {verdict_reason} Human mature means donors age 30 years and older in both the discovery fits and the tested MSSM-to-HBCC contrast. The verdict uses the unmoderated WLS replication files. The mature moderated intersection is {moderation_note}.

## Models and fixed families

MSSM and HBCC use donor-summed raw counts, log1p counts per 10,000, one-iteration voom WLS with age plus sex plus donor mean log UMI, and donor-unit inference. QC is primary normal data, at least 1500 UMIs per cell, at least 20 cells per donor and subtype, and at least 10 donors per unit. Depth-adjusted OLS uses the same design. The depth-unadjusted sensitivity uses age plus sex. BH correction is within cohort, subtype, and model over all tested genes. The seed is 20261108 and every block permutation uses 1000 iterations.

The SEA-AD fit is unchanged: 15 donors age 72 to 102, raw UMIs, released QC, Vip/Pvalb/Sst mappings, and WLS age plus sex plus mean log UMI plus multiome fraction plus Low ADNC indicator. Discovery is `{DISCOVERY_MODEL_LABEL}`. SEA-AD OLS is primary for genome-wide correlation and SEA-AD WLS is secondary. Targeted same-direction tests use WLS. Permutations and block bootstraps each use 1000 iterations with seeds 20261108 and 20261208.

## Donor retention

| Subtype | MSSM full N | MSSM mature N | MSSM mature ages | MSSM median | HBCC full N | HBCC mature N | HBCC mature ages | HBCC median |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(donor_lines)}

## FDR counts

| Subtype | Cohort | Full WLS | Full depth-adjusted OLS | Mature WLS | Mature depth-adjusted OLS | Mature depth-unadjusted OLS |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(fdr_lines)}

## MSSM-to-HBCC replication

Selected same-direction counts are among MSSM FDR genes. Both-FDR counts require the same sign in both cohorts. Hypergeometric and decile-matched values are expected counts. Block p is the coexpression-block permutation p for the both-FDR overlap.

| Subtype | Range | Shared genes | Genome-wide r | Selected r | Selected same/tested | Both FDR | Hypergeometric expected | Matched expected | Block p |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(replication_lines)}

The original full-range moderated three-subtype intersection is exactly FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16. DYNC1I1 occurs only in the unmoderated full-range intersection. The mature unmoderated three-subtype intersection is: {unmoderated_text}. The mature moderated three-subtype intersection is: {moderated_text}. The latter is {moderation_note}.

## Six-gene mature status

PASS means q < 0.05, NOT_FDR means tested with q >= 0.05, and NOT_TESTED means absent after the unchanged expression filter. Unmoderated statuses are the confirmed mature result. Moderated statuses are marked UNVERIFIED if the recovery gate failed and are not presented as confirmed replication.

| Gene | Cohort | Subtype | Unmoderated WLS | Moderated WLS |
|---|---|---|---|---|
{chr(10).join(six_lines)}

## Power and interpretation

Slope MDE80 is `(t_0.975 + t_0.80) * median WLS SE` on the log1p CP10K scale. Correlation MDE80 uses the Fisher z approximation with two-sided alpha 0.05 and the shared-gene count. Full-range median SE values come from the frozen all-cell-type audit, not a refit.

| Subtype | Cohort | Range | N | MDE80 per year | MDE80 per 10 years | Shared genes | Correlation MDE80 |
|---|---|---|---:|---:|---:|---:|---:|
{chr(10).join(power_lines)}

Versus the original full-range analysis, {comparison_text} HBCC donor loss therefore comes with a real slope-precision cost, while genome-wide correlation MDE80 remains near 0.065 to 0.066 because shared-gene counts stay near 1800. VIP and sst retain large both-FDR excesses despite fewer HBCC hits. A non-significant correlation is not evidence of no aging, and a small MDE80 is not an equivalence test.

## Mature-derived SEA-AD outcome

### Genome-wide effects

Intervals are 95 percent block-bootstrap percentile intervals. Block p is one-sided for positive OLS replication and BH q is across the three subtypes.

| Subtype | Shared genes | OLS Pearson r (95% CI) | WLS Pearson r (95% CI) | Block p | BH q |
|---|---:|---:|---:|---:|---:|
{chr(10).join(genome_lines)}

Genome-wide SEA-AD FDR gene counts are {integer(seaad_fdr['VIP'])} VIP, {integer(seaad_fdr['pvalb'])} pvalb, and {integer(seaad_fdr['sst'])} sst. These counts cover all filtered genes within each SEA-AD subtype.

### Target families

The three mature human families use mature unmoderated both-FDR and mature MSSM discovery sets. The mouse P1 families are young-derived from mouse 3 versus 18 months, not mature-derived targets. Same-direction SEA-AD FDR is the subset with SEA-AD q < 0.05. Only genome-wide correlation and mature both-FDR same-direction enrichment have primary three-subtype BH correction; other set-level p values are secondary.

| Subtype | Family | Tested/source | Same direction | Nominal | SEA-AD FDR | Block p | BH q |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(target_lines)}

### Sensitivities and donor influence

| Subtype | Arm | Status | N | Pearson r vs mature MSSM OLS | Both-FDR same/tested |
|---|---|---|---:|---:|---:|
{chr(10).join(sensitivity_lines)}

| Subtype | Successful drops | WLS correlation range | Sign flips | Both-FDR fraction range | Gene sign flips/evaluated |
|---|---:|---:|---:|---:|---:|
{chr(10).join(lodo_lines)}

This is a limited, model-sensitive result in 15 old donors, not evidence of no aging. The narrow age span, old-age range, MTG region, assay, donor count, residual pathology, and estimator sensitivity limit generalization. Not-AD-only fits remain untestable below the 10-donor floor. Mouse P1 sets are too sparse for a cross-species verdict when only a few genes pass the SEA-AD filter.

## Verification and deviations

The new complete MSSM/HBCC model families and SEA-AD gene effects were independently BH-recomputed from raw p-values. Maximum absolute q differences were {scientific(bh['max_abs_q_difference'].max())} for MSSM/HBCC and {scientific(seaad_bh['max_abs_q_difference'].max())} for SEA-AD, both below 1e-9. Selected replication-list q values retain their full-unit BH families; BH was not reapplied after selecting genes. Exact per-cohort, per-subtype family values are in the two BH audit CSVs.

| Model or correction family | Maximum absolute q difference |
|---|---:|
{chr(10).join(bh_family_lines)}

The frozen moderation recovery gate used every frozen unit and the exact recovered-s2 formula. It requires both relative differences to be at most 5 percent. The gate passed: {gate_pass}. Frozen gene files contain selected genes, so the recovered prior fit uses the available selected-gene rows rather than a complete tested-gene vector. This limitation and the observed discrepancies are reported rather than hidden.

| Cohort | Subtype | Recovered genes | Frozen s0_sq | Recovered s0_sq | Relative difference | Frozen d0 | Recovered d0 | Relative difference |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(gate_lines)}

The WLS copy returning final weighted residual s2 reproduced cached slopes, SEs, p values, and q values to 1e-9 before mature moderation. Mature priors were re-estimated separately for every mature unit. Because the gate failed, moderated mature intersections and six-gene statuses are unverified sensitivities; the unmoderated intersection is used for the confirmed result. WLS bootstrap intervals were added for the required SEA-AD report while the original primary OLS analysis was unchanged. No frozen original was recomputed or overwritten, no donor object was downloaded or re-hashed, and no original SEA-AD output was modified.
"""
    write_clean_markdown(text, REPORT)
    manifest_path = MANIFEST / "sea_ad_mature_interneuron_analysis_outputs.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["script_sha256"] = sea.sha256(Path(__file__))
    manifest["report_sha256"] = sea.sha256(REPORT)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["analyze", "report", "all"], default="all")
    parser.add_argument("--data-root", type=Path, default=REPO.parent / "agescape-revision-data" / "sea-ad")
    parser.add_argument("--discovery-source", type=Path, default=DISCOVERY_SOURCE)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstraps", type=int, default=1000)
    args = parser.parse_args()
    assert args.permutations == 1000 and args.bootstraps == 1000
    started = time.time()
    if args.stage in {"analyze", "all"}:
        analyze(args.data_root, args.discovery_source, args.permutations, args.bootstraps)
    if args.stage in {"report", "all"}:
        build_report()
    sea.log(f"mature SEA-AD stage {args.stage} completed in {time.time() - started:.1f} seconds")


if __name__ == "__main__":
    main()
