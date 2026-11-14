#!/usr/bin/env python3
"""Freeze the Step 7 ovary analysis plans before expression acquisition."""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "data" / "manifest"
HEAD = "0610c34bde6c176efc1c225cc8d817231e9ce54a"

COMMON = {
    "inferential_unit": "donor",
    "cell_level_inference": False,
    "cross_cohort_pooling": False,
    "premenopausal_claims": False,
    "causal_claims": False,
    "minimum_cells_per_donor_cell_type": 20,
    "released_annotations_only": True,
    "minimum_residual_degrees_of_freedom": 2,
    "confidence_interval": "two-sided donor-level OLS t 95 percent",
    "p1_signature": "Established human P1_Clean from src.aging_programs.load_program_tiers",
    "p1_summary": "Mean per-cell log1p(CPM) within donor and released broad cell type; detected P1 up mean minus detected P1 down mean",
    "p1_matched_null": "200 expression-decile-matched random-signature draws without replacement, excluding detected P1 genes, seed 20261105; centered two-sided empirical p=(1+extreme)/201",
    "pseudobulk": "Sum released corrected/raw UMI counts across all pieces and cells within donor-cell type, then log1p(CPM); filter to >=10 summed counts and detection in >=2 included donors without using age outcomes",
    "clr": "All released non-unknown coarse labels define the donor composition universe, including zero counts; 0.5-cell pseudocount",
    "effect_reporting": "Effect, standard error, 95% CI, nominal p and BH q; non-significance is not equivalence",
}

RESOURCES = {
    "script": "scripts/resource_guard.py",
    "minimum_free_gb": 30,
    "analysis_memory_gb": 8,
    "sequential_downloads_and_heavy_jobs": True,
    "backed_or_streamed_expression_reads": True,
    "raw_tars_permitted": False,
    "fastq_permitted": False,
}

JIN = {
    "created": "2026-11-14",
    "step": 7,
    "cohort": "GSE202601",
    "execution_order": 1,
    "frozen_at_repository_head": HEAD,
    "models_frozen_before_expression_download_or_outcomes": True,
    "current_plan_frozen_before_outcomes": True,
    "amendment_scope": "Only join/count-transformation gates were clarified after expression schema inspection and before outcomes; donor sets, endpoints, models, covariates, floors, and BH families were not changed",
    "resources": {**RESOURCES, "maximum_individual_download_gb": 1},
    "common": COMMON,
    "input": "Eight public per-donor Cell Ranger filtered feature-barcode H5 files and released snRNA metadata celltype2 broad annotations; no ATAC expression/peak download",
    "count_and_join_gate": "Stream every stored matrix value and require finite nonnegative integers; unique feature IDs and per-donor barcodes; exact one-to-one donor-specific H5/metadata barcode reconciliation after converting each matrix '-1' suffix to the released donor-specific metadata suffix '.1' through '.8'. Raw H5 matrices may contain additional unannotated barcodes; every released annotated barcode must match uniquely and only annotated nuclei enter biological summaries. This join clarification was frozen before outcomes after schema inspection.",
    "donors": {
        "young": ["ovary1", "ovary2", "ovary3", "ovary4"],
        "reproductively_aged": ["ovary5", "ovary6", "ovary7", "ovary8"],
        "each_donor_used_once": True,
        "continuous_age_model": "blocked_not_planned",
        "ovary5_exact_age": "UNRESOLVED_47_SNRNA_49_GEO_AND_SNATAC; group retained as reproductively aged; numeric age never used",
    },
    "eligibility": "Every one of four young and four reproductively-aged donors must have >=20 annotated nuclei in the broad cell type; otherwise all cell-type inferential endpoints are untestable",
    "models": {
        "primary_group_depth_adjusted": "intercept + reproductively_aged + mean_log_umi",
        "unadjusted_group_sensitivity": "intercept + reproductively_aged",
        "additional_covariates": "none; all donors female and no usable donor-level procurement/PMI measurements are released",
        "effect": "reproductively-aged minus young",
    },
    "endpoints": {
        "P1_transport": {"models": ["primary_group_depth_adjusted", "unadjusted_group_sensitivity"], "matched_null": "primary model; conventional and matched-null p-values receive separate BH families"},
        "donor_pseudobulk_gene": {"models": ["primary_group_depth_adjusted", "unadjusted_group_sensitivity"], "multiplicity": "BH across every tested cell-type-by-gene combination within model"},
        "donor_CLR_composition": {"models": ["primary_group_depth_adjusted", "unadjusted_group_sensitivity"], "depth": "mean log raw UMI over all annotated donor nuclei"},
        "Hellinger_displacement": {
            "models": ["primary_group_depth_adjusted", "unadjusted_group_sensitivity"],
            "reference": "Cross-fitted young-group raw-count reference within broad cell type after technical-gene exclusions; leave evaluated young donor out and use all eligible young donors for aged donors",
            "cell_filter": ">=1500 raw UMI and positive nontechnical UMI",
            "technical_exclusions": "MT-, RPL/RPS/MRPL/MRPS prefixes and XIST, TSIX, DDX3Y, EIF2S3Y, UTY, KDM5D",
            "floors": "20 qualifying nuclei per scored donor, >=50 cross-fitted young-reference nuclei per donor, and all four donors per group after filtering",
            "otherwise": "untestable; never substitute zero or a null result",
        },
        "same_nucleus_RNA_ATAC_coupling": {"status": "blocked_separate_assays_and_biosamples"},
        "donor_matched_RNA_ATAC_linkage": {"status": "blocked_no_valid_released_peak_gene_link_table_and_no_validated_EC_BEC_LEC_harmonization"},
    },
    "multiplicity": "BH separately by endpoint and model; P1 conventional and matched-null families separate; pseudobulk BH is across all cell-type-by-gene tests in each model; sensitivity families are not pooled with primary families",
}

POST = {
    "created": "2026-11-14",
    "step": 7,
    "cohort": "GSE329417",
    "execution_order": 2,
    "frozen_at_repository_head": HEAD,
    "models_frozen_before_expression_download_or_outcomes": True,
    "current_plan_frozen_before_outcomes": True,
    "amendment_scope": "Only join/count-transformation gates were clarified after expression schema inspection and before outcomes; donor sets, endpoints, models, covariates, floors, and BH families were not changed",
    "resources": {**RESOURCES, "maximum_compressed_atlas_download_gb": 7, "maximum_decompressed_atlas_gb": 40},
    "common": COMMON,
    "publication_status": "bioRxiv preprint; no peer-reviewed journal publication verified",
    "input": "Public final annotated atlas H5AD with released annotations and layers/raw_counts CellBender-corrected unnormalized UMI counts; public donor and processing metadata",
    "count_and_join_gate": "Stream every stored layers/raw_counts value and require finite nonnegative integers; unique feature IDs and observation IDs; exact donor/sample/plex reconciliation to the 64-row processing log and 28-donor clinical table; use only released annotations. Atlas unique_sample_ID values are reconciled one-to-one through each released per-sample CellBender H5AD filename stem and its _p1-_p4 suffix, not fuzzy joins. This filename-crosswalk clarification was frozen before outcomes after schema inspection.",
    "count_processing": {
        "primary": "layers/raw_counts CellBender-corrected counts",
        "x_check": "Require atlas X to equal layers/log1p_norm in streamed chunks; X is log-normalized and is never used as raw counts. Only layers/raw_counts enters endpoints. This transformation clarification was frozen before outcomes after direct HDF5 inspection.",
        "uncorrected_sensitivity": "blocked: original uncorrected Cell Ranger count matrices are not separately released",
    },
    "donor_aggregation": "Sum counts and average per-cell summaries over all annotated nuclei across repeated tissue pieces within donor-cell type; one donor row per endpoint. SN_0086 samples 3_4 and 16_4 are pooled as technical replicates within that donor aggregation, never treated as independent samples",
    "eligibility": ">=20 annotated nuclei per donor-cell type and >=20 eligible donors; each model must also be finite, full rank and have >=2 residual degrees of freedom",
    "technical_covariates": "mean_log_umi plus nuclei-weighted fractions of the donor-cell type in plexes 1, 2 and 3, with plex 4 as reference; composition uses all annotated donor nuclei for depth and plex fractions",
    "clinical_covariates": "BMI continuous and non_ovarian_abnormal indicator (released Pre-malignant or Malignant versus Benign); preserve original three-state labels in donor output",
    "models": {
        "primary_depth_plex_clinical": "intercept + exact_age_per_10y + mean_log_umi + plex_fraction_1 + plex_fraction_2 + plex_fraction_3 + BMI + non_ovarian_abnormal",
        "technical_depth_plex_sensitivity": "intercept + exact_age_per_10y + mean_log_umi + plex_fraction_1 + plex_fraction_2 + plex_fraction_3",
        "unadjusted_sensitivity": "intercept + exact_age_per_10y",
        "secondary_primary_plus_parity": "primary model + parity_live_births; age remains the inferential effect and parity adjustment is secondary",
        "effect": "within-postmenopausal exact-age association per 10 years",
    },
    "endpoints": {
        "P1_transport": {"models": ["primary_depth_plex_clinical", "technical_depth_plex_sensitivity", "unadjusted_sensitivity", "secondary_primary_plus_parity"], "matched_null": "primary model; conventional and matched-null p-values receive separate BH families"},
        "donor_pseudobulk_gene": {"models": ["primary_depth_plex_clinical", "technical_depth_plex_sensitivity", "unadjusted_sensitivity", "secondary_primary_plus_parity"], "multiplicity": "BH across every tested cell-type-by-gene combination within model"},
        "donor_CLR_composition": {"models": ["primary_depth_plex_clinical", "technical_depth_plex_sensitivity", "unadjusted_sensitivity", "secondary_primary_plus_parity"]},
        "Hellinger_displacement": {"status": "untestable_no_20_to_40_year_reference", "outcomes_examined": False},
        "premenopausal_or_Jin_pooled_analysis": {"status": "blocked_by_design", "outcomes_examined": False},
    },
    "multiplicity": "BH separately by endpoint and model; P1 conventional and matched-null families separate; pseudobulk BH is across all cell-type-by-gene tests in each model; sensitivity and secondary parity-adjusted families are not pooled with primary families",
}


def main():
    for name, plan in [("ovary_jin_analysis_plan.json", JIN), ("ovary_postmenopausal_analysis_plan.json", POST)]:
        path = MANIFEST / name
        path.write_text(json.dumps(plan, indent=2) + "\n")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
