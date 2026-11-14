#!/usr/bin/env python3
"""Render the Cui testis metadata/QC report from committed aggregate tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
PREFIX = "followup-cui-testis"
CAVEAT = (
    "**Biological caveat for every result:** changing germ-cell-stage composition can mimic "
    "broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls "
    "cannot support general whole-body aging claims."
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", default="PENDING")
    args = parser.parse_args()
    plan = json.loads((REPO / "data" / "manifest" / "cui_testis_analysis_plan.json").read_text())
    manifest = json.loads((REPO / "data" / "manifest" / "GSE254315.json").read_text())
    donors = pd.read_csv(RESULTS / f"{PREFIX}-donor-map.csv")
    included = donors[donors["included_in_inferential_n"] == 1].copy()
    qc = pd.read_csv(RESULTS / f"{PREFIX}-qc.csv")
    included = included.merge(qc[["donor", "raw_barcodes", "paper_qc_pass_cells"]], on="donor", validate="one_to_one")
    status = pd.read_csv(RESULTS / f"{PREFIX}-endpoint-status.csv")
    sperm = pd.read_csv(RESULTS / f"{PREFIX}-spermatogenesis-metadata.csv")
    validation = json.loads((RESULTS / f"{PREFIX}-validation.json").read_text())
    resources = json.loads((RESULTS / f"{PREFIX}-resource-use.json").read_text())
    inventory = pd.read_csv(RESULTS / f"{PREFIX}-file-inventory.csv")

    lines = [
        "# Cui testis GSE254315: mapping, frozen analysis plan, and descriptive QC",
        "",
        "Base revision: `8cd58bdfc0921636ba2846f898fc6d0ab96d3ccb`. Metadata and plan freeze: 2026-11-14.",
        "",
        CAVEAT,
        "",
        "## Decision",
        "",
        "The mandatory mapping is **resolved for GSE254315**: its 25 GEO records are 25 separate named donor libraries, not pools and not 35 independent measurements. The paper's 35-donor integrated atlas comprises these 25 plus ten reused external-only donor IDs. This step uses only the 25 GSE donors. No public donor-and-barcode cell-type crosswalk was found, so P1, donor/cell-type pseudobulk, CLR composition, and Hellinger analyses were **not run**. Raw-count validation and descriptive QC were completed; no expression outcome or age association was fitted.",
        "",
        "Primary sources: [Cui et al., Nature Aging (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12003174/), [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/soft/GSE254315_family.soft.gz), [GEO file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/suppl/filelist.txt), and the primary Fig. 1a/Extended Data Fig. 1/Fig. 6 source workbooks recorded in `data/manifest/GSE254315.json`.",
        "",
        "## The 25-versus-35 mapping",
        "",
        CAVEAT,
        "",
        "The primary Fig. 1a source is ordered Donor1 through Donor35. Its ages match GEO exactly for 24 of the 25 GSE donors. The sole discrepancy is **Donor8: 40 years in the primary source versus 39 in GEO**. The frozen primary uses 40/40s and declares a complete 39/30s sensitivity. The ten atlas-only donor IDs are Donor1–3, Donor26–29, and Donor31–33. The paper names GSE182786 and GSE215754 as external sources; those donors are excluded here and their ages are never added as independent observations.",
        "",
        "| Donor | GSM | Library | Primary age/decade | GEO age | BMI | Raw barcodes | Paper-QC pass |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in included.itertuples():
        lines.append(
            f"| {row.donor} | {row.gsm} | {row.library} | {row.primary_age}/{row.primary_decade} "
            f"| {row.geo_age} | {row.bmi:g} | {row.raw_barcodes} | {row.paper_qc_pass_cells} |"
        )
    decade_counts = included.groupby("primary_decade").size().to_dict()
    lines += [
        "",
        f"The primary decade distribution is 30s={decade_counts['30s']}, 40s={decade_counts['40s']}, 50s={decade_counts['50s']}, 60s={decade_counts['60s']}; N=25. The 20s donors in the integrated figure are external-only and are not in GSE254315.",
        "",
        "## Frozen pool/donor-level plan",
        "",
        CAVEAT,
        "",
        "The plan in `data/manifest/cui_testis_analysis_plan.json` was frozen after mapping and before expression download or outcomes. Its SHA-256 is recorded in the metadata manifest. If a valid annotation crosswalk is later released, the primary is an ordinal-decade effect (3/4/5/6) per decade with donor mean log actual raw UMI adjustment over all 25 donors. The 30s are young/reference, 40s are middle, and 50s/60s are older decades.",
        "",
        "Declared sensitivities are: unadjusted ordinal; BMI-plus-depth ordinal; categorical decades with 30s reference; old 50s/60s versus young 30s with middle 40s excluded only there; exact age per ten years; Donor8=39/30s versions of the applicable models; and incomplete-spermatogenesis-adjusted and complete-only models **only if** a released donor-linked status exists. No exact pool ages are invented and no cell-level p-values are permitted.",
        "",
        "Eligible released labels would require at least 20 cells per donor/cell type, six donors overall, three per retained group, full rank, and two residual degrees of freedom. Hellinger additionally requires at least two young reference donors, 20 qualifying cells per scored donor, 50 cross-fitted reference cells, and actual-UMI >=1500. Somatic labels and each comparable germ-cell stage are separate. P1 uses the established human P1_Clean and a 200-draw expression-decile-matched null (seed 20261114). Effects require standard errors, 95% CIs, and explicit BH families as frozen in the plan.",
        "",
        "## Annotation gate and blocked endpoints",
        "",
        CAVEAT,
        "",
        "GEO releases raw MatrixMarket counts, barcodes, and features only. Primary Fig. 1a releases donor age/BMI, Extended Data Fig. 1 releases age-aggregated cell totals, and Fig. 6a releases donor-linked elongating-spermatid frequencies. None is a per-cell donor/barcode annotation crosswalk. The released code repository contains notebooks that refer to local metadata and H5AD objects but no public annotation files. Marker-based relabeling, reclustering, or treating age-aggregated totals as donor observations would violate the input rules.",
        "",
        "| Endpoint | Somatic labels | Comparable germ-cell stages | Status | Effect/CI/p/q |",
        "|---|---:|---:|---|---|",
    ]
    for endpoint in ["P1", "donor_pseudobulk", "CLR", "Hellinger"]:
        part = status[status["endpoint"] == endpoint]
        lines.append(
            f"| {endpoint} | {(part['compartment'] == 'somatic').sum()} "
            f"| {(part['compartment'] == 'germ_cell_stage').sum()} | blocked; not run | none |"
        )
    lines += [
        "",
        "All 56 label-by-endpoint rows and their failed annotation reason are in `results/followup-cui-testis-endpoint-status.csv`. They are blocked, not null results, and enter no BH denominator. No whole-testis pseudobulk model was substituted because it could not separate cell-intrinsic changes from somatic/germ-stage composition.",
        "",
        "## Raw-count validation and descriptive QC",
        "",
        CAVEAT,
        "",
        f"All 75 per-donor files were downloaded sequentially under `scripts/resource_guard.py` with a 30 GB reserve and 0.2 GB per-file cap; SHA-256 and expected byte sizes passed. All stored count values in all 25 matrices are finite, nonnegative integers. Matrix dimensions match barcodes/features; barcodes and feature IDs are unique. Decompressed feature content is identical across all donors (canonical SHA-256 `{validation['feature_reference_sha256']}`), although gzip-container checksums can differ.",
        "",
        f"The matrices contain **{validation['raw_barcodes_total']:,} raw barcodes** and **{qc['raw_umi_sum'].sum():,} UMIs**, whereas the GEO summary reports 99,130 cells. Applying the paper's released thresholds descriptively (200–5500 detected features, 800–17000 UMI, mitochondrial fraction <=25%) retains **{qc['paper_qc_pass_cells'].sum():,}** barcodes ({100 * qc['paper_qc_pass_cells'].sum() / qc['raw_barcodes'].sum():.2f}%), **{99130 - qc['paper_qc_pass_cells'].sum():,} fewer** than the GEO summary. This difference is unresolved because final annotated barcodes are not public; it is not silently treated as a successful reproduction. These are not reclustered or annotated cells and were not used for age-associated QC tests. Per-donor UMI, feature, mitochondrial, and QC summaries are in `results/followup-cui-testis-qc.csv`.",
        "",
        f"The guarded QC process had maximum RSS {resources['max_rss_kib'] / 1024:.1f} MiB, below the 8 GB cap; one donor was processed at a time. External checkpoints are under `../agescape-revision-data/cui-testis/interim/qc/`.",
        "",
        "## Spermatogenesis information",
        "",
        CAVEAT,
        "",
        f"Primary Fig. 6a provides donor-linked elongating-spermatid `Freq` values for {sperm['released_freq_status'].eq('released').sum()} of 25 GSE donors; Donor34 is absent and remains missing in the aggregate CSV. The Methods define infertility as fewer than three sperm per seminiferous tubule, but the checked public source does not release donor-level PAS sperm counts or an incomplete-spermatogenesis indicator. Therefore the spermatogenesis-adjusted and complete-only sensitivity models are **untestable**, and no status was inferred from germ-cell frequencies. Released age-aggregated label totals are retained descriptively in `results/followup-cui-testis-released-age-cell-totals.csv`, not treated as independent donors.",
        "",
        "## What survived, died, and remains blocked",
        "",
        CAVEAT,
        "",
        "- Survived: one-to-one 25-donor GEO mapping; exact primary/GEO age reconciliation with one declared discrepancy; raw-count integrity; descriptive QC; and release of descriptive age/BMI and germ-stage source fields.",
        "- Died as an inferential shortcut: assigning all 35 atlas ages to GSE254315, treating cells as replicates, or treating the ten external donors as new libraries.",
        "- Blocked: every cell-resolved P1, gene-effect, CLR, and Hellinger estimand, including their matched-null and sensitivity models, because the required released annotation crosswalk is absent. Non-significance is neither reported nor implied.",
        "- Contradiction/complication: the 35-donor paper atlas is not a 35-library GSE cohort; Donor8 also has a one-year primary-versus-GEO discrepancy. These do not justify guessing labels or ages.",
        "- Not attempted: no FASTQ/RAW tar, external-donor analysis, clustering, marker imputation, expression outcome modeling, cell-level tests, or ovary work.",
        "",
        "Reproduce from the repository in the `agescape` environment:",
        "",
        "```bash",
        "python scripts/prepare_followup_cui_testis.py",
        "python scripts/fetch_followup_cui_testis.py",
        "python scripts/resource_guard.py run --path ../agescape-revision-data/cui-testis --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_cui_testis.py",
        "python scripts/report_followup_cui_testis.py --tests '<full-suite result>'",
        "```",
        "",
        "## Verification and required follow-up",
        "",
        f"Full-suite verification: `{args.tests}`. No commit or push was made by the testis agent.",
        "",
        "The remaining prerequisite is a **public, released donor-and-barcode-to-cell-type annotation crosswalk** (or annotated public object with raw counts) covering these libraries. If supplied, the frozen design can be evaluated without inventing labels. Any design change requires a new freeze before outcomes. Donor mapping itself does not require further user input.",
        "",
    ]
    (REPO / "reports" / "followup-cui-testis.md").write_text("\n".join(lines))
    print("Wrote reports/followup-cui-testis.md")


if __name__ == "__main__":
    main()
