#!/usr/bin/env python3
"""Render the PBMC analysis report from committed aggregate result tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
PREFIX = "followup-pbmc"
COHORTS = ["GSE213516", "GSE214546"]


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"{PREFIX}-{name}.csv")


def number(value, digits=4) -> str:
    if pd.isna(value):
        return "—"
    return f"{value:.{digits}g}"


def effect(row) -> str:
    if row.get("status") != "fit":
        return "—"
    return f"{number(row['effect'])} [{number(row['ci95_low'])}, {number(row['ci95_high'])}]"


def model_table(frame: pd.DataFrame, label_columns: list[str]) -> str:
    lines = [*label_columns, "Effect [95% CI]", "p", "q", "N", "Status/reason"]
    output = ["| " + " | ".join(lines) + " |", "|" + "|".join(["---"] * len(lines)) + "|"]
    for _, row in frame.iterrows():
        status = "fit" if row.get("status") == "fit" else f"{row.get('status', '')}: {row.get('reason', '')}"
        values = ["" if pd.isna(row.get(column, "")) else str(row.get(column, "")) for column in label_columns]
        values += [effect(row), number(row.get("pvalue")), number(row.get("q_bh")), number(row.get("n_donors")), status]
        output.append("| " + " | ".join(values) + " |")
    return "\n".join(output)


def family_counts(frame: pd.DataFrame, cohort: str, design: str, model: str) -> tuple[int, int]:
    part = frame[(frame["cohort"] == cohort) & (frame["design"] == design)
                 & (frame["model"] == model) & (frame["status"] == "fit")]
    return len(part), int((part["q_bh"] < 0.05).sum())


def main() -> None:
    p1 = read("p1-models")
    pseudo = read("pseudobulk-targeted")
    comp = read("clr-targeted")
    hell = read("hellinger-models")
    clr = read("clr-models")
    eligibility = read("eligibility-all")
    resources = json.loads((RESULTS / f"{PREFIX}-resource-use.json").read_text())
    validation = json.loads((RESULTS / f"{PREFIX}-validation-summary.json").read_text())
    labels = read("label-harmonization")
    identifiability = read("covariate-identifiability")
    strict_p1 = p1[(p1["design"] == "strict_group") & (p1["model"] == "primary")].sort_values(["cell_type", "cohort"])
    strict_hell = hell[(hell["design"] == "strict_group") & (hell["model"] == "primary")].sort_values(["cell_type", "cohort"])
    strict_pseudo = pseudo[(pseudo["design"] == "strict_group") & (pseudo["model"] == "primary")].sort_values(["gene", "cohort"])
    strict_comp = comp[(comp["design"] == "strict_group") & (comp["model"] == "primary")].sort_values(["common_cell_type", "cohort"])
    sensitivity = pd.concat([
        pseudo[pseudo["design"].isin(["atlas_group", "exact_age"])],
        comp[comp["design"].isin(["atlas_group", "exact_age"])],
    ], ignore_index=True).sort_values(["design", "common_cell_type", "gene", "cohort"], na_position="first")
    thomson_covariates = pseudo[(pseudo["cohort"] == "GSE214546") & pseudo["design"].eq("strict_group")
                                & pseudo["model"].isin(["depth", "cmv", "batch", "pool"])].sort_values(["gene", "model"])
    thomson_p1_strict = p1[(p1["cohort"] == "GSE214546") & p1["design"].eq("strict_group")]
    p1_survivors = {model: int((part["q_bh"] < 0.05).sum()) for model, part in thomson_p1_strict.groupby("model")}
    p1_batch_min_q = float(thomson_p1_strict.loc[thomson_p1_strict["model"].eq("batch"), "q_bh"].min())
    thomson_clr_covariates = clr[(clr["cohort"] == "GSE214546") & clr["design"].eq("strict_group")
                                 & clr["endpoint"].eq("CLR_composition")
                                 & clr["cell_type"].isin(["CD8 Naive", "CD16+ Monocytes"])
                                 & clr["model"].isin(["log_total_cells", "cmv", "batch", "pool"])].copy()
    thomson_clr_covariates = thomson_clr_covariates.sort_values(["cell_type", "model"])
    thomson_hell_strict = hell[(hell["cohort"] == "GSE214546") & hell["design"].eq("strict_group")]
    thomson_hell_batch_hits = thomson_hell_strict[(thomson_hell_strict["model"].eq("batch_depth"))
                                                  & (thomson_hell_strict["q_bh"] < 0.05)]
    thomson_hell_cmv_hits = thomson_hell_strict[(thomson_hell_strict["model"].eq("cmv_depth"))
                                                & (thomson_hell_strict["q_bh"] < 0.05)]
    exploratory_counts = []
    for cohort in COHORTS:
        path = RESULTS / f"{PREFIX}-{cohort}-pseudobulk-exploratory.csv.gz"
        for chunk in pd.read_csv(path, chunksize=300000):
            grouped = chunk.groupby(["design", "model"]).agg(tested=("pvalue", "size"), survivors=("q_bh", lambda x: int((x < 0.05).sum())))
            exploratory_counts.extend((cohort, design, model, *row) for (design, model), row in grouped.iterrows())
    exploratory = pd.DataFrame(exploratory_counts, columns=["cohort", "design", "model", "tested", "survivors"])
    exploratory = exploratory.groupby(["cohort", "design", "model"], as_index=False)[["tested", "survivors"]].sum()
    ineligible = eligibility[~eligibility["eligible"]].copy()
    untestable_models = pd.concat([
        p1[p1["status"] != "fit"], pseudo[pseudo["status"] != "fit"], comp[comp["status"] != "fit"],
        hell[hell["status"] != "fit"],
    ], ignore_index=True, sort=False)
    lines = [
        "# PBMC expression analysis",
        "",
        "Analysis base revision at plan freeze: `5e8bc6fff4fbb372c80c89338fe2c1ab6c929cb2`. The statistical plan was frozen in `data/manifest/pbmc_analysis_plan.json` before expression download.",
        "",
        "## Scope and main verdict",
        "",
        "Only **GSE213516 (Zhu)** and **GSE214546 (Thomson)** were analyzed, sequentially and independently. GSE157007 remained descriptive/excluded. No matrices were concatenated across cohorts, no p-values were pooled across studies, and no cell-level inferential test was used. Strict young `<=40` versus old `>60` is primary; atlas age `>=60` and exact-age models are sensitivities.",
        "",
        f"- **Robustness warning—Thomson strict P1:** primary {p1_survivors['primary']}/13, depth {p1_survivors['depth']}/13, CMV {p1_survivors['cmv']}/13, pool {p1_survivors['pool']}/13, but batch {p1_survivors['batch']}/13 cell types have q<0.05 (minimum batch q={p1_batch_min_q:.6f}). **Broad P1 is not robust to batch and must not be presented as confirmed aging biology.**",
        "- **Composition replicated in direction:** CD8-naive T-cell CLR decreased in both cohorts, with FDR below 0.05 in both. Non-classical-monocyte CLR increased in both, but only Thomson passed FDR. Targeted Thomson batch estimates also remain directionally consistent, but both have exploratory all-label batch-family q=0.139044; no favorable targeted batch q is claimed.",
        "- **CD8-naive gene replication was incomplete:** all six targeted effects were positive in Zhu; CLIC1, S100A4 and CCL5 passed its six-test BH family. Four were positive in Thomson, but none passed its BH family. PRF1 and GZMA were negative in Thomson.",
        "- **P1 was cohort-dependent:** no Zhu strict-primary OLS test passed BH, whereas 11 Thomson labels passed. Zhu matched-null P1 survived for four labels; Thomson matched-null survived for 12. These are not pooled or interpreted as a universal effect.",
        "- **Hellinger did not show a universal positive displacement:** three Zhu cell types passed strict depth-adjusted BH, all with negative old-minus-young effects; no Thomson cell type passed. Exploratory Thomson batch/CMV families do have survivors, but they do not reverse the primary conclusion or establish robust aging biology. Non-significance is not equivalence.",
        "- **MAIT was untestable** in both released annotation universes. No integrated labels were imported to create it.",
        "",
        "All numerical effects, confidence intervals and q values are in the linked aggregate CSVs. Tables below include every strict-primary P1 and Hellinger cell type, every targeted effect, all targeted atlas/exact-age sensitivities, Thomson targeted covariate sensitivities, and all ineligible endpoints. The full exploratory pseudobulk tables are gzip-compressed CSV rather than truncated hit lists.",
        "",
        "## Input validation and label joins",
        "",
        "| Cohort | Donors downloaded | Raw cells | Annotation rows | Matched rows | Unmatched annotation rows excluded | Counts | Joins |",
        "|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for cohort in COHORTS:
        value = validation[cohort]
        lines.append(f"| {cohort} | {value['donors']} | {value['raw_cells']:,} | {value['annotation_rows']:,} | {value['annotation_rows_matched']:,} | {value['unmatched_annotation_rows_excluded']:,} | finite nonnegative integers | one-to-one |")
    lines += [
        "",
        "All 17 Zhu and 32 Thomson per-donor files passed full stored-value validation and byte-size checks. SHA-256 values and per-file dimensions are in `results/followup-pbmc-GSE213516-validation.json` and `results/followup-pbmc-GSE214546-validation.json`. Duplicate gene-symbol rows were excluded from P1 and targeted gene tests rather than artificially renamed.",
        "",
        "Zhu annotation aliases were resolved by exact full barcode matches: HR→F40, HYD→F65, XMZ→F60 and ZY→F46. F55 uses `baseline_<barcode>-F55`; the other direct donors use `<donor>_<barcode>`. Across Zhu, 199 annotation-only rows were explicitly excluded; every retained annotation maps uniquely to a raw barcode.",
        "",
        "**Zhu genome-build adjudication:** GEO says hg19 and the primary Methods say GRCh38. The actual released feature files contain Ensembl IDs, symbols and feature types, but no reference-genome identity; MatrixMarket metadata identifies Cell Ranger 3.0.0 without resolving the assembly. The conflict is therefore unresolved from public records. Zhu gene results are **symbol-level only**, with no coordinate, locus or assembly-specific claim.",
        "",
        "Released-label inventories and the explicit common-label map are in `results/followup-pbmc-label-harmonization.csv`. Only CD8-naive and CD16/non-classical monocytes support the targeted cross-cohort comparisons. MAIT is absent/non-harmonizable in both.",
        "",
        "## Strict targeted donor-pseudobulk replication",
        "",
        "Effects are old-minus-young differences in donor log2CPM. BH is over six genes separately within each cohort and model variant.",
        "",
        model_table(strict_pseudo, ["cohort", "cell_type", "gene"]),
        "",
        "### Thomson CMV/batch/pool/depth targeted sensitivities",
        "",
        "All 23 strict donors supported full-rank sex, CMV, B-batch, P-pool and depth models at the donor level (9 residual df for the 12-level B-batch model). Each sensitivity is a separate six-test BH family.",
        "",
        model_table(thomson_covariates, ["cell_type", "gene", "model"]),
        "",
        "## Strict targeted CLR composition replication",
        "",
        "Effects are old-minus-young donor CLR changes with a 0.5-cell pseudocount. Only eligible tests enter the three-target BH family; MAIT remains an explicit untestable row.",
        "",
        model_table(strict_comp, ["cohort", "cell_type", "common_cell_type"]),
        "",
        "### Thomson targeted CLR covariate robustness",
        "",
        "These estimates come directly from `results/followup-pbmc-GSE214546-clr-models.csv`. Depth means log total annotated cells. The q values below are the existing exploratory all-label families, not a recomputed two-target or otherwise favorable targeted family. In particular, batch gives CD8 Naive -1.396839, 95% CI [-2.466337, -0.327341], and CD16+ Monocytes +1.040266, 95% CI [0.111162, 1.969369]; both have exploratory all-label batch-family q=0.139044. CMV, pool and depth results are shown without substituting their values into the primary targeted table.",
        "",
        model_table(thomson_clr_covariates, ["cell_type", "model"]),
        "",
        "### Thomson strict P1 robustness counts",
        "",
        f"Primary has {p1_survivors['primary']}/13, depth {p1_survivors['depth']}/13, CMV {p1_survivors['cmv']}/13 and pool {p1_survivors['pool']}/13 cell types with q<0.05, but batch has {p1_survivors['batch']}/13 (minimum batch q={p1_batch_min_q:.6f}). Thus the broad P1 pattern is not robust to batch and is not confirmed aging biology.",
        "",
        "| Model | Cell types with q<0.05 (of 13) | Minimum q |",
        "|---|---:|---:|",
        f"| primary | {p1_survivors['primary']} | {thomson_p1_strict.loc[thomson_p1_strict['model'].eq('primary'), 'q_bh'].min():.6f} |",
        f"| depth | {p1_survivors['depth']} | {thomson_p1_strict.loc[thomson_p1_strict['model'].eq('depth'), 'q_bh'].min():.6f} |",
        f"| CMV | {p1_survivors['cmv']} | {thomson_p1_strict.loc[thomson_p1_strict['model'].eq('cmv'), 'q_bh'].min():.6f} |",
        f"| pool | {p1_survivors['pool']} | {thomson_p1_strict.loc[thomson_p1_strict['model'].eq('pool'), 'q_bh'].min():.6f} |",
        f"| batch | {p1_survivors['batch']} | {p1_batch_min_q:.6f} |",
        "",
        "## Complete strict-primary P1 results",
        "",
        "P1 uses per-cell log1p(CPM), arithmetic donor means and old-minus-young OLS with sex. `q` is the OLS BH value. The matched-null q values are in the aggregate P1 table and use 200 expression-decile-matched draws, seed 20261114.",
        "",
        model_table(strict_p1, ["cohort", "cell_type"]),
        "",
        "### Strict P1 matched-null results",
        "",
        model_table(strict_p1.assign(q_bh=strict_p1["random_q_bh"], pvalue=strict_p1["random_pvalue"]), ["cohort", "cell_type"]),
        "",
        "## Complete strict depth-adjusted Hellinger results",
        "",
        "Distances use actual raw UMI >=1500, technical-gene exclusions, strict-young cross-fitted references and donor mean per-cell distances. Negative effects mean old cells were closer to the young reference; they do not indicate a universal increase in displacement.",
        "",
        model_table(strict_hell, ["cohort", "cell_type"]),
        "",
        "### Thomson Hellinger covariate sensitivity caveat",
        "",
        f"The strict depth-adjusted primary family has no survivors (minimum q={thomson_hell_strict.loc[thomson_hell_strict['model'].eq('primary'), 'q_bh'].min():.6f}). The exploratory batch-depth family has {len(thomson_hell_batch_hits)} survivors and the exploratory CMV-depth family has {len(thomson_hell_cmv_hits)}; pool-depth has none. These sensitivity-family survivors have negative effects and do not reverse the primary result or establish robust aging biology.",
        "",
        model_table(pd.concat([thomson_hell_batch_hits, thomson_hell_cmv_hits], ignore_index=True), ["cell_type", "model"]),
        "",
        "## Targeted atlas-age and exact-age sensitivities",
        "",
        "Atlas groups add age-60 donors; exact-age models use all 17 Zhu and all 32 Thomson adults with continuous age. These are separate BH families, not replacements for the strict primary contrast.",
        "",
        model_table(sensitivity, ["cohort", "design", "cell_type", "gene"]),
        "",
        "The complete P1, exploratory CLR and Hellinger atlas/exact-age tables are in `results/followup-pbmc-replication-sensitivities.csv`. Thomson P1 and Hellinger CMV/batch/pool/depth results are in the cohort-specific model CSVs.",
        "",
        "## Exploratory donor-pseudobulk families",
        "",
        "All features detected in at least 20% of eligible donors were tested without effect-based prefiltering. BH was applied separately by cohort, design, model and cell type. Full effects/CIs/q values are in `results/followup-pbmc-GSE213516-pseudobulk-exploratory.csv.gz` and `results/followup-pbmc-GSE214546-pseudobulk-exploratory.csv.gz`; all q<0.05 rows are repeated in `results/followup-pbmc-pseudobulk-exploratory-survivors.csv.gz` for convenience.",
        "",
        "| Cohort | Design | Model | Tested genes | q<0.05 |",
        "|---|---|---|---:|---:|",
    ]
    for _, row in exploratory.iterrows():
        lines.append(f"| {row['cohort']} | {row['design']} | {row['model']} | {int(row['tested']):,} | {int(row['survivors']):,} |")
    lines += [
        "",
        "Exploratory CLR results, including every eligible released label and raw-fraction sensitivities, are in `results/followup-pbmc-clr-models.csv`.",
        "",
        "## Ineligible and untestable endpoints",
        "",
        "Eligibility failures are not null effects. Generic failures below include cell types failing the 20-cell/three-donor/group/six-donor/reference floors. Model-specific untestable rows follow separately.",
        "",
    ]
    if ineligible.empty:
        lines.append("No generic eligibility failures.")
    else:
        lines.append(model_table(ineligible.assign(status="untestable"), ["cohort", "endpoint", "cell_type", "design", "reason"]))
    lines += [
        "",
        "### Model-specific untestable rows",
        "",
        model_table(untestable_models, ["cohort", "endpoint", "cell_type", "design", "model", "reason"]),
        "",
        "## Covariate identifiability",
        "",
        "The complete rank table is `results/followup-pbmc-covariate-identifiability.csv`. At the donor level, every strict/atlas/exact-age Thomson primary, depth, CMV, B-batch and P-pool candidate was full rank. Cell-type-specific rank failures remain untestable, notably the strict Thomson platelet Hellinger model. CMV and pool adjustments are reported separately rather than hidden in a pooled model.",
        "",
        "## Resource use",
        "",
        "Every expression download and analysis/validation/summarization job used `scripts/resource_guard.py` with a 30 GB disk reserve, 3 GB per-file cap, 8 GB memory cap, swap disabled in the job cgroup and one heavy-job lock. No RAW tar was downloaded. Zhu completed before any Thomson download.",
        "",
        "| Cohort | Download files | Download bytes | Analysis seconds | Peak RSS GiB | Donor-cache bytes |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for cohort in COHORTS:
        row = resources[cohort]
        lines.append(f"| {cohort} | {row['download_files']} | {row['download_bytes']:,} | {row['elapsed_seconds']:.1f} | {row['max_rss_kib'] / 1024**2:.3f} | {row['donor_aggregate_cache_bytes']:,} |")
    lines += [
        "",
        "## Interpretation and limitations",
        "",
        "The strict composition result supports a reproducible CD8-naive decrease and supports, less definitively across studies, increased non-classical-monocyte abundance. Their Thomson batch estimates remain directional but have exploratory all-label batch q=0.139044, so batch robustness is not established. The six-gene targeted program did not replicate uniformly at FDR. Broad Thomson P1 is not robust to batch (0/13 survivors) and is not confirmed aging biology. Hellinger effects were predominantly negative and did not establish greater old-cell displacement; exploratory batch/CMV survivors do not change that primary conclusion. Atlas-age and exact-age sensitivities must be read alongside the strict results rather than used to select a preferred design.",
        "",
        "Zhu has only five donors per strict group and an unresolved assembly record; its symbol-level results should not be used for locus-specific claims. Thomson CMV/batch/pool models have limited residual degrees of freedom and observational covariate confounding. Released label definitions differ, and MAIT cannot be tested. No result establishes longitudinal change, causality or equivalence from non-significance.",
        "",
        "## Reproduction",
        "",
        "```bash",
        "python scripts/fetch_followup_pbmc.py --cohort GSE213516",
        "python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_pbmc.py --cohort GSE213516",
        "python scripts/fetch_followup_pbmc.py --cohort GSE214546",
        "python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_pbmc.py --cohort GSE214546",
        "python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/summarize_followup_pbmc.py",
        "python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/report_followup_pbmc.py",
        "```",
        "",
        "The fetch script invokes `resource_guard.py download` separately for every per-donor file. Raw data and donor caches remain outside Git under `../agescape-revision-data/pbmc/`.",
    ]
    (REPO / "reports" / "followup-pbmc-analysis.md").write_text("\n".join(lines) + "\n")
    print("Wrote reports/followup-pbmc-analysis.md", flush=True)


if __name__ == "__main__":
    main()
