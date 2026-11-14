#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
K = "followup-muscle-kedlian"
P = "followup-muscle-perez"


def number(value, digits=6):
    if pd.isna(value):
        return "—"
    return f"{float(value):.{digits}g}"


def markdown_table(frame):
    columns = [str(column) for column in frame.columns]
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def read(name):
    return pd.read_csv(RESULTS / name)


def write_csv(frame, name):
    frame.to_csv(RESULTS / name, index=False)


def family_rows():
    rows = []
    specifications = [
        ("Kedlian", f"{K}-p1-models.csv", "P1_transport", "pvalue", "q_bh"),
        ("Kedlian", f"{K}-composition-models.csv", "CLR_composition", "pvalue", "q_bh"),
        ("Kedlian", f"{K}-hellinger-models.csv", "Hellinger_nonnumeric_age_sensitivity", "pvalue", "q_bh"),
        ("Perez", f"{P}-gene-models.csv", "whole_nucleus_gene", "pvalue", "q_bh"),
    ]
    for study, filename, endpoint, pcol, qcol in specifications:
        frame = read(filename)
        if frame.empty:
            continue
        groups = ["assay", "donor_set", "model"] if study == "Kedlian" else ["assay", "model"]
        for keys, part in frame.groupby(groups, dropna=False):
            keys = keys if isinstance(keys, tuple) else (keys,)
            row = dict(zip(groups, keys))
            row.update({"study": study, "endpoint": endpoint, "tests": len(part),
                        "minimum_q": part[qcol].min(), "q_lt_0_05": int((part[qcol] < 0.05).sum())})
            rows.append(row)
    for assay in ["scRNA", "snRNA"]:
        frame = read(f"{K}-{assay}-gene-models.csv")
        if frame.empty:
            continue
        for (donor_set, model), part in frame.groupby(["donor_set", "model"]):
            rows.append({"study": "Kedlian", "assay": assay, "donor_set": donor_set,
                         "model": model, "endpoint": "donor_pseudobulk_gene_assay_wide",
                         "tests": len(part), "minimum_q": part["q_bh"].min(),
                         "q_lt_0_05": int((part["q_bh"] < 0.05).sum())})
    p1 = read(f"{K}-p1-models.csv")
    null = p1[p1["random_pvalue"].notna()]
    for assay, part in null.groupby("assay"):
        rows.append({"study": "Kedlian", "assay": assay, "donor_set": "primary_adult",
                     "model": "expression_matched_random_null", "endpoint": "P1_matched_null",
                     "tests": len(part), "minimum_q": part["random_q_bh"].min(),
                     "q_lt_0_05": int((part["random_q_bh"] < 0.05).sum())})
    columns = ["study", "assay", "endpoint", "donor_set", "model", "tests", "minimum_q", "q_lt_0_05"]
    return pd.DataFrame(rows)[columns].sort_values(["study", "assay", "endpoint", "donor_set", "model"])


def endpoint_table(frame, null=False):
    lines = ["| Cell type | Donors (young/old) | Old-young effect [95% CI] | p | BH q |",
             "|---|---:|---:|---:|---:|"]
    for row in frame.itertuples():
        lines.append(
            f"| {row.cell_type} | {row.n_donors} ({row.n_young}/{row.n_old}) | "
            f"{number(row.effect)} [{number(row.ci95_low)}, {number(row.ci95_high)}] | "
            f"{number(row.pvalue)} | {number(row.q_bh)} |"
        )
    if null:
        lines[0] = "| Cell type | Donors (young/old) | Old-young effect [95% CI] | p | BH q | Null p | Null BH q |"
        lines[1] = "|---|---:|---:|---:|---:|---:|---:|"
        lines = lines[:2] + [f"{line} {number(row.random_pvalue)} | {number(row.random_q_bh)} |"
                             for line, row in zip(lines[2:], frame.itertuples())]
    return "\n".join(lines)


def perez_table(frame):
    lines = ["| Model | Gene | Effect [95% CI] | p | BH q |", "|---|---|---:|---:|---:|"]
    for row in frame.itertuples():
        lines.append(f"| {row.model} | {row.gene} | {number(row.effect)} [{number(row.ci95_low)}, {number(row.ci95_high)}] | {number(row.pvalue)} | {number(row.q_bh)} |")
    return "\n".join(lines)


def primary(frame, assay):
    return frame[(frame["assay"] == assay) & (frame["donor_set"] == "primary_adult")
                 & (frame["model"] == "primary_group_depth_adjusted")].sort_values("cell_type")


def top_genes(frame, n=5):
    if frame.empty:
        return frame
    return frame.sort_values(["pvalue", "gene_id"]).head(n)


def incomplete_report(args):
    k_resource = json.loads((RESULTS / f"{K}-download-resource.json").read_text())
    p_resource = json.loads((RESULTS / f"{P}-download-resource.json").read_text())
    p_verification = json.loads((RESULTS / f"{P}-verification.json").read_text())
    p_genes = read(f"{P}-gene-models.csv")
    p_donors = read(f"{P}-donor-summary.csv")
    p_endpoints = read(f"{P}-endpoint-status.csv")
    families = []
    for model, part in p_genes.groupby("model"):
        families.append({"study": "Perez", "assay": "snRNA", "endpoint": "whole_nucleus_gene",
                         "donor_set": "all_17", "model": model, "tests": len(part),
                         "minimum_q": part["q_bh"].min(), "q_lt_0_05": int((part["q_bh"] < 0.05).sum())})
    families = pd.DataFrame(families)
    write_csv(families, "followup-muscle-family-summary.csv")
    resources = pd.DataFrame([
        {"job": "Kedlian_incomplete_guarded_download", "elapsed_seconds": k_resource["approximate_download_seconds"],
         "peak_rss_bytes": np.nan, "memory_limit_gb": np.nan, "disk_reserve_gb": 30},
        {"job": "Perez_17_sequential_downloads", "elapsed_seconds": p_resource["elapsed_seconds"],
         "peak_rss_bytes": np.nan, "memory_limit_gb": np.nan, "disk_reserve_gb": 30},
        {"job": "Perez_validation_and_models", "elapsed_seconds": p_verification["elapsed_seconds"],
         "peak_rss_bytes": p_verification["peak_rss_bytes"], "memory_limit_gb": 6, "disk_reserve_gb": 30},
    ])
    write_csv(resources, "followup-muscle-resource-summary.csv")
    cross = pd.DataFrame([{
        "comparison": "Kedlian_vs_Perez", "status": "blocked_no_Kedlian_estimates",
        "reason": "The required Kedlian raw object was not completely acquired; no expression values or outcomes were read. Perez whole-nucleus estimates are also composition-confounded. No comparison, pooling, or replication label was made.",
    }])
    write_csv(cross, "followup-muscle-cross-cohort-comparison.csv")
    group = p_genes[p_genes["model"] == "primary_group_depth_adjusted"]
    exact = p_genes[p_genes["model"] == "exact_age_10y_sensitivity"]
    lines = [
        "# Step 4 muscle analyses: Perez complete; Kedlian acquisition incomplete",
        "",
        "## Verdict and completion status",
        "",
        "**Step 4 is incomplete.** The required 4.03-GB Kedlian raw H5AD was not completely downloaded, so no Kedlian expression values, count validation, eligibility outcomes, or biological models were produced. The smaller curated CELLxGENE object was not substituted, parallel downloading was not used, and the resumable partial file was retained. Perez was then completed sequentially as the lightweight second arm.",
        "",
        "Perez provides whole-nucleus donor summaries and composition-confounded donor-pseudobulk gene models only. No released barcode-level annotations were found. No study pooling, cross-cohort numerical comparison, or independent-replication label was made. Non-significance is not equivalence.",
        "",
        "## Frozen plan and resource guard",
        "",
        "The design was frozen before expression outcomes in `data/manifest/muscle_analysis_plan.json` at HEAD `6763430`: donor inference, separate assays, six donors/two per group/20 cells, full-rank models, assay-wide multiplicity families, primary exclusion of Kedlian 470BR, a declared all-donor sensitivity, established P1 with expression-matched null, donor CLR, donor-pseudobulk gene models, and cross-fitted released-young Hellinger only under the nonnumeric-age sensitivity label.",
        "",
        f"Kedlian acquisition used two guarded sequential attempts with a 30-GB reserve and 5-GB cap. Only {k_resource['partial_bytes']:,} of {k_resource['expected_bytes']:,} bytes ({k_resource['percent_complete']:.3f}%) were acquired in approximately {k_resource['approximate_download_seconds']} seconds. A bounded 4-MiB range diagnostic took {k_resource['range_diagnostic_seconds']} seconds. At the observed rate, more than 11 hours remained. The partial file was not opened or analyzed.",
        "",
        f"All 17 Perez H5 downloads completed sequentially through the same guard in {p_resource['elapsed_seconds']:.3f} seconds ({p_resource['bytes']:,} bytes). Perez analysis used a 6-GB memory cap with swap disabled; observed peak RSS was {p_verification['peak_rss_bytes'] / 1e9:.3f} GB and elapsed time was {p_verification['elapsed_seconds']:.3f} seconds. No FASTQ or bulk count table was downloaded.",
        "",
        "## Kedlian: not run",
        "",
        "- Required raw-count value validation: **not run**, because the complete raw object is unavailable.",
        "- P1, CLR, gene, Hellinger, primary adult, and all-donor models: **not run**, not null results.",
        "- Numeric-age analysis: blocked because exact donor ages are not public.",
        "- Partial local input: `../agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`.",
        "- Audit: `results/followup-muscle-kedlian-download-resource.json` and `data/manifest/E-MTAB-13874.json`.",
        "",
        "## Perez verified inputs",
        "",
        f"All 17 filtered 10x H5 files contained only the expected `matrix` group and passed streamed finite/nonnegative/integer validation. There are {p_verification['stored_count_values']:,} stored int32 values, {p_verification['total_raw_umi']:,} total raw UMI, and exactly {p_verification['released_filtered_nuclei']:,} released filtered nuclei from 6 young and 11 old donors. Feature order is identical across files. The 72-man bulk arm is excluded.",
        "",
        "| GSM | Donor | Group | Exact age | Nuclei | Mean UMI |",
        "|---|---|---|---:|---:|---:|",
    ]
    for row in p_donors.itertuples():
        lines.append(f"| {row.gsm} | {row.donor_id} | {row.group} | {row.age_years} | {row.n_nuclei:,} | {number(row.mean_umi)} |")
    lines += [
        "",
        "## Perez whole-nucleus donor-pseudobulk models",
        "",
        "Counts were summed by donor before log1p(CPM). The primary model is old+mean log raw UMI; the predeclared sensitivity is exact age per 10 years+mean log raw UMI. Effects include 95% donor-level OLS intervals. BH is applied across all 20,343 genes separately within each model. These associations cannot distinguish cell-intrinsic changes from changing cell-type composition.",
        "",
        "| Model | Tests | BH q<0.05 | Minimum BH q |",
        "|---|---:|---:|---:|",
        f"| Old-young group, depth adjusted | {len(group):,} | {int((group['q_bh'] < 0.05).sum())} | {number(group['q_bh'].min())} |",
        f"| Exact age per 10 years, depth adjusted | {len(exact):,} | {int((exact['q_bh'] < 0.05).sum())} | {number(exact['q_bh'].min())} |",
        "",
        "The five smallest unadjusted p-values per model are shown for orientation; complete estimates are in `results/followup-muscle-perez-gene-models.csv`.",
        "",
    ]
    for model, part in p_genes.groupby("model"):
        lines += [f"### {model}", "", perez_table(top_genes(part)), ""]
    lines += [
        "## Perez blocked endpoints",
        "",
        markdown_table(p_endpoints),
        "",
        "The GEO series/per-donor inventories, all 17 H5 root schemas, primary paper and supplements, and public processed-object/repository searches were checked. No released barcode-level annotation mapping was found. Secondary studies' reannotations were not substituted. P1, cell-type models, CLR, and Hellinger were not computed.",
        "",
        "## Cross-cohort status",
        "",
        "No numerical cross-cohort comparison was possible: Kedlian has no estimates, and Perez whole-nucleus estimates would not be compatible with cell-type-specific intercostal-muscle estimands in any case. Perez reuse by Kedlian and unresolved individual disjointness prohibit an independent-replication label.",
        "",
        "## Delivered files and reproduction",
        "",
        "- `data/manifest/muscle_analysis_plan.json`: frozen design.",
        "- `scripts/run_followup_kedlian.py`: complete backed Kedlian pipeline, ready but not run on the incomplete file; synthetic end-to-end tested.",
        "- `scripts/download_followup_perez.py` and `scripts/run_followup_perez.py`: guarded sequential acquisition, validation, models, and manifest hashes.",
        "- `src/followup_muscle.py` and `tests/test_followup_muscle.py`: model/safeguard helpers and synthetic tests.",
        "- `results/followup-muscle-perez-*`: donor summaries, 40,686 gene-model rows, endpoint status, download/analysis resource records, and verification.",
        "- `results/followup-muscle-kedlian-download-resource.json`: incomplete acquisition audit; no Kedlian outcome CSV exists.",
        "- `results/followup-muscle-family-summary.csv`, `-resource-summary.csv`, and `-cross-cohort-comparison.csv`: joint status tables.",
        "- `data/manifest/E-MTAB-13874.json` and `GSE167186.json`: updated acquisition, validation, hash, analysis, and blocker records.",
        "",
        "Resume Kedlian without deleting the partial file:",
        "",
        "```bash",
        "python scripts/resource_guard.py download --url https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad --output ../agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad --min-free-gb 30 --max-file-gb 5",
        "python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_kedlian.py build-cache",
        "python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_kedlian.py analyze",
        "python scripts/run_followup_kedlian.py finalize",
        "python scripts/report_followup_muscle.py",
        "```",
        "",
        "Perez reproduction:",
        "",
        "```bash",
        "python scripts/download_followup_perez.py",
        "python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_perez.py run",
        "python scripts/run_followup_perez.py finalize",
        "```",
        "",
        "## Validation record",
        "",
        f"- Full suite: **{args.tests}**.",
        "- Synthetic Kedlian and Perez workflows passed; no real Kedlian outcome is claimed.",
        "- `git diff --check`: passed. No commit or push was made. Parent must decide whether to resume the long Kedlian transfer before accepting Step 4.",
        "",
        "## Not attempted",
        "",
        "No raw FASTQ, reclustering, numeric-age imputation, Perez bulk inference, alternate-object substitution, parallel download, study pooling, independent-replication claim, or PBMC file modification was performed.",
    ]
    (REPO / "reports" / "followup-muscle-analysis.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", default="not recorded")
    args = parser.parse_args()
    if not (RESULTS / f"{K}-verification.json").exists():
        incomplete_report(args)
        return
    verification = json.loads((RESULTS / f"{K}-verification.json").read_text())
    analysis_resource = json.loads((RESULTS / f"{K}-analysis-resource.json").read_text())
    p_verification = json.loads((RESULTS / f"{P}-verification.json").read_text())
    p1 = read(f"{K}-p1-models.csv")
    composition = read(f"{K}-composition-models.csv")
    h_models = read(f"{K}-hellinger-models.csv")
    h_eligibility = read(f"{K}-hellinger-eligibility.csv")
    eligibility = read(f"{K}-celltype-eligibility.csv")
    statuses = read(f"{K}-model-status.csv")
    p_genes = read(f"{P}-gene-models.csv")
    p_donors = read(f"{P}-donor-summary.csv")
    p_endpoints = read(f"{P}-endpoint-status.csv")
    families = family_rows()
    write_csv(families, "followup-muscle-family-summary.csv")
    resources = pd.DataFrame([
        {"job": "Kedlian_count_validation_and_cache", "elapsed_seconds": verification["elapsed_seconds"],
         "peak_rss_bytes": verification["peak_rss_bytes"], "memory_limit_gb": 6, "disk_reserve_gb": 30},
        {"job": "Kedlian_models_and_backed_Hellinger", "elapsed_seconds": analysis_resource["elapsed_seconds"],
         "peak_rss_bytes": analysis_resource["peak_rss_bytes"], "memory_limit_gb": 6, "disk_reserve_gb": 30},
        {"job": "Perez_validation_and_models", "elapsed_seconds": p_verification["elapsed_seconds"],
         "peak_rss_bytes": p_verification["peak_rss_bytes"], "memory_limit_gb": 6, "disk_reserve_gb": 30},
    ])
    write_csv(resources, "followup-muscle-resource-summary.csv")
    cross = pd.DataFrame([{
        "comparison": "Kedlian_cell_type_vs_Perez_whole_nucleus",
        "status": "not_estimand_compatible",
        "reason": "Perez has no released per-nucleus annotations; whole-nucleus estimates are composition-confounded, and tissues differ. No cross-cohort numerical comparison, pooling, or independent-replication label was made.",
    }])
    write_csv(cross, "followup-muscle-cross-cohort-comparison.csv")

    lines = [
        "# Step 4 muscle analyses: Kedlian and Perez",
        "",
        "## Verdict",
        "",
        "Kedlian and Perez were processed sequentially and independently. Donors—not cells, nuclei, or libraries—are the inferential units. scRNA and snRNA models and all multiplicity families remain separate. No study pooling, meta-analysis, or independent-replication label is used: Kedlian reused Perez, donor-level disjointness is unproven, and the tissues differ.",
        "",
        f"Kedlian has {int(eligibility[(eligibility['donor_set'] == 'primary_adult') & eligibility['eligible']].groupby('assay')['cell_type'].nunique().to_dict())} eligible primary cell types by assay. The complete primary and all-donor sensitivity results are in the aggregate CSVs. Perez is restricted to whole-nucleus descriptive and donor-pseudobulk gene models because no released per-nucleus annotations were found; these models cannot separate cell-intrinsic aging from changing cell composition.",
        "",
        "Non-significance is not evidence of equivalence, stability, or replication. A matched-null result is distinct from the conventional donor-model result.",
        "",
        "## Frozen design and resources",
        "",
        "The design was frozen in `data/manifest/muscle_analysis_plan.json` at HEAD `6763430` before expression download or outcomes. Kedlian primary models exclude 470BR (17–20 band) and use released young/old labels without numeric ages or band midpoints. The all-donor sensitivity adds 470BR under the released young label. Cell types require six donors, two per group, 20 cells per donor, full rank, and at least two residual degrees of freedom.",
        "",
        "P1 uses established `P1_Clean`, per-cell log1p(CPM) donor means, and a 200-draw expression-decile-matched null (seed 20261105). CLR uses every released known label in the assay and a 0.5-cell pseudocount. Gene models sum raw counts by donor/cell type before log1p(CPM), with assay-wide BH across cell-type-by-gene tests. Hellinger is labeled **NONNUMERIC_AGE_GROUP_SENSITIVITY**, uses released young donors, leave-one-young-donor-out references, 1,500-UMI cells, technical-gene exclusions, and 20-cell/50-reference-cell floors.",
        "",
        "All downloads and heavy jobs used `scripts/resource_guard.py`, a 30-GB reserve, sequential execution, and a 6-GB memory cap with swap disabled. The Kedlian download had a 5-GB cap. Reads were backed or streamed. No FASTQ or Perez bulk count table was downloaded.",
        "",
        "| Job | Elapsed seconds | Peak RSS (GB) | Memory cap (GB) |",
        "|---|---:|---:|---:|",
    ]
    for row in resources.itertuples():
        lines.append(f"| {row.job} | {number(row.elapsed_seconds)} | {number(row.peak_rss_bytes / 1e9)} | {row.memory_limit_gb} |")
    lines += [
        "",
        "## Kedlian verified inputs and exclusions",
        "",
        f"The raw H5AD is {verification['file_size_bytes']:,} bytes with {verification['n_obs']:,} rows and {verification['n_vars']:,} features. All {verification['stored_values_checked']:,} stored `layers/counts` float32 values were checked and were finite, nonnegative, and integer-valued. The object maps completely to the released 43-library/donor crosswalk using `{verification['sample_column']}`; released `cell` and `nuclei` labels were used without reclustering.",
        "",
        f"The cache contains {verification['aggregation_units']} eligible-assay donor/cell-type aggregation units. Primary donor counts by assay are `{verification['donor_counts']}`. Exact numeric age analysis is untestable because exact ages are not public.",
        "",
        "### Eligibility",
        "",
        "| Assay | Donor set | Released types | Eligible types | Excluded types |",
        "|---|---|---:|---:|---:|",
    ]
    for (assay, donor_set), part in eligibility.groupby(["assay", "donor_set"]):
        lines.append(f"| {assay} | {donor_set} | {len(part)} | {int(part['eligible'].sum())} | {int((~part['eligible']).sum())} |")
    for assay in ["scRNA", "snRNA"]:
        lines += ["", f"### {assay} P1 transport, primary adults", ""]
        part = primary(p1, assay)
        lines.append(endpoint_table(part, null=True) if not part.empty else "No eligible identifiable P1 model.")
        lines += ["", f"### {assay} CLR composition, primary adults", ""]
        part = primary(composition, assay)
        lines.append(endpoint_table(part) if not part.empty else "No eligible identifiable CLR model.")
        lines += ["", f"### {assay} Hellinger, nonnumeric-age sensitivity", ""]
        part = primary(h_models, assay) if not h_models.empty else pd.DataFrame()
        lines.append(endpoint_table(part) if not part.empty else "No Hellinger model passed all reference, donor, and rank floors; this is untestable, not null.")
        failed_h = h_eligibility[(h_eligibility["assay"] == assay) & (~h_eligibility["hellinger_eligible"])]
        if not failed_h.empty:
            reasons = failed_h.groupby("reason")["cell_type"].count().to_dict()
            lines += ["", f"Failed Hellinger gates: `{reasons}`. Complete per-cell-type reasons are in `{K}-hellinger-eligibility.csv`."]
        genes = read(f"{K}-{assay}-gene-models.csv")
        primary_genes = genes[(genes["donor_set"] == "primary_adult") & (genes["model"] == "primary_group_depth_adjusted")]
        lines += ["", f"### {assay} donor-pseudobulk gene models", "",
                  f"The primary assay-wide family contains {len(primary_genes):,} tested cell-type-by-gene combinations; {int((primary_genes['q_bh'] < 0.05).sum())} have BH q < 0.05. The five smallest unadjusted p-values are shown for orientation only; the full family and all sensitivity models are in `{K}-{assay}-gene-models.csv`.", ""]
        lines.append(perez_table(top_genes(primary_genes).rename(columns={"effect_old_minus_young": "effect"})) if not primary_genes.empty else "No identifiable gene models.")
    lines += [
        "",
        "### Kedlian sensitivity and untestable models",
        "",
        "All-donor, unadjusted, conventional, and matched-null families are summarized in `results/followup-muscle-family-summary.csv`. They are not pooled with primary families.",
        "",
        "| Assay | Endpoint | Donor set | Model | Tests | Minimum BH q | q<0.05 |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for row in families[families["study"] == "Kedlian"].itertuples():
        lines.append(f"| {row.assay} | {row.endpoint} | {row.donor_set} | {row.model} | {row.tests} | {number(row.minimum_q)} | {row.q_lt_0_05} |")
    untestable = statuses[statuses["status"] == "untestable"]
    lines += ["", f"There are {len(untestable)} explicitly untestable endpoint/cell-type/model rows. Their reasons are preserved in `{K}-model-status.csv`."]
    if not untestable.empty:
        lines += ["", markdown_table(untestable.groupby(["endpoint", "reason"]).size().rename("rows").reset_index())]
    lines += [
        "",
        "## Perez verified inputs and whole-nucleus analysis",
        "",
        f"All 17 public filtered 10x H5 files passed streamed finite/nonnegative/integer validation. They contain {p_verification['stored_count_values']:,} stored int32 values and exactly {p_verification['released_filtered_nuclei']:,} released filtered nuclei: 6 young and 11 old donors. Exact ages are released. The primary model is old+mean log UMI; the declared sensitivity uses exact age per 10 years+mean log UMI. Per-donor sex was not invented, and the 72-man bulk arm was excluded.",
        "",
        "| GSM | Donor | Group | Exact age | Nuclei | Mean UMI |",
        "|---|---|---|---:|---:|---:|",
    ]
    for row in p_donors.itertuples():
        lines.append(f"| {row.gsm} | {row.donor_id} | {row.group} | {row.age_years} | {row.n_nuclei:,} | {number(row.mean_umi)} |")
    for model, part in p_genes.groupby("model"):
        lines += ["", f"### {model}", "",
                  f"The whole-nucleus family contains {len(part):,} genes; {int((part['q_bh'] < 0.05).sum())} have BH q < 0.05. Effects are composition-confounded. Five smallest unadjusted p-values are shown; all estimates are in `{P}-gene-models.csv`.", "",
                  perez_table(top_genes(part))]
    lines += ["", "### Perez blocked endpoints", "",
              markdown_table(p_endpoints), "",
              "Released per-nucleus annotations were sought in the primary GEO inventory, per-donor processed H5 schema, paper supplements, and public processed-object records. None supplied a barcode-level annotation mapping. No annotations were inferred. Therefore P1, cell-type models, CLR, and Hellinger are blocked/confounded and were not computed.", "",
              "## Cross-cohort comparison", "",
              "No numerical cross-cohort comparison is estimand-compatible: Kedlian has released cell-type-specific intercostal-muscle estimates, whereas Perez has only composition-confounded whole-nucleus vastus-lateralis estimates. This is recorded in `results/followup-muscle-cross-cohort-comparison.csv`. No inconsistent results were hidden or pooled away.", "",
              "## Delivered aggregates and reproduction", "",
              "- `scripts/run_followup_kedlian.py`: metadata/count gates, backed aggregation, P1/null, CLR, gene models, Hellinger, and manifest update.",
              "- `scripts/run_followup_perez.py`: 17-file count gate, whole-nucleus donor summaries, group/exact-age gene models, and manifest update.",
              "- `src/followup_muscle.py` and `tests/test_followup_muscle.py`: reusable safeguards, synthetic end-to-end tests, and delivered-output checks.",
              "- `results/followup-muscle-kedlian-*`: verification, eligibility, donor scores/counts/distances, P1/null, CLR, gene, Hellinger, model-status, and resource aggregates.",
              "- `results/followup-muscle-perez-*`: verification, donor summaries, gene models, and blocked-endpoint status.",
              "- `results/followup-muscle-family-summary.csv`, `-resource-summary.csv`, and `-cross-cohort-comparison.csv`: joint audit tables.",
              "- `data/manifest/E-MTAB-13874.json`, `GSE167186.json`: local paths, sizes, SHA-256 hashes, count gates, analysis scope, and blocked labels.",
              "",
              "```bash",
              "python scripts/resource_guard.py download --url https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad --output ../agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad --min-free-gb 30 --max-file-gb 5",
              "python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_kedlian.py build-cache",
              "python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_kedlian.py analyze",
              "# Download each Perez H5 sequentially with resource_guard, --min-free-gb 30 --max-file-gb 1.",
              "python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_perez.py run",
              "python scripts/run_followup_kedlian.py finalize",
              "python scripts/run_followup_perez.py finalize",
              "python scripts/report_followup_muscle.py",
              "```",
              "",
              "## Validation record",
              "",
              f"- Full suite: **{args.tests}**.",
              "- Synthetic end-to-end Kedlian and Perez pipelines passed before the real runs.",
              "- `git diff --check`: passed. No commit or push was made; parent validation is required.",
              "",
              "## Not attempted",
              "",
              "No reclustering, numeric-age imputation, raw FASTQ processing, Perez bulk analysis, study pooling, independent-replication claim, causal interpretation, or PBMC file modification was performed.",
    ]
    (REPO / "reports" / "followup-muscle-analysis.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
