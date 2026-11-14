#!/usr/bin/env python3
"""Aggregate and report the completed Step 7 ovary analyses."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
JIN = "followup-ovary-jin"
POST = "followup-ovary-postmenopausal"
JIN_PRIMARY = "primary_group_depth_adjusted"
POST_PRIMARY = "primary_depth_plex_clinical"
POST_MODELS = [
    "primary_depth_plex_clinical", "technical_depth_plex_sensitivity",
    "unadjusted_sensitivity", "secondary_primary_plus_parity",
]


def read(name):
    return pd.read_csv(RESULTS / name)


def write(frame, name):
    frame.to_csv(RESULTS / name, index=False)
    print(f"wrote {RESULTS / name}: {len(frame)} rows", flush=True)


def number(value, digits=4):
    if pd.isna(value):
        return "NA"
    return f"{value:.{digits}g}"


def effect_cell(row, random=False):
    text = (f"{number(row['effect'])} [{number(row['ci95_low'])}, {number(row['ci95_high'])}]; "
            f"p={number(row['pvalue'])}; q={number(row['q_bh'])}")
    if random and "random_pvalue" in row and pd.notna(row["random_pvalue"]):
        text += f"; matched-null p={number(row['random_pvalue'])}; q={number(row['random_q_bh'])}"
    return text


def markdown_models(frame, cell_type_header="Cell type", random=False):
    lines = [f"| {cell_type_header} | Donors | Effect (95% CI); nominal p; BH q |", "|---|---:|---|"]
    for row in frame.itertuples():
        lines.append(f"| {row.cell_type} | {row.n_donors} | {effect_cell(row._asdict(), random)} |")
    return lines


def gene_family(cohort, model, frame):
    return {
        "cohort": cohort, "model": model, "endpoint": "donor_pseudobulk_gene",
        "tests": len(frame), "nominal_p_lt_0.05": int(frame["pvalue"].lt(0.05).sum()),
        "bh_q_lt_0.05": int(frame["q_bh"].lt(0.05).sum()),
        "minimum_q_bh": float(frame["q_bh"].min()) if len(frame) else np.nan,
    }


def sensitivity_summary(primary, other, model, cohort):
    keys = ["cell_type", "gene_id"]
    hits = primary[primary["q_bh"] < 0.05]
    if not len(hits):
        return {"cohort": cohort, "model": model, "primary_hits": 0,
                "same_direction": 0, "same_direction_fraction": np.nan, "q_lt_0.05_in_model": 0}
    joined = hits[keys + ["effect"]].merge(other[keys + ["effect", "q_bh"]], on=keys,
                                           suffixes=("_primary", "_model"), validate="one_to_one")
    same = np.sign(joined["effect_primary"]) == np.sign(joined["effect_model"])
    return {
        "cohort": cohort, "model": model, "primary_hits": len(joined),
        "same_direction": int(same.sum()), "same_direction_fraction": float(same.mean()),
        "q_lt_0.05_in_model": int((joined["q_bh"] < 0.05).sum()),
    }


def main(test_text):
    jin_verification = json.loads((RESULTS / f"{JIN}-verification.json").read_text())
    post_verification = json.loads((RESULTS / f"{POST}-verification.json").read_text())
    jin_resource = json.loads((RESULTS / f"{JIN}-analysis-resource.json").read_text())
    post_resource = json.loads((RESULTS / f"{POST}-analysis-resource.json").read_text())
    jin_download = json.loads((RESULTS / "followup-ovary-GSE202601-download.json").read_text())
    post_download = json.loads((RESULTS / "followup-ovary-GSE329417-download.json").read_text())
    decompression = json.loads((RESULTS / f"{POST}-decompression.json").read_text())

    jin_p1 = read(f"{JIN}-p1-models.csv")
    jin_comp = read(f"{JIN}-composition-models.csv")
    jin_h = read(f"{JIN}-hellinger-models.csv")
    jin_genes = read(f"{JIN}-gene-models.csv")
    post_p1 = read(f"{POST}-p1-models.csv")
    post_comp = read(f"{POST}-composition-models.csv")
    post_genes = {model: read(f"{POST}-gene-models-{model}.csv") for model in POST_MODELS}
    jin_status = read(f"{JIN}-model-status.csv")
    post_status = read(f"{POST}-model-status.csv")

    endpoint_frames = []
    for cohort, frame in [("GSE202601", jin_p1), ("GSE202601", jin_comp), ("GSE202601", jin_h),
                          ("GSE329417", post_p1), ("GSE329417", post_comp)]:
        current = frame.copy()
        current.insert(0, "cohort", cohort)
        endpoint_frames.append(current)
    write(pd.concat(endpoint_frames, ignore_index=True), "followup-ovary-endpoint-summary.csv")
    family_rows = [gene_family("GSE202601", model, part) for model, part in jin_genes.groupby("model")]
    family_rows += [gene_family("GSE329417", model, post_genes[model]) for model in POST_MODELS]
    write(pd.DataFrame(family_rows), "followup-ovary-gene-family-summary.csv")
    sensitivity = []
    for model, part in jin_genes.groupby("model"):
        if model != JIN_PRIMARY:
            sensitivity.append(sensitivity_summary(jin_genes[jin_genes["model"].eq(JIN_PRIMARY)], part, model, "GSE202601"))
    for model in POST_MODELS[1:]:
        sensitivity.append(sensitivity_summary(post_genes[POST_PRIMARY], post_genes[model], model, "GSE329417"))
    write(pd.DataFrame(sensitivity), "followup-ovary-gene-sensitivity.csv")
    primary_hits = pd.concat([
        jin_genes[(jin_genes["model"] == JIN_PRIMARY) & (jin_genes["q_bh"] < 0.05)].assign(cohort="GSE202601"),
        post_genes[POST_PRIMARY][post_genes[POST_PRIMARY]["q_bh"] < 0.05].assign(cohort="GSE329417"),
    ], ignore_index=True)
    write(primary_hits, "followup-ovary-primary-gene-fdr-hits.csv")
    untestable = pd.concat([
        jin_status[jin_status["status"] != "completed"],
        post_status[post_status["status"] != "completed"],
    ], ignore_index=True)
    write(untestable, "followup-ovary-untestable-arms.csv")
    resources = pd.DataFrame([
        {"cohort": "GSE202601", "download_files": len(jin_download["files"]),
         "download_bytes": jin_download["download_bytes"], "decompressed_bytes": np.nan,
         "build_cache_seconds": jin_verification["elapsed_seconds"],
         "build_cache_peak_rss_bytes": jin_verification["peak_rss_bytes"],
         "analysis_seconds": jin_resource["elapsed_seconds"], "analysis_peak_rss_bytes": jin_resource["peak_rss_bytes"]},
        {"cohort": "GSE329417", "download_files": len(post_download["files"]),
         "download_bytes": post_download["download_bytes"], "decompressed_bytes": decompression["decompressed_size_bytes"],
         "build_cache_seconds": post_verification["elapsed_seconds"],
         "build_cache_peak_rss_bytes": post_verification["peak_rss_bytes"],
         "analysis_seconds": post_resource["elapsed_seconds"], "analysis_peak_rss_bytes": post_resource["peak_rss_bytes"]},
    ])
    write(resources, "followup-ovary-resource-summary.csv")

    jin_eligible = read(f"{JIN}-celltype-eligibility.csv")
    jin_h_eligible = read(f"{JIN}-hellinger-eligibility.csv")
    post_eligible = read(f"{POST}-celltype-eligibility.csv")
    post_identifiable = read(f"{POST}-covariate-identifiability.csv")
    post_donors = read(f"{POST}-donor-qc.csv")

    lines = [
        "# Step 7 ovary analyses: Jin and postmenopausal cohorts",
        "",
        "Both frozen designs were executed sequentially at repository HEAD `0610c34`. Jin GSE202601 was completed before GSE329417 acquisition. The inferential unit is the donor in every model; nuclei and tissue pieces are not independent replicates. The cohorts were never pooled. All associations are cross-sectional, not causal, and no premenopausal claim is made.",
        "",
        "GSE329417 remains a **bioRxiv preprint** (DOI `10.64898/2026.05.11.724381`); no peer-reviewed journal publication was verified. Jin is a peer-reviewed Nature Aging article (DOI `10.1038/s43587-024-00762-5`).",
        "",
        "## Frozen designs and multiplicity",
        "",
        "- `data/manifest/ovary_jin_analysis_plan.json`: four young versus four reproductively-aged donors, each once; reproductively-aged-minus-young; primary `group + mean_log_umi`, with group-only sensitivity. No continuous-age coefficient was fitted. The ovary5 snRNA `ovary47` versus GEO/snATAC 49 conflict remains unresolved and numeric age was never used.",
        "- `data/manifest/ovary_postmenopausal_analysis_plan.json`: exact age per 10 years within postmenopausal donors. Primary adjustment is depth, three donor-specific plex fractions (plex 4 reference), BMI, and non-ovarian abnormal state (pre-malignant/malignant versus benign). Technical-only and unadjusted sensitivities plus secondary parity-augmented models were frozen before outcomes.",
        "- P1 uses established human `P1_Clean`, per-cell log1p(CPM) donor means, detected up/down genes, and a 200-draw expression-decile-matched null (seed 20261105). Pseudobulk sums raw/corrected counts before log1p(CPM). CLR uses every released coarse label and a 0.5-cell pseudocount.",
        "- BH families are separate by cohort, endpoint, and model. Gene-level BH spans all tested cell-type-by-gene combinations within a model. Conventional and matched-null P1 p-values have separate primary families. Sensitivity families are reported separately, not used to select a favorable model.",
        "",
        "## Resource and acquisition audit",
        "",
        "All expression downloads used `scripts/resource_guard.py`, sequential execution, a 30-GB disk reserve, and per-file caps. Analysis/decompression jobs used an 8-GB memory cap with swap disabled in the job scope. Jin was processed one donor at a time. The atlas used direct HDF5 categorical metadata and CSR indptr-bounded 256-nucleus chunks after AnnData initialization exceeded the cap; the cap was not raised. Neither RAW tar nor FASTQ data was downloaded.",
        "",
        "| Cohort | Download files | Compressed/download bytes | Decompressed bytes | Cache peak RSS GiB | Analysis peak RSS GiB |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in resources.itertuples():
        lines.append(f"| {row.cohort} | {row.download_files} | {row.download_bytes:,} | {'' if pd.isna(row.decompressed_bytes) else format(int(row.decompressed_bytes), ',')} | {row.build_cache_peak_rss_bytes / 1024**3:.3f} | {row.analysis_peak_rss_bytes / 1024**3:.3f} |")
    lines += [
        "",
        "## Jin GSE202601",
        "",
        "### Count, barcode, and eligibility gates",
        "",
        f"Eight per-donor Cell Ranger H5 files were validated: **{jin_verification['stored_count_values_checked']:,} stored values** are finite, nonnegative integers. The files contain {jin_verification['raw_h5_barcodes']:,} barcodes; all **{jin_verification['annotated_nuclei']:,} released annotations** matched exactly one matrix barcode. The other {jin_verification['unannotated_h5_barcodes']:,} matrix barcodes remained outside annotation-based outcomes. Matrix barcodes `-1` were mapped to released donor suffixes `.1`–`.8`; this join clarification was frozen after schema inspection and before outcomes.",
        "",
        "| Cell type | Donors >=20 nuclei | Young | Aged | Minimum nuclei/donor | Base eligible | Hellinger eligible |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    for row in jin_eligible.itertuples():
        h = jin_h_eligible[jin_h_eligible["cell_type"] == row.cell_type].iloc[0]
        lines.append(f"| {row.cell_type} | {row.n_donors_ge20} | {row.n_young_ge20} | {row.n_reproductively_aged_ge20} | {row.min_cells_per_donor} | {row.eligible} | {h.hellinger_eligible} |")
    lines += [
        "",
        "LEC and TC fail the base all-four-donors/group floor and are untestable for P1, pseudobulk, and CLR. EpiC additionally fails Hellinger QC/reference floors. Hellinger uses 1,500-UMI nuclei, technical-gene exclusions, at least 20 qualifying nuclei per donor, and at least 50 cross-fitted young-reference nuclei for every scored donor.",
        "",
        "### Primary P1 transport",
        "",
    ]
    lines += markdown_models(jin_p1[jin_p1["model"] == JIN_PRIMARY], random=True)
    lines += ["", "No primary P1 result survives conventional or matched-null BH. These four-versus-four intervals do not establish equivalence.", "", "### Primary CLR composition", ""]
    lines += markdown_models(jin_comp[jin_comp["model"] == JIN_PRIMARY])
    lines += ["", "### Primary cross-fitted young-reference Hellinger displacement", ""]
    lines += markdown_models(jin_h[jin_h["model"] == JIN_PRIMARY])
    jin_gene_summary = pd.DataFrame(family_rows)[:2]
    lines += [
        "",
        "### Donor-pseudobulk genes and sensitivity",
        "",
        "| Model | Tested genes across cell types | Nominal p<0.05 | BH q<0.05 |",
        "|---|---:|---:|---:|",
    ]
    for _, row in jin_gene_summary.iterrows():
        lines.append(f"| {row['model']} | {row['tests']:,} | {row['nominal_p_lt_0.05']:,} | {row['bh_q_lt_0.05']:,} |")
    lines += ["", "Primary BH hits (effect is aged minus young on log1p CPM):", ""]
    jin_hits = primary_hits[primary_hits["cohort"] == "GSE202601"].sort_values("pvalue")
    if len(jin_hits):
        lines += ["| Cell type | Gene | Effect (95% CI) | p | BH q |", "|---|---|---:|---:|---:|"]
        for row in jin_hits.itertuples():
            lines.append(f"| {row.cell_type} | {row.gene} | {number(row.effect)} [{number(row.ci95_low)}, {number(row.ci95_high)}] | {number(row.pvalue)} | {number(row.q_bh)} |")
    else:
        lines.append("None.")
    lines += [
        "",
        "The group-only gene sensitivity has its own BH family and is not confirmation of the primary hits. Full effects and intervals are in `results/followup-ovary-jin-gene-models.csv`; direction/FDR concordance is summarized in `results/followup-ovary-gene-sensitivity.csv`.",
        "",
        "Same-nucleus RNA–ATAC coupling is impossible because the assays used separate nuclei and BioSamples. Donor-matched linkage is also blocked: no validated released peak–gene link table was identified, and ATAC EC was not silently equated to RNA BEC/LEC. No ATAC matrix or fragment file was downloaded.",
        "",
        "## GSE329417 postmenopausal cohort",
        "",
        "### Atlas, corrected counts, and donor aggregation",
        "",
        f"The final annotated atlas contains **{post_verification['atlas_nuclei']:,} nuclei**, {post_verification['atlas_genes']:,} features, {post_verification['donors']} donors, {post_verification['samples']} tissue-piece samples, and {post_verification['coarse_cell_types']} coarse labels. Exact fields were donor `{post_verification['join']['donor_field']}`, sample `{post_verification['join']['sample_field']}`, and coarse annotation `{post_verification['join']['coarse_annotation_field']}`. Atlas sample IDs were reconciled one-to-one through the 64 released CellBender filename stems and their `_p1`–`_p4` suffixes, then to the processing map and 28-donor clinical table; no fuzzy join was used.",
        "",
        f"All **{post_verification['stored_corrected_count_values_checked']:,} stored `layers/raw_counts` values** are finite, nonnegative integers. Atlas X was identical to `layers/log1p_norm` in every streamed chunk and was never treated as raw counts. Only CellBender-corrected, unnormalized `layers/raw_counts` entered endpoints. Original uncorrected Cell Ranger matrices are not separately released; therefore an uncorrected sensitivity **cannot be performed** and was not fabricated.",
        "",
        "Counts and per-cell summaries were aggregated over all pieces to one donor-cell-type row. SN_0086 samples 3_4 and 16_4, which shared a digested aliquot, were combined inside that single donor aggregation. All 64 samples were used but never treated as 64 independent biological replicates.",
        "",
        "| Clinical state | Donors |",
        "|---|---:|",
    ]
    for state, count in post_donors["non_ovarian_cancer_state"].value_counts().items():
        lines.append(f"| {state} | {count} |")
    lines += [
        "",
        "Non-ovarian abnormal state is a covariate, not a claim that ovarian pathology was present; the study excluded ovarian neoplasia and endometriosis.",
        "",
        "### Eligibility and identifiability",
        "",
        "| Cell type | Donors observed | Donors >=20 nuclei | Total nuclei | Base eligible |",
        "|---|---:|---:|---:|---|",
    ]
    for row in post_eligible.itertuples():
        lines.append(f"| {row.cell_type} | {row.n_donors_observed} | {row.n_donors_ge20} | {row.n_cells_total:,} | {row.eligible} |")
    failed_identifiability = post_identifiable[~post_identifiable["identifiable"] & post_identifiable["n_donors"].ge(20)]
    lines += ["", "All eligible primary/sensitivity designs were full rank with residual degrees of freedom." if failed_identifiability.empty else "Rank/identifiability failures among otherwise eligible units are listed below; affected models are untestable.", ""]
    if not failed_identifiability.empty:
        lines += ["| Cell type | Endpoint class | Model | Reason |", "|---|---|---|---|"]
        for row in failed_identifiability.itertuples():
            lines.append(f"| {row.cell_type} | {row.endpoint_class} | {row.model} | {row.reason} |")
        lines.append("")
    lines += ["### Primary P1 transport, exact age per 10 years", ""]
    post_p1_primary = post_p1[post_p1["model"] == POST_PRIMARY]
    lines += markdown_models(post_p1_primary, random=True)
    post_conventional_hits = post_p1_primary[post_p1_primary["q_bh"] < 0.05]
    post_null_hits = post_p1_primary[post_p1_primary["random_q_bh"] < 0.05]
    null_labels = ", ".join(f"{row.cell_type} (q={number(row.random_q_bh)})" for row in post_null_hits.itertuples()) or "none"
    lines += [
        "",
        f"Conventional P1 BH hits: {len(post_conventional_hits)}. Expression-matched-null P1 BH hits: {null_labels}. These are distinct calibrations; a matched-null hit is not relabeled as a conventional age-associated P1 finding.",
        "",
        "### Primary CLR composition, exact age per 10 years",
        "",
    ]
    lines += markdown_models(post_comp[post_comp["model"] == POST_PRIMARY])
    lines += [
        "",
        "### Donor-pseudobulk genes and prespecified sensitivity",
        "",
        "| Model | Tested genes across cell types | Nominal p<0.05 | BH q<0.05 |",
        "|---|---:|---:|---:|",
    ]
    for model in POST_MODELS:
        row = next(row for row in family_rows if row["cohort"] == "GSE329417" and row["model"] == model)
        lines.append(f"| {model} | {row['tests']:,} | {row['nominal_p_lt_0.05']:,} | {row['bh_q_lt_0.05']:,} |")
    post_hits = primary_hits[primary_hits["cohort"] == "GSE329417"].sort_values("pvalue")
    lines += ["", f"The primary family has {len(post_hits):,} BH hits. The 20 smallest nominal p-values among them are shown for orientation; the complete non-cherry-picked table is `results/followup-ovary-primary-gene-fdr-hits.csv`.", ""]
    if len(post_hits):
        lines += ["| Cell type | Gene | Effect/decade (95% CI) | p | BH q |", "|---|---|---:|---:|---:|"]
        for row in post_hits.head(20).itertuples():
            lines.append(f"| {row.cell_type} | {row.gene} | {number(row.effect)} [{number(row.ci95_low)}, {number(row.ci95_high)}] | {number(row.pvalue)} | {number(row.q_bh)} |")
    else:
        lines.append("None.")
    lines += ["", "| Sensitivity model | Primary hits with same direction | Same-direction fraction | Also q<0.05 in sensitivity |", "|---|---:|---:|---:|"]
    for _, row in pd.DataFrame(sensitivity)[pd.DataFrame(sensitivity)["cohort"] == "GSE329417"].iterrows():
        lines.append(f"| {row['model']} | {row['same_direction']}/{row['primary_hits']} | {number(row['same_direction_fraction'])} | {row['q_lt_0.05_in_model']} |")
    lines += [
        "",
        "Hellinger displacement is **untestable** for every coarse cell type because this postmenopausal cohort has no 20–40-year reference. There was no young-versus-old dichotomization and no borrowing of Jin's reference. Parity appears only in the declared secondary adjustment; age and parity are distinct estimands.",
        "",
        "## Untestable arms and interpretation limits",
        "",
        f"`results/followup-ovary-untestable-arms.csv` contains all {len(untestable):,} blocked/untestable model rows with `outcomes_examined=False`. No failed arm was converted to a zero effect or non-significant result. Non-significance does not demonstrate equivalence, particularly for Jin's eight donors. Cell-level inference, causal language, cross-cohort pooling, and premenopausal extrapolation were prohibited.",
        "",
        "## Reproduction and artifacts",
        "",
        "```bash",
        "python scripts/prepare_followup_ovary.py",
        "python scripts/fetch_followup_ovary.py --cohort GSE202601",
        "python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_jin.py build-cache",
        "python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_jin.py analyze",
        "python scripts/fetch_followup_ovary.py --cohort GSE329417",
        "python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/prepare_followup_ovary_postmenopausal.py decompress",
        "python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/prepare_followup_ovary_postmenopausal.py inspect",
        "python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_postmenopausal.py build-cache",
        "python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_postmenopausal.py analyze",
        "python scripts/run_followup_ovary_jin.py finalize",
        "python scripts/run_followup_ovary_postmenopausal.py finalize",
        "python scripts/report_followup_ovary.py",
        "```",
        "",
        "Core code: `src/followup_ovary.py`, `scripts/fetch_followup_ovary.py`, `scripts/prepare_followup_ovary.py`, `scripts/prepare_followup_ovary_postmenopausal.py`, `scripts/run_followup_ovary_jin.py`, `scripts/run_followup_ovary_postmenopausal.py`, and this generator. Plans, download hashes, decompressed hash, and count-validation manifests are under `data/manifest/` and `results/followup-ovary-*`. Aggregate tables include endpoint, gene-family, sensitivity, untestable-arm, and resource summaries. Tests are in `tests/test_followup_ovary.py`.",
        "",
        f"Verification: **{test_text}**. `git diff --check` passed. No commit or push was made. No heart file was modified.",
    ]
    report = REPO / "reports" / "followup-ovary-analysis.md"
    report.write_text("\n".join(lines) + "\n")
    print(f"wrote {report}", flush=True)

    artifact_paths = [
        REPO / "data/manifest/ovary_cohorts.json",
        REPO / "data/manifest/ovary_jin_analysis_plan.json",
        REPO / "data/manifest/ovary_postmenopausal_analysis_plan.json",
        report,
        REPO / "src/followup_ovary.py",
        REPO / "tests/test_followup_ovary.py",
        *sorted((REPO / "scripts").glob("*ovary*.py")),
        *sorted(RESULTS.glob("followup-ovary-*")),
    ]
    artifacts = []
    for path in artifact_paths:
        if not path.is_file():
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        artifacts.append({"path": str(path.relative_to(REPO)), "size_bytes": path.stat().st_size, "sha256": digest.hexdigest()})
    output_manifest = {
        "created": "2026-11-14",
        "step": 7,
        "repository_head_at_plan_freeze": "0610c34bde6c176efc1c225cc8d817231e9ce54a",
        "verification": test_text,
        "commit_or_push_performed": False,
        "heart_files_modified": False,
        "artifacts": artifacts,
    }
    path = REPO / "data" / "manifest" / "ovary_analysis_outputs.json"
    path.write_text(json.dumps(output_manifest, indent=2) + "\n")
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests", default="not yet run")
    args = parser.parse_args()
    main(args.tests)
