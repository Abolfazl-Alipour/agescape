# Jia 2026 lifespan heart atlas: access and metadata gate

Base revision: `0610c34bde6c176efc1c225cc8d817231e9ce54a`. Verification date: 2026-11-14.

## Scope and verdict

This was a **lightweight metadata and public-access check only**. No expression matrix, H5AD, loom, FASTQ, BAM, or count values were downloaded or read. No eligibility outcome, aging association, or biological result was computed. Small repository/API records, primary full text, supplementary metadata, and three short analysis scripts were read. Downloaded metadata are under `/home/reza/projects/agescape-revision-data/jia-heart-metadata/`; exact per-file sizes and source hashes are in `data/manifest/jia_heart.json`.

The publisher DOI page returned HTTP 403 in this environment. The same primary Science Advances article was therefore read as [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13274588/fullTextXML), including its actual data-availability statement. No accession was inferred from a secondary summary.

| Gate | Verdict | Reason |
|---|---|---|
| Human accession and public access | **READY** | The paper explicitly names Figshare DOI `10.6084/m9.figshare.31939227`; the actual v1 record is public, unembargoed, downloadable, and CC BY 4.0. |
| Complete 54-sample acquisition map | **READY WITH CONDITIONS** | Jia Figshare supplies only 44 AG1–AG5 matrices. The ten AG6 LV/RV matrices must come from the separately public Litviňuková/HCA release. |
| Adult-only exact continuous-age analysis, 14–75 years | **BLOCKED** | Only six internal postnatal snRNA donors have exact ages. The five AG6 snRNA donors have five-year bands. Six other internal donors with exact ages did **not** undergo snRNA-seq and cannot be substituted. |
| Adult group analysis | **BLOCKED as confirmatory evidence; descriptive sensitivity only after expression gates** | AG4/AG5 are all internal and AG6 is entirely external, so study/source is completely confounded with the old group. Only 6 versus 5 donors are available before cell-type exclusions. |
| Lifespan sensitivity analysis | **READY WITH CONDITIONS for stage-stratified/group sensitivity; BLOCKED as one exact continuous gestational-to-75-year slope** | Exact fetal ages and public matrices exist, but AG6 has bands and a different source. Count integrity, Jia barcode-level annotation, cross-release gene/barcode joins, and donor floors remain unchecked. Gestational weeks and postnatal years must not be silently placed on one numeric axis. |

“Ready with conditions” does not authorize expression acquisition here. Any future expression download or heavy job requires parent approval, a sequentially resource-guarded job, a 30-GB disk reserve, and at most 8 GB analysis memory.

## Primary accession and linked records

The Jia paper is Hao Jia et al., “Life-span–dependent transcriptional dynamics of the human heart,” *Science Advances* 12(25), eaeg2614 (2026), DOI `10.1126/sciadv.aeg2614`, PMID `42308288`, PMCID `PMC13274588`, first published 2026-06-17.

Its data-availability statement says: “Human heart snRNA-seq data have been deposited in the Figshare under the accession DOI: **10.6084/m9.figshare.31939227**.” It separately names mouse snRNA-seq `31939839`, hiPSC-CM bulk RNA-seq `31939743`, and analysis code `31941603`. The hiPSC bulk DOI is **not** the human heart atlas accession.

The actual [human Figshare API record](https://api.figshare.com/v2/articles/31939227) is version 1, DOI `10.6084/m9.figshare.31939227.v1`, posted 2026-04-05. It is public, not confidential, not embargoed, and downloads are enabled. Its license is [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). It contains 132 files totaling **3,431,641,059 bytes**: 44 `matrix.mtx.gz`, 44 `barcodes.tsv.gz`, and 44 `features.tsv.gz` files in 44 sample folders. The matrix files alone total 3,414,663,208 bytes. No expression file was downloaded.

The [code record](https://api.figshare.com/v2/articles/31941603) is likewise public, unembargoed, and CC BY 4.0. It contains 187 scripts totaling 590,883 bytes, not an annotated expression object. Three directly linked scripts were read: [public-data import](https://ndownloader.figshare.com/files/63467727), [DecontX processing](https://ndownloader.figshare.com/files/63467724), and [44+10-sample merging/QC](https://ndownloader.figshare.com/files/63467745).

## Verified 29-donor/54-tissue design

The paper reports 54 flash-frozen, transmural, nonfailing myocardial samples from 29 individuals, with 27 LV and 27 RV samples. Twenty-five donors contributed both chambers; four contributed one chamber. There are no atrial, apical, septal, or non-heart tissues among the 54 Jia analysis samples. The external HCA atlas contains other regions, but Jia explicitly subsets it to LV/RV nuclei.

| Group | Donors | LV/RV tissues | Age verified in table S1 | Source |
|---|---:|---:|---|---|
| AG1 | 5 | 5 / 5 | 13w+2d to 19w+5d | Jia |
| AG2 | 6 | 5 / 5 | 20w+3d to 27w+3d | Jia; PT45 LV only and PT42 RV only |
| AG3 | 7 | 6 / 6 | 29w+6d to 39w+3d | Jia; PT85 LV only and PT96 RV only |
| AG4 | 3 | 3 / 3 | 14, 20, 35 years | Jia |
| AG5 | 3 | 3 / 3 | 43, 48, 53 years | Jia |
| AG6 | 5 | 5 / 5 | 60–65, 60–65, 65–70, 70–75, 70–75 years | Litviňuková et al.; D2, D4, D5, D6, D7 |
| **Total** | **29** | **27 / 27** | Fetal gestational ages and postnatal 14–75-year range/bands | 44 Jia + 10 external tissues |

The primary [supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13274588/supplementaryFiles) was opened and only its index and table S1 were extracted for metadata review. Table S1 provides the following donor map. “Exact” below means the released value, not an invented midpoint. Fetal sex was inferred from sex-chromosome expression according to the paper, rather than independently clinically ascertained.

| Donor | Group/sample | Released age | Released sex | snRNA tissue |
|---|---|---|---|---|
| PT-13w | AG1-1 | 13w+2d | Male | LV+RV |
| PT-14w | AG1-2 | 14w+5d | Male | LV+RV |
| PT51 | AG1-4 | 17w+4d | Female | LV+RV |
| PT49 | AG1-5 | 18w+1d | Female | LV+RV |
| PT66 | AG1-6 | 19w+5d | Male | LV+RV |
| PT45 | AG2-1 | 20w+3d | Female | LV only |
| PT21 | AG2-2 | 22w+1d | Male | LV+RV |
| PT90 | AG2-3 | 23w+3d | Female | LV+RV |
| PT19 | AG2-4 | 25w | Female | LV+RV |
| PT42 | AG2-5 | 26w | Female | RV only |
| PT25 | AG2-6 | 27w+3d | Female | LV+RV |
| PT102 | AG3-1 | 29w+6d | Female | LV+RV |
| PT7 | AG3-2 | 30w+4d | Male | LV+RV |
| PT12 | AG3-3 | 30w+6d | Female | LV+RV |
| PT18 | AG3-4 | 32w+1d | Male | LV+RV |
| PT85 | AG3-5 | 33w+5d | Female | LV only |
| PT96 | AG3-6 | 37w | Male | RV only |
| PT50 | AG3-7 | 39w+3d | Male | LV+RV |
| HS316 | AG4-1 | 14 years | Male | LV+RV |
| HS257 | AG4-2 | 20 years | Male | LV+RV |
| HS235 | AG4-3 | 35 years | Male | LV+RV |
| HS197 | AG5-1 | 43 years | Male | LV+RV |
| HS173 | AG5-2 | 48 years | Male | LV+RV |
| HS246 | AG5-3 | 53 years | **Unresolved identifier conflict** | LV+RV |
| D2 | AG6-1 | 60–65 years | Male | LV+RV |
| D4 | AG6-2 | 70–75 years | Female | LV+RV |
| D5 | AG6-3 | 65–70 years | Female | LV+RV |
| D6 | AG6-4 | 70–75 years | Male | LV+RV |
| D7 | AG6-5 | 60–65 years | Male | LV+RV |

The sample-information sheet names **HS246**, but the gender sheet instead contains **HS247, Male** and has no HS246 row. The statistical sheet treats AG5 as three males. This is internally suggestive but not a clean donor-level mapping; HS246 sex must be reconciled before sex-adjusted modeling rather than silently relabeled.

Six additional internal AG6 donors—B2023-221 (73), B2023-275 (74), B2023-94 (62), B2023-223 (70), B2023-238 (66), and B2023-176 (71)—have `snRNA-seq = N`. Their exact ages must not be assigned to the five public AG6 snRNA donors or used to increase the inferential sample size.

### Health and procurement are not “uniformly healthy”

All 54 samples are described as nonfailing, and table S1 marks cardiac disease `N` for the 29 snRNA donors. Internal fetal samples were DCD; internal postnatal samples were DBD transplant-donor hearts excluded because of size or blood-type mismatch. External AG6 samples were DCD. The paper says cardiac developmental defects were excluded, but PT102 and PT18 have noncardiac fetal malformations recorded. D4 has hypertension and diabetes. For internal AG4/AG5 donors, hypertension, diabetes, cancer, and liver disease are `Unknown` (pulmonary disease is `N`). Thus “nonfailing” must not be rewritten as free of all disease or comorbidity.

## Assay, matrices, and released annotations

### Forty-four new Jia samples

The Methods specify 10x Chromium Single Cell 3′ v3.1 snRNA-seq, NovaSeq 6000, Cell Ranger 6.1.2, and `refdata-gex-GRCh38-2020-A`. The linked DecontX script expects `filtered_feature_bc_matrix` as the analysis input and separate `raw_feature_bc_matrix` as ambient-RNA background for exactly the 44 sample names represented in Figshare. The released record has one matrix/barcode/feature triplet per sample and no unfiltered background matrix.

Accordingly, the public triplets are consistent with the study’s filtered-count inputs, but their integer-valuedness, dimensions, and exact correspondence to the code inputs were not inspected. The unfiltered raw-feature matrices and DecontX-corrected matrices are not separately released in this record. Matrix-value validation remains a mandatory future gate; “raw counts available” must not be taken to mean raw sequencing reads or verified uncorrected integer matrices.

### Ten external AG6 samples

Jia’s Methods and code identify Litviňuková et al., [“Cells of the adult human heart”](https://www.nature.com/articles/s41586-020-2797-4), DOI `10.1038/s41586-020-2797-4`. The code reads `global_raw.h5ad` plus `Anno.csv`, selects `source == Nuclei`, LV/RV, and donors D2/D4/D5/D6/D7. This is a subset of an existing cohort, not five newly generated Jia donors and not an independent technical replicate of AG4/AG5.

The actual [HCA project record](https://service.azul.data.humancellatlas.org/index/projects/ad98d3cd-26fb-4ee3-99c9-8a2ab085e737) is project `ad98d3cd-26fb-4ee3-99c9-8a2ab085e737`, ENA `ERP123138` / `PRJEB39602`, with `accessible = true` and data-use restriction **`NRES`**. It lists 10x 3′ v2/v3 cell and nucleus assays, HiSeq 4000/NextSeq 500, six anatomical regions, and public raw reads/count objects. Relevant released object sizes are:

| External object | Bytes | Metadata status |
|---|---:|---|
| HCA contributor `hca_heart_global_ctl200723_freeze.h5ad` | 5,068,859,764 | Public contributor matrix; not downloaded |
| HCA `human-heart-10XV3.loom` | 3,255,905,369 | Public count matrix; not downloaded |
| HCA `human-heart-10XV2.loom` | 4,398,007,571 | Public count matrix; not downloaded |
| CELLxGENE all-cells/nuclei H5AD, dataset `d4e69e01-3ba2-4d6b-a15d-e7048f78f22e` | 2,953,098,273 | Public, raw counts in `raw.X`, annotations, 486,134 cells/nuclei; revised 2026-06-11; not downloaded |
| HCA raw FASTQ inventory, 768 files | 3,054,230,259,834 | Public but unnecessary for this planned processed-count analysis; not downloaded |

The [CELLxGENE collection API](https://api.cellxgene.cziscience.com/curation/v1/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72) confirms donor IDs D1–D7/D11/H2–H7, normal disease labels, both sexes, cells plus nuclei, six tissues, and raw counts in `raw.X`. Its current curated object is not automatically the same version Jia used; a future reproduction must select nuclei and document the object/version rather than treating every cell in the all-cells object as AG6.

### Annotation gap

For the 44 Jia samples, Figshare has only matrix/barcode/feature triplets; the code record has scripts; and supplementary tables provide marker gene sets and aggregate proportions. **No Jia barcode-level cell-type/state annotation object was found in these checked records.** The 442,239-nucleus, 11-major-cell-type labels are reported in the paper, but a per-barcode release mapping those labels to the 44 public matrices was not located. External HCA annotations are available, but they do not supply Jia’s labels for the new samples.

A future analysis must either locate a versioned Jia annotation object or predeclare and reproduce annotation from the released counts and marker tables. It must not import HCA labels onto Jia nuclei or treat supplementary aggregate proportions as per-cell annotations.

## Public terms and controlled access

- Jia human, mouse, bulk, and code Figshare records checked through their public APIs are public, unembargoed, downloadable, and CC BY 4.0.
- The external HCA project is `NRES`, accessible, and its raw reads and processed objects are publicly listed. CELLxGENE labels the collection `PUBLIC`.
- No dbGaP/EGA-controlled dataset, credential requirement, embargo, or data-use agreement was stated or encountered for the processed objects needed here. **No controlled-access component was verified.** This does not create exact AG6 ages: their five-year bands are a metadata limitation, not an access-control problem.
- No repository was contacted for authorization and no access request was made.

## Required future analysis boundary

1. **Inferential unit:** donor, not nucleus or chamber. LV and RV from one donor require donor-level aggregation, a donor-blocking term, or joint permutation. Do not count 54 tissues as 54 independent donors or add the six non-snRNA B2023 donors.
2. **Adult primary:** a full 14–75 exact-age model is blocked. A six-internal-donor exact-age sensitivity cannot represent ages 60–75. A 6-versus-5 AG4/5-versus-AG6 contrast is descriptive only because source and group are inseparable.
3. **Lifespan sensitivity:** analyze fetal and postnatal stages separately or use declared stage/group contrasts. Do not convert gestational weeks and postnatal years into an undocumented continuous lifespan variable. Preserve AG6 source, age-band, DCD/DBD, and ethnicity differences as limitations.
4. **Before any outcome:** validate finite, nonnegative integer filtered counts; establish whether uncorrected or DecontX counts are being used; verify unique barcodes and gene builds; obtain/reproduce Jia annotations; and verify D2/D4/D5/D6/D7 LV/RV nuclei joins in the selected HCA object/version.
5. **Eligibility:** repeat the standard gate per cell type before examining effects, including at least 20 annotated nuclei per donor-cell type, at least six eligible donors per model with residual degrees of freedom, and at least 50 cross-fitted young-reference nuclei plus 20 passing nuclei per scored donor for Hellinger. Ineligible units are descriptive/untestable, not null results.
6. **Endpoints and multiplicity:** freeze endpoints, age contrasts, identifiable covariates, and BH families before outcomes. Keep P1 transport, donor pseudobulk, CLR composition, and Hellinger separate. Non-significance is not equivalence.

## Not done and next decision

- No expression data or count values were acquired; no heavy job was run; no aging outcome was examined; no alternative heart dataset was substituted.
- No controlled-access request or repository contact was made. No ovary file was modified. No commit or push was performed.
- The accession blocker in Step 8 is resolved. The remaining blockers are exact AG6 ages, the HS246/HS247 sex conflict, absent Jia barcode-level annotations, unchecked matrices/joins, and complete AG6 source-age confounding.
- Parent review is required before any future expression acquisition or analysis. A future full adult continuous-age claim should remain blocked unless exact AG6 ages are lawfully available from a primary public metadata record; no midpoint imputation is permitted.

## Verification

- Parsed `data/manifest/jia_heart.json` and cross-checked every inventory tuple against the downloaded Figshare API record: 29 unique donors, 54 tissues, 44 Figshare folders, 132 exact file IDs/sizes, 3,431,641,059 total bytes, and matching supplied/computed MD5 fields in the source record.
- `git diff --check`: passed.
- Full suite with `/home/reza/miniconda3/envs/agescape/bin/python -m pytest -q`: **139 passed, 1 skipped, 1 warning in 6.45s**.
- Concurrent ovary-agent changes appeared in the worktree during this gate; no ovary file was read for modification or modified by this heart task.
