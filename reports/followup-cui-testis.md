# Cui testis GSE254315: mapping, frozen analysis plan, and descriptive QC

Base revision: `8cd58bdfc0921636ba2846f898fc6d0ab96d3ccb`. Metadata and plan freeze: 2026-11-14.

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

## Decision

The mandatory mapping is **resolved for GSE254315**: its 25 GEO records are 25 separate named donor libraries, not pools and not 35 independent measurements. The paper's 35-donor integrated atlas comprises these 25 plus ten reused external-only donor IDs. This step uses only the 25 GSE donors. No public donor-and-barcode cell-type crosswalk was found, so P1, donor/cell-type pseudobulk, CLR composition, and Hellinger analyses were **not run**. Raw-count validation and descriptive QC were completed; no expression outcome or age association was fitted.

Primary sources: [Cui et al., Nature Aging (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12003174/), [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/soft/GSE254315_family.soft.gz), [GEO file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/suppl/filelist.txt), and the primary Fig. 1a/Extended Data Fig. 1/Fig. 6 source workbooks recorded in `data/manifest/GSE254315.json`.

## The 25-versus-35 mapping

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

The primary Fig. 1a source is ordered Donor1 through Donor35. Its ages match GEO exactly for 24 of the 25 GSE donors. The sole discrepancy is **Donor8: 40 years in the primary source versus 39 in GEO**. The frozen primary uses 40/40s and declares a complete 39/30s sensitivity. The ten atlas-only donor IDs are Donor1–3, Donor26–29, and Donor31–33. The paper names GSE182786 and GSE215754 as external sources; those donors are excluded here and their ages are never added as independent observations.

| Donor | GSM | Library | Primary age/decade | GEO age | BMI | Raw barcodes | Paper-QC pass |
|---|---|---|---:|---:|---:|---:|---:|
| Donor4 | GSM8037606 | 20149X1 | 30/30s | 30.0 | 20.9 | 13426 | 10822 |
| Donor5 | GSM8037607 | 20188X1 | 32/30s | 32.0 | 25.6 | 5505 | 3698 |
| Donor6 | GSM8037621 | 20206X3 | 33/30s | 33.0 | 30.7 | 6208 | 3923 |
| Donor7 | GSM8037610 | 20188X4 | 38/30s | 38.0 | 30.8 | 4759 | 3307 |
| Donor8 | GSM8037622 | 20206X4 | 40/40s | 39.0 | 31.6 | 7171 | 3889 |
| Donor9 | GSM8037623 | 20206X5 | 42/40s | 42.0 | 31.9 | 7789 | 4658 |
| Donor10 | GSM8037624 | 20206X6 | 43/40s | 43.0 | 20 | 7390 | 4448 |
| Donor11 | GSM8037608 | 20188X2 | 43/40s | 43.0 | 30.1 | 5965 | 3371 |
| Donor12 | GSM8037611 | 20196X1 | 43/40s | 43.0 | 31.9 | 6594 | 3926 |
| Donor13 | GSM8037625 | 20207X1 | 43/40s | 43.0 | 32 | 2748 | 1803 |
| Donor14 | GSM8037626 | 20207X2 | 44/40s | 44.0 | 25.8 | 5135 | 3486 |
| Donor15 | GSM8037619 | 20206X1 | 44/40s | 44.0 | 44.9 | 6742 | 3801 |
| Donor16 | GSM8037615 | 20197X1 | 46/40s | 46.0 | 41.5 | 4074 | 2787 |
| Donor17 | GSM8037616 | 20197X2 | 48/40s | 48.0 | 32.4 | 3988 | 2547 |
| Donor18 | GSM8037612 | 20196X2 | 52/50s | 52.0 | 37.7 | 5826 | 3247 |
| Donor19 | GSM8037627 | 20207X3 | 53/50s | 53.0 | 34.3 | 7257 | 3029 |
| Donor20 | GSM8037613 | 20196X3 | 55/50s | 55.0 | 20.8 | 6966 | 5044 |
| Donor21 | GSM8037609 | 20188X3 | 55/50s | 55.0 | 32.6 | 5433 | 3594 |
| Donor22 | GSM8037617 | 20197X3 | 57/50s | 57.0 | 35.7 | 4123 | 3059 |
| Donor23 | GSM8037618 | 20197X4 | 57/50s | 57.0 | 46.5 | 5201 | 3604 |
| Donor24 | GSM8037614 | 20196X4 | 59/50s | 59.0 | 25.4 | 6320 | 3949 |
| Donor25 | GSM8037620 | 20206X2 | 59/50s | 59.0 | 34.5 | 6998 | 4344 |
| Donor30 | GSM8037629 | 20207X5 | 66/60s | 66.0 | 26.2 | 6495 | 4311 |
| Donor34 | GSM8037630 | 20207X6 | 67/60s | 67.0 | 30.8 | 1761 | 1272 |
| Donor35 | GSM8037628 | 20207X4 | 69/60s | 69.0 | 38.4 | 4866 | 3602 |

The primary decade distribution is 30s=4, 40s=10, 50s=8, 60s=3; N=25. The 20s donors in the integrated figure are external-only and are not in GSE254315.

## Frozen pool/donor-level plan

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

The plan in `data/manifest/cui_testis_analysis_plan.json` was frozen after mapping and before expression download or outcomes. Its SHA-256 is recorded in the metadata manifest. If a valid annotation crosswalk is later released, the primary is an ordinal-decade effect (3/4/5/6) per decade with donor mean log actual raw UMI adjustment over all 25 donors. The 30s are young/reference, 40s are middle, and 50s/60s are older decades.

Declared sensitivities are: unadjusted ordinal; BMI-plus-depth ordinal; categorical decades with 30s reference; old 50s/60s versus young 30s with middle 40s excluded only there; exact age per ten years; Donor8=39/30s versions of the applicable models; and incomplete-spermatogenesis-adjusted and complete-only models **only if** a released donor-linked status exists. No exact pool ages are invented and no cell-level p-values are permitted.

Eligible released labels would require at least 20 cells per donor/cell type, six donors overall, three per retained group, full rank, and two residual degrees of freedom. Hellinger additionally requires at least two young reference donors, 20 qualifying cells per scored donor, 50 cross-fitted reference cells, and actual-UMI >=1500. Somatic labels and each comparable germ-cell stage are separate. P1 uses the established human P1_Clean and a 200-draw expression-decile-matched null (seed 20261114). Effects require standard errors, 95% CIs, and explicit BH families as frozen in the plan.

## Annotation gate and blocked endpoints

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

GEO releases raw MatrixMarket counts, barcodes, and features only. Primary Fig. 1a releases donor age/BMI, Extended Data Fig. 1 releases age-aggregated cell totals, and Fig. 6a releases donor-linked elongating-spermatid frequencies. None is a per-cell donor/barcode annotation crosswalk. The released code repository contains notebooks that refer to local metadata and H5AD objects but no public annotation files. Marker-based relabeling, reclustering, or treating age-aggregated totals as donor observations would violate the input rules.

| Endpoint | Somatic labels | Comparable germ-cell stages | Status | Effect/CI/p/q |
|---|---:|---:|---|---|
| P1 | 8 | 6 | blocked; not run | none |
| donor_pseudobulk | 8 | 6 | blocked; not run | none |
| CLR | 8 | 6 | blocked; not run | none |
| Hellinger | 8 | 6 | blocked; not run | none |

All 56 label-by-endpoint rows and their failed annotation reason are in `results/followup-cui-testis-endpoint-status.csv`. They are blocked, not null results, and enter no BH denominator. No whole-testis pseudobulk model was substituted because it could not separate cell-intrinsic changes from somatic/germ-stage composition.

## Raw-count validation and descriptive QC

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

All 75 per-donor files were downloaded sequentially under `scripts/resource_guard.py` with a 30 GB reserve and 0.2 GB per-file cap; SHA-256 and expected byte sizes passed. All stored count values in all 25 matrices are finite, nonnegative integers. Matrix dimensions match barcodes/features; barcodes and feature IDs are unique. Decompressed feature content is identical across all donors (canonical SHA-256 `28fa7d6a0994e5fe6c1d2a39a3363ee1655916481e3ee240efd6cc8a361959c9`), although gzip-container checksums can differ.

The matrices contain **148,740 raw barcodes** and **805,486,323 UMIs**, whereas the GEO summary reports 99,130 cells. Applying the paper's released thresholds descriptively (200–5500 detected features, 800–17000 UMI, mitochondrial fraction <=25%) retains **95,521** barcodes (64.22%), **3,609 fewer** than the GEO summary. This difference is unresolved because final annotated barcodes are not public; it is not silently treated as a successful reproduction. These are not reclustered or annotated cells and were not used for age-associated QC tests. Per-donor UMI, feature, mitochondrial, and QC summaries are in `results/followup-cui-testis-qc.csv`.

The guarded QC process had maximum RSS 886.0 MiB, below the 8 GB cap; one donor was processed at a time. External checkpoints are under `../agescape-revision-data/cui-testis/interim/qc/`.

## Spermatogenesis information

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

Primary Fig. 6a provides donor-linked elongating-spermatid `Freq` values for 24 of 25 GSE donors; Donor34 is absent and remains missing in the aggregate CSV. The Methods define infertility as fewer than three sperm per seminiferous tubule, but the checked public source does not release donor-level PAS sperm counts or an incomplete-spermatogenesis indicator. Therefore the spermatogenesis-adjusted and complete-only sensitivity models are **untestable**, and no status was inferred from germ-cell frequencies. Released age-aggregated label totals are retained descriptively in `results/followup-cui-testis-released-age-cell-totals.csv`, not treated as independent donors.

## What survived, died, and remains blocked

**Biological caveat for every result:** changing germ-cell-stage composition can mimic broad expression displacement without cell-intrinsic aging drift. Testis hits or nulls cannot support general whole-body aging claims.

- Survived: one-to-one 25-donor GEO mapping; exact primary/GEO age reconciliation with one declared discrepancy; raw-count integrity; descriptive QC; and release of descriptive age/BMI and germ-stage source fields.
- Died as an inferential shortcut: assigning all 35 atlas ages to GSE254315, treating cells as replicates, or treating the ten external donors as new libraries.
- Blocked: every cell-resolved P1, gene-effect, CLR, and Hellinger estimand, including their matched-null and sensitivity models, because the required released annotation crosswalk is absent. Non-significance is neither reported nor implied.
- Contradiction/complication: the 35-donor paper atlas is not a 35-library GSE cohort; Donor8 also has a one-year primary-versus-GEO discrepancy. These do not justify guessing labels or ages.
- Not attempted: no FASTQ/RAW tar, external-donor analysis, clustering, marker imputation, expression outcome modeling, cell-level tests, or ovary work.

Reproduce from the repository in the `agescape` environment:

```bash
python scripts/prepare_followup_cui_testis.py
python scripts/fetch_followup_cui_testis.py
python scripts/resource_guard.py run --path ../agescape-revision-data/cui-testis --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_cui_testis.py
python scripts/report_followup_cui_testis.py --tests '<full-suite result>'
```

## Verification and required follow-up

Full-suite verification: `139 passed, 1 skipped, 1 warning in 6.60s`. No commit or push was made by the testis agent.

The remaining prerequisite is a **public, released donor-and-barcode-to-cell-type annotation crosswalk** (or annotated public object with raw counts) covering these libraries. If supplied, the frozen design can be evaluated without inventing labels. Any design change requires a new freeze before outcomes. Donor mapping itself does not require further user input.
