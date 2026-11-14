# PBMC metadata and public-access gate

Base revision: `67634307987dcf5f7761f39cd2c31576f8c16952`. Verification date: 2026-11-14.

## Scope and headline decisions

This was a lightweight metadata/access check of **GSE157007, GSE213516, and GSE214546 only**. No complete expression matrix, RAW tar, FASTQ, or integrated expression object was downloaded. GEO SOFT/file inventories, primary papers/supplements, small annotation/barcode probes, HTTP HEAD requests, and HDF5 schema-only range reads were used. No biological outcome was analyzed. GSA and Synapse arms were not accessed.

The exact URLs, companion barcode/feature files, byte sizes, and donor mappings are in `data/manifest/pbmc_geo_cohorts.json`.

| Cohort | Strict young <=40 / old >60 | Actual Filippov Table S1 | Metadata/access | Separate-cohort analysis verdict |
|---|---:|---:|---|---|
| GSE157007 | **8**: 2 young, 6 old; 3 F/5 M | 8 | Public raw counts pass; per-cell labels do not | **Blocked for donor-level inference**: only two young donors and no public per-cell annotation crosswalk found |
| GSE213516 | **10**: 5 young, 5 old; 6 F/4 M | **12**, because F60 and M60 are called old | Ready with conditions | P1, pseudobulk, CLR, and Hellinger ready after count/label/eligibility checks; resolve hg19/GRCh38 conflict for gene models |
| GSE214546 | **23**: 16 young adults, 7 old; 13 F/10 M | **25**, because PB00015 and PB00566 (both 60) are called old | Ready with conditions | All four endpoints ready after per-file count, label, donor/cell-floor, CMV, and pool checks |

The cutoff discrepancy is real. [Filippov et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11606963/) say young excludes donors older than 40 and old retains donors **over 60**, but their [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41514-024-00185-x/MediaObjects/41514_2024_185_MOESM2_ESM.xlsx) includes age-60 donors in the old group. Their [GSE214546 import code](https://raw.githubusercontent.com/f6v/PBMC_ageing_atlas/master/GSE214546/create_adatas.py) skips old samples only when `age < 60`, confirming the inclusive-age-60 implementation. The strict `>60` contrast is primary below; reproducing the atlas's `>=60` set is a separately labeled sensitivity, not a silent change.

## GSE157007 — Luo et al. 2022

Primary sources were the [paper](https://www.nature.com/articles/s43587-022-00198-9), [Extended Data Fig. 1](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs43587-022-00198-9/MediaObjects/43587_2022_198_Fig1_HTML.png), [primary supplementary workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-022-00198-9/MediaObjects/43587_2022_198_MOESM2_ESM.xlsx), [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157007/soft/GSE157007_family.soft.gz), and [exact file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157007/suppl/filelist.txt). The paper is Luo et al., *Nature Aging* 2, 348–364 (2022), DOI `10.1038/s43587-022-00198-9`.

### Donors, exclusions, and libraries

Extended Data Fig. 1a supplies exact ages and sex. Use only the original Young and non-frail Old controls. Exclude three cord-blood samples, five Frail samples, and F013 at age 42. The six Old controls are **non-frail rather than demonstrably disease-free**: primary Supplementary Table 1 records nonzero chronic-disease and disability scores for several of them. This cohort must not be redescribed as uniformly healthy.

| Donor | scRNA GSM | Age | Sex | Strict/atlas group | Matrix bytes |
|---|---|---:|---|---|---:|
| F012 | GSM4750303 | 27 | M | young | 36,878,290 |
| F013 | GSM4750304 | 42 | M | excluded: 41–60 | 57,020,995 |
| F014 | GSM4750305 | 23 | F | young | 39,459,879 |
| F020 | GSM4750309 | 78 | M | old | 45,446,507 |
| F021 | GSM4750310 | 72 | M | old | 36,719,329 |
| F023 | GSM4750311 | 80 | F | old | 53,233,979 |
| OH14 | GSM5684306 | 88 | F | old | 59,163,654 |
| OH15 | GSM5684307 | 97 | M | old | 65,791,732 |
| OH17 | GSM5684308 | 100 | M | old | 66,867,230 |

Each donor has one scRNA GSM. Separate V(D)J and surface-protein GSMs are modality-specific libraries, not extra donors or independent expression replicates. The assay is 10x 5' scRNA with separate TCR and CCR7/CD45RA/CD4/CD8 feature-barcode arms, NovaSeq 6000, Cell Ranger 3.1.0, GRCh38 v1.2.0.

### Counts, annotations, and public objects

- Per-donor MatrixMarket matrices, barcodes, and features are public. A 4,096-byte range probe of F012 confirmed a `matrix coordinate integer general` header and dimensions; no full-matrix integer-value validation was done. The series [RAW tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157007/suppl/GSE157007_RAW.tar) is 740,761,600 bytes and was not downloaded.
- The primary paper documents 17 broad clusters and T/monocyte subtypes. Primary Supplementary Table 4 provides donor-by-cell-type totals, but the checked releases do **not** provide a per-cell barcode-to-label crosswalk. Those totals can support descriptive comparison only; they cannot assign cells for P1, pseudobulk, or Hellinger.
- No cohort-released processed per-cell expression object was linked in GEO. The [Filippov code repository](https://github.com/f6v/PBMC_ageing_atlas) provides scripts, not a public integrated matrix/annotation download in the checked record. Integrated labels must not be imported without a barcode crosswalk.

**Endpoint verdicts:** separate-cohort P1, donor pseudobulk, CLR inference, and Hellinger inference are all **blocked**. Two young donors fail the prespecified three-donor/group floor, and the missing per-cell map independently blocks cell-resolved endpoints. Released composition totals may be shown descriptively, but must not become a new two-versus-six significance test.

## GSE213516 — Zhu et al. 2023

Primary sources were the [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213516/soft/GSE213516_family.soft.gz), [file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213516/suppl/filelist.txt), [paper full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10306289/), and [Figshare annotation API](https://api.figshare.com/v2/articles/23789844). The paper is Zhu et al., *Science Advances* 9, eabq7599 (2023), DOI `10.1126/sciadv.abq7599`.

### Exact donor map and filtering

All 17 volunteers were evaluated as generally healthy from medical history and vital signs. There is one GSM and one donor per sample; no additional disease exclusions are indicated. Middle-aged donors remain available for a separately declared exact-age sensitivity but are absent from the strict group contrast.

| Donor/GSM | Age | Sex | Strict group | Atlas group | Matrix bytes |
|---|---:|---|---|---|---:|
| F30 / GSM6588511 | 30 | F | young | young | 38,863,486 |
| F36 / GSM6588512 | 36 | F | young | young | 60,468,162 |
| F47 / GSM6588513 | 47 | F | excluded | excluded | 47,630,412 |
| F58 / GSM6588514 | 58 | F | excluded | excluded | 47,465,537 |
| F66 / GSM6588515 | 66 | F | old | old | 46,460,156 |
| F77 / GSM6588516 | 77 | F | old | old | 46,946,017 |
| F55 / GSM6588517 | 55 | F | excluded | excluded | 52,997,371 |
| F60 / GSM6588518 | 60 | F | **excluded: not >60** | **old** | 73,474,924 |
| F40 / GSM6588519 | 40 | F | young | young | 43,296,484 |
| F65 / GSM6588520 | 65 | F | old | old | 64,387,208 |
| F46 / GSM6588521 | 46 | F | excluded | excluded | 74,805,277 |
| M28 / GSM6588522 | 28 | M | young | young | 44,796,555 |
| M36 / GSM6588523 | 36 | M | young | young | 39,981,736 |
| M45 / GSM6588524 | 45 | M | excluded | excluded | 41,996,868 |
| M60 / GSM6588525 | 60 | M | **excluded: not >60** | **old** | 47,144,721 |
| M66 / GSM6588526 | 66 | M | old | old | 49,528,916 |
| M74 / GSM6588527 | 74 | M | old | old | 50,794,979 |

Thus the strict cohort is **10 donors, not 12**: five young (F30, F36, F40, M28, M36) and five old (F65, F66, F77, M66, M74). Strict sex counts are 6 F/4 M. The atlas's 12 donors are 7 F/5 M. Neither set is sex-restricted.

### Counts, annotations, and public objects

- GEO supplies per-donor MatrixMarket matrices, barcodes, and features. A 4,096-byte range probe of F30 confirmed an integer MatrixMarket header. The [RAW tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213516/suppl/GSE213516_RAW.tar) is 876,963,840 bytes and was not downloaded. Raw reads are withheld for privacy, but these processed raw-count matrices are public.
- The public Figshare v1 file [`SE_CYCT_JOCT_WHCT_SC_meta.csv`](https://ndownloader.figshare.com/files/41726226) is 18,499,390 bytes, MD5 `f258e3ce234de79118821119925c2efd`. Its header and initial rows verify `X`, `Age`, `sampleID`, `Cohort`, `Celltype`, `UMAP1`, and `UMAP2` fields. This is annotation/UMAP metadata, not an expression matrix; other cohorts in the CSV are excluded.
- The barcode join is deterministic but must be validated: for `Cohort == "SE"`, `X` equals `sampleID + "_"` plus the 16-base GEO barcode without `-1`. For example, GEO `AAACCTGAGAGACGAA-1` maps to `F30_AAACCTGAGAGACGAA`. Require unique, one-to-one, full-coverage matches per selected sample before using labels.
- The assay is 10x 3' scRNA, Cell Ranger 3.0.0, HiSeq PE150/GEO GPL20301. **Genome-build provenance conflicts:** GEO says `hg19`, while the primary Methods say GRCh38. Resolve feature IDs/reference provenance before P1 or all-gene pseudobulk analysis; do not silently choose one statement.

**Endpoint verdicts:** P1, donor pseudobulk, CLR, and Hellinger are **ready conditionally**. Five young and five old donors pass the cohort-level floor. Each cell type still needs >=20 cells/donor, >=3 donors/group, >=6 scored donors, full rank, and the Hellinger reference/depth checks. Gene-level work additionally needs the genome-build reconciliation.

## GSE214546 — Thomson et al. 2023

Primary sources were the [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214546/soft/GSE214546_family.soft.gz), [file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214546/suppl/filelist.txt), [paper full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10602854/), and [primary Supplementary Table 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-023-01641-8/MediaObjects/41590_2023_1641_MOESM3_ESM.xlsx). The paper is Thomson et al., *Nature Immunology* 24, 1947–1959 (2023), DOI `10.1038/s41590-023-01641-8`.

### Arm selection, sex, and donor map

Use the healthy **confirmatory scRNA PBMC arm**, which has 16 pediatric, 16 young-adult, and 16 older-adult donors with balanced sex design. Exclude pediatric donors from the Filippov-replication young arm even though their numeric ages are <=40; otherwise age filtering would introduce a developmental contrast absent from the atlas subset. Also exclude ages 41–60 from the strict contrast, all TEA-seq samples, IMM19 technical controls, stimulated/fixed-RNA arms, and pooled/unknown-age controls. CMV-positive healthy donors remain eligible; CMV is a covariate/sensitivity factor, not a post hoc exclusion.

The **TEA-seq discovery arm is sex-restricted to females** (eight pediatric and eight older donors) and is excluded. The selected scRNA arm is mixed sex. Strict selection has 13 F/10 M; atlas selection has 13 F/12 M.

The confirmatory assay is modified 10x 3' scRNA with cell hashing, eight donors pooled per library plus a common control, NovaSeq, and Ensembl GRCh38-3.0.0. GSMs below are per-donor demultiplexed labeled H5 objects; batch/pool identity remains in the object and must not be mistaken for donor replication.

| Donor / subject | GSM | Age | Sex / CMV | Strict group | H5 bytes |
|---|---|---:|---|---|---:|
| PB00002 / BR1002 | GSM6929079 | 29 | M / Neg | young adult | 61,725,418 |
| PB00003 / BR1003 | GSM6929080 | 31 | F / Neg | young adult | 60,470,432 |
| PB00004 / BR1004 | GSM6929082 | 31 | M / Neg | young adult | 63,955,392 |
| PB00006 / BR1005 | GSM6929083 | 28 | F / Neg | young adult | 74,300,299 |
| PB00010 / BR1006 | GSM6929111 | 32 | M / Neg | young adult | 61,380,274 |
| PB00015 / BR2003 | GSM6929084 | 60 | M / Neg | **excluded; atlas old** | 66,912,167 |
| PB00022 / BR1010 | GSM6929081 | 30 | F / Pos | young adult | 61,072,574 |
| PB00023 / BR1009 | GSM6929112 | 27 | F / Neg | young adult | 68,991,299 |
| PB00026 / BR1012 | GSM6929113 | 32 | M / Pos | young adult | 56,905,201 |
| PB00030 / BR1013 | GSM6929085 | 32 | F / Neg | young adult | 82,892,346 |
| PB00031 / BR1014 | GSM6929114 | 27 | F / Pos | young adult | 65,027,099 |
| PB00041 / BR1018 | GSM6929086 | 31 | M / Neg | young adult | 71,661,908 |
| PB00352 / BR2010 | GSM6929110 | 61 | F / Pos | old | 88,142,620 |
| PB00363 / BR1040 | GSM6929089 | 29 | M / Neg | young adult | 77,844,808 |
| PB00366 / BR1019 | GSM6611359 | 35 | F / Neg | young adult | 83,772,948 |
| PB00384 / BR1016 | GSM6929087 | 34 | M / Neg | young adult | 67,974,137 |
| PB00386 / BR1043 | GSM6929091 | 35 | M / Pos | young adult | 62,306,351 |
| PB00388 / BR2016 | GSM6929088 | 62 | F / Pos | old | 78,578,850 |
| PB00507 / BR1017 | GSM6929090 | 30 | F / Neg | young adult | 71,636,240 |
| PB00537 / BR2014 | GSM6929095 | 62 | F / Pos | old | 70,094,258 |
| PB00539 / BR2020 | GSM6929101 | 63 | M / Neg | old | 70,658,900 |
| PB00552 / BR2024 | GSM6929094 | 62 | F / Pos | old | 65,244,971 |
| PB00562 / BR2023 | GSM6929096 | 61 | M / Neg | old | 67,107,812 |
| PB00566 / BR2026 | GSM6929103 | 60 | M / Neg | **excluded; atlas old** | 68,406,896 |
| PB00569 / BR2018 | GSM6929098 | 63 | F / Pos | old | 66,637,524 |

Strict exclusion of PB00015 and PB00566 gives **23 donors, not 25**. Older donors aged 55–59 are also excluded; pediatric, TEA, stimulation, and technical-control records are outside this 25-row atlas candidate set.

### Counts, annotations, and public objects

- Every listed donor has a public `*_labeled.h5`. A schema-only check of PB00002 used 34 range requests totaling 2,218,730 bytes. It found `/matrix/data` dtype `uint16`, 28,471,702 stored entries, 33,538 features, 22,592 cells, and per-cell `seurat_pbmc_type`, score, barcode, batch/pool, and QC fields. No expression values were extracted. All other files still need equivalent schema and count-value validation after guarded acquisition.
- These labeled H5 files are processed public objects containing the raw count matrix and cohort-released Seurat PBMC annotations. Integrated Filippov labels are unnecessary. ATAC fragments and TEA objects are outside this plan.
- The series [RAW tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214546/suppl/GSE214546_RAW.tar) is **32,088,821,760 bytes** and was not downloaded. Use small per-donor H5 files, not this tar.

**Endpoint verdicts:** P1, donor pseudobulk, CLR, and Hellinger are **ready conditionally**. Apply donor/cell floors and test sex, CMV, pool/batch, and depth identifiability. Do not add pediatric or female-only TEA cells to enlarge the young reference or common cell types.

## Prespecified separate-cohort replication plan

1. **Freeze inputs before outcomes.** Use the strict 8/10/23 adjudication above. Validate finite, nonnegative integer counts; unique GSM/donor mapping; barcode-label joins; and released annotation universes. Keep GSE157007 descriptive. Exact-age models and age>=60 atlas-reproduction models are separate sensitivities, with no midpoint imputation.
2. **Common eligibility.** Per cohort and harmonized cell type require >=20 annotated cells per donor, >=3 donors per age group, >=6 scored donors, and residual degrees of freedom. Record every ineligible cell type. P1 requires at least one detected up- and down-program gene. Hellinger additionally requires >=2 young-reference donors, >=50 cross-fitted young-reference cells, >=20 passing cells per scored donor, and actual-UMI adjustment where identifiable.
3. **P1 transport.** Compute donor-by-cell-type P1 from cohort raw counts using the established library normalization and up-minus-down score. Fit old-minus-young effects with sex where identifiable; add depth as a sensitivity. Report effects and 95% CIs per cohort and BH within each cohort/model family.
4. **Donor pseudobulk.** Sum raw counts within donor and cell type; never sum across donors before modeling. The targeted CD8-naive replication set is the six published Filippov-supported increases: `PTGER2`, `CLIC1`, `S100A4`, `PRF1`, `GZMA`, and `CCL5`. Test these within each eligible cohort with a six-target multiplicity family. Any all-gene arm is exploratory, with BH separately by cohort and cell type. No cell-level Wilcoxon p-values and no cross-study p-value pooling.
5. **CLR composition.** Use one declared released-label universe per cohort, a 0.5-cell pseudocount, donor normalization, and centered log-ratio transformation. Targeted replication concerns CD8-naive T-cell decrease, MAIT decrease, and non-classical-monocyte increase. Other harmonized types are exploratory. Raw-fraction effects are descriptive sensitivities, not replacements for CLR.
6. **Hellinger.** Compute donor mean per-cell Hellinger distance to a cross-fitted strict-young reference independently in each cohort. Subtract a donor's own cells when that donor belongs to the reference. Apply the established technical-gene exclusions and depth adjustment; do not relax floors for GSE157007.
7. **Replication reporting.** Harmonize only defensible cohort-released labels. Produce forest plots of old-minus-young effects and 95% CIs for GSE213516 and GSE214546, plus direction concordance and explicit untestable rows. Do not concatenate matrices, pool donors to increase N, combine p-values, or call non-significance equivalence. CMV/pool effects and the age-60 sensitivity remain visible.

## Not done and next gate

- No expression analysis, count-value outcome check, QC threshold selection, differential expression, P1 score, CLR fraction, Hellinger distance, or replication estimate was computed.
- Only representative MatrixMarket headers and one H5 schema were inspected. Full per-file count validation and per-cell label coverage remain prerequisites.
- No integrated labels were imported. No GSA/Synapse data, muscle files, authorization channels, or external repositories for access requests were used.
- Future acquisition must be sequential under `scripts/resource_guard.py`, with 30 GB free reserve, per-file limits, <=8 GB analysis memory, and small donor batches. Parent review is required before expression work.
