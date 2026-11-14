# Muscle metadata and public-access gate

Base revision reviewed: `a33576495be8fa4b4c53701c2418f86eb85e140b`. Verification date: 2026-11-14.

## Scope and decisions

This was a lightweight metadata/access check only. No expression matrix was downloaded, no count values were analyzed, and no biological outcome was computed. Small text/API requests, supplementary metadata, HTTP HEAD requests, and HDF5 schema-only byte-range reads were used. HDF5 checks were limited to group names, dimensions, data types, and categorical labels; no H5/H5AD expression values were read.

| Study | Public metadata/access | Count-layer status | Analysis-age decision | Verdict |
|---|---|---|---|---|
| Kedlian 2024, E-MTAB-13874 | Ready: public BioStudies, supplement, and atlas objects verified | `layers/counts` exists but is float32; integer-valuedness not tested | Exact ages are not public in the checked supplement/SDRF; use young/old labels, not invented midpoints | Ready for later guarded count-value validation; continuous-age analysis blocked |
| Perez 2022, GSE167186 | Ready: public GEO sample records and per-donor files verified | All 17 filtered H5 `/matrix/data` arrays are int32 | Exact numeric ages available for all 17 snRNA donors | Ready for later guarded count-file validation |
| Cross-study replication label | Blocked | Perez is explicitly reused in Kedlian's comparisons; individual donor disjointness is not proven by public anonymized IDs | Do not pool donors or double-count Perez evidence | Do not call one study an independent replication of the other yet |

“Ready” does not mean that an expression download or analysis was performed. No Lake file was changed, and no external repository was contacted for authorization.

## Kedlian 2024 — E-MTAB-13874

Primary sources read were the [paper](https://www.nature.com/articles/s43587-024-00613-3), [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-024-00613-3/MediaObjects/43587_2024_613_MOESM3_ESM.xlsx), [BioStudies study JSON](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-13874), [public-release/file information](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-13874/info), [file inventory](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-13874/files), [SDRF](https://www.ebi.ac.uk/biostudies/files/E-MTAB-13874/E-MTAB-13874.sdrf.txt), [IDF](https://www.ebi.ac.uk/biostudies/files/E-MTAB-13874/E-MTAB-13874.idf.txt), and [atlas download page](https://www.muscleageingcellatlas.org/). The paper's [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11108788/fullTextXML) was also read to check the documented Perez reuse.

### Donors, ages, sex, assays, and health fields

- The supplement and SDRF agree on **17 unique donors and 43 human libraries**: 21 cell/scRNA libraries and 22 nuclei/snRNA libraries.
- scRNA covers 12 donors (5 young, 7 old); snRNA covers 13 donors (6 young, 7 old). These assay-specific counts overlap and must not be added to create 25 donors.
- Overall there are 8 young and 9 old donors, 10 male and 7 female. Age is released only as bands. In particular, the public young arm includes `17-20`, and the old arm includes `50-55` and `55-60`; the prose summaries (“approximately 20–40” and “approximately 60–80”) are not exact per-donor ages.
- These are intercostal-muscle samples from deceased transplant organ donors. `Donor_type` is DBD or DCD. The SDRF `injury` field records terminal injury/cause-of-death information, not a clean healthy/disease flag. It must not be reinterpreted as muscle disease or discarded as irrelevant without a prespecified rule.
- Supplementary Table 1 contains technical replicates and multiple assays per donor. The inferential unit remains the donor.

| Donor | Age band | Group | Sex | Type | SDRF injury field | scRNA SampleID(s) | snRNA SampleID(s) |
|---|---|---|---|---|---|---|---|
| 470BR | 17-20 | young | M | DBD | hypoxic brain damage | 5386STDY8090404, 5386STDY8090405 | 5386STDY8090406, 5386STDY8090407, 5386STDY8552709 |
| 398B | 20-25 | young | F | DBD | Trauma RTA - car | — | 5386STDY7796286, 5386STDY7796287 |
| 583B | 25-30 | young | M | DBD | Intracranial haemorrhage | WS_A_SKM9845427, WS_A_SKM9845428 | WS_A_SKM10691778 |
| 566B | 25-30 | young | F | DBD | Intracranial haemorrhage | — | WS_A_SKM9426499, WS_A_SKM9426500 |
| 411C | 25-30 | young | M | DCD | Hypoxic brain damage | 5386STDY7835292, 5386STDY7835293 | — |
| 464C | 35-40 | young | M | DCD | Hypoxic brain damage | 5386STDY8047212 | 5386STDY8552613, 5386STDY8552614 |
| 529C | 35-40 | young | M | DCD | Intracranial thrombosis | — | WS_A_SKM9426501, WS_A_SKM9426502 |
| 591C | 35-40 | young | M | DCD | Intracranial haemorrhage | WS_A_SKM9845429, WS_A_SKM9845430 | — |
| 621B | 50-55 | old | M | DBD | Intracranial haemorrhage | WS_A_SKM10490125, WS_A_SKM10490126 | — |
| 582C | 55-60 | old | F | DCD | hypoxic brain damage | WS_A_SKM9845425, WS_A_SKM9845426 | WS_A_SKM10691779 |
| 362C | 60-65 | old | M | DCD | Other | 5386STDY7600838, 5386STDY7600839 | 5386STDY7600836, 5386STDY7600837, WS_A_SKM9426496 |
| 502B | 60-65 | old | F | DBD | Intracranial haemorrhage | — | WS_A_SKM9426495 |
| 500C | 60-65 | old | M | DCD | Intracranial haemorrhage | — | 5386STDY8493510 |
| 343B | 70-75 | old | F | DBD | Intracranial haemorrhage | mus_SNuc7511880 | WS_A_SKM9426497, WS_A_SKM9426498 |
| 367C | 65-70 | old | M | DCD | Intracranial haemorrhage | 5386STDY7645353, 5386STDY7645354 | 5386STDY7645355 |
| 339C | 70-75 | old | F | DCD | Intracranial haemorrhage | mus_SNuc7468112 | mus_SNuc7468111 |
| 640C | 70-75 | old | F | DCD | Intracranial haemorrhage | WS_A_SKM10490127, WS_A_SKM10490128 | — |

The complete 43-row SampleID-to-ENA mapping is in `data/manifest/E-MTAB-13874.json`.

### Raw-count objects, sizes, and access

All atlas H5AD links returned HTTP 200 to HEAD requests with byte sizes and range support. The full human raw object was additionally inspected with 55 schema-only range requests totaling 3,544,752 bytes. It has `X` and `layers/counts` sparse matrices, both with float32 `data`. Thus a raw-count layer is present, but integer-valued counts are **not** established by dtype and were not tested by reading values.

| Object | Bytes | Public URL |
|---|---:|---|
| Human raw cells+nuclei H5AD | 4,025,846,616 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad |
| Human processed cells+nuclei H5AD | 2,022,042,678 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_human_pp_cells2nuclei_2023-06-22.h5ad |
| MuSC H5AD | 235,029,240 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_MuSC_human_2023-06-22.h5ad |
| Myonuclei-to-myofiber H5AD | 766,049,888 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_myonuclei2myofiber_human_2023-06-22.h5ad |
| Myonuclei H5AD | 649,019,946 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_myonuclei_human_2024-04-14.h5ad |
| Immune H5AD | 224,028,501 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_immune_human_2023-06-22.h5ad |
| Fibroblasts/Schwann H5AD | 325,448,911 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_fibroblasts_Schwann_human_2023-06-22.h5ad |
| Endothelium/SMC H5AD | 282,584,411 | https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/SKM_Endothelium_SMC_human_2023-06-22.h5ad |

BioStudies itself lists the 14,339-byte IDF and 132,328-byte SDRF and marks the study public (release 2024-03-07); raw sequencing is linked through ENA `ERP158366`. The 301,754-byte Supplementary Table 1 workbook was read in memory. No raw sequencing file was downloaded.

**Kedlian verdict:** metadata and public access pass. A donor-level young/old design is supported, subject to cell/donor eligibility checks and count-value validation. A continuous-age model using the checked public metadata is blocked. The 4.03 GB raw H5AD must use the resource guard if later downloaded; no download was made here.

## Perez 2022 — GSE167186

Primary records read were the [GEO text record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167186&targ=self&form=text&view=full), [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/soft/GSE167186_family.soft.gz), [supplementary-file directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/suppl/), [file inventory with exact sizes](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/suppl/filelist.txt), [simplified snRNA sample sheet](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/suppl/GSE167186_SimplifiedMetadataSheet.xlsx), [paper full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9792217/fullTextXML), and [primary supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9792217/supplementaryFiles). Only Supplementary Materials and Methods and Supplementary Tables 1–2 were extracted from the 2,388,235-byte archive.

### snRNA donor mapping and exact ages

GEO contains 17 snRNA GSMs (`GSM5098737`–`GSM5098753`), each with one donor-like HM label, exact age, Young/Old group, gene-by-barcode CSV, and filtered 10x H5. The simplified sheet independently gives the same GSM/HM/H5 mapping. There are 6 young and 11 old donors; no snRNA donor is labeled sarcopenic in GEO, so bulk frail labels must not be imported into this arm.

The filtered H5 barcode dimensions sum to **134,486 released filtered nuclei**. This total comes from the current GEO files, not a later secondary analysis. All 17 `/matrix/data` arrays have dtype int32 and `/matrix/barcodes` dimensions matching the second `/matrix/shape` dimension. The checks used 55 range requests totaling 3,589,111 bytes and read no count values.

| GSM | Donor label | Group | Exact age | Filtered nuclei | CSV bytes | H5 bytes |
|---|---|---|---:|---:|---:|---:|
| GSM5098737 | HM1 / OldHM1 | Old | 75 | 6,587 | 9,077,523 | 3,016,762 |
| GSM5098738 | HM2 / OldHM2 | Old | 76 | 7,626 | 11,950,749 | 3,130,748 |
| GSM5098739 | HM3 / OldHM3 | Old | 75 | 4,206 | 6,477,572 | 2,861,643 |
| GSM5098740 | HM4 / OldHM4 | Old | 79 | 11,512 | 9,426,872 | 2,957,839 |
| GSM5098741 | HM5 / YouHM5 | Young | 20 | 5,877 | 7,808,502 | 3,003,523 |
| GSM5098742 | HM6 / OldHM6 | Old | 82 | 7,311 | 7,372,139 | 3,155,099 |
| GSM5098743 | HM7 / OldHM7 | Old | 72 | 4,478 | 6,021,744 | 2,723,950 |
| GSM5098744 | HM8 / YouHM8 | Young | 21 | 5,510 | 5,034,109 | 2,855,478 |
| GSM5098745 | HM9 / OldHM9 | Old | 70 | 9,552 | 14,068,517 | 3,670,266 |
| GSM5098746 | HM10 / YouHM10 | Young | 21 | 4,739 | 7,820,931 | 2,999,165 |
| GSM5098747 | HM11 / OldHM11 | Old | 75 | 3,991 | 6,038,507 | 2,675,123 |
| GSM5098748 | HM12 / OldHM12 | Old | 72 | 5,342 | 8,464,293 | 3,147,965 |
| GSM5098749 | HM13 / YouHM13 | Young | 19 | 3,623 | 6,384,026 | 2,467,164 |
| GSM5098750 | HM15 / OldHM15 | Old | 79 | 3,537 | 4,670,888 | 2,343,927 |
| GSM5098751 | HM21 / YouHM21 | Young | 22 | 14,685 | 21,072,996 | 26,004,630 |
| GSM5098752 | HM22 / YouHM22 | Young | 21 | 13,829 | 15,836,519 | 19,284,317 |
| GSM5098753 | HM23 / OldHM23 | Old | 90 | 22,081 | 26,626,284 | 32,048,595 |
| **Total** | **17 donors** | **6 Young, 11 Old** | — | **134,486** | — | — |

The exact H5 URLs are recorded in `data/manifest/GSE167186.json`; each follows `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5098nnn/<GSM>/suppl/<GSM>_<HM>_filtered_feature_bc_matrix.h5`.

### Sex and the separate bulk arm

- The primary paper states that biopsies came from **72 adult men**. Therefore a mixed-sex analysis is not supported for the snRNA subset. However, GEO has no per-snRNA-GSM sex field and no public HM-to-bulk-X crosswalk; per-donor sex is a study-level source-cohort assignment, not an independently released per-sample field.
- The bulk arm is **72 distinct bulk RNA-seq samples**, `GSM5098661`–`GSM5098736`, on GPL20301/HiSeq 4000. It is not single-cell or single-nucleus evidence.
- The primary supplement divides the bulk cohort into 19 young, 29 old healthy, and 24 frail men. These clinical groups and exact bulk ages must remain separate from the 17-snRNA-donor table above.
- The 1,335,795-byte `GSE167186_counts.csv.gz` header contains `Symbol` plus the 72 bulk X sample labels. It is the bulk gene-count table, not an snRNA object.

### Which files contain snRNA raw counts

- Each `GSM…_HM…_filtered_feature_bc_matrix.h5` contains a filtered 10x sparse count matrix with int32 counts. These are the preferred verified raw-count inputs.
- Each corresponding `GSM…_HM….csv.gz` has a gene-by-barcode layout and also represents a per-donor nuclei count matrix. Only a small layout/header probe was made; these CSV matrices were not downloaded or compared value-by-value with H5.
- `GSE167186_RAW.tar` (292,526,080 bytes) bundles the 17 CSV and 17 H5 files. It was HEAD-checked but not downloaded.
- `GSE167186_SimplifiedMetadataSheet.xlsx` (10,124 bytes) is mapping metadata, not counts.

**Perez verdict:** donor mapping, exact ages, public raw-count H5 files, and the released filtered-nuclei total pass the metadata gate. The 72-man bulk arm is excluded from any single-cell inferential N. Later analysis still requires guarded acquisition, per-donor eligibility checks, and reconciliation of any downstream QC with the released 134,486-nuclei starting point.

## Donor overlap and replication status

The studies have different documented recruitment frameworks: Kedlian used UK CBTM deceased transplant organ donors and intercostal muscle; Perez used McMaster/community vastus-lateralis biopsies. Their released donor identifiers do not overlap. Nevertheless, anonymized public IDs alone do not prove individual-level disjointness.

More importantly, Kedlian explicitly states that the authors obtained and reprocessed the Perez quadriceps data to compare myonuclei populations between intercostal and leg muscle. That is documented **dataset reuse**, not a second independent Perez cohort. The 17 new Kedlian donors must therefore be distinguished from Perez-derived results in Kedlian figures and text.

**Decision:** no shared human donor is documented, but donor-level independence is not proven and Perez data reuse is certain. Do not label Perez an independent replication of Kedlian (or vice versa), do not pool their donors, and do not count Perez twice when summarizing evidence. A future replication label requires an explicit overlap adjudication; compatible reanalysis can otherwise be described only as a cross-cohort comparison.

## Not done and next gate

- No expression matrix, raw FASTQ, or raw-count archive was downloaded; no expression job was run.
- No biological outcomes, harmonized cell types, QC exclusions, differential expression, composition, transport, Hellinger, or replication estimates were examined.
- Kedlian float32 count values were not tested for integer-valuedness; Perez CSV and H5 matrices were not value-compared.
- Per-cell/nucleus donor eligibility and downstream QC remain unevaluated.
- Before analysis, obtain parent review, resolve the replication-label rule, freeze group/age and assay contrasts, and use sequential guarded downloads with the 30 GB disk reserve and documented memory limits.
