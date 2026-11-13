# Follow-up Step 3: Lake et al. 2023 kidney metadata/access gate

## Verdict

- **Full 45-reference-donor analysis: BLOCKED.** Public exact ages and a common single-cell/nucleus analysis mapping were not verified for all 45 reference donors. Do not impute ages or substitute the secondary 20–80-year summary.
- **Restricted public scCv3/snCv3 signature and donor-level CLR arms: READY to advance to count-matrix validation.** There are 18 scCv3 and 10 snCv3 healthy reference donors with public exact ages and openly downloadable processed count matrices. Analyze the assays separately, keep donors as the inferential units, and freeze contrasts, covariates, eligible cell types, and multiplicity before outcomes. Readiness here is an access/metadata verdict, not an analysis result.
- **Hellinger: BLOCKED pending a separate eligibility check.** The young-reference definition, per-donor cell floors, and young-reference cell floors have not been evaluated against the count matrices and reconciled metadata. The access gate does not establish Hellinger eligibility.
- **Raw-sequencing arm: BLOCKED by authorization.** KPMP raw sequencing requires a KPMP data use agreement, and HuBMAP raw sequencing requires authorized dbGaP access. No request or external contact was made.

The important distinction is that **processed count matrices do not require KPMP/NIDDK authorization**. Controlled access applies to raw sequencing and other restricted human data, not to the public GEO `Counts` RDS files used for the proposed reduced arm.

## Primary sources and repository records actually read

- Lake et al., *Nature* 619, 585–594 (2023), DOI [10.1038/s41586-023-05769-3](https://www.nature.com/articles/s41586-023-05769-3): abstract, methods, and complete data-availability statement.
- [GEO SuperSeries GSE183279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183279), the full-study processed-data accession.
- [GEO GSE183276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276), the **scCv3 subseries only**, and [GEO GSE183277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183277), the snCv3 subseries. Their public metadata and field-description files were read locally.
- Nature [Supplementary Tables 1–37 ZIP](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05769-3/MediaObjects/41586_2023_5769_MOESM6_ESM.zip); Tables 1 and 2 were read from the archive.
- Public GEO file inventories for [SNARE2/GSE183273](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183273/suppl/), [Slide-seq2/GSE183274](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183274/suppl/), [Visium/GSE183275](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183275/suppl/), and the [full SuperSeries file list](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183279/suppl/filelist.txt).
- [CELLxGENE collection API](https://api.cellxgene.cziscience.com/curation/v1/collections/bcb61471-2a44-4d00-a0af-ff085512674c), including dataset versions, donor IDs, development-stage labels, and `raw.X` declarations. Metadata-only CELLxGENE Census `latest` queries were also made on 2026-11-13 for `donor_id`, `development_stage`, `disease`, and `suspension_type`; no expression values were read.
- [KPMP compiled dataset collection, DOI 10.48698/3z31-8924](https://doi.org/10.48698/3z31-8924), and [dbGaP phs002249](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002249).

## Full-study data locations

The paper reports more than 400,000 cells/nuclei from 45 healthy reference donors and 48 diseased patients across KPMP, HuBMAP, and HCA. Its data-availability statement names **GSE183279**, not GSE183276, as the processed-data SuperSeries.

| Data | Verified location | Access |
|---|---|---|
| Processed snCv3, scCv3, SNARE2, Slide-seq2, and Visium | GEO SuperSeries `GSE183279`; subseries `GSE183277`, `GSE183276`, `GSE183273`, `GSE183274`, and `GSE183275`, respectively | Public |
| Processed scCv3 and snCv3 | CELLxGENE collection `bcb61471-2a44-4d00-a0af-ff085512674c`; public API declares counts in `raw.X` | Public; distinct processed release from GEO |
| KPMP raw snCv3, scCv3, and Visium sequencing | KPMP repository and compiled DOI `10.48698/3z31-8924` | Restricted; KPMP DUA required |
| HuBMAP raw snCv3, SNARE2, and Slide-seq sequencing | HuBMAP compiled collection and dbGaP `phs002249` | Restricted; NIH DAC authorization required |
| KPMP pilot-nephrectomy Slide-seq raw sequencing | GEO SuperSeries `GSE183279` | Public exception stated by the paper |
| CZI/HCA living-donor scCv3 raw sequencing | `GSE169285` | Named as publicly available in the paper |
| Additional Visium not in KPMP | `GSE171406` | Named as available in the paper |
| COVID-AKI raw sequencing | Washington University KTRC | Request-required; not contacted |

The paper also states that 3D-imaging raw data are openly downloadable. Those imaging data are not inputs to the planned single-cell signature, CLR, or Hellinger analyses.

### Public count-file check

Unauthenticated HTTP HEAD requests returned `200 OK` for all three RNA count files checked; no matrix bytes were downloaded:

- scCv3 `GSE183276_..._Counts_03282022.RDS.gz`: 1,409,085,514 bytes.
- snCv3 `GSE183277_..._Counts_03282022.RDS.gz`: 910,246,297 bytes.
- SNARE2 RNA `GSE183273_..._Counts_03282022.RDS.gz`: 263,011,754 bytes.

GEO labels these as processed data and explicitly says raw data are not provided on the scCv3/snCv3 series records. In this context, “raw data not provided” must not be misread as absence of the released count matrices.

## Healthy donor and sample mapping

Supplementary Table 1 contains 110 omic experiment rows from 35 reference and 36 disease participants. Supplementary Table 2 contains 105 spatial/imaging experiment rows from 11 reference and 30 disease participants. Their union is exactly **45 reference and 48 disease participants**, with no participant assigned to both reference and disease. Donor 3593 is the one reference donor shared between the two tables.

The 45 reference donors comprise:

- 31 represented in the public scCv3/snCv3 metadata: 18 scCv3 and 13 snCv3;
- four additional SNARE2-only reference donors: 3409, 3494, 3505, and 3510;
- ten additional spatial-only reference donors: 20-576, 21-015, 21-019, 21-020, 3579, 3611, IU-F52, IU-F59, IU-M32, and IU-M61.

SNARE2 donors 3499, 3504, and 3535 overlap the snCv3 reference set; spatial donor 3593 also overlaps snCv3. The supplements provide experiment, specimen, patient, condition, assay, sex, and repository fields needed for the full row-level crosswalk. These modalities must not be concatenated to manufacture a single-cell donor set.

### GEO scCv3/snCv3 reconciliation

- GSE183276 metadata contain 109,741 cells, 45 donors, and 49 libraries. The healthy subset is **18 donors, 20 specimens, and 20 libraries**; the disease subsets are 15 CKD donors/17 libraries and 12 AKI donors/12 libraries. Thus, the GEO overall-design phrase “20 healthy reference individuals” is inconsistent with the released patient field: 20 is the healthy **library/specimen count**, not the donor count.
- GSE183277 metadata contain 203,702 cells, 36 donors, and 44 libraries. The healthy subset is **13 donors, 18 specimens, and 21 libraries**; disease subsets are 13 CKD donors/13 libraries and 10 AKI donors/10 libraries.
- Multiple libraries from one donor are not independent inferential units. In particular, scCv3 donor PRE019 has three libraries, and snCv3 donor 3535 has two specimens and five libraries.
- GEO metadata have patient, specimen, library, condition, sex, region, assay, and annotation fields but no age field. Ages below come from the public CELLxGENE donor metadata, not from GEO.
- The CELLxGENE releases have 107,344 sc cells and 172,847 sn nuclei, so they are not cell-for-cell identical to the GEO metadata releases. The crosswalk below uses only unique donor-level age mappings; it does not assume cell-level equivalence between releases.

## Are exact donor ages openly available?

**Yes for 28 of the 31 public scCv3/snCv3 reference donors; no for the other three, which have decade-only labels.** Public exact ages were not verified for the 14 additional SNARE2-only/spatial-only reference donors in this gate. The public GEO metadata alone cannot supply these ages.

The CELLxGENE Census query returned one donor/development-stage/disease combination per donor in each dataset. Healthy status (`disease = normal`) and the following IDs were reconciled to GEO `condition.l1 = Ref` using released library/specimen identifiers. In particular, GEO `PRE018`, `PRE055`, and `PRE062` map to CELLxGENE `PRE018-1`, `PRE055-1`, and `PRE062-1`; GEO `PRE098` maps through specimen/library `PRE98sc` to CELLxGENE `PRE98sc`; GEO `1153EO1`-style IDs map to `Sample1153-EO1`-style IDs. These formatting differences do not create extra donors.

### scCv3: 18 reference donors, all with exact public ages

| GEO patient | CELLxGENE donor | Specimens | Libraries | Age (years) |
|---|---|---:|---:|---:|
| 1153EO1 | Sample1153-EO1 | 1 | 1 | 35 |
| 1153EO2 | Sample1153-EO2 | 1 | 1 | 35 |
| 1153EO3 | Sample1153-EO3 | 1 | 1 | 30 |
| 1157EO1 | Sample1157-EO1 | 1 | 1 | 49 |
| 1157EO2 | Sample1157-EO2 | 1 | 1 | 41 |
| 1157EO3 | Sample1157-EO3 | 1 | 1 | 47 |
| 1158EO1 | Sample1158-EO1 | 1 | 1 | 48 |
| 1158EO2 | Sample1158-EO2 | 1 | 1 | 56 |
| 1158EO3 | Sample1158-EO3 | 1 | 1 | 35 |
| 1162EO1 | Sample1162-EO1 | 1 | 1 | 39 |
| 1162EO2 | Sample1162-EO2 | 1 | 1 | 52 |
| PRE018 | PRE018-1 | 1 | 1 | 56 |
| PRE019 | PRE019 | 3 | 3 | 53 |
| PRE027 | PRE027 | 1 | 1 | 31 |
| PRE038 | PRE038 | 1 | 1 | 42 |
| PRE055 | PRE055-1 | 1 | 1 | 66 |
| PRE062 | PRE062-1 | 1 | 1 | 40 |
| PRE098 | PRE98sc | 1 | 1 | 57 |

### snCv3: 10 exact-age and three decade-only reference donors

| GEO patient / CELLxGENE donor | Specimens | Libraries | Public development stage |
|---|---:|---:|---|
| 18-142 | 1 | 1 | seventh decade; exact age unavailable |
| 18-162 | 1 | 1 | seventh decade; exact age unavailable |
| 18-312 | 1 | 1 | eighth decade; exact age unavailable |
| 3490 | 1 | 1 | 14 years |
| 3499 | 1 | 1 | 57 years |
| 3504 | 1 | 1 | 45 years |
| 3535 | 2 | 5 | 76 years |
| 3593 | 2 | 2 | 66 years |
| 3613 | 1 | 1 | 59 years |
| KRP446 | 2 | 2 | 29 years |
| KRP460 | 2 | 2 | 63 years |
| KRP461 | 1 | 1 | 54 years |
| KRP462 | 2 | 2 | 69 years |

The public 14-year value for donor 3490 is a direct reason not to impose a 20–80-year range from a secondary summary. Any adult-only eligibility rule must be declared before outcomes rather than silently deleting this donor. Decade labels must not be converted to numeric midpoints.

## Unresolved facts and required next checks

- Exact public ages for the 14 SNARE2-only/spatial-only reference donors remain unresolved, as do exact ages for snCv3 donors 18-142, 18-162, and 18-312.
- Integer/raw-count contents, feature sets, matrix dimensions, barcode concordance, and count-layer behavior remain unverified because expression matrices were deliberately not downloaded or opened.
- GEO-to-CELLxGENE age mappings are donor-level only. The releases' different cell totals preclude assuming identical cell membership without a separate barcode/QC reconciliation.
- A restricted scCv3/snCv3 analysis still needs a frozen analysis plan, explicit assay handling, eligible cell types, donor/cell floors, and a raw-integer matrix check. The 18- and 10-donor assay-specific sets must not be pooled casually to inflate inference.
- Hellinger remains untestable until its young-reference rule and both donor and young-reference cell floors pass. No non-significance or equivalence claim is available from this gate.

## Local metadata provenance and work not done

Small public metadata are under `/home/reza/projects/agescape-revision-data/lake_kidney/`; filenames, byte sizes, and SHA-256 hashes are recorded in `data/manifest/lake_kidney.json`. Total local download size was 36,690,103 bytes (about 35 MiB). The count-file requests were headers only.

No expression matrix was downloaded or analyzed; no signature, CLR, Hellinger, or aging model was run. No KPMP/NIDDK/dbGaP/KTRC contact or access request was made. No TS2 file was modified, no heavy job was started, and no commit or push was performed.
