# Ovary metadata and public-access gate

Base revision: `8cd58bdfc0921636ba2846f898fc6d0ab96d3ccb`. Verification date: 2026-11-14.

## Scope and verdicts

This gate covers **GSE202601 and GSE329417 only**. No expression matrix, RAW tar, FASTQ file, or large processed atlas was downloaded; no biological outcome or aging association was computed. Checks used GEO SOFT/file inventories, primary text, small public metadata files, and HTTP HEAD requests. Downloaded metadata and article XML are under `/home/reza/projects/agescape-revision-data/ovary-metadata/`. Exact URLs, hashes for downloaded metadata, all 64 GSE329417 sample mappings/files, and byte sizes are in `data/manifest/ovary_cohorts.json`.

| Cohort | Metadata/access verdict | Permitted future analysis | Blocked or conditional arms |
|---|---|---|---|
| Jin, GSE202601 | **Ready with conditions** | Eight-donor young versus reproductively-aged contrast only, after count/barcode/eligibility checks | No powered continuous-age analysis. Same-nucleus RNA–chromatin coupling is blocked; donor-matched exploratory linkage remains conditional on valid peak–gene links and cell-type harmonization. |
| Postmenopausal, GSE329417 | **Ready with conditions** | Within-postmenopausal exact-age associations after corrected-count, atlas-join, donor-eligibility, and covariate checks; descriptive QC otherwise | No premenopausal comparison and no pooling with Jin. The 64 tissue-piece samples are not 64 independent donors. |

## GSE202601 — Jin et al.

Primary sources: [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202601/soft/GSE202601_family.soft.gz), [exact file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202601/suppl/filelist.txt), [snRNA metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202601/suppl/GSE202601_human_ovary_snRNA-seq_metadata.txt.gz), [snATAC metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202601/suppl/GSE202601_human_ovary_snATAC-seq_metadata.txt.gz), and the [primary paper/full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11839473/). The peer-reviewed article is Jin et al., *Nature Aging* 5, 275–290 (2025), first published 2024-11-22, DOI `10.1038/s43587-024-00762-5`, PMID `39578560`. GEO became public on 2023-03-12.

### Eight-donor map, exact ages, and groups

GEO gives exact ages rather than only the paper's ranges. There is one snRNA and one snATAC GSM per donor. All 16 GSM records say female and ovary; the paper states that all eight samples had **normal ovarian histology**, with H&E performed on all eight (3–5 slides per sample). Tissues were fresh-frozen normal ovaries from BioIVT/Cureline, collected after sudden death or bilateral salpingo-oophorectomy and frozen within 16 hours post-mortem/post-surgery.

| Donor | snRNA GSM | snATAC GSM | GEO age | snRNA metadata label | Group |
|---|---|---|---:|---|---|
| ovary1 | GSM6125533 | GSM6125541 | 23 | ovary23 | young |
| ovary2 | GSM6125534 | GSM6125542 | 27 | ovary27 | young |
| ovary3 | GSM6125535 | GSM6125543 | 28 | ovary28 | young |
| ovary4 | GSM6125536 | GSM6125544 | 29 | ovary29 | young |
| ovary5 | GSM6125537 | GSM6125545 | 49 | **ovary47** | reproductively aged |
| ovary6 | GSM6125538 | GSM6125546 | 51 | ovary51 | reproductively aged |
| ovary7 | GSM6125539 | GSM6125547 | 52 | ovary52 | reproductively aged |
| ovary8 | GSM6125540 | GSM6125548 | 54 | ovary54 | reproductively aged |

The paper describes four young donors aged **23–29** and four reproductively-aged donors aged **49–54**. Those ranges and all group labels are consistent with GEO. However, the released snRNA nucleus metadata calls ovary5 `ovary47`, whereas GEO and the snATAC metadata call it 49. This is an unresolved exact-age discrepancy, not a license to choose 47 or 49 silently. It does not change the four-versus-four group contrast, but ovary5 must not contribute an exact numeric age until reconciled. No continuous-age analysis is planned.

### RNA/ATAC pairing and annotations

The modalities assay the **same donors and flash-frozen tissues, not the same nuclei**. RNA and ATAC have separate extraction protocols, library kits, GSMs, and BioSample accessions. Thus the word “multiome” in the study description must not be interpreted as same-nucleus paired RNA+ATAC measurements.

- snRNA metadata contain 42,568 rows and fields `age`, `group`, `celltype1`, and `celltype2`. The eight broad `celltype2` labels are BEC, EpiC, GC, IC, LEC, SC, SMC, and TC. The `_Y`/`_O` suffixes in `celltype1` are group-specific labels, not additional independent biological replicates or cell types.
- snATAC metadata contain 41,550 rows and fields `age`, `group`, and `celltype`. Its seven transferred labels are EC, EpiC, GC, IC, SC, SMC, and TC. EC is not directly interchangeable with the separate RNA BEC/LEC labels; harmonization is required before cross-modality analysis.
- The primary Methods describe Cicero co-accessibility links within 500 kb and ATAC gene activity from gene bodies plus 2 kb upstream. The GEO inventory has no standalone validated peak–gene link table. Any links used later need exact GRCh38 coordinates, gene identifiers, version, and cell-type mapping; inferred co-accessibility is not experimental or causal evidence.

### Counts and public objects

Each RNA donor has a public Cell Ranger filtered feature-barcode H5; each ATAC donor has a filtered peak-barcode H5 plus fragments and index. Exact per-file URLs and byte sizes are in the manifest. The [series RAW tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202601/suppl/GSE202601_RAW.tar) is **12,368,220,160 bytes** and was not downloaded. Count integrity, dimensions, and barcode joins remain unchecked.

| Public processed object | Bytes | Gate disposition |
|---|---:|---|
| `GSE202601_human_ovary_snRNA-seq_count.rds.gz` | 154,322,091 | Public; not downloaded or inspected |
| `GSE202601_human_ovary_snATAC-seq_count.rds.gz` | 1,072,649,963 | Public; not downloaded or inspected |
| `GSE202601_human_ovary_snRNA-seq_metadata.txt.gz` | 238,912 | Downloaded; schema/donor/group labels checked |
| `GSE202601_human_ovary_snATAC-seq_metadata.txt.gz` | 210,481 | Downloaded; schema/donor/group labels checked |

### Prespecified analysis boundary

1. **Primary:** reproductively-aged-minus-young donor-level effects, four donors per group, using each donor once per endpoint. No nuclei-as-replicates tests and no chronological-age slope.
2. **Before outcomes:** validate finite, nonnegative integer raw counts, unique matrix/metadata barcode joins, and released annotations. Require at least 20 annotated nuclei per donor-cell type and all four donors per group for the eight-donor contrast; otherwise mark that cell type descriptive/untestable. P1 requires detected up- and down-program genes. Hellinger additionally requires at least 50 cross-fitted young-reference nuclei and at least 20 passing nuclei per scored donor.
3. **Endpoints:** keep P1 transport, donor pseudobulk, CLR composition, and Hellinger displacement separate. Report effects with uncertainty and BH correction within each prespecified endpoint/cell-type family. Do not interpret non-significance in four-versus-four comparisons as equivalence.
4. **Sensitivity:** inspect procurement/PMI and depth metadata before modeling. Do not overfit covariates in this tiny design. Preserve the ovary5 age conflict explicitly.
5. **RNA–chromatin:** same-nucleus coupling is **blocked**. Only if valid links and harmonized cell types are later established may exploratory donor-matched pseudobulk RNA/peak associations be attempted. These are cross-sectional associations, **not chronological causal coupling**.

## GSE329417 — Byrne et al. postmenopausal cohort

Primary sources: [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329417/soft/GSE329417_family.soft.gz), [exact file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329417/suppl/filelist.txt), [28-donor metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329417/suppl/GSE329417_human_ovary_snucRNAseq_public_donor_metadata.tsv.gz), [64-sample processing log](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329417/suppl/GSE329417_human_ovary_snucRNAseq_10xFlex_sample_processing_metadata.xlsx), and the [primary preprint](https://www.biorxiv.org/content/10.64898/2026.05.11.724381v2.full). GEO became public on 2026-05-11.

**Publication status:** this is a bioRxiv preprint, DOI `10.64898/2026.05.11.724381`, first posted 2026-05-14; version 2 is available. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13193008/) labels it a preprint. [Crossref](https://api.crossref.org/works/10.64898/2026.05.11.724381) still classifies it as posted content/subtype preprint with no journal container, although it records an accepted date of 2026-11-01. No peer-reviewed journal publication was verified.

### Twenty-eight donors, 64 samples, and four GEO libraries

The public donor table provides **exact integer ages for all 28 donors**, spanning 50–84 years. The processing workbook has 64 rows and maps each tissue-piece sample to donor, probe barcode, plex, loading target, mass, and independence notes. All 64 rows matched one unique released CellBender H5AD filename. Donor-reported tissue-piece counts sum to 64.

| GEO library | Multiplexing |
|---|---|
| GSM9703623 | Plex 1, 16 demultiplexed samples |
| GSM9703624 | Plex 2, 16 demultiplexed samples |
| GSM9703625 | Plex 3, 16 demultiplexed samples |
| GSM9703626 | Plex 4, 16 demultiplexed samples |

| Donor | Age | Tissue pieces | Donor | Age | Tissue pieces |
|---|---:|---:|---|---:|---:|
| SN_0076 | 56 | 3 | SN_0153 | 84 | 1 |
| SN_0079 | 63 | 2 | SN_0155 | 58 | 2 |
| SN_0080 | 57 | 2 | SN_0161 | 74 | 2 |
| SN_0082 | 54 | 2 | SN_0163 | 50 | 2 |
| SN_0083 | 63 | 2 | SN_0164 | 54 | 2 |
| SN_0085 | 70 | 4 | SN_0168 | 54 | 3 |
| SN_0086 | 75 | 6 | SN_0169 | 71 | 2 |
| SN_0087 | 68 | 2 | SN_0170 | 62 | 1 |
| SN_0088 | 60 | 2 | SN_0238 | 67 | 2 |
| SN_0092 | 63 | 2 | SN_0244 | 65 | 2 |
| SN_0093 | 65 | 2 | SN_0245 | 75 | 2 |
| SN_0094 | 63 | 2 | SN_0249 | 75 | 4 |
| SN_0095 | 69 | 2 | SN_0253 | 61 | 2 |
| SN_0106 | 63 | 2 | SN_0263 | 54 | 2 |

The assay is multiplexed 10x Fixed RNA Profiling/Flex snucRNA-seq with 16 probe barcodes per plex, Cell Ranger multi 8.0.0, GRCh38-2024-A, and Element AVITI sequencing. Tissue pieces include cortex and medulla. Donors are women with postmenopausal ovaries; ovarian neoplasia and endometriosis were excluded, and a pathologist confirmed no premalignant/malignant ovarian pathology. Non-ovarian benign, premalignant, and malignant surgical indications remain in the public covariates and must not be erased by calling the cohort uniformly healthy.

**Independence caveat:** SN_0086 processing samples `3_4` and `16_4` shared one digested tissue/nuclei isolation aliquot. They are technical replicates and must be combined or selected before biological modeling. The workbook says all other samples used independent tissue chunks/dissociations. Repeated pieces from every donor still require donor-level aggregation or donor-aware models.

### Counts, annotations, and public objects

The 64 public per-sample H5ADs contain **unnormalized but CellBender-corrected UMI counts** in `X` and `layers["raw_counts"]`, restricted to Cell Ranger filtered features/barcodes. The final atlas uses the same corrected-count definition. These are not separately released original, uncorrected Cell Ranger raw matrices. The distinction and a corrected-count sensitivity analysis are mandatory; integer/value checks and atlas joins remain future gates.

The preprint reports 21 coarse cell types and 37 fine cell states. Final annotations are public in the atlas H5AD/SCE, with a separate immune-compartment object. Annotation fields and barcode uniqueness were not inspected because no expression object was downloaded.

| Public processed object | Bytes | Gate disposition |
|---|---:|---|
| `GSE329417_human_ovary_snucRNAseq_final_annotated_atlas.h5ad.gz` | 5,746,830,392 | Public; not downloaded |
| `GSE329417_human_ovary_snucRNAseq_final_annotated_atlas_SCE.rds.gz` | 6,009,522,227 | Public; not downloaded |
| `GSE329417_human_ovary_snucRNAseq_immune_compartment_annotated_scVI_milo.h5ad.gz` | 60,181,789 | Public; not downloaded |
| `GSE329417_human_ovary_snucRNAseq_public_donor_metadata.tsv.gz` | 809 | Downloaded; 28 rows checked |
| `GSE329417_human_ovary_snucRNAseq_10xFlex_sample_processing_metadata.xlsx` | 14,658 | Downloaded; 64 rows checked |

The [series RAW tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329417/suppl/GSE329417_RAW.tar) is **15,488,409,600 bytes** and was not downloaded. Per-sample file sizes and URLs are recorded individually in the manifest.

### Prespecified analysis boundary

1. **Primary:** within-postmenopausal associations with exact continuous age only, reported per 10 years with 95% intervals. No young-versus-old contrast, no premenopausal arm, and no pooling with Jin.
2. **Replication unit:** donor. Aggregate to one donor value per endpoint or use donor random effects/donor-clustered uncertainty for repeated pieces. Handle the SN_0086 shared-aliquot pair as technical replication. Never treat 64 samples or 439,011 nuclei as independent biological replicates.
3. **Before outcomes:** validate corrected-count integrity, sample/donor/atlas joins, and annotation coverage. Require at least 20 annotated nuclei per donor-cell type and at least 20 eligible donors with residual degrees of freedom for a parsimonious primary model; otherwise use descriptive QC only.
4. **Covariates/multiplicity:** freeze a small identifiable set before outcomes, considering plex and prespecified BMI, parity, and non-ovarian cancer state. Keep parity secondary and distinct from chronological age. Separate donor-pseudobulk expression, CLR composition, and any eligible vector endpoint; apply BH within each declared endpoint/cell-type family. Include a corrected-count sensitivity analysis.
5. **Fallback:** if future donor joins, counts, or annotations fail these checks, stop at descriptive annotation and QC. No association has been estimated in this gate.

## Not done and next decision

- No expression data, matrix values, differential expression, P1, composition, Hellinger, entropy, RNA–chromatin association, or other aging outcome was analyzed.
- No large object was acquired and no heavy job was run. Future expression work needs parent approval and sequential `scripts/resource_guard.py` execution with a 30 GB disk reserve and at most 8 GB analysis memory; use per-donor/sample files rather than either RAW tar.
- The Jin ovary5 47/49 discrepancy remains unresolved. Same-nucleus Jin coupling is not possible from the released design, and donor-matched linkage is not ready until valid peak–gene links are verified.
- GSE329417 supplies corrected rather than original uncorrected raw counts; this limitation remains explicit. Its publication is currently verified as a preprint, not a journal article.
- No testis file was modified, and no commit or push was performed. Parent review is required before any expression acquisition or analysis.
