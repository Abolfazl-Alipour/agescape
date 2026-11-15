# Dataset 12: Human Subcutaneous Adipose Tissue Aging Access Check (Zhou et al., `HRA001279`)

## 1. Access-Check Disposition
- **Status**: **BLOCKED (Controlled Access, Access Check Only; Do Not Submit DAC Application)**.
- **Primary Accession**: [`HRA001279`](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001279) in NGDC/CNCB GSA-Human (plus corresponding-author request clause and code-only repository `https://github.com/TGXBio/huACA` / `colatutu/huACA`).
- **Primary Publication**: Zhou et al. (2023), *Signal Transduction and Targeted Therapy*, "High-resolution aging niche of human adipose tissues", `DOI: 10.1038/s41392-023-01315-9`, `PMID: 36918557`, `PMCID: PMC10014933`.
- **Observed Access Terms**: Anonymous HTTPS inspection of `https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001279` returned an automated `SecurityCheck` HTML gate (`Content-Length: 1014` bytes). GSA-Human controlled human sequencing datasets require registered investigator authentication, institutional verification, and formal Data Access Committee (DAC) approval before download. Per project rules, no DAC application was submitted and no gate bypass was attempted. Inspection of the linked GitHub repository (`TGXBio/huACA` / `colatutu/huACA`) confirmed it contains R scripts only (`ASAT_code.R`, `GSAT_code.R`, `scRNA_seq_functions.R`, `README.md`) and zero processed count matrices or Seurat `.Rdata` files.

## 2. Verification of Donor Counts, Exact Ages, Sex, and Assay
Exact donor metadata, previously listed as unresolved in [`followup-dataset-checklist.md`](../notes/followup-dataset-checklist.md) (#32), are now **fully verified** from Supplementary Table 1 (`41392_2023_1315_MOESM2_ESM.xlsx`) and Supplementary Methods (`41392_2023_1315_MOESM1_ESM.docx`) of `PMC10014933`:

- **Total donors**: **`14` healthy surgical subcutaneous adipose tissue (SAT) donors** across two cohorts (`Cohort 1`: `6` donors, `3` young / `3` old; `Cohort 2`: `8` donors, `4` young / `4` old), **`10` males (`M`) and `4` females (`F`)**.
- **Exact ages exist**: Yes, exact chronological integer years exist for all `14` donors in [verified donor CSV](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv), forming an extreme **bimodal young-vs-old design** with zero middle-aged donors (`30-67` yr empty):
  - **Young (`7` donors, ages `16-29` yr)**: Cohort 1 `YSAT-P1..3` (`22, 29, 28`), Cohort 2 `YSAT-P4..7` (`21, 16, 29, 27`).
  - **Old (`7` donors, ages `68-87` yr)**: Cohort 1 `OSAT-P1..3` (`71, 78, 74`), Cohort 2 `OSAT-P4..7` (`72, 87, 80, 68`).
- **Assay and anatomical confounding limits**: 10x Genomics Chromium scRNA-seq of stromal vascular fraction (SVF) cells (`Cell Ranger v2.1.0`, GRCh38, Illumina HiSeq X Ten) with `souporcell v2.4` genotype-based donor demultiplexing (`MOESM1`). Critically, **anatomical harvest site separates cohorts almost cleanly** (`Cohort 1`: `5` Flank + `1` Abdomen; `Cohort 2`: `8/8` Thigh), and SVF dissociation excludes mature adipocytes.

| Cohort | Individual | Group | Age (yr) | BMI (`kg/m^2`) | Sex | Anatomical Site |
|---|---|---|---:|---:|---|---|
| Cohort 1 | `YSAT-P1` | Young | [22](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [23.62](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Flank |
| Cohort 1 | `YSAT-P2` | Young | [29](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [25.25](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Flank |
| Cohort 1 | `YSAT-P3` | Young | [28](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [20.40](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | F | Abdomen |
| Cohort 1 | `OSAT-P1` | Old | [71](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [30.00](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Flank |
| Cohort 1 | `OSAT-P2` | Old | [78](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [20.20](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Flank |
| Cohort 1 | `OSAT-P3` | Old | [74](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [25.88](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | F | Flank |
| Cohort 2 | `YSAT-P4` | Young | [21](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [18.59](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Thigh |
| Cohort 2 | `YSAT-P5` | Young | [16](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [18.75](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Thigh |
| Cohort 2 | `YSAT-P6` | Young | [29](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [27.78](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Thigh |
| Cohort 2 | `YSAT-P7` | Young | [27](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [30.00](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | F | Thigh |
| Cohort 2 | `OSAT-P4` | Old | [72](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [24.73](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Thigh |
| Cohort 2 | `OSAT-P5` | Old | [87](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [20.52](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Thigh |
| Cohort 2 | `OSAT-P6` | Old | [80](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [20.75](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | F | Thigh |
| Cohort 2 | `OSAT-P7` | Old | [68](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | [25.23](../results/followup-organ-12-adipose-zhou-hra001279-donors.csv) | M | Thigh |

## 3. Analyses Run and Next Steps
- **None run**: No controlled-access download or expression analysis performed. If Reza decides to pursue `HRA001279`, formal GSA-Human DAC authorization is required, after which `Cohort 1` (`Flank/Abdomen`) and `Cohort 2` (`Thigh`) should be analyzed with cohort/anatomical-site stratification under binary `old_vs_young`.
