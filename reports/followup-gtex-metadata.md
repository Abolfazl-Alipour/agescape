# Follow-up Step 10: GTEx v8 public metadata and access gate

Base revision: `2bb297e43f55b954d31ec1d93bf0dffc85a86eef`. Verification date: 2026-11-14.

## Scope and verdict

This was a **lightweight public-metadata/access check only**. No gene-count or TPM matrix bytes, matrix headers, range requests, raw sequences, protected phenotypes, or exact ages were downloaded or read. No aging outcome, expression summary, eligibility result based on expression, or biological claim was computed. No heavy job or access request was run.

The official GTEx records establish that the v8 gene-count and TPM objects, de-identified sample annotations, donor IDs, sex, and decade age bins are openly available. Public sample-to-donor mapping is complete for all 17,382 RNA-seq samples. All six requested tissues are **metadata-ready with conditions** for a future tissue-specific bulk ordinal age-bin analysis. None is expression-validated or authorized for analysis by this gate. Kidney has particularly sparse extreme age bins; ovary and testis are single-sex and have few 70–79 donors.

Exact ages remain protected. The [GTEx FAQ](https://gtexportal.org/home/faq) states that exact donor age is protected because of donor consent, public ages are 10-year brackets, and exact age requires dbGaP access. No age midpoint was substituted and no protected record was sought.

Downloaded public metadata and source records are under `/home/reza/projects/agescape-revision-data/gtex-metadata/`. Exact URLs, sizes, hashes, covariate details, per-tissue blockers, and the frozen design are in `data/manifest/gtex_v8.json`.

## Primary records and actual public files

The [GTEx bulk-expression download record](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression), [metadata record](https://www.gtexportal.org/home/downloads/adult-gtex/metadata), and official [open-access file-inventory API](https://gtexportal.org/api/v2/dataset/openAccessFilesMetadata?project_id=adult-gtex) were checked. The API identifies the files below as GTEx Analysis V8. Unauthenticated HTTP HEAD returned `200 OK` for both expression objects and gave exactly the inventory sizes; this verifies availability, not matrix contents.

| Public v8 input | Exact official URL | Bytes | What was verified |
|---|---|---:|---|
| Gene read counts, `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz` | https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz | 917,091,994 | Portal description “Gene read counts”; HEAD 200 and matching size; not downloaded |
| Gene TPM, `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` | https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | 1,625,807,147 | Portal description “Gene TPMs”; HEAD 200 and matching size; not downloaded |
| Sample attributes | https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | 11,512,258 | Downloaded; 22,951 rows, 63 columns; SHA-256 `74f6ab4c34ed2648d708a0ae6e6dff324f6c86ea723ae7d1c37d76f5221148f0` |
| Sample-attribute dictionary | https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx | 15,760 | Downloaded and read; field definitions |
| Subject phenotypes | https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt | 20,271 | Downloaded; 980 unique subjects; SHA-256 `821bdaff39e7a9a1d166919b3c786724c2b79c2861aeb936a2537a0f59b066f7` |
| Subject-phenotype dictionary | https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx | 10,605 | Downloaded and read; age/sex/Hardy definitions |

The expression objects have no local SHA-256 because they were deliberately not downloaded. Their object-header MD5 values, ETags, and modification timestamps are recorded in the manifest rather than represented as content hashes computed locally.

## Donor IDs, mapping, age bins, and sex

The analysis-freeze assay code is `SMAFRZE == "RNASEQ"`, not the prose spelling `RNA-Seq`. This selects **17,382 samples from 948 donors**. The subject table has 980 unique donors across assays, so its 980 must not be reported as the RNA-seq donor count.

For each `SAMPID`, the first two hyphen-delimited fields give the public donor ID: for example, `GTEX-1117F-0226-SM-5GZZ7` maps to `GTEX-1117F`. All 17,382 RNA-seq samples matched exactly one `SUBJID` in the public subject table; there were zero unmatched samples and IDs were unique within their respective tables. Within each requested tissue, sample count equals distinct donor count and there are no duplicate donors.

Public `AGE` has exactly six labels: `20-29`, `30-39`, `40-49`, `50-59`, `60-69`, and `70-79`. Across all 980 subject records their counts are 84, 78, 153, 315, 317, and 33, respectively. Age and sex are nonmissing in those public subject rows. `SEX` is coded `1=Male`, `2=Female`; the dictionary describes donor sex identification from self-report, family/next of kin, or medical-record abstraction.

## Tissue choice and public age spread

The set is motivated by completed or newly planned cohorts: whole blood for the PBMC arm, skeletal muscle for Kedlian/Perez, lung for De Man and existing integrated-lung work, kidney cortex for existing kidney cohorts and the restricted/blocked Lake arms, ovary for Jin/postmenopausal cohorts, and testis for Cui. These are **bulk gene-level contextual analyses**, not direct single-cell or cell-composition replications.

Counts below are distinct donors computed directly from the joined public metadata before any expression/QC exclusion. Each tissue has one RNA-seq sample per donor.

| GTEx tissue | Donors (= samples) | Male / female | 20–29 | 30–39 | 40–49 | 50–59 | 60–69 | 70–79 | Metadata readiness |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Whole Blood | 755 | 501 / 254 | 68 | 68 | 113 | 234 | 249 | 23 | Ready with conditions |
| Muscle - Skeletal | 803 | 543 / 260 | 67 | 65 | 124 | 255 | 264 | 28 | Ready with conditions |
| Lung | 578 | 395 / 183 | 33 | 40 | 93 | 200 | 190 | 22 | Ready with conditions |
| Kidney - Cortex | 85 | 66 / 19 | 1 | 8 | 8 | 29 | 35 | 4 | Limited power; sparse extremes |
| Ovary | 180 | 0 / 180 | 22 | 14 | 37 | 61 | 42 | 4 | Single-sex; sparse oldest bin |
| Testis | 361 | 361 / 0 | 36 | 31 | 52 | 117 | 116 | 9 | Single-sex; sparse oldest bin |

Readiness does not mean a gene model was run. Every tissue still requires GCT header/metadata reconciliation, count-integrity and gene-build checks, and model-rank/dispersion diagnostics before outcomes. Kidney's single 20–29 donor and four 70–79 donors make its trend potentially tail-driven. Ovary and testis require omission of the invariant sex coefficient. Testis bulk expression also cannot distinguish changing germ-cell-stage abundance from cell-intrinsic change.

## Usable public technical covariates

The sample dictionary was read, and actual field coverage/categories were tabulated for each requested tissue.

| Field | Meaning and observed usability | Frozen use |
|---|---|---|
| `SMRIN` | RNA Integrity Number; no missing values in the six tissues | Continuous primary covariate |
| `SMTSISCH` | Total sample ischemic time; 9 missing in Whole Blood, otherwise complete here | Continuous primary covariate in hours; complete-case |
| `SMEXNCRT` | RNA-SeQC exonic mapping fraction; complete in all six | Continuous primary covariate |
| `SMCENTER` | Collection site B1/C1/D1; one missing Muscle value | Categorical primary covariate |
| `SMNABTCHT` | Isolation protocol; one level in Whole Blood, two elsewhere | Categorical primary where estimable |
| `SMRDTTL` | Total filtered reads; complete | QC/depth sensitivity; library size is handled by the count-model offset |
| `SMGEBTCHT` | Expression batch type | Invariant `TruSeq.v1`; no coefficient |
| `SMRDLGTH` | Read length | 76 bp except five 101-bp samples; 76-bp-only sensitivity |
| `SMNABTCH`, `SMGEBTCH` | Isolation and expression batch IDs | Available but high-cardinality (73–344 and 68–208 levels per tissue); not primary fixed effects |
| `SMATSSCR` | Autolysis score | Missing in all 755 blood samples, 4 muscle and 5 lung; nonblood sensitivity only |
| `DTHHRDY` | Hardy death classification | Optional complete-case sensitivity, not a technical primary covariate |
| `SMGAPPCT`, `SM350NRM`, `SM550NRM`, `SMNUMGPS` | Selected released QC fields | Entirely missing in all six tissues; unusable |

Whole Blood needs a specific timing caveat: 329 released `SMTSISCH` values are negative, with observed range -1,287 to 1,401 minutes. The other five tissues have no negative values. The primary design preserves the released signed values rather than silently truncating them; a nonnegative-blood sensitivity is frozen. Nine missing blood values yield 746 primary core-covariate complete donors before expression gates. One missing muscle collection site yields 802; the other tissue complete-case counts currently equal their donor counts (578 lung, 85 kidney, 180 ovary, 361 testis).

## Access and redistribution terms

The [GTEx Portal Data License](https://gtexportal.org/home/license) was read from the official Portal application text. It states that downloading indicates acceptance of the terms, with no charges, usage fees, or royalties paid to the Portal. Users must acknowledge the Portal as the source using its requested wording and must not imply endorsement. Redistributors must maintain the current version or clearly disclose that their product does not reflect the most current/accurate Portal data. The license also contains warranty/liability disclaimers and notes that data formats may change. The full license, not this summary, governs reuse.

The [protected-access record](https://gtexportal.org/home/protectedDataAccess) states that raw sequences, expanded sample attributes, and expanded subject phenotypes including age are protected; access requires an approved dbGaP application and v8 protected data are hosted in AnVIL. The v8 study accession is `phs000424.v8.p2`. No authorization was requested, no controlled record was accessed, and Reza alone would decide and complete any future application.

## Frozen future design: bulk ordinal age-bin gene models only

This design is frozen before any GTEx expression or outcome access.

1. **Unit and tissues.** Fit separate models in Whole Blood, Muscle - Skeletal, Lung, Kidney - Cortex, Ovary, and Testis. One tissue RNA-seq sample is one donor measurement. Do not pool tissues or treat donors contributing to several tissues as independent in a cross-tissue test.
2. **Input gate.** Use the gene-read-count GCT, not TPM, for the primary model. Before examining age coefficients, verify unique sample columns, exact correspondence to the 17,382 public `RNASEQ` IDs, one selected sample per tissue donor, finite nonnegative integer counts, unique gene identifiers, and the released gene build. Stop affected tissues on mismatch.
3. **Fixed samples and genes.** Use complete cases for age, sex where estimable, RIN, ischemic time, exonic rate, collection site, and isolation protocol. Do not impute. Within each tissue, retain genes with raw count at least 10 in at least `ceil(0.10 × analysis samples)` samples; this filter uses counts only, not age associations.
4. **Primary model.** Per tissue and gene, fit a negative-binomial GLM with empirical-Bayes dispersion sharing and median-of-ratio size-factor offset:
   `expression ~ age_bin_order + SEX + SMRIN + SMTSISCH/60 + SMEXNCRT + categorical(SMCENTER) + categorical(SMNABTCHT) + offset(log size factor)`.
   Center continuous covariates within tissue for numerical stability. Omit `SEX` in ovary/testis and omit invariant protocol terms.
5. **Ordinal estimand.** Encode the six age bins as 1 through 6 in chronological order. Report the age coefficient as log2 fold-change per one-decade-bin step, 95% Wald interval, and adjusted p-value. This is not an exact-age slope in years, a longitudinal effect, or proof of stability when nonsignificant.
6. **Eligibility.** Require at least 50 complete-case donors, at least three occupied bins, positive residual degrees of freedom, full rank, and finite dispersion/model fit. Report sparse-level/leverage diagnostics. Failed models are untestable, not null.
7. **Primary multiplicity.** Apply Benjamini–Hochberg at q < 0.05 in **one family containing every retained gene-by-tissue primary age Wald test across all six tissues**. Do not replace this with within-tissue correction after seeing results.
8. **Secondary sensitivity families.** Separately correct all gene-by-tissue tests within each frozen arm: (a) categorical age with levels 20–29, 30–39, 40–49, 50–59, and pooled 60–79, joint likelihood-ratio test, at least 10 donors per level (kidney currently ineligible); (b) primary ordinal model excluding 70–79 donors; (c) extended complete-case autolysis/Hardy adjustment, with no blood autolysis model; (d) Whole Blood restricted to nonnegative ischemic time; and (e) 76-bp-only samples. Label all secondary and do not choose a favorable correction/model post hoc.
9. **Scientific boundary.** No cell fractions, per-cell drift, Hellinger displacement, entropy, single-cell composition, or deconvolution-as-measured-composition claims are included. Bulk associations cannot identify which cell type changed or distinguish abundance from within-cell regulation. No causal or direct single-cell replication claim is allowed.

Any future expression acquisition or analysis requires parent approval and sequential resource-guarded execution with a 30-GB disk reserve and at most 8 GB analysis memory. This metadata gate did not authorize or perform it.

## Not done and verification

- No expression matrix or outcomes were accessed; no exact ages were requested or imputed; no heavy job was run; no SEA-AD file was modified; no commit or push was performed.
- Public metadata computations and HTTP HEAD evidence are saved outside Git under `../agescape-revision-data/gtex-metadata/` (14 MiB total); source sizes and SHA-256 hashes are in the manifest.
- Independent reconciliation passed for JSON syntax, all four local annotation/dictionary sizes and SHA-256 hashes, both expression HEAD records, the complete 17,382-sample donor join, and every tissue's donor, sex, age-bin, complete-case, and batch counts.
- `git diff --check`: passed.
- Full suite with `/home/reza/miniconda3/envs/agescape/bin/python -m pytest -q`: **149 passed, 1 skipped, 1 warning in 10.40s**. The warning is the existing AnnData `ImplicitModificationWarning` in `tests/test_followup_muscle.py::test_kedlian_pipeline_synthetic_end_to_end`.
- The pre-existing untracked `notes/followup-sea-ad-interneurons-plan.md` remained untouched. This task created only `data/manifest/gtex_v8.json` and `reports/followup-gtex-metadata.md` in the repository.
