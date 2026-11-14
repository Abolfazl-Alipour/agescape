# Step 9: SEA-AD donor-pseudobulk interneuron extension

## Outcome

The SEA-AD MTG control cohort did **not** independently replicate the previously reported MSSM/HBCC interneuron aging pattern. No gene passed SEA-AD subtype-wide FDR, the frozen PsychAD+HBCC both-FDR sets showed no same-direction enrichment in gene-block permutations, and genome-wide correlations were near zero or negative depending on subtype and model. This is a limited, model-sensitive non-replication in 15 old donors, **not evidence of equivalence or no aging**. It does not overturn the original DLPFC result, but it does narrow its external support.

Only the new Step 9 donor-pseudobulk gene analysis was run. No P1 transport, signature, composition, Hellinger, alternative-reference, or GTEx expression analysis was performed, and no MSSM/HBCC expression was refitted.

## Verified inputs and exclusions

- Public release: Allen Institute SEA-AD MTG RNA-seq final nuclei, release `2024-02-13`, acquired as per-donor H5AD objects from the public archive. Source URLs, sizes, and SHA-256 hashes for all 15 objects and the donor workbook are in `data/manifest/sea-ad-mtg-2024-02-13.json`.
- Donor metadata: `sea-ad_cohort_donor_metadata_072524.xlsx`, sheet `SEA-AD_Cohort_Metadata`. The public source is the [Allen SEA-AD data page](https://brain-map.org/consortia/sea-ad/our-data). The associated atlas paper is Gabitto et al., *Nature Neuroscience* 27, 2366-2383 (2024), [doi:10.1038/s41593-024-01774-5](https://pmc.ncbi.nlm.nih.gov/articles/PMC11614693/).
- Of 84 donors, 42 had dementia, 10 no-dementia donors lacked checked consensus Control, 17 no-dementia clinical controls had Intermediate/High ADNC, and 15 no-dementia clinical controls had Not AD/Low ADNC and entered the cell-type gate.
- Exact age and sex were required to agree between each donor object and workbook. Diagnosis and ADNC were also checked in both. The included range is 72-102 years (donor SD 9.35), with 8 female and 7 male donors, 8 Not AD and 7 Low ADNC. This fails the project's general M4 age-SD >=12-year power floor and is explicitly a limited old-age replication cohort.
- Raw counts are the nonnegative integer `UMIs` layer; normalized non-integer `X` was ignored. For every analyzed donor-subtype, UMI row sums agreed with released `Number of UMIs`. Nuclei had to be in the released final analysis set, have passing preparation flags and no failed experiment component, and have >=1,500 UMIs. No new doublet threshold was imposed on the released QC-filtered object.
- Assays were verified as `10Xv3.1` singleome and `10xMulti` multiome. Five donors contributed multiome nuclei; the primary model adjusts for each donor's multiome fraction. Singleome-only sensitivity retains 15 donors.
- Released annotation mappings were exact: `Vip -> VIP`, `Pvalb -> pvalb`, and `Sst -> sst`. Every mapped nucleus also had Class `Neuronal: GABAergic`. Supertypes were not pooled.

All three subtypes passed the frozen floor of >=10 donors and >=20 nuclei per donor-subtype with 15 donors each:

| Subtype | Nuclei | Median per donor | Range per donor |
|---|---:|---:|---:|
| VIP | 18,699 | 1,222 | 515-2,505 |
| pvalb | 15,688 | 1,012 | 459-2,288 |
| sst | 9,914 | 427 | 241-1,495 |

The 50-nucleus sensitivity therefore uses the same 15 donors. The Not-AD-only sensitivity has eight donors and was declared untestable without relaxing the floor.

## Frozen models and multiplicity families

The pre-outcome plan and hash are in `notes/followup-sea-ad-interneurons-plan.md` and `data/manifest/sea_ad_interneuron_analysis_plan.json`. Donor pseudobulk was summed raw UMIs. The primary response was `log1p(CP10K)` and the primary model was Rev2 one-iteration voom-style WLS:

`expression ~ exact age + sex + mean log UMI + multiome fraction + Low ADNC indicator`

Every primary subtype had 15 donors and nine residual degrees of freedom. Gene effects, standard errors, t-based 95% CIs, and BH q values are in `results/followup-sea-ad-gene-effects.csv`. BH correction covers all filtered genes separately within each SEA-AD subtype. Positive slopes mean increased expression with age; the table includes effects per year and per 10 years.

Two pre-outcome source corrections are documented in the frozen plan:

1. The Rev2 all-gene figure table was not retained. Genome-wide comparisons therefore use the already-computed original MSSM `age_plus_sex_depth` OLS aggregate from the original project, with SEA-AD OLS as the like-for-like primary comparison and SEA-AD WLS as secondary. No MSSM/HBCC expression was read or refit. These are not the manuscript's Rev2 WLS correlations.
2. `P1_Clean` has 139 mouse genes but 118 unique mapped human symbols (76 up, 42 down) in the frozen one-to-one ortholog table. Twenty-one unmapped mouse genes are unavailable, not counted as failures.

The six targeted families are fully enumerated in `results/followup-sea-ad-replication-summary.csv`: PsychAD+HBCC both-FDR; all PsychAD discovery genes; PsychAD discovery without HBCC both-FDR support; all mapped mouse P1; mouse-P1/PsychAD overlap; and mouse P1 outside PsychAD discovery. The frozen unmoderated both-FDR sets have 76 VIP, 66 pvalb, and 249 sst genes. Gene-level q values always use the genome-wide SEA-AD subtype family. Primary set-level BH correction is across three subtypes separately for genome-wide correlation and both-FDR same-direction enrichment. Other set-level tests are secondary.

## Primary results

### Genome-wide effects

Block-bootstrap intervals are 95% percentile intervals over the SEA-AD coexpression blocks. Permutation p values are one-sided for positive replication, using 1,000 within-block permutations; 14 blocks were available from 15 donors.

| Subtype | Shared genes | OLS Pearson r (95% CI) | OLS Spearman rho | Block p | BH q | Secondary WLS Pearson r |
|---|---:|---:|---:|---:|---:|---:|
| VIP | 1,840 | 0.0054 (-0.0777 to 0.1003) | 0.0043 | 0.4166 | 0.6249 | -0.0356 |
| pvalb | 1,756 | -0.1775 (-0.3078 to -0.0440) | -0.1752 | 1.0000 | 1.0000 | -0.1917 |
| sst | 1,815 | 0.0390 (-0.0607 to 0.1179) | 0.0041 | 0.0969 | 0.2907 | -0.1674 |

There is no positive genome-wide replication. The negative pvalb point estimate does not establish a true negative cross-cohort relationship; it is model-sensitive and the one-sided replication test is not a test of a negative effect.

### Frozen PsychAD+HBCC both-FDR sets

The primary targeted estimand uses SEA-AD WLS slopes. "SEA-AD FDR" below means same-direction genes that also pass FDR in SEA-AD; these would be triple-cohort FDR genes, not a relabeling of the original two-cohort sets.

| Subtype | Tested/source | Same direction (fraction; block-bootstrap 95% CI) | Same-direction nominal | SEA-AD FDR | Block p | BH q |
|---|---:|---:|---:|---:|---:|---:|
| VIP | 71/76 | 29 (40.8%; 33.3-65.1%) | 4 | 0 | 0.8921 | 1.0000 |
| pvalb | 64/66 | 19 (29.7%; 14.3-32.9%) | 3 | 0 | 0.9620 | 1.0000 |
| sst | 240/249 | 86 (35.8%; 39.8-52.0%) | 9 | 0 | 1.0000 | 1.0000 |

No gene passed FDR among 1,923 VIP, 1,831 pvalb, or 1,890 sst genes tested genome-wide in SEA-AD. Missing or expression-filtered targeted genes remain in the denominators above as unavailable, not as opposite-direction evidence. The cluster bootstrap interval need not contain the observed fraction when block resampling is unstable; the intervals and observed values are reported without smoothing over that behavior.

### PsychAD discovery versus discovery-only support

| Subtype | All PsychAD discovery, tested/source | Same direction | Block p | Discovery without HBCC both-FDR, tested/source | Same direction | Block p |
|---|---:|---:|---:|---:|---:|---:|
| VIP | 183/202 | 95 (51.9%) | 0.2088 | 112/126 | 66 (58.9%) | 0.0240 |
| pvalb | 135/143 | 49 (36.3%) | 0.9500 | 71/77 | 30 (42.3%) | 0.7413 |
| sst | 354/375 | 148 (41.8%) | 1.0000 | 114/126 | 62 (54.4%) | 0.2318 |

The VIP discovery-only block p = 0.0240 is a **secondary, uncorrected set-level result** outside the two declared primary three-subtype correction families. It has zero same-direction SEA-AD-FDR genes and is not promoted as a robust replication. All PsychAD discovery subsets have zero same-direction SEA-AD-FDR genes.

### Mouse-derived P1 sets

Only 4 VIP, 3 pvalb, and 6 sst genes from the 118 mapped P1 symbols passed the SEA-AD expression filter. Same-direction counts were 2, 3, and 3, respectively, with no SEA-AD-FDR genes. PsychAD overlap was only 2 VIP, 1 pvalb, and 1 sst gene. These sparse sets cannot support a mouse-P1 replication verdict or any mouse-exclusive/interneuron-exclusive biological claim. Full and disjoint family counts are in the aggregate table.

## Sensitivities and donor influence

- Removing Low-ADNC adjustment left WLS genome-wide Pearson correlations at -0.0441 VIP, -0.1963 pvalb, and -0.1659 sst; both-FDR same-direction fractions were unchanged for VIP/pvalb and changed from 86 to 85 of 240 for sst. Pathology adjustment therefore does not explain the primary non-replication, but residual AD/co-pathology confounding remains possible.
- Singleome-only WLS correlations were -0.0216 VIP, -0.1634 pvalb, and -0.0703 sst. Same-direction counts were 29/71, 18/64, and 114/240. The sst change shows technical/model sensitivity rather than stable replication.
- OLS with the same primary design gave same-direction counts of 34/71 VIP, 16/64 pvalb, and 111/240 sst. Its genome-wide sst correlation was +0.0390 whereas WLS was -0.1674. This estimator dependence is a central limitation.
- The 50-nucleus floor gave identical results because every included donor-subtype already exceeded it.
- Not-AD-only models were untestable with eight donors, below the frozen ten-donor floor.
- Leave-one-donor-out WLS correlations ranged from -0.1416 to +0.0184 for VIP, -0.3007 to -0.0913 for pvalb, and -0.3339 to +0.2230 for sst. VIP had three sign flips and sst two. Both-FDR same-direction fractions ranged widely: 0.338-0.577 VIP, 0.203-0.391 pvalb, and 0.278-0.531 sst. Gene-level sign flips occurred in 34/71 VIP, 23/64 pvalb, and 136/240 sst covered both-FDR genes. These are not stable donor-level replications.

Per-donor expression and metadata rows remain only in external checkpoints. Public CSVs contain gene-level effects or cohort/subtype aggregates.

## What survived, died, or remains untestable

- **Survived:** source, diagnosis/pathology, exact-age/sex, raw-count, assay, annotation, QC, and donor-floor gates; complete reporting of null and model-sensitive results.
- **Did not replicate:** genome-wide positive correspondence, the 76/66/249 both-FDR gene programs, and SEA-AD gene-level FDR support in all three subtypes.
- **Untestable or inadequate:** Not-AD-only age models; a meaningful mouse-P1 set replication because only 3-6 genes were testable; any claim of no aging, equivalence, stability, AD-adjusted independence, or interneuron exclusivity.
- **Relationship to the manuscript:** SEA-AD does not add independent support to the DLPFC interneuron findings. It also does not invalidate those larger MSSM/HBCC analyses: region, age span, donor number, broad-subclass composition, assay, and estimator choice differ. The appropriate revision is to describe the interneuron program as cohort-specific with limited external generalizability, not universal human interneuron aging.
- **AD-aging boundary:** Low-ADNC adjustment and control selection cannot remove all AD or co-pathology. Age and residual pathology can remain entangled in this narrow old cohort; no causal or AD-independent aging claim is made.

## Not attempted

No Hellinger distance or young-reference relaxation; no alternative-reference workaround; no P1/signature/composition rerun; no 21-control-donor DLPFC rerun; no moderated-set substitution; no all-cell-type claim; no MSSM/HBCC expression job; no GTEx files or analyses touched; no commit or push.

## Reproduction and verification

From the repository at starting HEAD `2bb297e43f55b954d31ec1d93bf0dffc85a86eef`:

```bash
conda run -n agescape python scripts/resource_guard.py run \
  --path ../agescape-revision-data/sea-ad --min-free-gb 30 --memory-gb 8 -- \
  /home/reza/miniconda3/envs/agescape/bin/python \
  scripts/run_followup_sea_ad_interneurons.py --stage all
```

Acquisition is separately available through `--stage acquire`; it downloads only the public workbook and the 15 gated per-donor objects sequentially through `resource_guard.py`. No login or DUA was required.

Targeted tests: **6 passed**. Full suite: **155 passed, 1 skipped** in 9.47 seconds (the original 149 passing tests plus six Step 9 tests); one unrelated pandas `ImplicitModificationWarning` was emitted by the muscle test. Output SHA-256 hashes are in `data/manifest/sea_ad_interneuron_analysis_outputs.json`.

Reza review is required before proceeding to another step. No next-step analysis is authorized by this report.
