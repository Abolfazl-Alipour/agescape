# Block 9 mature-only interneuron refit

## Verdict

INCONCLUSIVE. MSSM was refit for donors age 30 years and older, but every subtype is UNTESTABLE for mature-only MSSM-to-HBCC replication because `/home/reza/projects/agescape-revision/scratch/hbcc.h5ad` is absent. Under the frozen rule, any subtype UNTESTABLE in either cohort makes the block INCONCLUSIVE. No cross-cohort null is claimed. No gene among FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16 can be confirmed as a surviving mature-only cross-cohort shared gene. All six are UNTESTABLE in HBCC for VIP, pvalb, and sst. In the supplemental MSSM-only WLS fit, GABRB1, IGF1R, and TRMT11 pass FDR in all three subtypes; FBN1 and PHF3 pass in pvalb and sst but not VIP; ZBTB16 passes only in sst. These MSSM-only patterns are not substitutes for replication.

## Frozen specification

The only change to the original donor selection is restriction to donors age 30 years and older in both cohorts. The source of this cutoff is `scripts/run_mature_only_composition.py`: Human mature = 30 years and older, Reza's maturity definition, based on brain development completing around age 30. Mouse maturity is not used.

The original `run_unit_de` implementation is reused for donor-pseudobulk voom-WLS and OLS. QC remains primary data, disease normal, at least 1500 raw UMIs per cell, at least 20 eligible cells per donor and cell type, and at least 10 donors per unit. Pseudobulk is summed raw counts. The response is log1p counts per 10,000. The design is age in years plus sex plus donor mean log UMI. Low-expression filtering, technical-gene masking, one-iteration voom precision weighting, and donor unit inference are unchanged. There is no downsampling.

The technical mask excludes uppercase MT-, RPL, RPS, MRPL, and MRPS prefixes, plus XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D. The original human model has no study, pool, or cohort indicator within either separately fit cohort. The shared-rule unadjusted sensitivity therefore removes depth and fits the original OLS function with age plus sex. It is labeled depth-unadjusted and does not replace the depth-adjusted primary.

The required replication statistics are MSSM-to-HBCC effect-size correlations, same-direction counts among MSSM FDR genes, both-FDR same-direction counts against hypergeometric expectation, and the moderated three-subtype shared-gene intersection. These require both cohorts and are UNTESTABLE in the mature-only arm. The original full-range rows are copied from the specified targets, not recomputed.

The frozen verdict rule is: PROVES only if at least one subtype retains significant MSSM-to-HBCC replication, with the effect-size correlation passing the original replication criterion and both-FDR same-direction count above hypergeometric expectation. DOES NOT PROVE requires no replication with adequate donor counts and per-subtype power reasoning. INCONCLUSIVE applies if donor counts are too small or any subtype is UNTESTABLE in either cohort. The last condition applies here.

## Families fixed before correction

The original report states: Benjamini-Hochberg correction was applied within each cohort and cell type. The original implementation further defines this as BH q within each cohort by cell type over the filtered tested genes. The mature WLS and depth-adjusted OLS families retain this definition separately for each model. The depth-unadjusted OLS sensitivity is BH corrected separately within each cohort by subtype and is never pooled with the primary family. No cross-cohort mature family was created after seeing results.

The fixed seed inherited from the original analysis is 20261108. The MSSM fits are deterministic. No replication permutation was run because HBCC is absent.

## Inputs and blockers

MSSM input: `/home/reza/workspace/rev2-mssm-null-scratch/mssm.h5ad`, present. Its four mature pseudobulk caches were built without copying the h5ad. Each backed h5ad batch was gated on at least 20 GiB available memory.

HBCC input: `/home/reza/projects/agescape-revision/scratch/hbcc.h5ad`, absent. The HBCC arm stopped at this path. No download or source substitution was attempted.

The MSSM source-cohort QC route retains 213 donors at age 20 and older and 209 donors at age 30 and older before cell-type eligibility. The mature source-cohort age range is 31 to 89, median 67.0. The original comparison target describes 212 MSSM and 88 HBCC adult donors. Across the union of tested MSSM interneuron units, 195 mature donors met the 20-cell floor in at least one unit, age 31 to 89, median 66.0. The table below uses original Rev2 per-unit eligible counts for like-for-like donor comparison.

## Donor retention

| Subtype | MSSM full N | MSSM mature N | Mature age range | Mature median | HBCC full N | HBCC mature N |
|---|---:|---:|---:|---:|---:|---:|
| VIP | 197 | 193 | 31 to 89 | 67.0 | 80 | UNTESTABLE |
| pvalb | 189 | 185 | 31 to 89 | 67.0 | 80 | UNTESTABLE |
| sst | 177 | 174 | 31 to 89 | 66.0 | 82 | UNTESTABLE |
| GABA | 155 | 153 | 31 to 89 | 66.0 | 71 | UNTESTABLE |

All MSSM units remain above the original 10-donor floor. No MSSM subtype is UNTESTABLE on donor count.

## Full-range and mature-only FDR counts

Full-range values are the specified original Phase 1 OLS comparison targets. Mature values are the required Rev2 WLS primary and its original depth-adjusted OLS arm. This distinction is retained rather than presenting different estimators as identical. The mature tested-gene column is the BH family size for each MSSM subtype and model.

| Subtype | MSSM full OLS FDR | MSSM mature tested genes | MSSM mature WLS FDR | MSSM mature OLS FDR | HBCC full OLS FDR | HBCC mature FDR |
|---|---:|---:|---:|---:|---:|---:|
| VIP | 294 | 1932 | 148 | 145 | 86 | UNTESTABLE |
| pvalb | 333 | 1873 | 133 | 147 | 79 | UNTESTABLE |
| sst | 123 | 1906 | 366 | 204 | 704 | UNTESTABLE |
| GABA | 5 | 1883 | 39 | 32 | 4 | UNTESTABLE |

Depth-unadjusted OLS sensitivity FDR counts, each using its own cohort-by-subtype family:

| Subtype | MSSM mature age plus sex OLS FDR |
|---|---:|
| VIP | 385 |
| pvalb | 353 |
| sst | 359 |
| GABA | 153 |

## Replication comparison

Full-range Pearson r and same-direction counts below are among original MSSM FDR genes. Both-FDR counts and expectations use all original shared tested genes. Every mature entry is UNTESTABLE because HBCC was not available.

| Subtype | Full r | Mature r | Full same direction | Mature same direction | Full both FDR | Full expected | Mature both FDR |
|---|---:|---:|---:|---:|---:|---:|---:|
| VIP | 0.8946 | UNTESTABLE | 272 | UNTESTABLE | 24 | 0.425 | UNTESTABLE |
| pvalb | 0.8562 | UNTESTABLE | 295 | UNTESTABLE | 20 | 0.442 | UNTESTABLE |
| sst | 0.9295 | UNTESTABLE | 118 | UNTESTABLE | 62 | 1.455 | UNTESTABLE |
| GABA | 0.7249 | UNTESTABLE | 4 | UNTESTABLE | 0 | 0.000 | UNTESTABLE |

For additional context, the original full-range Rev2 genome-wide statistics were:

| Subtype | Genome-wide r | Unmoderated both FDR | Hypergeometric expected | Block permutation p |
|---|---:|---:|---:|---:|
| VIP | 0.7217 | 76 | 23.547 | 0.0010 |
| pvalb | 0.6839 | 66 | 18.215 | 0.0010 |
| sst | 0.7614 | 249 | 121.543 | 0.0010 |

The original full-range moderated both-FDR counts were 77 VIP, 65 pvalb, and 239 sst. The original Phase 1 targets were 24 VIP, 20 pvalb, and 62 sst. These are different original analysis versions and are not merged.

## Six moderated shared genes

The full-range moderated intersection is exactly 6 genes: FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16. DYNC1I1 is not in this intersection. It appears only when the unmoderated file is used. The mature moderated intersection itself is UNTESTABLE. The table gives supplemental MSSM mature WLS FDR status beside the required HBCC status. NOT_FDR means tested but q >= 0.05. NOT_TESTED means absent after the unchanged expression filter.

| Gene | MSSM VIP | HBCC VIP | MSSM pvalb | HBCC pvalb | MSSM sst | HBCC sst |
|---|---|---|---|---|---|---|
| FBN1 | NOT_FDR | UNTESTABLE | PASS | UNTESTABLE | PASS | UNTESTABLE |
| GABRB1 | PASS | UNTESTABLE | PASS | UNTESTABLE | PASS | UNTESTABLE |
| IGF1R | PASS | UNTESTABLE | PASS | UNTESTABLE | PASS | UNTESTABLE |
| PHF3 | NOT_FDR | UNTESTABLE | PASS | UNTESTABLE | PASS | UNTESTABLE |
| TRMT11 | PASS | UNTESTABLE | PASS | UNTESTABLE | PASS | UNTESTABLE |
| ZBTB16 | NOT_FDR | UNTESTABLE | NOT_FDR | UNTESTABLE | PASS | UNTESTABLE |

## Power and interpretation

- VIP: MSSM retained 193 donors, above the original 10-donor floor. HBCC has no mature refit because the input file is absent. Cross-cohort correlation, overlap power, and MDE80 are therefore not estimable. This is UNTESTABLE, not a null.
- pvalb: MSSM retained 185 donors, above the original 10-donor floor. HBCC has no mature refit because the input file is absent. Cross-cohort correlation, overlap power, and MDE80 are therefore not estimable. This is UNTESTABLE, not a null.
- sst: MSSM retained 174 donors, above the original 10-donor floor. HBCC has no mature refit because the input file is absent. Cross-cohort correlation, overlap power, and MDE80 are therefore not estimable. This is UNTESTABLE, not a null.
- GABA: MSSM retained 153 donors, above the original 10-donor floor. HBCC has no mature refit because the input file is absent. Cross-cohort correlation, overlap power, and MDE80 are therefore not estimable. This is UNTESTABLE, not a null.

The high MSSM retention does not rescue cross-cohort replication. With no HBCC mature model there are zero eligible HBCC donor inputs for every subtype, no shared mature effect vectors, and no meaningful correlation or overlap MDE80. Accordingly, these results neither prove survival from maturity nor prove loss of the program. The MSSM-only FDR counts describe one cohort and cannot satisfy the two-cohort verdict rule.

## Outputs and reproduction

- `results/followup-11pt-block9-mature-interneuron.csv`: mature per-gene MSSM WLS, depth-adjusted OLS, and depth-unadjusted OLS results.
- `results/followup-11pt-block9-mature-replication.csv`: original columns plus `full_range`, `status`, and `status_reason`. Full-range reference rows are paired with UNTESTABLE mature rows.
- `scripts/run_followup_11pt_block9_mature_interneuron.py`: original model reuse, age 30 donor filter, memory guard, validation, and report generation.

Run with `/home/reza/miniconda3/envs/agescape/bin/python scripts/run_followup_11pt_block9_mature_interneuron.py`. Resume uses the existing pseudobulk caches. No manuscript or `report/` file was changed.
