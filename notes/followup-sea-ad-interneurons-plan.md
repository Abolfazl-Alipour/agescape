# Step 9 SEA-AD interneuron extension: frozen analysis plan

Frozen on 2026-11-14 before fitting any age-expression model. This is the only new Step 9 analysis. It does not rerun signature transport, composition, Hellinger displacement, or the existing 21-control-donor SEA-AD DLPFC analyses.

## Inputs and cohort gate

- Expression: Allen Institute public archived `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`, acquired as the released per-donor objects under `MTG/RNAseq/previous_objects/donor_objects/`. The raw-count layer is `UMIs`; normalized `X` is not used.
- Donor metadata: public `sea-ad_cohort_donor_metadata_072524.xlsx`, sheet `SEA-AD_Cohort_Metadata`, from the Allen SEA-AD data page.
- Metadata-only inspection found 84 donors. The primary cohort requires `Cognitive Status == No dementia`, checked consensus clinical Control, and `Overall AD neuropathological Change` in `Not AD` or `Low`. Fifteen donors pass before cell-type floors (age range 72-102 years). Dementia, Intermediate/High ADNC, missing diagnosis/pathology, and non-control consensus diagnoses are excluded.
- Use only released final nuclei with `Used in analysis == True`, passing released RNA/library preparation flags, no failed experiment component, and `Number of UMIs >= 1500`. The UMI layer must be nonnegative integer-valued and agree with released per-nucleus UMI totals. No new doublet cutoff is imposed on the already QC-filtered final-nuclei release.
- Annotation mapping is fixed at the released subclass level: `Vip` to VIP, `Pvalb` to pvalb, and `Sst` to sst. Supertypes are not pooled. Every mapped nucleus must also carry released Class `Neuronal: GABAergic`.
- A subtype is testable with at least 10 eligible donors, at least 20 QC-passing nuclei in each included donor-subtype, and at least five residual degrees of freedom in the primary design. Donors below the nucleus floor are excluded from that subtype only. All exclusions are tabulated.

## Primary gene model

- Sum raw UMIs per donor and subclass. This donor pseudobulk is the inferential unit; nuclei are not independent replicates.
- Filter genes separately within subtype using the existing Rev2 rule: remove `MT-`, `RPL/RPS`, `MRPL/MRPS`, and the same sex-linked technical genes as `run_rev2_de_interneuron.py`; require total pseudobulk counts >=10, detection in >=10% of donors, and mean CP10K >=1.
- Response is `log1p(CP10K)`, using the all-gene donor-subtype library denominator. The primary continuous-age effect is the slope per year on this scale; reports also give the slope per 10 years. Positive means higher expression with older age.
- Primary model is the existing Rev2 one-iteration voom-style precision-weighted least squares: `expression ~ age_years + sex_male + mean_log_umi + multiome_fraction + adnc_low`. `mean_log_umi` is the donor mean log nucleus UMI total within that subtype. `multiome_fraction` is the fraction of eligible nuclei from a released method containing `multiome` (case-insensitive); all other released methods must be recognized singleome 10x v3.1. `adnc_low` contrasts Low with Not AD. Constant covariates are omitted and recorded. If age is aliased with another term or the full nonconstant design is rank-deficient, the affected primary subtype is untestable rather than silently relabeled.
- Gene standard errors and 95% confidence intervals use the fitted t model and residual degrees of freedom, as in Rev2. Benjamini-Hochberg FDR is applied across all filtered genes within each SEA-AD subtype. The small, old cohort does not meet the general M4 age-SD >=12-year power floor; this is a limited replication, not a newly powered lifespan cohort.

## Frozen replication estimands and families

Existing effects and gene sets are read without refitting. The unmoderated Rev2 targeted sets come from `reviewer-rev2-interneuron-genes.csv`. Its all-gene `reviewer-rev2-interneuron-figure-data.csv` was not retained in either repository, so it cannot supply genome-wide effects without a prohibited non-SEA-AD expression job. Genome-wide effects instead come from the already-computed, verified local aggregate `../agescape/results/h3-human-neuronal-gene-replication-all.csv` (original donor-pseudobulk OLS, model `age_plus_sex_depth`). No MSSM/HBCC expression data are read or refit. MSSM is the PsychAD discovery cohort and HBCC is the independent replication cohort. The unmoderated, same-direction MSSM-and-HBCC FDR sets are fixed at 76 VIP, 66 pvalb, and 249 sst genes; the moderated sensitivity sets are not substituted.

1. **Genome-wide shared-gene family:** all unique gene symbols present in the retained original OLS aggregate and passing the SEA-AD subtype filter. Primary Pearson and secondary Spearman correlations compare original MSSM OLS and SEA-AD OLS age slopes on the same response scale. Correlations with SEA-AD primary WLS slopes are secondary cross-method estimates. This source-model distinction is reported and is not conflated with the Rev2 WLS correlation printed in the manuscript. Same-direction counts are descriptive; non-significance is not equivalence.
2. **PsychAD+HBCC both-FDR family:** the 76/66/249 same-direction both-FDR sets. Report SEA-AD coverage, same-direction fraction, same-direction nominal fraction, and same-direction SEA-AD-FDR fraction (the last is triple-cohort FDR, not a relabeled two-cohort result).
3. **PsychAD discovery family:** all MSSM-FDR genes in the corresponding subtype (202 VIP, 143 pvalb, 375 sst), with the both-FDR subset above and the disjoint `PsychAD discovery without HBCC both-FDR support` subset reported separately. This avoids treating discovery-only and independently supported genes as one family.
4. **Mouse-derived P1 family:** `P1_Clean` contains 139 mouse genes; the frozen one-to-one ortholog table maps these to 118 unique human symbols (76 up and 42 down), with no up/down overlap. Evaluate those 118 human symbols within each subtype and report the 21 unmapped mouse genes as unavailable rather than silently counting them. Report the full mapped P1 set, its overlap with the PsychAD discovery set, and the disjoint mapped P1 genes outside that subtype's PsychAD discovery set. "Mouse-derived" describes provenance only; it is not a claim of mouse-exclusive biology.

For every targeted set, denominators are genes detectable and testable in SEA-AD; missing and filtered-out genes remain explicit. Gene-level q values always come from the genome-wide SEA-AD subtype family above. No additional selective threshold is used to define a successful replication.

## Coexpression-block calibration

Use the existing Rev2 method: 30 or fewer dominant-PC coexpression blocks from SEA-AD donor-by-gene logCP10K pseudobulk, seed 20261108, and 1,000 permutations. Within each block, permute SEA-AD gene-level FDR status and effect sign while holding discovery gene-set membership and MSSM effects fixed. Calculate nulls for genome-wide Pearson correlation and for same-direction and same-direction-SEA-AD-FDR overlap in each targeted family. Empirical p is `(1 + null >= observed) / 1001`; report null mean and observed-minus-null. Separately, use 1,000 coexpression-block bootstrap resamples (seed 20261208) for percentile 95% intervals for genome-wide correlations and targeted same-direction fractions. Set-level block tests form separate three-subtype families for (a) genome-wide correlation and (b) PsychAD+HBCC both-FDR replication; report raw and BH-adjusted p across those three tests. Other set-level tests remain secondary and are not promoted after seeing results.

## Declared sensitivities

- Unweighted OLS with the same nonconstant primary design.
- Primary WLS without `adnc_low`, to expose rather than erase pathology adjustment effects.
- Singleome-only pseudobulk with `age + sex + mean_log_umi + adnc_low`, only if donor/nucleus/residual-df floors still pass.
- Not-AD-only pseudobulk with the same primary terms, only if floors pass; otherwise explicitly untestable.
- A 50-nuclei-per-donor-subtype floor, only if at least 10 donors and residual-df requirements still pass.
- Leave-one-donor-out refits for genome-wide correlation and the both-FDR sets, including sign-flip counts and slope ranges for covered both-FDR genes.

No alternative young reference, Hellinger analysis, signature score, cell-fraction model, or relaxed donor floor is permitted.

## Interpretation limits fixed in advance

This is a cross-sectional, middle-temporal-gyrus analysis in a narrow 72-102-year control range with at most 15 donors before subtype floors. It cannot establish equivalence, an absence of aging, an interneuron-exclusive mechanism, or independence from residual AD/co-pathology. Low-ADNC adjustment is not proof that AD-aging confounding is absent. Effects are compared with DLPFC PsychAD/HBCC results across regions and assays, so differences can reflect region, composition within broad subclasses, technical design, and limited power as well as biology.
