# Agescape follow-up execution plan

## Status and boundary

The manuscript review copy is complete. No follow-up dataset analysis starts until Reza finishes reviewing it. The work below is a sequential plan, not a claim that any new result exists.

## Rules for every step

1. Verify the accession, release, file inventory, raw-count layer, and sample-to-donor mapping before expression analysis. Catalog values are leads, not verified inputs.
2. Freeze the primary contrast, eligible cell types, reference ages, covariates, multiplicity family, and sensitivity analyses before examining outcomes.
3. Use donors, animals, fish, or independently measured pools as inferential units, never cells as independent replicates. Use continuous age only where numeric ages exist; otherwise use documented groups.
4. Keep composition, Hellinger displacement, Shannon entropy, off-program expression, and signature transport as separate endpoints. Do not treat non-significance as equivalence.
5. Run one memory-intensive job at a time with a process memory cap and no competing downloads. Stream backed H5AD reads, use small donor batches, and checkpoint aggregate tables. Stop rather than risk another machine-level memory failure.
6. Preserve study identity. Do not concatenate cohorts to manufacture donor numbers or cross-study p-values. Joint age permutations must follow each donor across all contributed cell types and organs.
7. Commit and push each completed metadata gate, analysis, and report with a clear message. Every reported number must trace to a committed aggregate result file.
8. Do not contact external services or request access directly. Where authorization is required, Reza decides and acts. Record blocked arms instead of bypassing them.
9. After each step, Reza reviews the result before the next step starts. Record explicitly what was done, what was untestable, and what was not attempted.

## Step 1: De Man human lung aging and entropy replication

### Metadata gate

- Use GSE281219 for the 32 newly profiled Baylor lungs and GSE136831 for the 28 IPF Cell Atlas control subjects. The published 60-donor analysis cannot be reproduced from GSE281219 alone.
- Build a crosswalk of accession, GEO sample, source, subject, age, sex, smoking, lobe/site, assay, barcode, released annotation, and QC disposition. Verify one age and sex per subject and raw integer counts.
- Reconcile the released QC and annotations. Require at least 20 cells per donor-cell type and six eligible donors per model, with adequate age spread. Require at least 50 remaining young-reference cells for Hellinger.
- Analyze Baylor and IPF controls separately for Agescape inference. Any integrated model must include source and remain identifiable.

### Analysis

- Reproduce the published KL-based entropy method for AT1, AT2, gCap, and aerocyte: released cell caps, age groups, zero-gene removal, thinning, and square-root transformation. Reproduce estimates, but use donor means and donor-level continuous-age models for inference. Cell-level Wilcoxon tests are descriptive only.
- Separately compute direct within-cell Shannon entropy after the Agescape technical-gene exclusions. Do not use the existing `run_h3_entropy.py`, which implements dispersion rather than Shannon entropy.
- Primary donor model per cell type: entropy ~ age + sex + mean log actual UMI. Add cell-type-intercept global models with donor-clustered errors and BH correction within the new external-validation family.
- Include fixed-depth, exact-multinomial, realized-binomial-thinning, no-square-root, full-depth, stringent-exclusion, donor-balanced, and repeated-seed sensitivity arms. Record realized UMIs because binomial thinning does not guarantee exact realized depth.
- Secondary endpoints if eligible: P1 transport with an expression-matched random null, cross-fitted Hellinger, CLR composition, and off-program expression. Do not merge entropy with off-program expression or Hellinger.

### Deliverables and stop rules

- Planned outputs: `data/manifest/GSE281219.json`, `GSE136831_controls.json`, `scripts/run_followup_gse281219.py`, metadata/eligibility tables, donor endpoints, model tables, robustness tables, and `reports/followup-gse281219.md`.
- Stop affected tests if subject linkage or raw counts are missing, QC cannot be reconciled, source/site is inseparable from age, P1 lacks detected up or down genes, or donor/reference floors fail. Mark them untestable rather than reporting a null.

## Step 2: Tabula Sapiens 2.0 within-person organ and skin analyses

### Metadata gate

- Use processed public counts and metadata from GSE306755 and the corresponding release, not the old four-donor Tabula Sapiens 1.0 skin object. Raw reads and their transfer agreement are unnecessary.
- Build the complete TSP donor/sample/organ crosswalk: age, sex, disease, assay, chemistry, anatomical site, enrichment, cell annotation, raw counts, and primary-data flags. Remove duplicate cells across organ files.
- Use 10x droplet data as primary and analyze Smart-seq separately. Verify the 24-donor age distribution, including only six donors over 60, before choosing models.

### Analysis

- Fit organ-by-cell-type P1, Hellinger, and secondary CLR models with continuous age, sex, depth, and identifiable assay/site covariates. Use donor-clustered standard errors and organ-interaction tests.
- For paired organ contrasts, restrict to shared donors and use donor-block bootstrap or joint permutation. One age assignment follows each donor across every organ. Do not use donor fixed effects in the chronological-age model because age is constant within donor.
- Run a skin-specific P1, Hellinger, and direct Shannon arm after checking anatomy and enrichment. Do not call Tabula 1.0 to 2.0 an independent replication because the releases overlap.
- Do not describe this as the manuscript's proteomic organ-age test. Tabula Sapiens 2.0 supplies transcriptomes, not plasma-proteomic organ ages.
- Optional RNA organ-age models must be exploratory, nested leave-one-donor-out, tuned only in inner folds, and excluded if they do not beat an age-mean baseline. They are not substitutes for plasma-proteomic ages.

### Deliverables and stop rules

- Planned outputs: `data/manifest/tabula_sapiens_2_0.json`, `scripts/run_followup_ts2.py`, donor/organ eligibility, organ-interaction models, skin models, robustness tables, and `reports/followup-tabula-sapiens-2.md`.
- Stop confirmatory paired-organ claims when fewer than six donors share the required organs/cell types. Report within-stratum estimates or untestable when age is inseparable from assay, anatomy, release, or enrichment.

## Step 3: Lake kidney

- Verify the full Lake dataset, healthy/control subset, participant mapping, and numeric age fields in primary metadata. Do not treat GSE183276 as the full-study accession or rely on the secondary 20-80 summary without checking it.
- Raw-data authorization is the bottleneck. Reza must decide and complete any KPMP/NIDDK request; the agent will not contact the repository.
- If authorized metadata are available, run signature scoring and donor-level CLR composition, then Hellinger only where young-reference and donor floors pass. Keep this separate from the four existing kidney cohorts for replication comparison.
- If access or age metadata remain unresolved, stop at the documented blocker and proceed to Step 4.

## Step 4: Muscle, Kedlian then Perez

- Kedlian (E-MTAB-13874): verify donor age and sample mapping, then use continuous age if numeric or the documented young/aged groups otherwise. Run the appropriate group contrast with 17 donors rather than an overconfident continuous-age model.
- Perez (GSE167186): verify exact ages and the released filtered-nuclei count. Use the count from the downloaded version, not a conflicting secondary total. Remember the bulk cohort is male-only and the snRNA arm is small.
- Analyze the studies independently and compare only harmonized cell types and compatible estimands. Check donor overlap before calling one a replication of the other.
- Deliver a joint muscle replication report with separate cohort tables and any inconsistent results shown rather than pooled away.

## Step 5: PBMC, open GEO cohorts only

- Start with GSE157007, GSE213516, and GSE214546. Do not use GSA or Synapse cohorts without separate authorization decisions.
- Verify exact age, sex, donor mapping, and raw counts. Analyze each cohort separately using continuous age where numeric and group contrasts where only labels exist.
- Replicate common immune-cell effects across studies without concatenating matrices or pooling study-specific p-values. State sex restrictions by cohort.
- Planned output: a three-cohort PBMC replication table with effect sizes, confidence intervals, and a clear account of any cell types that cannot be harmonized.

## Step 6: Cui testis

- Resolve the 25-GEO-sample versus 35-donor mapping before modeling. The maximum inferential N is 25 independent libraries/pools unless independent donor-level measurements are demonstrated.
- Use numeric age only where each pool has a defensible value; otherwise use pool-level age groups. Do not assign 35 individual ages as 35 independent observations.
- Put the biological caveat next to every result: changing germ-cell-stage composition can create broad expression displacement without cell-intrinsic aging drift. Analyze comparable germ-cell stages separately and model incomplete spermatogenesis where documented.
- If the mapping cannot be resolved from public supplements, stop at descriptive summaries and ask Reza for the mapping rather than guessing.

## Step 7: Ovary, Jin then postmenopausal cohort

- Jin (GSE202601): verify eight-donor mapping and exact ages. Treat it as a young/aged group contrast, not a powered continuous-age study. Multiome coupling can be exploratory and is distinct from chronological aging.
- Postmenopausal cohort (GSE329417): confirm demultiplexed per-donor metadata and exact ages. Analyze only within-postmenopausal age associations if the 28-donor age spread supports them. Do not create a premenopausal comparison or pool it with Jin to manufacture one.
- If exact ages are unavailable, restrict the second cohort to descriptive annotation and QC.

## Step 8: Jia heart, conditional and blocked

- Do not acquire or analyze Jia data until its accession is verified in the paper's data-availability statement.
- If verified, choose continuous or group analysis from actual donor metadata and repeat the standard eligibility gate.
- If the accession remains unresolved, mark this arm blocked and proceed. Do not substitute an unverified heart dataset.

## Step 9: SEA-AD interneuron extension

- Do not duplicate the existing 21-control-donor signature and composition analyses.
- A new donor-pseudobulk gene-level analysis is permissible because it does not require a young reference. Include diagnosis/pathology filters and assay covariates where identifiable.
- SEA-AD has zero eligible vector units under the existing young-reference rule. Do not silently relax that rule. Any alternative-reference Hellinger analysis requires Reza's explicit approval and a separate sensitivity label.

## Step 10: GTEx bulk, last

- Keep GTEx outside the single-cell H2/H3 evidence chain. Bulk data cannot estimate observed cell-type fractions or per-cell displacement.
- Public decade bins support group or ordinal-bin analyses only. Exact continuous-age modeling requires authorized exact ages; Reza must obtain any required dbGaP authorization.
- Run only tissue-specific bulk gene-level models with declared technical covariates and multiplicity control. Do not present deconvolution as measured single-cell composition.

## Explicitly not planned

- No unverified heart, liver, SenNet, gut, retinal, HRCA, or superior-frontal-gyrus accession will be used without primary-source verification.
- No GSA-Human data will be accessed without terms and authorization being resolved.
- Tiny skin, adipose, and other C-grade datasets will not be used to manufacture a general human conclusion.
- No longitudinal claim will be made: every planned analysis is cross-sectional.
- No full-pipeline rerun will overwrite established mouse, human, or naked-mole-rat outputs.

## Completion report for each step

Each step ends with a short report containing: verified inputs and exclusions; primary and sensitivity results with effect sizes and intervals; what survived; what died or was untestable; what contradicted the manuscript; what was not done and why; commit hash and test output; and the decision required before the next step.
