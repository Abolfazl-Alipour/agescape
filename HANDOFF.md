# Current Project Handoff

Last updated: 2026-11-02. This is the operating handoff for the scientific artifact. The public website is not the current development priority.

## Repository and current state

- Repository: https://github.com/Abolfazl-Alipour/agescape
- Branch: `master`
- Main manuscript: `report/manuscript-aging-mechanisms.md`
- Reviewer guide: `PROJECT_OVERVIEW.md`
- Reproduction commands: `README.md`
- Latest verification: 48 tests pass and all 146 app-data anchors pass. One Bayesian-only test is skipped in the core environment. All raw-data checksums pass.
- Controlled-access donor-level DLPFC data, credentials, and derived donor-level values are not committed.

## Read in this order

1. `README.md`: current findings and reproduction.
2. `PROJECT_OVERVIEW.md`: scientific design, results, and limitations.
3. `report/manuscript-aging-mechanisms.md`: current manuscript. The new vector analysis has not yet been added.
4. `report/extended-cross-species-aging.md`: integrated fly, worm, zebrafish, mouse, human, and NMR interpretation.
5. `report/cross-species-expression-vector.md`: direct whole-vector analysis in mouse, human, and NMR.
6. `report/human-expression-vector-interpretation.md`: depth, assay, power, gene drivers, and HBCC replication.
7. `report/cross-species-bayesian-vector.md`: Bayesian effects, prior sensitivity, and diagnostics.
8. `report/organ-cell-composition.md`: direct three-species cell-fraction analysis.
9. `report/transcriptional-disorder.md`: entropy and off-program definitions and frequentist results.
10. `report/transcriptional-disorder-bayesian.md`: matching frequentist and Bayesian global comparison.
11. `notes/extended-analysis-methods.md` and `notes/extended-dataset-inventory.md`: extended calculations, designs, and reproduction.
12. `notes/cross-species-expression-vector-methods.md` and `notes/cross-species-bayesian-vector-methods.md`: original vector methods.
13. `notes/transcriptional-disorder-methods.md`: detailed calculations, priors, diagnostics, and output map.
14. `PLAN.md`, `KICKOFF.md`, and `notes/crucibles/`: historical material, not current evidence or preregistration.

## Non-negotiable rules

1. Never invent data, statistics, citations, or availability claims. Quantitative claims must trace to committed result tables.
2. Do not treat non-significance as equivalence or biological stability.
3. Keep donor or animal replication distinct from cell-level measurements.
4. Never commit credentials, raw controlled-access data, or per-donor DLPFC values. Only approved aggregate results may be public.
5. Preserve the distinction between nominal and multiplicity-adjusted evidence, global and cell-type estimands, and exploratory and confirmatory analyses.
6. Do not describe transcriptomic entropy or off-program expression as direct evidence of DNA damage, chromatin restoration, or epigenetic causation.
7. Run tests and relevant verification gates after changing results or claims.
8. Do not deploy or modify the personal website unless Reza explicitly asks.

## Current scientific results

### Signature transport

The primary signature is `P1_Clean`, containing 139 genes: 89 up-regulated and 50 down-regulated. `P0_Raw257` is the 257-gene sensitivity score. P1 applies technical, biotype, ribosomal, mitochondrial, and ambient-RNA filters. It does not remove sex-dimorphic or immediate-early genes; those belong to P2 and P3.

One BH family covers 252 primary human cell-type tests. Survivors are MSSM L2/3 neurons, MSSM oligodendrocytes, and PsychAD EN_L3_5_IT_3 neurons. No OPC or microglial result survives. Glial transport is a candidate, not an established mechanism.

Primary PsychAD whole-cortex slopes are +0.000186 before adjustment, p = 0.3494, and -0.000209 after adjustment, p = 0.4004. Exact M9/M10 donor overlap is unresolved, so M10 is not independent replication.

### Organ cell fractions

The direct question was whether cell-type fractions within organs change with age. Eight of 13 adequately covered mouse organs and 3 of 11 human cohort datasets had at least one FDR-significant fraction change. NMR skin was inconclusive because five animals permit only ten age-label assignments. These are compositional associations, not evidence of fibrosis or organ destruction.

See `results/organ-composition-age-tests.csv`, `results/organ-composition-summary.csv`, and `notes/organ-composition-methods.md`.

### Transcriptional disorder

The third hypothesis was separated into within-cell Shannon entropy, off-program expression, and the distinct earlier measurement of cell-to-cell variability.

Cells were downsampled without replacement to 1,500 UMIs. Entropy was calculated per cell as `-sum(p_g * log(p_g))` after selected gene exclusions. Off-program genes were uncommon in young target cells but characteristic of another lung cell type. Young references were cross-fitted by leaving out the evaluated young animal. Animal means were the inferential units.

No frequentist cell-type analysis detected an FDR-significant primary increase. Entropy estimates were sensitive to balancing at 20 cells per animal. The frequentist global entropy effect was -0.079 standard deviations over 15 months, p = 0.699. The Bayesian posterior mean was -0.060, with 94% HDI -0.508 to +0.356 and posterior probability of an increase 0.395.

Off-program expression showed a weak positive average trend. The frequentist global effect was +0.244 standard deviations over 15 months, p = 0.278. The Bayesian posterior mean was +0.243, with 94% HDI -0.224 to +0.702 and posterior probability of an increase 0.841. Across sensitivity analyses, increase probabilities ranged from 0.731 to 0.857, and every interval included zero. No individual cell-type off-program result survived FDR correction.

All 18 Bayesian models converged with zero divergences, but all had influential-observation warnings. The primary off-program model also had imperfect calibration. Frequentist analysis remains the confirmatory reference; Bayesian hierarchical estimates are exploratory. The third hypothesis is not established.

### Whole-expression-vector displacement across species

The direct vector analysis uses no UMI downsampling. It compares each complete within-cell expression-proportion vector with a cross-fitted young reference using Hellinger distance. Raw-count Euclidean distance is a separate depth-sensitive endpoint. Tissues and cohorts are analyzed independently, and species summaries use standardized effects per 10 percent of maximum lifespan.

Mouse shows an overall positive association: effect +0.1028 SD, 95% CI +0.0678 to +0.1378, p = 0.00000398 across 62 tissue-by-cell-type units and 23 animals. It remains positive after UMI-depth adjustment and after excluding pooled identifiers. Five cell types pass FDR, all positive: bladder cells, kidney fenestrated cells, podocytes, limb-muscle mesenchymal stem cells, and lung fibroblasts.

Human shows no overall association: effect +0.0365 SD, 95% CI -0.0221 to +0.0951, p = 0.2218 across 135 cohort-by-cell-type units and 479 cohort-specific adult donors. None of seven cohort summaries is significant. Six cell types pass within-cohort FDR in the unadjusted analysis, but all six lose significance after UMI-depth adjustment. Age-associated depth has the sign needed to create positive MSSM neuronal and negative integrated-lung epithelial associations. Lung hits also remain null after assay adjustment.

Normalized donor-pseudobulk analysis gives a more specific result. No gene survives depth-adjusted FDR in the two lung epithelial units. In MSSM interneurons, 294 VIP, 333 pvalb, and 123 sst genes survive depth-adjusted FDR. These changes replicate in HBCC: 272 of 294 VIP, 295 of 333 pvalb, and 118 of 123 sst genes have the same direction, with 24, 20, and 62 passing FDR in both cohorts. Seven genes pass FDR in all three MSSM subtypes, and all seven have the same direction in HBCC. RNLS passes FDR in all three HBCC subtypes. This supports reproducible interneuron transcriptional aging, not a universal whole-vector disruption claim.

NMR shows a negative point estimate of -0.4744 SD per 10 percent of maximum lifespan, exact p = 0.2, across basal and corneum skin compartments in five animals. Both compartment effects are negative. This does not establish stability or a significant between-species difference. The simple hypothesis of positive mouse and human effects with a weaker NMR effect is not observed because the overall human effect is null.

Bayesian analysis is consistent with these qualifications. Mouse posterior means are +0.1495 unadjusted and +0.0782 depth-adjusted, with both 94% HDIs above zero. The human global unadjusted posterior is positive (+0.0893), but an alternative hierarchical model is near zero and depth adjustment gives +0.0133 with every prior-sensitivity interval crossing zero. NMR has posterior mean -0.1581 with 94% HDI -0.7210 to +0.4798. Human global R-hat values are 1.0115 to 1.0152, all models have influence warnings, and NMR calibration is poor. Bayesian estimates are exploratory.

See `report/cross-species-expression-vector.md`, `report/human-expression-vector-interpretation.md`, `report/cross-species-bayesian-vector.md`, and the corresponding `notes/` and `results/` files. The manuscript has not yet been updated with this analysis.

### Extended fly, worm, and zebrafish analyses

The extended analysis tests the same cell-fraction and whole-expression-vector questions in three new model systems. It supports cell-type-specific aging, not a universal increase in transcriptomic displacement. Full interpretation is in `report/extended-cross-species-aging.md`, with methods in `notes/extended-analysis-methods.md` and dataset decisions in `notes/extended-dataset-inventory.md`.

Fruit fly AFCA shows FDR-significant fraction changes in 36 of 162 tested cell types. All 36 retain direction in observed-record sensitivity analysis, 30 retain FDR, and 31 retain FDR after total-nuclei adjustment. Thirty-eight cell types pass depth-adjusted expression-vector FDR in both directions, but the global effect is -0.0467 SD over 40 days, p = 0.7958. Bayesian global HDIs cross zero at every tested prior.

Gao worm shows ARSC and spermatheca fraction changes after platform adjustment. Germline also passes FDR in the whole-animal model. With a day-1 expression reference, the platform-adjusted global vector effect is +1.6755 SD over 14 days, with platform-stratified permutation p = 0.0025. With the more balanced day-6 reference, the effect is +0.5319, permutation p = 0.0535, and no individual cell group passes FDR. Somatic day-6 analysis is not significant. The independent Roux atlas gives positive descriptive slopes in 32 of 33 cell groups but has technical rather than biological replication. The analyses support worm expression displacement, but its magnitude and somatic evidence are reference-dependent.

Zebrafish testis shows an immune-fraction increase of +0.0926 CLR units per month, q = 0.0203. It remains significant after total-cell and incomplete-spermatogenesis adjustment. No fine germ-cell fraction change passes FDR. Fine expression-vector hits do not survive exact age-permutation FDR, and Bayesian global HDIs cross zero. Immune-fraction change is stronger evidence than cell-intrinsic drift in this dataset.

The rat accession is verified as `GSE137869`, but one library per tissue, sex, and condition does not support animal-level significance tests. The primate cardiopulmonary accessions are `CRA002577` and `CRA002689`; they are not yet ingested. The 55-animal rhesus preprint still lacks a verified public accession. Genome-only and bulk long-lived species are not used as H2 or H3 evidence.

### Methylation clock

Predictive R-squared was 0.3169 in mouse liver and 0.6383 in the human multi-tissue evaluation. MAE was 0.36 years and 5.06 years, respectively. Mouse scoring used 724 of 816 probes. These were benchmark reproductions, not independent validations.

## Pipeline map

| Stage | Entry point | Main outputs |
|---|---|---|
| Mouse QC and atlas | `scripts/run_m1_pipeline.py` | `results/m1-*.csv` |
| Pseudobulk DE and ranking | `scripts/run_m2_pipeline.py`, `scripts/run_m2_audit.py` | `results/m2-*.csv` |
| Trajectory and TMS agreement | `scripts/run_m3_pipeline.py`, `scripts/run_m3_agreement.py` | `results/m3-*.csv` |
| Signature and heterogeneity | `scripts/run_m5_pipeline.py`, `scripts/run_m5_heterogeneity.py` | `results/m5-*.csv` |
| Mouse, human skin, and NMR scores | `scripts/run_m6_pipeline.py` | `results/m6-*.csv` |
| Public human cohorts | `scripts/run_m9_pipeline.py` | `results/m9-*.csv` |
| DLPFC lifespan atlas | `scripts/run_m10_pipeline.py` | `results/m10-*.csv` |
| Organ cell fractions | `scripts/run_organ_composition.py` | `results/organ-composition-*.csv` |
| Cell-to-cell variability | `scripts/run_h3_entropy.py`, `scripts/run_d2_noise_sensitivity.py` | `results/h3-entropy-mouse-lung.csv`, `results/d2-noise-sensitivity.csv` |
| Entropy and off-program expression | `scripts/run_h3_transcriptional_disorder.py` | `results/h3-transcriptional-disorder*.csv` |
| Bayesian hierarchical analysis | `scripts/run_h3_bayesian_disorder.py` | `results/h3-bayesian-*.csv`, `results/h3-frequentist-global.csv` |
| Mouse-lung direct vector analysis | `scripts/run_h3_expression_vector.py` | `results/h3-expression-vector*.csv` |
| Cross-species direct vector analysis | `scripts/run_cross_species_expression_vector.py` | `results/h3-cross-species-vector-*.csv` |
| Human depth and assay audits | `scripts/audit_human_expression_vector.py`, `scripts/audit_human_vector_assay.py` | `results/h3-human-*-audit.csv` |
| Human gene drivers and HBCC replication | `scripts/run_human_vector_gene_drivers.py`, `scripts/run_human_neuronal_gene_replication.py` | `results/h3-human-*gene*.csv` |
| Bayesian vector models | `scripts/run_cross_species_bayesian_vector.py`, `scripts/run_human_bayesian_global_vector.py` | `results/h3-*bayesian*.csv` |
| Extended fly, worm, and zebrafish analyses | `scripts/run_afca_aging.py`, `scripts/run_gao_worm_aging.py`, `scripts/run_roux_worm_aging.py`, `scripts/run_zebrafish_testis_aging.py` | `results/extended-*.csv` |

`scripts/run_all.py` covers the core pipeline, direct organ-fraction and transcriptional-disorder stages, and the mouse-lung vector analysis. M9, M10, the cross-species vector analysis, human follow-up analyses, extended atlas analyses, and Bayesian analyses have separate environments or inputs. Extended reproduction commands are in `notes/extended-analysis-methods.md`.

## Reproduction and verification

Core analysis:

```bash
conda create -n agescape --file environment.lock
conda activate agescape
python -m pip install -r requirements-supplemental.txt
python scripts/fetch_data.py
python scripts/run_all.py
pytest -v
python scripts/check_app_data.py
```

Cross-species whole-vector analysis:

```bash
python scripts/run_cross_species_expression_vector.py
python scripts/run_cross_species_expression_vector.py --summarize-only
```

The first command downloads public full-transcriptome human data and recomputes all species. The second recomputes statistics from committed public donor-level distances without downloading cells again.

Human follow-up:

```bash
python scripts/audit_human_expression_vector.py
python scripts/audit_human_vector_assay.py
python scripts/run_human_vector_gene_drivers.py
python scripts/run_human_neuronal_gene_replication.py
python scripts/summarize_human_neuronal_programs.py
```

Bayesian analysis:

```bash
conda env create -f environment-bayesian.yml
conda activate agescape-bayesian
python scripts/run_h3_bayesian_disorder.py
python scripts/run_cross_species_bayesian_vector.py
python scripts/run_human_bayesian_global_vector.py
```

The Bayesian stage uses `bayesian-workflow` revision `18ab0de9b63b5bb03caa96146235e8f1534d45d4`. Large NetCDF files and machine-specific diagnostic artifacts stay local under `results/h3-bayesian/` and `results/h3-cross-species-bayesian/` and are gitignored. Aggregate posterior and diagnostic tables are committed.

M9 downloads public CELLxGENE data. M10 requires separately authorized Synapse inputs under `data/external/dlpfc-atlas/`. Never commit Synapse credentials or controlled donor-level data.

## Remaining scientific work

- Increase animal-level replication for transcriptional-disorder tests.
- Add technical covariates and direct ambient-RNA estimates.
- Predeclare future effect definitions and sensitivity analyses before evaluating new cohorts.
- Use chromatin measurements and experimental perturbation to test the proposed DNA-repair mechanism; transcriptomics alone cannot do so.
- Resolve cross-species comparisons with harmonized tissues, orthologs, age coverage, and expression summaries.
- Obtain independent scientific review before presenting any exploratory off-program trend as a finding.

## Repository hygiene

Before committing, inspect `git status`, keep data and credentials out of the commit, run `git diff --check`, and rerun relevant gates. Historical audits are records of corrections, not current instructions.
