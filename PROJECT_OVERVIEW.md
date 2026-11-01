# Project Overview: Conserved Aging Explorer

**Subtitle:** within mouse; cross-species conservation not established

## Scientific question

Agescape asks when a mouse-derived aging-associated expression signature shows comparable age associations in other cohorts and species, and when atlas design prevents a defensible answer. It separates expression within cell types from cell-type abundance, UMI sampling noise, and DNA methylation clock performance.

The project is cross-sectional and observational. It is not a perturbation study, a longitudinal study, proof that naked mole-rats do not age, or evidence that a candidate glial association is a causal mechanism.

## Data and analysis design

| Arm | Data | Inferential unit | Purpose |
|---|---|---|---|
| Mouse | Tabula Muris Senis droplet atlas, 23 mice | Mouse | Differential expression, signature derivation, composition, noise, and targeted markers |
| Public human cohorts | 11 CELLxGENE cohorts, 621 scored control donors | Donor | Library-normalized signature transport and organ cell-fraction analysis |
| Human lifespan atlas | PsychAD DLPFC, 284 control donors, 27 cell types | Donor | Covariate-adjusted lifespan and adult-stratum sensitivity models |
| Naked mole-rat | GSE152886, five animals | Animal | Descriptive skin comparison and equivalence bounds |
| Methylation | GSE199979 and GSE184211, 1,000 array profiles | Array profile | Published Clock 2 benchmark with per-sample predictive metrics |
| Longevity context | AnAge Build 14, 4,141 species | Species record | Lifespan context, not single-cell evidence for every species |

The primary score is the 139-gene `P1_Clean` signature, with 89 up-regulated and 50 down-regulated mouse genes. The 257-gene `P0_Raw257` score is a sensitivity analysis. Public human cells are library-size normalized before donor-cell-type aggregation. PsychAD donor-cell-type pseudobulk is CPM normalized before log transformation.

One BH family covers 225 public-cohort cell-type tests and 27 primary covariate-adjusted PsychAD cell-type tests. PsychAD adult models are a separate secondary sensitivity family.

## Main results

1. **No primary whole-cortex association was detected.** PsychAD slopes were +0.000186 per year before adjustment (p = 0.3494) and -0.000209 after adjustment (p = 0.4004).
2. **Three primary cell-type tests survived global correction.** These were MSSM L2/3 neurons, MSSM oligodendrocytes, and PsychAD EN_L3_5_IT_3 neurons. The PsychAD neuronal association attenuated in adults.
3. **Glial transport is a candidate, not an established mechanism.** No OPC or microglial result survived the 252-test primary correction. PsychAD adult microglia survived only the secondary adult-only family.
4. **Cell-type fractions change in some organs, not universally.** Eight of 13 adequately covered mouse organs and 3 of 11 human cohort datasets had age-associated fraction changes. Integrated human lung and living-donor kidney results were sensitive to observed-record filtering. Five NMR animals cannot establish stability. See `report/organ-cell-composition.md`.
5. **Cross-species remodeling concordance was not demonstrated.** Three of seven homologous pairs had the same slope sign, binomial p = 1.000, and Spearman rho = 0.000.
6. **Mouse fibroblast Cdkn1a increased with age.** The slope was +0.1080 per month with q = 0.0322. The nominal transcriptional-noise increase depended on the UMI sampling model.
7. **The NMR comparison is inconclusive.** The slope was +0.0058 per normalized lifespan with p = 0.9332 in five animals. The post hoc equivalence margin exceeded the mouse whole-body point estimate.
8. **Clock correlation and prediction are distinct.** Mouse liver predictive R-squared was 0.3169, compared with squared Pearson correlation of 0.851. Human multi-tissue predictive R-squared was 0.6383. Consortium overlap prevents an independent-validation claim.

## Donor overlap and data use

PsychAD contains MSSM and HBCC source-bank donors. Direct identifiers are not shared with the public cohort exports. Age, sex, and bank matching gives possible and unique-candidate counts, not confirmed individual identities. PsychAD is therefore a second analysis of overlapping source populations, not independent replication.

Only aggregate PsychAD regression results are public. Controlled donor-level data and derived donor-level scores remain outside Git. The exact AD Knowledge Portal terms are recorded in `notes/adkp-data-use-verification.md`.

## Reviewer guide

Start with `README.md` for reproduction and `report/manuscript-aging-mechanisms.md` for the current claims. `HANDOFF.md` describes operations. `notes/crucibles/` contains historical audits whose numerical claims may predate the current pipelines. `PLAN.md` and `KICKOFF.md` are historical planning documents, not evidence of preregistration.

## Principal limitations

- Cross-sectional sampling and partial sex-by-age confounding in TMS.
- Public human cohorts have fewer technical covariates than PsychAD.
- Source-population overlap between PsychAD and public MSSM/HBCC cohorts.
- Detection and study-origin effects in integrated lung atlases.
- Different tissues and expression summaries across species, with no universal biological effect scale.
- Five NMR animals, group-mean ages, skin-only sampling, and incomplete ortholog coverage.
- Methylation evaluation data are not demonstrated to be independent of clock development data.
- App-data anchors check consistency with committed results, not the correctness of the scientific design.
