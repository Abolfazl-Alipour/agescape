# M9 Addendum — The Human Cohort Arm: Does the Mouse Aging Program Travel?

Date: 2026-10-21. Pipeline: `scripts/run_m9_pipeline.py` (CELLxGENE census, stable 2025-11-08).
Companion note: `notes/m9-data-scoping.md`. Data: `results/m9-donor-scores.csv`,
`results/m9-cohort-summary.csv`, `results/m9-regressions.csv`, `figures/m9-*.png`,
`results/m9-provenance.log`. Web: graph G15, export `humancohorts.json`.

## Question

Chapter 5 previously reported the human arm as untestable with public data (Tabula
Sapiens skin, 4 donors). M9 asks the question directly: does the 261-gene mouse aging
program (167 up / 94 down; the exact set scored in the app) track donor age in human
single-cell data that carries real per-donor ages?

## Data

Eleven public CELLxGENE cohorts, control donors only (disease == "normal") with numeric
ages from HsapDv `development_stage`: MSSM DLPFC, HBCC DLPFC, RADC DLPFC, SEA-AD DLPFC,
normal adult kidney, living-donor kidney, adult kidney v1.5 (snRNA and scRNA), lung
CellRef 1.0, integrated human lung atlas, CZI immune-aging BAL. 624 control donors
selected; 621 scored after per-cohort filters; 8,280 donor by cell-type units; 98.4%
program-gene coverage in every cohort (188 human orthologs of the 261 mouse genes;
3 genes map to both directions in human space and net out, matching M6).

## Methods

Per cell: program score = mean log1p expression of up genes minus down genes.
Per donor by cell type (>=20 cells): mean score. OLS per cohort by cell type and per
whole cohort: score ~ age + sex. Pre-registered M4 floor for the power verdict
(>=6 donors and donor-age SD >=12 years). No new statistical inventions beyond M2-M6.

## Results (whole-cohort, control donors)

| cohort | n | ages | slope/yr (x10^-3) | 95% CI (x10^-3) | p | verdict |
|---|---|---|---|---|---|---|
| MSSM DLPFC | 213 | 19-89 | -0.193 | [-0.366, -0.021] | 0.028 | POWERED |
| HBCC DLPFC | 134 | 2-85 | +0.168 | [-0.032, +0.369] | 0.100 | POWERED |
| RADC DLPFC | 18 | 64-89 | +0.902 | [-1.222, +3.027] | 0.380 | LOW |
| SEA-AD DLPFC | 21 | 29-89 | +0.148 | [-0.206, +0.501] | 0.392 | POWERED |
| Kidney (4 cohorts) | 65 pooled | 14-76 | all CIs span 0 | | 0.29-0.78 | mixed |
| Lung (3 cohorts) | 170 pooled | 2-81 | all CIs span 0 | | 0.49-0.89 | mixed |

Headline: the mouse program does NOT rise with age in any human cohort. In the largest
cohort (MSSM, n=213) the whole-donor score significantly declines with age
(p = 0.028) — human DLPFC ages away from the mouse lung/kidney old-age direction.

## Results (cell type: the direction flips)

The signal is cell-type-private, in opposite directions by compartment:

- UP with age (meets floor): SEA-AD microglia (+1.16x10^-3/yr, p=0.005), SEA-AD OPC
  (+0.83, p=0.0001), SEA-AD SST interneurons (+0.43, p=0.007), SEA-AD oligodendrocyte
  (+0.62, p=0.013), lung integrated acinar cells (+8.47, p=0.0015), basal and goblet cells.
- DOWN with age (meets floor): MSSM L2/3 IT excitatory neurons (-0.49, p=0.0002),
  MSSM VIP interneurons (-0.33, p=0.007), MSSM GABAergic neurons (-0.27, p=0.011).

Glia and secretory epithelia move toward the mouse old-age program with age; human
excitatory and inhibitory cortical neurons move against it. This agrees directionally
with the DLPFC lifespan atlas (Yang et al., Nature 2026), which reports late-life glial
immune reactivation and cell-type-specific trajectories — and it independently
reproduces that study's central architecture in control donors only, with a mouse
program as the readout.

## Interpretation and limits

1. Cross-tissue transport: the program was derived in mouse lung/kidney. Its failure
   to rise in human brain, kidney, and lung is evidence AGAINST one universal program
   and FOR the project's cell-type-private thesis. It is not evidence that human
   tissues do not age — the per-cell-type table shows they do, in opposite directions.
2. Whole-donor scores weight cell types by abundance (oligodendrocyte-heavy in DLPFC);
   per-cell-type rows are the correct read.
3. Disease cohorts contribute control donors only; composition differs across cohorts
   and could confound slopes. Replication across the four DLPFC cohorts (direction of
   glial rise is consistent in SEA-AD and MSSM OPC/oligo) mitigates but does not
   eliminate this.
4. Exploratory multiplicity: cell-type p-values are not FDR-corrected in this pass;
   the G15 table is a ranking, not a confirmed hit list.
5. The 3-gene human-space direction conflict (CKB, ID2, STMN1) nets out in the score;
   documented in the provenance log.
6. The DLPFC Lifespan Atlas (syn52396927) remains requested under DUA; once granted it
   becomes the replication cohort in a single neurotypical lifespan design.

## Bottom line for the app and report

The human question is no longer "untestable." With public data alone: the mouse aging
program does not travel to humans as a whole-organism signal; where it does appear, its
direction is cell-type-private — up in glia and secretory lung cells, down in cortical
neurons. Aging is not one program, in a second species, with a second method.
