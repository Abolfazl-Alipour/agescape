# Depth versus age: Gao heart and Yazar blood

Date: 2026-11-17. Donor-level analysis using only the committed donor tables. No count matrices were loaded. Read-only on the repo.

## Inputs

- `results/followup-between-subject-gao-donors.csv`: Gao heart, GSE290367, 117 donors, ages 20 to 91 years, 7 studies.
- `results/followup-between-subject-yazar-donors.csv`: Yazar blood, GSE196830, 981 donors, ages 19 to 97 years, 75 pools.

Each table has one row per donor per cell type. This analysis uses only the single pooled row per donor (`All_LV_Nuclei` for Gao, `All_Cells` for Yazar), so the donor is the unit and N equals the donor count. Age, sex, and study or pool are constant within donor in both tables.

Depth means the total number of RNA molecules captured per cell. The tables carry it as mean log UMI (the average of the log molecule count) per donor.

## Method

Ordinary least squares (OLS), fit with statsmodels 0.15.0, classical standard errors. Slopes are per decade of age.

- Primary outcome: `mean_log_umi_ge1500` (mean log UMI among cells with at least 1500 molecules).
- Sensitivity outcome: `mean_log_umi_all` (mean log UMI among all cells).
- Secondary outcome: log of `n_cells_all` (number of captured cells per donor).
- Exploratory outcome: fraction of cells with at least 1500 molecules (`n_cells_ge1500` divided by `n_cells_all`). Added because the primary outcome is computed on a selected subset: cells that fall below 1500 molecules leave the average, which can hide a real decline.
- Model A: outcome on age only.
- Model B: outcome on age plus sex plus batch. Gao uses 6 study indicator columns (reference study ENCODE v4 Snyder). Yazar uses 75 pool indicators (74 dummies plus reference). Pool dummies were feasible: 981 donors, 9 to 14 donors per pool.

## Results

Slope is change per decade. CI is the 95 percent confidence interval.

### Gao heart, N = 117 donors

| Outcome | Model | Slope | CI low | CI high | p | R2 |
|---|---|---:|---:|---:|---:|---:|
| mean_log_umi_ge1500 | A unadjusted | -0.0390 | -0.0790 | 0.0009 | 0.055 | 0.032 |
| mean_log_umi_ge1500 | B adjusted | -0.0021 | -0.0192 | 0.0150 | 0.81 | 0.85 |
| mean_log_umi_all | A unadjusted | -0.0912 | -0.1704 | -0.0120 | 0.024 | 0.043 |
| mean_log_umi_all | B adjusted | -0.0268 | -0.0601 | 0.0065 | 0.11 | 0.85 |
| log_n_cells_all | A unadjusted | 0.0153 | -0.0431 | 0.0736 | 0.61 | 0.002 |
| log_n_cells_all | B adjusted | -0.0073 | -0.0515 | 0.0370 | 0.75 | 0.50 |
| frac_cells_ge1500 (exploratory) | A unadjusted | -0.0357 | -0.0656 | -0.0057 | 0.020 | 0.046 |
| frac_cells_ge1500 (exploratory) | B adjusted | -0.0172 | -0.0336 | -0.0008 | 0.040 | 0.75 |

### Yazar blood, N = 981 donors

| Outcome | Model | Slope | CI low | CI high | p | R2 |
|---|---|---:|---:|---:|---:|---:|
| mean_log_umi_ge1500 | A unadjusted | 0.0054 | 0.0023 | 0.0084 | 0.0006 | 0.012 |
| mean_log_umi_ge1500 | B adjusted | -0.0003 | -0.0024 | 0.0019 | 0.81 | 0.66 |
| mean_log_umi_all | A unadjusted | 0.0345 | 0.0246 | 0.0444 | 1.7e-11 | 0.045 |
| mean_log_umi_all | B adjusted | 0.0030 | -0.0001 | 0.0061 | 0.055 | 0.93 |
| log_n_cells_all | A unadjusted | -0.0159 | -0.0253 | -0.0064 | 0.0010 | 0.011 |
| log_n_cells_all | B adjusted | -0.0272 | -0.0370 | -0.0174 | 6.9e-08 | 0.25 |
| frac_cells_ge1500 (exploratory) | A unadjusted | 0.0348 | 0.0255 | 0.0442 | 5.8e-13 | 0.052 |
| frac_cells_ge1500 (exploratory) | B adjusted | 0.0024 | -0.0003 | 0.0051 | 0.083 | 0.94 |

Full numbers are in `depth-vs-age-results.csv` in this directory.

## Verdicts in plain words

Gao heart: depth does not decline significantly with age once study is accounted for. The unadjusted decline is study confounding: adding study indicators moves R2 from 0.03 to 0.85, meaning depth differs far more between studies than with age. The exploratory share of cells clearing 1500 molecules does decline slightly after adjustment, minus 1.7 percentage points per decade, p 0.040, but that is one exploratory outcome with no correction for multiple testing, and the magnitude is small. It is a flag for a closer look, not a finding. The number of captured cells per donor does not change with age.

Yazar blood: depth does not decline with age. Unadjusted, depth slightly increases with age, and that too is pool confounding: adding pool indicators moves R2 from 0.01 to 0.66 for the primary outcome and from 0.05 to 0.93 for the all-cells outcome, and the slopes go to zero. After adjustment there is no depth change in either direction. The number of captured cells per donor does decline with age, about 2.7 percent fewer cells per decade, and that survives pool adjustment (p 6.9e-08). That is a statement about how many cells were captured per donor, not about depth per cell, and cell counts also reflect sampling and loading.

Overall: in these two large cohorts, per-cell depth shows no robust decline with age after batch adjustment. The earlier MSSM observation that depth falls with age does not generalize to Gao heart or Yazar blood at the donor level. This does not prove depth never changes with age anywhere; it says these two cohorts show no adjusted decline.

## Limits, stated plainly

1. Depth is confounded with technical factors: tissue dissociation, capture efficiency, sequencing run, and postmortem interval. A depth change would not by itself prove a biological mechanism, and no mechanism is claimed here.
2. The primary outcome conditions on cells with at least 1500 molecules, which selects against low-depth cells. The all-cells outcome and the exploratory fraction outcome were included to check that selection. The fraction result in Gao is the only adjusted signal and it is exploratory.
3. Age enters as a straight-line slope. A nonlinear age pattern would not be captured.
4. Gao has only 117 donors across 7 studies, and study is strongly tangled with depth. The adjusted Gao estimates are precise about the study-adjusted slope but rest on few donors in some studies (3 to 54 donors per study).
5. Yazar pool indicators absorb most depth variation (R2 0.93). If pool assignment correlates with a real biological depth change, adjustment could remove real signal along with batch signal. Pool explains part of age variation too, which is why the unadjusted and adjusted Yazar slopes differ in sign.
6. These are donor-level means. Nothing here describes how depth varies across cells within a donor, or within a cell type.
7. No correction for multiple testing was applied across the outcomes. The Gao exploratory p value of 0.040 would not survive one.

## Data problems

None blocking. Notes: the donor tables are donor-by-cell-type rows, not one row per donor; only the pooled row per donor was used. The ge1500 outcome has missing values in some small cell-type rows (134 in Gao, 1083 in Yazar) but zero missing values in the pooled rows used here. One Gao donor has a fraction of cells above 1500 molecules of 0.001, far below the median of 0.62; it was kept in all models.
