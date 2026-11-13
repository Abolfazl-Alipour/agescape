# Follow-up Step 1: De Man lung aging

## Status

IPF Cell Atlas control analysis is complete. The Baylor-only analysis is not complete because the GEO raw matrices do not include cell annotations or the integrated object. No labels were invented, and no full 60-donor reproduction is claimed.

## Verified inputs

- Supplementary Data 1 contains 60 subjects with exact ages and smoking status: 32 Baylor subjects and 28 IPF Cell Atlas controls.
- GSE281219 contains 32 Baylor raw matrices. GSE136831 contains 96,303 control cells from 28 subjects and 38 libraries. All 28 controls map to Supplementary Data 1.
- IPF libraries are pooled by subject for inference. Gene filtering passes 92,049 cells; regression eligibility additionally requires QC and the donor/model floors.
- Target labels are ATI, ATII, aerocyte (`VE_Capillary_A`), and gCap (`VE_Capillary_B`).
- Source provenance and hashes are in `data/manifest/de_man_sources.json`; donor and cell-count tables are in `results/followup-de-man-donor-crosswalk.csv` and `results/followup-de-man-celltype-counts.csv`.

## Methods

The published-KL arm follows the released entropy code: zero-gene removal, independent binomial thinning toward the minimum library, square-root transformation, and KL divergence from a uniform gene distribution. The analysis uses donor means and continuous-age models with sex and mean log actual UMI as covariates. Cell-level Wilcoxon tests were not used for inference.

The Shannon arm computes direct within-cell Shannon entropy after mitochondrial, ribosomal, and sex-linked exclusions. Fixed-depth analysis uses exact multinomial sampling to 1,500 UMIs for cells whose post-exclusion library is at least 1,500 UMIs. Full-depth analysis uses no downsampling. The full-depth result is identical across random seeds.

The primary seed is 20261113; 20 consecutive seeds are reported for sampling robustness. BH correction is applied across the four target cell types separately for each method, using the primary seed only. Seed replicates are not treated as independent significance tests.

## Primary results

No result had q below 0.05 in the primary donor-level models.

| Method | Cell type | Donors | Slope per year | 95% CI | p | q |
|---|---|---:|---:|---:|---:|---:|
| Published KL | AT1 | 19 | -0.000382 | -0.001899 to +0.001134 | 0.599 | 0.637 |
| Published KL | AT2 | 18 | -0.000559 | -0.002377 to +0.001259 | 0.521 | 0.637 |
| Published KL | Aerocyte | 17 | -0.000672 | -0.003675 to +0.002331 | 0.637 | 0.637 |
| Published KL | gCap | 22 | -0.002005 | -0.004181 to +0.000170 | 0.0687 | 0.275 |
| Shannon, 1,500 UMIs | AT1 | 24 | +0.000629 | -0.003589 to +0.004846 | 0.759 | 0.759 |
| Shannon, 1,500 UMIs | AT2 | 26 | -0.000869 | -0.004241 to +0.002504 | 0.599 | 0.759 |
| Shannon, 1,500 UMIs | Aerocyte | 17 | +0.001752 | -0.006173 to +0.009678 | 0.641 | 0.759 |
| Shannon, 1,500 UMIs | gCap | 22 | +0.004459 | -0.002688 to +0.011605 | 0.206 | 0.759 |
| Shannon, full depth | AT1 | 25 | +0.004575 | -0.002869 to +0.012018 | 0.215 | 0.308 |
| Shannon, full depth | AT2 | 26 | -0.001983 | -0.006642 to +0.002676 | 0.387 | 0.387 |
| Shannon, full depth | Aerocyte | 17 | +0.005159 | -0.003716 to +0.014034 | 0.231 | 0.308 |
| Shannon, full depth | gCap | 22 | +0.005619 | -0.002300 to +0.013538 | 0.153 | 0.308 |

Negative published-KL slopes mean increasing entropy with age, because the published score is divergence from a uniform distribution. Positive Shannon slopes mean increasing entropy with age. These are not interchangeable with off-program expression or Hellinger displacement.

## Robustness across 20 seeds

- Published-KL gCap was negative in every seed and significant at p < 0.05 in 6 of 20 seeds. This is sampling-sensitive evidence, not a confirmed entropy increase.
- Published-KL AT2 was negative in 19 of 20 seeds and significant in 2. AT1 and aerocyte were negative in all 20 seeds and significant in none.
- Full-depth gCap, aerocyte, and AT1 had positive Shannon slopes; AT2 was negative. None was significant.
- At 1,500 UMIs, gCap, aerocyte, and AT1 were positive in every seed; AT2 was negative in every seed. None was significant.
- Per-seed effects, intervals, and p-values are in `results/followup-de-man-entropy-models.csv`; the compact summary is `results/followup-de-man-entropy-robustness.csv`.

## Interpretation and limits

The primary IPF-control analysis does not reproduce an FDR-significant age-associated entropy increase in the four published target cell types. gCap shows a directionally consistent published-KL increase, but its significance depends on cell sampling. This is not evidence that entropy is stable: confidence intervals and cohort composition limit what can be concluded.

The published paper used cell-level comparisons, which give cells disproportionate weight relative to donors. This analysis instead makes donors the inferential unit. Different estimands, a smaller IPF-only cohort, QC, and depth handling can therefore produce different results.

## Not done

- Baylor cell annotation and Baylor entropy models were not run because no annotation or integrated object was released in the raw GEO package. A new annotation was not substituted for the published labels.
- No full 60-donor published-method reproduction is claimed.
- P1 transport, Hellinger displacement, composition, and off-program expression were not run because the Baylor annotation is missing for the new dataset, and existing integrated-lung analyses already cover the overlapping HLCA resource.
- No somatic-mutation analysis was attempted.
- No journal submission or external contact was made.

Peak memory for the completed entropy run was 0.85 GB. Tests for entropy functions and downsampling passed.
