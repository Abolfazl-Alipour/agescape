# Dataset 07: Human Retina (Yan et al., 2020, GSE148077)

## Verdict
**No global gene-level aging signature (`0 / 120604` tests), cross-cell-type directional concordance (`r = +0.1037, p = 0.9321`), or joint analytic + empirical `P1_Clean` transport signal was detected across the unbiased foveal retina cohort (`n = 6` donors, ages 53 to 78 yr) or the 7-donor sensitivity (`53 to 78 yr`).** One major-class fraction shifted nominally and after Family A BH correction (`RGC` CLR slope `+0.9635` per decade, `q_BH = 0.0052`), reflecting wide inter-donor foveal layer/eccentricity variation (`3.1% to 18.2%` RGC, `0.5% to 38.3%` Rod) across `n = 6` postmortem terminal-disease donors (`PMI 3 to 14 h`) lacking anyone under age `53 yr`.

## Primary Source Verification
- **Publication**: Yan et al. (2020) *Scientific Reports* 10:9802 (PMID 32555229, PMCID PMC7299956).
- **Exact donor reconciliation**: GEO `GSE148077` GSM records (`GSM4453751` to `GSM4453769`) omit donor ages and sexes. Primary **Table S1** in the published Springer supplementary PDF (`41598_2020_66092_MOESM1_ESM.pdf`) establishes the exact 7-donor table (`H9 53F`, `H3 60M`, `H4 64M`, `H11 65M`, `H5 69M`, `H1 74M`, `H2 78M`), correcting secondary notes (`H2 53M / H3 78F / H4 78M` was scrambled; only `53 to 78` matched the primary span).
- **Raw count verification**: Extracted from `Sanes_GSE148077` inside the primary HRCA scRNA-seq integration (`9847a0f8-18a0-411c-bff9-827a607e7fcf.h5ad`), yielding **93,842 cells across all 19 libraries** mapped 1-to-1 to Table S3 (`H11` coded `Hu11` with `H11FoveaS1_*` barcodes; `H2` coded `Pt2` with `H2Fovea1S1_*` and `H2Fovea2S1_*` barcodes). `raw/X` contains 100% finite integer UMI counts (`max = 2,399`), collapsed across 35,475 Ensembl IDs to 35,466 unique HGNC symbols.
- **Unbiased vs. sorted tissue**: Only the `6` foveal donors (`H9, H3, H4, H11, H5, H2`, `63,438` cells) are unbiased whole-fovea dissociations (`suspension_enrichment_factors = na`, `suspension_depletion_factors = na`). Peripheral libraries (`H1` peripheral-only `74M` and peripheral `H3 60M`) used `CD73` depletion (removing rods) or `CD90` enrichment (isolating RGCs/amacrine/astrocytes), so `All7_MajorClass_SortedIncluded` excludes Family A composition.

## Cohort Summary (`n = 7` total, `n = 6` unbiased fovea)
| Donor (HRCA / Paper) | Age (yr) | Sex | 10x Kit | Sampling Protocol | Total Cells | Unbiased Fovea Cells | Cause of Death (Table S1) |
|---|---:|:---:|:---:|---|---:|---:|---|
| H9 (H9) | 53 | F | V2 | fovea_only | 9187 | 9187 | Interstitial Lung Disease |
| H3 (H3) | 60 | M | V2 | fovea_and_sorted_periphery | 27514 | 8831 | Metastatic tonsillar SCC |
| H4 (H4) | 64 | M | V2 | fovea_only | 3807 | 3807 | Diffuse B cell lymphoma |
| Hu11 (H11) | 65 | M | V3 | fovea_only | 1453 | 1453 | Metastatic melanoma |
| H5 (H5) | 69 | M | V2 | fovea_only | 32059 | 32059 | Metastatic melanoma |
| H1 (H1) | 74 | M | V2 | sorted_periphery_only | 11721 | 0 | Lung cancer |
| Pt2 (H2) | 78 | M | V2 | fovea_only | 8101 | 8101 | Metastatic melanoma |

## Family A: Compositional Shifts (CLR OLS, Unbiased Fovea `n = 6`)
- Primary `Fovea_MajorClass_Primary` (`7` classes `>=5` cells in `>=4` donors; `~ age_per_decade`): **`1` / 7 classes at `q_BH < 0.05`** (`RGC`: `effect_clr = +0.9635`, `p = 0.000736`, `q_BH = 0.0052`, driven by lower RGC % in `H9` `3.29%` and `H3` `3.11%` vs `H5` `18.22%` and `Pt2/H2` `14.01%`).
- `Fovea_MajorClass_SexAdj` (`~ age_per_decade + sex_male`): **`1` / 7 class at `q_BH < 0.05`** (`RGC`: `effect_clr = +0.7904`, `q_BH = 0.0138`).
- Finer `Fovea_CellType` (`15` cell ontology classes): **`0` classes at `q_BH < 0.05`** (`midget ganglion cell of retina` reaches `p = 0.0080, q_BH = 0.1194`).

| Major Class | Donors Present | Fovea Cells | Median Fraction | Effect (CLR/dec) | p-value | q_BH (Family A) |
|---|---:|---:|---:|---:|---:|---:|
| RGC | 6 | 8148 | 0.1125 | +0.9635 | 0.0007 | 0.0052 |
| Cone | 6 | 2142 | 0.0253 | +0.3067 | 0.0367 | 0.1285 |
| MG | 6 | 12382 | 0.1920 | +0.5564 | 0.1656 | 0.3864 |
| Rod | 6 | 7435 | 0.0558 | -0.7906 | 0.4188 | 0.6204 |
| AC | 6 | 5114 | 0.0919 | -0.1748 | 0.4431 | 0.6204 |
| BC | 6 | 24717 | 0.3704 | +0.3241 | 0.5873 | 0.6852 |
| HC | 6 | 3295 | 0.0405 | +0.2090 | 0.7454 | 0.7454 |

## Family B: Donor-Level Gene Age Tests and Directional Concordance
- Primary global Family B (`Fovea_MajorClass_Primary`, `M_B = 120604` gene x major-class tests across `7` classes, model `~ age_per_decade + log_lib`): **`0` hits at `q_BH < 0.05`** (minimum global `q_BH = 1.0000`).
- Zero Family B global hits also occurred in `Fovea_CellType` (`M_B = 224,091`), `Fovea_Pooled` (`M_B = 18,327`), and `All7_MajorClass_SortedIncluded` (`M_B = 130,244`). In the underpowered `Fovea_MajorClass_SexAdj` sensitivity (`n = 6`, `residual_df = 2` with 4 parameters), `1` singleton (`AC: FANCF`, `r2 = 1.0`, `t = -1975.4`, `q_BH = 0.0309`) is an exact overfitting artifact of boundary female `H9`.
- Cross-cell-type directional concordance (1,000 permutations): primary `r = +0.1037` (`p_one_sided = 0.9321`).

| Major Class | n_donors | Genes Tested | Within q<0.05 | Global q<0.05 | Min p-value | Min q_BH Global |
|---|---:|---:|---:|---:|---:|---:|
| AC | 6 | 18542 | 0 | 0 | 1.88e-04 | 1.0000 |
| BC | 6 | 18027 | 0 | 0 | 1.42e-04 | 1.0000 |
| Cone | 6 | 16387 | 0 | 0 | 1.74e-04 | 1.0000 |
| HC | 6 | 17563 | 0 | 0 | 7.35e-05 | 1.0000 |
| MG | 6 | 17797 | 0 | 0 | 4.70e-05 | 1.0000 |
| RGC | 6 | 16566 | 0 | 0 | 7.45e-04 | 1.0000 |
| Rod | 6 | 15722 | 0 | 0 | 1.92e-03 | 1.0000 |

### Cross-Cell-Type Concordance Summary
| Suite | n_cell_types | Common Genes | Mean Pairwise Pearson r | Perm p (1-sided) | Perm p (2-sided) |
|---|---:|---:|---:|---:|---:|
| Fovea_MajorClass_Primary | 7 | 13425 | +0.1037 | 0.9321 | 0.9321 |
| Fovea_MajorClass_SexAdj | 7 | 13425 | +0.0547 | 0.8751 | 0.8751 |
| Fovea_CellType | 15 | 6123 | +0.0675 | 0.9760 | 0.9760 |
| All7_MajorClass_SortedIncluded | 8 | 11448 | +0.1188 | 0.9221 | 0.9221 |

## Family C: P1_Clean Mouse Signature Transport
- **Zero cell types or suites reached `q_BH < 0.05` on BOTH analytic OLS and the 200-draw empirical null** (minimum analytic `q_BH_comp = 0.4245`; `min q_BH_emp = 0.0348`).
- Specifically, Muller glia (`MG`, `Fovea_MajorClass_Primary`, `n = 6`) showed a small positive point estimate (`composite_effect = +0.2053`) that ranked above most 200 random gene sets (`empirical_null_pvalue = 0.0050`, `q_BH_emp = 0.0348`) solely under unadjusted variance across `n = 6` donors, while its actual donor-level OLS was completely non-significant (`composite_pvalue = 0.1376`, `q_BH_comp = 0.6750`) and dropped above `0.05` in `Fovea_MajorClass_SexAdj` (`q_BH_emp = 0.0697, q_BH_comp = 0.8407`), `Fovea_CellType` (`Mueller cell`, `q_BH_emp = 0.0746`), and `All7_MajorClass_SortedIncluded` (`n = 7`, `q_BH_emp = 0.1194, q_BH_comp = 0.9217`).

| Suite | Group | n_donors | P1 Up/Down | Composite Effect | Analytic p | Empirical p (200 draws) | q_BH (Comp) | q_BH (Emp) |
|---|---|---:|:---:|---:|---:|---:|---:|---:|
| Fovea_MajorClass_Primary | MG | 6 | 66/32 | +0.2053 | 0.1376 | 0.0050 | 0.6750 | 0.0348 |
| Fovea_MajorClass_Primary | AC | 6 | 62/33 | +0.1784 | 0.2893 | 0.0498 | 0.6750 | 0.1741 |
| Fovea_MajorClass_Primary | HC | 6 | 61/30 | +0.1517 | 0.2186 | 0.1741 | 0.6750 | 0.4063 |
| Fovea_MajorClass_Primary | BC | 6 | 59/31 | -0.0484 | 0.6642 | 0.2438 | 0.7749 | 0.4249 |
| Fovea_MajorClass_Primary | RGC | 6 | 59/27 | -0.0491 | 0.6636 | 0.3035 | 0.7749 | 0.4249 |
| Fovea_MajorClass_Primary | Rod | 6 | 53/24 | +0.1280 | 0.4279 | 0.4577 | 0.7488 | 0.5340 |
| Fovea_MajorClass_Primary | Cone | 6 | 57/24 | -0.0012 | 0.9974 | 0.9751 | 0.9974 | 0.9751 |
| Fovea_Pooled | All_Retina_Pooled | 6 | 64/33 | +0.1223 | 0.4245 | 0.0995 | 0.4245 | 0.0995 |
| All7_MajorClass_SortedIncluded | MG | 7 | 65/32 | +0.1316 | 0.5192 | 0.0149 | 0.9217 | 0.1194 |
| All7_MajorClass_SortedIncluded | HC | 7 | 61/30 | +0.2049 | 0.1170 | 0.0348 | 0.9217 | 0.1393 |
| All7_MajorClass_SortedIncluded | AC | 7 | 62/32 | +0.1612 | 0.2395 | 0.0647 | 0.9217 | 0.1493 |
| All7_MajorClass_SortedIncluded | Astrocyte | 4 | 58/24 | +0.2267 | 0.6709 | 0.0746 | 0.9217 | 0.1493 |
| All7_MajorClass_SortedIncluded | BC | 7 | 60/31 | -0.0382 | 0.7755 | 0.5423 | 0.9217 | 0.7629 |
| All7_MajorClass_SortedIncluded | Rod | 7 | 51/21 | -0.0797 | 0.5991 | 0.5721 | 0.9217 | 0.7629 |
| All7_MajorClass_SortedIncluded | Cone | 7 | 56/24 | +0.0297 | 0.9217 | 0.8507 | 0.9217 | 0.9502 |
| All7_MajorClass_SortedIncluded | RGC | 7 | 58/26 | +0.0116 | 0.9217 | 0.9502 | 0.9217 | 0.9502 |

## Key Limitations
1. **No young donors**: all 7 donors are ages `53 to 78 yr` (`25 yr` range), testing only late-midlife-to-old differences rather than young-vs-aged aging.
2. **Small `n` and sex collinearity**: unbiased fovea has `n = 6` (`df = 3` unadjusted for sex; single female `H9` is also the youngest donor at `53 yr`).
3. **Postmortem terminal illness & 10x kit heterogeneity**: 6/7 donors died of advanced malignancies (`H9` had ILD), PMI spans `3 to 14 h` (`Pt2/H2 78M` longest at `14 h`), foveal dissection eccentricity varies widely across donors (`Rod 0.5% to 38.3%`), and `H11` used 10x 3' v3 while the other 6 used 10x 3' v2.
4. **Source provenance**: counts were analyzed from the HRCA CellxGene v1-curated UMI extraction of `GSE148077` rather than re-processing FASTQs.

## Output Files
- `results/followup-organ-07-eye-yan-gse148077-donors.csv`
- `results/followup-organ-07-eye-yan-gse148077-fractions.csv`
- `results/followup-organ-07-eye-yan-gse148077-gene-summary.csv`
- `results/followup-organ-07-eye-yan-gse148077-top-genes.csv`
- `results/followup-organ-07-eye-yan-gse148077-concordance.csv`
- `results/followup-organ-07-eye-yan-gse148077-p1-signature.csv`
