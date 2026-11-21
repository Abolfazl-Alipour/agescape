# Item 7 mature interneuron replication and mature-derived SEA-AD test

## Verdict

PROVES. Significant positive genome-wide replication and both-FDR overlap above hypergeometric expectation survives in VIP, pvalb, sst. Human mature means donors age 30 years and older in both the discovery fits and the tested MSSM-to-HBCC contrast. The verdict uses the unmoderated WLS replication files. The mature moderated intersection is UNVERIFIED because the frozen recovery gate failed.

## Models and fixed families

MSSM and HBCC use donor-summed raw counts, log1p counts per 10,000, one-iteration voom WLS with age plus sex plus donor mean log UMI, and donor-unit inference. QC is primary normal data, at least 1500 UMIs per cell, at least 20 cells per donor and subtype, and at least 10 donors per unit. Depth-adjusted OLS uses the same design. The depth-unadjusted sensitivity uses age plus sex. BH correction is within cohort, subtype, and model over all tested genes. The seed is 20261108 and every block permutation uses 1000 iterations.

The SEA-AD fit is unchanged: 15 donors age 72 to 102, raw UMIs, released QC, Vip/Pvalb/Sst mappings, and WLS age plus sex plus mean log UMI plus multiome fraction plus Low ADNC indicator. Discovery is `mature MSSM age_plus_sex_depth OLS (30-plus)`. SEA-AD OLS is primary for genome-wide correlation and SEA-AD WLS is secondary. Targeted same-direction tests use WLS. Permutations and block bootstraps each use 1000 iterations with seeds 20261108 and 20261208.

## Donor retention

| Subtype | MSSM full N | MSSM mature N | MSSM mature ages | MSSM median | HBCC full N | HBCC mature N | HBCC mature ages | HBCC median |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| VIP | 197 | 193 | 31 to 89 | 67.0 | 80 | 59 | 30 to 85 | 48.0 |
| pvalb | 189 | 185 | 31 to 89 | 67.0 | 80 | 59 | 30 to 85 | 48.0 |
| sst | 177 | 174 | 31 to 89 | 66.0 | 82 | 61 | 30 to 85 | 48.0 |
| GABA | 155 | 153 | 31 to 89 | 66.0 | 71 | 51 | 30 to 85 | 48.0 |

## FDR counts

| Subtype | Cohort | Full WLS | Full depth-adjusted OLS | Mature WLS | Mature depth-adjusted OLS | Mature depth-unadjusted OLS |
|---|---|---:|---:|---:|---:|---:|
| VIP | MSSM_DLPFC | 218 | 235 | 148 | 145 | 385 |
| VIP | HBCC_DLPFC | 224 | 236 | 8 | 0 | 0 |
| pvalb | MSSM_DLPFC | 152 | 237 | 133 | 147 | 353 |
| pvalb | HBCC_DLPFC | 232 | 219 | 14 | 4 | 7 |
| sst | MSSM_DLPFC | 400 | 230 | 366 | 204 | 359 |
| sst | HBCC_DLPFC | 604 | 540 | 61 | 27 | 27 |
| GABA | MSSM_DLPFC | 75 | 53 | 39 | 32 | 153 |
| GABA | HBCC_DLPFC | 45 | 14 | 0 | 0 | 0 |

## MSSM-to-HBCC replication

Selected same-direction counts are among MSSM FDR genes. Both-FDR counts require the same sign in both cohorts. Hypergeometric and decile-matched values are expected counts. Block p is the coexpression-block permutation p for the both-FDR overlap.

| Subtype | Range | Shared genes | Genome-wide r | Selected r | Selected same/tested | Both FDR | Hypergeometric expected | Matched expected | Block p |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| VIP | Full range | 1853 | 0.7217 | 0.8551 | 190/202 | 76 | 23.547 | 24.180 | 0.0010 |
| VIP | Mature 30-plus | 1862 | 0.7523 | 0.8998 | 136/139 | 4 | 0.597 | 0.574 | 0.0010 |
| pvalb | Full range | 1790 | 0.6839 | 0.8507 | 141/143 | 66 | 18.215 | 19.587 | 0.0010 |
| pvalb | Mature 30-plus | 1797 | 0.6299 | 0.7738 | 120/128 | 8 | 0.997 | 1.146 | 0.0010 |
| sst | Full range | 1808 | 0.7614 | 0.8894 | 363/375 | 249 | 121.543 | 125.962 | 0.0010 |
| sst | Mature 30-plus | 1826 | 0.6996 | 0.8634 | 335/355 | 41 | 11.665 | 12.365 | 0.0010 |
| GABA | Full range | 1781 | 0.3315 | 0.4927 | 57/69 | 0 | 1.743 | 1.886 | 1.0000 |
| GABA | Mature 30-plus | 1792 | 0.4085 | 0.5770 | 31/36 | 0 | 0.000 | 0.000 | 1.0000 |

The original full-range moderated three-subtype intersection is exactly FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16. DYNC1I1 occurs only in the unmoderated full-range intersection. The mature unmoderated three-subtype intersection is: TRMT11. The mature moderated three-subtype intersection is: TRMT11. The latter is UNVERIFIED because the frozen recovery gate failed.

## Six-gene mature status

PASS means q < 0.05, NOT_FDR means tested with q >= 0.05, and NOT_TESTED means absent after the unchanged expression filter. Unmoderated statuses are the confirmed mature result. Moderated statuses are marked UNVERIFIED if the recovery gate failed and are not presented as confirmed replication.

| Gene | Cohort | Subtype | Unmoderated WLS | Moderated WLS |
|---|---|---|---|---|
| FBN1 | MSSM_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| FBN1 | MSSM_DLPFC | pvalb | PASS | PASS (UNVERIFIED) |
| FBN1 | MSSM_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| FBN1 | HBCC_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| FBN1 | HBCC_DLPFC | pvalb | NOT_FDR | NOT_FDR (UNVERIFIED) |
| FBN1 | HBCC_DLPFC | sst | NOT_FDR | NOT_FDR (UNVERIFIED) |
| GABRB1 | MSSM_DLPFC | VIP | PASS | PASS (UNVERIFIED) |
| GABRB1 | MSSM_DLPFC | pvalb | PASS | PASS (UNVERIFIED) |
| GABRB1 | MSSM_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| GABRB1 | HBCC_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| GABRB1 | HBCC_DLPFC | pvalb | NOT_FDR | NOT_FDR (UNVERIFIED) |
| GABRB1 | HBCC_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| IGF1R | MSSM_DLPFC | VIP | PASS | PASS (UNVERIFIED) |
| IGF1R | MSSM_DLPFC | pvalb | PASS | PASS (UNVERIFIED) |
| IGF1R | MSSM_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| IGF1R | HBCC_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| IGF1R | HBCC_DLPFC | pvalb | NOT_FDR | NOT_FDR (UNVERIFIED) |
| IGF1R | HBCC_DLPFC | sst | NOT_FDR | NOT_FDR (UNVERIFIED) |
| PHF3 | MSSM_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| PHF3 | MSSM_DLPFC | pvalb | PASS | PASS (UNVERIFIED) |
| PHF3 | MSSM_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| PHF3 | HBCC_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| PHF3 | HBCC_DLPFC | pvalb | NOT_FDR | NOT_FDR (UNVERIFIED) |
| PHF3 | HBCC_DLPFC | sst | NOT_FDR | NOT_FDR (UNVERIFIED) |
| TRMT11 | MSSM_DLPFC | VIP | PASS | PASS (UNVERIFIED) |
| TRMT11 | MSSM_DLPFC | pvalb | PASS | PASS (UNVERIFIED) |
| TRMT11 | MSSM_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| TRMT11 | HBCC_DLPFC | VIP | PASS | PASS (UNVERIFIED) |
| TRMT11 | HBCC_DLPFC | pvalb | PASS | PASS (UNVERIFIED) |
| TRMT11 | HBCC_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| ZBTB16 | MSSM_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| ZBTB16 | MSSM_DLPFC | pvalb | NOT_FDR | NOT_FDR (UNVERIFIED) |
| ZBTB16 | MSSM_DLPFC | sst | PASS | PASS (UNVERIFIED) |
| ZBTB16 | HBCC_DLPFC | VIP | NOT_FDR | NOT_FDR (UNVERIFIED) |
| ZBTB16 | HBCC_DLPFC | pvalb | NOT_FDR | NOT_FDR (UNVERIFIED) |
| ZBTB16 | HBCC_DLPFC | sst | NOT_FDR | NOT_FDR (UNVERIFIED) |

## Power and interpretation

Slope MDE80 is `(t_0.975 + t_0.80) * median WLS SE` on the log1p CP10K scale. Correlation MDE80 uses the Fisher z approximation with two-sided alpha 0.05 and the shared-gene count. Full-range median SE values come from the frozen all-cell-type audit, not a refit.

| Subtype | Cohort | Range | N | MDE80 per year | MDE80 per 10 years | Shared genes | Correlation MDE80 |
|---|---|---|---:|---:|---:|---:|---:|
| VIP | MSSM_DLPFC | full_range | 197 | 0.001280 | 0.012799 | 1853 | 0.0650 |
| VIP | MSSM_DLPFC | mature30 | 193 | 0.001364 | 0.013635 | 1862 | 0.0649 |
| VIP | HBCC_DLPFC | full_range | 80 | 0.001725 | 0.017248 | 1853 | 0.0650 |
| VIP | HBCC_DLPFC | mature30 | 59 | 0.002380 | 0.023799 | 1862 | 0.0649 |
| pvalb | MSSM_DLPFC | full_range | 189 | 0.001410 | 0.014097 | 1790 | 0.0662 |
| pvalb | MSSM_DLPFC | mature30 | 185 | 0.001487 | 0.014872 | 1797 | 0.0660 |
| pvalb | HBCC_DLPFC | full_range | 80 | 0.001777 | 0.017767 | 1790 | 0.0662 |
| pvalb | HBCC_DLPFC | mature30 | 59 | 0.002404 | 0.024040 | 1797 | 0.0660 |
| sst | MSSM_DLPFC | full_range | 177 | 0.001417 | 0.014168 | 1808 | 0.0658 |
| sst | MSSM_DLPFC | mature30 | 174 | 0.001478 | 0.014779 | 1826 | 0.0655 |
| sst | HBCC_DLPFC | full_range | 82 | 0.001938 | 0.019380 | 1808 | 0.0658 |
| sst | HBCC_DLPFC | mature30 | 61 | 0.002616 | 0.026160 | 1826 | 0.0655 |
| GABA | MSSM_DLPFC | full_range | 155 | 0.001634 | 0.016342 | 1781 | 0.0663 |
| GABA | MSSM_DLPFC | mature30 | 153 | 0.001714 | 0.017139 | 1792 | 0.0661 |
| GABA | HBCC_DLPFC | full_range | 71 | 0.001994 | 0.019936 | 1781 | 0.0663 |
| GABA | HBCC_DLPFC | mature30 | 51 | 0.002835 | 0.028355 | 1792 | 0.0661 |

Versus the original full-range analysis, VIP: genome-wide r 0.7217 to 0.7523; both-FDR 76 to 4; mature/full-range slope MDE80 ratio 1.065 in MSSM and 1.380 in HBCC. pvalb: genome-wide r 0.6839 to 0.6299; both-FDR 66 to 8; mature/full-range slope MDE80 ratio 1.055 in MSSM and 1.353 in HBCC. sst: genome-wide r 0.7614 to 0.6996; both-FDR 249 to 41; mature/full-range slope MDE80 ratio 1.043 in MSSM and 1.350 in HBCC. GABA: genome-wide r 0.3315 to 0.4085; both-FDR 0 to 0; mature/full-range slope MDE80 ratio 1.049 in MSSM and 1.422 in HBCC. HBCC donor loss therefore comes with a real slope-precision cost, while genome-wide correlation MDE80 remains near 0.065 to 0.066 because shared-gene counts stay near 1800. VIP and sst retain large both-FDR excesses despite fewer HBCC hits. A non-significant correlation is not evidence of no aging, and a small MDE80 is not an equivalence test.

## Mature-derived SEA-AD outcome

### Genome-wide effects

Intervals are 95 percent block-bootstrap percentile intervals. Block p is one-sided for positive OLS replication and BH q is across the three subtypes.

| Subtype | Shared genes | OLS Pearson r (95% CI) | WLS Pearson r (95% CI) | Block p | BH q |
|---|---:|---:|---:|---:|---:|
| VIP | 1634 | 0.0308 (-0.0751 to 0.1261) | 0.0028 (-0.1010 to 0.1011) | 0.0749 | 0.1124 |
| pvalb | 1566 | -0.1686 (-0.2689 to -0.0609) | -0.1567 (-0.2540 to -0.0476) | 1.0000 | 1.0000 |
| sst | 1630 | 0.0478 (-0.0326 to 0.1259) | -0.1692 (-0.2365 to -0.1107) | 0.0749 | 0.1124 |

Genome-wide SEA-AD FDR gene counts are 0 VIP, 0 pvalb, and 0 sst. These counts cover all filtered genes within each SEA-AD subtype.

### Target families

The three mature human families use mature unmoderated both-FDR and mature MSSM discovery sets. The mouse P1 families are young-derived from mouse 3 versus 18 months, not mature-derived targets. Same-direction SEA-AD FDR is the subset with SEA-AD q < 0.05. Only genome-wide correlation and mature both-FDR same-direction enrichment have primary three-subtype BH correction; other set-level p values are secondary.

| Subtype | Family | Tested/source | Same direction | Nominal | SEA-AD FDR | Block p | BH q |
|---|---|---:|---:|---:|---:|---:|---:|
| VIP | mature_both_fdr | 4/4 | 1 (25.0%) | 1 | 0 | 0.8691 | 0.9760 |
| VIP | mature_discovery_all | 124/139 | 64 (51.6%) | 13 | 0 | 0.2188 | secondary |
| VIP | mature_discovery_without_both_fdr | 120/135 | 63 (52.5%) | 12 | 0 | 0.1908 | secondary |
| VIP | mouse_p1_all_young_derived | 4/118 | 2 (50.0%) | 1 | 0 | 0.7662 | secondary |
| VIP | mouse_p1_mature_discovery_overlap_young_derived | 1/1 | 0 (0.0%) | 0 | 0 | 1.0000 | secondary |
| VIP | mouse_p1_outside_mature_discovery_young_derived | 3/117 | 2 (66.7%) | 1 | 0 | 0.5355 | secondary |
| pvalb | mature_both_fdr | 7/8 | 2 (28.6%) | 1 | 0 | 0.8482 | 0.9760 |
| pvalb | mature_discovery_all | 121/128 | 46 (38.0%) | 3 | 0 | 0.8881 | secondary |
| pvalb | mature_discovery_without_both_fdr | 114/120 | 44 (38.6%) | 2 | 0 | 0.8611 | secondary |
| pvalb | mouse_p1_all_young_derived | 3/118 | 3 (100.0%) | 0 | 0 | 0.1528 | secondary |
| pvalb | mouse_p1_mature_discovery_overlap_young_derived | 1/1 | 1 (100.0%) | 0 | 0 | 0.6903 | secondary |
| pvalb | mouse_p1_outside_mature_discovery_young_derived | 2/117 | 2 (100.0%) | 0 | 0 | 0.2108 | secondary |
| sst | mature_both_fdr | 38/41 | 16 (42.1%) | 4 | 0 | 0.9760 | 0.9760 |
| sst | mature_discovery_all | 333/355 | 142 (42.6%) | 24 | 0 | 1.0000 | secondary |
| sst | mature_discovery_without_both_fdr | 295/314 | 126 (42.7%) | 20 | 0 | 0.9990 | secondary |
| sst | mouse_p1_all_young_derived | 5/118 | 3 (60.0%) | 0 | 0 | 0.6364 | secondary |
| sst | mouse_p1_mature_discovery_overlap_young_derived | 1/1 | 1 (100.0%) | 0 | 0 | 0.4426 | secondary |
| sst | mouse_p1_outside_mature_discovery_young_derived | 4/117 | 2 (50.0%) | 0 | 0 | 0.8312 | secondary |

### Sensitivities and donor influence

| Subtype | Arm | Status | N | Pearson r vs mature MSSM OLS | Both-FDR same/tested |
|---|---|---|---:|---:|---:|
| VIP | primary_wls | tested | 15 | 0.0028 | 1/4 |
| VIP | ols_same_design | tested | 15 | 0.0308 | 1/4 |
| VIP | wls_without_adnc | tested | 15 | -0.0090 | 1/4 |
| VIP | singleome_wls | tested | 15 | 0.0173 | 1/4 |
| VIP | floor_50_nuclei_wls | tested | 15 | 0.0028 | 1/4 |
| VIP | not_ad_only_wls | untestable | 8 | UNTESTABLE | UNTESTABLE |
| pvalb | primary_wls | tested | 15 | -0.1567 | 2/7 |
| pvalb | ols_same_design | tested | 15 | -0.1686 | 2/7 |
| pvalb | wls_without_adnc | tested | 15 | -0.1627 | 2/7 |
| pvalb | singleome_wls | tested | 15 | -0.1264 | 2/7 |
| pvalb | floor_50_nuclei_wls | tested | 15 | -0.1567 | 2/7 |
| pvalb | not_ad_only_wls | untestable | 8 | UNTESTABLE | UNTESTABLE |
| sst | primary_wls | tested | 15 | -0.1692 | 16/38 |
| sst | ols_same_design | tested | 15 | 0.0478 | 14/38 |
| sst | wls_without_adnc | tested | 15 | -0.1670 | 16/38 |
| sst | singleome_wls | tested | 15 | -0.0707 | 21/38 |
| sst | floor_50_nuclei_wls | tested | 15 | -0.1692 | 16/38 |
| sst | not_ad_only_wls | untestable | 8 | UNTESTABLE | UNTESTABLE |

| Subtype | Successful drops | WLS correlation range | Sign flips | Both-FDR fraction range | Gene sign flips/evaluated |
|---|---:|---:|---:|---:|---:|
| VIP | 15 | -0.1017 to 0.0676 | 6 | 0.250 to 0.500 | 2/4 |
| pvalb | 15 | -0.2614 to -0.0710 | 0 | 0.286 to 0.714 | 3/7 |
| sst | 15 | -0.3481 to 0.2362 | 3 | 0.395 to 0.553 | 18/38 |

This is a limited, model-sensitive result in 15 old donors, not evidence of no aging. The narrow age span, old-age range, MTG region, assay, donor count, residual pathology, and estimator sensitivity limit generalization. Not-AD-only fits remain untestable below the 10-donor floor. Mouse P1 sets are too sparse for a cross-species verdict when only a few genes pass the SEA-AD filter.

## Verification and deviations

The new complete MSSM/HBCC model families and SEA-AD gene effects were independently BH-recomputed from raw p-values. Maximum absolute q differences were 1.149e-13 for MSSM/HBCC and 0.000e+00 for SEA-AD, both below 1e-9. Selected replication-list q values retain their full-unit BH families; BH was not reapplied after selecting genes. Exact per-cohort, per-subtype family values are in the two BH audit CSVs.

| Model or correction family | Maximum absolute q difference |
|---|---:|
| MSSM/HBCC rev2_voom_wls | 1.149e-13 |
| MSSM/HBCC age_plus_sex_depth_ols | 3.536e-14 |
| MSSM/HBCC age_plus_sex_ols_depth_unadjusted | 1.589e-14 |
| MSSM/HBCC rev2_voom_wls_moderated | 0.000e+00 |
| SEA-AD primary_voom_wls | 0.000e+00 |
| SEA-AD ols_same_design | 0.000e+00 |
| SEA-AD genome_wide_block_permutation | 0.000e+00 |
| SEA-AD mature_both_fdr_block_permutation | 0.000e+00 |

The frozen moderation recovery gate used every frozen unit and the exact recovered-s2 formula. It requires both relative differences to be at most 5 percent. The gate passed: False. Frozen gene files contain selected genes, so the recovered prior fit uses the available selected-gene rows rather than a complete tested-gene vector. This limitation and the observed discrepancies are reported rather than hidden.

| Cohort | Subtype | Recovered genes | Frozen s0_sq | Recovered s0_sq | Relative difference | Frozen d0 | Recovered d0 | Relative difference |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| MSSM_DLPFC | VIP | 202 | 1.413055 | 0.960124 | 32.05% | 4.053685 | 5.074197 | 25.17% |
| MSSM_DLPFC | pvalb | 137 | 1.458388 | 0.958068 | 34.31% | 3.949667 | 5.385622 | 36.36% |
| MSSM_DLPFC | sst | 362 | 1.404325 | 0.997894 | 28.94% | 4.326758 | 4.722075 | 9.14% |
| MSSM_DLPFC | GABA | 65 | 1.309581 | 0.797468 | 39.11% | 4.743125 | 4.931786 | 3.98% |
| HBCC_DLPFC | VIP | 202 | 1.672292 | 1.024637 | 38.73% | 3.364052 | 5.161732 | 53.44% |
| HBCC_DLPFC | pvalb | 137 | 1.706020 | 0.905314 | 46.93% | 3.195126 | 4.166436 | 30.40% |
| HBCC_DLPFC | sst | 362 | 1.414791 | 0.878628 | 37.90% | 3.339865 | 3.753006 | 12.37% |
| HBCC_DLPFC | GABA | 65 | 1.553144 | 0.887372 | 42.87% | 3.869196 | 4.817075 | 24.50% |

The WLS copy returning final weighted residual s2 reproduced cached slopes, SEs, p values, and q values to 1e-9 before mature moderation. Mature priors were re-estimated separately for every mature unit. Because the gate failed, moderated mature intersections and six-gene statuses are unverified sensitivities; the unmoderated intersection is used for the confirmed result. WLS bootstrap intervals were added for the required SEA-AD report while the original primary OLS analysis was unchanged. No frozen original was recomputed or overwritten, no donor object was downloaded or re-hashed, and no original SEA-AD output was modified.
