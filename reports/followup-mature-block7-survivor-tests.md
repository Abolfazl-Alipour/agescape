# Mature-to-Old (30-plus) Block 7: Survivor-Bias Within-Type Tests

## Status and scope

Phase 2 Item 13 of the Agescape mature-to-old rebuild. This report reruns 11pt Block 7 strictly on human mature donors aged 30 years and older (`age_years >= 30`), reusing the frozen protocol, frozen 271-gene death/stress scoring list, and frozen verdict rules from `reports/followup-11pt-block7-survivor-within-type.md` verbatim. Mature signature P1m does not enter Block 7 because Block 7 scores the frozen Hallmark apoptosis plus unfolded protein response gene list.

## Frozen protocol as applied at 30-plus

### Frozen death/stress gene list and citation

Source: Liberzon et al. 2015, Cell Systems, The Molecular Signatures Database Hallmark Gene Set Collection, DOI `10.1016/j.cels.2015.12.004`, PMID `26771021`. Official MSigDB sets: `HALLMARK_APOPTOSIS` (systematic ID `M5902`, 161 genes) and `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` (systematic ID `M5922`, 113 genes). Their union is 272 unique symbols (`ATF3` and `DNAJC3` occur in both). The shared human technical mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. `RPS14` is the sole excluded source-list gene, leaving 271 fixed scoring genes. Absent dataset symbols contribute zero counts and remain in the fixed denominator of 271 (Gao has 248 present and 23 absent; Yazar has 271 present and 0 absent).

Apoptosis source genes (161): ADD1, AIFM3, ANKH, ANXA1, APP, ATF3, AVPR1A, BAX, BCAP31, BCL10, BCL2L1, BCL2L10, BCL2L11, BCL2L2, BGN, BID, BIK, BIRC3, BMF, BMP2, BNIP3L, BRCA1, BTG2, BTG3, CASP1, CASP2, CASP3, CASP4, CASP6, CASP7, CASP8, CASP9, CAV1, CCNA1, CCND1, CCND2, CD14, CD2, CD38, CD44, CD69, CDC25B, CDK2, CDKN1A, CDKN1B, CFLAR, CLU, CREBBP, CTH, CTNNB1, CYLD, DAP, DAP3, DCN, DDIT3, DFFA, DIABLO, DNAJA1, DNAJC3, DNM1L, DPYD, EBP, EGR3, EMP1, ENO2, ERBB2, ERBB3, EREG, ETF1, F2, F2R, FAS, FASLG, FDXR, FEZ1, GADD45A, GADD45B, GCH1, GNA15, GPX1, GPX3, GPX4, GSN, GSR, GSTM1, GUCY2D, H1-0, HGF, HMGB2, HMOX1, HSPB1, IER3, IFITM3, IFNB1, IFNGR1, IGF2R, IGFBP6, IL18, IL1A, IL1B, IL6, IRF1, ISG20, JUN, KRT18, LEF1, LGALS3, LMNA, LUM, MADD, MCL1, MGMT, MMP2, NEDD9, NEFH, PAK1, PDCD4, PDGFRB, PEA15, PLAT, PLCB2, PLPPR4, PMAIP1, PPP2R5B, PPP3R1, PPT1, PRF1, PSEN1, PSEN2, PTK2, RARA, RELA, RETSAT, RHOB, RHOT2, RNASEL, ROCK1, SAT1, SATB1, SC5D, SLC20A1, SMAD7, SOD1, SOD2, SPTAN1, SQSTM1, TAP1, TGFB2, TGFBR3, TIMP1, TIMP2, TIMP3, TNF, TNFRSF12A, TNFSF10, TOP2A, TSPO, TXNIP, VDAC2, WEE1, XIAP.

Unfolded protein response source genes (113): ALDH18A1, ARFGAP1, ASNS, ATF3, ATF4, ATF6, ATP6V0D1, BAG3, BANF1, CALR, CCL2, CEBPB, CEBPG, CHAC1, CKS1B, CNOT2, CNOT4, CNOT6, CXXC1, DCP1A, DCP2, DCTN1, DDIT4, DDX10, DKC1, DNAJA4, DNAJB9, DNAJC3, EDC4, EDEM1, EEF2, EIF2AK3, EIF2S1, EIF4A1, EIF4A2, EIF4A3, EIF4E, EIF4EBP1, EIF4G1, ERN1, ERO1A, EXOC2, EXOSC1, EXOSC10, EXOSC2, EXOSC4, EXOSC5, EXOSC9, FKBP14, FUS, GEMIN4, GOSR2, H2AX, HERPUD1, HSP90B1, HSPA5, HSPA9, HYOU1, IARS1, IFIT1, IGFBP1, IMP3, KDELR3, KHSRP, KIF5B, LSM1, LSM4, MTHFD2, MTREX, NABP1, NFYA, NFYB, NHP2, NOLC1, NOP14, NOP56, NPM1, PAIP1, PARN, PDIA5, PDIA6, POP4, PREB, PSAT1, RPS14, RRP9, SDAD1, SEC11A, SEC31A, SERP1, SHC1, SKIC3, SLC1A4, SLC30A5, SLC7A5, SPCS1, SPCS3, SRPRA, SRPRB, SSR1, STC2, TARS1, TATDN2, TSPYL2, TUBB2A, VEGFA, WFS1, WIPI1, XBP1, XPOT, YIF1A, YWHAZ, ZBTB17.

Final scoring symbols (271): ADD1, AIFM3, ALDH18A1, ANKH, ANXA1, APP, ARFGAP1, ASNS, ATF3, ATF4, ATF6, ATP6V0D1, AVPR1A, BAG3, BANF1, BAX, BCAP31, BCL10, BCL2L1, BCL2L10, BCL2L11, BCL2L2, BGN, BID, BIK, BIRC3, BMF, BMP2, BNIP3L, BRCA1, BTG2, BTG3, CALR, CASP1, CASP2, CASP3, CASP4, CASP6, CASP7, CASP8, CASP9, CAV1, CCL2, CCNA1, CCND1, CCND2, CD14, CD2, CD38, CD44, CD69, CDC25B, CDK2, CDKN1A, CDKN1B, CEBPB, CEBPG, CFLAR, CHAC1, CKS1B, CLU, CNOT2, CNOT4, CNOT6, CREBBP, CTH, CTNNB1, CXXC1, CYLD, DAP, DAP3, DCN, DCP1A, DCP2, DCTN1, DDIT3, DDIT4, DDX10, DFFA, DIABLO, DKC1, DNAJA1, DNAJA4, DNAJB9, DNAJC3, DNM1L, DPYD, EBP, EDC4, EDEM1, EEF2, EGR3, EIF2AK3, EIF2S1, EIF4A1, EIF4A2, EIF4A3, EIF4E, EIF4EBP1, EIF4G1, EMP1, ENO2, ERBB2, ERBB3, EREG, ERN1, ERO1A, ETF1, EXOC2, EXOSC1, EXOSC10, EXOSC2, EXOSC4, EXOSC5, EXOSC9, F2, F2R, FAS, FASLG, FDXR, FEZ1, FKBP14, FUS, GADD45A, GADD45B, GCH1, GEMIN4, GNA15, GOSR2, GPX1, GPX3, GPX4, GSN, GSR, GSTM1, GUCY2D, H1-0, H2AX, HERPUD1, HGF, HMGB2, HMOX1, HSP90B1, HSPA5, HSPA9, HSPB1, HYOU1, IARS1, IER3, IFIT1, IFITM3, IFNB1, IFNGR1, IGF2R, IGFBP1, IGFBP6, IL18, IL1A, IL1B, IL6, IMP3, IRF1, ISG20, JUN, KDELR3, KHSRP, KIF5B, KRT18, LEF1, LGALS3, LMNA, LSM1, LSM4, LUM, MADD, MCL1, MGMT, MMP2, MTHFD2, MTREX, NABP1, NEDD9, NEFH, NFYA, NFYB, NHP2, NOLC1, NOP14, NOP56, NPM1, PAIP1, PAK1, PARN, PDCD4, PDGFRB, PDIA5, PDIA6, PEA15, PLAT, PLCB2, PLPPR4, PMAIP1, POP4, PPP2R5B, PPP3R1, PPT1, PREB, PRF1, PSAT1, PSEN1, PSEN2, PTK2, RARA, RELA, RETSAT, RHOB, RHOT2, RNASEL, ROCK1, RRP9, SAT1, SATB1, SC5D, SDAD1, SEC11A, SEC31A, SERP1, SHC1, SKIC3, SLC1A4, SLC20A1, SLC30A5, SLC7A5, SMAD7, SOD1, SOD2, SPCS1, SPCS3, SPTAN1, SQSTM1, SRPRA, SRPRB, SSR1, STC2, TAP1, TARS1, TATDN2, TGFB2, TGFBR3, TIMP1, TIMP2, TIMP3, TNF, TNFRSF12A, TNFSF10, TOP2A, TSPO, TSPYL2, TUBB2A, TXNIP, VDAC2, VEGFA, WEE1, WFS1, WIPI1, XBP1, XIAP, XPOT, YIF1A, YWHAZ, ZBTB17.

### Part A specification (restricted to age 30-plus)

- Per-cell program score is the mean within-cell proportion across the fixed 271 scoring genes: sum of post-mask raw counts divided by `271 * post-mask UMI`. Primary analyses use full-depth counts without downsampling. Cells with post-mask UMI below 1500 are excluded as an eligibility filter only.
- All steps are restricted strictly to donors aged 30 years and older (`age_years >= 30.0`): cell eligibility, donor-cell eligibility (at least 20 eligible cells per mature donor per cell type), the per-cell-type 90th-percentile high-score threshold (pooled across eligible mature donor-cell-type pairs only, numpy linear quantile interpolation, high defined as `score >= threshold` including ties), and every regression. No full-range or young-derived threshold is reused.
- Test 1 outcome is each mature donor's fraction of high-scoring cells among that donor's eligible cells in the cell type. Primary OLS fits `fraction_high ~ age_per_decade + sex_male + batch` (Gao 7-study fixed effects with `ENCODE v4 (Snyder)` reference; Yazar 75-pool fixed effects) requiring at least 6 mature donors, singleton-batch pruning (`min_batch_donors = 2`), rank-checked covariate pruning, and residual df at least 5 (Gao) or 10 (Yazar). Unadjusted sensitivity retains `age_per_decade + sex_male` without batch (never substituted for a failed batch-adjusted fit). Models without sufficient residual df or donor/batch floors are recorded as UNTESTABLE, never null.
- Labeled 1500 sensitivity repeats test 1 only using exact marginal `Binomial(original_list_count, 1500 / original_UMI)` draws from `SeedSequence([20261117, 7, dataset_code, donor_idx])` with mature cells sorted by original h5ad row index, and recomputes the 90th-percentile threshold on mature rarefied scores only.
- Test 1 BH families are each `(dataset, test)` slice across evaluable non-pooled cell types at 30-plus, separately for `test1_primary_full`, `test1_unadjusted_full`, and `test1_rarefied_1500`. Every BH family was independently recomputed using a manual step-up implementation and cross-checked against `statsmodels.stats.multitest.multipletests`.
- Test 2 runs only for cell types passing primary test 1 BH `q < 0.05` at 30-plus. Unremoved and removed pseudobulk gene OLS models (`log2CPM ~ age_per_decade + sex_male + log_lib + batch`) use the same common mature donor set (at least 20 cells per donor in both arms, singleton batches pruned, at least 6 common mature donors), gene expression filter `CPM >= 1` in at least half of donors (minimum 3 donors), and within-arm gene BH (`q_bh_within_group < 0.05`). A predicted pass-to-fail change requires an unremoved pass becoming a removed fail with `|slope_removed| < |slope_unremoved|`. Slope change is `removed - unremoved`. Removing captured high-scoring cells is a lower bound because dying cells are under-captured.
- MDE80 is `(t_0.975 + t_0.80) * SE` using residual df, reported on both the raw outcome scale and standardized by the donor-outcome SD (`mde80_sd_per_decade`). Small MDE80 is predefined as at most `0.50` outcome SD per decade.

### Part B specification (mature slopes only)

- Joins ONLY mature-recomputed per-cell-type batch-adjusted CLR fraction slopes per decade and mature-recomputed full-depth primary Shannon spread slopes per decade (`metric=shannon`, `depth_mode=full_depth`, `analysis_role=primary_full_depth`, `model_spec=primary`, non-pooled, non-stratified-global). Legacy full-range slopes are never joined.
- Fraction slope negative means losing share; spread slope is Shannon change per decade (negative means narrowing). Frozen prediction: survivor bias predicts a positive cross-type association between fraction slope and spread slope, with more negative spread slopes in fraction-losing cell types.
- Weighted least squares (`spread_slope ~ 1 + fraction_slope`) with inverse-variance weights `1 / SE_spread^2`, fitted separately per dataset, one test per dataset, no multiplicity correction across the two datasets, never pooled. Zero fraction slope counts as gaining.

### Frozen verdict rules

- Part A per dataset: `PROVES` only if test 1 passes BH correction (`q < 0.05`) and test 2 produces at least one predicted pass-to-fail gene change with absolute-slope attenuation. `DOES NOT PROVE` if no test 1 result passes BH and every evaluable primary null has MDE80 at most 0.50 SD per decade. `INCONCLUSIVE` if any primary null has MDE80 above 0.50 SD per decade, or test 1 is significant but removal produces zero predicted changes, or the required removal comparison is not evaluable. Overall Part A is `PROVES` if either dataset meets `PROVES`, `DOES NOT PROVE` only if both meet `DOES NOT PROVE`, and otherwise `INCONCLUSIVE`.
- Part B: `PROVES` only if the cross-type association is positive and significant (`p < 0.05`) in at least one dataset and the other dataset is not significantly opposite. `DOES NOT PROVE` if the association is null or opposite in both datasets. `INCONCLUSIVE` if one is significantly positive and the other significantly opposite, or if a dataset is UNTESTABLE.
- Neither Part A nor Part B alone proves survivor bias as a mechanism. The Nov 17 cross-type correlation is not rerun here.

## Input provenance and mature cohort audit

| Input | Path | Provenance / SHA-256 prefix | Verification |
| --- | --- | --- | --- |
| Gao counts cache (Block 1) | `/home/reza/projects/agescape-revision-data/tmp_organs/block1_cache/gao` | Reused Block 1 donor-sharded post-mask UMI >= 1500 cache | 109 of 117 donors >= 30 (8 excluded; ages 30 to 91; 482902 mature cells scored, 39859 under-30 cells excluded) |
| Yazar counts cache (Block 1) | `/home/reza/projects/agescape-revision-data/tmp_organs/block1_cache/yazar` | Reused Block 1 donor-sharded post-mask UMI >= 1500 cache | 932 of 981 donors >= 30 (49 excluded; ages 30 to 97; 592250 mature cells scored, 25510 under-30 cells excluded) |
| Core mature human composition | `results/followup-mature-composition-human.csv` | committed c0a034a (2026-11-21); sha256 `250179f836a2` | Checked first: 285 rows across 11 core brain/kidney/lung cohorts; 0 Gao and 0 Yazar rows (not joined) |
| Gao heart mature spread | `results/followup-mature-fulldepth-gao-spread.csv` | untracked working-tree output from parallel Phase 2 runner; sha256 `8386cf72f383` | Verified mature primary full-depth Shannon rows: 13, max N donors = 101 |
| Gao heart mature CLR fractions | `results/followup-mature-organ-transports-a-gao.csv` | untracked working-tree output from parallel Phase 2 runner; sha256 `ba6b107d2aa6` | Verified mature `fractions_clr` (`arm=pooled_7study_primary`): 13 cell types, N = 109 mature donors (8 under-30 excluded) |
| Yazar blood mature spread | `results/followup-mature-fulldepth-yazar-spread.csv` | untracked working-tree output from parallel Phase 2 runner; sha256 `c112de004a48` | Verified mature primary full-depth Shannon rows: 17, max N donors = 909 |
| Yazar blood mature CLR fractions | `results/followup-mature-organ-transports-a-yazar.csv` | untracked working-tree output from parallel Phase 2 runner; sha256 `7a191ff66c5b` | Verified mature `fractions_clr` (`arm=pool_fe`): 28 cell types, N = 932 mature donors (49 under-30 excluded) |

## Results

### Part A test 1: high death/stress-score cell fractions at age 30-plus

#### Gao (primary full depth, batch-adjusted)

Evaluated 13 of 13 cell types (0 UNTESTABLE). BH family size 13; 2 passed BH q < 0.05. Evaluated primary nulls with MDE80 > 0.50 SD/decade: 0 (max MDE80 = 0.411 SD/decade).

| Cell type | Status | N orig | N mature | Resid df | Slope per decade | 95 percent CI | p | BH q | MDE80 SD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Adipocyte | EVALUATED | 37 | 35 | 28 | +0.029407 | -0.001347 to 0.060160 | 6.018e-02 | 2.117e-01 | 0.365 |
| Cardiomyocyte | EVALUATED | 108 | 101 | 93 | -0.012452 | -0.029095 to 0.004192 | 1.408e-01 | 3.050e-01 | 0.223 |
| Endocardial | EVALUATED | 72 | 66 | 58 | -0.012678 | -0.030910 to 0.005554 | 1.692e-01 | 3.143e-01 | 0.308 |
| Endothelial | EVALUATED | 86 | 79 | 72 | -0.011547 | -0.030710 to 0.007615 | 2.336e-01 | 3.631e-01 | 0.275 |
| Epicardial | EVALUATED | 19 | 17 | 13 | -0.024718 | -0.072942 to 0.023506 | 2.882e-01 | 3.631e-01 | 0.353 |
| Fibroblast | EVALUATED | 105 | 98 | 90 | -0.011095 | -0.032022 to 0.009831 | 2.950e-01 | 3.631e-01 | 0.225 |
| LEC | EVALUATED | 36 | 34 | 29 | -0.012980 | -0.049712 to 0.023753 | 4.756e-01 | 4.756e-01 | 0.411 |
| Lymphoid | EVALUATED | 59 | 56 | 49 | -0.029237 | -0.047625 to -0.010849 | 2.446e-03 | 3.179e-02 | 0.259 |
| Mast | EVALUATED | 17 | 16 | 10 | -0.018086 | -0.072315 to 0.036143 | 4.745e-01 | 4.756e-01 | 0.401 |
| Myeloid | EVALUATED | 106 | 99 | 91 | -0.006494 | -0.019058 to 0.006070 | 3.072e-01 | 3.631e-01 | 0.257 |
| Neuronal | EVALUATED | 60 | 56 | 49 | -0.021248 | -0.045248 to 0.002751 | 8.141e-02 | 2.117e-01 | 0.304 |
| Pericyte | EVALUATED | 95 | 88 | 81 | -0.024178 | -0.041156 to -0.007200 | 5.811e-03 | 3.777e-02 | 0.209 |
| vSMC | EVALUATED | 86 | 79 | 72 | -0.015068 | -0.031204 to 0.001067 | 6.674e-02 | 2.117e-01 | 0.227 |

- `gao test1_rarefied_1500`: 13 evaluated of 13, 2 BH hit(s): Lymphoid (slope -0.030565, p=4.360e-03, q=4.369e-02), Pericyte (slope -0.023265, p=6.721e-03, q=4.369e-02).
- `gao test1_unadjusted_full`: 13 evaluated of 13, 0 BH hit(s).

#### Yazar (primary full depth, batch-adjusted)

Evaluated 16 of 28 cell types (12 UNTESTABLE). BH family size 16; 3 passed BH q < 0.05. Evaluated primary nulls with MDE80 > 0.50 SD/decade: 1 (max MDE80 = 0.653 SD/decade).

| Cell type | Status | N orig | N mature | Resid df | Slope per decade | 95 percent CI | p | BH q | MDE80 SD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | EVALUATED | 252 | 245 | 192 | -0.006669 | -0.017576 to 0.004237 | 2.293e-01 | 6.114e-01 | 0.167 |
| CD14-positive monocyte | EVALUATED | 423 | 402 | 342 | -0.001013 | -0.006545 to 0.004519 | 7.190e-01 | 9.684e-01 | 0.105 |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_FLOOR | 6 | 5 | NA | NA | NA to NA | NA | NA | NA |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| CD4-positive, alpha-beta cytotoxic T cell | EVALUATED | 75 | 74 | 48 | +0.006590 | -0.013861 to 0.027040 | 5.202e-01 | 9.324e-01 | 0.303 |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| central memory CD4-positive, alpha-beta T cell | EVALUATED | 955 | 909 | 832 | +0.006713 | 0.003744 to 0.009683 | 1.034e-05 | 8.269e-05 | 0.051 |
| central memory CD8-positive, alpha-beta T cell | EVALUATED | 91 | 85 | 58 | +0.010165 | -0.016855 to 0.037185 | 4.545e-01 | 9.324e-01 | 0.329 |
| conventional dendritic cell | UNTESTABLE_DONOR_FLOOR | 5 | 5 | NA | NA | NA to NA | NA | NA | NA |
| dendritic cell | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| double negative thymocyte | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| effector memory CD4-positive, alpha-beta T cell | EVALUATED | 261 | 250 | 199 | +0.000206 | -0.009233 to 0.009645 | 9.657e-01 | 9.684e-01 | 0.128 |
| effector memory CD8-positive, alpha-beta T cell | EVALUATED | 746 | 721 | 645 | +0.001256 | -0.003231 to 0.005743 | 5.827e-01 | 9.324e-01 | 0.067 |
| erythrocyte | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| gamma-delta T cell | EVALUATED | 120 | 111 | 75 | +0.000326 | -0.015996 to 0.016648 | 9.684e-01 | 9.684e-01 | 0.207 |
| hematopoietic precursor cell | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| innate lymphoid cell | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| memory B cell | EVALUATED | 244 | 227 | 178 | +0.001618 | -0.008255 to 0.011492 | 7.467e-01 | 9.684e-01 | 0.132 |
| mucosal-associated invariant T cell | EVALUATED | 30 | 25 | 14 | +0.006241 | -0.056242 to 0.068724 | 8.335e-01 | 9.684e-01 | 0.653 |
| naive B cell | EVALUATED | 319 | 307 | 265 | +0.009607 | 0.001324 to 0.017890 | 2.318e-02 | 9.272e-02 | 0.108 |
| naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 821 | 777 | 701 | +0.008188 | 0.004576 to 0.011800 | 9.947e-06 | 8.269e-05 | 0.064 |
| naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 384 | 344 | 284 | +0.012181 | 0.005405 to 0.018956 | 4.696e-04 | 2.504e-03 | 0.101 |
| natural killer cell | EVALUATED | 815 | 783 | 707 | +0.002343 | -0.001013 to 0.005699 | 1.710e-01 | 5.471e-01 | 0.061 |
| plasmablast | UNTESTABLE_NOT_EVALUABLE_WITHIN_BATCH_DF | 14 | 14 | NA | NA | NA to NA | NA | NA | NA |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_FLOOR | 3 | 3 | NA | NA | NA to NA | NA | NA | NA |
| platelet | UNTESTABLE_DONOR_FLOOR | 0 | 0 | NA | NA | NA to NA | NA | NA | NA |
| regulatory T cell | EVALUATED | 54 | 54 | 36 | +0.006851 | -0.015423 to 0.029126 | 5.367e-01 | 9.324e-01 | 0.484 |
| transitional stage B cell | EVALUATED | 201 | 185 | 139 | +0.000444 | -0.011782 to 0.012669 | 9.429e-01 | 9.684e-01 | 0.147 |

UNTESTABLE primary cell types in Yazar:

| Cell type | Status / Reason | N mature donors (>= 20 cells) | N mature donors with >= 1 cell | Covariates omitted |
| --- | --- | --- | --- | --- |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_FLOOR | 5 | 777 | none |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_FLOOR | 0 | 414 | none |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_FLOOR | 0 | 194 | none |
| conventional dendritic cell | UNTESTABLE_DONOR_FLOOR | 5 | 829 | none |
| dendritic cell | UNTESTABLE_DONOR_FLOOR | 0 | 126 | none |
| double negative thymocyte | UNTESTABLE_DONOR_FLOOR | 0 | 525 | none |
| erythrocyte | UNTESTABLE_DONOR_FLOOR | 0 | 144 | none |
| hematopoietic precursor cell | UNTESTABLE_DONOR_FLOOR | 0 | 640 | none |
| innate lymphoid cell | UNTESTABLE_DONOR_FLOOR | 0 | 140 | none |
| plasmablast | UNTESTABLE_NOT_EVALUABLE_WITHIN_BATCH_DF | 14 | 748 | batch_block(df_limit) |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_FLOOR | 3 | 589 | none |
| platelet | UNTESTABLE_DONOR_FLOOR | 0 | 432 | none |

- `yazar test1_rarefied_1500`: 16 evaluated of 28, 2 BH hit(s): naive thymus-derived CD4-positive, alpha-beta T cell (slope +0.009357, p=5.670e-07, q=9.072e-06), central memory CD4-positive, alpha-beta T cell (slope +0.006310, p=1.723e-05, q=1.378e-04).
- `yazar test1_unadjusted_full`: 17 evaluated of 28, 2 BH hit(s): naive thymus-derived CD4-positive, alpha-beta T cell (slope +0.007493, p=2.184e-04, q=3.713e-03), central memory CD4-positive, alpha-beta T cell (slope +0.005669, p=3.850e-03, q=3.273e-02).

### Independent BH recomputation verification

Every BH family in Part A test 1 and test 2 was independently recomputed from raw p-values using a standalone step-up implementation (`p_sorted * m / rank`, reverse cumulative minimum, clipped to `[0, 1]`) and checked against `statsmodels.stats.multitest.multipletests(..., method='fdr_bh')`.

| Stage | Dataset | BH family | Family size | Passing q < 0.05 | Min raw p | Min BH q | Max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| test1 | gao | test1_primary_full | 13 | 2 | 2.446e-03 | 3.179e-02 | 5.551e-17 |
| test1 | gao | test1_rarefied_1500 | 13 | 2 | 4.360e-03 | 4.369e-02 | 5.551e-17 |
| test1 | gao | test1_unadjusted_full | 13 | 0 | 2.871e-02 | 3.732e-01 | 5.551e-17 |
| test1 | yazar | test1_primary_full | 16 | 3 | 9.947e-06 | 8.269e-05 | 0.000e+00 |
| test1 | yazar | test1_rarefied_1500 | 16 | 2 | 5.670e-07 | 9.072e-06 | 0.000e+00 |
| test1 | yazar | test1_unadjusted_full | 17 | 2 | 2.184e-04 | 3.713e-03 | 1.110e-16 |
| test2 | gao | Lymphoid (removed) | 10551 | 0 | 3.960e-05 | 4.178e-01 | 2.220e-16 |
| test2 | gao | Lymphoid (unremoved) | 10775 | 0 | 1.846e-05 | 1.989e-01 | 2.220e-16 |
| test2 | gao | Pericyte (removed) | 12007 | 1 | 3.353e-06 | 4.026e-02 | 2.220e-16 |
| test2 | gao | Pericyte (unremoved) | 12076 | 1 | 2.258e-06 | 2.726e-02 | 2.220e-16 |
| test2 | yazar | central memory CD4-positive, alpha-beta T cell (removed) | 10800 | 115 | 8.463e-12 | 9.140e-08 | 2.220e-16 |
| test2 | yazar | central memory CD4-positive, alpha-beta T cell (unremoved) | 10971 | 148 | 1.425e-13 | 1.563e-09 | 2.220e-16 |
| test2 | yazar | naive thymus-derived CD4-positive, alpha-beta T cell (removed) | 10088 | 135 | 1.780e-15 | 1.796e-11 | 2.220e-16 |
| test2 | yazar | naive thymus-derived CD4-positive, alpha-beta T cell (unremoved) | 10319 | 194 | 6.519e-17 | 6.727e-13 | 2.220e-16 |
| test2 | yazar | naive thymus-derived CD8-positive, alpha-beta T cell (removed) | 8488 | 0 | 1.459e-05 | 8.561e-02 | 2.220e-16 |
| test2 | yazar | naive thymus-derived CD8-positive, alpha-beta T cell (unremoved) | 8746 | 4 | 3.890e-07 | 3.402e-03 | 2.220e-16 |

### Part A test 2: removal sensitivity at age 30-plus

| Dataset | Cell type | Status | Common mature donors | Unremoved cells | Removed cells | High cells removed | Common genes | Unremoved pass | Removed pass | Predicted pass to fail | Other pass to fail | Fail to pass | Pass to pass | Fail to fail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Lymphoid | EVALUATED | 54 | 6082 | 5483 | 599 | 10549 | 0 | 0 | 0 | 0 | 0 | 0 | 10549 |
| gao | Pericyte | EVALUATED | 87 | 30395 | 27366 | 3029 | 11987 | 1 | 1 | 0 | 0 | 0 | 1 | 11986 |
| yazar | central memory CD4-positive, alpha-beta T cell | EVALUATED | 902 | 139734 | 125778 | 13956 | 10800 | 146 | 115 | 37 | 4 | 10 | 105 | 10644 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 749 | 105941 | 95405 | 10536 | 10088 | 190 | 135 | 56 | 9 | 10 | 125 | 9888 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 310 | 17101 | 15486 | 1615 | 8488 | 4 | 0 | 3 | 1 | 0 | 0 | 8484 |

- yazar `central memory CD4-positive, alpha-beta T cell`: 37 predicted pass-to-fail genes with absolute-slope attenuation (median |slope| reduction = 0.0043 log2CPM/decade; top attenuated genes: PAQR8, CEP164, UVSSA, ENSG00000247121, ENSG00000274092).
- yazar `naive thymus-derived CD4-positive, alpha-beta T cell`: 56 predicted pass-to-fail genes with absolute-slope attenuation (median |slope| reduction = 0.0108 log2CPM/decade; top attenuated genes: PPP1CB-DT, ENSG00000273445, FAM13A, MYB, ARHGAP35).
- yazar `naive thymus-derived CD8-positive, alpha-beta T cell`: 3 predicted pass-to-fail genes with absolute-slope attenuation (median |slope| reduction = 0.0293 log2CPM/decade; top attenuated genes: RASGRF2, ANXA1, NR4A2).

Pass/fail change counts use genes evaluated in both arms (`EVALUATED_BOTH`) on the exact same common mature donor set. Slope change is `removed - unremoved`. Caveat: removing captured high-scoring cells provides a lower bound on stress/dying-cell contribution because dying cells are under-captured in droplet scRNA-seq/snRNA-seq.

### Part B: narrowing localization (mature slopes only)

| Dataset | Status | N types | Association slope | 95 percent CI | p | MDE80 | Losing N | Gaining N | Mean spread losing | Mean spread gaining | Reason / Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | EVALUATED | 13 | -0.00537 | -0.12660 to 0.11586 | 9.241e-01 | 0.16946 | 7 | 6 | -0.01559 | -0.01287 | evaluated on mature slopes |
| yazar | EVALUATED | 17 | -0.01384 | -0.03194 to 0.00427 | 1.241e-01 | 0.02546 | 10 | 7 | -0.00690 | -0.00262 | evaluated on mature slopes |

The association slope is the change in mature Shannon spread slope per unit mature CLR fraction slope per decade (WLS with weights `1 / SE_spread^2`). Fraction slope negative means losing share; zero fraction slope counts as gaining. Frozen prediction: positive association. One test per dataset, no multiplicity correction across datasets, never pooled.

#### Part B matched per-cell-type mature slopes

| Dataset | Cell type | Spread N | Mature spread slope | Spread SE | Fraction N | Mature CLR fraction slope | Fraction SE | Direction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 35 | -0.00547 | 0.02494 | 109 | -0.10640 | 0.13039 | losing |
| gao | Cardiomyocyte | 101 | -0.01003 | 0.02042 | 109 | -0.02697 | 0.06165 | losing |
| gao | Endocardial | 66 | -0.01466 | 0.00875 | 109 | -0.01393 | 0.12347 | losing |
| gao | Endothelial | 79 | -0.00471 | 0.00586 | 109 | -0.07854 | 0.08771 | losing |
| gao | Epicardial | 17 | +0.02586 | 0.03154 | 109 | +0.24603 | 0.13406 | gaining |
| gao | Fibroblast | 98 | -0.01157 | 0.00833 | 109 | +0.01851 | 0.03252 | gaining |
| gao | LEC | 34 | -0.00106 | 0.01014 | 109 | +0.07530 | 0.08666 | gaining |
| gao | Lymphoid | 56 | -0.03530 | 0.01699 | 109 | +0.14289 | 0.06020 | gaining |
| gao | Mast | 16 | -0.01644 | 0.05109 | 109 | -0.18524 | 0.09796 | losing |
| gao | Myeloid | 99 | -0.03450 | 0.01015 | 109 | +0.00514 | 0.03639 | gaining |
| gao | Neuronal | 56 | -0.02703 | 0.01101 | 109 | -0.07931 | 0.05311 | losing |
| gao | Pericyte | 88 | -0.03076 | 0.01255 | 109 | -0.04405 | 0.03887 | losing |
| gao | vSMC | 79 | -0.02065 | 0.01774 | 109 | +0.04658 | 0.03853 | gaining |
| yazar | CD14-low, CD16-positive monocyte | 245 | -0.00760 | 0.00232 | 932 | +0.14529 | 0.02535 | gaining |
| yazar | CD14-positive monocyte | 402 | -0.00616 | 0.00232 | 932 | +0.09453 | 0.02380 | gaining |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 74 | +0.00104 | 0.00933 | 932 | +0.15503 | 0.02826 | gaining |
| yazar | central memory CD4-positive, alpha-beta T cell | 909 | +0.00110 | 0.00128 | 932 | +0.02330 | 0.00986 | gaining |
| yazar | central memory CD8-positive, alpha-beta T cell | 85 | +0.00893 | 0.00733 | 932 | -0.08123 | 0.01791 | losing |
| yazar | effector memory CD4-positive, alpha-beta T cell | 250 | -0.00498 | 0.00388 | 932 | -0.04672 | 0.01482 | losing |
| yazar | effector memory CD8-positive, alpha-beta T cell | 721 | -0.00187 | 0.00166 | 932 | +0.09465 | 0.02023 | gaining |
| yazar | gamma-delta T cell | 111 | +0.00558 | 0.00468 | 932 | -0.02591 | 0.02276 | losing |
| yazar | memory B cell | 227 | -0.00822 | 0.00401 | 932 | -0.05696 | 0.01601 | losing |
| yazar | mucosal-associated invariant T cell | 25 | +0.00186 | 0.01413 | 932 | -0.14055 | 0.02543 | losing |
| yazar | naive B cell | 307 | -0.00225 | 0.00286 | 932 | -0.00334 | 0.02153 | losing |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 777 | -0.00110 | 0.00165 | 932 | -0.08920 | 0.01720 | losing |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 344 | +0.00283 | 0.00309 | 932 | -0.40476 | 0.02024 | losing |
| yazar | natural killer cell | 783 | -0.00464 | 0.00156 | 932 | +0.11766 | 0.01408 | gaining |
| yazar | plasmablast | 14 | -0.04974 | 0.06391 | 932 | -0.08214 | 0.02597 | losing |
| yazar | regulatory T cell | 54 | -0.00024 | 0.01338 | 932 | +0.06007 | 0.01212 | gaining |
| yazar | transitional stage B cell | 185 | -0.02189 | 0.00508 | 932 | -0.03614 | 0.02120 | losing |

## Verdicts under the frozen rules

- Part A Gao: **INCONCLUSIVE** (2 primary test 1 BH hit(s) but 0 predicted pass-to-fail gene changes with absolute-slope attenuation (evaluated removal types: 2)).
- Part A Yazar: **PROVES** (3 primary test 1 BH hit(s) and 96 predicted pass-to-fail gene change(s) with absolute-slope attenuation across 3 evaluated type(s)).
- Part A overall: **PROVES**.
- Part B: **DOES NOT PROVE** (association is null or opposite in both mature datasets).

Neither Part A nor Part B alone proves survivor bias as a mechanism. These tests evaluate within-type high-stress cell dependence and cross-type narrowing localization only. The Nov 17 cross-type correlation was not rerun here.

## Comparison vs original full-range Block 7

| Component | Original full-range | Mature 30-plus | Mature summary |
| --- | --- | --- | --- |
| Part A Gao | DOES NOT PROVE | INCONCLUSIVE | 2 primary test 1 BH hit(s) but 0 predicted pass-to-fail gene changes with absolute-slope attenuation (evaluated removal types: 2) |
| Part A Yazar | PROVES | PROVES | 3 primary test 1 BH hit(s) and 96 predicted pass-to-fail gene change(s) with absolute-slope attenuation across 3 evaluated type(s) |
| Part A overall | PROVES | PROVES | PROVES if either dataset PROVES |
| Part B | DOES NOT PROVE | DOES NOT PROVE | association is null or opposite in both mature datasets |

### Original test 1 BH hits and any new mature test 1 BH hits

| Dataset | Cell type | Status vs original | N orig | N mature | Orig slope | Orig BH q | Mature slope | Mature BH q | Orig predicted pass-to-fail | Mature predicted pass-to-fail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Lymphoid | NEW_AT_30_PLUS | 59 | 56 | -0.020645 | 1.369e-01 | -0.029237 | 3.179e-02 | NA | 0 |
| gao | Pericyte | NEW_AT_30_PLUS | 95 | 88 | -0.013091 | 2.211e-01 | -0.024178 | 3.777e-02 | NA | 0 |
| yazar | central memory CD4-positive, alpha-beta T cell | SURVIVED | 955 | 909 | +0.006311 | 7.484e-06 | +0.006713 | 8.269e-05 | 77 | 37 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | SURVIVED | 821 | 777 | +0.007257 | 1.217e-05 | +0.008188 | 8.269e-05 | 83 | 56 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | SURVIVED | 384 | 344 | +0.010752 | 1.675e-04 | +0.012181 | 2.504e-03 | 8 | 3 |

### Part B comparison vs original full-range Block 7

| Dataset | Orig N types | Orig slope | Orig p | Mature status | Mature N types | Mature slope | Mature p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gao | 13 | -0.03739 | 3.138e-01 | EVALUATED | 13 | -0.00537 | 9.241e-01 |
| yazar | 17 | -0.00901 | 2.224e-01 | EVALUATED | 17 | -0.01384 | 1.241e-01 |

## Caveats and power notes

- Removing captured high-scoring cells in Part A test 2 is a lower bound because dying or severely stressed cells are under-captured in single-cell and single-nucleus RNA-seq.
- Any primary test 1 null with MDE80 exceeding 0.50 outcome SD per decade is treated as underpowered (`INCONCLUSIVE`), not as a powered null.
- Neither Part A nor Part B alone establishes survivor bias as a causal mechanism.
- All randomness uses fixed master seed `20261117` and deterministic per-donor derived seeds (`SeedSequence([20261117, 7, dataset_code, donor_idx])`). Donor is the regression unit in every Part A model, and cell type is the regression unit in Part B.
