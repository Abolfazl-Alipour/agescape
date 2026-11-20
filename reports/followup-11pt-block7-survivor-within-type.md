# 11pt Follow-Up Block 7: Survivor-Bias Within-Type Tests

## Frozen protocol and gene list

This protocol and gene list were frozen before any Block 7 score, model, or result was computed. The source reference lookup did not download expression data. No gene symbol, threshold, model, family, or verdict rule may be edited after seeing results.

### Inputs and exact source files

- Gao h5ad: `/home/reza/projects/agescape-revision-data/tmp_organs/organ17_gao_heart/GSE290367_integrated_RNA_adata.h5ad`.
- Yazar h5ad: `/home/reza/projects/agescape-revision-data/tmp_organs/hcata/yazar_1e44db10.h5ad`, with `YAZAR_CACHE` and exclusion of rows labeled `peripheral blood mononuclear cell`.
- Donor tables: `results/followup-between-subject-gao-donors.csv` and `results/followup-between-subject-yazar-donors.csv`.
- Full-depth Part 2 spread inputs, identified by listing `results/`: `results/followup-fulldepth-gao-spread.csv` and `results/followup-fulldepth-yazar-spread.csv`.
- Paper-cohort fraction input: `results/reviewer-phase1-composition-global-bh.csv`, 431 tests, `effect` is CLR slope per year. It is recorded as required, but is not joined to Gao or Yazar because it contains no Gao or Yazar rows and no matching full-depth Part 2 spread slopes.
- Gao fraction slopes: `results/followup-organ-17-heart-gao-fractions.csv`, `subset=pooled_7study_primary`, `effect_clr` per decade.
- Yazar fraction slopes: `results/followup-hcata-pool-sensitivity-fractions.csv`, `study=yazar`, `arm=pool_fe`, `effect_clr` per decade.
- The existing Block 1 donor-sharded cache is reused. It contains raw integer counts after duplicate-symbol collapse and the human technical mask, restricted to post-mask UMI >= 1500 cells. Thus this block does not reread the large h5ad matrices and does not load them into memory.

### Frozen death/stress gene list and citation

Source: Liberzon et al. 2015, Cell Systems, The Molecular Signatures Database Hallmark Gene Set Collection, DOI `10.1016/j.cels.2015.12.004`, PMID `26771021`. Official MSigDB sets: `HALLMARK_APOPTOSIS`, systematic ID `M5902`, 161 genes, https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_APOPTOSIS.html; and `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`, systematic ID `M5922`, 113 genes, https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_UNFOLDED_PROTEIN_RESPONSE.html. Their union is 272 genes because `ATF3` and `DNAJC3` occur in both.

Apoptosis source genes: ADD1, AIFM3, ANKH, ANXA1, APP, ATF3, AVPR1A, BAX, BCAP31, BCL10, BCL2L1, BCL2L10, BCL2L11, BCL2L2, BGN, BID, BIK, BIRC3, BMF, BMP2, BNIP3L, BRCA1, BTG2, BTG3, CASP1, CASP2, CASP3, CASP4, CASP6, CASP7, CASP8, CASP9, CAV1, CCNA1, CCND1, CCND2, CD14, CD2, CD38, CD44, CD69, CDC25B, CDK2, CDKN1A, CDKN1B, CFLAR, CLU, CREBBP, CTH, CTNNB1, CYLD, DAP, DAP3, DCN, DDIT3, DFFA, DIABLO, DNAJA1, DNAJC3, DNM1L, DPYD, EBP, EGR3, EMP1, ENO2, ERBB2, ERBB3, EREG, ETF1, F2, F2R, FAS, FASLG, FDXR, FEZ1, GADD45A, GADD45B, GCH1, GNA15, GPX1, GPX3, GPX4, GSN, GSR, GSTM1, GUCY2D, H1-0, HGF, HMGB2, HMOX1, HSPB1, IER3, IFITM3, IFNB1, IFNGR1, IGF2R, IGFBP6, IL18, IL1A, IL1B, IL6, IRF1, ISG20, JUN, KRT18, LEF1, LGALS3, LMNA, LUM, MADD, MCL1, MGMT, MMP2, NEDD9, NEFH, PAK1, PDCD4, PDGFRB, PEA15, PLAT, PLCB2, PLPPR4, PMAIP1, PPP2R5B, PPP3R1, PPT1, PRF1, PSEN1, PSEN2, PTK2, RARA, RELA, RETSAT, RHOB, RHOT2, RNASEL, ROCK1, SAT1, SATB1, SC5D, SLC20A1, SMAD7, SOD1, SOD2, SPTAN1, SQSTM1, TAP1, TGFB2, TGFBR3, TIMP1, TIMP2, TIMP3, TNF, TNFRSF12A, TNFSF10, TOP2A, TSPO, TXNIP, VDAC2, WEE1, XIAP.

Unfolded protein response source genes: ALDH18A1, ARFGAP1, ASNS, ATF3, ATF4, ATF6, ATP6V0D1, BAG3, BANF1, CALR, CCL2, CEBPB, CEBPG, CHAC1, CKS1B, CNOT2, CNOT4, CNOT6, CXXC1, DCP1A, DCP2, DCTN1, DDIT4, DDX10, DKC1, DNAJA4, DNAJB9, DNAJC3, EDC4, EDEM1, EEF2, EIF2AK3, EIF2S1, EIF4A1, EIF4A2, EIF4A3, EIF4E, EIF4EBP1, EIF4G1, ERN1, ERO1A, EXOC2, EXOSC1, EXOSC10, EXOSC2, EXOSC4, EXOSC5, EXOSC9, FKBP14, FUS, GEMIN4, GOSR2, H2AX, HERPUD1, HSP90B1, HSPA5, HSPA9, HYOU1, IARS1, IFIT1, IGFBP1, IMP3, KDELR3, KHSRP, KIF5B, LSM1, LSM4, MTHFD2, MTREX, NABP1, NFYA, NFYB, NHP2, NOLC1, NOP14, NOP56, NPM1, PAIP1, PARN, PDIA5, PDIA6, POP4, PREB, PSAT1, RPS14, RRP9, SDAD1, SEC11A, SEC31A, SERP1, SHC1, SKIC3, SLC1A4, SLC30A5, SLC7A5, SPCS1, SPCS3, SRPRA, SRPRB, SSR1, STC2, TARS1, TATDN2, TSPYL2, TUBB2A, VEGFA, WFS1, WIPI1, XBP1, XPOT, YIF1A, YWHAZ, ZBTB17.

The shared human technical mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. `RPS14` is the one source-list gene excluded by this mask. It is explicitly recorded here and is not silently dropped. The frozen scoring list is the remaining 271 unique symbols. Symbols absent from a dataset contribute zero counts and remain in the fixed denominator.

Final scoring symbols: ADD1, AIFM3, ALDH18A1, ANKH, ANXA1, APP, ARFGAP1, ASNS, ATF3, ATF4, ATF6, ATP6V0D1, AVPR1A, BAG3, BANF1, BAX, BCAP31, BCL10, BCL2L1, BCL2L10, BCL2L11, BCL2L2, BGN, BID, BIK, BIRC3, BMF, BMP2, BNIP3L, BRCA1, BTG2, BTG3, CALR, CASP1, CASP2, CASP3, CASP4, CASP6, CASP7, CASP8, CASP9, CAV1, CCL2, CCNA1, CCND1, CCND2, CD14, CD2, CD38, CD44, CD69, CDC25B, CDK2, CDKN1A, CDKN1B, CEBPB, CEBPG, CFLAR, CHAC1, CKS1B, CLU, CNOT2, CNOT4, CNOT6, CREBBP, CTH, CTNNB1, CXXC1, CYLD, DAP, DAP3, DCN, DCP1A, DCP2, DCTN1, DDIT3, DDIT4, DDX10, DFFA, DIABLO, DKC1, DNAJA1, DNAJA4, DNAJB9, DNAJC3, DNM1L, DPYD, EBP, EDC4, EDEM1, EEF2, EGR3, EIF2AK3, EIF2S1, EIF4A1, EIF4A2, EIF4A3, EIF4E, EIF4EBP1, EIF4G1, EMP1, ENO2, ERBB2, ERBB3, EREG, ERN1, ERO1A, ETF1, EXOC2, EXOSC1, EXOSC10, EXOSC2, EXOSC4, EXOSC5, EXOSC9, F2, F2R, FAS, FASLG, FDXR, FEZ1, FKBP14, FUS, GADD45A, GADD45B, GCH1, GEMIN4, GNA15, GOSR2, GPX1, GPX3, GPX4, GSN, GSR, GSTM1, GUCY2D, H1-0, H2AX, HERPUD1, HGF, HMGB2, HMOX1, HSP90B1, HSPA5, HSPA9, HSPB1, HYOU1, IARS1, IER3, IFIT1, IFITM3, IFNB1, IFNGR1, IGF2R, IGFBP1, IGFBP6, IL18, IL1A, IL1B, IL6, IMP3, IRF1, ISG20, JUN, KDELR3, KHSRP, KIF5B, KRT18, LEF1, LGALS3, LMNA, LSM1, LSM4, LUM, MADD, MCL1, MGMT, MMP2, MTHFD2, MTREX, NABP1, NEDD9, NEFH, NFYA, NFYB, NHP2, NOLC1, NOP14, NOP56, NPM1, PAIP1, PAK1, PARN, PDCD4, PDGFRB, PDIA5, PDIA6, PEA15, PLAT, PLCB2, PLPPR4, PMAIP1, POP4, PPP2R5B, PPP3R1, PPT1, PREB, PRF1, PSAT1, PSEN1, PSEN2, PTK2, RARA, RELA, RETSAT, RHOB, RHOT2, RNASEL, ROCK1, RRP9, SAT1, SATB1, SC5D, SDAD1, SEC11A, SEC31A, SERP1, SHC1, SKIC3, SLC1A4, SLC20A1, SLC30A5, SLC7A5, SMAD7, SOD1, SOD2, SPCS1, SPCS3, SPTAN1, SQSTM1, SRPRA, SRPRB, SSR1, STC2, TAP1, TARS1, TATDN2, TGFB2, TGFBR3, TIMP1, TIMP2, TIMP3, TNF, TNFRSF12A, TNFSF10, TOP2A, TSPO, TSPYL2, TUBB2A, TXNIP, VDAC2, VEGFA, WEE1, WFS1, WIPI1, XBP1, XIAP, XPOT, YIF1A, YWHAZ, ZBTB17.

### Part A specification

- Per cell, program score is the mean within-cell proportion across the fixed 271 scoring genes: the sum of their post-mask raw counts divided by `271 * post-mask UMI`. Primary analyses use full-depth counts and no downsampling. Cells with post-mask UMI < 1500 were excluded as an eligibility filter only.
- For each dataset and cell type, donor-cell eligibility is at least 20 eligible cells. The high-score threshold is the 90th percentile, computed with numpy linear quantile interpolation, among cells from eligible donor-cell-type pairs pooled across donors. High means `score >= threshold`, including ties. The threshold and tie convention are fixed here. The same full-depth high-score set is removed in test 2.
- Test 1 outcome is each donor's fraction of high-scoring cells among that donor's eligible cells in the cell type. Primary OLS is outcome on `age_per_decade + sex_male + batch`: Gao 7-study fixed effects and Yazar 75-pool fixed effects. At least 6 donors and the repository rank and residual-df gates are required. The unadjusted sensitivity retains age and sex but removes batch. No unadjusted fallback is substituted for a failed batch-adjusted model.
- The labeled 1500 sensitivity repeats test 1 only. For each cell, aggregate source-list counts after multinomial rarefaction to 1500 UMIs have the exact marginal distribution `Binomial(original_list_count, 1500 / original_UMI)`. Independent fixed-seed draws across cells therefore simulate the exact aggregate-list score without changing any primary data. The 90th-percentile rule is recomputed on rarefied scores within each cell type. Seeds derive from `SeedSequence([20261117, 7, dataset_code, donor_idx])`, with cells sorted by original h5ad row index.
- Test 2 is run only for cell types with primary test 1 BH q < 0.05. Donor pseudobulk gene models mirror the follow-up gene tests: log2CPM, age per decade, sex, log library size, and the same study or pool fixed effects; genes pass the original expression filter of CPM >= 1 in at least half of donors, with a minimum of 3 donors. Both arms require at least 20 cells per donor after their respective cell filters and at least 6 common donors. Unremoved and removed fits use the same common donor set. The gene family is the original `q_bh_within_group` family: BH separately within each dataset, cell type, and removal arm across that arm's expressed genes. Pass means q < 0.05. The descriptive predicted change is an unremoved pass becoming a removed fail while the absolute age slope becomes smaller. Slope change is removed minus unremoved. All other pass/fail changes are also counted.
- Removing captured high-scoring cells is a lower bound. Dying cells are under-captured, so removal cannot remove all dying or stressed cells.
- Test 1 BH families are stated before correction as each `(dataset, test)` slice across evaluable non-pooled cell types, separately for primary full depth, unadjusted full depth, and rarefied 1500. Test 2 uses the gene-level family stated above. No BH correction is applied across slope changes.
- MDE80 is `(t_0.975 + t_0.80) * SE`, using residual df, on the raw outcome scale and standardized by the donor-outcome SD. A small MDE80 for test 1 is predefined as at most 0.50 outcome SD per decade, matching the repository convention.

### Part B specification

- Spread slopes are read, not recomputed, from the two exact full-depth CSVs above. Required rows are `metric=shannon`, `depth_mode=full_depth`, `analysis_role=primary_full_depth`, `model_spec=primary`, non-pooled, and non-stratified-global. These are the batch-adjusted and depth-adjusted per-cell-type primary rows needed for a cell-type join. A `stratified_global_primary` row is one cross-type summary and cannot provide per-cell-type slopes.
- Gao and Yazar are joined separately on exact cell-type keys to their batch-adjusted CLR fraction slopes above. Fraction slope negative means losing share. Spread slope is Shannon change per decade, so a decrease is negative.
- Prediction frozen before running: survivor bias predicts a positive association between fraction slope and spread slope, with more negative spread slopes concentrated in fraction-losing types.
- For each dataset, fit weighted least squares of spread slope on fraction slope across cell types with inverse variance weights `1 / SE_spread^2`. Report the association slope, 95 percent CI, two-sided p, and MDE80. Also report mean spread slope and cell-type counts for fraction-losing versus fraction-gaining types. Zero fraction slope is included in the gaining group.
- There is one test per dataset and no multiplicity correction across the two datasets. The datasets are never pooled.

### Verdict rules

- Part A, per dataset: PROVES only if test 1 passes BH correction and test 2 produces at least one predicted pass-to-fail change with absolute-slope attenuation. DOES NOT PROVE if no test 1 result passes BH and every evaluable primary null has MDE80 at most 0.50 SD per decade. INCONCLUSIVE if a null has larger MDE80, or test 1 is significant but removal changes nothing in the predicted direction, or the required removal comparison is not evaluable. The overall Part A summary is PROVES if either dataset meets PROVES, DOES NOT PROVE only if both datasets meet DOES NOT PROVE, and otherwise INCONCLUSIVE.
- Part B: PROVES only if the cross-type association is positive and significant in at least one dataset and the other dataset is not significantly opposite. DOES NOT PROVE if the association is null or opposite in both datasets. A significant positive result in one dataset with a significant opposite result in the other is INCONCLUSIVE.
- Neither Part A nor Part B alone proves survivor bias as a mechanism. The Nov 17 cross-type correlation is not rerun.

### Outputs and reproducibility

- `results/followup-11pt-block7-death-stress.csv`
- `results/followup-11pt-block7-removal-sensitivity.csv`
- `results/followup-11pt-block7-narrowing-localization.csv`
- `reports/followup-11pt-block7-survivor-within-type.md`
- `scripts/run_followup_11pt_block7_survivor_within_type.py`

All randomness uses fixed seed `20261117` and deterministic derived seeds. Donor is the unit in every regression. No permutation or bootstrap is used. No split-half technical null is required by this design. No primary analysis uses downsampling. Exactly one deep review round follows completion of the outputs.

## Results

All three result CSVs were generated from the frozen protocol above. Primary analyses use full-depth counts. Donor is the unit in every regression.

### Part A test 1: high death/stress-score cell fractions

#### Gao

Evaluated 13 of 13 cell types. BH family size 13; 0 passed q < 0.05.

| Cell type | Status | N donors | Slope per decade | 95 percent CI | p | BH q | MDE80 SD |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Adipocyte | EVALUATED | 37 | +0.013710 | -0.010417 to 0.037836 | 2.550e-01 | 5.823e-01 | 0.302 |
| Cardiomyocyte | EVALUATED | 108 | -0.006765 | -0.020016 to 0.006485 | 3.135e-01 | 5.823e-01 | 0.180 |
| Endocardial | EVALUATED | 72 | -0.004586 | -0.018737 to 0.009565 | 5.197e-01 | 6.093e-01 | 0.243 |
| Endothelial | EVALUATED | 86 | -0.005588 | -0.019805 to 0.008629 | 4.363e-01 | 6.093e-01 | 0.211 |
| Epicardial | EVALUATED | 19 | -0.028168 | -0.066396 to 0.010061 | 1.372e-01 | 4.457e-01 | 0.272 |
| Fibroblast | EVALUATED | 105 | -0.004689 | -0.020700 to 0.011323 | 5.625e-01 | 6.093e-01 | 0.181 |
| LEC | EVALUATED | 36 | -0.005896 | -0.035632 to 0.023841 | 6.887e-01 | 6.887e-01 | 0.356 |
| Lymphoid | EVALUATED | 59 | -0.020645 | -0.036255 to -0.005036 | 1.053e-02 | 1.369e-01 | 0.227 |
| Mast | EVALUATED | 17 | -0.013870 | -0.053829 to 0.026089 | 4.609e-01 | 6.093e-01 | 0.303 |
| Myeloid | EVALUATED | 106 | -0.005448 | -0.015318 to 0.004423 | 2.761e-01 | 5.823e-01 | 0.206 |
| Neuronal | EVALUATED | 60 | -0.019357 | -0.038726 to 0.000012 | 5.013e-02 | 2.211e-01 | 0.255 |
| Pericyte | EVALUATED | 95 | -0.013091 | -0.026242 to 0.000059 | 5.102e-02 | 2.211e-01 | 0.168 |
| vSMC | EVALUATED | 86 | -0.003799 | -0.016462 to 0.008863 | 5.521e-01 | 6.093e-01 | 0.185 |

gao test1_rarefied_1500: 13 evaluated, 0 BH hits.
gao test1_unadjusted_full: 13 evaluated, 0 BH hits.

#### Yazar

Evaluated 16 of 28 cell types. BH family size 16; 3 passed q < 0.05.

| Cell type | Status | N donors | Slope per decade | 95 percent CI | p | BH q | MDE80 SD |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | EVALUATED | 252 | -0.006569 | -0.015724 to 0.002587 | 1.587e-01 | 3.596e-01 | 0.139 |
| CD14-positive monocyte | EVALUATED | 423 | -0.000631 | -0.005033 to 0.003771 | 7.782e-01 | 8.301e-01 | 0.085 |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_NOT_EVALUABLE_BATCH_SINGLE_LEVEL | 6 | NA | NA to NA | NA | NA | NA |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| CD4-positive, alpha-beta cytotoxic T cell | EVALUATED | 75 | +0.006987 | -0.009226 to 0.023200 | 3.907e-01 | 6.251e-01 | 0.241 |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| central memory CD4-positive, alpha-beta T cell | EVALUATED | 955 | +0.006311 | 0.003872 to 0.008751 | 4.677e-07 | 7.484e-06 | 0.042 |
| central memory CD8-positive, alpha-beta T cell | EVALUATED | 91 | +0.010315 | -0.012515 to 0.033146 | 3.699e-01 | 6.251e-01 | 0.281 |
| conventional dendritic cell | UNTESTABLE_DONOR_FLOOR | 5 | NA | NA to NA | NA | NA | NA |
| dendritic cell | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| double negative thymocyte | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| effector memory CD4-positive, alpha-beta T cell | EVALUATED | 261 | +0.002385 | -0.005186 to 0.009955 | 5.353e-01 | 6.605e-01 | 0.103 |
| effector memory CD8-positive, alpha-beta T cell | EVALUATED | 746 | +0.002580 | -0.001140 to 0.006301 | 1.737e-01 | 3.596e-01 | 0.056 |
| erythrocyte | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| gamma-delta T cell | EVALUATED | 120 | +0.000299 | -0.011410 to 0.012008 | 9.596e-01 | 9.596e-01 | 0.150 |
| hematopoietic precursor cell | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| innate lymphoid cell | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| memory B cell | EVALUATED | 244 | +0.005167 | -0.002403 to 0.012738 | 1.798e-01 | 3.596e-01 | 0.101 |
| mucosal-associated invariant T cell | EVALUATED | 30 | +0.010348 | -0.031479 to 0.052175 | 6.084e-01 | 6.953e-01 | 0.480 |
| naive B cell | EVALUATED | 319 | +0.008657 | 0.001824 to 0.015489 | 1.321e-02 | 5.283e-02 | 0.089 |
| naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 821 | +0.007257 | 0.004318 to 0.010195 | 1.521e-06 | 1.217e-05 | 0.052 |
| naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 384 | +0.010752 | 0.005743 to 0.015761 | 3.141e-05 | 1.675e-04 | 0.075 |
| natural killer cell | EVALUATED | 815 | +0.000932 | -0.001874 to 0.003737 | 5.146e-01 | 6.605e-01 | 0.051 |
| plasmablast | UNTESTABLE_NOT_EVALUABLE_WITHIN_BATCH_DF | 14 | NA | NA to NA | NA | NA | NA |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_FLOOR | 3 | NA | NA to NA | NA | NA | NA |
| platelet | UNTESTABLE_DONOR_FLOOR | 0 | NA | NA to NA | NA | NA | NA |
| regulatory T cell | EVALUATED | 54 | +0.006851 | -0.015423 to 0.029126 | 5.367e-01 | 6.605e-01 | 0.484 |
| transitional stage B cell | EVALUATED | 201 | +0.005931 | -0.002762 to 0.014623 | 1.797e-01 | 3.596e-01 | 0.104 |

yazar test1_rarefied_1500: 16 evaluated, 3 BH hits.
yazar test1_unadjusted_full: 17 evaluated, 2 BH hits.

### Part A test 2: removal sensitivity

| Dataset | Cell type | Status | N donors | Common genes | Unremoved pass | Removed pass | Predicted pass to fail | Other pass to fail | Fail to pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| yazar | central memory CD4-positive, alpha-beta T cell | EVALUATED | 948 | 10777 | 310 | 232 | 77 | 17 | 16 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 791 | 10072 | 319 | 240 | 83 | 12 | 16 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 349 | 8552 | 21 | 12 | 8 | 1 | 0 |

Pass/fail change counts use genes evaluated in both arms. Slope change is removed minus unremoved. Removal is a lower bound because dying cells are under-captured.

### Part B: narrowing localization

| Dataset | N types | Association slope | 95 percent CI | p | MDE80 | Losing N | Gaining N | Mean spread losing | Mean spread gaining |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | 13 | -0.03739 | -0.11536 to 0.04057 | 3.138e-01 | 0.10898 | 7 | 6 | -0.01430 | -0.01307 |
| yazar | 17 | -0.00901 | -0.02408 to 0.00607 | 2.224e-01 | 0.02121 | 10 | 7 | -0.00687 | -0.00194 |

The association slope is the change in Shannon spread slope per unit CLR fraction slope per decade. Fraction slope negative means losing share. The prediction is a positive association. There is one test per dataset and no multiplicity correction across datasets.

## Verdicts

- Part A Gao: DOES NOT PROVE.
- Part A Yazar: PROVES.
- Part A overall: PROVES.
- Part B: DOES NOT PROVE.

For Part A, PROVES requires both a BH-passing test 1 result and at least one predicted pass-to-fail gene change with absolute-slope attenuation. A no-hit DOES NOT PROVE verdict requires every evaluated primary null to have MDE80 at most 0.50 outcome SD per decade. Larger MDE80 or a required but unevaluable removal comparison is INCONCLUSIVE.

Neither Part A nor Part B alone proves survivor bias as a mechanism. These tests address within-type stress-cell dependence and cross-type localization only. The Nov 17 cross-type correlation was not rerun.

## Reproducibility and input audit

The paper composition file contains 431 tests as required. It is not part of either cross-type model because it has no Gao or Yazar cohort rows. Gao joins 13 cell types and Yazar joins 17 cell types. Exact spread and fraction filenames, model rows, mask, gene source, fixed seeds, families, and eligibility rules are recorded in the frozen protocol.
