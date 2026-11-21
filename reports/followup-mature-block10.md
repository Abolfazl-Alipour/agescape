# Block 10 mature rerun: Published-signature benchmark at 30-plus

## Pre-analysis freeze

This registry, the exact gene lists, the analysis specification, the coverage rule, and all multiplicity families were frozen before any composite slope or matched-null model was fit. Four published human signatures were selected before running. No candidate required substitution because all four lists and directions were retrievable. No P1 or P0 genes were used as a benchmark signature.

## Mature sensitivity: what changed and what did not

Scoring is restricted to donors aged 30-plus in both primary cohorts (Gao heart, Yazar blood) and in both arms. Matched nulls are recomputed on the 30-plus donors with 2000 draws and the same fixed seed 20261105. Every other frozen rule is unchanged: the four published gene lists and directions, the 50 percent coverage rule, the eligibility floors (at least 20 cells per donor per unit and at least 6 donors per unit), the Block 5 composite and empirical tests, and BH correction within each signature, cohort, and arm across the frozen units, separately for composite and empirical tests.

Limitation, stated explicitly: the four external published signatures were NOT rederived from mature donors. Their gene lists and directions come from the original publications, whose derivations used younger donors. This rerun is a mature-TESTED sensitivity only: it asks whether the published signatures still pass when scored on donors aged 30-plus. It is not a mature-derived benchmark. Full compliance with a mature-to-old derivation would require rederiving each signature from its source study restricted to mature donors, which is outside this sensitivity.

Unique donors scored: Gao heart 109 of 117 full-range donors at 30-plus; Yazar blood 932 of 981 full-range donors at 30-plus. Donor restriction audit: `results/followup-mature-block10-donor-counts.csv` records, for every cohort, arm, and unit, the donor count at full range and at 30-plus. The unit is the donor, never the cell.

### Frozen registry and citations

### Harries2011_6gene

Harries LW, Hernandez D, Henley W, Wood AR, Holly AC, Bradley-Smith RM, Yaghootkar H, Dutta A, Murray A, Frayling TM, Guralnik JM, Bandinelli S, Singleton A, Ferrucci L, Melzer D. 2011. Human aging is characterized by focused changes in gene expression and deregulation of alternative splicing. Aging Cell 10(5):868-878. DOI: 10.1111/j.1474-9726.2011.00726.x.

Source: Main article section 'Six of the most deregulated 25 probes provide good discrimination between young and old', with directions from Table 1a and Table 1b. Exact model genes: LRRN3, CD248, CCR6, GRAP, VAMP5, CD27.

Retrieval URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC3173580/

Published content: 6 genes, 1 up and 5 down. Scoring uses 6 unique human symbols, 1 up and 5 down.

**Exact up list used (1):** VAMP5

**Exact down list used (5):** LRRN3, CD248, CCR6, GRAP, CD27

### deMagalhaes2009_73gene

de Magalhaes JP, Curado J, Church GM. 2009. Meta-analysis of age-related gene expression profiles identifies common signatures of aging. Bioinformatics 25(7):875-881. DOI: 10.1093/bioinformatics/btp073.

Source: MSigDB C2 CGP gene sets DEMAGALHAES_AGING_UP and DEMAGALHAES_AGING_DN, contributed by the paper author and sourced to the paper's published meta-signature. The 56 up and 17 down source identifiers map to 54 and 16 current unique human symbols in MSigDB.

Retrieval URL: https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=DEMAGALHAES_AGING_UP&fileType=gmt ; https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=DEMAGALHAES_AGING_DN&fileType=gmt

Published content: 73 source identifiers, 56 up and 17 down. Scoring uses 70 unique human symbols, 54 up and 16 down.

**Exact up list used (54):** ADIPOR2, ANXA3, ANXA5, APOD, B2M, C1QA, C1QB, C1QC, C3, C4A, CLIC4, CLU, CTSS, DCLK1, DERL1, EFCAB14, EFEMP1, FCGR2A, FCGR2B, GBP2, GFAP, GNS, GPNMB, GSTA1, H1-2, HBA1, HCST, HLA-G, IL33, JCHAIN, LAPTM5, LGALS3, LITAF, LYZ, MGST1, MPEG1, MSN, MT1F, NDRG1, NPC2, PCSK6, PSMD11, PTGES3, RASA3, RNF213, S100A4, S100A6, SERPING1, SGK1, SPP1, TMED10, TXNIP, VAT1, VWF

**Exact down list used (16):** ACSS2, ATP5MC3, CA4, CALB1, COL1A1, COL3A1, COL4A5, CX3CL1, DIABLO, FABP3, GHITM, NDUFB11, NREP, TFRC, UQCRFS1, UQCRQ

### Sood2015_150probe

Sood S, Gallagher IJ, Lunnon K, Rullman E, Keohane A, Crossland H, Phillips BE, Cederholm T, Jensen T, van Loon LJC, Lannfelt L, Kraus WE, Atherton PJ, Howard R, Gustafsson T, Hodges A, Timmons JA. 2015. A novel multi-tissue RNA diagnostic of healthy ageing relates to cognitive health status. Genome Biology 16:185. DOI: 10.1186/s13059-015-0750-x.

Source: Additional file 1, worksheet 'Top 150', columns 'HU133 Probeset_ID', 'Gene Symbol', and 'Ratio of Y:0 muscle'. The 150 probe sets comprise 148 unique symbols because RUNX1 and FAM20C each occur twice.

Retrieval URL: https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0750-x/MediaObjects/13059_2015_750_MOESM1_ESM.xlsx

Published content: 150 probe sets, 12 up and 138 down. Scoring uses 148 unique human symbols, 12 up and 136 down.

**Exact up list used (12):** MBNL1, CFLAR, RBMS3, DHX29, PCDH9, PNISR, SKAP2, RP11-271C24.3, ALDH6A1, TGFBR3, UNC13C, CD36

**Exact down list used (136):** HIST1H3E, RAB3A, MAPK8IP2, ZNF205, SMTNL2, SHISA4, KLK15, TFRC, SARDH, DCBLD2, LOC100128008, ARMC5, BAIAP2L2, CNTROB, MAPRE3, 217079_at, FUT7, CNPY4, RBM17, PCYT2, HMGB2, CARM1, HCN4, 215488_at, AK1, TNPO2, TNK2, SIN3A, C7orf50, 243906_at, ARF1, TNFRSF10C, SIGLEC8, AGER, SIRT5, C6orf136, FAM57B, U2AF2, RUNX1, 213690_s_at, CDK13, CDC42EP4, ADAP1, LOC100652912, ZNF628, SLC25A5-AS1, GPC2, ADAM11, RIF1, 244182_at, NPEPL1, COL1A1, NRG1, MRC2, SERGEF, ATG4D, FAM223B, RNF207, PTCRA, LMNA, CALR, RPP25, RASAL3, ZNF503, TRIM62, PWWP2B, TFR2, NRXN2, SCN2B, MRM1, AMN, CHRAC1, EMID2, ENOX2, 235671_at, GALNT1, PLCH2, ARL3, ABCB9, MVK, ATP2A1, SLC35F3, TUSC5, 240116_at, FBXW7, SPN, ZNF574, POF1B, CCDC149, PLEKHB1, GTSE1, IL34, 239060_at, 240241_at, EVX1, PRKAR2A, AGRN, ZCCHC3, 234400_at, ARHGEF1, ICAM5, FLJ11710, BHLHE41, FAM20C, KIF20B, SOX30P1, EIF3H, MMP11, ZDHHC4, EPB41L4B, FAM27E3, 231402_at, ADAMTS16, JAG1, NENF, SCNN1A, SMOX, RNASET2, CRYGEP, PHF19, ANKRD13B, SEMA7A, SURF2, ATXN7L3B, LOC79999, TMUB1, BLOC1S3, SEZ6L2, JMJD8, CAMTA1, CITED2, IL12RB1, DNAJC22, PRX, 241563_at, OTUB2

### Peters2015_1497gene

Peters MJ, Joehanes R, Pilling LC, Schurmann C, Conneely KN, Powell J, Reinmaa E, et al. 2015. The transcriptional landscape of age in human peripheral blood. Nature Communications 6:8570. DOI: 10.1038/ncomms9570.

Source: Supplementary Data 1, worksheet 'Data.1'. Select rows with sign-DISC = y and sign-REPL = y, then use NEW-Gene-ID and Direction. This yields exactly 1497 unique genes, 600 positive and 897 negative, matching the paper.

Retrieval URL: https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9570/MediaObjects/41467_2015_BFncomms9570_MOESM436_ESM.xlsx

Published content: 1497 genes, 600 up and 897 down. Scoring uses 1497 unique human symbols, 600 up and 897 down.

**Exact up list used (600):** GZMH, GPR56, SYT11, FGFBP2, TGFBR3, FCRL6, CX3CR1, RHOC, PRSS23, NKG7, PPP2R2B, GZMB, TTC38, PRR5L, ADRB2, OSBPL5, C1orf21, LGALS1, NUAK1, GTSCR1, PDGFD, ANXA1, APOBEC3H, MCOLN2, KYNU, F2R, PRF1, CACNA2D2, ATP2B4, CD84, GZMA, PDGFRB, IGFBP3, SUSD1, KLF6, KLRD1, ZEB2, VAMP5, IGFBP7, F5, B3GAT1, VCAN, CRIP1, HAVCR2, GNLY, GK5, ERBB2, PYHIN1, NME8, FAM13A, AGPAT4, STOM, AHNAK, S100A10, CLDND2, PAPSS2, KLF11, GPR141, ZBP1, CHN1, SLFN11, CD244, ILK, CST7, PTGDS, KLRF1, PDK4, NMUR1, PTGDR, LGR6, GBP5, BHLHE40, APOBEC3G, LPP, C1orf162, P2RX7, SLC25A20, GAB3, EOMES, KLF10, MTSS1, PTPN7, TSHZ2, SNTB2, CTNNA1, IL10RA, SORT1, PTGER2, APH1B, MARCH1, ITGB1, CD70, ZNF365, MT1A, PDZD4, CKAP4, FASLG, BATF, GPR137B, IRAK3, ASCL2, LPCAT1, CYP4F22, MAFB, FAM214B, LINC00152, CMKLR1, SH2D1B, ANXA2, TMX4, MT1E, PLBD1, KIF1B, GPR114, NCAM1, SLC31A1, MAP3K5, ARG1, TTLL3, CPT1A, LMNA, FCGR1A, IL1R2, PARP3, ASGR2, CAPN2, STAT1, CPEB3, RRAS, LRRFIP2, ROPN1L, GFI1, PLEK, CCL5, VAMP3, GLIPR2, FBXO6, IDS, MEGF9, TSPAN2, YPEL1, MBOAT1, RCAN2, FRMD3, B4GALT5, SLC26A8, WARS, SLC22A15, ABI3, ANKRD22, CD58, AGPAT4-IT1, RNASE2, TPST1, TIGIT, RAP1GAP2, CD163, C3AR1, SRPK1, LHFPL2, EIF4G3, CDKN1C, SLAMF6, MS4A6A, ACVR2A, WIPI1, NBEAL2, RALGDS, UPP1, APMAP, ALOX5AP, SERPINB8, ANTXR2, GDPD5, LYPD2, SLC31A2, VSIG4, GABARAPL1, CD63, MGC24103, RASSF1, SLCO4C1, SERPINB6, ENTPD1, PRDM1, DSC2, E2F2, EFHD2, ARHGEF40, TP53I3, ATL3, RAP2A, CTSW, ANG, AKR1C3, SLC27A3, ATXN1, UBR2, PALLD, NPC2, TNF, SLFN13, OGDH, STX11, RAB11FIP5, CTSL1, MGST3, C19orf59, ANKRD13A, GALC, NDST1, TMTC1, CDC42EP4, MS4A4A, PARP8, ADAP2, TRAFD1, NPC1, MSC, HK3, CRISPLD2, LILRB1, TIPARP, LINGO2, GBP1, SFTPD, IFNG, SYNE2, FAM46A, FAR2, ACTA2, SMARCD3, DKFZP586I1420, HMOX1, PCNX, NRG1, CYB561, CCDC170, SLC7A7, ZNF683, TNFAIP6, DMXL2, CD86, ATG7, KIF13B, VIM, CAPNS1, MR1, PRDM8, PLA2G16, GNS, VNN3, MMP9, RBPJ, SLCO3A1, PROS1, ARHGAP26, METRNL, MAPK14, PLOD1, LGALSL, GOLM1, GBA, ITGA2B, HERPUD2, CTSG, DUSP18, HMGB2, MPO, NOD2, GNPTAB, TCN2, PHC2, USP3, GSTO1, FGR, LPAR6, FOXD2, TNIP3, EHD4, PRSS30P, SLC4A4, VCL, PPP2R5C, IL12RB1, SMAD7, TSPAN17, CSRP1, ANXA4, MVP, CEBPA, MFSD10, TTYH3, PFKFB2, DOK2, CLSTN3, SLAMF7, TCN1, ACAA2, DUSP3, JAKMIP1, RAB31, CD151, SLC1A7, CAST, SGPL1, RAB37, ID2, FGL2, GLT25D2, PADI4, MYOF, CFLAR, DUSP5, S1PR5, CD97, F13A1, PAQR4, NACC2, SRBD1, SERPINB1, SESTD1, DIAPH1, RAB27A, MXRA7, CREB3L2, PLAUR, ITGAL, LILRA5, TBXAS1, ARHGAP19, CYFIP1, DIRC2, ABCC3, APOBEC3C, ITGAM, GYG1, DPYSL2, EDARADD, ADAM17, CXCL16, PLEKHF1, UBE2F, C12orf75, ADCY9, FAM126B, VNN1, CPM, ADAMTSL4, ADHFE1, WDFY2, AGFG1, KLF4, LOXL3, PLXDC2, DEFA4, CYP1B1, PCBD1, PGD, YIPF1, CD300C, NCEH1, S100A8, LINC00260, TAGAP, MAP3K3, EMP1, TSPO, TCTA, FEZ1, LINC00324, WBP1L, GALNT10, MAF, ANXA5, VAV3, LAIR2, SOX13, SLC9A9, STX6, CLEC4E, ST6GALNAC6, PARP9, TIMP1, TTC16, SH3BP2, CLIC1, S100A12, GNG2, FAM53B, QPCT, SH3BGRL3, KIAA0319L, BNC2, NLRP12, PIK3CG, C3orf62, THBS1, DGKG, SLC22A4, DENND1B , APLP2, FLVCR2, VOPP1, DYSF, RASSF4, CDC14B, ZNF385A, MARCO, DLGAP4, AKIRIN2, TRIM69, AKAP13, CSF1R, PHTF1, IL18BP, BTN3A1, HCG27, FBXL5, TP53I11, ABR, TNFSF13B, GCNT1, DRAXIN, DAAM2, BCL6, CTSD, AGPAT2, IL15, ABCA6, GAS7, RAB24, MAPK1, P2RX1, LPGAT1, MYL6, DYNLT1, INPPL1, MAPKAPK3, ATP1B3, NCR1, GRB10, HP, GDE1, C1orf216, TRIM21, TCIRG1, SORCS3, LRRC33, SIAH2, FAM134A, ERGIC1, ACADVL, MFN2, PRRG4, FAM160A2, IRF9, FUCA2, DHRS9, TAGLN, PNPLA6, IFI30, GPX1, REPS2, CFH, RAB7L1, PTP4A2, CCR2, SPON2, BLOC1S1, DEDD, PDZK1IP1, PLIN2, STAT3, HNMT, CPNE8, CPQ, IRF1, GLA, TADA3, BPI, PSMB9, DENND1A, OASL, CTSC, EPG5, MYD88, ACOT9, AGPAT3, CDK14, MEOX1, CALCOCO1, PAM, GPR68, SH2D2A, TRIM38, ADAR, AGAP1, ABHD16A, WDR1, IQCG, PTK2B, SPPL2A, GNGT2, IDH1, NEDD9, UBAC1, TRIB1, SLC39A1, TLN1, ACTN4, GPR97, VPS33A, BMX, EPB41L4A, PLOD3, SLC20A1, IFIT3, TNFRSF1B, ENPP4, ARRB1, SSH1, NLRP3, FAM129B, SETBP1, TOX, CALML4, PLCB2, ELOVL6, CYBB, IL10RB, C4orf34, ABHD3, STAT2, AOAH, ICAM4, NAGA, SAMD9L, RARRES3, STARD3NL, INPP5D, SLITRK4, KIAA0913, ARHGAP24, ARF3, CDKN1A, ATF3, FCGR2B, SSBP4, C11orf75, TBCB, RIPK3, MAD2L2, FAM70A, MGLL, ITPRIP, ELF4, RUNX3, OAS1, EIF2C4, QSOX1, CYSTM1, MGAM, ETV6, C20orf197, EBP, SH3BP5L, HLA-B, ABHD5, ORM1, PFKFB4, GUSB, FHOD1, ETFDH, PARP14, GP6, ARSG, CMTM3, TK2, SLC26A6, RPS6KA1, DUSP10, TBC1D2, FANCA, TMBIM1, RIPK2, HLA-E, ACSL1, SNX11, EMP3, EIF4EBP2, ENDOD1, USF1, SCIMP, USP4, PTTG1, TAP1, ANKS1A, HSPA6, NLRC5

**Exact down list used (897):** CD248, LRRN3, NELL2, LEF1, CCR7, ABLIM1, MYC, CD27, FAM102A, SERPINE2, SLC16A10, FCGBP, BACH2, PDE9A, NOG, FLNB, NT5E, ITM2C, ATF7IP2, CR2, FAIM3, PHGDH, LDHB, SIRPG, PDE7A, NOSIP, PAICS, BZW2, OXNAD1, SCML1, RPL22, LDLRAP1, LTB, FAM134B, LBH, SUSD3, PIK3IP1, MFGE8, AGMAT, SLC7A6, IMPDH2, GRAP, NUCB2, FBL, SATB1, IL7R, ZSCAN18, KLRB1, HSP90AB1, VPREB3, TCF7, IGJ, BCL11B, NMT2, TNFRSF25, PTPRK, BTBD11, EDAR, ELOVL4, KLHL3, CD8B, NBEA, ATM, ID3, PCED1B, COBLL1, MZB1, GPR18, CEP41, ETS1, FAM117B, EPHX2, CAMK4, CCR6, TXK, LEPROTL1, DKC1, SPRY1, FOXO1, MMP28, LGMN, DENND2D, CLUAP1, OCIAD2, CD79A, RCAN3, TSPAN13, HOOK1, KIAA1147, THNSL1, ATIC, NAA16, CD28, RORC, RAD54B, MAN1C1, STRBP, TARBP1, BANK1, TRA@, FAM159A, NPM3, ZNF154, DGKA, KAT2A, POU2AF1, EPHA1, TPM2, C1QBP, LTBP3, PLCG1, ITGA6, TMEM204, CD96, RRS1, HSPH1, SESN1, TBC1D4, ANXA2R, CD7, STMN1, APEX1, ICOS, SAE1, FLT3LG, PMEPA1, RAPGEF6, ZNHIT6, PASK, TRIB2, DDX47, CDCA7L, GNG7, GEMIN4, TRAF5, CD3E, CXCR5, FBLN2, RNF144A, CD22, AK5, EIF2S3, HSPD1, DNAJA3, TXNDC5, BCAS4, IARS, C6orf48, RSL1D1, VEGFB, IMP3, SEPT6, POLR1C, PLAG1, AXIN2, TTC3, ZNF135, FAM129C, SPOCK2, CNN3, ZNF439, MTA3, C2orf40, EPB41L4A-AS1, EFNA1, ZNF256, NIPAL3, EXOSC5, PAQR8, ENOSF1, EIF3E, IGLL1, TEX10, MAGED1, UBASH3A, SPTBN1, HSF2, ARHGEF18, TRAP1, EIF3H, ZNF550, BCL11A, HLA-DOA, HSP90B1, SCD, POLR1E, SNRPN, ADSL, LARS, SYPL1, SLC38A1, CIRH1A, GSPT2, E2F5, SERINC5, EIF3L, CHMP7, DDX10, LAGE3, WDR75, CCDC109B, PFAS, FAM213A, TIGD7, PKIA, ARHGAP32, EEF2, DPH5, OBSCN, DOCK10, NAT10, ZNF544, TC2N, BTLA, TLE2, CCT2, C10orf2, EVL, GNL3, RPS6, IRF8, CD79B, SRPRB, GPA33, NOL11, TESPA1, PUS1, SET, MKI67IP, ITK, NRIP1, PDE3B, GRPEL2, C9orf123, SEH1L, TNFRSF17, ABCE1, MORC2, TRAF3IP2, ZNF671, DBP, ALKBH2, ATP6V0E2, MAPKAPK5-AS1, RPS4X, PPAT, KIAA0125, TNFRSF13B, LY9, WDR4, N6AMT1, CDK6, MAL, UBA52, PRMT3, UXT, RPL18, SPIB, BNIP3, BEX2, BLK, TTC27, GIMAP5, RPL13, WDR12, EIF3M, KIF5C, CCT8, DSC1, TMEM97, KLHL6, PACSIN1, KIAA0020, DPP4, MTUS1, DYRK2, RPS23, GPATCH4, PNPO, DPH2, URI1, LRPPRC, TRMT13, PPRC1, PDCD2L, BCKDHB, AFF3, PARP1, RPL8, CD19, EBF1, NOP2, CHI3L2, EIF2D, MARS2, NME1, CBLL1, CPA5, HSD17B8, C12orf57, PPA1, RPS13, RPL5, RBM17, SH3YL1, RASGRP2, SRM, FAM50B, CCDC50, TCEA3, CYP2J2, FAM153B, AKTIP, ALG8, NREP, NCL, RPUSD2, HLA-DOB, CERS6, GTPBP4, NPAT, RPLP0, MRPL45, HADH, HMGN1, CSTF3, RALA, ZNF529, ZNF22, GAR1, ADPRM, DDX55, BMS1, TMEM66, METTL1, ZNF30, LTV1, GAMT, NUP107, CCT3, CDK4, TMEM99, FCRLA, TCEA2, PHF14, AARS, NARS2, ZNF302, DNASE1L3, EIF4B, TOP1MT, C6orf108, RBM26, SRSF6, FOXJ3, H2AFY2, WRN, PUS7, C12orf23, SPON1, RBL2, ITFG2, KIAA0355, GLTSCR2, CHCHD10, LSM2, EBAG9, ZCCHC11, ADK, TSEN2, GPX7, CCT4, NOA1, NOP58, SREBF1, EIF4A2, EIF2A, MFNG, TMEM156, APBA2, SERBP1, ZBTB25, MAML2, AHCY, RGCC, RRP15, RPSAP58, NUP160, TSPAN3, DNMT1, OSBPL10, C16orf58, SMARCA1, SLC38A11, RPS2, KDM1A, RPL13P5, SHMT2, SERPINF1, HNRPDL, NOC3L, WDR59, C10orf35, CMTM8, NOL8, TIMM21, C5orf51, CYBASC3, ZDHHC9, TCF3, SERTAD2, CD3D, GIMAP7, NOLC1, C22orf29, EIF3G, MRPS27, DIMT1, TYSND1, PITPNC1, TNFRSF21, ZNF101, ANAPC1, ZNF480, BRIX1, RCL1, SLAIN1, RABGGTB, TMEM8B, WDR67, PRKCA, ESYT1, ZNF329, ZZZ3, ARMCX2, TRMT112, DDX18, YTHDC1, SDCCAG3, LGALS3BP, VPS51, TMEM123, GAL3ST4, ANKRD46, PDIA5, HDAC2, ST6GAL1, SH3PXD2A, ZCCHC7, KCNG1, TRIM28, RPL19, GTPBP3, HVCN1, MAP3K4, RPP40, PDCD4, LTA, EPRS, PAG1, TTC39B, MRI1, C11orf1, RRP1B, ZNF146, RPS16, CIAPIN1, HNRNPD, TOP2B, NHP2L1, DCXR, MIR600HG, PSIP1, PIK3C2B, C3orf37, MIF, LSM5, TCF4, ELP2, PAFAH1B3, NONO, GOLGA3, LAT, CCNB1IP1, USP11, NSMCE1, CRLF3, RPF2, UBIAD1, OCIAD1, ZNF548, CCDC41, AKR1B1, LCMT2, MTHFD1L, RPL4, MRPL3, P2RX5, RASGRP3, PRPSAP2, HAPLN3, HNRNPM, EIF3D, TRAT1, CMSS1, LUC7L, DDX21, DNMT3A, ZFAND1, LRRC8C, PLAC8, ABCB1, MRPL24, TIMD4, PVRIG, NXT1, C7orf23, THEM6, C1orf109, PLEKHB1, MRPS26, EEF1B2, CXCR7, LCK, SEPHS1, CLEC2D, TPP2, APRT, QSOX2, USP47, THOC1, RPAP2, TRAPPC6A, EML4, ZNF573, HLTF, ZNF827, LAMP5, MRPS31, TNFAIP8, FKBP4, LYRM4, IGHV5-78 , BIN1, AUTS2, ATP5L, NDFIP1, SEPT1, RPL36, SALL2, CD5, ZNF239, ORC5, MRPS18B, CCT7, C12orf65, MSH6, RNMT, ITM2A, PPP1R3E, TSPYL2, NET1, IGIP, PPAN, TOMM20, RCN2, RPS20, PNN, NUP88, LIG1, C22orf32, TRMT5, GRAMD3, LUC7L3, SNRPF, PARM1, UBAP2, WDR74, HNRNPH3, PPAPDC1B, NFKB1, NSUN6, CXCR6, C21orf7, NARG2, PEBP1, MARCH3, NKRF, C21orf2, ST13, NAP1L1, TLR10, HNRNPA0, B3GALT2, RPUSD4, CECR5, DDX31, MST4, BYSL, TIMM9, PRMT6, EPHA4, RDH14, CNOT7, RPS14, TUT1, EIF2AK3, KDSR, SFXN2, FARSB, PIGN, ATP5G2, ZNF470, TCL1A, TFAM, ZMYND19, VIPR1, WDR89, NGDN, CTPS2, POLR3H, ZBED5, UBE2D2, DOLPP1, TSR1, MCCC1, SEL1L3, TMEM42, KLHL34, VRK1, DAP3, BLMH, HIST1H3H, MRPS33, NSG1, RNMTL1, FUT8, NAA25, BCL9, PTPLAD1, CLPP, ERCC6L2, TMEM220, CUTA, UTP14A, NOP56, EXOSC2, ZEB1, MS4A1, MPHOSPH10, NOL6, RPL13AP3, ADORA2A, DPPA4, ZNF428, SRSF5, RIOK1, ZNF830, SLC24A6, FXR1, AFG3L2, MRTO4, P2RY10, KIAA0922, KLF12, ANP32B, PNOC, MAK16, TASP1, RPL28, PSAT1, YTHDF1, MIS12, EFHD1, PDCL3, NUDT9, UQCC, IL27RA, PAK1IP1, EIF4A3, FGF9, BIVM, ATP5G1, FAM165B, U2AF1L4, LARP1B, DENR, MRPS9, C11orf31, CYP2E1, IPCEF1, RRN3P1, DCK, SRSF1, TGDS, HSPE1, NACA, CXXC5, PPIL1, CEP135, PIGP, PPWD1, NSUN2, NCBP2, KRI1, USP44, NIN, ZNF175, GRWD1, LOC338799, THAP11, LARP4, RPS29, TDP1, VARS, CYCS, IFT74, HNRNPH1, RPL17, BAG3, TOMM70A, SPECC1L, HSPA9, ERGIC3, CRELD1, TBRG4, SEPT9, LEO1, HIBADH, MTERFD1, ANAPC16, AHSA1, MZF1, SLC25A15, ABCF1, QTRT1, SIRT4, GNB2L1, CBY1, GYPE, CTPS1, XPO4, IPO4, LPHN1, CACYBP, SUPT3H, OPN3, FKBP11, MTR, ZNF721, ZNF260, AKAP7, ILF2, UBTF, NOB1, MPHOSPH9, LAPTM4B, BOD1, ARHGAP15, RBM14, ITPKB, RPSA, EWSR1, RPS6KB1, BCCIP, TFB2M, DHX15, C1QTNF6, GORASP2, FAM43A, PRKRA, MBTPS1, MGC57346, RUVBL1, STIP1, CD37, POLE3, FAM220A, SFXN4, RPP21, PRKACB, LCLAT1, MATN2, PJA1, ATAD1, CCDC57, POLRMT, TCERG1, PODXL2, DCAF4, RBFA, ALG13, CHORDC1, PRMT1, TMEM5, RPL11, ZNF581, DZIP3, CLNS1A, SURF6, AGPAT5, RPAIN, PTPN13, PRPF19, EBPL, TCTN1, SMARCC1, UFM1, FBXO15, EXOSC8, TMEM230, ZNF121, METTL25, RPL27A, ZNF234, UBE2G2, NSUN5, CYB5A, KIAA0586, USP13, CCR9, HSPA8, RASGRP1, SOX8, GLMN, METAP2, DHPS, CELA1, PPIL3, TGFBR2, SRRM1, TAF4, CEPT1, MTRF1, GTF3A, MSH2, TRIM44, TSPYL5, KANSL2, PCSK5, RIOK2, TMEM209, CSNK2A2, ACADM, DDX1, EPM2A, PGRMC2, RPL30, TSGA10, ZNF706, SMARCAD1, ZNF10, NGLY1, RPL35, TRA2B, EI24, N6AMT2, ZNF800, CD72, ZNF189, MORC2-AS1, COLQ, GART, TPT1, SUMF2, SRFBP1, CDC16, RBM11, HDAC1, TCTN3, MOAP1, SAAL1, ZRANB2, NUP62, MIR17HG, CCDC34, ABCF2, RPL23AP64, LARS2, OXCT1, MAGEE1, ZNF337, CCT6A, BCS1L, RPS18, SH2D3A, HIVEP2, ZNF567, CCNC, FBXO10, IPO11, PRIM1, HIC2, UBXN1, MRPL35, POLG2, LSM4

## Exact analysis specification frozen before running

The unit is the donor, never the cell. The analysis reuses `run_family_c_p1_signature` from the Block 5 pipeline without changing its scoring or random-set code. Full-depth raw counts are aggregated to donor pseudobulk, the human technical mask is applied before normalization, and expression is log2(CPM + 1). There is no downsampling or rarefaction. The mask excludes uppercase MT-, RPL/RPS, MRPL/MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D, plus ERCC features through the Block 5 extended mask.

Only donors aged 30-plus are scored, in both cohorts and both arms. For each signature, cohort, arm, and cell-type unit, the composite score is mean published age-up log2CPM minus mean published age-down log2CPM. The OLS predictor is age_per_decade. Both arms include sex_male and log_lib. The batch-adjusted arm adds Gao study fixed effects for seven studies (six indicator columns) or Yazar pool fixed effects, using the Block 5 pool construction and minimum residual df of 10. The unadjusted arm removes only batch covariates. Gao uses minimum residual df 2. The expression filter is the Block 5 filter: at least 1 CPM in at least half of donors, at least 200 genes in the expression pool, and at least one detected up and down signature gene.

Shared eligibility floors are at least 20 cells per donor per unit and at least 6 donors per unit. Units below either floor are excluded rather than treated as nulls. For each signature and cohort, only units eligible in both arms enter the batch contrast, so arm comparisons contain the same units. Yazar uses pool fixed effects and excludes the peripheral blood mononuclear cell doublet label as in Block 5.

The matched random-set comparison uses exactly 2000 iterations, matched separately on up and down set sizes and expression deciles, with fixed seed 20261105 for every fit. Composite empirical p-values use the Block 5 centered two-sided matched-null statistic. A unit passes only when composite BH q < 0.05 and empirical BH q < 0.05. MDE80 is (t(0.975, df) + t(0.80, df)) times the composite standard error and is reported for every tested unit, including every null.

Coverage denominator is the number of unique published human symbols used. N detected is the number present in the cohort post-mask gene universe. N scored is the number detected in at least one eligible unit expression pool. Per-unit scored counts are also reported. Coverage below 50 percent makes that signature UNTESTABLE in that cohort, not null. All four signatures are human, so mouse one-to-one ortholog mapping is not applicable; n_genes_mapped equals the unique symbols used.

The optional paper human cohorts were not added. The complete MSSM/HBCC brain, kidney, and lung extension does not have one already cached batch-contrast pipeline with these exact covariates across all those cohorts. No partial secondary arm or new download was substituted. The two required primary cohorts use only existing local data.

A split-half technical null is not part of this transport design; the required design-specific null is the 2000-iteration expression-matched random-set comparison. Age bins are not used because the prespecified model is continuous age per decade.

## Multiplicity families frozen before correction

BH correction is within each signature, cohort, and arm across the frozen units, separately for composite and empirical tests. There are 32 families: 4 signatures x 2 cohorts x 2 arms x 2 test types. An UNTESTABLE coverage slice has family size zero and no q-values. Family sizes are identical for composite and empirical tests within a slice. Family sizes at 30-plus are in the table below and in the registry CSV; a family can shrink where units become ineligible at 30-plus.

| Signature | Cohort | Coverage | Detected | Scored | Status | Units per arm and test family |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| Harries2011_6gene | Gao_Heart | 100.00% | 6 | 6 | TESTABLE | 13 |
| Harries2011_6gene | Yazar_Blood | 100.00% | 6 | 6 | TESTABLE | 17 |
| deMagalhaes2009_73gene | Gao_Heart | 90.00% | 63 | 61 | TESTABLE | 14 |
| deMagalhaes2009_73gene | Yazar_Blood | 100.00% | 70 | 49 | TESTABLE | 17 |
| Sood2015_150probe | Gao_Heart | 78.38% | 116 | 113 | TESTABLE | 14 |
| Sood2015_150probe | Yazar_Blood | 84.46% | 125 | 94 | TESTABLE | 17 |
| Peters2015_1497gene | Gao_Heart | 85.30% | 1277 | 1253 | TESTABLE | 14 |
| Peters2015_1497gene | Yazar_Blood | 89.38% | 1338 | 1296 | TESTABLE | 17 |

## Verdict rule frozen before running

A signature SURVIVES if at least one unit passes both tests after batch adjustment in either testable primary cohort. It is LOST_WITH_BATCH only if it has at least one unadjusted both-test pass in each primary cohort and zero adjusted both-test passes in both cohorts. Other outcomes are recorded as NO_MAJORITY_CATEGORY or PARTIALLY_UNTESTABLE and are not mislabeled as batch losses. PROVES requires a majority of the four benchmarked signatures to be LOST_WITH_BATCH. DOES NOT PROVE requires a majority to SURVIVE. INCONCLUSIVE applies if most signatures are coverage-untestable, results split evenly, or neither majority rule is met. The verdict names every surviving and non-surviving signature. Four signatures cannot prove a claim about ALL published signatures; the verdict is about this frozen benchmarked set only.

## Results

The pre-analysis section above was written before model fitting. This section contains the completed batch contrast.

### Batch-contrast summary

Both pass means composite BH q < 0.05 and empirical BH q < 0.05 in the frozen slice. N is the number of tested units.

| Signature | Cohort | Arm | N units | Composite passes | Empirical passes | Both pass |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Harries2011_6gene | Gao_Heart | unadjusted | 13 | 0 | 0 | 0 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | 13 | 0 | 0 | 0 |
| Harries2011_6gene | Yazar_Blood | unadjusted | 17 | 8 | 7 | 7 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | 17 | 5 | 7 | 5 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | 14 | 0 | 0 | 0 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | 14 | 0 | 1 | 0 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | 17 | 1 | 0 | 0 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | 17 | 1 | 0 | 0 |
| Sood2015_150probe | Gao_Heart | unadjusted | 14 | 3 | 3 | 3 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | 14 | 0 | 2 | 0 |
| Sood2015_150probe | Yazar_Blood | unadjusted | 17 | 2 | 2 | 1 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | 17 | 1 | 0 | 0 |
| Peters2015_1497gene | Gao_Heart | unadjusted | 14 | 0 | 7 | 0 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | 14 | 0 | 8 | 0 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | 17 | 16 | 17 | 16 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | 17 | 16 | 17 | 16 |

### Per-signature outcome

| Signature | Unadjusted both-pass units, Gao | Adjusted, Gao | Unadjusted both-pass units, Yazar | Adjusted, Yazar | Classification |
| --- | ---: | ---: | ---: | ---: | --- |
| Harries2011_6gene | 0 | 0 | 7 | 5 | SURVIVES |
| deMagalhaes2009_73gene | 0 | 0 | 0 | 0 | NO_MAJORITY_CATEGORY |
| Sood2015_150probe | 3 | 0 | 1 | 0 | LOST_WITH_BATCH |
| Peters2015_1497gene | 0 | 0 | 16 | 16 | SURVIVES |

### Mature vs original: survived vs died vs untestable

Original counts are both-pass counts from the full-range benchmark (`results/followup-11pt-block10-benchmark-results.csv`, 246 tested rows). Mature counts are both-pass counts at 30-plus from this rerun. Died means the signature SURVIVED originally and has zero adjusted both-pass units at 30-plus in both cohorts while remaining testable. Untestable means coverage below 50 percent at 30-plus.

| Signature | Original Gao unadj | Original Gao adj | Original Yazar unadj | Original Yazar adj | Original class | Mature Gao unadj | Mature Gao adj | Mature Yazar unadj | Mature Yazar adj | Mature class | Change |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | --- |
| Harries2011_6gene | 0 | 0 | 7 | 4 | SURVIVES | 0 | 0 | 7 | 5 | SURVIVES | unchanged |
| deMagalhaes2009_73gene | 0 | 0 | 0 | 0 | NO_MAJORITY_CATEGORY | 0 | 0 | 0 | 0 | NO_MAJORITY_CATEGORY | unchanged |
| Sood2015_150probe | 1 | 1 | 0 | 0 | SURVIVES | 3 | 0 | 1 | 0 | LOST_WITH_BATCH | no longer survives at 30-plus |
| Peters2015_1497gene | 0 | 0 | 17 | 16 | SURVIVES | 0 | 0 | 16 | 16 | SURVIVES | unchanged |

### Units passing both tests after batch adjustment

| Signature | Cohort | Unit | Slope per decade | Composite q | Empirical q | N donors | MDE80 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| Harries2011_6gene | Yazar_Blood | All_Cells | 0.1985 | 7.511e-47 | 0.004248 | 932 | 0.03585 |
| Harries2011_6gene | Yazar_Blood | central memory CD4-positive, alpha-beta T cell | 0.09951 | 9.054e-13 | 0.004248 | 931 | 0.03721 |
| Harries2011_6gene | Yazar_Blood | effector memory CD8-positive, alpha-beta T cell | 0.1179 | 3.65e-05 | 0.004248 | 917 | 0.07387 |
| Harries2011_6gene | Yazar_Blood | naive thymus-derived CD4-positive, alpha-beta T cell | 0.1423 | 7.95e-18 | 0.004248 | 925 | 0.04409 |
| Harries2011_6gene | Yazar_Blood | naive thymus-derived CD8-positive, alpha-beta T cell | 0.09531 | 0.03276 | 0.04248 | 623 | 0.103 |
| Peters2015_1497gene | Yazar_Blood | All_Cells | 0.09675 | 2.599e-33 | 0.0006535 | 932 | 0.02117 |
| Peters2015_1497gene | Yazar_Blood | CD14-low, CD16-positive monocyte | 0.02265 | 0.003238 | 0.001699 | 277 | 0.02065 |
| Peters2015_1497gene | Yazar_Blood | CD14-positive monocyte | 0.01424 | 0.006072 | 0.0006535 | 464 | 0.01404 |
| Peters2015_1497gene | Yazar_Blood | CD4-positive, alpha-beta cytotoxic T cell | 0.03616 | 0.001376 | 0.0006535 | 261 | 0.03017 |
| Peters2015_1497gene | Yazar_Blood | central memory CD4-positive, alpha-beta T cell | 0.02243 | 2.879e-16 | 0.0006535 | 931 | 0.007347 |
| Peters2015_1497gene | Yazar_Blood | central memory CD8-positive, alpha-beta T cell | 0.02354 | 0.02494 | 0.002655 | 261 | 0.02903 |
| Peters2015_1497gene | Yazar_Blood | effector memory CD4-positive, alpha-beta T cell | 0.02513 | 1.533e-05 | 0.0006535 | 658 | 0.01557 |
| Peters2015_1497gene | Yazar_Blood | effector memory CD8-positive, alpha-beta T cell | 0.03297 | 6.994e-11 | 0.0006535 | 917 | 0.01361 |
| Peters2015_1497gene | Yazar_Blood | gamma-delta T cell | 0.05949 | 0.0003528 | 0.0006535 | 292 | 0.04438 |
| Peters2015_1497gene | Yazar_Blood | memory B cell | 0.01776 | 0.00741 | 0.0006535 | 573 | 0.0181 |
| Peters2015_1497gene | Yazar_Blood | mucosal-associated invariant T cell | 0.05238 | 0.01415 | 0.001214 | 74 | 0.05753 |
| Peters2015_1497gene | Yazar_Blood | naive thymus-derived CD4-positive, alpha-beta T cell | 0.02484 | 2.808e-27 | 0.0006535 | 925 | 0.006098 |
| Peters2015_1497gene | Yazar_Blood | naive thymus-derived CD8-positive, alpha-beta T cell | 0.03071 | 3.321e-07 | 0.0006535 | 623 | 0.01614 |
| Peters2015_1497gene | Yazar_Blood | natural killer cell | 0.019 | 2.403e-11 | 0.0006535 | 928 | 0.007625 |
| Peters2015_1497gene | Yazar_Blood | regulatory T cell | 0.06206 | 4.556e-09 | 0.0006535 | 574 | 0.02834 |
| Peters2015_1497gene | Yazar_Blood | transitional stage B cell | 0.03099 | 0.0001038 | 0.0006535 | 456 | 0.02133 |

### Complete per-unit results and MDE80

Every tested unit is shown so each null has its MDE80 alongside it. UNTESTABLE coverage slices have no model or MDE80 and remain listed in the frozen coverage table.

| Signature | Cohort | Arm | Unit | Slope per decade | Composite q | Empirical q | Both pass | N donors | MDE80 |
| --- | --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Adipocyte | 0.3797 | 0.6425 | 0.7081 | False | 44 | 0.7587 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | All_LV_Nuclei | -0.06365 | 0.6425 | 0.7081 | False | 109 | 0.1389 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Cardiomyocyte | -0.008305 | 0.9816 | 0.978 | False | 107 | 0.2205 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Endocardial | -0.2174 | 0.4914 | 0.7081 | False | 80 | 0.3424 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Endothelial | -0.07138 | 0.8693 | 0.7546 | False | 103 | 0.2497 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Epicardial | -0.4077 | 0.8693 | 0.7546 | False | 22 | 1.254 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Fibroblast | -0.05125 | 0.8693 | 0.7546 | False | 109 | 0.2398 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | LEC | -0.07126 | 0.9816 | 0.9771 | False | 59 | 0.7533 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Lymphoid | -0.03197 | 0.9816 | 0.978 | False | 105 | 0.4892 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Myeloid | 0.001747 | 0.9816 | 0.978 | False | 109 | 0.2138 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Neuronal | -0.1374 | 0.8693 | 0.8323 | False | 95 | 0.72 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | Pericyte | 0.2646 | 0.06568 | 0.2014 | False | 109 | 0.2611 |
| Harries2011_6gene | Gao_Heart | batch_adjusted | vSMC | 0.07546 | 0.8693 | 0.8323 | False | 107 | 0.4079 |
| Harries2011_6gene | Gao_Heart | unadjusted | Adipocyte | 0.3109 | 0.9455 | 0.9965 | False | 44 | 0.8816 |
| Harries2011_6gene | Gao_Heart | unadjusted | All_LV_Nuclei | -0.05217 | 0.9455 | 0.9965 | False | 109 | 0.216 |
| Harries2011_6gene | Gao_Heart | unadjusted | Cardiomyocyte | -0.06567 | 0.9455 | 0.9965 | False | 107 | 0.3626 |
| Harries2011_6gene | Gao_Heart | unadjusted | Endocardial | -0.1943 | 0.9455 | 0.9965 | False | 80 | 0.334 |
| Harries2011_6gene | Gao_Heart | unadjusted | Endothelial | -0.03337 | 0.9455 | 0.9965 | False | 103 | 0.2425 |
| Harries2011_6gene | Gao_Heart | unadjusted | Epicardial | -0.02586 | 0.9455 | 0.9965 | False | 22 | 1.098 |
| Harries2011_6gene | Gao_Heart | unadjusted | Fibroblast | -0.01606 | 0.9455 | 0.9965 | False | 109 | 0.3654 |
| Harries2011_6gene | Gao_Heart | unadjusted | LEC | 0.01678 | 0.9455 | 0.9965 | False | 59 | 0.6973 |
| Harries2011_6gene | Gao_Heart | unadjusted | Lymphoid | 0.02761 | 0.9455 | 0.9965 | False | 105 | 0.5432 |
| Harries2011_6gene | Gao_Heart | unadjusted | Myeloid | 0.005481 | 0.9455 | 0.9965 | False | 109 | 0.2181 |
| Harries2011_6gene | Gao_Heart | unadjusted | Neuronal | -0.07983 | 0.9455 | 0.9965 | False | 95 | 0.7623 |
| Harries2011_6gene | Gao_Heart | unadjusted | Pericyte | 0.2136 | 0.9455 | 0.4678 | False | 109 | 0.4538 |
| Harries2011_6gene | Gao_Heart | unadjusted | vSMC | 0.04417 | 0.9455 | 0.9965 | False | 107 | 0.4865 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | All_Cells | 0.1985 | 7.511e-47 | 0.004248 | True | 932 | 0.03585 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | CD14-low, CD16-positive monocyte | 0.07576 | 0.671 | 0.762 | False | 277 | 0.3255 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | CD14-positive monocyte | 0.08218 | 0.4977 | 0.469 | False | 464 | 0.2061 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.1036 | 0.6114 | 0.762 | False | 261 | 0.3701 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | central memory CD4-positive, alpha-beta T cell | 0.09951 | 9.054e-13 | 0.004248 | True | 931 | 0.03721 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | central memory CD8-positive, alpha-beta T cell | 0.1402 | 0.1461 | 0.2899 | False | 261 | 0.2157 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | effector memory CD4-positive, alpha-beta T cell | 0.02375 | 0.6714 | 0.7975 | False | 658 | 0.1122 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | effector memory CD8-positive, alpha-beta T cell | 0.1179 | 3.65e-05 | 0.004248 | True | 917 | 0.07387 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | gamma-delta T cell | 0.1708 | 0.1228 | 0.04248 | False | 292 | 0.2445 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | memory B cell | 0.1432 | 0.05625 | 0.04976 | False | 573 | 0.172 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | mucosal-associated invariant T cell | 0.1202 | 0.4977 | 0.469 | False | 74 | 0.3247 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | naive B cell | 0.001865 | 0.9688 | 0.9015 | False | 796 | 0.1336 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | naive thymus-derived CD4-positive, alpha-beta T cell | 0.1423 | 7.95e-18 | 0.004248 | True | 925 | 0.04409 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.09531 | 0.03276 | 0.04248 | True | 623 | 0.103 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | natural killer cell | 0.02589 | 0.4977 | 0.762 | False | 928 | 0.07328 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | regulatory T cell | 0.01992 | 0.7676 | 0.7975 | False | 574 | 0.1574 |
| Harries2011_6gene | Yazar_Blood | batch_adjusted | transitional stage B cell | 0.02844 | 0.7537 | 0.7975 | False | 456 | 0.1843 |
| Harries2011_6gene | Yazar_Blood | unadjusted | All_Cells | 0.2277 | 3.172e-65 | 0.006372 | True | 932 | 0.03416 |
| Harries2011_6gene | Yazar_Blood | unadjusted | CD14-low, CD16-positive monocyte | 0.09297 | 0.4251 | 0.6096 | False | 287 | 0.2677 |
| Harries2011_6gene | Yazar_Blood | unadjusted | CD14-positive monocyte | 0.01495 | 0.8123 | 0.9735 | False | 468 | 0.1766 |
| Harries2011_6gene | Yazar_Blood | unadjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.0612 | 0.6019 | 0.6096 | False | 275 | 0.2998 |
| Harries2011_6gene | Yazar_Blood | unadjusted | central memory CD4-positive, alpha-beta T cell | 0.1275 | 7.051e-24 | 0.006372 | True | 931 | 0.0339 |
| Harries2011_6gene | Yazar_Blood | unadjusted | central memory CD8-positive, alpha-beta T cell | 0.2078 | 0.0007202 | 0.04127 | True | 270 | 0.1576 |
| Harries2011_6gene | Yazar_Blood | unadjusted | effector memory CD4-positive, alpha-beta T cell | 0.0382 | 0.4043 | 0.6096 | False | 658 | 0.1003 |
| Harries2011_6gene | Yazar_Blood | unadjusted | effector memory CD8-positive, alpha-beta T cell | 0.1211 | 1.315e-06 | 0.008496 | True | 917 | 0.06589 |
| Harries2011_6gene | Yazar_Blood | unadjusted | gamma-delta T cell | 0.2945 | 0.0001155 | 0.006372 | True | 298 | 0.1966 |
| Harries2011_6gene | Yazar_Blood | unadjusted | memory B cell | 0.1543 | 0.005902 | 0.03115 | True | 573 | 0.1442 |
| Harries2011_6gene | Yazar_Blood | unadjusted | mucosal-associated invariant T cell | 0.1066 | 0.3321 | 0.32 | False | 93 | 0.2314 |
| Harries2011_6gene | Yazar_Blood | unadjusted | naive B cell | 0.04618 | 0.4043 | 0.6096 | False | 796 | 0.1191 |
| Harries2011_6gene | Yazar_Blood | unadjusted | naive thymus-derived CD4-positive, alpha-beta T cell | 0.1888 | 2.296e-35 | 0.006372 | True | 925 | 0.04028 |
| Harries2011_6gene | Yazar_Blood | unadjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.1076 | 0.001041 | 0.12 | False | 623 | 0.08529 |
| Harries2011_6gene | Yazar_Blood | unadjusted | natural killer cell | 0.05254 | 0.05352 | 0.4146 | False | 928 | 0.0671 |
| Harries2011_6gene | Yazar_Blood | unadjusted | regulatory T cell | 0.04536 | 0.4251 | 0.6096 | False | 576 | 0.1361 |
| Harries2011_6gene | Yazar_Blood | unadjusted | transitional stage B cell | 0.03878 | 0.55 | 0.6096 | False | 456 | 0.1559 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Adipocyte | -0.0374 | 0.4184 | 0.004664 | False | 44 | 0.06437 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | All_LV_Nuclei | -0.0141 | 0.4477 | 0.001749 | False | 109 | 0.03312 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Cardiomyocyte | -0.007341 | 0.5155 | 0.3371 | False | 107 | 0.02581 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Endocardial | -0.02645 | 0.3531 | 0.01799 | False | 80 | 0.03776 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Endothelial | -0.001035 | 0.9233 | 0.5812 | False | 103 | 0.02492 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Epicardial | -0.003248 | 0.9233 | 0.94 | False | 22 | 0.09905 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Fibroblast | -0.01297 | 0.4184 | 0.001749 | False | 109 | 0.02419 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | LEC | -0.02317 | 0.4893 | 0.08956 | False | 59 | 0.06779 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Lymphoid | 0.01652 | 0.4893 | 0.03323 | False | 105 | 0.04973 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Mast | -0.01548 | 0.5155 | 0.3866 | False | 69 | 0.05695 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Myeloid | -0.01624 | 0.4477 | 0.001749 | False | 109 | 0.03927 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Neuronal | -0.03638 | 0.2481 | 0.001749 | False | 95 | 0.04263 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | Pericyte | -0.01704 | 0.4184 | 0.002799 | False | 109 | 0.03319 |
| Peters2015_1497gene | Gao_Heart | batch_adjusted | vSMC | -0.01511 | 0.4477 | 0.08956 | False | 107 | 0.03741 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Adipocyte | -0.03774 | 0.4007 | 0.005597 | False | 44 | 0.06155 |
| Peters2015_1497gene | Gao_Heart | unadjusted | All_LV_Nuclei | -0.01267 | 0.5861 | 0.001749 | False | 109 | 0.04036 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Cardiomyocyte | -0.01395 | 0.4552 | 0.05685 | False | 107 | 0.0316 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Endocardial | -0.0164 | 0.5729 | 0.1018 | False | 80 | 0.04722 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Endothelial | -0.007208 | 0.6121 | 0.09025 | False | 103 | 0.02865 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Epicardial | -0.03644 | 0.4552 | 0.07541 | False | 22 | 0.08643 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Fibroblast | -0.01618 | 0.4007 | 0.001749 | False | 109 | 0.02487 |
| Peters2015_1497gene | Gao_Heart | unadjusted | LEC | -0.03148 | 0.4552 | 0.01866 | False | 59 | 0.07334 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Lymphoid | -0.005705 | 0.8388 | 0.4763 | False | 105 | 0.05734 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Mast | 0.001255 | 0.9638 | 0.9605 | False | 69 | 0.07841 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Myeloid | -0.009088 | 0.7869 | 0.03098 | False | 109 | 0.06101 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Neuronal | -0.05212 | 0.1397 | 0.001749 | False | 95 | 0.05609 |
| Peters2015_1497gene | Gao_Heart | unadjusted | Pericyte | -0.02826 | 0.4474 | 0.001749 | False | 109 | 0.05207 |
| Peters2015_1497gene | Gao_Heart | unadjusted | vSMC | -0.009617 | 0.6121 | 0.2659 | False | 107 | 0.03846 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | All_Cells | 0.09675 | 2.599e-33 | 0.0006535 | True | 932 | 0.02117 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | CD14-low, CD16-positive monocyte | 0.02265 | 0.003238 | 0.001699 | True | 277 | 0.02065 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | CD14-positive monocyte | 0.01424 | 0.006072 | 0.0006535 | True | 464 | 0.01404 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.03616 | 0.001376 | 0.0006535 | True | 261 | 0.03017 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | central memory CD4-positive, alpha-beta T cell | 0.02243 | 2.879e-16 | 0.0006535 | True | 931 | 0.007347 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | central memory CD8-positive, alpha-beta T cell | 0.02354 | 0.02494 | 0.002655 | True | 261 | 0.02903 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | effector memory CD4-positive, alpha-beta T cell | 0.02513 | 1.533e-05 | 0.0006535 | True | 658 | 0.01557 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | effector memory CD8-positive, alpha-beta T cell | 0.03297 | 6.994e-11 | 0.0006535 | True | 917 | 0.01361 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | gamma-delta T cell | 0.05949 | 0.0003528 | 0.0006535 | True | 292 | 0.04438 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | memory B cell | 0.01776 | 0.00741 | 0.0006535 | True | 573 | 0.0181 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | mucosal-associated invariant T cell | 0.05238 | 0.01415 | 0.001214 | True | 74 | 0.05753 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | naive B cell | 0.008058 | 0.06906 | 0.01099 | False | 796 | 0.01242 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | naive thymus-derived CD4-positive, alpha-beta T cell | 0.02484 | 2.808e-27 | 0.0006535 | True | 925 | 0.006098 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.03071 | 3.321e-07 | 0.0006535 | True | 623 | 0.01614 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | natural killer cell | 0.019 | 2.403e-11 | 0.0006535 | True | 928 | 0.007625 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | regulatory T cell | 0.06206 | 4.556e-09 | 0.0006535 | True | 574 | 0.02834 |
| Peters2015_1497gene | Yazar_Blood | batch_adjusted | transitional stage B cell | 0.03099 | 0.0001038 | 0.0006535 | True | 456 | 0.02133 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | All_Cells | 0.09594 | 6.977e-31 | 0.0005664 | True | 932 | 0.02207 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | CD14-low, CD16-positive monocyte | 0.01829 | 0.006518 | 0.001062 | True | 287 | 0.01862 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | CD14-positive monocyte | 0.02145 | 3.112e-06 | 0.0005664 | True | 468 | 0.0125 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.03379 | 0.0003526 | 0.0005664 | True | 275 | 0.02573 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | central memory CD4-positive, alpha-beta T cell | 0.02627 | 1.002e-20 | 0.0005664 | True | 931 | 0.007549 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | central memory CD8-positive, alpha-beta T cell | 0.02797 | 0.0008004 | 0.0005664 | True | 270 | 0.02295 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | effector memory CD4-positive, alpha-beta T cell | 0.02638 | 3.446e-07 | 0.0005664 | True | 658 | 0.01403 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | effector memory CD8-positive, alpha-beta T cell | 0.03118 | 1.463e-11 | 0.0005664 | True | 917 | 0.0125 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | gamma-delta T cell | 0.07442 | 6.568e-08 | 0.0005664 | True | 298 | 0.03677 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | memory B cell | 0.01864 | 0.0006207 | 0.0005664 | True | 573 | 0.01497 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | mucosal-associated invariant T cell | 0.06059 | 4.76e-06 | 0.0005664 | True | 93 | 0.0346 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | naive B cell | 0.007381 | 0.06829 | 0.003998 | False | 796 | 0.01134 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | naive thymus-derived CD4-positive, alpha-beta T cell | 0.0306 | 5.382e-37 | 0.0005664 | True | 925 | 0.006339 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.03222 | 2.611e-10 | 0.0005664 | True | 623 | 0.01376 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | natural killer cell | 0.01991 | 2.573e-13 | 0.0005664 | True | 928 | 0.007323 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | regulatory T cell | 0.06895 | 1.463e-11 | 0.0005664 | True | 576 | 0.02739 |
| Peters2015_1497gene | Yazar_Blood | unadjusted | transitional stage B cell | 0.03077 | 3.112e-06 | 0.0005664 | True | 456 | 0.01794 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Adipocyte | 0.1587 | 0.06576 | 0.04198 | False | 44 | 0.1516 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | All_LV_Nuclei | 0.0409 | 0.08839 | 0.06997 | False | 109 | 0.04851 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Cardiomyocyte | 0.04614 | 0.07461 | 0.1707 | False | 107 | 0.05014 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Endocardial | 0.01979 | 0.7353 | 0.5695 | False | 80 | 0.1053 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Endothelial | 0.00242 | 0.9769 | 0.6636 | False | 103 | 0.08533 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Epicardial | -0.06289 | 0.4905 | 0.5845 | False | 22 | 0.1811 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Fibroblast | 0.02503 | 0.4984 | 0.3647 | False | 109 | 0.07636 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | LEC | 0.08531 | 0.3016 | 0.1707 | False | 59 | 0.1488 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Lymphoid | 0.04927 | 0.4516 | 0.5695 | False | 105 | 0.1156 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Mast | 0.08393 | 0.4516 | 0.3647 | False | 69 | 0.1926 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Myeloid | 0.03185 | 0.4516 | 0.3067 | False | 109 | 0.07923 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Neuronal | -0.001077 | 0.9769 | 0.8281 | False | 95 | 0.105 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | Pericyte | 0.01445 | 0.7353 | 0.611 | False | 109 | 0.08466 |
| Sood2015_150probe | Gao_Heart | batch_adjusted | vSMC | 0.08253 | 0.1288 | 0.04198 | False | 107 | 0.1103 |
| Sood2015_150probe | Gao_Heart | unadjusted | Adipocyte | 0.1776 | 0.02819 | 0.009329 | True | 44 | 0.1555 |
| Sood2015_150probe | Gao_Heart | unadjusted | All_LV_Nuclei | 0.0561 | 0.02819 | 0.006997 | True | 109 | 0.05643 |
| Sood2015_150probe | Gao_Heart | unadjusted | Cardiomyocyte | 0.0485 | 0.2986 | 0.1729 | False | 107 | 0.07895 |
| Sood2015_150probe | Gao_Heart | unadjusted | Endocardial | 0.03678 | 0.4788 | 0.4249 | False | 80 | 0.1017 |
| Sood2015_150probe | Gao_Heart | unadjusted | Endothelial | 0.01787 | 0.6273 | 0.681 | False | 103 | 0.08175 |
| Sood2015_150probe | Gao_Heart | unadjusted | Epicardial | -0.01244 | 0.8345 | 0.7578 | False | 22 | 0.1739 |
| Sood2015_150probe | Gao_Heart | unadjusted | Fibroblast | 0.03954 | 0.3773 | 0.07171 | False | 109 | 0.08109 |
| Sood2015_150probe | Gao_Heart | unadjusted | LEC | 0.06708 | 0.3773 | 0.2737 | False | 59 | 0.1437 |
| Sood2015_150probe | Gao_Heart | unadjusted | Lymphoid | 0.08167 | 0.3498 | 0.1203 | False | 105 | 0.1493 |
| Sood2015_150probe | Gao_Heart | unadjusted | Mast | 0.06371 | 0.4788 | 0.4249 | False | 69 | 0.1985 |
| Sood2015_150probe | Gao_Heart | unadjusted | Myeloid | 0.03362 | 0.4788 | 0.1729 | False | 109 | 0.107 |
| Sood2015_150probe | Gao_Heart | unadjusted | Neuronal | -0.01658 | 0.767 | 0.8196 | False | 95 | 0.1269 |
| Sood2015_150probe | Gao_Heart | unadjusted | Pericyte | 0.03397 | 0.4788 | 0.2923 | False | 109 | 0.1021 |
| Sood2015_150probe | Gao_Heart | unadjusted | vSMC | 0.1097 | 0.02819 | 0.009329 | True | 107 | 0.1107 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | All_Cells | 0.01773 | 0.0005743 | 0.5369 | False | 932 | 0.01193 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | CD14-low, CD16-positive monocyte | -0.01332 | 0.9546 | 0.9745 | False | 277 | 0.06563 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | CD14-positive monocyte | -0.02782 | 0.2257 | 0.4333 | False | 464 | 0.03948 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.005702 | 0.9546 | 0.9745 | False | 261 | 0.1021 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | central memory CD4-positive, alpha-beta T cell | -0.008997 | 0.4484 | 0.6117 | False | 931 | 0.01787 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | central memory CD8-positive, alpha-beta T cell | 0.02656 | 0.9546 | 0.9745 | False | 261 | 0.1214 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | effector memory CD4-positive, alpha-beta T cell | -0.04367 | 0.2257 | 0.4333 | False | 658 | 0.06324 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | effector memory CD8-positive, alpha-beta T cell | 0.02554 | 0.05941 | 0.4333 | False | 917 | 0.02649 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | gamma-delta T cell | 0.01662 | 0.9546 | 0.9745 | False | 292 | 0.08481 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | memory B cell | -0.01022 | 0.9546 | 0.9745 | False | 573 | 0.06073 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | mucosal-associated invariant T cell | -0.008647 | 0.9546 | 0.9745 | False | 74 | 0.1933 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | naive B cell | 0.005204 | 0.9546 | 0.9745 | False | 796 | 0.0423 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | naive thymus-derived CD4-positive, alpha-beta T cell | -0.0001296 | 0.9894 | 0.9745 | False | 925 | 0.02737 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.004018 | 0.9546 | 0.9745 | False | 623 | 0.06374 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | natural killer cell | -0.01269 | 0.3708 | 0.5369 | False | 928 | 0.02219 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | regulatory T cell | -0.004221 | 0.9546 | 0.9745 | False | 574 | 0.07251 |
| Sood2015_150probe | Yazar_Blood | batch_adjusted | transitional stage B cell | 0.01718 | 0.9546 | 0.9745 | False | 456 | 0.0672 |
| Sood2015_150probe | Yazar_Blood | unadjusted | All_Cells | 0.00723 | 0.4005 | 0.931 | False | 932 | 0.01459 |
| Sood2015_150probe | Yazar_Blood | unadjusted | CD14-low, CD16-positive monocyte | -0.04842 | 0.06537 | 0.04673 | False | 287 | 0.05598 |
| Sood2015_150probe | Yazar_Blood | unadjusted | CD14-positive monocyte | -0.0385 | 0.01044 | 0.0623 | False | 468 | 0.03323 |
| Sood2015_150probe | Yazar_Blood | unadjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.03472 | 0.4829 | 0.6637 | False | 275 | 0.08066 |
| Sood2015_150probe | Yazar_Blood | unadjusted | central memory CD4-positive, alpha-beta T cell | -0.01473 | 0.06537 | 0.3738 | False | 931 | 0.01647 |
| Sood2015_150probe | Yazar_Blood | unadjusted | central memory CD8-positive, alpha-beta T cell | -0.009664 | 0.99 | 0.931 | False | 270 | 0.08774 |
| Sood2015_150probe | Yazar_Blood | unadjusted | effector memory CD4-positive, alpha-beta T cell | -0.06512 | 0.01044 | 0.04673 | True | 658 | 0.05469 |
| Sood2015_150probe | Yazar_Blood | unadjusted | effector memory CD8-positive, alpha-beta T cell | 0.01992 | 0.06537 | 0.6637 | False | 917 | 0.02398 |
| Sood2015_150probe | Yazar_Blood | unadjusted | gamma-delta T cell | 0.01512 | 0.9017 | 0.931 | False | 298 | 0.06765 |
| Sood2015_150probe | Yazar_Blood | unadjusted | memory B cell | 0.00854 | 0.9197 | 0.931 | False | 573 | 0.05167 |
| Sood2015_150probe | Yazar_Blood | unadjusted | mucosal-associated invariant T cell | 0.003826 | 0.999 | 0.931 | False | 93 | 0.1215 |
| Sood2015_150probe | Yazar_Blood | unadjusted | naive B cell | -0.006053 | 0.9197 | 0.725 | False | 796 | 0.03731 |
| Sood2015_150probe | Yazar_Blood | unadjusted | naive thymus-derived CD4-positive, alpha-beta T cell | -0.0005878 | 0.999 | 0.931 | False | 925 | 0.02358 |
| Sood2015_150probe | Yazar_Blood | unadjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 2.187e-05 | 0.999 | 0.931 | False | 623 | 0.05102 |
| Sood2015_150probe | Yazar_Blood | unadjusted | natural killer cell | -0.00668 | 0.6886 | 0.6637 | False | 928 | 0.02065 |
| Sood2015_150probe | Yazar_Blood | unadjusted | regulatory T cell | 0.003342 | 0.999 | 0.931 | False | 576 | 0.06269 |
| Sood2015_150probe | Yazar_Blood | unadjusted | transitional stage B cell | 0.04572 | 0.06537 | 0.3806 | False | 456 | 0.0563 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Adipocyte | -0.1255 | 0.3121 | 0.4285 | False | 44 | 0.1736 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | All_LV_Nuclei | 0.007474 | 0.8741 | 0.7238 | False | 109 | 0.08852 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Cardiomyocyte | -0.01034 | 0.8741 | 0.8859 | False | 107 | 0.06568 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Endocardial | -0.02746 | 0.8741 | 0.7638 | False | 80 | 0.1221 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Endothelial | -0.01552 | 0.8741 | 0.7238 | False | 103 | 0.092 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Epicardial | -0.1532 | 0.5261 | 0.3533 | False | 22 | 0.2724 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Fibroblast | 0.01096 | 0.8741 | 0.7238 | False | 109 | 0.107 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | LEC | -0.04009 | 0.8741 | 0.7238 | False | 59 | 0.1959 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Lymphoid | 0.04639 | 0.8741 | 0.7238 | False | 105 | 0.1486 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Mast | -0.0491 | 0.8741 | 0.7238 | False | 69 | 0.1993 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Myeloid | -0.01034 | 0.8741 | 0.7238 | False | 109 | 0.1085 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Neuronal | -0.07313 | 0.8433 | 0.4285 | False | 95 | 0.1755 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | Pericyte | -0.06971 | 0.3121 | 0.04198 | False | 109 | 0.08887 |
| deMagalhaes2009_73gene | Gao_Heart | batch_adjusted | vSMC | -0.002375 | 0.9581 | 0.9255 | False | 107 | 0.1277 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Adipocyte | -0.04192 | 0.8291 | 0.639 | False | 44 | 0.3503 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | All_LV_Nuclei | 0.03151 | 0.8291 | 0.1119 | False | 109 | 0.1502 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Cardiomyocyte | -0.00771 | 0.8291 | 0.9985 | False | 107 | 0.06564 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Endocardial | 0.07854 | 0.8291 | 0.2041 | False | 80 | 0.2623 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Endothelial | 0.08138 | 0.8291 | 0.1119 | False | 103 | 0.2286 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Epicardial | -0.06812 | 0.8291 | 0.5156 | False | 22 | 0.3544 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Fibroblast | 0.02613 | 0.8291 | 0.2457 | False | 109 | 0.2353 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | LEC | 0.02971 | 0.8683 | 0.7825 | False | 59 | 0.5085 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Lymphoid | 0.06281 | 0.8291 | 0.2457 | False | 105 | 0.2854 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Mast | 0.1791 | 0.8291 | 0.1399 | False | 69 | 0.4789 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Myeloid | 0.02888 | 0.8291 | 0.2457 | False | 109 | 0.2785 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Neuronal | -0.1234 | 0.8291 | 0.1119 | False | 95 | 0.3017 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | Pericyte | -0.04575 | 0.8291 | 0.2029 | False | 109 | 0.1474 |
| deMagalhaes2009_73gene | Gao_Heart | unadjusted | vSMC | 0.02825 | 0.8291 | 0.5425 | False | 107 | 0.1713 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | All_Cells | 0.03287 | 5.523e-08 | 0.1869 | False | 932 | 0.01541 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | CD14-low, CD16-positive monocyte | 0.03428 | 0.5401 | 0.5834 | False | 277 | 0.06265 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | CD14-positive monocyte | 0.0004667 | 0.9765 | 0.9319 | False | 464 | 0.04444 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | CD4-positive, alpha-beta cytotoxic T cell | -0.004611 | 0.9405 | 0.945 | False | 261 | 0.0898 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | central memory CD4-positive, alpha-beta T cell | 0.001415 | 0.9365 | 0.9319 | False | 931 | 0.01808 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | central memory CD8-positive, alpha-beta T cell | 0.009964 | 0.9029 | 0.8611 | False | 261 | 0.08563 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | effector memory CD4-positive, alpha-beta T cell | 0.01172 | 0.8928 | 0.8611 | False | 658 | 0.04644 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | effector memory CD8-positive, alpha-beta T cell | -0.01218 | 0.5401 | 0.8611 | False | 917 | 0.02471 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | gamma-delta T cell | 0.00965 | 0.9029 | 0.8611 | False | 292 | 0.0718 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | memory B cell | 0.01262 | 0.8928 | 0.8611 | False | 573 | 0.04722 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | mucosal-associated invariant T cell | 0.01849 | 0.8928 | 0.8611 | False | 74 | 0.08677 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | naive B cell | -0.01108 | 0.8928 | 0.8611 | False | 796 | 0.03973 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | naive thymus-derived CD4-positive, alpha-beta T cell | 0.01096 | 0.5401 | 0.8611 | False | 925 | 0.02338 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.009072 | 0.9029 | 0.8611 | False | 623 | 0.05405 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | natural killer cell | 0.01901 | 0.5401 | 0.5225 | False | 928 | 0.02879 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | regulatory T cell | 0.02653 | 0.5401 | 0.8611 | False | 574 | 0.05683 |
| deMagalhaes2009_73gene | Yazar_Blood | batch_adjusted | transitional stage B cell | -0.01176 | 0.8928 | 0.8611 | False | 456 | 0.05929 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | All_Cells | 0.03129 | 7.498e-07 | 0.5579 | False | 932 | 0.0159 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | CD14-low, CD16-positive monocyte | 0.0264 | 0.2677 | 0.5579 | False | 287 | 0.04921 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | CD14-positive monocyte | 0.0008633 | 0.9494 | 0.925 | False | 468 | 0.03817 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | CD4-positive, alpha-beta cytotoxic T cell | 0.01312 | 0.772 | 0.8313 | False | 275 | 0.06895 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | central memory CD4-positive, alpha-beta T cell | 0.002149 | 0.8153 | 0.925 | False | 931 | 0.01677 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | central memory CD8-positive, alpha-beta T cell | 0.03008 | 0.3214 | 0.6266 | False | 270 | 0.06423 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | effector memory CD4-positive, alpha-beta T cell | 0.009958 | 0.7061 | 0.8313 | False | 658 | 0.04125 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | effector memory CD8-positive, alpha-beta T cell | -0.01172 | 0.2677 | 0.785 | False | 917 | 0.02234 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | gamma-delta T cell | 0.001757 | 0.9494 | 0.945 | False | 298 | 0.05807 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | memory B cell | 0.02099 | 0.2677 | 0.5579 | False | 573 | 0.03938 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | mucosal-associated invariant T cell | 0.04536 | 0.2641 | 0.5579 | False | 93 | 0.06801 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | naive B cell | -0.02 | 0.2677 | 0.6362 | False | 796 | 0.03469 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | naive thymus-derived CD4-positive, alpha-beta T cell | 0.01885 | 0.08779 | 0.5579 | False | 925 | 0.02057 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | naive thymus-derived CD8-positive, alpha-beta T cell | 0.0162 | 0.4673 | 0.7901 | False | 623 | 0.04404 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | natural killer cell | 0.0195 | 0.2014 | 0.5579 | False | 928 | 0.02598 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | regulatory T cell | 0.02638 | 0.2677 | 0.6266 | False | 576 | 0.04877 |
| deMagalhaes2009_73gene | Yazar_Blood | unadjusted | transitional stage B cell | -0.008491 | 0.772 | 0.925 | False | 456 | 0.0503 |

## Verdict

**INCONCLUSIVE**. Neither majority rule was met for the four frozen signatures.

Surviving signatures: Harries2011_6gene, Peters2015_1497gene.

Non-surviving signatures and their categories: deMagalhaes2009_73gene = NO_MAJORITY_CATEGORY; Sood2015_150probe = LOST_WITH_BATCH.

Counts: 2 SURVIVES, 1 LOST_WITH_BATCH, 0 PARTIALLY_UNTESTABLE, and 1 NO_MAJORITY_CATEGORY.

These four signatures cannot prove a claim about ALL published signatures. This verdict is only about the frozen benchmarked set and its tested Gao and Yazar units. Absence of a passing unit is not proof of no biological age association; per-unit MDE80 values are supplied in the result CSV so nulls can be judged for power.

## Reproducibility and outputs

- `results/followup-mature-block10-benchmark-signatures.csv`: frozen citations, sources, exact up and down lists, coverage at 30-plus, and family sizes at 30-plus.
- `results/followup-mature-block10-benchmark-results.csv`: one row per signature, cohort, arm, and tested unit at 30-plus, plus cohort-level UNTESTABLE rows where required.
- `results/followup-mature-block10-donor-counts.csv`: donor counts at full range and at 30-plus for every cohort, arm, and unit.
- `scripts/run_mature_block10_published_benchmark.py`: mature restriction, preflight freeze, Block 5 transport, matched nulls recomputed at 30-plus, BH correction, and report generation.
- Fixed seed: 20261105. Random iterations: 2000 per tested signature, cohort, arm, and unit, recomputed on the 30-plus donors.
- Original full-range benchmark for comparison: `results/followup-11pt-block10-benchmark-results.csv` and `reports/followup-11pt-block10-benchmark.md`.

## Independent verification

BH q-values were recomputed from the raw p-values in results/followup-mature-block10-benchmark-results.csv with a from-scratch Benjamini-Hochberg step-up implementation, separately for all 32 families (4 signatures x 2 cohorts x 2 arms x 2 test types). Per-family results are in results/followup-mature-block10-bh-check.csv (32 rows). Maximum absolute q difference across all 32 families is 1.39e-15; every family matches and every pass count matches. Both-pass counts, per-signature classifications, and the verdict were recomputed from the p-values alone and reproduce the report exactly: INCONCLUSIVE with 2 SURVIVES, 1 LOST_WITH_BATCH, 0 PARTIALLY_UNTESTABLE, and 1 NO_MAJORITY_CATEGORY.

The donor restriction was checked by reloading both cohort inputs independently: Gao heart 109 of 117 unique donors at 30-plus and Yazar blood 932 of 981 unique donors at 30-plus, in both arms, with every per-unit count matching results/followup-mature-block10-donor-counts.csv (65 rows) and every 30-plus minimum age at least 30.0 years.

Decisive unit-level changes behind the verdict flip, all read from the mature and original result CSVs: Sood Gao Adipocyte was the single original adjusted both-pass unit (composite q 0.00398, empirical q 0.00700) and at 30-plus its composite q is 0.0658 while its empirical q is 0.0420, so it fails on the composite test and Sood has zero adjusted passes. Sood Gao unadjusted passes rose from 1 to 3 (Adipocyte, All_LV_Nuclei, vSMC) and Sood Yazar unadjusted gained 1 pass (effector memory CD4-positive T cell, slope minus 0.0651 per decade, a negative-direction pass under the two-sided test), which together with zero adjusted passes in both cohorts makes Sood LOST_WITH_BATCH under the frozen rule. Harries Yazar adjusted passes rose from 4 to 5 with naive thymus-derived CD8-positive T cell added (composite q 0.0328, empirical q 0.0425). Peters Yazar unadjusted passes fell from 17 to 16 with naive B cell dropping out. No signature-cohort became untestable at 30-plus; all 32 family sizes are unchanged from the original.
