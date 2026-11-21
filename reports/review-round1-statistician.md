## 1. Verdict: MAJOR REVISION

The manuscript is not statistically reliable in its current form. Two abstract-level problems are fatal to the claims they support: the 624-test composition family was chosen after both correction outcomes were available and promotes two fragile associations, and the categorical statement that P1 exceeded no matched-random null is contradicted by the committed 2,000-draw analysis. The stored BH arithmetic is mostly correct, but correct arithmetic cannot rescue outcome-dependent family selection, selective sensitivity reporting, or verdict rules that confuse non-significance with absence. The required changes affect the abstract, primary family definitions, transport conclusions, several block verdicts, and the theoretical Discussion. This is not a cosmetic revision.

## 2. FATAL issues

### F1. The blanket matched-random-null claim is false under the committed 2,000-draw P1 analysis

The abstract says that none of the 252-test hits exceeded matched random-signature nulls, with lowest q = 0.052 (`report/plos-one/final-plos-one.tex:37`). The Results repeat that no MSSM cell type had q below 0.05 and no PsychAD cell type survived (`report/plos-one/final-plos-one.tex:235`). Those statements rely on 200 random signatures, for which MSSM oligodendrocytes have empirical p = 0.0049751 and q = 0.0522388 (`results/reviewer-rev2-transport-null-mssm.csv:2`) and PsychAD Micro has q = 0.1007463 (`results/reviewer-phase1-signature-random-null.csv:22`).

The later Block 5 analysis evaluates P0 and P1 with the same primary models and 2,000 random sets (`reports/followup-11pt-block5-p0-transport.md:9-11`). Its committed P1 rows directly contradict the manuscript:

- MSSM oligodendrocytes have composite q = 0.0345577, empirical p = 0.0004997501, empirical q = 0.01049475, and `both_pass=True` (`results/followup-11pt-block5-p0-transport.csv:475`). This is not merely another marginal random-test result. It passes both transport tests.
- PsychAD Micro has empirical q = 0.0404798 (`results/followup-11pt-block5-p0-transport.csv:43`).
- MSSM endothelial cells and vascular leptomeningeal cells have empirical q = 0.0454773 and 0.0314843, respectively (`results/followup-11pt-block5-p0-transport.csv:469` and `:489`).

Thus four core-cohort units exceed the larger matched-random null, and one passes both corrected tests. The manuscript discloses the 2,000-draw sensitivity in Methods (`report/plos-one/final-plos-one.tex:120`) but suppresses its contrary P1 outcomes when stating the transport verdict. Choosing the 200-draw analysis because its Monte Carlo floor leaves oligodendrocytes just above 0.05 is not defensible.

**Concrete fix:** Make the 2,000-draw P1 results the primary matched-null analysis, or report both analyses side by side with equal prominence. Correct the abstract, Results, Discussion, Conclusions, and the scoped transport claim to include MSSM oligodendrocyte transport and the other random-test passes. Explain the Monte Carlo sensitivity rather than declaring that no core-cohort unit exceeded the null.

### F2. The 624-test family is outcome-dependent family-shopping, and the two promoted associations fail the committed robustness checks

The manuscript calls the 624-test family primary (`report/plos-one/final-plos-one.tex:132-134`) and promotes mouse lung vein endothelial cells and human kidney cortex collecting duct cells into the main findings (`report/plos-one/final-plos-one.tex:218` and `:37`). It also says exploratory fly, worm, zebrafish, and naked mole rat arms do not drive the main claims (`report/plos-one/final-plos-one.tex:60`). These positions are incompatible.

The mouse row has p = 0.00271978, core-only q = 0.0689545, and 624-family q = 0.0390405 (`results/rev3-composition-global-bh-624.csv:83`). The human row has p = 0.00247762, core-only q = 0.0667410, and 624-family q = 0.0377082 (`results/rev3-composition-global-bh-624.csv:198`). The raw p-values do not change. The promotions occur only because 193 non-mammalian tests alter the BH ranks. The project record explicitly says both 431-test and 624-test outputs had been produced without selecting one, that choosing 624 rescues these two borderline mammalian rows, and that the family policy was deferred for a later decision (`reports/followup-overnight-summary.md:77-99` and `:122-127`). That is family selection after seeing the correction outcomes, not a pre-stated confirmatory family in any meaningful sense.

The promoted results are also not robust within the committed sensitivity tables:

- Mouse vein endothelial cells have total-cell-adjusted q = 0.0627264 and propeller p = 0.0527594 (`results/reviewer-phase1-composition-total-adjusted.csv:218` and `results/reviewer-phase1-composition-propeller.csv:83`). The primary row contains only 320 cells and 14 present donors (`results/rev3-composition-global-bh-624.csv:83`).
- Human collecting duct cells have total-cell-adjusted q = 0.4231248 and propeller p = 0.0588314 (`results/reviewer-phase1-composition-total-adjusted.csv:597` and `results/reviewer-phase1-composition-propeller.csv:198`). The primary row contains only 45 cells and 11 present donors (`results/rev3-composition-global-bh-624.csv:198`).

The abstract therefore elevates two associations that fail the original mammalian family and both sensitivity checks, while claiming that the added exploratory species do not drive the main findings.

**Concrete fix:** Restore the 431-test mammalian analysis as primary and demote both promoted tests to exploratory sensitivity results. If 624 is retained, provide documentary evidence that its scope was fixed before inspecting the mammalian p-values, remove the claim that exploratory arms do not drive the main claims, and disclose the total-cell and propeller failures next to both associations. Do not call either result established without independent replication.

## 3. MAJOR issues

### M1. Non-significant mature drift is repeatedly misrepresented as proof that drift stops after maturity

The abstract says mouse profiles drifted from young references "but not past maturity" (`report/plos-one/final-plos-one.tex:37`). The Results say the displacement was concentrated in the young-to-mature transition (`report/plos-one/final-plos-one.tex:257`) and that the mature Bayesian refit confirms this (`report/plos-one/final-plos-one.tex:264`). The Conclusions repeat the absence claim (`report/plos-one/final-plos-one.tex:457`).

The 18-month-reference estimate is +0.007229 with 95 percent CI from -0.026995 to +0.041454 and p = 0.6656 (`results/reviewer-rev2-hellinger.csv:135`). The mature Bayesian estimate is +0.0389439 with 94 percent interval from -0.0997789 to +0.177334 and posterior probability positive 0.7068 (`results/rev3-bayesian-vector-mature-only-posterior.csv:183`). That interval includes meaningful continuation in either direction. A significant full-range slope and a non-significant mature slope do not establish a significant difference between those regimes. This is the classic error of treating "significant versus not significant" as evidence of an interaction, and it violates the manuscript's own statement that non-significance is not equivalence (`report/plos-one/final-plos-one.tex:195` and `:422`).

**Concrete fix:** Replace every claim that drift stops or is confirmed to be transition-concentrated with "no post-maturity continuation was detected, but the mature estimate is imprecise and does not establish absence." A trajectory claim requires a direct regime contrast or change-point model with an interval for the difference.

### M2. Block 8 calls two steady, same-direction trends "transition-concentrated" without testing the step difference

The pre-stated rule gives transition classifications precedence when the young-to-middle step passes BH and the middle-to-old step does not, even when an ordinal trend also passes (`report/plos-one/final-plos-one.tex:179` and `reports/followup-11pt-block8-categorical-signature.md:44-50`). The manuscript then claims transition-concentrated effects rather than a universally steady trend (`report/plos-one/final-plos-one.tex:355`). The rule is followed mechanically, but it does not support that scientific classification.

For Yazar central memory CD4 cells, both steps are positive: 0.0000124867 and 0.00000578552. The ordinal trend is highly significant at q = 0.000000172362, and the second step itself has p = 0.0189947 but q = 0.123466 (`results/followup-11pt-block8-categorical-signature.csv:154-157`). For naive CD4 cells, both steps are again positive, 0.00000890608 and 0.00000562485, with trend q = 0.000000649631 and second-step p = 0.0106050 but q = 0.123466 (`results/followup-11pt-block8-categorical-signature.csv:266-269`). In both cases, the spline-versus-linear sensitivity is non-significant, p = 0.314712 and 0.428429 (`results/followup-11pt-block8-categorical-signature.csv:158` and `:270`). Failure of the second step to clear BH is not evidence that it is zero, especially when its nominal interval excludes zero and the ordinal trend passes overwhelmingly.

**Concrete fix:** Reclassify these two cell types as steady trend or inconclusive, not transition-concentrated. Test the contrast between the two step coefficients directly and disclose the ordinal and spline results next to any transition claim. Restrict the narrow Block 8 finding to Gao myeloid cells if it survives such a direct contrast.

### M3. Block 7 does not identify stress-cell dependence because it has no removal control

The stress analysis removes cells above the 90th percentile and declares dependence when a significant high-stress-fraction slope is followed by at least one predicted pass-to-fail gene with attenuated absolute slope (`report/plos-one/final-plos-one.tex:181` and `reports/followup-11pt-block7-survivor-within-type.md:30-38` and `:49-53`). The Results use 77, 83, and 8 such changes to claim stress-cell dependence (`report/plos-one/final-plos-one.tex:353`; `reports/followup-11pt-block7-survivor-within-type.md:132-140`).

There is no matched random-cell removal, no placebo score, no alternative percentile, and no formal test of the slope change. Removing an age-enriched subset mechanically changes the age contrast, the remaining cell mix, and measurement precision. Thousands of gene-level BH decisions will then change without establishing a specific stress mechanism. A rule requiring only one predicted change across roughly 8,552 to 10,777 common genes is particularly weak. The reported net pass-count changes also include other pass-to-fail and fail-to-pass changes, so the analysis is not a controlled intervention with a defined error rate.

**Concrete fix:** Add matched random-cell and non-stress-score removal controls, repeat several thresholds, and test the distribution of slope attenuation against those controls. Define a family-level criterion for reproducible attenuation rather than allowing one gene to prove dependence. Until then, relabel Part A as suggestive stress-score stratification, not PROVES.

### M4. The mouse concordance "detection" is a nominal one-sided result with no cross-species multiplicity control

The abstract and Conclusions say concordance was detected in mouse at p = 0.044 (`report/plos-one/final-plos-one.tex:37` and `:457`). Methods define one permutation test per species and scope but no family spanning those tests (`report/plos-one/final-plos-one.tex:195`). The underlying mouse result is observed concordance 0.185485 against null mean 0.137770, with one-sided `joint_p_greater = 0.0439560` from 1,000 permutations (`results/reviewer-rev2-concordance.csv:2`). The other within-tissue p-values are human 0.700300, fly 0.000999001, worm 0.00299700, zebrafish 0.360639, and naked mole rat 0.70 (`results/reviewer-rev2-concordance.csv:8-12`). A six-test BH check gives mouse q = 0.0879121, so the mouse detection does not survive the most obvious cross-species family.

The human value itself is copied correctly: +0.189108 against null mean +0.207635, p = 0.700300, 1,000 donor-joint permutations (`results/reviewer-rev2-concordance.csv:9`). The problem is not the null construction. It is promoting a marginal, unadjusted mouse p-value alongside the much smaller fly and worm values.

**Concrete fix:** State the one-sided test and its multiplicity scope explicitly. Correct across the predeclared species and scopes, or downgrade mouse concordance to nominal exploratory evidence.

### M5. Block 3 turns a bidirectional dispersion pattern into a dataset-level PROVES while omitting six significant contractions

The dataset rule declares PROVES if at least one cell type spreads (`reports/followup-11pt-block3-profile-cloud.md:31-35`). Two Yazar types meet that rule: memory B cells, difference +0.00890075 and q = 0.0155478 (`results/followup-11pt-block3-profile-cloud.csv:184`), and naive CD8 T cells, difference +0.0219792 and q = 0.000999500 (`results/followup-11pt-block3-profile-cloud.csv:208`). The local numerical claims in the manuscript are correct.

However, six other evaluated types have significant negative batch-residualized differences:

- CD14-positive monocytes: -0.0162773, q = 0.000999500 (`results/followup-11pt-block3-profile-cloud.csv:88`).
- Central memory CD4 cells: -0.0241450, q = 0.000999500 (`results/followup-11pt-block3-profile-cloud.csv:118`).
- Effector memory CD8 cells: -0.0197407, q = 0.000999500 (`results/followup-11pt-block3-profile-cloud.csv:154`).
- Gamma-delta T cells: -0.0163732, q = 0.0113693 (`results/followup-11pt-block3-profile-cloud.csv:166`).
- Naive CD4 cells: -0.0184350, q = 0.000999500 (`results/followup-11pt-block3-profile-cloud.csv:202`).
- Natural killer cells: -0.0286362, q = 0.000999500 (`results/followup-11pt-block3-profile-cloud.csv:214`).

The manuscript reduces all of these to "other Yazar types remained inconclusive" (`report/plos-one/final-plos-one.tex:347`). They are inconclusive only for the one-sided spreading verdict; they are statistically significant contractions. A dataset-level PROVES based on two positives while six significant negatives are suppressed is scientifically misleading.

**Concrete fix:** Replace the single dataset verdict with explicit bidirectional heterogeneity: two spreading types and six contracting types. Report the negative q-values in the main text and do not use the at-least-one-positive rule to summarize the dataset.

### M6. Block 1 retains a PROVES label after its own diagnostics invalidate the measurement interpretation

The pre-stated Block 1 rule requires a positive BH-significant primary slope and agreement in evaluable composition arms (`reports/followup-11pt-block1-residual-dispersion.md:27`). Naive CD4 T cells formally meet it, with primary slope 0.0005 and q = 0.00048 and reweighted slope 0.0007 (`reports/followup-11pt-block1-residual-dispersion.md:71`). The manuscript calls this the sole proving type while acknowledging fragility (`report/plos-one/final-plos-one.tex:343`).

The block's own diagnostics attribute 66.4 percent of the slope to finite-cell effects, show that adding log eligible cell count collapses the slope to 0.0001 with p = 0.1802, show severe residual imbalance after reweighting, and show that reweighting steepens the cell-count decline that generates the artifact (`reports/followup-11pt-block1-residual-dispersion.md:300-302`). Unadjusted, rarefied, gene-fraction, shared-trend, and cell-count-adjusted alternative estimators do not corroborate the signal (`reports/followup-11pt-block1-residual-dispersion.md:301-302`). A verdict rule that ignores these diagnostics is not measuring biological residual dispersion robustly.

**Concrete fix:** Downgrade naive CD4 T cells to INCONCLUSIVE or rename the outcome "nominal rule pass, not robust to finite-cell calibration." Incorporate the finite-cell calibration and cell-count-adjusted model into the verdict rule rather than quarantining them as non-binding diagnostics.

### M7. The Discussion assigns the mouse young-to-mature result to humans

The Discussion says the weak continuation model "fits human, where drift is young-to-mature and concordance is null" (`report/plos-one/final-plos-one.tex:387`). The young-to-mature drift result is explicitly mouse: the 18-month-reference estimate is +0.0072 with p = 0.67 (`report/plos-one/final-plos-one.tex:257`). Human overall displacement is +0.0365 with p = 0.2218 and is described as small and imprecise, not as an established young-to-mature trajectory (`report/plos-one/final-plos-one.tex:262`). This is a direct species mix-up in the theoretical interpretation, not a wording nuance.

**Concrete fix:** Replace "human" with "mouse" for the drift clause, or remove the trajectory argument entirely. Do not use the imprecise human estimate to support a developmental-continuation theory.

## 4. MINOR issues

### m1. Batch attenuation is overinterpreted as proof that unadjusted depth slopes are not age effects

The manuscript says that because unadjusted depth slopes vanish after adjustment, they were "batch confounding, not age effects" (`report/plos-one/final-plos-one.tex:333`). The committed adjusted primary slopes are -0.002107 with p = 0.8071 in heart and -0.0002662 with p = 0.8068 in blood (`results/csv-tests-2026-11-17/depth-vs-age/depth-vs-age-results.csv:3` and `:9`). Adjustment attenuates the associations, but it cannot prove that the unadjusted associations contain no age effect, particularly when pool is strongly associated with age (`report/plos-one/final-plos-one.tex:167`). The related statement that unadjusted public-cohort results are upper bounds on aging signal (`report/plos-one/final-plos-one.tex:397`) is also not generally valid when adjustment can remove true signal or reverse signs.

The separate Gao/Yazar transport wording is appropriately limited: "consistent with a power explanation" is explicitly labeled an interpretation, not a proven result, and GTEx, HRCA, Perez, and the 1 of 13 and 5 of 16 pass counts are disclosed (`report/plos-one/final-plos-one.tex:246`). That sentence should not be strengthened without a formal cross-cohort power analysis.

**Concrete fix:** Use "attenuated after batch adjustment" throughout. Retain the existing power-explanation caveat and add that age-pool collinearity limits separation of technical and biological effects.

### m2. The MSSM depth result lacks the technical covariates used elsewhere in the manuscript

The four MSSM neuronal depth models include only age and sex (`reports/followup-11pt-block6-depth-extension.md:27-33`; `results/followup-11pt-block6-mssm-depth.csv:3`, `:6`, `:9`, and `:12`). The reported primary slopes are accurate: GABAergic -0.0985208, VIP -0.0662748, parvalbumin -0.0598932, and somatostatin -0.0560026 log UMI per decade, with q from 0.0000260080 to 0.00499545 (`results/followup-11pt-block6-mssm-depth.csv:3`, `:6`, `:9`, and `:12`). The mature sensitivities also remain negative (`results/followup-11pt-block6-mssm-depth.csv:2`, `:5`, `:8`, and `:11`). However, the manuscript's primary PsychAD models include brain bank, postmortem interval, and three technical principal components (`report/plos-one/final-plos-one.tex:120`). Without comparable technical adjustment, "direct neuronal depth decline" remains an age association with capture depth, not evidence that the decline is independent of RNA quality or processing.

**Concrete fix:** Add postmortem and available technical covariates to the MSSM depth models, or state explicitly that the result is adjusted only for sex and remains vulnerable to technical confounding.

### m3. Several per-organ "significant" drift claims use unadjusted p-values

The manuscript calls lung significant at 3,000 UMIs with p = 0.030 and bladder significant at full depth with p = 0.048 (`report/plos-one/final-plos-one.tex:259`). Those are the raw values in `results/reviewer-rev2-hellinger.csv:117` and `:52`. Multiple tissues and depths were examined, yet no q-values accompany these local significance claims. Methods say whole-expression vectors are corrected within species and cohort (`report/plos-one/final-plos-one.tex:195`), but the per-organ paragraph does not report such corrections.

**Concrete fix:** Report the relevant adjusted q-values or replace "significant" with "nominally associated" and identify these organ-level results as exploratory.

### m4. Figure 1's caption omits the model that generated all 25 displayed results

The caption says only that the 25 associations pass within-stratum FDR at pseudocount 0.5 (`report/plos-one/final-plos-one.tex:211`). Every figure-data row uses `age_sex_platform_log_total_cells`, not the primary age-plus-sex model (`results/reviewer-phase3-composition-figure-data.csv:1-4`). The distinction matters beside the separate 624-test primary-family claim.

**Concrete fix:** State "log-total-recovered-cell and platform adjusted sensitivity model" in the caption and distinguish these 25 results from the 45 primary-family passes.

### m5. The abstract's 381 interneuron associations should be identified as moderated

The abstract reports 381 replicated associations and six shared genes (`report/plos-one/final-plos-one.tex:37`). The number is the moderated sum 77 + 65 + 239, not the unmoderated sum 76 + 66 + 249 = 391 (`report/plos-one/final-plos-one.tex:280-282`; `reports/followup-11pt-block9-mature-interneuron.md:75-87`). The six-gene intersection is correct in the moderated CSV.

**Concrete fix:** Write "381 associations under empirical-Bayes moderation" so the abstract does not appear inconsistent with the unmoderated counts in the Results.

### Independent BH recomputation checks

I recomputed BH values from raw p-values with an independent script rather than using the stored q-values. "Match" below means exact to CSV floating-point precision, with maximum absolute discrepancy at most 2.22e-16.

- **Composition families:** The 431-row core family gives 15 passes and matches `q_global_bh` exactly (`results/reviewer-phase1-composition-global-bh.csv:1-432`). The 624-row family gives 45 passes and matches `q_global_bh_624` exactly (`results/rev3-composition-global-bh-624.csv:1-625`). The embedded 431-row recalculation also matches `q_global_bh_431` exactly. The two promotion q-values are arithmetically correct; F2 concerns their inferential use, not their computation.
- **Primary signature family:** Combining the 225 M9 rows using `p` with the 27 M10 rows using the fully adjusted `p_covar_adj` gives 252 tests and three passes, matching `q_global_bh` exactly (`results/m9-regressions.csv:1-226` and `results/m10-celltype-regressions.csv:1-28`). The three hits are MSSM L2/3 neurons q = 0.0207907, MSSM oligodendrocytes q = 0.0345577, and PsychAD EN_L3_5_IT_3 q = 0.0207713 (`results/m9-regressions.csv:23-24` and `results/m10-celltype-regressions.csv:6`).
- **Original matched-random families:** The 21 MSSM and 27 PsychAD empirical-p families match exactly, with zero passes in each (`results/reviewer-rev2-transport-null-mssm.csv:1-22` and `results/reviewer-phase1-signature-random-null.csv:1-28`). This verifies only the 200-draw arithmetic; F1 documents the contradictory 2,000-draw results.
- **Block 5 sensitivity families:** The P0 and P1 252-test composite corrections and the empirical corrections for MSSM 21, PsychAD 27, Gao 14, and Yazar 17 all match exactly (`results/followup-11pt-block5-p0-transport.csv:1-635`). These recomputations reproduce three P1 empirical passes in MSSM and one in PsychAD.
- **Block 10 families:** For each of Harries, de Magalhaes, Sood, and Peters, I recomputed all four Gao/Yazar by adjusted/unadjusted slices, giving 16 slices. I independently recomputed composite and empirical BH within every slice, hence 32 correction sequences with sizes 13, 14, or 17 (`results/followup-11pt-block10-benchmark-results.csv:1-247`). Every q-value and every pass count matches exactly. The manuscript's "three of four" verdict follows the frozen majority rule, and the scope limitation to four signatures is explicit (`reports/followup-11pt-block10-benchmark.md:96-99` and `:415-425`).

### Other verdict, power, and inferential-unit checks

- **Donor-as-unit:** No pseudoreplication was found in the reviewed primary claims. The manuscript defines donors, animals, fish, or pools as inferential units (`report/plos-one/final-plos-one.tex:70`), and Blocks 1 to 10 use donor-level models or donor-label permutations, including profile-cloud permutations (`reports/followup-11pt-block3-profile-cloud.md:19-28` and `:274-278`).
- **Powered-null thresholds:** The 0.50 SD rule is applied consistently in Blocks 1, 6, and 7, while Blocks 2 and 3 use their pre-stated split-half technical floors (`reports/followup-11pt-block1-residual-dispersion.md:27`; `reports/followup-11pt-block6-depth-extension.md:35-41`; `reports/followup-11pt-block7-survivor-within-type.md:38-39`; `reports/followup-11pt-block2-gcl.md:28-36`; `reports/followup-11pt-block3-profile-cloud.md:27-35`). The off-program and entropy nulls correctly report MDE80 and failed equivalence rather than claiming absence (`report/plos-one/final-plos-one.tex:313`). M1 is the substantive exception in interpretation.
- **Block 2:** Yazar natural killer cells remain INCONCLUSIVE because the significant slope magnitude 0.0002595 is below the 0.004101 noise floor, not because the q-value is ignored (`report/plos-one/final-plos-one.tex:345`; `reports/followup-11pt-block2-gcl.md:79`). This follows the stated rule.
- **Block 4:** Yazar reverse transport has zero both-test passes in 51 units, but its MDE80 values extend from 0.3736 to 10.3 (`reports/followup-11pt-combined-summary.md:47`). The manuscript correctly calls this DOES NOT PROVE rather than evidence of absence and identifies Gao as INCONCLUSIVE because its composite lacks an up component (`report/plos-one/final-plos-one.tex:349` and `:359`).
- **Block 5:** The strict P0 verdict follows its rule: P0 passes both tests in zero of 59 units where P1 fails both (`report/plos-one/final-plos-one.tex:351`; `reports/followup-11pt-block5-p0-transport.md:137-141`). This does not cure F1's separate P1 reporting failure.
- **Block 6:** Part A Yazar is honestly INCONCLUSIVE because two members of the fixed 18-type candidate family are untestable and ten additional rows fall below donor floors (`reports/followup-11pt-combined-summary.md:81-83`). These are not presented as powered nulls.
- **Block 7 Part B:** The cross-type associations are negative and null in both Gao, -0.03739 with p = 0.3138, and Yazar, -0.00901 with p = 0.2224 (`reports/followup-11pt-block7-survivor-within-type.md:142-149`). DOES NOT PROVE follows the stated localization rule, and the manuscript correctly avoids claiming a complete survivor-bias mechanism (`report/plos-one/final-plos-one.tex:353` and `:399`). This is non-proof, not a powered equivalence result.
- **Block 9:** Mature interneuron replication is honestly INCONCLUSIVE because HBCC is unavailable, and the manuscript does not claim a mature replication null (`report/plos-one/final-plos-one.tex:282`; `reports/followup-11pt-block9-mature-interneuron.md:98-105`).
- **Block 10:** Three survivors among four frozen signatures satisfy DOES NOT PROVE under the written majority rule, and the manuscript correctly refuses to extrapolate to all published signatures (`report/plos-one/final-plos-one.tex:357`; `reports/followup-11pt-block10-benchmark.md:96-99`).
- **Checker scope:** The checker enforces abstract length, citations, figures, section order, and process residue, not statistical correctness (`scripts/check_plos_manuscript.py:13-60`). Its pass provides no validation of the numbers or verdicts.

### Abstract-to-source audit

Every numeric string in the abstract was traced below. The values that are accurate as copies do not escape the inferential objections above.

- **"624-test", "45 passes", and two promoted associations:** The independent recomputation confirms 45 passes in 624 rows. The promoted mouse and human rows are `results/rev3-composition-global-bh-624.csv:83` and `:198`. Their inferential promotion is fatal under F2.
- **"Five of the original 15" mature fraction hits:** The source contains nine mouse DIES rows, one mouse UNTESTABLE row, and five human SURVIVES rows (`results/mature-only-composition-verdicts.csv:2-16`). The study-aware lung attenuation is documented in `results/d3-human-composition-study-sensitivity.csv:16` and `:28`. The abstract statement matches these sources.
- **"139-gene", "three hits in 252 tests", and "lowest q = 0.052":** P1 has 139 rows (`results/m5-program-hierarchy.csv:259-397`). The three 252-test hits match `results/m9-regressions.csv:23-24` and `results/m10-celltype-regressions.csv:6`. The original lowest MSSM random-family q is 0.0522388 (`results/reviewer-rev2-transport-null-mssm.csv:2`). The blanket "none exceeded" statement is false under the later 2,000-draw rows cited in F1.
- **Heart endocardium and five blood cell types:** Gao endocardium has slope +0.0879173, composite q = 0.0150077, and empirical q = 0.00932867 (`results/followup-organ-17-heart-gao-p1-signature.csv:4`). The five Yazar both-test passes are central memory CD4, naive CD4, memory B, naive CD8, and natural killer cells (`reports/followup-11pt-block5-p0-transport.md:108-114`). These counts match.
- **Mouse drift p = 0.000004:** The source is effect +0.102797 with p = 0.00000398288 (`results/h3-cross-species-vector-global.csv:114`). The rounded abstract value matches.
- **Mature p = 0.67 and Bayesian probabilities 0.999 and 0.71:** The mature-reference p is 0.665610 (`results/reviewer-rev2-hellinger.csv:135`). The full-range depth-adjusted posterior probability positive is 0.999 (`results/h3-cross-species-bayesian-posterior.csv:75`), and the mature probability is 0.7068 (`results/rev3-bayesian-vector-mature-only-posterior.csv:183`). The numbers match; the absence interpretation fails under M1.
- **Human displacement p = 0.22:** The source is +0.0364830 with p = 0.221842 (`results/h3-cross-species-vector-global.csv:110`). The abstract correctly calls it small and imprecise.
- **"381 associations" and "six genes":** The moderated counts are 77 VIP, 65 parvalbumin, and 239 somatostatin, totaling 381, and the intersection is FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16 (`results/reviewer-rev2-interneuron-genes-moderated.csv:1-769`; `reports/followup-11pt-block9-mature-interneuron.md:83-96`). These match.
- **Concordance numbers:** Mouse p = 0.0439560, fly p = 0.000999001, worm p = 0.00299700, human +0.189108 with p = 0.700300, zebrafish p = 0.360639, and naked mole rat exact p = 0.70 all match `results/reviewer-rev2-concordance.csv:2` and `:8-12`. Human, mouse, fly, worm, and zebrafish used 1,000 permutations; naked mole rat used all ten assignments (`report/plos-one/final-plos-one.tex:147`). The human null is donor-joint and is reported with the correct null direction. Multiplicity remains a problem under M4.
- **Memory B and naive CD8 cloud spreading:** The two positive differences and q-values match `results/followup-11pt-block3-profile-cloud.csv:184` and `:208`. The six omitted contractions are documented under M5.
- **Four MSSM neuronal units:** The four primary slopes range from -0.0560026 to -0.0985208 log UMI per decade and all four primary q-values pass (`results/followup-11pt-block6-mssm-depth.csv:3`, `:6`, `:9`, and `:12`). The numbers match; the covariate caveat is m2.
- **P0 and reverse transport:** The strict verdict statements match `reports/followup-11pt-block5-p0-transport.md:137-141` and `reports/followup-11pt-block4-reverse-transport.md:36-45`. These are scoped non-proof verdicts, not evidence of universal non-transport.
- **Three of four published signatures:** Harries, Sood, and Peters survive, while de Magalhaes has no unadjusted both-test pass to lose (`reports/followup-11pt-block10-benchmark.md:127-134` and `:415-425`). The abstract count and its scope are correct.
