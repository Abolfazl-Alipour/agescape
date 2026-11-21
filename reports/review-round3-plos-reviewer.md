# Hostile PLOS ONE review of Agescape

## Verdict

**MAJOR REVISION.** The manuscript contains an unusually large amount of donor-level analysis, and much of the arithmetic checks out, but the current presentation is not technically reliable enough for publication. A central categorical statement about P1 failing every matched-random test is contradicted by the repository's newer, larger analysis of the same primary MSSM model. The two promoted composition findings become significant only after exploratory tests are added to the correction family, despite the promise that those tests do not drive the main claims. Several Discussion statements go beyond the measured evidence, and the advertised reproduction route does not regenerate this manuscript. The 292-word abstract passes the length check, the four included figures show the basic quantities named in their captions, and the ten follow-up verdict labels mostly match their source rows, but those strengths do not repair selective reporting. The single change that would most improve the paper is to freeze one primary analysis manifest and regenerate every claim from it, explicitly including the 2,000-draw P1 null and the core-only composition results.

## FATAL issues

### F1. The categorical P1 matched-null claim is contradicted by the repository's own larger primary analysis

The Abstract states that none of the three global-screen associations exceeded matched random-signature nulls, with lowest q = 0.052 (`report/plos-one/final-plos-one.tex:37`). The Results repeat that no MSSM cell type had q below 0.05 (`report/plos-one/final-plos-one.tex:235`). That statement rests on 200 random signatures per cell type (`report/plos-one/final-plos-one.tex:120`). However, the newer Block 5 comparison evaluates the same `MSSM_DLPFC`, `oligodendrocyte`, `P1_Clean`, `primary`, `age_plus_sex` model with 2,000 random sets and reports empirical p = 0.00049975, empirical q = 0.01049475, and `both_pass=True` (`results/followup-11pt-block5-p0-transport.csv:475`). The Block 5 report itself marks this P1 row as both-pass (`reports/followup-11pt-block5-p0-transport.md:103`). This is not an unadjusted sensitivity or a different biological model. It is the same primary comparison with more Monte Carlo draws. The manuscript uses the newer Block 5 analysis while suppressing its contrary P1 result.

**Required fix:** Make the 2,000-draw analysis primary, or report both resolutions prominently wherever the null claim appears. Replace the categorical "none exceeded" claim with the actual result: MSSM oligodendrocytes pass the 2,000-draw matched-null test. Reassess the Abstract, Results, Discussion, Conclusions, and any claim that P1 transport failed in all core brain, kidney, and lung cohorts. Do not hide behind the lower-resolution 200-draw run.

## MAJOR issues

### M1. The 624-test family manufactures two headline promotions from exploratory arms that supposedly do not drive the paper

The Introduction says exploratory fly, worm, zebrafish, and naked mole rat analyses "do not drive the title or the main claims" (`report/plos-one/final-plos-one.tex:60`). Yet the primary family adds 193 fly, worm, and zebrafish tests to 431 core tests (`report/plos-one/final-plos-one.tex:134`), and the two promoted core associations are significant only in that enlarged family (`report/plos-one/final-plos-one.tex:218`). Mouse lung vein endothelial cells have core-only q = 0.0689545 but 624-family q = 0.0390405 (`results/rev3-composition-global-bh-624.csv:83`). Human kidney collecting duct cells have core-only q = 0.0667410 but 624-family q = 0.0377082 (`results/rev3-composition-global-bh-624.csv:198`). Thus exploratory tests directly determine two Abstract-level findings (`report/plos-one/final-plos-one.tex:37`). This is not independent replication or additional biological support. It is a change in BH ranks caused by expanding the family.

The underlying rows also show fragility that the headline omits: the mouse row has 320 total cells and `n_present=14` (`results/rev3-composition-global-bh-624.csv:83`), while the human row has only 45 total cells and `n_present=11` (`results/rev3-composition-global-bh-624.csv:198`). The manuscript gives no independent replication for either promoted result.

**Fix:** Lead with the core-only 431-test analysis and demote both promoted tests to exploratory sensitivity results. If the 624-family is retained as primary, justify its scope and timing, remove the claim that exploratory arms do not drive main findings, and show donor-level cell counts, leverage diagnostics, pseudocount and total-cell sensitivities, and independent replication for both tests. Add this dependence to Limitations rather than celebrating it as promotion.

### M2. "Maturation-concentrated" is an overclaim when both Yazar cell types also have continuing old-versus-mid increases and strong ordinal trends

The Results say categorical models prove transition-concentrated effects and label Yazar central memory CD4 and naive CD4 cells maturation-concentrated (`report/plos-one/final-plos-one.tex:355`). The classification rule gives transition labels precedence whenever the young-to-mid step passes correction and old-to-mid does not, even if the ordinal trend also passes (`reports/followup-11pt-block8-categorical-signature.md:46`). That rule does not establish concentration in maturation. For central memory CD4 cells, old-to-mid is positive at 0.000005786 with nominal p = 0.018995, and the ordinal trend has q = 1.724e-07 (`results/followup-11pt-block8-categorical-signature.csv, rows 154 to 157`). For naive CD4 cells, old-to-mid is positive at 0.000005625 with nominal p = 0.010605, and trend q = 6.496e-07 (`results/followup-11pt-block8-categorical-signature.csv, rows 266 to 269`). A nonsignificant corrected old-to-mid contrast is not evidence that the later effect is absent, especially when its point estimate is positive and the continuous trend is decisive.

**Fix:** Rename these outcomes as significant young-to-mid effects with unresolved continuing increases, not maturation-concentrated effects. Report all three contrasts together and stop using the precedence rule to erase the ordinal evidence. If concentration is a scientific claim, test the difference between early and late slopes and require evidence against a meaningful late effect, not merely q >= 0.05.

### M3. The Discussion assigns the mouse young-to-mature drift result to humans

The Discussion says the weak continuation model "fits human, where drift is young-to-mature and concordance is null" (`report/plos-one/final-plos-one.tex:387`). The young-to-mature result is explicitly mouse: the 18-month reference gives effect +0.0072 with p = 0.67 (`report/plos-one/final-plos-one.tex:257`). The human overall estimate is +0.0365 with p = 0.2218 and is described as small and imprecise, not as established young-to-mature drift (`report/plos-one/final-plos-one.tex:262`). This is a direct species mix-up in the theoretical interpretation.

**Fix:** Rewrite the sentence so that young-to-mature drift refers to mouse, while human drift remains imprecise and human concordance remains null. Remove any implication that the human data establish a developmental-versus-aging transition.

### M4. The interneuron Discussion invents temporal ordering and disease progression that the analysis cannot measure

The Discussion calls VIP, parvalbumin, and somatostatin interneurons "the subtypes that change most and earliest" and presents their expression changes as a plausible preclinical footprint of an ordered vulnerability sequence (`report/plos-one/final-plos-one.tex:391`). The cited result file contains genome-wide correlations and FDR overlap counts, not temporal onset, longitudinal ordering, disease conversion, or a comparison establishing "most and earliest" (`results/reviewer-rev2-replication-summary.csv, rows 2 to 5`). The manuscript itself concedes that other cell types show age-associated expression and that SEA-AD did not replicate the program (`report/plos-one/final-plos-one.tex:282`). Cross-sectional age slopes in two brain banks cannot establish a preclinical disease sequence.

**Fix:** Delete "most and earliest" and the claimed ordered-vulnerability footprint. Retain the narrower finding of replicated age-associated expression in three subtypes, and label any Alzheimer comparison as speculative motivation for a longitudinal study rather than interpretation of demonstrated ordering.

### M5. The manuscript does not define the decision rules behind its repeated "proves", "does not prove", and "inconclusive" language

The Methods say verdicts were pre-stated but do not provide the dataset-level rules (`report/plos-one/final-plos-one.tex, lines 169 to 171`). The multiplicity section lists correction families, not verdict thresholds (`report/plos-one/final-plos-one.tex, lines 193 to 195`). Results then assert that Part B proves a depth decline, profile-cloud dispersion proves spreading, and other tests do not prove effects (`report/plos-one/final-plos-one.tex, lines 335 to 347`). The actual rules are materially different across blocks: Block 1 uses cell-type proportions and powered-null thresholds (`reports/followup-11pt-block1-residual-dispersion.md:27`), whereas Block 3 allows a dataset-level PROVES verdict from a single cell type (`reports/followup-11pt-block3-profile-cloud.md:35`). A reader cannot reproduce the manuscript's inferential labels from the Methods alone, and "proof" is much stronger than the underlying observational evidence warrants.

**Fix:** Add a table giving every block's unit-level and dataset-level verdict rule, MDE80 threshold, technical-floor requirement, treatment of untestable rows, and handling of conflicting arms. Replace "prove" with specific statistical statements such as "met the prespecified spreading criterion". Link each rule to its exact report and script in Data Availability.

### M6. The advertised reproduction route does not reproduce the submitted manuscript or its follow-up analyses

Data Availability points readers to the repository (`report/plos-one/final-plos-one.tex:461`). The README presents `python scripts/run_all.py` as the core pipeline (`README.md, lines 44 to 48`), but the runner's documented stages contain no M9, M10, revision, or 11-point follow-up stages (`scripts/run_all.py, lines 9 to 20`), and its implemented later stages stop at earlier H2, H3, and fibroblast analyses (`scripts/run_all.py, lines 117 to 163`). The README lists separate M9 and PsychAD commands (`README.md, lines 60 to 79`) but no commands for the ten follow-up blocks. Therefore the advertised route cannot regenerate the paper's decisive new claims.

Source identification is also incomplete in the manuscript. The Methods name the Census release but not the 11 dataset IDs (`report/plos-one/final-plos-one.tex:70`), which exist only in code (`scripts/run_m9_pipeline.py, lines 32 to 44`). The follow-up table gives the retina atlas accession merely as "Integration" and Tabula Sapiens merely as "CELLxGENE" (`report/plos-one/final-plos-one.tex, lines 498 to 509`). Finally, no immutable manuscript commit or archived release is specified in Data Availability (`report/plos-one/final-plos-one.tex:461`).

**Fix:** Supply one workflow manifest with ordered commands, inputs, outputs, seeds, and environment specifications for every number in the submitted paper. List all Census dataset IDs, exact retina and Tabula Sapiens source versions, and an immutable commit plus archived release. Verify the workflow from a clean checkout rather than citing a legacy runner.

### M7. The highlighted MSSM neuronal depth decline lacks the technical adjustment needed for a credible interpretation

The Abstract and Results emphasize depth declines in all four MSSM neuronal units (`report/plos-one/final-plos-one.tex:37,335`). The MSSM model uses only age and sex (`report/plos-one/final-plos-one.tex:187`; `reports/followup-11pt-block6-depth-extension.md, lines 29 to 31`). It does not include postmortem interval or the technical principal components used in the primary PsychAD model (`report/plos-one/final-plos-one.tex:120`). The four primary MSSM slopes are indeed negative, from -0.0560 to -0.0985 log UMI per decade (`results/followup-11pt-block6-mssm-depth.csv:3,6,9,12`), but this establishes an association with captured library size, not its origin. Limitations discuss depth generically (`report/plos-one/final-plos-one.tex:439`) without acknowledging that this highlighted MSSM analysis lacks available postmortem and technical adjustment and conditions on cells passing a 1,500-count eligibility filter (`reports/followup-11pt-block6-depth-extension.md:29`).

**Fix:** Add PMI and technical covariates if available, and report sensitivity to the eligibility threshold and cell composition. State explicitly in Results and Limitations that this is a technical measurement association with unresolved postmortem and selection explanations, not evidence of biological change in transcriptional output.

## MINOR issues

### m1. The Abstract is within the limit but overpacked and hides important analysis scopes

The Abstract is 292 words, so length is compliant, but it compresses composition, transport, drift, concordance, interneurons, and all ten follow-ups into one paragraph (`report/plos-one/final-plos-one.tex:37`). It does not say that the original matched-null claim concerns 200 draws in MSSM and PsychAD (`report/plos-one/final-plos-one.tex:120`). Its "381 associations" is the moderated total, 77 + 65 + 239, whereas the preceding Results paragraph gives the unmoderated total of 391 (`report/plos-one/final-plos-one.tex, lines 280 to 282`). These are not interchangeable analyses.

**Fix:** Shorten the Abstract to the primary questions and strongest supported results. Identify the null cohorts and draw count, and write "381 associations under empirical-Bayes moderation" if that number remains after F1 is corrected.

### m2. Figure 1 shows a sensitivity model but its caption does not identify it as such

The caption describes 25 within-stratum FDR associations (`report/plos-one/final-plos-one.tex:211`), and the image does show 25 slopes with intervals. However, every source row uses `age_sex_platform_log_total_cells`, not the primary age-plus-sex model (`results/reviewer-phase3-composition-figure-data.csv, rows 1 to 26`). This distinction matters because the adjacent text emphasizes the separate 624-test primary family (`report/plos-one/final-plos-one.tex:218`). The two promoted primary tests are not the subject of this figure.

**Fix:** State "log-total-recovered-cell and platform adjusted sensitivity model" in the caption and explicitly distinguish these 25 results from the 45 primary-family passes. Consider adding the two promoted tests in a separate panel with their core-only and full-family q values.

### m3. Appendix A falsely says the main text leads with mature-reference results

Appendix A claims, "The main text leads with mature-reference results" (`report/plos-one/final-plos-one.tex:479`). The vector Results instead lead with young-reference effect +0.1027 and extensive young-reference sensitivities (`report/plos-one/final-plos-one.tex, lines 248 to 255`), introducing the 18-month result only afterward (`report/plos-one/final-plos-one.tex:257`).

**Fix:** Replace the sentence with an accurate statement that both references are reported, or reorder the Results if mature-reference results are genuinely intended to lead.

### m4. Common abbreviations are not defined at first use, contrary to generalist readability

UMI is used without expansion in the cohort definition and mouse processing (`report/plos-one/final-plos-one.tex:70,115`). MSSM and HBCC appear as table labels without expansion (`report/plos-one/final-plos-one.tex, lines 99 to 100`). FDR appears in the Figure 1 caption (`report/plos-one/final-plos-one.tex:211`) before "false discovery rate" is spelled out much later (`report/plos-one/final-plos-one.tex:363`). By contrast, cell fractions and off-program expression are correctly scoped as recovered-cell shares and as distinct from Shannon entropy (`report/plos-one/final-plos-one.tex:125,151`).

**Fix:** Expand UMI, MSSM, HBCC, FDR, and other cohort acronyms at first use. Retain the correct fraction and off-program definitions.

### m5. The controlled-access statement is broader than the analyses actually performed

The follow-up Methods state, "controlled-access data were not analyzed" (`report/plos-one/final-plos-one.tex:167`). Yet PsychAD controlled-access donor-level data were analyzed (`report/plos-one/final-plos-one.tex:70`), and Data Availability correctly says only aggregate outputs were shared (`report/plos-one/final-plos-one.tex:461`). The intended scope is presumably new follow-up cohorts, but the blanket sentence is misleading.

**Fix:** Write that no new controlled-access follow-up cohorts were analyzed and that the previously authorized PsychAD analysis remained controlled-access with aggregate outputs only.

### m6. The Census citation does not document the Census resource actually used

The Methods cite reference 2 for CZ CELLxGENE Census release 2025-11-08 (`report/plos-one/final-plos-one.tex:70`), but reference 2 is the Tabula Sapiens atlas paper (`report/plos-one/final-plos-one.tex:523`), not documentation for the Census software or release. This matters because exact release and dataset identity determine the downloaded observations.

**Fix:** Add the appropriate CELLxGENE Census citation and persistent release URL, while retaining Tabula Sapiens only where that atlas itself is discussed.

## Audit boundary

All 61 listed references are cited in the body, and 13 source spot checks found matching bibliographic records and relevant claims. The reality and exact claim matching of the other 48 references are **UNVERIFIED**, not alleged to be defective. The zebrafish limitation is present and correctly restricts inference to testis and 97.23 percent germline composition (`report/plos-one/final-plos-one.tex:134,427`). The reported human concordance value +0.1891 with p = 0.7003 is correct (`results/reviewer-rev2-joint-concordance-human.csv:2`). Block 4 has zero of 51 Yazar both-test passes and no estimable Gao arm, Block 5 has zero P0 rescues among 59 strict P1 failures, Block 6 has 13 Gao fits with minimum negative q = 0.2474 and 16 Yazar fits with 12 untestable rows, and Blocks 7, 9, and 10 otherwise match their cited CSV summaries. These verified points do not remove F1 or M2.

FATAL 1, MAJOR 7, MINOR 6.
