# Removing zebrafish and naked mole rat from the Agescape paper

Date: 2026-11-25. Branch: `remove-zebrafish-nmr`, from master 43f9734. Ordered by Reza 2026-11-25 01:08 PT.

## What was removed

1. All zebrafish and naked mole rat results, nulls, caveats, dataset rows, and references were removed from `report/plos-one/final-plos-one.tex`, except one Methods paragraph (see item 4).
2. Dataset table: deleted the naked mole rat skin row, the zebrafish row, and the naked mole rat spleen row.
3. Results: deleted the zebrafish broad/fine vector paragraph, the naked mole rat skin vector sentence, the naked mole rat spleen T-cell sentence, and the naked mole rat mortality sentence with its citation.
4. Methods: deleted the zebrafish and naked mole rat sentences from the drift, concordance, and Bayesian Methods subsections, the 31-year naked mole rat lifespan from the AnAge scaling list, and the generic "fish" from the inferential-unit definitions. The single permitted mention stays: the Primary composition family paragraph states both species were examined and excluded, with reasons (zebrafish: 4 mature fish at 12, 20, and 22 months, below the 6-fish floor, testis only; naked mole rat: no mature window, 3 animals at 1.1 years vs 2 at 11.3 years, no intermediate ages).
5. Limitations: deleted the naked mole rat arm sentence, the zebrafish testis/germline paragraph, the naked mole rat spleen colony-confound clause, the LACA (Longevity and Aging Cell Atlas) clauses in Results and Limitations, and the clause saying the exclusions qualify the primary fraction count. The three-most-constrain list now reads: pooling species on different mature clocks, the all-male old mouse group, and the missing progeroid control.
6. Appendix: deleted the naked mole rat spleen T-cell row from the follow-up table and the exclusion restatement from Appendix C.
7. Data Availability: removed GEO GSE214390 (naked mole rat spleen) from the follow-up accession list.
8. References: removed entries 6 (Sposato, zebrafish), 7 (Savina, naked mole rat), 8 (Hilton, naked mole rat), and 18 (Ruby, naked mole rat mortality). All later entries and all in-text citations were renumbered by script: entries 9 to 17 shift down by 3, entries 19 to 81 shift down by 4. Verified: cited set equals listed set, 1 to 77, no gaps.
9. The 624-test full-range family became the 605-test family everywhere (see below). The substring 431 now appears nowhere in the tex. The substring 624 remains only where it is not the family: the CELLxGENE collection UUID containing a624, the pooled q value 0.624 for marrow precursor B cells, and the interval bound 2.624.

## The 605-test family recompute

1. Method: Benjamini-Hochberg correction over the p-values in `results/rev3-composition-global-bh-624.csv`, dropping the 15 zebrafish and 4 naked mole rat rows. Script: `/tmp/recompute_605.py` (run on Beast in the worktree).
2. Verification: the same BH code reproduces the committed `q_global_bh_624` column to a maximum absolute difference of 2.2e-16 with 45 passes, and reproduces the committed `q_global_bh_431` column on the 431 core rows to the same tolerance. The new file `results/rev3-composition-global-bh-605.csv` (605 rows, new `q_global_bh_605` column) has 45 passes, the same 45 rows, zero gained, zero lost. Passes by species: fruit fly 26, mouse 11, human 6, worm 2. The 624 file is untouched as the historical record.
3. Core mammalian family: 431 tests become 427 tests (136 mouse + 291 human). The 427-family BH was recomputed the same way and has the same 15 passes as the 431-family (all 15 are mouse/human rows; naked mole rat contributed 0 passes). The 32 retained / 12 not retained / 1 untestable (mouse bladder cell) reconciliation is unchanged.
4. Values used in the tex:
   - Line 262 rows (macrophage, mesenchymal stem cell): 605-family q 0.847 and 0.848; 427-family q 0.865 for both (exact 0.8647695733928269). Their 624-family q values were 0.855 and 0.856; their 431-family q was 0.870.
   - Promoted tests: mouse lung vein endothelial 605-family q 0.0379 (was 0.0390); human kidney cortex collecting duct 605-family q 0.0366 (was 0.0377).
   - Dopaminergic PAM neuron row: 605-family q 0.0482 (was 0.0497). Spermatheca row: 605-family q 0.00169 (was 0.00174).
   - Concordance: full-range six-species family mouse q 0.0879 becomes four-species family mouse q 0.0586 (BH over mouse p 0.043956, human p 0.700300, fly p 0.000999, worm p 0.002997; recomputed independently and cross-checked against `results/followup-mature-concordance-bh-check.csv`). The verdict is unchanged: still does not survive correction. Fly and worm full-range q values (0.00400, 0.00599) are not cited in the tex.

## Not cleanly removed (none)

1. Everything listed above was removed cleanly. No clause was dropped for lack of a recompute.
2. The LACA atlas clauses were deleted outright (they named a naked mole rat atlas). The record that the LACA atlas could not be verified lives in the repo reports, not in the paper.
3. The 624 file, the hierarchy Level 0/Level 2 CSVs, and the concordance CSVs still contain zebrafish and naked mole rat rows as historical records. They are not cited in the tex. No tex change was needed for the four main figures (all mouse/human only).

## Verification

1. `scripts/check_plos_manuscript.py`: PASS (see summary file).
2. Tectonic compile: PDF page count in the summary file.
3. Grep: zebrafish / naked mole rat appear only in the permitted Methods paragraph; 431 appears nowhere; 624 appears only as the UUID, 0.624, and 2.624; cited set equals listed set (77 references).
