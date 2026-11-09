# PLOS ONE manuscript

Source of truth: `final-plos-one.tex`.

Compiled deliverable: `final-plos-one.pdf`, built with Tectonic 0.15.0.

The manuscript uses the middle-claim conclusion: no uniform transcriptomic program was detected in the surveyed datasets, without settling upstream causes or excluding shared mechanisms with cell-type-specific consequences. It contains four data figures and two cohort tables. No AI-generated figure or six-species meta-analysis is included.

Affiliation remains blank. Funding, competing interests, and ethics remain pending author confirmation. The manuscript has not been submitted.

To rebuild from this directory:

```bash
tectonic final-plos-one.tex
```

The older `report/manuscript-aging-mechanisms.md` is not part of this manuscript build.

## Rev2 stage 3 (final refresh)

Two commits landed after stage 2 and made the concordance section stale; this stage fixes that plus two small cleanups. It was rebased on top of the parallel revision (merge cdf4a4d) that rewrote the Introduction, Methods cohorts, interneuron Results, and Discussion; the stage 3 changes below were re-applied onto that rewritten text. The headline middle claim is unchanged.

- Joint donor-level concordance nulls REDONE for fly, worm, and zebrafish (`scripts/run_rev2_joint_concordance.py`, commit ed1577a; 1,000 permutations each, one age assignment per donor or pool shared across cell types, permuted within cohort, sex, and platform strata). Fruit fly: observed plus 0.2191 against a joint null of plus 0.0384 (95 percent interval plus 0.0126 to plus 0.0874), p = 0.001, 99 units, 72 donors, 2 strata; head plus 0.2581 vs null plus 0.0396 and body plus 0.1336 vs null plus 0.0360, both p = 0.001. C. elegans: observed plus 0.3718 against a joint null of plus 0.1841 (95 percent interval plus 0.0723 to plus 0.3271), p = 0.003, 13 units, 12 donors, single whole-animal stratum (small sample, wide null interval). Zebrafish: observed plus 0.3415 against a joint null of plus 0.3254 (95 percent interval plus 0.2400 to plus 0.4365), p = 0.361, 14 units, 6 donors, NOT significant; the null itself is large, so the observed value barely exceeds it. Mouse is unchanged (within-tissue plus 0.1868 vs null plus 0.1378, p = 0.0430; cross-tissue not significant). Human was NOT redone (no joint human run exists) and stays reference-only. The concordance Results paragraph, the concordance Methods paragraph (which had wrongly implied human used 1,000 permutations; corrected), the exploratory-arms Discussion paragraphs, the Discussion coordination sentence, and the Conclusions now report the redone numbers. Within-species concordance in fly and worm is shared direction within those species, not evidence of a cross-species uniform program, so the middle claim stands.
- Figure 3 REBUILT from the new joint-null data (`results/reviewer-rev2-concordance-figure-data.csv`, same design and filenames as the stage 1 build): joint null squares with 95 percent intervals now shown for mouse, naked mole rat, fly, worm, and zebrafish; the human row stays an open circle marked reference only, with no null square. The redone observed values use the same detection filter as the mouse run, so they differ slightly from the earlier Phase 1 reference values (fly 0.2191 vs 0.1936, worm 0.3718 vs 0.3537, zebrafish 0.3415 vs 0.3393). Note: the parallel revision added `scripts/build_rev2_manuscript_figures.py`, whose concordance builder reads `results/reviewer-rev2-joint-concordance-human.csv`, which does not exist, so that builder currently crashes; Figure 3 was rebuilt with the stage 1 design instead, which handles the human reference-only row.
- Transport null files DEDUPLICATED. `results/reviewer-rev2-transport-null.csv` (stage 2) and `results/reviewer-rev2-transport-null-mssm.csv` (commit 38e7740) were the same null run, differing only in last-digit float rounding on four rows (L6, L6b, T cell, vascular leptomeningeal); all headline numbers identical. The un-suffixed file was deleted with `git rm`; the `-mssm.csv` file is kept because the merged script `scripts/run_rev2_transport_null_h5ad.py` writes that filename (line 294), so the committed file now matches what the script reproduces. The manuscript never cited the deleted filename. Note: the older superseded script `scripts/run_rev2_transport_null.py` still writes the deleted filename; it was not touched in this stage.
- Leave-one-donor-out wording PINNED to the unmoderated set, and the sentence RESTORED: the parallel interneuron rewrite had dropped the leave-one-donor-out stability sentence entirely. `results/reviewer-rev2-lodo-summary.csv` covers 391 genes x 2 cohorts = 782 refits, and those 391 entries (76 VIP, 66 parvalbumin, 249 somatostatin gene-by-cell-type entries) are the unmoderated both-FDR set, not the moderated 381 (77, 65, 239). Verified from the file: 0 rows with a sign flip, median largest single-donor shift 0.48 standard errors, 3 cases above one standard error, 4 cases losing nominal significance. The manuscript now says "leave-one-donor-out refits of the 391 unmoderated replicated genes".
- Abstract: 238 words after the parallel rewrite (this stage did not change the abstract), no statistical results. Affiliation stays blank; Funding, Competing Interests, and Ethics stay marked pending author confirmation.

See `~/workspace/agescape-final/plos-one/REV2_STAGE3_NOTES.md` for the full change list, verified numbers, figure decision, and anything not done.
