# Other-species matched per-gene replication and mature concordance

Analysis ordered by Reza on 2026-11-24; executed with `kimi-code/k3.1-pro` on the `other-species-pergene` branch created from master `8266f51`. No manuscript or `.tex` file was changed.

## Headline

**No matched cross-species per-gene test could run.** Under the existing mouse-human test's exact normalized same-tissue cell-type rule, every candidate species has **0 matched pairs against mouse and 0 against human**. A second audit that ignored tissue also found 0 name matches in every source cell type. Thus there are no sign-agreement, Spearman, Pearson, or correlation p-values to estimate. This is an applicability blockage, not a negative replication result.

The separate **within-species mature-to-old concordance refits did run** for fly and worm. Fly has joint permutation p = 0.000999 and worm has p = 0.011988, each with 1,000 joint pool-level permutations. These are concordances among cell types within one species, not matched cross-species per-gene replication tests.

## Fixed reference design and matching rule

The reference implementation is `scripts/run_mature_matched_p1m.py`; its committed CSV and report were read before this audit.

- Mouse inputs are the existing 36 P1m-eligible units in lung, kidney, and marrow, with the committed unit-specific 18m-to-30m effects and gene list. They were not re-derived.
- Human inputs are the existing matched-test roster: kidney, lung, and blood cohorts, including the existing donor caches. The present-source roster contains 249 cohort-cell-type rows and 177 unique cell-type labels.
- A match requires identical cell-type names after lowercasing and removing punctuation and whitespace, within the same tissue. As in the reference test, mouse marrow targets human blood. No synonyms, fuzzy matches, broad cross-tissue classes, or new ontology mappings were introduced.
- A candidate unit is slope-fittable in the unfloored arm with at least 2 inferential units and at least 2 ages. The floored arm requires at least 3 inferential units and at least 50 total cells. Pools, fish, or animals are the inferential units, never cells.
- Shared or similar cell-type names would not prove evolutionary homology. For these atlases even that conservative name criterion produced no pairs. Ortholog mapping cannot repair the absence of a matched cell-type pair, so no expression slopes were fabricated and no per-pair result CSV was written.

The reproducible audit is `scripts/audit_mature_other_species_pergene.py`. Its 232 cell-type-scope rows are in `results/followup-mature-other-species-pergene-cell-types.csv`; the eight source-scope summaries are in `results/followup-mature-other-species-pergene-sufficiency.csv`.

## Per-species sufficiency and decisions

“Fittable” and “floored” below are counts of source cell types under the matched per-gene rules above, not the stricter six-unit concordance rule. “Matches” are exact normalized same-tissue matches to the fixed reference rosters.

| Species / source | Analysis scope | Inferential units and ages | Cell types | Fittable | Floored | Mouse matches | Human matches | Decision |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| Fly AFCA head | 30+ days | 36 pools; 30/50/70 days | 90 | 90 | 88 | 0 | 0 | Neither test can run |
| Fly AFCA body | 30+ days | 36 pools; 30/50/70 days | 72 | 72 | 63 | 0 | 0 | Neither test can run |
| Worm Gao | 6+ days | 9 pools; 6/12/14 days | 16 | 16 | 15 | 0 | 0 | Neither test can run |
| Worm Gao, full-range check | 1+ days | 12 pools; 1/6/12/14 days | 16 | 16 | 16 | 0 | 0 | Adding day 1 creates no match |
| Zebrafish fine labels | 12+ months | 4 fish; 12/20/22 months | 15 | 15 | 15 | 0 | 0 | Neither test can run |
| Zebrafish fine labels, full-range check | 5–22 months | 6 fish; 5/12/20/22 months | 16 | 15 | 15 | 0 | 0 | Adding 5-month fish creates no match |
| NMR skin | All available ages; no defensible mature window | 5 animals; three at 1.1 and two at 11.3 years | 2 | 2 | 2 | 0 | 0 | Neither test can run |
| NMR spleen | All available ages; no defensible mature window | 7 animals total; three at 2 and four at 24 years; 5–7 represented per subset | 5 | 5 | 5 | 0 | 0 | Neither test can run; raw inputs also absent |

### Fly

The resolved AFCA head and body H5AD files exist and contain only the established AFCA 30/50/70-day design after the source filters. There are 72 pools and 162 annotated organ-cell-type combinations. Power is not the blocker: all 162 are fittable under the two-pool rule and 151 pass the three-pool/50-cell floor. Nevertheless, none of the 162 labels matches any mouse or human reference label, even after discarding the tissue requirement. Fly-vs-mouse and fly-vs-human therefore both have **0 matched pairs and were not run**.

### Worm

The resolved Gao integrated H5AD exists. Restricting normal N2 singlets to 6+ days leaves 58,236 annotated cells in 9 pools at days 6, 12, and 14. All 16 labels are fittable and 15 pass the matched-test floor. None matches a reference label, with or without the tissue restriction. Worm-vs-mouse and worm-vs-human therefore both have **0 matched pairs and were not run**. This conclusion does not rely on the older feasibility audit.

### Zebrafish

The fine-label source contains 6 fish at 5, 12, 20, and 22 months; the separate 27-month broad-label atlas was not substituted for comparable fine labels. At the established 12-month mature cutoff there are 4 fish, 15 represented fine labels, and all 15 pass the matched per-gene floor. Thus the earlier six-fish composition floor is not the reason for exclusion here. Instead, none of the testis labels matches a mouse or human reference label, even ignoring tissue. Zebrafish-vs-mouse and zebrafish-vs-human both have **0 matched pairs and were not run**. The germline skew remains a further interpretation limit: 35,511 of 36,523 broad-known cells in the 5–22-month source are germ cells (97.229143%).

### Naked mole-rat

For GSE152886 skin, the raw count file is present on Beast at `/home/reza/projects/agescape/data/raw/GSE152886_nmr_skin_raw.csv.gz`, although it is not under `agescape-revision-data`. The committed donor table documents two analyzed compartments: basal (5 animals, 5,182 cells) and corneum (5 animals, 5,023 cells). Both pass the matched-test floor using all available ages, but neither label has an exact reference match. There is also no defensible mature-only window from three 1.1-year and two 11.3-year animals. NMR-skin-vs-mouse and NMR-skin-vs-human therefore both have **0 matched pairs and were not run**.

For GSE214390 spleen, the expected raw directory `/home/reza/projects/agescape-revision-data/tmp_organs/nmr_gse214390` is still absent. No files matching the seven GSM identifiers or the publication supplement number were found under the project data root. The committed donor table documents 7 animals and five annotated T/ILC subsets; all five are slope-fittable and pass the minimal matched-test floor using all ages. However, none has an exact mouse or human reference match, and raw matrices are unavailable for a new expression refit. NMR-spleen-vs-mouse and NMR-spleen-vs-human therefore both have **0 matched pairs and were not run**. Existing committed result tables were used only to audit counts, never to invent new slopes.

### Matched cell-type pairs and per-gene outputs

There are **no matched cell-type pairs used**, because there were no eligible pairs. Consequently there are no per-species-pair CSVs with the reference output columns, no pooled pair-gene rows, and no floored/unfloored sign-agreement or correlation estimates. Assigning NA to those statistics is intentional; zero is not substituted for an unrun test.

## Secondary: mature-to-old concordance refits

### Feasibility and methods

Both refits are possible under the rules used for the mouse/human mature concordance analyses:

- Fly uses 30+ days. `run_reviewer_phase1_gene_concordance.fly_strata` already applies this restriction before pseudobulk aggregation. Its existing cache was verified to contain only 30/50/70-day pools and was reused.
- Worm uses 6+ days. `worm_strata(mature_only=True)` now removes day-1 cells before pseudobulk aggregation. A separate external cache, `concordance-cache/worm/mature-6plus`, prevents reuse or overwrite of the full-range cache. The cache was rebuilt from the resolved raw H5AD.
- Raw counts are aggregated by pool and cell type, requiring at least 20 cells in a pool pseudobulk and at least 5 pools to construct a unit; concordance retains units with at least 6 pools. Expression is log1p counts per 10,000. Detection filtering is recomputed after the age restriction: detection in at least 50% of pool pseudobulks in at least half of a stratum's units.
- Per-gene age slopes use continuous age plus sex and mean log UMI where variable. The statistic is mean pairwise Pearson correlation among unit slope vectors. Each of 1,000 permutations assigns one permuted age per pool jointly across all its cell types, within cohort, sex, and platform strata. P values are one-sided greater-tail values with the plus-one correction. Seeds are 20261993 for fly and 20261988 for worm.

### Results read from the output CSVs

| Species | Mature scope | Units | Pools | Unit pairs | Detected genes | Observed mean Pearson | Joint null mean (SD) | Null 95% interval | Joint p |
| --- | --- | ---: | ---: | ---: | --- | ---: | --- | --- | ---: |
| Fruit fly | 30+ days | 99 | 72 | 2,491 | 3,414 head; 4,707 body | 0.2191073747 | 0.0381143239 (0.0192844551) | 0.0129765332 to 0.0863783341 | 0.0009990010 |
| C. elegans | 6+ days | 10 | 9 | 45 | 12,723 | 0.3620947004 | 0.1600886478 (0.0743010658) | 0.0421104308 to 0.3127656922 | 0.0119880120 |

Fly strata are head: 59 units, 36 pools, 1,711 pairs, observed 0.2580904213, null mean 0.0394030949 (SD 0.0265476624), p = 0.0009990010; and body: 40 units, 36 pools, 780 pairs, observed 0.1335945634, null mean 0.0352872891 (SD 0.0192827451), p = 0.0009990010. Both p values correspond to zero of 1,000 null values at or above the observed statistic. Worm p = 0.0119880120 corresponds to 11 of 1,000.

Worm produced 12 units with at least 5 eligible pool pseudobulks; 10 had at least 6 and entered the test. The small 9-pool design and platform-stratified null remain important power limits. No new multiple-testing family or cross-species evolutionary claim is introduced for these secondary within-species refits.

The previously committed worm mature CSV contained observed 0.2760037370, 13 units, 78 pairs, and p = 0.0149850150, but the loader on master had no explicit worm mature filter and the earlier cache/provenance was not reproducible from that code. The fresh pre-aggregation 6+ rebuild gives the different result above and supersedes that orphaned output. Fly's observed and null values reproduce the prior mature CSV exactly; its new `notes` field now records the 30+ scope explicitly.

Outputs:

- `results/followup-mature-concordance-fly.csv`
- `results/followup-mature-concordance-fly-strata.csv`
- `results/followup-mature-concordance-worm.csv`
- `results/followup-mature-concordance-worm-strata.csv`

## Rat-data check

No rat single-cell aging data were found anywhere under `/home/reza/projects/agescape-revision-data/`: recursive name checks found 0 rat/Rattus/norvegicus directories and 0 matching data filenames, and the three inventory/manifest files are human PBMC resources. A rat paper citation inside an XML literature file is not a local dataset. No external rat search was performed.

## Verification and reproducibility

- The audit was rerun after correcting NMR donor-table cell counting. All 232 saved cell-type rows were re-read; every same-tissue and any-tissue mouse/human match field is empty, and all eight summary rows report 0 matches.
- Fly and worm output CSVs were re-read after completion. From the mature caches, independent pseudobulk-slope and pairwise-correlation calculations reproduced fly head 0.2580904212821753, fly body 0.13359456343275292, their pair-weighted pooled value 0.2191073746653349, and worm 0.36209470042800107 to absolute tolerance 1e-12. Unit, pool, pair, and detection-gene counts also match.
- The four original full-range fly/worm concordance CSVs remained byte-identical by SHA-256. The only changed concordance results are the requested mature outputs.
- The complete project test suite passed: 214 passed, 1 skipped. `git diff --check` also passed.
- Run commands used the `agescape` conda environment. Fly: `python scripts/run_rev2_joint_concordance.py --species 'Fruit fly' --mature-only --permutations 1000`. Worm: the same command with `--species 'C. elegans' --refresh-cache`.
