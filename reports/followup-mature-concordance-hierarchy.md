# Hierarchical concordance: Level 0, Level 1, Level 2, Level 3

Branch: `followup-mature-concordance-hierarchy` (base `0df5da2`).
Script: `scripts/run_mature_concordance_hierarchy.py`.
Outputs: `results/followup-mature-hierarchy-level0.csv` (311 rows),
`results/followup-mature-hierarchy-level2.csv` (6 rows),
`results/followup-mature-hierarchy-level2-details.csv` (21 rows),
`results/followup-mature-hierarchy-level2-human-pairs.csv` (6,675 pairs),
`results/followup-mature-hierarchy-level2-human-cohort-pairs.csv` (21 cohort pairs),
`results/followup-mature-hierarchy-level2-fly-pairs.csv` (2,016 pairs),
`results/followup-mature-hierarchy-level0-gene-set-sensitivity.csv`,
`results/followup-mature-hierarchy-level0-fidelity.csv`.
Tests: `tests/test_mature_concordance_hierarchy.py` (16 tests).
Full suite: `python -m pytest tests/ -q` gives 230 passed, 1 skipped.

No tex file was edited. No existing result CSV was overwritten.
Level 1 and Level 3 reuse committed results; they were not rerun.
Level 0 was computed fresh; Level 2 was computed fresh only where missing
(human cross-cohort, fly head vs body).

Terms, defined on first use. A unit is one cell type in one tissue or cohort,
one separate test. A slope vector is the list of per-gene age slopes for one
unit. Concordance is the mean pairwise Pearson correlation of slope vectors
between units. BH is Benjamini-Hochberg correction, which turns p-values into
q-values that control the false discovery rate; q below 0.05 is the pass bar.
A pseudobulk is the summed counts of one donor's cells of one type, so the
donor (or animal, or pool) is the inferential unit, never the cell.

## Level 0: per-cell-type age-slope significance

Method. For every unit, the per-gene OLS (ordinary least squares) age model
of the concordance pipeline was refit on the same pseudobulks, same
covariates (sex where variable, mean log UMI where genuine), same mature
windows (mouse 18-plus months, human 30-plus years, fly 30-plus days, worm
6-plus days; zebrafish and naked mole rat as in their original runs). BH
correction was applied WITHIN each unit over the gene set the pipeline used
for that unit's slope vector. The gene sets: mouse uses one species-wide
primary mask (10,633 genes); human, fly, worm, and zebrafish use per-stratum
primary masks; naked mole rat uses the unfiltered set (20,273 genes), which
is the set its committed concordance run used.

The fit is a safe OLS (Frisch-Waugh-Lovell via QR), not the pipeline's
pinv(X.T @ X) fit. The pipeline fit squares the condition number and is
unstable on ill-conditioned designs. Verified on the fly data: the fly
mean_log_umi column is log(10000) plus float32 rounding noise, range
3.8e-06, a near-constant column. On the spot-checked fly unit (head Johnston
organ neuron), safe vs pipeline p-values differ by up to 0.418 and BH hit
counts differ (109 vs 94). Level 0 inference therefore uses the safe fit;
the pipeline fit is used only for the Level 1 fidelity check below.

Two guards, both verified in the code and the outputs. First, a unit whose
age slope is not identified (rank-deficient design, one age value, or age
collinear with sex) gets no p-values and is reported as NA, never as hits.
Second, a gene with zero standard error (constant expression, or an exact
fit where the residual sum of squares is at most 1e-15 of the total sum of
squares, which happens when a gene is detected in only one donor) gets
p-value 1. Without this guard those genes produce fake p-values near 1e-30.
Zero-SE gene counts on the primary gene set: mouse 8,346 genes in 57 of 59
units; human 24,390 in 85 of 127; fly 1,208 in 37 of 99; worm 2,923 in 9 of
10; zebrafish 10,905 in all 14 units; naked mole rat 4,386 in both units.

Verdict rule for this level, stated before the numbers: PASS means at least
one unit has at least one gene at q below 0.05. The median, total, and
concentration are reported alongside, because one hit and ten thousand hits
are not the same strength. Zebrafish and naked mole rat NULLs are unpowered,
not evidence of absence.

Results, primary gene set only, from
`results/followup-mature-hierarchy-level0.csv`:

1. Mouse: 59 units, 13 with at least 1 significant gene, median 0,
   total 2,884 significant genes. Top units are all limb muscle:
   mesenchymal stem cell 1,916 of 10,633 genes; endothelial cell 642;
   macrophage 211. Those three units hold 2,769 of 2,884 hits (96.0%).
   Verdict: PASS, concentrated in limb muscle.
2. Human: 127 units, 35 with at least 1 significant gene, median 0,
   total 8,771. Top units are all Mount Sinai School of Medicine (MSSM)
   brain: L2/3 intratelencephalic projecting glutamatergic neuron 1,676
   of 22,286 genes; L2/3-6 intratelencephalic projecting glutamatergic
   neuron 1,551; microglial cell 1,131. Verdict: PASS, concentrated in
   MSSM brain.
3. Fruit fly: 99 units, 87 with at least 1 significant gene,
   median 141.5, total 30,418. Top units: outer photoreceptor cell
   2,392 of 3,414 genes; adult fat body 2,204 of 4,707; follicle cell
   2,001 of 4,707. Seven units have no identifiable age slope and are NA
   (see Limitations). Verdict: PASS, broad across units.
4. C. elegans: 10 units, 5 with at least 1 significant gene, median 2,
   total 2,868. Embryonic_cells 2,358 of 12,723 genes and Germline 463
   hold 2,821 of 2,868 hits (98.4%). The remaining hits are Hypodermis 24,
   Intestine 19, Muscle 4. Verdict: PASS, concentrated in embryonic and
   germline cells.
5. Zebrafish: 14 units, 0 with at least 1 significant gene, median 0,
   total 0. Verdict: NULL, unpowered (6 fish).
6. Naked mole rat: 2 units, 0 with at least 1 significant gene, median 0,
   total 0. Verdict: NULL, unpowered (5 animals).

## Level 1: within tissue or cohort

Existing committed results, not rerun. Sources:
`results/followup-mature-concordance-mouse.csv`,
`results/followup-mature-concordance-human.csv`,
`results/followup-mature-concordance-fly.csv`,
`results/followup-mature-concordance-worm.csv`,
`results/reviewer-rev2-joint-concordance-zebrafish.csv`,
`results/reviewer-rev2-concordance.csv`.

1. Mouse within-tissue: observed 0.2105 vs null 0.1465, p 0.0220,
   q 0.0440, 59 units, 14 mice. Verdict: PASS.
2. Human within-cohort: observed 0.1793 vs null 0.2126, p 0.8152,
   127 units, 416 donors. Verdict: NULL.
3. Fruit fly within-stratum: observed 0.2191 vs null 0.0381,
   p 0.000999, 99 units, 72 pools. Verdict: PASS.
4. C. elegans: observed 0.3621 vs null 0.1601, p 0.011988, 10 units,
   9 pools. Verdict: PASS.
5. Zebrafish: p 0.3606, testis only. Verdict: NULL.
6. Naked mole rat: exact p 0.70, 5 animals. Verdict: NULL.

Fidelity check. The Level 0 pipeline refit reproduces these committed
Level 1 observed statistics from the same slopes; every check passes
(`results/followup-mature-hierarchy-level0-fidelity.csv`). Human, fly,
worm, zebrafish, and naked mole rat match to 1e-16 or better. Mouse
within-tissue mean recomputes as 0.2098 vs committed 0.2105 (difference
0.000667) and cross-tissue as 0.0605 vs 0.0604, inside the 1e-3 tolerance
for the mouse rebuild: the mouse strata are rebuilt from the raw Tabula
Muris Senis file and the rebuild no longer reproduces the original
doublet calls bit-for-bit (Fat 20 vs 15 committed, Kidney 77 vs 75,
the other 14 tissues exact). Units (59), primary genes (10,633), and
within-tissue pairs (260) match exactly.

## Level 2: whole species

Primary rows are in `results/followup-mature-hierarchy-level2.csv`.
Decompositions and sensitivities are in
`results/followup-mature-hierarchy-level2-details.csv`; they are computed
from the same null draws as their primary row, not independent tests.

1. Mouse cross-tissue (existing): observed 0.0604 vs null 0.0297,
   p 0.1029, 1,451 pairs, 59 units, 14 mice. Verdict: NULL.
2. Human cross-cohort (new): pairs of units whose cohorts differ, mean
   pairwise Pearson correlation of slope vectors, null by permuting ages
   independently within each cohort (and sex where variable), one shuffled
   age per donor shared across its cell types, 1,000 permutations, seed
   20261620. Pair-weighted primary: observed 0.0139 vs null 0.0001,
   p 0.0729, 6,675 pairs, 127 units, 416 donors. Verdict: NULL.
   Decompositions from the same draws: equal cohort-pair weighting
   (each of the 21 cohort pairs counts 1/21) gives 0.0161, p 0.0260;
   strict cross-organ (brain vs lung vs kidney, excluding the multi-organ
   ImmuneAging_BAL cohort) gives 0.0040, p 0.2747, 3,130 pairs;
   same-organ cross-cohort gives 0.0496, p 0.0519, 919 pairs; same-organ
   same-cell-type gives 0.2986, p 0.0010, 28 pairs; same-organ
   different-cell-type gives 0.0418, p 0.0759, 891 pairs; pairs involving
   ImmuneAging_BAL give 0.0132, p 0.2228, 2,626 pairs; leave-one-cohort-out
   is null for all cohorts except leaving Lung_Integrated out, 0.0219,
   p 0.0180, 3,699 pairs. The 21 cohort-pair rows are in
   `results/followup-mature-hierarchy-level2-human-cohort-pairs.csv`.
   The verdict depends on weighting: the pair-weighted primary is null,
   the equal-weighted version passes. The strongest human cross-cohort
   signal is the same cell type in different cohorts (28 pairs).
3. Fruit fly head vs body (new): 59 head units vs 40 body units. The head
   and body pools are the same 36 pools sampled in both organs: pool
   identities were inferred by stripping the head/body tag from donor IDs,
   giving a complete 1:1 match with identical age and sex in all 36 pairs
   (0 identical raw IDs). The null therefore permutes one age per
   head/body pool pair, shared across both organs and all cell types,
   within sex, 1,000 permutations, seed 20261995. Primary (92 identifiable
   units, 2,016 pairs): observed 0.0867 vs null 0.0122, p 0.0020.
   Verdict: PASS. Sensitivities all pass: all-strata gene set 0.0952,
   p 0.0020; technical-masked (mitochondrial and ribosomal genes dropped)
   0.0691, p 0.0020; independent-per-organ null p 0.0010. The primary uses
   92 units and 2,016 pairs, not the 2,360 pairs in the brief, because 7
   units have no identifiable age slope (see Limitations).
4. C. elegans, zebrafish, naked mole rat: single stratum each, so Level 1
   IS Level 2. Their Level 1 rows are carried verbatim: worm 0.3621,
   p 0.0120, 45 pairs, 10 units, 9 pools; zebrafish 0.3415, p 0.3606,
   91 pairs, 14 units, 6 fish; naked mole rat 0.6146, p 0.70, 1 pair,
   2 units, 5 animals. The naked mole rat p is exact: 10 assignments of
   2 aged among 5 animals, 6 of the 9 null draws at or above the observed
   value, (6+1)/10 = 0.70. Verdict: NOT APPLICABLE as separate Level 2
   tests.

## Level 3: across species

Existing committed results, not rerun.

1. Mouse vs human matched per-gene test: 28 matched cell-type pairs,
   pooled sign agreement 0.391 (floored) / 0.388 (unfloored), no
   correlation. Verdict: NULL (negative).
   Source: `reports/followup-mature-matched-p1m.md`.
2. Mouse vs human pathway test: formal verdicts 0 agree, 1 disagree,
   1,422 unclear of 1,423 function rows; floor-free sign agreement
   19.7% vs 38.6% expected under independence, so functions go
   systematically opposite. Verdict: NULL (opposite).
   Source: `reports/followup-mature-matched-p1m-pathways.md`.
3. Fruit fly, worm, zebrafish, naked mole rat: 0 matched cell-type pairs
   with mouse or human under the exact-name rule. Verdict: NOT
   APPLICABLE (applicability block, not a power problem).
   Source: `reports/followup-mature-other-species-pergene.md`.

## Full hierarchy verdict table

| Species | Level 0 per cell type | Level 1 within tissue | Level 2 whole species | Level 3 across species |
|---|---|---|---|---|
| Mouse | PASS (13/59 units, median 0, limb muscle) | PASS (p 0.0220, q 0.0440) | NULL (p 0.1029) | NULL (matched test negative) |
| Human | PASS (35/127 units, median 0, MSSM brain) | NULL (p 0.8152) | NULL (p 0.0729) | NULL (matched test negative) |
| Fruit fly | PASS (87/99 units, median 141.5) | PASS (p 0.0010) | PASS (p 0.0020) | NOT APPLICABLE (0 matched pairs) |
| C. elegans | PASS (5/10 units, median 2, embryonic/germline) | PASS (p 0.0120) | NOT APPLICABLE (single stratum) | NOT APPLICABLE (0 matched pairs) |
| Zebrafish | NULL (0/14 units, 6 fish) | NULL (p 0.3606) | NOT APPLICABLE (single stratum) | NOT APPLICABLE (0 matched pairs) |
| Naked mole rat | NULL (0/2 units, 5 animals) | NULL (p 0.70) | NOT APPLICABLE (single stratum) | NOT APPLICABLE (0 matched pairs) |

Reading the table: per-cell-type age signal exists in mouse, human, fly,
and worm, but it is concentrated (limb muscle in mouse, MSSM brain in
human, embryonic and germline cells in worm) and broad only in fly.
Within-tissue agreement holds in mouse, fly, and worm, not in human.
Whole-species agreement holds only in fly. Cross-species agreement holds
nowhere it could be tested.

## Methods (detail)

Level 0. `scripts/run_mature_concordance_hierarchy.py level0 --species
<name>` fits every unit once with the safe OLS and once with the pipeline
OLS (fidelity only), applies within-unit BH over each candidate gene set,
and writes per-species tables plus fidelity checks to
`scratch/hierarchy/level0-<slug>.csv`. The primary row per unit is the
pipeline gene set for its species. The 7 unidentifiable fly units stay in
the final CSV as NA rows so it has exactly one row per unit (311 total).
Level 0 is its own level: concordance was never recomputed on
significant-genes-only.

Level 2. `level2 --species Human|Fruit fly --permutations 1000` builds
slope vectors on the common-gene set, computes the mean cross-group
Pearson correlation, and runs the joint donor-level permutation null.
Human uses the independent-by-cohort scheme (donors are disjoint across
cohorts). Fly uses the paired_shared scheme (36 head/body pool pairs).
Per-iteration identifiability masking is counted; no fly iteration needed
masking beyond the 7 dropped units. The human null draws are stored as a
full pair matrix so all decompositions come from the same draws.

Assembly. `assemble` validates every fidelity file, the 311-unit counts,
and that each fraction equals its count divided by genes tested, then
writes the final CSVs. It raises rather than writing on any failure.

## Limitations

1. Fly Level 2 primary is 92 units and 2,016 pairs, not the 2,360 pairs
   the brief expected (59 head x 40 body). Seven units have no
   identifiable age slope and were dropped before the run: head medullary
   intrinsic neuron Mi9, neuron 5, neuron 6; body adult midgut-hindgut
   hybrid zone, adult reticular neuropil associated glial cell_body,
   gustatory receptor neuron, stalk follicle cell. Three body units have
   a single age value; four have age collinear with sex. A 99-unit run
   would need slope vectors for units with no identified slope, which do
   not exist. The 7 units are NA rows in the Level 0 CSV.
2. Mouse Level 0 signal is 96.0% limb muscle. The 30-month mice are a
   separate sequencing batch from the 18- to 24-month mice in the source
   data, and the cleaned files carry no chemistry or batch field, so a
   batch confound at 30 months cannot be checked from this data. This is
   an open question, not a finding.
3. Worm Level 0 signal is 98.4% embryonic and germline cells. The somatic
   signal is 47 genes across 3 units.
4. Zebrafish Level 0 (6 fish) and naked mole rat Level 0 (5 animals) are
   unpowered nulls. The naked mole rat concordance p has a floor of 0.10
   (10 assignments); its Level 0 null is likewise a power statement, not
   evidence of no age effect.
5. The human cross-cohort verdict depends on pair weighting. The
   pair-weighted primary is null (p 0.0729); equal cohort-pair weighting
   passes (p 0.0260). The 21 cohort-pair p-values and the decompositions
   are uncorrected slices of one null test, not independent findings.
6. Level 0 BH correction is within each unit only. No correction was
   applied across the 311 units, and the per-species PASS verdicts are
   descriptive (at least one unit with at least one significant gene),
   not family-wise claims.
7. Fly and worm mature concordance refits, and the fly/worm Level 1 rows
   used here, are exploratory: they were computed after the primary
   two-test mammalian family and are not part of it.
