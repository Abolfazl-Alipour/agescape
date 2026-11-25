# SAMap homology matching: mouse vs human validation (Stage 1)

**Branch:** `samap-homology-matching` (from master 43f9734)
**Date:** 2026-11-25
**Status:** STAGE 1 GATE FAILED. 5 of 28 name-matched pairs recovered (gate: >= 14).
Per the pre-stated protocol, work stops here. No fly, worm, zebrafish, or
naked mole rat arms were run.

## Scope change (Reza, Nov 25 01:08 PT)

Zebrafish and naked mole rat were removed from the Agescape paper entirely.
Their SAMap arms were dropped before any BLAST or mapping ran for them.
Their proteomes, input files, and prepped files remain on disk untouched
(complete files another step may reuse); nothing was deleted for them.
The default BLAST pair list now covers only hu, ms, dm, ce, and the driver
refuses dr/hg pairs.

## Inputs (Step 0, all 6 species built before the scope change)

Per-species mature-only single-cell AnnData, concordance-pipeline unit
labels (`tissue | cohort | unit`), subsampled to at most 3,000 cells per
unit, seed 20261125. All 311 hierarchy Level 0 units are present
(mouse 59, human 127, fly 99, worm 10, zebrafish 14, NMR 2).
Human input uses the 7 M9 cohorts only (MSSM, HBCC, Kidney x2, Lung x3);
Yazar blood is NOT in the human input, which matters for validation below.

| Species | Cells | Units | Genes in input | Matched in FASTA | Prepped genes |
|---------|------:|------:|---------------:|-----------------:|--------------:|
| Human | 306,050 | 127 | 58,140 | 20,331 (35.0%) | 20,331 |
| Mouse | 81,903 | 59 | 20,138 | 15,792 (78.4%) | 15,792 |
| Fly | 142,130 | 99 | 15,992 | 12,660 (79.2%) | 12,660 |
| Worm | 20,117 | 10 | 26,454 | 19,143 (72.4%) | 19,143 |
| Zebrafish | 27,501 | 14 | 30,390 | 19,741 (65.0%) | 19,741 |
| NMR | 6,000 | 2 | 20,195 | 11,541 (57.1%) | 11,541 |

Notes on the match fractions. Human 35.0% is a denominator effect: the
58,140 input features include ~35,000 non-protein-coding features, and the
protein FASTA covers 23,869 genes; 20,331 of 23,869 protein-coding genes
(85%) matched. Zebrafish rose from 55.1% to 65.0% after a bug fix: the FASTA
parser split gene symbols on ":", which corrupted zebrafish symbols that
themselves contain colons (`si:ch211-...`, `zgc:...`); the split is now
fly-only. The remaining zebrafish unmatched genes are mostly GRCz10 to
GRCz11 symbol renames. Mouse unmatched genes are mostly Rik/Gm/LOC predicted
labels. NMR unmatched genes (8,654) are not rescuable by symbol: the
expression IDs come from the HetGla_female_1.0 annotation while the only
available Ensembl proteome is the renamed maternal assembly; the release-102
old-assembly FASTA returns 404.

Two build-script bugs were found and fixed (scripts, not results): the
manifest grouped mouse units by cell-type name only (45 rows instead of 59
tissue-by-unit rows), and `build_samap_inputs.py` was once overwritten empty
by a double-`cat` transfer and restored from the local copy (408 lines,
verified).

## BLAST (reciprocal, e-value < 1e-6)

**Method change, documented (thresholds untouched).** The first full-isoform
hu-to-ms run was measured at ~1,000 queries/min: 245,535 human proteins
would take ~4 hours for one direction and ~13 hours for all remaining arms.
It was killed at 22% (54,565 of 245,535 queries, 266 MB partial output
deleted as incomplete). The driver now blasts ONE representative protein per
expressed gene (the longest isoform). SAMap collapses protein hits to gene
level by max bit score, so the gene-level graph it consumes is the same up
to rare isoform-specific edge cases (a shorter isoform aligning better than
the longest). E-value threshold, match rule, and gates are unchanged. The
full-isoform run remains available as a fallback if validation had failed
for that reason; it did not (see attribution below).

Representative FASTAs (`data/samap/maps/rep_<sid>.fa`): hu 20,331, ms 15,792,
dm 12,660, ce 19,143 proteins (one per expressed gene, all genes covered).
BLAST hu-ms completed in ~23 min total (vs ~5 h projected full-isoform).

| Direction | Rows | Unique query genes with hits |
|-----------|-----:|-----------------------------:|
| hu to ms | 778,360 | 17,962 of 20,331 (88%) |
| ms to hu | 689,650 | 15,603 of 15,792 (99%) |

Sanity checks on the graph: reciprocal best-hit gene pairs 14,445
(literature one-to-one human-mouse ortholog count is ~15,134). Spot checks
all correct: human ACTB to mouse Actb, GAPDH to Gapdh, TP53 to Trp53, MSN
to Msn (ENSMUSG00000031207; an earlier recalled ID, 31200, is Mtcp1 and was
wrong, the graph was right).

## SAMap mouse vs human (Stage 1)

`scripts/run_samap_pair.py hu ms`, SAMap v3 (sc-samap 3.0.1), keys = the
concordance units, `save_processed=False`. The run took ~92 min wall
(388k cells; memory peaked at 94% of Beast's 38 GB twice and survived).
Note for later runs: the script now calls `sm.run(umap=False)` because UMAP
is only a visualization projection and the mapping scores do not use it.

| Quantity | Value |
|----------|------:|
| Human cells / clusters | 306,050 / 127 |
| Mouse cells / clusters | 81,903 / 59 |
| Genes matched (hu / ms) | 20,331 / 15,792 |
| Pairs scored (unordered) | 7,493 |
| Pairs passing the full match rule | 27 |

Match rule (pre-stated, all must hold): alignment score >= 0.30;
reciprocal best hit at cluster level; >= 40 enriched gene pairs
(GenePairFinder, align_thr 0.1); score in the top 1% of pairs in that
comparison. Top-1% threshold in this comparison: 0.2076. The lowest passing
score is 0.3146 (vein endothelial cell); 45 pairs score >= 0.30, 33 are
reciprocal-best, 163 carry >= 40 gene pairs, 75 are top-1%, and 27 pass all
four rules at once.

The 27 passing pairs are in `results/samap-hums/matches.csv` (passing rows)
and `results/samap-hums/gene_pairs.csv` (enriched gene pairs). The strongest
are exact-name kidney matches (fenestrated cell to endothelial cell 0.894;
loop of Henle thick ascending limb 0.859; distal convoluted tubule 0.844;
non-classical monocyte 0.759). Several passing pairs are cross-tissue
(e.g. lung progenitor cell to marrow naive T cell 0.787; lung mast cell to
marrow naive B cell 0.604), which the reciprocal-best rule allowed because
each side's best partner was the other.

## Validation against the 28 name-matched pairs: GATE FAIL

`scripts/validate_samap_hums.py` compares the 27 passing pairs against the
28 name-matched pairs in `results/followup-mature-matched-p1m.csv`
(28 distinct pairs; 55 rows = 27 floored + 28 unfloored arms).
Full per-pair diagnostics in `results/samap-hums/validation.csv`.

**Recovered: 5 of 28. Gate (>= 14) FAILS. Work stops here per protocol.**

The 5 recovered: kidney distal convoluted tubule (SCR), kidney loop of
Henle thick ascending limb (SCR), lung alveolar macrophage (Integrated),
lung non-classical monocyte (Integrated), lung vein endothelial cell
(Integrated). All five are exact-name, same-tissue matches.

### Attribution of the 23 misses

1. **10 pairs are untestable: one or both clusters are absent from the
   SAMap input.** The input was built from the concordance-pipeline roster
   (311 units), while the validation pairs come from the matched test's
   broader human roster. Absent human clusters: Kidney B cell (both kidney
   cohorts), BAL and CellRef B cell, CellRef classical and non-classical
   monocyte, LivingDonor distal convoluted tubule, Lung_Integrated T cell
   (also absent on the mouse side: no Lung T cell unit). Absent entirely:
   both Yazar blood pairs (hematopoietic precursor cell, naive B cell),
   because Yazar blood is not among the 7 human input cohorts. That is 10
   of 28 untestable by construction, not by SAMap.
2. **Of the 18 testable pairs, 5 pass.** Among those 18: 8 score >= 0.30,
   5 are reciprocal-best, 10 have >= 40 gene pairs, 9 are top-1%.
3. **The binding sub-rule is reciprocal best hit.** Three testable pairs
   score >= 0.30, carry 660 to 970 enriched gene pairs, and sit in the
   top 1%, but fail ONLY on reciprocity because SAMap matched each to a
   near-synonym unit instead of the exact-name unit: kidney epithelial
   cell of proximal tubule (SCR, 0.692) matched to kidney proximal
   convoluted tubule epithelial cell (0.749, a passing pair); kidney
   collecting duct principal cell (SCR, 0.524) matched elsewhere; lung
   CD8-positive alpha-beta T cell (Integrated, 0.375) matched to spleen
   mature NK T cell (0.550, a passing pair). The remaining testable misses
   score below 0.30 (kidney T cell 0.025 to 0.088; lung B cell 0.072;
   alveolar macrophage CellRef 0.002; classical monocyte 0.097 to 0.202;
   non-classical monocyte BAL 0.026; vein endothelial CellRef 0.080;
   proximal tubule LivingDonor 0.004).
4. **Even on the testable denominator, 5 of 18 is below the gate.**
   Relaxing the rule would not fix this: the three reciprocity-only
   failures would bring it to 8, still below 14, and the rule was
   pre-stated and is not changed after seeing results.

### What this means (interpretation, not a re-litigation)

SAMap does find real mouse-human cell-type matches (27 passing pairs, top
scores 0.89, 0.86, 0.84 on exact-name kidney types, 14,445 reciprocal
best-hit gene pairs in the underlying graph). What it does not do is
reproduce the exact-name pair list: broad immune types (B cell, T cell,
classical monocyte) score near zero between the same-named units, and
several epithelial/T-cell pairs land on near-synonym neighbours instead.
Two structural causes are visible in the data: (a) the validation pairs
use a broader human roster than the concordance roster the SAMap input was
built from, so 10 pairs were never in the run; (b) same-name does not mean
same-cluster under SAMap when annotations differ in granularity across
cohorts (proximal tubule vs proximal convoluted tubule; CD8 T cell vs
mature NK T cell). Whether a looser rule, a broader input roster, or a
different validation set would recover more is a new decision for Reza,
not a change made here.

## What was NOT run (stopped by the gate)

- Fly arms (mouse vs fly, human vs fly) and the label-permutation negative
  control: no BLAST tables were built for dm pairs, no SAMap runs started.
- Worm arms (mouse vs worm, human vs worm): not started.
- Matched slope test (Stage 3): not started; its gate (>= 3 fly pairs AND
  >= 150 pooled P1m-ortholog rows) was never reached.
- MetaNeighbor AUROC corroboration: not run; it belongs to Stage 2, which
  the failed gate never opened.
- Zebrafish and NMR arms: dropped Nov 25 01:08 PT at Reza's order before
  any BLAST or mapping; their files remain on disk.

## Deliverables on branch `samap-homology-matching`

- `scripts/build_samap_inputs.py` (input builder, 311 units, seed 20261125)
- `scripts/build_samap_gene_maps.py` (FASTA harmonization + matched-ID counts)
- `scripts/run_samap_blast.py` (reciprocal BLAST, representative FASTAs)
- `scripts/run_samap_pair.py` (SAMap pair run + match rule, umap=False)
- `scripts/run_samap_permutation.py` (label-permutation control driver, unused)
- `scripts/validate_samap_hums.py` (Stage 1 validation + per-rule diagnostics)
- `tests/test_samap_match_rule.py` (6 tests, all pass)
- `results/samap-hums/` (alignment_scores.csv, gene_pairs.csv,
  gene_pair_counts.csv, matches.csv, run_summary.json, validation.csv)
- `results/samap-gene-map-counts.csv` (per-species matched-ID counts; the
  working copy is `data/samap/gene_map_counts.csv`, gitignored with the
  rest of `data/`)
- `data/samap/` (proteomes, inputs, prepped, maps incl. rep_*.fa; dr/hg files
  present but unused)
- This report: `reports/followup-samap-homology-matching.md`

## Environment and cost

- Beast, conda env `samap` (sc-samap 3.0.1, BLAST+ 2.5.0+), agescape env
  for builds (Python 3.11.17). MetaNeighbor installed in R (unused; gate
  stopped Stage 2).
- Wall time to this point: ~11 h from task start (00:35 PT) to gate verdict
  (11:40 PT), dominated by the SAMap hu-ms run (~92 min) and the killed
  full-isoform BLAST (~55 min) plus input builds.
- Beast disk: 14 GB free at 98% (was 17 to 18 GB at start). A 1 TB ZFS disk
  was added Nov 25 01:01 PT but is not yet mounted (needs Reza's sudo);
  no action taken here.
