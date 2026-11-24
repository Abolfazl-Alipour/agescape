# Mature HBCC interneuron refit: provenance

## Verdict
TRACEABLE. The mature HBCC per-gene age slopes are a genuine mature-only
(donors 30 years and older) refit by a named script from an existing input
file. They were not copied from the full-range run.

## The conflict
- reports/followup-11pt-block9-mature-interneuron.md line 5: INCONCLUSIVE.
  Every subtype UNTESTABLE for mature-only MSSM-to-HBCC replication because
  /home/reza/projects/agescape-revision/scratch/hbcc.h5ad is absent.
  Line 29: "The HBCC arm stopped at this path. No download or source
  substitution was attempted." Donor table: HBCC mature N UNTESTABLE for
  all four subtypes.
- reports/followup-mature-interneuron.md line 5: PROVES. Mature replication
  survives in VIP, pvalb, sst. Donor table: HBCC mature N 59 (VIP),
  59 (pvalb), 61 (sst), 51 (GABA); ages 30 to 85, median 48.0. Replication
  table (mature rows): VIP r 0.7523, both-FDR 4; pvalb r 0.6299, both-FDR 8;
  sst r 0.6996, both-FDR 41. Tex Methods line 180 and Results line 413
  report the same numbers.

## The script
scripts/run_followup_mature_interneuron_replication.py (Item 7 mature MSSM
and HBCC interneuron replication). Its fit_hbcc_unit() calls
block9.fetch_mature_unit_h5ad() and run_rev2_de_interneuron.run_unit_de(),
the same donor-pseudobulk voom-WLS machinery as the original, with the sole
selection change age 30 to 120. MSSM units were validated against the
block-9 CSV cache and only refit on cache failure; HBCC units were fit fresh.

## The input file
HBCC_PATH in the Item-7 script is
/home/reza/projects/agescape-revision-data/human-h5ads/5c97eeeb-7e52-44b3-b010-b832b1f5424c.h5ad.
It exists: 14,150,652,636 bytes, dated Nov 9 19:04. The script's
validate_hbcc() asserts shape (1,486,324 cells, 34,176 genes), 300 donors,
ages 2 to 85, and per-type cell counts VIP 74,520, pvalb 67,895, sst 71,641,
GABA 31,968. The run record run_summary_mature.json (Nov 21 20:45) carries
the same validation payload, including the identical size and
disease-normal count 722,982.

## Why block 9 said UNTESTABLE
scripts/run_followup_11pt_block9_mature_interneuron.py sets HBCC_PATH to
REPO / "scratch" / "hbcc.h5ad", i.e.
/home/reza/projects/agescape-revision/scratch/hbcc.h5ad. That path was never
populated and is still absent (checked Nov 24). Block 9 stopped there.
Item 7 used the real HBCC object at the UUID path instead. The two scripts
share the fetch and fit machinery; only the HBCC_PATH constant differs.
Block 9's verdict describes its own missing input path, not the absence of
HBCC data.

## Corroboration
- Scratch dir
  /home/reza/projects/agescape-revision-data/tmp_organs/followup-11pt-block9-mature-interneuron/
  holds per-unit caches counts_hbcc_mature30_*.npz and
  de_hbcc_mature30_*.npz for all four subtypes, plus
  item7_replication_run.log and run_summary_mature.json.
- results/followup-mature-interneuron-replication.csv (mature rows) matches
  the tex exactly: VIP r 0.7522709920913488, both-FDR 4, hypergeometric
  expected 0.5972073039742213; pvalb r 0.6298883175716962, both-FDR 8,
  expected 0.9972175848636616; sst r 0.6996293500214019, both-FDR 41,
  expected 11.664841182913472; block-permutation p 0.000999 for each.
- The moderated mature intersection stays UNVERIFIED (recovery gate failed,
  disclosed in the report). That is a methods caveat, not a provenance gap.

## Remaining gap (not blocking)
Block 9's report was never corrected to point at the UUID path, so the two
reports still contradict each other on the page. The fix is a one-line
correction in the block-9 report (or a pointer to this note), not a
re-analysis. No tex change is needed: the tex numbers are backed by the
Item-7 files above.
