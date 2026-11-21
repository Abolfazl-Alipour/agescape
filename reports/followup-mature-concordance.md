# Follow-up: mature-to-old concordance (audit item 6)

Scope: mouse gene-effect concordance refit at 18-plus months and human
seven-cohort concordance refit at 30-plus years, with the six-species
correction summary recomputed only from the mature inputs. Mouse mature
means 18 months and older; human mature means 30 years and older. An
analysis counts as mature-to-old only when both its derivation/selection
and its tested contrast use mature-and-older animals/donors only. The
restriction was applied before pseudobulk aggregation and before detection
filtering in both reruns.

## Verdicts

- Mouse within-tissue, primary gene set: observed mean Pearson 0.2105,
  joint null mean 0.1465 (SD 0.0272), p = 0.0220, 1000 joint permutations,
  59 units, 14 mice. BH q = 0.0440 in the mature two-test family and
  q = 0.0460 in the mature six-row mouse family. The nominal signal
  strengthened versus the original (p = 0.0440, q = 0.0879) and passes
  0.05 under both mature BH families.
- Mouse cross-tissue, primary gene set: observed 0.0604, joint null mean
  0.0297 (SD 0.0214), p = 0.1029, q = 0.1049 in the six-row mouse family.
  Still null, as in the original (p = 0.1528).
- Human seven-cohort, primary gene set: observed 0.1793, joint null mean
  0.2126 (SD 0.0344), p = 0.8152, q = 0.8152 in the mature two-test
  family, 416 donors, 127 units, 7 strata. Still null, as in the original
  (p = 0.7003, 480 donors).
- Six-species summary: the only mature-to-old members are mouse and
  human. BH over those two tests gives mouse q = 0.0440 and human
  q = 0.8152. Fly, worm, zebrafish, and naked mole-rat are NOT
  APPLICABLE as mature-to-old and are retained as exploratory only with
  their original p values.
- No test in this item came out UNTESTABLE. All 7 human strata and all
  7 retained mouse tissues were evaluable.

## Methods

Mouse script: scripts/run_rev2_concordance.py gained --mature-only and
--output-suffix flags (default behavior unchanged). In mature mode, cells
from mice younger than 18 months are removed before donor aggregation,
units are defined from the mature pseudobulks (at least 6 mature donors
per unit), the Phase 1 78-unit list is not forced, detection masks are
recomputed from the mature units, and the joint donor-level permutation
(one age assignment per mouse shared across all its cell types, permuted
within sex) is unchanged at 1000 permutations, seed 20261108. Scrublet
was ON (M1 settings, expected rate 0.06); zero scrublet failures.
Checkpoint is separate:
scratch/rev2_concordance_checkpoint_mature.npz. Outputs:
results/followup-mature-concordance-mouse.csv and
results/followup-mature-concordance-mouse-strata.csv. Run time 12:08
wall, peak RSS 4.1 GB, exit status 0.

Human script: scripts/run_rev2_joint_concordance.py gained --mature-only
and --output-suffix flags (default behavior unchanged). In mature mode
the donor age filter changes from 20-plus to 30-plus before pseudobulk
and detection filtering. The existing 20-plus cache was backed up first
to concordance-cache/human-20plus-backup-20261121/ and the cache was
rebuilt with --refresh-cache. Same seed (20261613) and 1000 permutations
per stratum. Outputs: results/followup-mature-concordance-human.csv and
results/followup-mature-concordance-human-strata.csv. Run time 56:44
wall, peak RSS 10.8 GB, exit status 0, no errors.

The four original files were verified byte-identical by sha256 before
and after (reviewer-rev2-concordance.csv,
reviewer-rev2-concordance-strata.csv,
reviewer-rev2-joint-concordance-human.csv,
reviewer-rev2-joint-concordance-human-strata.csv). Nothing under
report/ and no manuscript tex file was touched.

## Mature mouse results

| Scope | Gene set | Observed | Null mean (SD) | p | q, 6-row mouse family | Genes |
| --- | --- | --- | --- | --- | --- | --- |
| within-tissue | primary | 0.2105 | 0.1465 (0.0272) | 0.0220 | 0.0460 | 10633 |
| within-tissue | lenient | 0.2039 | 0.1445 (0.0261) | 0.0230 | 0.0460 | 14283 |
| within-tissue | unfiltered | 0.2028 | 0.1439 (0.0261) | 0.0230 | 0.0460 | 19955 |
| cross-tissue | primary | 0.0604 | 0.0297 (0.0214) | 0.1029 | 0.1049 | 10633 |
| cross-tissue | lenient | 0.0554 | 0.0275 (0.0197) | 0.1049 | 0.1049 | 14283 |
| cross-tissue | unfiltered | 0.0550 | 0.0271 (0.0195) | 0.1029 | 0.1049 | 19955 |

All six rows were reproduced exactly from the raw null draws in the
mature checkpoint (p, null mean, and null SD all match the CSV).
Within-tissue pairs 260, cross-tissue pairs 1451, units 59, donors 14.
Split-half reliability: 59 of 59 units estimated, mean Spearman-Brown
0.3042, no missing units. Noise ceiling within 0.2798 (observed over
ceiling 0.7523), cross 0.2190 (ratio 0.2759).

Donor pool: 14 mature mice, 4 female and 10 male. Female ages are 18
(x2) and 21 (x2) only, giving 6 distinct joint assignments, so the
female stratum null is coarse; the male stratum (ages 18 x2, 24 x4,
30 x4) gives 3150 assignments, not coarse. Total distinct joint
assignments across sex strata: 18900. The within-tissue result rests
on this small pool and the coarse female null; that is a limitation of
the mature mouse result, stated here rather than hidden.

Per-tissue mature donors (run log, post-QC, before the unit floor):
Fat 6, Heart_and_Aorta 9, Kidney 11, Limb_Muscle 12, Lung 10,
Marrow 13, Spleen 10. The same counts were recomputed independently
from the interim M1 QC files with the run filters and match exactly
on all 7 tissues; the recomputation is saved in
results/followup-mature-concordance-mouse-donors.json.

Retained tissues and units (mature units vs original units):
Fat 4 (was 4), Heart_and_Aorta 4 (was 5), Kidney 12 (was 13),
Limb_Muscle 7 (was 7), Lung 11 (was 12), Marrow 13 (was 14),
Spleen 8 (was 8). Total 59 of the original 78; 19 units dropped,
0 added. Tissues lost entirely because they have fewer than 6 mature
donors in every cell type: Bladder (4 mature donors), Mammary_Gland
(4), Thymus (5), Tongue (4), Trachea (0). Dropped units:
Bladder bladder cell, bladder urothelial cell, endothelial cell,
leukocyte; Heart_and_Aorta cardiomyocyte; Kidney brush cell; Lung
fibroblast of lung; Mammary_Gland B cell, T cell, basal cell,
endothelial cell, luminal epithelial cell of mammary gland,
macrophage, stromal cell; Marrow plasma cell; Thymus DN4 thymocyte,
professional antigen presenting cell; Tongue basal cell of
epidermis, keratinocyte.

Per-tissue within-tissue joint p, descriptive only (no BH family was
declared across the 7 strata): Fat 0.918, Heart_and_Aorta 0.467,
Kidney 0.409, Limb_Muscle 0.001, Lung 0.144, Marrow 0.019, Spleen
0.411. Limb_Muscle (observed 0.4284 vs null 0.1902) and Marrow
(observed 0.2923 vs null 0.1887) carry most of the mature
within-tissue signal; no per-tissue claim is made beyond that
description.

## Mature human results

Summary: observed 0.1793, joint null mean 0.2126 (SD 0.0344, 95 pct
interval 0.1538 to 0.2861), p = 0.8152, 1000 permutations, 1326 pairs,
127 units, 416 donors, 7 strata.

Per-cohort table (mature donors vs original donors; mature units vs
original units; mature p vs original p):

| Cohort | Donors | Units | p (mature) | p (original) |
| --- | --- | --- | --- | --- |
| HBCC_DLPFC | 64 vs 88 | 20 vs 20 | 0.648 | 0.459 |
| ImmuneAging_BAL | 8 vs 12 | 26 vs 29 | 0.825 | 0.440 |
| Kidney_LivingDonor | 19 vs 19 | 6 vs 6 | 0.060 | 0.060 |
| Kidney_SCR_v15 | 18 vs 18 | 9 vs 9 | 0.143 | 0.143 |
| Lung_CellRef | 13 vs 17 | 15 vs 16 | 0.406 | 0.326 |
| Lung_Integrated | 86 vs 114 | 31 vs 31 | 0.417 | 0.741 |
| MSSM_DLPFC | 208 vs 212 | 20 vs 20 | 0.924 | 0.873 |

Kidney_LivingDonor and Kidney_SCR_v15 are unchanged because all of
their donors were already 30-plus; their rows are identical to the
original. ImmuneAging_BAL keeps only 8 donors at 30-plus, the
smallest stratum, and still evaluates (26 units, p = 0.825). No
stratum was UNTESTABLE. Total donors fell 480 to 416 (64 lost, all
under 30); units fell 131 to 127; pairs fell 1422 to 1326.

## Six-species correction summary

Original family verification (recomputed from
results/reviewer-rev2-concordance.csv, primary within-tissue joint p,
one per species): sorted p values are fly 0.000999, worm 0.002997,
mouse 0.043956, zebrafish 0.360639, naked mole-rat 0.700, human
0.700300. Benjamini-Hochberg q values are fly 0.005994, worm
0.008991, mouse 0.087912, zebrafish 0.540959, naked mole-rat 0.700300,
human 0.700300. This confirms the locked original statement, mouse
q = 0.0879.

Mature family (mature-to-old members only, 2 tests): mouse
within-tissue primary p = 0.021978, q = 0.043956; human
p = 0.815185, q = 0.815185. Full table in
results/followup-mature-concordance-six-species.csv.

Secondary sensitivity, BH over the six mature mouse rows: within-tissue
q = 0.045954 for all three gene sets; cross-tissue q = 0.104895 for
all three gene sets. For comparison, the same family on the original
rows gives within-tissue q = 0.091908 and cross-tissue q = 0.152847.

Not applicable as mature-to-old (retained as exploratory only, original
p values, not rerun, not mature proof): fruit fly p = 0.001, worm
p = 0.003, zebrafish p = 0.361 (testis germline caveat remains),
naked mole-rat exact p = 0.70 (5 animals, 1.1 vs 11.3 years only; do
not claim a null beyond this tiny exploratory contrast).

All q values were recomputed from raw p values with statsmodels
multipletests and are recorded in
results/followup-mature-concordance-bh-check.csv; every q in this
report matches that file.

## Power loss (MDE80)

MDE80 = 2.4865 times null SD (one-sided alpha 0.05, 80 pct power,
normal approximation on the mean-Pearson statistic). Detectable
threshold = null mean plus MDE80. Full table in
results/followup-mature-concordance-mde80.csv.

| Analysis | MDE80 original | MDE80 mature | Observed minus null, original | Observed minus null, mature |
| --- | --- | --- | --- | --- |
| Mouse within | 0.0625 | 0.0677 | 0.0477 | 0.0640 |
| Mouse cross | 0.0370 | 0.0532 | 0.0131 | 0.0307 |
| Human | 0.0762 | 0.0855 | -0.0185 | -0.0334 |

- Mouse within: MDE80 rose 8 pct (0.0625 to 0.0677) while the
  observed-minus-null effect rose 34 pct (0.0477 to 0.0640). The
  mature observed value 0.2105 sits 0.0037 below the normal-approx
  detectable threshold 0.2142, yet the permutation p = 0.022, so the
  normal MDE80 is slightly conservative here. Net: power slightly
  reduced, effect slightly larger, verdict strengthened.
- Mouse cross: MDE80 rose 44 pct (0.0370 to 0.0532); the effect rose
  135 pct but remains below MDE80 (0.0307 vs 0.0532). p moved 0.153
  to 0.103 and stays null. The restriction cost real power here.
- Human: MDE80 rose 12 pct (0.0762 to 0.0855). The observed value
  sits below the null mean in both the original and the mature run,
  so no positive effect of the original size is detectable and none
  is present in the data; p moved 0.700 to 0.815 and stays null.

## What changed versus the original verdicts

- Mouse within-tissue: p 0.0440 to 0.0220; q 0.0879 (six-species
  family) to 0.0440 (mature two-test family) and 0.0460 (mature
  six-row mouse family). The nominal signal is now significant under
  both mature BH families, where the original was nominal-only.
- Mouse cross-tissue: p 0.1528 to 0.1029; still null.
- Human: p 0.7003 to 0.8152; still null, on 416 instead of 480
  donors.
- The mouse strengthening is real in the numbers but rests on 14
  mice with a coarse female null; do not overstate it beyond the
  p and q values above.

## Deviations, handover, and unverified items

- Kimi (kimi-code/k3.1-pro) implemented both scripts, ran the mature
  mouse job, and launched the mature human job, then died on its
  5-hour usage quota at about 13:13 PST Nov 21. No Kimi work was
  lost: both heavy jobs were running under nohup with log files.
- Antigravity Claude Opus 5.5 High: 429 on launch, individual quota
  reached, resets in about 140 hours. Never ran.
- Antigravity Gemini 4 Neon XHigh: verified the mouse checkpoint
  reproduction, the original-file checksums, and the script diff,
  then polled the human run and died on its own 429 at about 13:42
  PST Nov 21 (resets in about 4.5 hours).
- The runner (Muse Spark subagent) finished the remaining work
  directly: BH families, MDE80, six-species table, this report, and
  the commit. No heavy compute was rerun; every number above was
  read from the output CSVs or recomputed from the saved null draws.
- Unverified items: none material. The female-stratum coarse null
  is a stated limitation, not a gap. Per-tissue mouse donor counts
  were recomputed from the interim M1 QC files and match the run log
  exactly on all 7 tissues. The mature human cache was rebuilt from
  the same local h5ad files as the original run, not re-fetched
  from Census.

## Files

- scripts/run_rev2_concordance.py (added --mature-only,
  --output-suffix)
- scripts/run_rev2_joint_concordance.py (added --mature-only,
  --output-suffix)
- results/followup-mature-concordance-mouse.csv
- results/followup-mature-concordance-mouse-strata.csv
- results/followup-mature-concordance-human.csv
- results/followup-mature-concordance-human-strata.csv
- results/followup-mature-concordance-six-species.csv
- results/followup-mature-concordance-bh-check.csv
- results/followup-mature-concordance-mde80.csv
- results/followup-mature-concordance-mouse-donors.json
- reports/followup-mature-concordance.md (this report)
- scratch/mature-concordance/mouse.log (mouse run log, 12:08 wall)
- scratch/mature-concordance/human.log (human run log, 56:44 wall)
- scratch/rev2_concordance_checkpoint_mature.npz (mouse null draws)
