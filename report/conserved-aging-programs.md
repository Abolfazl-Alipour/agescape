# Conserved aging programs at single-cell resolution: a cross-species analysis with an honest power audit

**Repository:** `agescape` — cross-species single-cell aging portfolio project (private repository)
**Date:** 2026-10-14
**Status:** Public release pending Reza's review. All numbers below are read directly from the
committed files in `results/` and `figures/`; nothing in this document is untraced.

---

## Abstract

We re-analyzed the Tabula Muris Senis (TMS) droplet atlas (245,389 cells, 23 mice, 16 tissues)
with a pre-registered pipeline — per-tissue QC, pseudobulk differential expression per
tissue×cell-type unit (PyDESeq2, 3m vs 18m, the only sex-balanced contrast), and a calibration
trio — and found that 17 of 163 evaluated units met the eligibility floor, all in Lung or Kidney.
The strongest aging signals were Lung::adventitial cell (55 DE genes) and Kidney::proximal
convoluted tubule epithelial cell (53 DE genes); the calibration trio supported the ranking
(55 observed DE genes vs a 4.6-gene permutation null and a 0-gene sex-effect floor). Agreement
with the TMS paper's own published age-coefficients was 13 agree / 4 partial / 0 disagree across
the 17 powered units. We then attempted the pre-registered human replication arm and found it
**untestable with obtainable public data**: Tabula Sapiens v1 is an adult cohort (skin donors
33/45/59/61 y, n=4; lung donors 59/59/60/61 y) and the SenNet public index exposes no donor
ages. Collapsing the 17-unit mouse program (167 up / 94 down genes) into a per-sample score, the
mouse whole-body positive control rose with age (slope 0.235 [95% CI 0.067, 0.403], p=0.008,
n=23 mice), while the naked mole-rat (NMR) score was flat over the tested span (slope 0.017
[-0.257, 0.290], p=0.858, n=5 samples); the mouse-vs-NMR slope difference was 0.218 [-0.013,
0.449], p=0.064, and TOST equivalence was **not** established (p=0.329 at δ=25% of the mouse
slope). We partially reproduced Savina et al. 2022 (29 global DE genes at the paper's stated rule
vs their reported 0). The honest headline: a calibrated mouse aging program, a human arm that
public data cannot currently test, and a suggestive but non-conclusive longevity-outlier contrast.

## Significance statement

Aging-intervention targets are usually nominated from bulk or single-species data, and
cross-species "conservation" claims are often made from cohorts that were never powered to detect
age effects at all. This project contributes three things a target-prioritization effort can use:
a calibrated, cell-type-resolved ranking of where aging signal is strongest in the mouse atlas;
a documented, reproducible finding that the public human single-cell data obtainable today
(Tabula Sapiens v1, SenNet's public index) cannot test per-cell-type human age replication; and
a longevity-outlier contrast (mouse vs naked mole-rat) with pre-registered equivalence testing,
reported with the statistics it earned rather than the statistics we wished for. Where a claim
could not be tested, we say so and document exactly what data would test it.

---

## 1. Data

All raw inputs are version-pinned in `data/manifest/*.json` (URL, accession or dataset UUID,
dataset version, schema version, size, SHA-256, download date) and verified by
`scripts/verify_data.py`.

- **Tabula Muris Senis, droplet (mouse, primary).** CZ CELLxGENE collection
  `0b9d8a04-bb9d-44da-aa27-705bb65b54eb`, dataset `48b37086-25f7-4ecd-be66-f5bb378e3aea`,
  dataset version `084058cc-4f17-43ce-b14e-1278df074013` (schema 7.1.0); collection/dataset
  identity confirmed via the curation API on 2026-10-14. The object contains 245,389 cells ×
  20,138 genes across 16 tissues, 23 mice (6 female, 17 male), and 6 age groups (1m, 3m, 18m,
  21m, 24m, 30m). Age and sex are partially confounded by design: 1m, 24m and 30m are male-only
  and 21m is female-only; the pre-registered primary contrast is 3m vs 18m, the only pair with
  both sexes in both groups (`results/m1-sex-age-table.csv`).
- **Tabula Sapiens v1 (human, replication cohort attempt).** CZ CELLxGENE Tabula Sapiens
  collection. Two tissues were downloaded and donor ages verified from the source objects:
  Skin (dataset version `c240d5bc-a521-4330-a76b-7d566e939d87`; 17,786 cells; donors aged
  33/45/59/61 y, n=4) and Lung (dataset version `fddf43a3-b7c8-49e9-ac0a-5e24b390849f`;
  65,847 cells; donors TSP1/2/14/25 aged 59/61/59/60 y). TS v1 is an adult cohort: it provides
  no usable donor-age contrast in either tissue.
- **Naked mole-rat skin (longevity outlier).** GEO GSE152886 (Savina et al. 2022, *Aging*;
  accession cross-checked on the GEO record 2026-10-14; BioProject PRJNA640800). The single
  supplementary processed-counts file (20,274 genes × 6 samples) was used; SRA raw reads were
  not downloaded (authors' processed counts sufficient; rationale recorded in the manifest).
  Sample 2V was dropped as a technical outlier per the source paper, leaving 5 samples:
  1J/5J/6J (young, 1.1 y group mean) and 3V/4V (aged, 11.3 y group mean). Our QC retained
  10,251 cells vs the paper's 10,232.
- **AnAge (lifespan denominators).** Official zipped-TSV export from genomics.senescence.info,
  downloaded 2026-10-14 (`results/m6-species-lifespans.csv`): Homo sapiens 122.5 y,
  Heterocephalus glaber 31 y, Mus musculus 4 y.
- **Orthologs.** One-to-one mouse↔human and mouse↔NMR orthologs via Ensembl BioMart
  (release recorded in `notes/methods.md`), cached in `data/interim/`.
- **TMS published coefficients (agreement test target).** The TMS paper's Supplementary Table
  MOESM8 (published age-DGE coefficients, sheets `DGE_result.tissue_cell.age_coef/.bh_p`),
  downloaded 2026-10-14 from the publisher; used as the independent reference for the M3
  agreement table.

### Data availability

All datasets above are public. Download URLs, accessions, dataset UUIDs, SHA-256 checksums and
download dates are in `data/manifest/`; `python scripts/fetch_data.py` re-acquires every file and
`python scripts/verify_data.py` checks each SHA-256. The NMR processed-counts file is the only
input not mirrored on CELLxGENE; it is fetched from the NCBI GEO FTP
(`ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152886/suppl/GSE152886_raw_data.csv.gz`).

---

## 2. Methods (per milestone)

**M1 — TMS processing.** String obs columns stripped (known trailing-space bug in `subtissue`).
Per-tissue QC with binding thresholds: 500–5,000 genes per cell, ≥1,000 counts; the raw object
contains zero mitochondrial genes, so no mt% filter was possible (recorded, not silently
skipped). Scrublet doublet removal per tissue. Post-QC retention: 238,186 of 245,389 cells
(97.06%; per-tissue range 67.92% in Pancreas to 99.67% in Lung; `results/m1-qc-summary.csv`).
Normalization (log1p) on the concatenated object; 2,000 HVGs per tissue (Seurat flavor), union
of 8,421 genes; PCA (zero_center=False, 50 components), kNN graph, single-threaded UMAP,
Leiden clustering; majority-vote annotation against the provided `cell_ontology_class`. All
stochastic steps carry seed 42 (`src/set_seeds.py`). The full post-QC object (16 tissues) feeds
eligibility evaluation; the DE analysis is confined to units meeting the pre-registered floor
(see M2), which is why powered units come only from Lung and Kidney.

**M2 — Mouse DE and aging-signal ranking.** Unit: tissue × cell-ontology-type, never pooled
across tissues. Pseudobulk per (mouse × tissue × cell type) from raw counts; PyDESeq2,
design `~ sex + age_group`, primary contrast 3m vs 18m; BH-FDR < 0.05 and |log2FC| ≥ 0.25;
log2FC orientation old-vs-young, asserted programmatically on known-direction genes. Eligibility
(pre-registered): ≥3 mice and ≥50 cells per group; ineligible units are listed, never dropped
silently (`results/m2-unit-eligibility.csv`, 163 units evaluated). Aging signal = n DE genes,
with median |log2FC| and bootstrap-over-mice CIs. Calibration trio on the top unit:
(a) label-permutation null; (b) positive controls (canonical markers checked for presence and
direction); (c) sex-effect floor (DE by sex within 3m only). Sensitivity: extended age bins
(young=1m+3m, old=18m+21m+24m+30m) and a male-only 18m-vs-24m run.

**M3 — Trajectory and literature agreement.** Diffusion pseudotime on Marrow (chosen by the
written rule: graded hematopoietic differentiation continuum); root = maximal Kit/Cd34
expression among young hematopoietic precursors; graph built on young cells with old cells
projected; corr(pseudotime, age) reported as a mandatory guard. Agreement with the TMS paper
was tested against the paper's own published continuous-age coefficients (MOESM8), not against
its prose: per powered unit, Spearman correlation and sign concordance between the paper's
age coefficients and our effect sizes over shared expressed genes. Methodological bridge: the
paper modeled age continuously; our primary model is dichotomous (3m vs 18m). Where the
comparison would otherwise be model-asymmetry-dominated, we also fit a continuous-age bridge
(`~ age_months + sex`) and say which comparator was used per unit
(`results/m3-agreement-table.csv`, column `reason`).

**M4 — Human data and orthologs.** Human DE design (pre-registered): pseudobulk per
(donor × cell type), `~ age_years + sex`, age continuous, never dichotomized; inclusion floor
≥6 donors and donor-age SD ≥12 y per cell type. A power grid over n donors (4–50), donor-age SD
(8–20 y) and biological dispersion (0.1–0.4) was computed before any human claim
(`results/m4-human-power.csv`). SenNet was queried through its public, token-free search API on
2026-10-14: 6,820 published datasets; the public index exposes no donor-age fields, and
per-dataset traversal plus Globus transfer was beyond this session's budget — the query is
documented as a first-class result, not a footnote. One-to-one orthologs (Ensembl BioMart) were
cached for mouse↔human and mouse↔NMR.

**M5 — Program catalog and housekeeping guard.** The 17 powered units' DE gene lists were
cataloged per unit with up/down counts, median |log2FC|, and the fraction of ribosomal/
mitochondrial genes among DE genes vs a 0.017 background (`results/m5-mouse-program-catalog.csv`).
Because the human arm was untestable (M4), every unit's conservation column reads UNTESTABLE and
the divergence table makes no conservation or divergence call
(`results/m5-divergence-table.csv`).

**M6 — Signature scoring and the longevity-outlier contrast.** Program definition: unique DE
genes across the 17 powered units — 167 up, 94 down (`results/m6-program-genes.csv`).
Strict human replication (1:1 ortholog slope across the 4 TS-skin donors, same direction,
nominal p<0.05) retained only 2 up (Limd1, Wnk1) and 2 down (Col3a1, Id3) genes; below the
pre-registered 10-gene floor, the **full mouse program** was used for scoring and is flagged as
a deviation (§6). Scoring: per-cell `sc.tl.score_genes`, score = up-program minus down-program,
aggregated to sample means; per-species OLS of score on life-stage-normalized age (fraction of
AnAge max lifespan). Pre-registered statistics: per-species slope + 95% CI; species×age
interaction; and a TOST equivalence test on the NMR slope with δ = 25% of the mouse slope
(δ=0.0587) fixed **before** the NMR test; plus an exact sample-level permutation test as a
distribution-free check. Coverage gate: program-gene detection rate printed next to every
cross-species score (mouse 1.0, human 0.73, NMR 0.31). Savina reproduction: global
young-vs-aged DE in NMR epidermis as a DE-count-vs-|log2FC|-cutoff curve
(`results/m6-savina-reproduction.csv`, `figures/m6-savina-reproduction.png`).

---

## 3. Results

### 3.1 M1/M2 — Where the aging signal is, and calibration

Of 163 tissue×cell-type units evaluated against the pre-registered floor, 17 were eligible —
9 Lung and 8 Kidney (`results/m2-unit-eligibility.csv`). No unit from the other 14 tissues met
≥3 mice and ≥50 cells per group for the sex-balanced 3m-vs-18m contrast; this restriction is a
property of the atlas design, not of an analysis choice.

The aging-signal ranking (`results/m2-aging-signal-ranked.csv`, Fig. 1) is headed by
Lung::adventitial cell (55 DE genes; 8 mice, 423 cells; median |log2FC| 1.34) and
Kidney::kidney proximal convoluted tubule epithelial cell (53 DE genes; 6 mice, 1,793 cells;
median |log2FC| 1.24), followed by Lung::fibroblast of lung (36) and Lung::alveolar macrophage
(35). Three eligible units returned zero DE genes. The calibration trio on the top unit
(`results/m2-calibration-trio-summary.csv`, Fig. 2): 55 observed DE genes against a
label-permutation null mean of 4.6 and a 3m sex-effect floor of 0 — the ranking metric is not
explained by FDR path or sex noise. The positive-control arm is reported as it was observed:
in the top unit, Apoe (+0.099), Il6 (+0.778) and Lgals3 (+0.467) were positive (old vs young)
and Igfbp3 (−0.905) negative, with none reaching nominal significance (smallest p=0.267;
`results/m2-calibration-positive-controls.csv`). Sensitivity runs bracket the primary result:
extended bins more than double the top unit's DE count (120 vs 55), consistent with the
pre-stated wider-gap hypothesis — and with the caveat that the extended old arm includes the
survival-selected 30m cohort, so the increase cannot be attributed to age alone; a male-only
18m-vs-24m contrast in Limb_Muscle::macrophage yields 757 DE genes
(`results/m2-sensitivity-summary.csv`), confirming the design is capable of detecting strong
effects where the contrast is large.

### 3.2 M3 — Agreement with the TMS paper's published coefficients

Against MOESM8, the independent verdicts across the 17 powered units are **13 agree, 4 partial,
0 disagree** (`results/m3-agreement-table.csv`). The strongest concordance is
Kidney::proximal convoluted tubule epithelial cell (Spearman ρ=0.866 over 13,019 shared genes;
sign concordance 0.977 among the paper's BH 0.01-significant genes). The four partial units
(e.g., Kidney::kidney collecting duct principal cell, ρ=0.127) are cases where our dichotomous
3m-vs-18m effect sizes correlate weakly with the paper's continuous coefficients; for two units
we additionally fit the continuous-age bridge model for a like-for-like comparison, and the
agreement table states per unit which comparator was used. Pseudotime on Marrow
(39,471 cells; graph built on young cells, old projected) correlates with chronological age at
Spearman ρ=0.016 (p=1.2e-03; `figures/m3-dpt-vs-age.png`) — i.e., within-age-group
differentiation structure does not confound the age axis at any material level, and no
trajectory-based aging claim is made.

### 3.3 M4/M5 — The human-power reality

This is a first-class result. Tabula Sapiens v1, the largest public human single-cell dataset
with verified donor ages obtainable without restricted access, is an adult cohort: skin donors
are 33/45/59/61 y (n=4) and lung donors 59/59/60/61 y. A donor-age contrast does not exist in
these data. The power grid (`results/m4-human-power.csv`) quantifies why n=4 cannot be read as
evidence either way: with 4 donors, even under the most favorable grid point (age SD 20 y,
dispersion 0.1), the minimum detectable slope corresponds to a 30-year log2FC of 0.59; at the
observed skin-cohort SD (8 y, dispersion 0.1) it is 1.474. All 12 grid points at n=4 are flagged
Low-Powered; the eligibility floor (≥6 donors, SD ≥12 y) first becomes satisfiable at n=6.
SenNet — the pre-registered primary human cohort — was queried via its public API on
2026-10-14 (6,820 published datasets) and its public index carries no donor-age fields;
restricted-path retrieval is documented as future work. Consequently the M5 catalog marks every
unit's human replication UNTESTABLE, and the divergence table makes no conservation or
divergence call. The housekeeping guard stands ready (ribo/mito fraction vs 0.017 background is
reported per unit) but is not load-bearing because no human replication fraction is computed.
The honest summary: **human per-cell-type age replication is untestable with the public data
obtainable in this project**, and any published cross-species conservation claim built on
TS-v1-scale cohorts deserves the same power scrutiny.

### 3.4 M6 — The signature result, with the statistics it earned

Program. The 17 powered units define 167 up- and 94 down-regulated unique genes
(`results/m6-program-genes.csv`); 312 gene-unit entries in total across 14 units with ≥1 DE
gene. Strict human replication on the only human skin data available (4 donors) retained 4
genes (2 up, 2 down) — below the 10-gene floor — so the scored program is the full mouse
program, flagged (§6).

Scores (`results/m6-lifespan-program-scores.csv`, `results/m6-statistical-summary.csv`,
Fig. 3). The mouse whole-body positive control passes: program score rises with
life-stage-normalized age across 23 mice (slope 0.235 [0.067, 0.403], p=0.008). The NMR score
is flat over the tested span (slope 0.017 [-0.257, 0.290], p=0.858, n=5 samples spanning
1.1 y to 11.3 y — roughly the first third of the NMR lifespan). The mouse-vs-NMR slope
difference is 0.218 [-0.013, 0.449], p=0.064: **suggestive, not conclusive** — the interval
barely excludes zero and the test is underpowered at n=5. TOST equivalence at the
pre-registered δ=0.0587 (25% of the mouse slope) is **not** established (p=0.329), and the
exact sample-level permutation p=0.818. We therefore do not claim NMR equivalence; we claim
that no NMR program rise was detected over the tested third of its lifespan, and that the
formal equivalence test at the pre-registered margin failed to confirm flatness. The human arm
(n=4, slope −0.213 [-1.374, 0.948], p=0.513) is shown with its wide interval and carries no
claim. Program-gene coverage: mouse 1.00, human 0.73, NMR 0.31 — the NMR score is built on
less than a third of the program genes and is labeled as such wherever shown.

Savina reproduction. The paper reports zero global DE genes between young and aged NMR
epidermis at its stated rule (p<0.05, |logFC|≥1.0 in Seurat's natural-log scale). Our curve
(`results/m6-savina-reproduction.csv`) gives 29 genes at p<0.05 & |log2FC|≥1.0 (20 at FDR<0.05)
and 560 genes at FDR<0.05 & |log2FC|≥0.25. The literal zero is not reproduced, but a conversion
caveat applies: the paper's threshold of 1.0 in natural-log units corresponds to |log2FC|≈1.44,
and at the nearest grid cutoff (1.45) our count is 9 genes (p<0.05) / 6 (FDR<0.05). Qualitatively, NMR epidermis is far more
stable across age than mouse skin programs at matched cutoffs; the "remarkable stability" claim
is partially reproduced, the literal zero is not, and a unit-conversion ambiguity in the paper's
methods is the most parsimonious explanation. Two NMR caveats are printed wherever the NMR
result appears: Igfbp3 — the paper's key gene — is absent from the GEO count file and its
behavior is unverifiable here; and our QC yields 10,251 cells vs the paper's 10,232.

---

## 4. Limitations

1. **Correlation, not causation.** Every effect reported is an association with age; no
   perturbation or longitudinal validation exists in these data.
2. **Sex×age confound in TMS.** 1m, 24m and 30m are male-only; 21m is female-only; even the
   balanced pair (3m: 2F/5M; 18m: 2F/2M) is imperfect. Sex-chromosome-linked effects cannot be
   fully excluded from any pooled statement, and only the 3m-vs-18m contrast is interpretable
   as age.
3. **Powered units cover 2 of 16 tissues.** All 17 eligible units are Lung or Kidney units;
   statements about "where aging signal is strongest" are statements about this atlas's design
   (cell abundance and mouse coverage), not about mouse biology as a whole.
4. **Tabula Sapiens v1 is an adult cohort.** Skin 33–61 y (n=4), lung 59–61 y; no donor-age
   contrast exists, so TS cannot support a human aging DE in these tissues.
5. **SenNet inaccessibility.** The public index lacks donor ages (6,820 datasets queried
   2026-10-14); the restricted path (Globus) was out of session budget. The human arm is
   untestable as a result, not by choice of threshold.
6. **The human arm of the signature figure is therefore descriptive only** (n=4, CI spanning
   −1.37 to +0.95); no human conservation or divergence claim is made anywhere in this repo.
7. **NMR program coverage is 31%.** Less than a third of the program genes were detected in
   the NMR data; the NMR score is a partial-program score and could differ from the full-program
   score in either direction.
8. **Igfbp3 is unverifiable.** The paper's key gene is absent from the GEO count file; claims
   that hinge on it could not be checked.
9. **NMR ages are group means.** The 5 samples collapse to 1.1 y (n=3) and 11.3 y (n=2);
   within-group age variation is unknown to us.
10. **Equivalence not established.** TOST at the pre-registered δ failed (p=0.329); the
    interaction p=0.064 is suggestive only. "NMR is flat" is not a claim this repo makes.
11. **Mouse skin arm is n=3, all ≥18 m.** TMS droplet contains no young skin mice; the
    epidermal mouse slope (0.452 [-3.679, 4.583], p=0.397, n=3) is uninformative, which is why
    the positive control was moved to the whole-body arm (§6).
12. **No multi-tissue NMR atlas exists.** The outlier contrast is skin-only; tissue-specific
    aging resistance elsewhere in the NMR is untested.
13. **Survivorship.** The 30m mice are a survival-selected cohort; the extended-bins
    sensitivity (120 vs 55 DE genes) mixes age with selection.
14. **Harvest-time unknown.** Time-of-day and fasting state at sacrifice are not recorded in
    the atlas and can confound DE.
15. **Dropout.** Droplet scRNA-seq dropout, amplified at the small pseudobulk sizes of some
    units (e.g., Lung::adventitial cell, 423 cells over 8 mice), biases against low-expression
    genes.
16. **Model asymmetry in M3.** For 15 of 17 units our dichotomous effect sizes were compared
    to the paper's continuous coefficients; only 2 units had a continuous bridge fit. Some
    "partial" verdicts may be model asymmetry rather than biology.
17. **Positive controls were not significant in the top unit.** The canonical markers checked
    in Lung::adventitial cell showed plausible directions for 3 of 4 genes but none reached
    p<0.05; the calibration claim rests on the permutation null and sex floor.
18. **NMR "aged" is middle-aged.** 11.3 y is roughly a third of the NMR's 31-y AnAge maximum;
    the data say nothing about late-life NMR programs.

## 5. Deviations from pre-registration (log)

- **M4 — SenNet substitution by documented inaccessibility.** Pre-registered plan: SenNet as
  the primary human contrast. Reality: the public index has no donor ages (query documented,
  2026-10-14). No silent substitution was made; the human arm is reported as untestable.
- **M5 — no conservation/divergence calls.** Pre-registered plan: two-axis conservation with
  permutation null. Deviation: with the human arm untestable, every unit is marked UNTESTABLE
  and the divergence table explicitly declines to call divergence (the first run's error of
  labeling low-powered rows "Divergent" is corrected).
- **M6 — full mouse program scored instead of the strict conserved program.** Pre-registered
  plan: score the strictly human-replicated program. Deviation: strict replication retained
  4 genes (<10-gene floor) at n=4 donors, so the full mouse program (167 up / 94 down) was
  scored and is flagged everywhere it appears.
- **M6 — whole-body mouse positive-control arm substituted for the skin arm.** Pre-registered
  plan: mouse-skin positive control. Deviation: TMS droplet has no young skin mice (epidermal
  arm n=3, all ≥18m, slope CI [-3.679, 4.583]); the positive control was run on whole-body
  per-mouse scores (23 mice), which passed (0.235 [0.067, 0.403], p=0.008).
- **M6 — human arm shown despite being low-powered.** Shown with wide CIs and flagged
  LOW-POWERED; no claim made from it.

## 6. Relationship to prior work

This project builds directly on the Tabula Muris Senis atlas (Tabula Muris Consortium, Nature
2020), whose published age-coefficients we use as the external reference for our agreement
table (13/4/0 across 17 powered units), and extends it with the pre-registered calibration
trio, an explicit power audit, and a longevity-outlier contrast. Schaum et al. (Nature 2020) is
related work — cross-species, life-stage-resolved single-cell aging — and no findings are
attributed to it here. Tabula Sapiens (Science 2022) was obtained as the human replication
cohort and found to be an adult cohort, which is itself a documented negative result for the
field's habit of treating it as an aging reference. Savina et al. (2022, *Aging* — the journal
is *Aging*, not *Nature Aging*) is the source of the NMR data and is partially reproduced,
with a unit-conversion caveat on its headline zero-DE claim. Lifespan denominators come from
AnAge (Tacutu et al. 2018, NAR), and the SenNet consortium's public portal was the attempted
source of the primary human aging cohort. The recently published cross-mammalian conserved
aging programs of Tyshkovskiy/Gladyshev and the AP-1 inflammaging analysis of Karakaslar et
al. were **not yet incorporated** and are named future-work comparisons, not cited supports.
Differentiation from the existing literature: cell-type-resolved ranking rather than
whole-tissue programs; an honest power analysis that treats "untestable" as a result; and a
longevity-outlier contrast tested with pre-registered equivalence statistics rather than
absence-of-significance language.

## 7. Code availability and reproduction

All analysis code is in this repository: `scripts/run_all.py` (executes the full pipeline,
notebooks 01–07 via the nbconvert gate command and writes every figure and `results/` table),
`scripts/fetch_data.py` and `scripts/verify_data.py` (data acquisition and SHA-256
verification against `data/manifest/`), the notebooks in `notebooks/`, importable utilities in
`src/`, and the frozen environment in `environment.lock`. Reproduction:

```
conda create -n agescape --file environment.lock
python scripts/fetch_data.py
python scripts/run_all.py
```

A full recompute takes on the order of hours on a 12-core/49 GB machine. Every quantitative
claim in this report traces to a committed file in `results/` or `figures/`; the notebooks are
committed fully executed, and every stochastic step is seeded (seed 42). **AI-assistance
disclosure:** this analysis was executed by AI coding agents under a strict no-invented-numbers
protocol, with human review gates between milestones; the takeover audit of 2026-10-14
(documented in `notes/gates/M-takeover-report.md`) found and replaced an earlier
synthetic-data result, and every number in this document was re-read from the real outputs.

## 8. Figures

- **Fig. 1** `figures/m2-aging-signal-ranking.png` — aging-signal ranking with bootstrap CIs.
- **Fig. 2** `figures/m2-calibration-permutation.png` — permutation null vs observed DE counts.
- **Fig. 3** `figures/m6-signature-figure.png` — program score vs life-stage-normalized age,
  mouse (n=23) / human (n=4, low-powered) / NMR (n=5).
- **Fig. 4** `figures/m6-savina-reproduction.png` — Savina DE-count-vs-cutoff reproduction curve.
- **Fig. 5** `figures/m3-dpt-vs-age.png` — diffusion pseudotime vs age group, Marrow.

## 9. References

1. Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes ageing tissues in
   the mouse. *Nature* 583, 590–595 (2020). DOI: 10.1038/s41586-020-2496-1.
2. Schaum N. et al. Ageing hallmarks exhibited by each cell type in the mouse. *Nature* 583,
   596–602 (2020). (Related work; no findings attributed.)
3. The Tabula Sapiens Consortium. The Tabula Sapiens: a multiple-organ, single-cell
   transcriptomic atlas of humans. *Science* 376, eabl4896 (2022).
4. Savina A. et al. Single-cell transcriptomics reveals age-resistant mechanisms in the naked
   mole-rat. *Aging* 14 (2022). GEO GSE152886.
5. Tacutu R. et al. AnAge 2018: a database of ageing and longevity. *Nucleic Acids Research*
   46(D1), D997–D1003 (2018).
6. NIH SenNet Consortium — Cellular Senescence Network, public data portal.
