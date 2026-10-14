# Methodological Limits — Naked Mole-Rat (NMR) Comparative Analysis

Every limit below is stated with the file it is read from. This document replaces the
first-run version, which described a cohort design (cultured fibroblasts, a δ = 0.15
TOST "success," a p = 2.2e-16 species×age interaction) that does not correspond to any
executed analysis.

## 1. Tissue and cohort coverage

- **Epidermis only.** The only obtainable NMR single-cell dataset is GSE152886
  (Savina et al. 2022, *Aging*): epidermal cells, 5 samples after the paper's own 2V
  exclusion (1J/5J/6J young, 3V/4V aged). There is **no multi-tissue NMR atlas**; the
  program score cannot be tested beyond skin.
- **Two age groups as group means.** The paper reports group-mean ages only — young
  1.1y, aged 11.3y; individual ages are not in GEO. Normalized to the AnAge maximum
  (31y), the observed span is **≈0.03–0.36 of maximum lifespan** — a two-point
  contrast, not a lifespan continuum (`results/m6-lifespan-program-scores.csv`).
- **Cross-sectional**, not longitudinal: young and aged animals are different
  individuals; animal-level n = 5 with 3 young / 2 aged.

## 2. Program-gene coverage in the NMR matrix

- **31% coverage.** Only 0.314 of the 257-gene mouse program has detectable NMR
  orthologs in the count matrix (`coverage_nmr` gate, `results/m6-gates.csv`) — a
  consequence of HetGla assembly/annotation gaps. The program score is computed on
  the surviving ~31% of genes; the slope and equivalence tests inherit that sparsity.
- **Human arm coverage** (for the TS-skin replication context) is 0.732
  (`coverage_human` gate) — better, but the human arm has only n = 4 donors.

## 3. What could NOT be verified

- **Igfbp3 is unverifiable.** The paper's key within-cluster finding centers on Igfbp3,
  whose Ensembl ID `ENSHGLG00000037774` is **absent from the GEO count file** (checked
  directly against `data/raw/GSE152886_nmr_skin_raw.csv.gz` in Notebook 06). The
  paper's central gene-level claim cannot be confirmed or refuted with obtainable data.
- **The paper's "no DE genes" result was only partially reproduced.** On our QC object
  (10,251 cells across the 5 post-filter samples), the DE-count-vs-cutoff curve
  (`results/m6-savina-reproduction.csv`) gives **29 genes at the paper's stated rule**
  (p < 0.05 and |log2FC| ≥ 1.0) vs the paper's 0; 560 genes at our standard rule
  (FDR < 0.05, |log2FC| ≥ 0.25). **Caveat:** the paper's fold-changes are reported on
  the Seurat natural-log scale; converting thresholds between ln and log2 shifts the
  comparison. We flag the discrepancy rather than adjudicate it.

## 4. Statistical interpretation of "flat"

- **Equivalence is NOT established.** The pre-registered TOST margin was δ = 25% of
  the mouse whole-body slope (δ = 0.0587); TOST gives **p = 0.329**
  (`results/m6-statistical-summary.csv`). The exact sample-level permutation test on
  the aged-vs-young score difference gives p = 0.818. The mouse-vs-NMR slope
  difference is borderline (0.218 [−0.013, 0.449], p = 0.064).
- **What we can say:** with 5 samples over ~1/3 of normalized lifespan and 31% program
  coverage, there is **no detectable change** in the program score with age in NMR
  epidermis. **What we cannot say:** that NMR aging is "proven flat" — the data are
  underpowered to establish equivalence, and the human/mouse comparison arms carry
  their own power limits (mouse skin n = 3; human n = 4).

## 5. Interaction with the rest of the project

- The NMR contrast is anchored to the **mouse program** (M2/M5). Because strict human
  replication in TS skin (n = 4) intersected only 4 program genes, the full mouse
  program was scored (deviation D3, `notes/methods.md` §7) — the NMR result is
  therefore a statement about the *mouse* aging program in NMR epidermis, not about a
  three-species conserved program.
- Lifespan denominators are AnAge species maxima (mouse 4.0y, NMR 31y, human 122.5y);
  normalized ages are only as good as those denominators (`results/m6-species-lifespans.csv`).
