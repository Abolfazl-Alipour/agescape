# Reviewer response Phase 3: manuscript package inputs

## Scope and provenance

Phase 3 started from commit `868ca108d8e14ed89ceca73c77b7086dfd54faf9`. It produced package inputs only. The repository manuscript was not edited. No submission package was rebuilt, and no submission, publication, or external contact occurred.

## 1. Figures and figure data

| Figure | Files | Underlying data | Caption |
|---|---|---|---|
| Mouse vector drift by cell type | `figures/reviewer-phase3-mouse-vector-drift.png` and `.pdf` | `results/reviewer-phase3-mouse-drift-figure-data.csv` | `reports/reviewer_phase3_figure_captions.md` |
| Interneuron gene changes | `figures/reviewer-phase3-interneuron-gene-changes.png` and `.pdf` | `results/reviewer-phase3-interneuron-figure-data.csv` | Same caption file |
| Composition changes with CIs | `figures/reviewer-phase3-composition-changes.png` and `.pdf` | `results/reviewer-phase3-composition-figure-data.csv` | Same caption file |
| Gene-effect concordance | `figures/reviewer-phase3-universality-concordance.png` and `.pdf` | `results/reviewer-phase3-universality-figure-data.csv` | Same caption file |

All four images were generated and inspected. The captions state what each estimate measures and what it does not establish.

## 2. References

`reports/reviewer_phase3_references.md` contains 32 references verified against Crossref metadata with title similarity at least 0.85. The full audit, including requested and matched titles, is `results/reviewer-phase3-reference-metadata.csv`.

The existing Zhang citation, `Characterization of the conserved mammalian aging signature through single-cell cross-species meta-analysis`, could not be verified. It is excluded from the verified list and marked for author correction. No DOI was invented for it.

## 3. Methods draft

`reports/reviewer_phase3_methods_draft.md` contains the complete author-facing Methods draft: dataset accessions and inferential units, normalization, regression models, composition transformations, young-reference construction, permutations, cell-type harmonization, Bayesian priors, depth sensitivity, annotation stability, matched signature nulls, and clock benchmarking.

## 4. Residue removal

`reports/reviewer_phase3_residue_and_scope.md` lists exact repository-manuscript locations for repository and explorer links and records whether Source lines, commit hashes, test counts, author-decision appendices, and draft banners were found. The external draft was unavailable, so the author must repeat the residue check against the actual submission file.

## 5. Scope and title recommendation

The main claims should use mouse and human. Fly, worm, zebrafish, NMR, proliferation moderation, and clock benchmarking should be exploratory or supplementary for the design reasons documented in the report.

Recommended title: **Cell-type-specific compositional and transcriptional changes in mouse and human aging**.

The author decides the final title and scope.

## Exact commands

Commands ran from `/home/reza/projects/agescape` using the core environment unless noted:

```bash
/home/reza/miniconda3/envs/agescape/bin/python scripts/run_reviewer_phase3_figures.py
/home/reza/miniconda3/envs/agescape/bin/python scripts/build_reviewer_phase3_references.py
```

Reference verification used read-only Crossref metadata requests. Earlier source-specific verification used the public atlas, PubMed, publisher, and GEO records cited in the metadata and Phase 2 report.

## Stop boundary

Phase 3 verification gave 68 passing tests and one Bayesian-only skip in the core environment; all 146 app-data anchors passed. Figure data, captions, references, compilation, and `git diff --check` passed. The author writes the manuscript revision. No manuscript edit was made by this analysis.
