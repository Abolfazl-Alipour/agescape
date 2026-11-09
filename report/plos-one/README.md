# PLOS ONE submission version, Rev2 stage 1

Manuscript: "Cell-type-specific compositional and transcriptional changes in mouse and human aging", sole author Abolfazl Alipour.

- `final-plos-one.tex` — LaTeX source, standard IMRaD order (Materials and Methods after the Introduction). Abstract 230 words, no statistical results.
- `final-plos-one.pdf` — compiled PDF (Tectonic).
- `figures-v2/` — data figures 1 to 4. Figure 3 was rebuilt in Rev2 stage 1 from the joint donor-level concordance null; see `results/reviewer-rev2-concordance-figure-data.csv`.

Status: draft, NOT submitted anywhere. Affiliation is blank and the Funding, Competing Interests, and Ethics Statement sections are marked pending author confirmation. The older manuscript at `report/manuscript-aging-mechanisms.md` is untouched.

Target journal: PLOS ONE (impact factor 2.8, soundness-only review), chosen 2026-11-08. The cover letter draft is not in the repo; it lives with the working copy at `~/workspace/agescape-final/plos-one/cover-letter-draft.md`.

## What changed in Rev2 stage 1

Headline is now the middle claim only: no uniform transcriptomic program was detected in the surveyed datasets, with upstream causes not settled and shared mechanisms with cell-type-specific consequences not excluded. The Introduction adds an operational definition of "universal" (same direction in every cell type, overlapping gene sets, comparable lifespan-scaled magnitude, present in each species tested).

- Methylation clock (H1, Universal Clock 2) deleted entirely, including its Methods subsection, Results paragraph, Discussion and Limitations mentions, and references. Reference list renumbered to 24 entries, all cited in text.
- Proliferation analysis deleted from the main text entirely.
- Off-program conclusion rewritten consistently in the abstract, Results, Discussion, and Conclusions: no statistically significant increase was detected in mouse lung at the tested depths, but confidence intervals were wide and biologically relevant increases were not excluded. Scoped to mouse lung only.
- Results updated with the verified Rev2 reanalysis numbers: joint donor-level concordance null (mouse within-tissue plus 0.187 vs null plus 0.138, p = 0.043; cross-tissue not significant); vector depth survival with dispersion dead at every fixed depth; the 18-month reference arm showing drift concentrated in the young-to-mature transition; human vector interval reported honestly with no equivalence claim; exact sign test replaced by a donor-cluster bootstrap; human lung club and multiciliated findings downgraded to between-study associations; signature method sensitivity stated (moderated WLS 151 calls vs DESeq2 312 vs unmoderated WLS 1; 310 of 765 P1 instances flip sign on leave-one-donor-out).
- Theory cut to a short predictions-only section; no verdicts on whole theories. Only the same-cell multiome test is kept as the next step.
- Meta-analysis removed from the main text; kept as a short text note in Supporting Information.
- Figure 5 (AI-generated illustration) and the schematic summary subsection removed. Supporting Figures S1, S2, and S3 removed: S1 was AI-generated, S2 lists the deleted proliferation and clock arms and gives theory verdicts, and S3 states the human vector null is evidence of absence, a framing now dropped.
- Appendix: Evidence and Caveats deleted. All process residue removed (reviewer-response language, phase labels, commit hashes, file paths, submission status). Acknowledgements dropped. A cohort table (Table 1) was added to Methods; cells that could not be verified are marked "not reported". The Bayesian Methods subsection was deleted (no Bayesian results paragraph added). Figure 1 caption now states its selection rule (25 results passing false discovery rate correction within their own strata at pseudocount 0.5, sorted by slope; 21 of the 25 also pass the wider global correction).
- Data Availability gives the repository URL (https://github.com/Abolfazl-Alipour/agescape) and notes a Zenodo archive is pending with no DOI claimed.

## Pending refresh (stage 2)

- Transport matched random null in MSSM still running. When it lands, update the signature section verdict on the three survivors (MSSM L2/3 minus 0.001303 q = 0.0208; MSSM oligodendrocytes plus 0.001152 q = 0.0346; PsychAD EN_L3_5_IT_3 minus 0.001383 q = 0.0208) and the Discussion validation paragraph.
- Human interneuron genome-wide moderated redo pending. When it lands, refresh the interneuron section baseline numbers (MSSM FDR VIP 294, pvalb 333, sst 123; both-cohort FDR 24, 20, 62; effect correlations r 0.89, 0.86, 0.93; same-direction 272 of 294, 295 of 333, 118 of 123) and Figure 4 if the gene sets change.

See `~/workspace/agescape-final/plos-one/REV2_STAGE1_NOTES.md` for the full per-fix change list, refresh list, and figure decisions.
