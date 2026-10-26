# Crucible Log — 10 Rounds (5 perspectives × 2 rounds)

Every perspective audited the plan against live files, data and endpoints. All findings below were folded into Plan v3.

| Round | Perspective | Verdict on draft | Headline findings (all fixed in v3) |
|---|---|---|---|
| 1, 2 | Biostatistician | Reject v1; v2 had 3 blockers | `immune-endo-depleted` 3m male Lung pools; PyDESeq2 `independent_filter` flips the #1 unit (Kidney PCT vs Lung adventitial); 1-cell pseudobulks; exact-permutation floors; zero-power TOST; 4 M6 code bugs; Q* additivity bug (908/4815 negative between-organ values); CR2 absent in statsmodels; 51-level subtissue crash; floor survival 17→12→9; 15 / 6 valid global permutations |
| 3, 4 | Single-cell biology | Reject v1; v2 had 3 blockers | Fabricated `orthologs-mouse-nmr.csv`; Igfbp3 IS in GSE152886 (`ENSHGLG00000018612`); Ins1/Ins2/Rn45s soup; `Gm2a` regex false positive; Mgp/Col3a1 ambient in depleted Lung; M9 missing CP10K; M10 Brain_bank × Age collapse; celltype-mapping.md fails 11 of 17 units; Ransick table unobtainable (HTTP 403) |
| 5, 6 | Reproducibility / CI | Reject v1; v2 had 5 blockers | Wave 1 race conditions; dead BioMart endpoint; manifest placement crashes fetch/verify; notebook 07 legacy bomb; `set_seeds` thread vars are a no-op after imports; `/tmp` and `/home/reza` path lint |
| 7, 8 | Security / DUA | Reject v1; v2 had 4 blockers | Tracked per-donor DUA scatter PNG; `Developmental` stratum singular matrix; CI test contradicts audit schema; unconditional richText tooltip breaks 5 charts; favicon 404 / CSP |
| 9, 10 | Hiring lead and subagent architect | Reject v1; v2 had 3 blockers | Split-personality title; G0–G16 primary-vs-audit contract; G9 geometry; G16 redesign; ECharts `cssText` vs `style-src`; Worker E overload |

**Unresolved caveat:** the Single-Cell Biology subagent hit a quota error (429, resets in about 167h) right after delivering its Round 2 report. Its review is complete. Subagent-based execution may hit the same quota, so see the fallback in the plan.

**Ground-truth items verified live:**
- Zhang 2021 *eLife* DOI `10.7554/eLife.62293`; supplementary 3 gives 330 FACS GAGs and 59 droplet GAGs.
- SenMayo: 125 genes in human, 117 in the mouse panel.
- Savina 2022 *Aging* DOI `10.18632/aging.204054`.
- Karakaslar 2023 *Aging Cell* DOI `10.1111/acel.13792`.
- Tyshkovskiy 2026 *Nature* DOI `10.1038/s41586-026-10542-3`.
- Yang 2026 *Nature* DOI `10.1038/s41586-026-10271-7`.
- Ensembl 106 `HetGla_female_1.0`: Igfbp3 maps to `ENSHGLG00000018612`, with 29,793 UMIs in 5,841 cells.
- Ensembl 106 FTP URLs for GTF and Compara homologies.
