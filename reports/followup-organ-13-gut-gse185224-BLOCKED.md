# Follow-up Blocker Report: Healthy Adult Intestinal Epithelium (`GSE185224`)

**Status**: BLOCKED / NOT ANALYZED (`notes/followup-dataset-checklist.md` #34: `SKIP`; severe sample-size underpowering `n=3`, `residual_df = 0` with `log_lib`, narrow `29-53` yr adult span)
**Primary accession**: [`GSE185224`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185224) (`GSM5655495` to `GSM5655500`, `PRJNA768275`)
**Primary publication**: Burclaff et al. (2022), *Cellular and Molecular Gastroenterology and Hepatology*, "A Proximal-to-Distal Survey of Healthy Adult Human Small Intestine and Colon Epithelium by Single-Cell Transcriptomics", `DOI: 10.1016/j.jcmgh.2022.02.007`, `PMID: 35176508`, `PMCID: PMC9043569`

## 1. Primary-Source Metadata Verification

- **Accession resolved**: Previously flagged in [`followup-dataset-checklist.md`](../notes/followup-dataset-checklist.md) (#34) with unresolved accession, `GSE185224` (`BioProject PRJNA768275`) is now verified as the open GEO series for Burclaff et al. (2022), containing 10x Chromium v3.1 gene-expression (`GEX`) and cell-surface oligo hashtag (`CSP`, `TotalSeq-B0251` to `B0256`) libraries for `12,590` filtered intestinal epithelial cells across 6 gut regions (`Duodenum`, `Jejunum`, `Ileum`, `Ascending Colon`, `Transverse Colon`, `Descending Colon`).
- **Donors, ages, sex, and disease/metabolic notes**: **`n = 3` brain-dead male organ donors (`0` females)** recruited under donor acceptance criteria `age <= 65` years, with verified chronological ages **`29`, `45`, and `53` years**. Per `PMC9043569`, `Donor 3` differs from the other two donors on BMI, ancestry, and comorbidity (`only donor with type II diabetes`).

| Cohort Attribute | Verified Value | Primary Source |
|---|---|---|
| Total biological donors | `3` | `PMC9043569` / `GSE185224` SOFT |
| Sex | `3` male, `0` female (`100%` male) | `PMC9043569` |
| Exact chronological ages | `29`, `45`, `53` yr (`0` donors `>= 60` yr) | `PMC9043569` |
| Anatomical coverage | `6` intestinal regions hashed per donor (`12,590` cells total) | `GSE185224` SOFT / `PMC9043569` |
| Assay | 10x Chromium Next GEM Single Cell 3' v3.1 + TotalSeq-B hashtags | `PMC9043569` / `GSE185224` |

## 2. Documented Blockers (Why No Expression Matrices Were Downloaded or Modeled)

1. **Zero degrees of freedom with standard covariates (`n = 3`)**: At the biological donor unit (`n = 3`), linear modeling with `~ intercept + age + log_lib` consumes all `3` parameters, leaving **`residual_df = 0`** (even an unadjusted slope `~ intercept + age` has only `residual_df = 1`). Multiple gut regions (`n = 6` hashed regions per donor) represent correlated within-donor anatomical replicates, not independent aging subjects.
2. **Non-aging adult span (`29-53` yr) and single-sex confounding**: All `3` donors are young-to-middle-aged males (`29, 45, 53`), lacking any elderly (`>= 60` yr) or female representation, while `1` of the `3` donors (`Donor 3`) carries type II diabetes alongside distinct BMI and ancestry.
3. **Conclusion**: `GSE185224` cannot support donor-level aging inference (`Family A`, `Family B`, or `Family C`) under Agescape's `residual_df >= 3` requirement and remains excluded as documented underpowered data.
