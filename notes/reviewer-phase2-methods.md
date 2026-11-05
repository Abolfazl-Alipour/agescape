# Reviewer-response Phase 2 methods

## Scope

Phase 2 performs new runs on existing datasets. No new single-cell dataset was collected. The manuscript was not edited. Outputs use the `reviewer-phase2-` prefix.

## Entropy and off-program depth sweep

Mouse lung singlet cells were evaluated at 500, 1,500, and 3,000 UMIs per cell and at full depth. Each fixed-depth arm included cells with at least that many UMIs and used without-replacement downsampling with random seed 42. The full-depth arm used no downsampling and retained the original 1,500-UMI inclusion threshold.

Young-reference off-program marker definitions were fixed across depths using the original young cells with at least 1,500 UMIs. References remained cross-fitted by animal. Off-program expression was the fraction of actual total UMIs assigned to those markers. Entropy was calculated over the same retained gene classes as the original analysis.

Animal means were the inferential units, with at least 20 cells per animal-cell-type unit. Models included age and sex. The library-size arm added the animal's mean log actual pre-downsampling UMI total. Global models standardized outcomes within cell type, used cell-type intercepts, and clustered standard errors by animal. Effects are reported per 15 months with t-based intervals. The fixed-depth 1,500 effect is expected to reproduce the original point estimate; its t-based interval and p-value can differ slightly from the earlier normal-approximation global output.

The verdict rule was fixed before interpretation: a significant positive full-depth off-program result would require withdrawal of the downsampling-null interpretation. A null full-depth result supports the null but does not establish equivalence. The separate direct vector arm was rerun without downsampling using `scripts/run_h3_expression_vector.py`.

## Human adult and nonlinear models

The public human vector and composition pipelines already restrict donors to age 20 years or older. Adult-only status and age distributions were verified from the donor-level inputs. PsychAD developmental donors are not part of these public vector and composition inputs.

Nonlinear age models used a cubic B-spline with three degrees of freedom, fixed before examining outcomes. Linear and nonlinear models were fitted to the same donor-cell-type units. Vector models included sex and mean log UMI; composition models used the existing CLR outcome with sex when variable. F tests compared nonlinear and linear fits. BH correction was applied within each cohort separately for linear and nonlinear p-values.

## Annotation stability

No per-cell annotation-confidence probability was available in the TMS or public human inputs. Doublet scores were not relabeled as annotation confidence.

Excluded-cell and unassigned-cell counts were tabulated by age and tissue or cohort. Mouse exclusions compare the raw droplet object with the cleaned object. Human exclusions use primary normal-disease cells and the existing 1,500-UMI threshold. Unassigned labels include unknown, missing, or equivalent labels.

Marker stability used a fixed canonical marker panel selected by cell identity. Donor-cell-type marker expression was aggregated as log CP10K and regressed on age, sex, and mean log UMI. At least six donors and 20 cells per donor-cell-type unit were required. BH correction was applied across the targeted marker tests. These tests can detect marker change but cannot prove annotation identity loss or exclude every annotation error.

## Clock overlap

GSE199979 metadata and the Universal Clock publication were checked. The publisher's complete Supplementary Data 1 to 14 workbook was downloaded from `https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00462-6/MediaObjects/43587_2023_462_MOESM3_ESM.xlsx` and searched for the accession, a tested GSM identifier, and a tested BXD sample title. Absence of these identifiers leaves sample-level membership unresolved. GSE184211 remains a consortium-development dataset as recorded in `results/h1-training-overlap-audit.csv`.

## Reproduction

```bash
python scripts/run_reviewer_phase2_depth.py
python scripts/run_h3_expression_vector.py
python scripts/run_reviewer_phase2_human_age.py
python scripts/run_reviewer_phase2_annotation.py
python scripts/run_reviewer_phase2_clock.py
```
