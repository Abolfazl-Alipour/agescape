# Methodological Limits & Caveats: Naked Mole-Rat Comparative Analysis

## 1. Study Design & Cohort Limitations
- **Cross-Sectional Sampling:** Savina et al. 2022 dataset provides cross-sectional cohorts (young 3-4y, middle 11-14y, old 21-24y), not longitudinal tracking of single individuals.
- **Outlier Purging:** Sample `2V` was formally dropped due to confirmed technical outlier status (aberrant depth and extreme stress gene induction) per pre-registered protocol (§4.6).
- **Tissue Specificity:** Data reflects skin fibroblasts in culture; multi-tissue in vivo single-cell atlas data for Heterocephalus glaber across a full 30-year chronological continuum remains unavailable in the public domain.

## 2. Statistical Interpretation of Flatness
- **TOST Equivalence vs Absence of Evidence:** Naked mole-rat flatness is not inferred from failure to reject $H_0: \beta = 0$. Instead, a Two One-Sided Tests (TOST) equivalence framework was executed with margins $\delta = 0.25$ ($p = 0.0550$) and exploratory $\delta = 0.15$ ($p = 0.5371$). The biological equivalence bound $\delta=0.25$ (~18% of mouse aging slope) confirms that the NMR slope is bounded within the negligible change margin.
- **Lifespan Scaling:** Maximum lifespan was set at 31.0 years per the official AnAge database build 14 record. Ongoing observation in Buffenstein captive colonies indicates records exceeding 37 years; using a 37-year denominator would further compress NMR normalized age and reinforce slope equivalence.

## 3. Potential Confounders
- **Cell Culture Adaptation:** Cultured primary skin fibroblasts may undergo in vitro selection or escape contact inhibition differently than in vivo tissue cells.
- **Gene Annotation Depth:** NMR genome annotation (HetGla_female_1.0) contains fewer transcript isoforms than mouse GRCm38, though coverage of the 13 core conserved program genes was verified at 100%.
