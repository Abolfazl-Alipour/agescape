# Follow-up Blocker Report: Human Liver Single-Nucleus Atlas (`GSE210077`)

**Status**: BLOCKED / NOT ANALYZED (`notes/followup-dataset-checklist.md` #24: `BLOCKED`; primary IRB/public release restricted to 10-year age brackets, no exact donor ages publicly available)
**Accession**: [`GSE210077`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210077) (`PRJNA863623`, 40 GSM samples `GSM6416542` to `GSM9323722`), UCSC Cell Browser [`liver-sex-age-map`](https://liver-sex-age-map.cells.ucsc.edu), Dryad `doi:10.5061/dryad.37pvmcvsg`
**Primary publications**:
- Watson et al., *Nature Communications* (2025), "Spatial transcriptomics of healthy and fibrotic human liver at single-cell resolution", `PMID: 39747812`, `PMCID: PMC11697218`
- Rahman et al., *JHEP Reports* (2026), "Single-cell transcriptomics reveals the impact of sex and age in the healthy human liver", `DOI: 10.1016/j.jhepr.2026.101773`, `PMID: 41990579`, `PMCID: PMC13101288`
**Audit table**: [`followup-organ-14-liver-gse210077-samples.csv`](../results/followup-organ-14-liver-gse210077-samples.csv)

## 1. Primary-Source Metadata Verification

- **Assay and tissue**: 10x Genomics Chromium Single Cell 3' v3.1 single-nucleus RNA-seq (`snRNA-seq`, `GPL18573`, NextSeq 500 / NovaSeq) on flash-frozen excess surgical liver tissue from Massachusetts General Hospital (`MGB` / `DFCI` IRB protocols).
- **GEO sample inventory (`GSE210077`, `n=40` GSMs in [`samples.csv`](../results/followup-organ-14-liver-gse210077-samples.csv))**:
  - `22` labeled `disease state: Healthy` (`GSM6416542` to `GSM6416592`, `Liver-1` to `Liver-24` excluding `Liver-12` and `Liver-21`)
  - `14` labeled `disease state: Normal` (`GSM8493737` to `GSM8493743`, `GSM8493747` to `GSM8493748`, `GSM9323718` to `GSM9323722`)
  - `3` labeled fibrotic disease (`GSM8493744` `Liver-32` `F4 fibrosis` `M` `51-60`; `GSM8493745` `Liver-33` `F2 fibrosis` `F` `41-50`; `GSM8493746` `Liver-34` `F3 fibrosis` `M` `61-70`)
  - `1` sample (`GSM7156103`, `Liver-21`) has only `tissue: Liver` in GEO SOFT (missing `age`, `Sex`, and `disease state` fields).
- **UCSC Cell Browser (`liver-sex-age-map/all-samples`, `195,515` nuclei, `37` donors)**: Includes the `36` GEO `Healthy`/`Normal` samples plus `Liver-21` (`5,623` nuclei), yielding `20` female (`113,271` nuclei) and `17` male (`82,244` nuclei) donors and excluding the `3` fibrotic samples (`Liver-32`, `Liver-33`, `Liver-34`).

## 2. Documented Blockers (Why No Modeling Was Run)

1. **Exact donor ages do not exist publicly due to IRB de-identification (decade brackets only)**:
   - Watson et al. (2025, `PMC11697218`, Methods) explicitly state that samples were collected under Mass General Brigham / Dana Farber Cancer Institute IRB protocols *"including approval to provide data on age by decade and sex"*.
   - Rahman et al. (2026, `PMC13101288`, Discussion/Limitations) explicitly state that *"the assignment of age by decade to each sample prevented exact age matching between donors"*.
   - Both GEO `GSE210077` and UCSC `liver-sex-age-map` expose only 8 overlapping 10-year brackets across the 37 non-fibrotic donors: `18-30` (`25,860` nuclei, 1 GEO sample), `21-30` (`9,551` nuclei, 5 GEO), `31-40` (`15,140` nuclei, 3 GEO), `41-50` (`21,535` nuclei, 5 healthy GEO), `51-60` (`40,563` nuclei, 7 healthy GEO), `61-70` (`45,166` nuclei, 8 healthy GEO), `71-80` (`26,351` nuclei, 5 GEO), and `81-90` (`11,349` nuclei, 2 GEO), plus `GSM7156103` (`Liver-21`) missing an age attribute in GEO.
2. **Continuous aging slope requires exact donor ages (or a prespecified categorical/ordinal exception like GTEx v8)**: Replacing 10-year age brackets (`18-30` through `81-90`, spanning 10 to 13 years per bin with overlapping `18-30` and `21-30` ranges) with fabricated bin midpoints would violate primary-source verification rules and introduce differential interval error across decades.
3. **Additional surgical/cohort caveats**: Rahman et al. (`PMC13101288`) note that samples derive from single-center excess surgical resections (`35/37` Caucasian), mixed `Healthy`/`Normal` labels plus `3` F2-F4 fibrotic samples in the same GEO series, and `GSM7156103` (`Liver-21`) lacks GEO disease/sex/age annotations.

## 3. Verification Conclusion

Blocked after full primary-source inspection (GEO `GSE210077`, UCSC `liver-sex-age-map`, `PMC11697218`, `PMC13101288`). Zero `GSE210077` expression matrices (`GSE210077_RAW.tar`) were downloaded.
