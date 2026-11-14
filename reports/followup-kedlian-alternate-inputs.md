# Kedlian alternate-input schema check

Date: 2026-11-14. Scope: bounded HTTP range/schema probes only; no whole-object download, expression analysis, or biological outcomes.

## Conclusion

**No equivalent raw-count object was found among the official processed alternatives.** The 4.03-GB raw H5AD remains required for the frozen raw-count pipeline. A complete copy of the exact browser-download file is needed at:

- URL: `https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`
- Destination: `/home/reza/projects/agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`

The existing partial file is not usable. The download was not resumed.

## Probe findings

Every object below had a float32 sparse `X`, **no `raw` group, and no layers**. Up to 3,072 values per object (small chunks at the start, midpoint, and end) were finite and nonnegative but **not integer-valued**. These are tiny value probes, not full count validation.

| Official object | Rows x genes | Count status | Annotation/sample scope |
|---|---:|---|---|
| `SKM_human_pp_cells2nuclei_2023-06-22.h5ad` | 183,161 x 29,400 | No raw counts; noninteger `X` sample | Complete cells+nuclei object; 17 donors and all 43 SampleIDs |
| `SKM_MuSC_human_2023-06-22.h5ad` | 17,528 x 19,668 | No raw counts; noninteger `X` sample | Targeted scRNA MuSC subset |
| `SKM_myonuclei2myofiber_human_2023-06-22.h5ad` | 87,522 x 26,689 | No raw counts; noninteger `X` sample | Targeted myofiber cells+nuclei |
| `SKM_myonuclei_human_2024-04-14.h5ad` | 74,626 x 25,900 | No raw counts; noninteger `X` sample | Targeted snRNA myonuclei; overlaps the preceding object |
| `SKM_immune_human_2023-06-22.h5ad` | 20,321 x 20,220 | No raw counts; noninteger `X` sample | Targeted scRNA immune subset |
| `SKM_fibroblasts_Schwann_human_2023-06-22.h5ad` | 20,611 x 20,772 | No raw counts; noninteger `X` sample | Targeted scRNA fibroblast/Schwann subset |
| `SKM_Endothelium_SMC_human_2023-06-22.h5ad` | 19,317 x 20,221 | No raw counts; noninteger `X` sample | Targeted scRNA endothelial/mural subset |

## Metadata and annotation sufficiency

The full processed object has `SampleID`, `DonorID`, `Age_group`, `Age_bin`, `Sex`, `batch` (`cells`/`nuclei`), `10X_version`, and `annotation_level0/1/2`. All 183,161 rows have level-0/1/2 annotations, with no missing category observed. Its metadata would support donor, young/old, sex, assay, and released-label mapping, but its noninteger `X` cannot substitute for raw counts.

Two released band-label discrepancies require reconciliation if this object is ever used for metadata: 470BR is `15-20` here versus `17-20` in the checked SDRF/supplement, and 367C is `70-75` here versus `65-70`. Their `Age_bin` young/old assignments agree. The smaller objects carry comparable donor/age/sex/sample metadata for their respective subsets.

The smaller-object union is **not the full eligible cell-type universe**. At levels 1 and 2, exact-label comparison is missing `Adipocyte`, `Mesothelium`, and `RBC`; level 0 is absent from several targeted files. Targeted files also use different annotation versions/granularity, and the two myonuclei objects overlap. They cannot be concatenated or deduplicated into an equivalent full object from these schema checks.

## Resource use and untouched artifacts

The final seven-object probe set read 13,884,919 range bytes in 222 requests including HEAD requests; development repeats brought the total to about 19 MB in 302 requests. Requests were capped at 64 KiB each, and servers returning anything other than HTTP 206 were rejected. Detailed schema evidence is outside the repository at `/home/reza/projects/agescape-revision-data/muscle/kedlian-alternate-schema-probes.json`.

No existing result or manifest was changed. No PBMC file was modified, and no commit or push was made.
