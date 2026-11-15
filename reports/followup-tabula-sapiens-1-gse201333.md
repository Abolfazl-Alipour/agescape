# Follow-up Audit Report: Tabula Sapiens 1.0 (`GSE201333`) Accession Verification and 100% Containment in TS2 (`GSE306755`)

**Status**: SKIP / NOT RE-ANALYZED (`notes/followup-dataset-checklist.md` #2: `SKIP`, superseded by Tabula Sapiens 2.0 `GSE306755`; verified 100% donor and donor-tissue containment)
**Primary Accessions**: GEO [`GSE201333`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201333) (`BioProject PRJNA830893`, `GSM6058681`); The Tabula Sapiens Consortium, *Science* (2022), DOI `10.1126/science.abl4896`, `PMID: 35549404`
**Superseding Release**: Tabula Sapiens 2.0, GEO [`GSE306755`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306755) (`BioProject PRJNA1312052`) and CZ CELLxGENE Collection [`e5f58829-1a66-40b5-a624-9046778e74f5`](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) (`35` datasets, updated in place from the original `2021-07-14` collection)
**Containment Table**: [`followup-tabula-sapiens-1-gse201333-containment.csv`](../results/followup-tabula-sapiens-1-gse201333-containment.csv)

## 1. Primary-Source Verification

1. **GEO `GSE201333` (`GSM6058681`) vs TS2 `GSE306755`**:
   - Primary GEO record `GSE201333` (`GSM6058681`) releases `GSM6058681_TabulaSapiens.h5ad.gz` and `GSM6058681_TabulaSapiens_metadata.csv.gz` (`483,152` cells across `24` tissues and `15` donors: `TSP1, TSP2, TSP3, TSP4, TSP5, TSP6, TSP7, TSP8, TSP9, TSP10, TSP11, TSP12, TSP13, TSP14, TSP15`; `456,101` `10X` cells and `27,051` `smartseq2` cells).
   - `GSM6058681_TabulaSapiens_metadata.csv.gz` contains `12` columns (`cell_id`, `organ_tissue`, `method`, `donor`, `anatomical_information`, `n_counts_UMIs`, `n_genes`, `cell_ontology_class`, `free_annotation`, `manually_annotated`, `compartment`, `gender`) and **zero donor age column**; GEO `GSM6058681` sample characteristics record only `tissue: Various tissues` and `genotype: mixed sample (non-demultiplexed)`.
   - CZ CELLxGENE collection `e5f58829-1a66-40b5-a624-9046778e74f5` (first published `2021-07-14`) has been updated in place to Tabula Sapiens 2.0 (`1,136,218` all-cells object, `35` datasets), and GEO `GSE306755` explicitly documents that *"Nine new donors and four new tissues were collected and analyzed together with the original Tabula Sapiens 1.0 dataset"*.

2. **Exact `100%` Donor and Donor-Tissue Containment Verification**:
   - **Donors (`15/15 = 100%`)**: Every one of the `15` TS1 donors (`TSP1, TSP2, TSP3, TSP4, TSP5, TSP6, TSP7, TSP8, TSP9, TSP10, TSP11, TSP12, TSP13, TSP14, TSP15`) appears in TS2 ([`followup-ts2-donor-summary.csv`](../results/followup-ts2-donor-summary.csv), `24` donors total), alongside the `9` newly added TS2 donors (`TSP17, TSP19, TSP20, TSP21, TSP25, TSP26, TSP27, TSP28, TSP30`).
   - **Donor-tissue pairs (`58/58 = 100%`)**: Every observed `(donor, organ_tissue)` pair in `GSM6058681_TabulaSapiens_metadata.csv.gz` (`58` pairs across `24` tissues) appears in the TS2 crosswalk ([`followup-ts2-donor-organ-crosswalk.csv`](../results/followup-ts2-donor-organ-crosswalk.csv), `123` pairs across `28` tissues). All `58` pairs differ in cell count between TS1 and TS2 due to TS2 re-QC/ambient-RNA processing (`decontx`/`min200-2500`). The `4` tissues unique to TS2 are `Ear`, `Ovary`, `Stomach`, and `Testis`.
   - **Sex & age reconciliation**: `14/14` TS1 donors present in TS2 `media-24.xlsx` match curated sex (`TSP9` is absent from `media-24.xlsx`, leaving `14` age-verified TS1 donors with curated ages `22..74` yr in TS2).

## 2. Disposition & Scientific Rationale

- **Not analyzed separately (`SKIP`, already represented / superseded by TS2)**: Re-analyzing TS1 (`GSE201333`) would strictly subset `15` donors whose cells are already re-integrated into Tabula Sapiens 2.0 (`GSE306755`), violating cohort non-redundancy while providing zero independent donors beyond TS2. No TS1 H5AD was downloaded; `GSM6058681_TabulaSapiens_metadata.csv.gz` was removed after generating the containment table.
