# Follow-up Blocker Report: Adult Human Heart Cell Atlas v2 (`Kanemaru 2023`, `E-MTAB-12916` / `E-MTAB-12919`)

**Status**: BLOCKED / NOT ANALYZED (`notes/followup-dataset-checklist.md` #22: `BLOCKED`; primary donor metadata provides 5-year age brackets only, zero exact continuous ages publicly available)
**Accessions**: ArrayExpress [`E-MTAB-12916`](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12916) (multiome snRNA-seq), [`E-MTAB-12919`](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12919) (multiome snATAC-seq), [`E-MTAB-12975`](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12975) (Visium); CZ CELLxGENE Collection [`3116d060-0a8e-4767-99bb-e866badea1ed`](https://cellxgene.cziscience.com/collections/3116d060-0a8e-4767-99bb-e866badea1ed) (`Heart Global` dataset `1e1e07c3-bfcb-4a0d-91de-c2614a891409`, `704,296` cells, `4.97 GB` NOT downloaded)
**Primary publication**: Kanemaru et al., *Nature* (2023), "Spatially resolved multiomics of human cardiac niches", DOI `10.1038/s41586-023-06311-1`, `PMID: 37438528`, `PMCID: PMC10371870`
**Audit table**: [`followup-organ-16-heart-kanemaru-donors.csv`](../results/followup-organ-16-heart-kanemaru-donors.csv)

## 1. Primary-Source Metadata Verification & `25` vs `22` Donor Reconciliation

1. **Five-year age brackets only across ALL `25` donors**:
   - Inspection of `41586_2023_6311_MOESM3_ESM.xlsx` (`Supplementary Table 1: Anonymized donor information`, exported in [`donors.csv`](../results/followup-organ-16-heart-kanemaru-donors.csv)) and ArrayExpress SDRFs (`E-MTAB-12916.sdrf.txt`, `E-MTAB-12919.sdrf.txt`) confirms that **all `25` donors (`14` Female, `11` Male) have ages reported strictly in 5-year brackets**: `20-25` (`2`), `40-45` (`1`), `45-50` (`6`), `50-55` (`5`), `55-60` (`2`), `60-65` (`4`), `65-70` (`1`), and `70-75` (`4`). Zero exact donor ages exist in any public Kanemaru repository.
   - CZ CELLxGENE (`1e1e07c3-bfcb-4a0d-91de-c2614a891409`) collapses these further into 5 10-year `HsapDv` decade ontologies (`third`, `fifth`, `sixth`, `seventh`, and `eighth decade stage`).
2. **Exact reconciliation of `25` donors (`MOESM3`) vs `22` donors in CZ `Heart Global`**:
   - **`15` re-analyzed donors from Litvinukova et al. (2020, Heart Cell Atlas v1)**: `D1..D8, D11, H2..H7`, all present in CZ `Heart Global` (brackets `40-45` to `70-75`).
   - **`10` newly recruited Kanemaru et al. (2023) donors**: `7` (`AH1, A61, AV3, AV10, AV13, AV14, AH2`) have sc/snRNA-seq / Multiome and are present in CZ `Heart Global` (`15 + 7 = 22` donors, `704,296` cells).
   - **Why `3` donors (`AV1`, `AH5`, `AH6`) are absent from CZ `Heart Global`**: Inspection of `41586_2023_6311_MOESM10_ESM.xlsx` (`Table S8: sequenced_samples_inAnalysis_20`) verifies that `AV1` (`45-50`, `3` Visium-FFPE SAN), `AH5` (`20-25`, `4` Visium SAN), and `AH6` (`50-55`, `1` Visium AVN) were profiled exclusively by Visium spatial transcriptomics with zero scRNA/snRNA libraries.
3. **Additional ArrayExpress / assay issues**:
   - `E-MTAB-12916.sdrf.txt` / `E-MTAB-12919.sdrf.txt` include a multi-donor pooled SAN library (`Mixture of AH1 and A61`) and `E-MTAB-12916` carries an uncorrected string-replacement corruption (`atrioventriculeft cardiac atriumr node`).
   - CZ `Heart Global` spans `8` distinct anatomical compartments (`LV, RV, LA, RA, AX, SP, SAN, AVN`) and mixes fresh dissociated cells (`scRNA` `3' v2/v3`), `snRNA`, and `10x multiome` nuclei.

## 2. Documented Blockers (Why No Modeling Was Run)

1. **No exact continuous donor ages exist publicly**: All `25` primary donor records (`MOESM3`) and ArrayExpress SDRFs provide 5-year anonymized brackets (`20-25` to `70-75`) only (and CZ CELLxGENE gives 10-year `HsapDv` bins), making continuous aging slope estimation impossible without fabricating bin midpoints (e.g., `52.5` as used in secondary re-analyses like Gao et al.).
2. **Assay/region fragmentation on the small new subset**: Only `7` new sc/snRNA donors (`AH1, A61, AV3, AV10, AV13, AV14, AH2`) are independent of Litvinukova 2020, distributed heterogeneously across `SAN`, `AVN`, `LV`, `RV`, `LA`, and `RA` with mixed multi-donor `SAN` pooling (`Mixture of AH1 and A61`).

## 3. Verification Conclusion

Blocked after full primary-source verification (`MOESM3`, `MOESM10`, `E-MTAB-12916.sdrf.txt`, `E-MTAB-12919.sdrf.txt`, CZ CELLxGENE API). The `4.97 GB` `Heart Global` H5AD was NOT downloaded.
