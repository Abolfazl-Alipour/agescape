# Follow-up Access & Confounding Audit: SenNet Human Liver Multiome (`PMC12903365`)

**Status**: ACCESS & CONFOUNDING AUDIT ONLY (`notes/followup-dataset-checklist.md` #25: `BLOCKED`; no DAC application submitted, zero SenNet/HTAN files downloaded)
**Primary publication**: Karpova et al., *Cell Genomics* (2026), "Cellular senescence in human liver under normal aging and cancer", DOI `10.1016/j.xgen.2025.101133`, `PMID: 41576948`, `PMCID: PMC12903365`
**Accessions / Portals**:
- **Human normal liver**: [SenNet Data Portal](https://sennetconsortium.org/) (`TMC - Washington University`, organ `UBERON:0002107`; code on [GitHub `ding-lab/SenNet_liver`](https://github.com/ding-lab/SenNet_liver) and Zenodo DOI `10.5281/zenodo.17584291`)
- **Human mCRC cohort**: [HTAN DCC Data Portal](https://data.humantumoratlas.org/) (`HTAN WUSTL Atlas`)
- **Mouse only (GEO)**: `GSE293958` (2-week snRNA-seq), `GSE310392` (Xenium), `GSE311064` (snRNA-seq); ZERO human GEO RNA-seq accessions exist in `PMC12903365`
**Audit table**: [`followup-organ-15-liver-sennet-pmc12903365-audit.csv`](../results/followup-organ-15-liver-sennet-pmc12903365-audit.csv)

## 1. Controlled-Access Terms (Verified Without Application)

1. **Governance (`https://sennetconsortium.org/external-data-use/`, approved June 15, 2022)**:
   - SenNet follows NIH Genomic Data Sharing (`GDS`) Policy (`US 45 CFR 46`).
   - **Controlled-access sequencing**: DNA and RNA raw sequencing data (`data_access_level = protected`) require permission reviewed and granted by a designated **NIH Data Access Committee (`DAC`)** under **NIH Security Best Practices** (General Research Use `GRU` or unrestricted consent groups, with explicit prohibition against participant re-identification).
   - Processed / imaging data tiers are designated `public` (`CC-BY 4.0` after consortium verification/release).
2. **Inventory on SenNet Public Search API (`TMC - Washington University`, `UBERON:0002107`, [`audit.csv`](../results/followup-organ-15-liver-sennet-pmc12903365-audit.csv))**:
   - `164` published Liver dataset records across **`43` unique human donors** (ages `22` to `88` yr; `<40`: `15`, `40-59`: `8`, `>=60`: `20`; `26` Female, `16` Male, `1` Unknown in SenNet source metadata):
     - `39` `RNAseq` datasets across `34` unique donors (`<40`: `15`, `40-59`: `0`, `>=60`: `19`; `22` to `88` yr): **`protected`**
     - `34` `ATACseq` datasets across `31` unique donors (`<40`: `14`, `40-59`: `0`, `>=60`: `17`; `22` to `88` yr): **`protected`**
     - `34` `10X Multiome` datasets across `31` unique donors (`<40`: `14`, `40-59`: `0`, `>=60`: `17`; `22` to `88` yr): **`protected`**
     - `27` `Xenium` spatial imaging datasets (`25` donors): `public`
     - `10` `PhenoCycler` protein imaging datasets (`10` donors): `public`
     - `13` `10x Multiome [Salmon + ArchR + Muon]` (`10` donors, `~3.5-7.4 GB` `.h5mu` per dataset) + `7` `RNAseq [Salmon]` (`5` donors) derived datasets spanning **`14` unique donors** (`[22, 30, 33, 34, 37, 37, 37]` vs `[61, 63, 69, 71, 72, 74, 75]`, zero middle-aged `40-59`): **`public`**
   - Note that the `PMC12903365` `Table S1` (`1-s2.0-S2666979X26000017` Elsevier supplemental spreadsheets) was not retrieved automatically because Elsevier's supplemental CDN returned a Google reCAPTCHA challenge; no bypass was attempted, and clinical Metavir fibrosis/steatosis/ALT scores per SenNet ID are absent from the public SenNet API donor metadata.

## 2. Cohort Design, Fibrosis Confounding, and Clinical Selection Limits

1. **High baseline liver pathology (`60%` fibrosis, `3` cirrhosis)**:
   - `PMC12903365` Results (`p[18]`, STAR Methods `p[146]`) report that the `95` consented donors (`89` with parenchyma, `43` profiled multi-omically) were surgical patients undergoing **cholecystectomy, bariatric surgery, or open abdominal surgery** at Barnes-Jewish Hospital (`IRB 202108156`).
   - Across profiled "normal" liver donors: `79%` White, `60%` Female, median BMI `29 kg/m2`, **`60%` had at least mild histological fibrosis (`Metavir >= 1`)**, and **`3` donors had cirrhosis** (Figure S1A). Pathologist fibrosis/steatosis/inflammation scores are stored in `Table S1`, not in the open SenNet donor JSON.
   - `PMC12903365` Results (`p[48]`) further document acute liver-injury outliers (`ALT = 242, AST = 71` in a young donor with `>40%` `CRP+` hepatocytes) and report directly that `CDKN1A+` hepatocytes associate with fibrosis (`Metavir >= 1`, univariate `p = 0.007`), demonstrating that liver pathology directly confounds the aging signal unless strictly stratified/filtered.
2. **Bimodal selection on Multiome (`n = 34`) and Public Subset (`n = 14`)**:
   - Multiome snRNA-seq + snATAC-seq (`n = 34`) intentionally prioritized extreme ages (`15` `<40` yr and `19` `>=60` yr, **`0` donors aged `40-59` yr**), and the `14` donors with public derived `Salmon`/`Muon` objects are likewise bimodal (`7` aged `22-37` vs `7` aged `61-75`, without public Metavir fibrosis labels to exclude the `60%` fibrotic / `3` cirrhotic donors).
   - Meanwhile, the `HTAN` `mCRC` cohort (`26` patients: `6` untreated, `20` chemotherapy-treated; chemotherapy expanded senescent hepatocytes `5`-fold relative to aging) is cancer/chemotherapy-confounded.

## 3. Verification Conclusion

No NIH DAC application or HTAN controlled-access request was submitted, and no SenNet `h5mu` / raw sequencing files were downloaded. Primary sequencing (`39` RNAseq, `34` ATACseq, `34` 10X Multiome) is controlled-access (`protected`, NIH DAC / GDS Policy), and even the partial `14`-donor public derived tier lacks open Metavir fibrosis annotations in a cohort where `60%` have `Metavir >= 1` fibrosis (including `3` cirrhotic surgical patients) and `0` donors span ages `40-59`.
