# Item 2 Report: Verification of the LACA Multi-Tissue Naked Mole Rat Aging Atlas Listing

Date: 2026-11-14
Status: **UNVERIFIED / UNRELEASED PORTAL PLACEHOLDER (NOT A REAL DOWNLOADABLE PUBLIC DATASET)**
Action taken on expression data: None (no download attempted, per safety and verification rules).

## 1. Summary Verdict

The "NMR Multi-tissue Lifespan" entry (`LACA-2024-002`, claiming `1.2M` cells and `380GB`) and the companion package listing "NMR Full Atlas v2.1" (claiming `660K` cells, `18 tissues`, `114 annotated cell types`, and `8.2 GB h5ad + metadata`) on the Longevity & Aging Cell Atlas portal (`https://www.longevityaging.net/LACA/`) **cannot be verified as a real, accessible, or published primary dataset**. Primary-source inspection shows:

1. No `Heterocephalus glaber` (`SP239`) or `LACA-2024-002` file exists on the LACA download server.
2. Direct `SP239` portal and data endpoints return `HTTP 403 Forbidden`.
3. The portal's interactive viewer (`atlas.php?dataset=LACA-2024-002`) explicitly contains a developer comment labeling `SP239` an `"unpublished self-test"` and silently redirects `LACA-2024-002` to the public human atlas (`data/SP121/umap`).
4. Exhaustive searches across NCBI GEO, NCBI BioProject, NCBI PubMed, ENA, and CELLxGENE identify zero 1.2M-cell or 660K-cell multi-tissue naked mole rat single-cell aging datasets or accessions.

## 2. Direct Portal Inspection (`https://www.longevityaging.net/LACA/`)

All checks were conducted read-only on `2026-11-14`:

| Endpoint Checked | HTTP Status | Observed Finding |
| --- | --- | --- |
| `https://www.longevityaging.net/LACA/datasets.php` | `200 OK` | Displays card `LACA-2024-002`, title `"NMR Multi-tissue Lifespan"`, description `"First complete multi-tissue atlas of naked mole rat aging from young to geriatric ages."`, metadata `"Naked Mole Rat | Multi-tissue | Cells: 1.2M | Size: 380GB | snRNA-seq | Multi-tissue | Primary"`, and button `download.php?dataset=LACA-2024-002`. Lists no GEO/BioProject/SRA/GSA accession or DOI. |
| `https://www.longevityaging.net/LACA/packages.php` | `200 OK` | Displays an internally contradictory card `"NMR Full Atlas v2.1"` claiming `"Complete naked mole-rat multi-tissue atlas with 660K cells, 18 tissues, and 114 annotated cell types"`, `"8.2 GB"`, `"h5ad + metadata"`, with a non-functional dead link `<a href="#" class="pkg-download" aria-label="Download NMR Full Atlas (8.2 GB)">`. Lists 0 DOIs and 0 accessions. |
| `https://www.longevityaging.net/LACA/download.php` and `download.php?dataset=LACA-2024-002` | `200 OK` | States `"Every file listed here exists on this server"`, yet contains **zero** matches for `SP239`, `LACA-2024-002`, `Naked`, `Heterocephalus`, `glaber`, or `NMR`. The 1.2M-cell `380GB` dataset and 8.2 GB package do not appear in the server package table. |
| `https://www.longevityaging.net/LACA/species.php?SP_ID=SP239` | `200 OK` | Lists `Heterocephalus glaber` (`SP239`) linking to `content.php?SP_ID=SP239` and 17 tissue subpages (`content.php?tisID=...&SP_ID=SP239`). |
| `https://www.longevityaging.net/LACA/content.php?SP_ID=SP239` and `content.php?tisID=sys09tis02&SP_ID=SP239` | `403 Forbidden` | Both species and tissue pages return `HTTP Error 403: Forbidden`. |
| `https://www.longevityaging.net/LACA/data/SP239/umap` | `403 Forbidden` | Underlying UMAP data path is blocked (`HTTP Error 403: Forbidden`). |
| `https://www.longevityaging.net/LACA/atlas.php?dataset=LACA-2024-002` | `200 OK` | Hardcoded JavaScript explicitly reveals the dataset is unpublished and redirected to Human (`SP121`): `// laca_atlas_default_hidden 2026-10-17: SP239 unpublished self-test; path switched to public human atlas` and `'LACA-2024-002': 'data/SP121/umap'`. |
| `https://www.longevityaging.net/LACA/index.php` and `publications.php` | `200 OK` | Zero DOIs and zero public accessions (`GSE*`, `PRJNA*`, `PRJEB*`, `CRA*`, `HRA*`, `OMIX*`) found on either page. |

## 3. Primary Public Repository Verification

To verify whether a multi-tissue single-cell aging atlas of *Heterocephalus glaber* exists independently under any primary accession:

1. **NCBI GEO DataSets (`gds`)**: Querying `("Heterocephalus glaber"[Organism] OR "naked mole rat" OR "naked mole-rat") AND ("single cell" OR "single-cell" OR "snRNA" OR "scRNA")` returned 48 records. All single-cell datasets are tissue-focused, not multi-tissue aging atlases:
   - `GSE152886` (Savina et al. 2022, 6 samples; NMR skin aging, already included in Agescape).
   - `GSE185720`, `GSE185721`, `GSE185722`, `GSE185723`, `GSE185724`, `GSE202903` (hematopoiesis / bone marrow / spleen / thymus / blood / lymph node).
   - `GSE132642` (Hilton et al. 2019, spleen and circulating immune cells).
   - `GSE214390` (NMR T cells across tissues/thymus, 66 samples) and `GSE305565` (autochthonous lung cancer model, 3 samples).
2. **NCBI BioProject and ENA (`tax_id=10181`)**: Querying all *Heterocephalus glaber* projects returned 43 BioProjects (`PRJNA1510787` reference genome, `PRJNA1306540` lung cancer, `PRJNA1058087` cross-species spleen study on `Myotis pilosus` reusing Hilton et al. 2019 NMR spleen data, `PRJNA826455` / `GSE152886` skin aging, `PRJEB47073`/`58293`/`58299`/`58852` skin and liver DMBA/TPA, `PRJNA386219` / `GSE98746` bulk liver RNA-seq, `PRJDB7685` bulk tissues). Zero projects match a multi-tissue 660K-cell or 1.2M-cell snRNA-seq aging atlas.
3. **NCBI PubMed**: Querying `Heterocephalus glaber` / `naked mole-rat` with `single-cell` / `snRNA` / `atlas` returned 17 articles (including Savina 2022 skin PMID `35507806`, Emmrich 2022 hematopoiesis PMID `35694726`, Hilton 2019 immunity PMID `31751331`, ovary PMID `39126475`, T-cell evolution PMID `38605005`, spleen PMID `37681451`). None reports a 1.2M-cell or 660K-cell multi-tissue NMR aging atlas.
4. **CELLxGENE Curation API**: Across all `406` public collections (`https://api.cellxgene.cziscience.com/curation/v1/collections?visibility=PUBLIC`), **`0`** collections reference *Heterocephalus glaber* or naked mole-rat.

## 3.5 Direct Record Verification of `PRJNA1058087` and Scope Note on `GSE98746` (No Download)

1. **`PRJNA1058087` is NOT an aging naked mole-rat single-cell spleen cohort**:
   - **NCBI BioProject (`uid=1058087`, registered `2023-12-27`, Beijing Institute of Genomics CAS)**: Registered under `Myotis pilosus` (`taxID=203696`) with `0` NCBI SRA runs and external link to CNCB-NGDC BioProject `PRJCA014972` / GSA `CRA009824` (registered species: `Mus musculus` and `Myotis ricketti`).
   - **Accompanying paper (*Aging Cell* 2023, `PMID: 37681451`, `PMCID: PMC10652311`, `10.1111/acel.13982`)**: Methods Sections `2.1`, `2.3`, and `3.1` report generating new 10x single-cell splenic RNA-seq solely from `Myotis` bats (`2` adult males, `2` adult females; alongside bat genome assembly and mouse/cell-line assays), while naked mole-rat splenic immune cells were reused from Hilton et al. 2019 (`GSE132642`, adult NMR vs mouse natural/LPS/saline immune repertoire comparison, no young-vs-old NMR aging arm). `PRJNA1058087` appears in text searches for `"naked mole rat"` solely because the abstract compares bat, naked mole-rat, human, and mouse splenic cells; it adds no new NMR aging single-cell data.
2. **Bulk liver `GSE98746` (`PRJNA386219`)**:
   - Comprises `18` samples total (`Rattus norvegicus` and `Heterocephalus glaber`, `"Aging in naked mole-rats and rats"`), including `9` naked mole-rat bulk liver RNA-seq samples (`3` old). As a small bulk tissue comparison, a per-sample 139-gene score cannot separate within-cell expression changes from cell-type composition shifts and could at most serve as a suggestive side check, not a standalone paper claim.

## 4. Conclusion and Recommendation

`LACA-2024-002` (`1.2M` cells) and `"NMR Full Atlas v2.1"` (`660K` cells) on `longevityaging.net/LACA/` are unreleased/internal placeholders (`SP239 unpublished self-test`) with no downloadable artifact, no primary repository accession, and a silent frontend redirect to human data (`SP121`). Therefore, Agescape's NMR comparative evidence remains strictly restricted to the verified public epidermal cohort `GSE152886` (Savina et al. 2022, `n = 5` animals after excluding `2V`), and no download or analysis of LACA NMR should be attempted.
