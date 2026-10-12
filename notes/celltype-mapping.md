# Cross-Species Cell Ontology Mapping: Mouse ↔ Human ↔ Naked Mole-Rat

## Pre-registered Methodological Decision (§4.5.2)
- Authority: Cell Ontology (CL) terms.
- Mapping standard: Curated one-to-one or one-to-many ontology mappings with written biological rationale.
- Constraint: Conservation analysis is computed on **tissue-matched pairs only**.
- Exclusion rule: Unmatched cell types are formally classified as `excluded, reported` (never scored as zero, never dropped silently).

| Mouse Tissue | Mouse Cell Type (`cell_ontology_class`) | Mouse CL ID | Human Mapped Cell Type | Human CL ID | Mapping Disposition | Rationale / Marker Verification |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Lung** | fibroblast of lung | CL:0002553 | lung fibroblast | CL:0002553 | Matched | 1:1 identity; verified Col1a1/COL1A1, Pdgfra/PDGFRA |
| **Lung** | adventitial cell | CL:0000000 | lung adventitial fibroblast | CL:0002553 | Matched | Subtype collapsed to perivascular adventitial fibroblast (Pi16+) |
| **Lung** | alveolar macrophage | CL:0000583 | alveolar macrophage | CL:0000583 | Matched | 1:1 identity; verified Siglecf/SIGLEC8, Pparg/PPARG |
| **Lung** | classical monocyte | CL:0000860 | classical monocyte | CL:0000860 | Matched | 1:1 identity; verified CD14+/Ly6c2+ |
| **Lung** | non-classical monocyte | CL:0000875 | non-classical monocyte | CL:0000875 | Matched | 1:1 identity; verified CD16+/Itgal+ |
| **Lung** | NK cell | CL:0000623 | natural killer cell | CL:0000623 | Matched | 1:1 identity; verified Ncr1/NCR1, Klrb1c |
| **Lung** | bronchial smooth muscle cell | CL:0002598 | bronchial smooth muscle cell | CL:0002598 | Matched | 1:1 identity; verified Acta2/ACTA2, Tagln/TAGLN |
| **Lung** | B cell | CL:0000236 | B cell | CL:0000236 | Matched | 1:1 identity; verified Cd19/CD19, Ms4a1/MS4A1 |
| **Lung** | CD4-positive, alpha-beta T cell | CL:0000624 | CD4-positive, alpha-beta T cell | CL:0000624 | Matched | 1:1 identity; verified Cd4/CD4, Cd3e/CD3E |
| **Kidney** | kidney proximal convoluted tubule epithelial cell | CL:0002306 | proximal tubule epithelial cell | CL:0002306 | Matched | 1:1 identity; verified Slc34a1/SLC34A1 |
| **Kidney** | kidney loop of Henle thick ascending limb epithelial cell | CL:0002341 | thick ascending limb cell | CL:0002341 | Matched | 1:1 identity; verified Slc12a1/SLC12A1 |
| **Kidney** | epithelial cell of proximal tubule | CL:0002306 | proximal tubule epithelial cell | CL:0002306 | Matched | S3 segment mapped to human proximal tubule |
| **Kidney** | fenestrated cell | CL:0000650 | kidney fenestrated endothelial cell | CL:0000650 | Matched | Glomerular/peritubular endothelium (Pecam1/PECAM1) |
| **Kidney** | macrophage | CL:0000235 | kidney macrophage | CL:0000235 | Matched | 1:1 identity; verified Cd68/CD68 |
| **Kidney** | kidney distal convoluted tubule epithelial cell | CL:0002340 | distal convoluted tubule cell | CL:0002340 | Matched | 1:1 identity; verified Slc12a3/SLC12A3 |
| **Kidney** | kidney collecting duct principal cell | CL:0005011 | collecting duct principal cell | CL:0005011 | Matched | 1:1 identity; verified Aqp2/AQP2 |
| **Kidney** | kidney cortex artery cell | CL:0002543 | renal arterial vascular smooth muscle cell | CL:0002543 | Matched | Artery wall smooth muscle / endothelial combo |
| **Skin (NMR)** | skin fibroblast | CL:0002620 | skin fibroblast | CL:0002620 | Matched | Cross-species skin aging comparison (Savina 2022 GEO GSE152886) |

### Formal Exclusion Log (Unmatched Types)
- **Tongue basal keratinocytes:** No tissue-matched human tongue cohort with ≥6 donors and age SD ≥12y in primary SenNet/TS v1 release. Disposition: `excluded, reported`.
- **Trachea chondrocytes:** Sampled in mouse 10x droplet but absent in human adult SenNet healthy aging core panel. Disposition: `excluded, reported`.
- **Mammary gland stroma:** Mouse-specific parity/endocrine variance without matched human donor cohort. Disposition: `excluded, reported`.
