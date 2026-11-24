# Pathway-level uniformity test, matched mouse-human P1m pairs

## Numbers

The requested row-weighted fraction is AGREE divided by all testable function-pair rows; UNCLEAR rows remain in that denominator. Directional agreement, AGREE/(AGREE+DISAGREE), excludes UNCLEAR and is shown alongside it. Null p values are upper-tail, plus-one corrected.

| set | functions | testable rows | AGREE | DISAGREE | UNCLEAR | AGREE/testable | null mean | null 95% range | p | directional agreement | null mean | p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | 987 | 1423 | 0 | 1 | 1422 | 0 | 0 | 0–0 | 1 | 0 | 0 | 1 |
| repeatable | 108 | 502 | 0 | 0 | 502 | 0 | 0 | 0–0 | 1 | NA | 0 | NA |

Noise floors are gene-count matched: for a row with n genes on a side, its threshold is the 95th percentile of absolute permuted means among all side-testable rows with that n (nearest available n if absent). Across observed testable rows, mouse thresholds range 1.87658–2.73525 log2FC and human thresholds range 0.000282529–0.00478299 slope units. At n=5 the thresholds are 2.73525 and 0.00478299, respectively. Every row's exact thresholds are in the main CSV; estimates use 2,000 shuffles.

### Repeatable functions with at least three AGREE pairs

The repeatable analysis contains functions testable in at least three pairs. The table below further requires at least three AGREE pairs; directions are mouse/human.

None.

### Top disagreeing functions

| GO id | name | tested pairs | AGREE | DISAGREE | UNCLEAR | directions |
| --- | --- | --- | --- | --- | --- | --- |
| GO:0045321 | leukocyte activation | 1 | 0 | 1 | 0 | up/down:1 |

### Per-pair summary

| mouse unit | human cohort | human unit | pair genes | testable functions | AGREE | DISAGREE | UNCLEAR | AGREE/testable | directional agreement |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Kidney::B cell | Kidney_LivingDonor | B cell | 4 | 0 | 0 | 0 | 0 | NA | NA |
| Kidney::B cell | Kidney_SCR_v15 | B cell | 4 | 0 | 0 | 0 | 0 | NA | NA |
| Kidney::T cell | Kidney_LivingDonor | T cell | 0 | 0 | 0 | 0 | 0 | NA | NA |
| Kidney::T cell | Kidney_SCR_v15 | T cell | 0 | 0 | 0 | 0 | 0 | NA | NA |
| Kidney::epithelial cell of proximal tubule | Kidney_LivingDonor | epithelial cell of proximal tubule | 15 | 94 | 0 | 0 | 94 | 0 | NA |
| Kidney::epithelial cell of proximal tubule | Kidney_SCR_v15 | epithelial cell of proximal tubule | 14 | 93 | 0 | 0 | 93 | 0 | NA |
| Kidney::kidney collecting duct principal cell | Kidney_SCR_v15 | kidney collecting duct principal cell | 17 | 94 | 0 | 0 | 94 | 0 | NA |
| Kidney::kidney distal convoluted tubule epithelial cell | Kidney_LivingDonor | kidney distal convoluted tubule epithelial cell | 7 | 26 | 0 | 0 | 26 | 0 | NA |
| Kidney::kidney distal convoluted tubule epithelial cell | Kidney_SCR_v15 | kidney distal convoluted tubule epithelial cell | 6 | 18 | 0 | 0 | 18 | 0 | NA |
| Kidney::kidney loop of Henle thick ascending limb epithelial cell | Kidney_SCR_v15 | kidney loop of Henle thick ascending limb epithelial cell | 11 | 105 | 0 | 0 | 105 | 0 | NA |
| Lung::B cell | ImmuneAging_BAL | B cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::B cell | Lung_CellRef | B cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::B cell | Lung_Integrated | B cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::CD4-positive, alpha-beta T cell | Lung_Integrated | CD4-positive, alpha-beta T cell | 0 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::CD8-positive, alpha-beta T cell | Lung_Integrated | CD8-positive, alpha-beta T cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::T cell | Lung_Integrated | T cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::alveolar macrophage | Lung_CellRef | alveolar macrophage | 0 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::alveolar macrophage | Lung_Integrated | alveolar macrophage | 0 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::classical monocyte | ImmuneAging_BAL | classical monocyte | 4 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::classical monocyte | Lung_CellRef | classical monocyte | 4 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::classical monocyte | Lung_Integrated | classical monocyte | 4 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::non-classical monocyte | ImmuneAging_BAL | non-classical monocyte | 5 | 2 | 0 | 0 | 2 | 0 | NA |
| Lung::non-classical monocyte | Lung_CellRef | non-classical monocyte | 5 | 2 | 0 | 0 | 2 | 0 | NA |
| Lung::non-classical monocyte | Lung_Integrated | non-classical monocyte | 5 | 2 | 0 | 0 | 2 | 0 | NA |
| Lung::vein endothelial cell | Lung_CellRef | vein endothelial cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Lung::vein endothelial cell | Lung_Integrated | vein endothelial cell | 1 | 0 | 0 | 0 | 0 | NA | NA |
| Marrow::hematopoietic precursor cell | Yazar_Blood | hematopoietic precursor cell | 187 | 987 | 0 | 1 | 986 | 0 | 0 |
| Marrow::naive B cell | Yazar_Blood | naive B cell | 0 | 0 | 0 | 0 | 0 | NA | NA |

## Interpretation

Observed AGREE/testable is 0.000000, below or equal to the permutation mean of 0.000000; the upper-tail p value is 1. The corresponding repeatable-set values are 0.000000, 0.000000, and 1. These are descriptive function-level direction comparisons, not evidence that GO annotations themselves are causal aging mechanisms.

## Design and grouping

The committed P1m derivation supplied 1,914 genes and 36 eligible mouse units; only the 28 matched pairs in `results/followup-mature-matched-p1m.csv` were analyzed, without floored/unfloored arms. Mouse effects are per-unit 30-month versus 18-month log2FC values from the existing DE CSVs. The P1m list was not re-derived. The cached 299 pair-gene rows cover every one of the 22 pairs with tested human genes; the other six pairs have zero tested genes. All 299 rows were checked against freshly computed slopes from the owner-approved human cache.

GO grouping used direct database mapping through `org.Mm.eg.db` `GOALL`, filtered to `ONTOLOGYALL == BP`, with names from `GO.db`; this includes direct and inherited GO Biological Process annotations. clusterProfiler was not used. The mapping contains 1,851 annotated mouse genes, 9,550 GO terms, and 193,603 gene-term memberships. A gene can belong to many terms. The CSV is the complete 9,550 terms × 28 pairs grid; absent directions are NA, never zero-filled.

Human orthologs and collision handling are exactly those of the matched test: use `human_ortholog` from the P1m gene CSV, sort within each mouse unit by human ortholog, unit FDR, and mouse symbol, and keep the first mouse gene per human ortholog. Human slopes are continuous-age OLS fits for donors aged 30+ on donor-cell-type log1p counts-per-10,000 pseudobulks, with sex included where variable and identifiable. Mouse n counts all annotated P1m DE genes in that unit; human n counts annotated, ortholog-mapped genes present with finite slopes in that human unit.

A row is testable only with at least 5 genes on each side. Otherwise it is UNTESTABLE. For testable rows, either absolute mean below its noise floor makes the row UNCLEAR; otherwise equal mean signs give AGREE and opposite signs DISAGREE. Up/down and positive/negative counts and shares are recorded in the main CSV.

The null uses 2,000 deterministic shuffles (seed 20261126) of gene labels relative to the complete GO membership table. Each shuffle preserves every gene's number of GO memberships, every GO term's global size, and the linkage between each gene's mouse and human measurements; the same shuffle is applied across all 28 pairs. Counts, means, noise-floor classifications, and both agreement fractions are recomputed. The repeatable null reselects functions testable in at least three pairs within each shuffle. Noise floors are fixed at the gene-count-matched pooled null 95th percentiles for observed and null classifications.

## Verification

The match check found 28 unique matched pairs. The independent spot-check selected `GO:0045321` (leukocyte activation) in `Marrow::hematopoietic precursor cell` versus `Yazar_Blood::hematopoietic precursor cell`. Recomputing from the unit DE CSV, P1m gene list, GO mapping, and pair-gene slopes gave mouse n=41, mean=2.26977161879; human n=37, mean=-0.000519883457624; verdict=DISAGREE. Counts, means, and direction counts exactly match the serialized main CSV.

After writing, the script re-read both result CSVs and recomputed the pooled and per-pair tables used in this report from those files. Assertions cover 28-pair completeness, unique function-pair keys, NA means at zero counts, cached pair counts, and exact observed/null summary consistency.

## Limits

- GO terms overlap extensively and are not independent tests. Row weighting gives more influence to terms and pairs with many annotations, and the pooled fractions are not gene- or pair-independent estimates.
- The 95th-percentile noise floors adapt to gene count but pool across functions and pairs; they do not adapt to function identity, tissue, effect distribution, or donor number. UNCLEAR is a threshold classification, not a confidence interval or a multiple-testing-adjusted pathway p value.
- Permutation tests exchangeability of gene-to-function labels while preserving the measured values and annotation multiplicity. It does not model donor sampling, cell-type mismatches, mouse gene-selection uncertainty, GO annotation error, or dependence among terms.
- Only 22 pairs have any tested human genes, and gene coverage is highly unequal. Six zero-gene pairs cannot contribute testable functions. The repeatable set reduces sparsity but does not remove unequal tissue, cohort, donor, or gene-count representation.
- Mouse effects are selected 18m-to-30m two-group log2FC values; human effects are cross-sectional continuous-age slopes on a different scale. Equal signs are not equivalent effect sizes or proof of conserved causal programs. The mouse 30m group is all male, so mouse age and sex remain confounded.
- Same-named cell types are not proven homologous across species. Existing ortholog and exact-name matching choices propagate into every function result.

## Reproducibility and manuscript boundary

Ran with Python 3.11.17, pandas 2.3.3, NumPy 2.4.6, SciPy 1.17.1; R version 4.3.3 (2024-02-29); org.Mm.eg.db 3.18.0; GO.db 3.18.0; AnnotationDbi 1.64.1. Database source metadata follows.

| metadata key | value |
| --- | --- |
| R_version | R version 4.3.3 (2024-02-29) |
| org.Mm.eg.db_package_version | 3.18.0 |
| GO.db_package_version | 3.18.0 |
| AnnotationDbi_package_version | 1.64.1 |
| BiocManager_package_version | 1.30.27 |
| org.Mm.eg.db:DBSCHEMAVERSION | 2.1 |
| org.Mm.eg.db:Db type | OrgDb |
| org.Mm.eg.db:Supporting package | AnnotationDbi |
| org.Mm.eg.db:DBSCHEMA | MOUSE_DB |
| org.Mm.eg.db:ORGANISM | Mus musculus |
| org.Mm.eg.db:SPECIES | Mouse |
| org.Mm.eg.db:EGSOURCEDATE | 2023-Sep11 |
| org.Mm.eg.db:EGSOURCENAME | Entrez Gene |
| org.Mm.eg.db:EGSOURCEURL | ftp://ftp.ncbi.nlm.nih.gov/gene/DATA |
| org.Mm.eg.db:CENTRALID | EG |
| org.Mm.eg.db:TAXID | 10090 |
| org.Mm.eg.db:GOSOURCENAME | Gene Ontology |
| org.Mm.eg.db:GOSOURCEURL | http://current.geneontology.org/ontology/go-basic.obo |
| org.Mm.eg.db:GOSOURCEDATE | 2023-07-27 |
| org.Mm.eg.db:GOEGSOURCEDATE | 2023-Sep11 |
| org.Mm.eg.db:GOEGSOURCENAME | Entrez Gene |
| org.Mm.eg.db:GOEGSOURCEURL | ftp://ftp.ncbi.nlm.nih.gov/gene/DATA |
| org.Mm.eg.db:KEGGSOURCENAME | KEGG GENOME |
| org.Mm.eg.db:KEGGSOURCEURL | ftp://ftp.genome.jp/pub/kegg/genomes |
| org.Mm.eg.db:KEGGSOURCEDATE | 2011-Mar15 |
| org.Mm.eg.db:GPSOURCENAME | UCSC Genome Bioinformatics (Mus musculus) |
| org.Mm.eg.db:GPSOURCEURL | NA |
| org.Mm.eg.db:GPSOURCEDATE | 2021-Nov23 |
| org.Mm.eg.db:ENSOURCEDATE | 2023-May10 |
| org.Mm.eg.db:ENSOURCENAME | Ensembl |
| org.Mm.eg.db:ENSOURCEURL | ftp://ftp.ensembl.org/pub/current_fasta |
| org.Mm.eg.db:UPSOURCENAME | Uniprot |
| org.Mm.eg.db:UPSOURCEURL | http://www.UniProt.org/ |
| org.Mm.eg.db:UPSOURCEDATE | Mon Sep 18 16:19:56 2023 |
| GO.db:GOSOURCENAME | Gene Ontology |
| GO.db:GOSOURCEURL | http://current.geneontology.org/ontology/go-basic.obo |
| GO.db:GOSOURCEDATE | 2023-07-27 |
| GO.db:Db type | GODb |
| GO.db:package | AnnotationDbi |
| GO.db:DBSCHEMA | GO_DB |
| GO.db:GOEGSOURCEDATE | 2023-Sep11 |
| GO.db:GOEGSOURCENAME | Entrez Gene |
| GO.db:GOEGSOURCEURL | ftp://ftp.ncbi.nlm.nih.gov/gene/DATA |
| GO.db:DBSCHEMAVERSION | 2.1 |

Script: `scripts/run_mature_matched_p1m_pathways.py`. Use `--refresh-human-cache` to rebuild the external cache and pair-gene rows. Started `2026-11-24T17:33:32.752971+00:00`; runtime 76.518 seconds. Main data: `results/followup-mature-matched-p1m-pathways.csv`; pooled data: `results/followup-mature-matched-p1m-pathways-summary.csv`. Human cache and GO mapping are under `/home/reza/projects/agescape-revision-data/matched-p1m`.

No manuscript or submission file was changed. Any TeX change is only a proposal for later owner review: if these descriptive results are cited, report both the observed and null fractions, the noise-floor rule, and the limits above rather than claiming pathway conservation from agreement alone.
