# Cross-species whole-expression-vector analysis

## Question

The analysis tests whether aging is associated with displacement of complete single-cell expression vectors within each species, tissue, cohort, and cell type. It does not require matched tissues across species. It also does not test DNA damage, chromatin restoration, or any other upstream mechanism.

A whole-vector age association can reflect coordinated mean-expression changes, increased dispersion, cell-state transitions, technical effects, or combinations of these processes. It is not synonymous with increased Shannon entropy or random transcriptional noise.

## Distance calculations

No UMI downsampling is used. Mitochondrial, ribosomal, and sex-linked marker genes are excluded using species-appropriate gene symbols. Each cell must have at least 1,500 UMI in mouse and human data. The existing NMR processing requires at least 2,000 UMI and 700 detected genes.

For normalized profile distance, each cell's raw UMI vector is converted to within-cell proportions. A young reference probability vector is constructed from pooled counts of the same cell type. Hellinger distance is calculated as `sqrt(1 - sum_g sqrt(p_cell,g * p_reference,g))`. It ranges from 0 to 1.

Normalized Hellinger distance without depth adjustment is the primary endpoint. Depth-adjusted Hellinger models are sensitivity analyses. For the literal raw-vector sensitivity analysis, Euclidean distance is calculated between each cell's raw UMI vector and the young mean raw-count vector. This endpoint is expected to depend strongly on library size and is not a pure biological measurement.

Young animals or donors are excluded from their own reference. Every evaluated animal or donor must contribute at least 20 cells. At least 50 young reference cells must remain after this exclusion. Animal or donor means, not individual cells, are the inferential units.

## Species designs

Mouse uses Tabula Muris Senis cells from animals aged 1 or 3 months as the young reference. Each tissue and cell type is tested separately. Continuous age models include sex, with a second model adding mean log UMI depth. At least six animals are required for a unit-level regression. Pooled animal identifiers are retained in the primary atlas-wide analysis because they are part of the experimental design. Atlas-wide and tissue summaries are repeated after excluding pooled identifiers.

Human uses public CELLxGENE Census release 2025-11-08 and the 11 cohorts already analyzed in this project. Only primary, normal-disease cells from donors aged 20 years or older are included. Donors aged 20 to 40 years form the young reference. Each cohort and cell type is tested separately with the same donor and reference thresholds as mouse. Continuous age models include sex, with a second model adding mean log UMI depth. Cohorts or cell types without adequate adult and young-reference coverage are not tested.

NMR uses GSE152886 skin after the project's existing annotation and exclusion of sample 2V. Animals aged 1.1 years form the young reference and animals aged 11.3 years form the aged group. Compartments with all five animals evaluable are tested using all ten assignments of two aged animals among five. The minimum attainable two-sided p value is 0.1. Sex adjustment and continuous-age modeling are not used for this five-animal design.

## Multiplicity and species summaries

Benjamini-Hochberg correction is applied across cell types within each mouse tissue or human cohort, separately for each endpoint and depth model. NMR correction is applied across evaluable compartments separately for each endpoint.

For mouse and human, species and tissue or cohort summaries standardize the outcome within each tissue-by-cohort-by-cell-type unit. Models contain unit fixed intercepts, age, and centered sex. Depth-adjusted versions add centered mean log UMI. Standard errors are clustered by cohort-specific donor identifiers, with t-based inference using the cluster count. Mouse summaries are repeated after excluding pooled animal identifiers.

Species-wide continuous-age effects are reported per 10 percent of maximum lifespan. The committed AnAge values are 4 years for mouse, 122.5 years for human, and 31 years for NMR, giving age scales of 0.4, 12.25, and 3.1 years. The NMR global standardized aged-minus-young contrast is rescaled to the same 3.1-year unit using the observed 10.2-year age gap. These standardized effects support descriptive comparison but do not erase differences in tissue, assay, cell-type composition, or age design.

## Interpretation

The proposed pattern is a positive mouse and human age association accompanied by a smaller NMR association. A non-significant NMR result alone does not demonstrate stability or a between-species difference. Five NMR animals, skin-only sampling, and exact-permutation resolution remain fundamental limitations. Cross-species effect-size comparisons are descriptive unless a direct contrast and adequate replication support stronger inference.

## Outputs

- `results/h3-cross-species-vector-donors.csv`: public animal- or donor-level distances and metadata.
- `results/h3-cross-species-vector-coverage.csv`: evaluable donor-unit counts by tissue and cohort.
- `results/h3-cross-species-vector-tests.csv`: unit-level continuous-age tests and exact NMR contrasts.
- `results/h3-cross-species-vector-global.csv`: tissue, cohort, and species summaries on the common lifespan scale.
