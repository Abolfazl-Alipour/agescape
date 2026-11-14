# Step 7 ovary analyses: Jin and postmenopausal cohorts

Both frozen designs were executed sequentially at repository HEAD `0610c34`. Jin GSE202601 was completed before GSE329417 acquisition. The inferential unit is the donor in every model; nuclei and tissue pieces are not independent replicates. The cohorts were never pooled. All associations are cross-sectional, not causal, and no premenopausal claim is made.

GSE329417 remains a **bioRxiv preprint** (DOI `10.64898/2026.05.11.724381`); no peer-reviewed journal publication was verified. Jin is a peer-reviewed Nature Aging article (DOI `10.1038/s43587-024-00762-5`).

## Frozen designs and multiplicity

- `data/manifest/ovary_jin_analysis_plan.json`: four young versus four reproductively-aged donors, each once; reproductively-aged-minus-young; primary `group + mean_log_umi`, with group-only sensitivity. No continuous-age coefficient was fitted. The ovary5 snRNA `ovary47` versus GEO/snATAC 49 conflict remains unresolved and numeric age was never used.
- `data/manifest/ovary_postmenopausal_analysis_plan.json`: exact age per 10 years within postmenopausal donors. Primary adjustment is depth, three donor-specific plex fractions (plex 4 reference), BMI, and non-ovarian abnormal state (pre-malignant/malignant versus benign). Technical-only and unadjusted sensitivities plus secondary parity-augmented models were frozen before outcomes.
- P1 uses established human `P1_Clean`, per-cell log1p(CPM) donor means, detected up/down genes, and a 200-draw expression-decile-matched null (seed 20261105). Pseudobulk sums raw/corrected counts before log1p(CPM). CLR uses every released coarse label and a 0.5-cell pseudocount.
- BH families are separate by cohort, endpoint, and model. Gene-level BH spans all tested cell-type-by-gene combinations within a model. Conventional and matched-null P1 p-values have separate primary families. Sensitivity families are reported separately, not used to select a favorable model.

## Resource and acquisition audit

All expression downloads used `scripts/resource_guard.py`, sequential execution, a 30-GB disk reserve, and per-file caps. Analysis/decompression jobs used an 8-GB memory cap with swap disabled in the job scope. Jin was processed one donor at a time. The atlas used direct HDF5 categorical metadata and CSR indptr-bounded 256-nucleus chunks after AnnData initialization exceeded the cap; the cap was not raised. Neither RAW tar nor FASTQ data was downloaded.

| Cohort | Download files | Compressed/download bytes | Decompressed bytes | Cache peak RSS GiB | Analysis peak RSS GiB |
|---|---:|---:|---:|---:|---:|
| GSE202601 | 8 | 118,288,235 |  | 0.901 | 1.023 |
| GSE329417 | 1 | 5,746,830,392 | 21,429,914,240 | 0.806 | 0.636 |

## Jin GSE202601

### Count, barcode, and eligibility gates

Eight per-donor Cell Ranger H5 files were validated: **69,548,923 stored values** are finite, nonnegative integers. The files contain 43,494 barcodes; all **42,568 released annotations** matched exactly one matrix barcode. The other 926 matrix barcodes remained outside annotation-based outcomes. Matrix barcodes `-1` were mapped to released donor suffixes `.1`–`.8`; this join clarification was frozen after schema inspection and before outcomes.

| Cell type | Donors >=20 nuclei | Young | Aged | Minimum nuclei/donor | Base eligible | Hellinger eligible |
|---|---:|---:|---:|---:|---|---|
| BEC | 8 | 4 | 4 | 79 | True | True |
| EpiC | 8 | 4 | 4 | 21 | True | False |
| GC | 8 | 4 | 4 | 25 | True | True |
| IC | 8 | 4 | 4 | 51 | True | True |
| LEC | 7 | 4 | 3 | 10 | False | False |
| SC | 8 | 4 | 4 | 2983 | True | True |
| SMC | 8 | 4 | 4 | 70 | True | True |
| TC | 6 | 3 | 3 | 2 | False | False |

LEC and TC fail the base all-four-donors/group floor and are untestable for P1, pseudobulk, and CLR. EpiC additionally fails Hellinger QC/reference floors. Hellinger uses 1,500-UMI nuclei, technical-gene exclusions, at least 20 qualifying nuclei per donor, and at least 50 cross-fitted young-reference nuclei for every scored donor.

### Primary P1 transport

| Cell type | Donors | Effect (95% CI); nominal p; BH q |
|---|---:|---|
| BEC | 8 | -0.0003492 [-0.1407, 0.14]; p=0.9951; q=0.9951; matched-null p=0.9005; q=0.9005 |
| EpiC | 8 | 0.03077 [-0.03354, 0.09509]; p=0.2734; q=0.8365; matched-null p=0.607; q=0.9005 |
| GC | 8 | 0.0159 [-0.3089, 0.3407]; p=0.9048; q=0.9951; matched-null p=0.8607; q=0.9005 |
| IC | 8 | 0.01552 [-0.2288, 0.2599]; p=0.8767; q=0.9951; matched-null p=0.7264; q=0.9005 |
| SC | 8 | 0.1151 [-0.1286, 0.3587]; p=0.2788; q=0.8365; matched-null p=0.0199; q=0.1194 |
| SMC | 8 | 0.05347 [-0.1632, 0.2701]; p=0.5538; q=0.9951; matched-null p=0.2736; q=0.8209 |

No primary P1 result survives conventional or matched-null BH. These four-versus-four intervals do not establish equivalence.

### Primary CLR composition

| Cell type | Donors | Effect (95% CI); nominal p; BH q |
|---|---:|---|
| BEC | 8 | -0.4732 [-1.656, 0.7099]; p=0.351; q=0.5265 |
| EpiC | 8 | 1.316 [-0.05809, 2.689]; p=0.0571; q=0.1713 |
| GC | 8 | -0.6141 [-2.491, 1.263]; p=0.4387; q=0.5265 |
| IC | 8 | 0.3079 [-0.1979, 0.8137]; p=0.1784; q=0.3568 |
| SC | 8 | 0.692 [0.009436, 1.374]; p=0.04789; q=0.1713 |
| SMC | 8 | 0.2799 [-0.8156, 1.375]; p=0.5403; q=0.5403 |

### Primary cross-fitted young-reference Hellinger displacement

| Cell type | Donors | Effect (95% CI); nominal p; BH q |
|---|---:|---|
| BEC | 8 | -0.005744 [-0.0389, 0.02742]; p=0.6748; q=0.8434 |
| GC | 8 | 0.01838 [-0.01231, 0.04908]; p=0.1843; q=0.8434 |
| IC | 8 | -0.005066 [-0.01824, 0.008108]; p=0.3683; q=0.8434 |
| SC | 8 | 0.00556 [-0.02574, 0.03686]; p=0.6671; q=0.8434 |
| SMC | 8 | 0.001247 [-0.02125, 0.02374]; p=0.8923; q=0.8923 |

### Donor-pseudobulk genes and sensitivity

| Model | Tested genes across cell types | Nominal p<0.05 | BH q<0.05 |
|---|---:|---:|---:|
| primary_group_depth_adjusted | 100,009 | 6,741 | 3 |
| unadjusted_group_sensitivity | 100,009 | 6,775 | 7 |

Primary BH hits (effect is aged minus young on log1p CPM):

| Cell type | Gene | Effect (95% CI) | p | BH q |
|---|---|---:|---:|---:|
| SMC | CCDC86 | -1.552 [-1.692, -1.413] | 9.829e-07 | 0.04474 |
| GC | TMEM229B | -3.677 [-4.028, -3.326] | 1.32e-06 | 0.04474 |
| IC | AC005702.1 | -2.004 [-2.196, -1.812] | 1.342e-06 | 0.04474 |

The group-only gene sensitivity has its own BH family and is not confirmation of the primary hits. Full effects and intervals are in `results/followup-ovary-jin-gene-models.csv`; direction/FDR concordance is summarized in `results/followup-ovary-gene-sensitivity.csv`.

Same-nucleus RNA–ATAC coupling is impossible because the assays used separate nuclei and BioSamples. Donor-matched linkage is also blocked: no validated released peak–gene link table was identified, and ATAC EC was not silently equated to RNA BEC/LEC. No ATAC matrix or fragment file was downloaded.

## GSE329417 postmenopausal cohort

### Atlas, corrected counts, and donor aggregation

The final annotated atlas contains **439,011 nuclei**, 18,099 features, 28 donors, 64 tissue-piece samples, and 21 coarse labels. Exact fields were donor `SenNet_ID`, sample `unique_sample_ID`, and coarse annotation `cell_type_coarse1`. Atlas sample IDs were reconciled one-to-one through the 64 released CellBender filename stems and their `_p1`–`_p4` suffixes, then to the processing map and 28-donor clinical table; no fuzzy join was used.

All **877,413,549 stored `layers/raw_counts` values** are finite, nonnegative integers. Atlas X was identical to `layers/log1p_norm` in every streamed chunk and was never treated as raw counts. Only CellBender-corrected, unnormalized `layers/raw_counts` entered endpoints. Original uncorrected Cell Ranger matrices are not separately released; therefore an uncorrected sensitivity **cannot be performed** and was not fabricated.

Counts and per-cell summaries were aggregated over all pieces to one donor-cell-type row. SN_0086 samples 3_4 and 16_4, which shared a digested aliquot, were combined inside that single donor aggregation. All 64 samples were used but never treated as 64 independent biological replicates.

| Clinical state | Donors |
|---|---:|
| Benign | 18 |
| Malignant | 6 |
| Pre-malignant | 4 |

Non-ovarian abnormal state is a covariate, not a claim that ovarian pathology was present; the study excluded ovarian neoplasia and endometriosis.

### Eligibility and identifiability

| Cell type | Donors observed | Donors >=20 nuclei | Total nuclei | Base eligible |
|---|---:|---:|---:|---|
| ArtEnd | 28 | 28 | 4,007 | True |
| ArtSMC | 28 | 28 | 28,547 | True |
| B Cell | 10 | 0 | 11 | False |
| CapEnd | 28 | 28 | 4,474 | True |
| Epithelium | 23 | 15 | 1,384 | False |
| LympEnd | 28 | 26 | 1,438 | True |
| Macrophage | 28 | 28 | 2,590 | True |
| Mast | 5 | 0 | 19 | False |
| Mono | 14 | 0 | 32 | False |
| Myofibroblast | 26 | 18 | 4,854 | False |
| NK Cell | 27 | 5 | 354 | False |
| Pericyte | 28 | 28 | 1,981 | True |
| Plasma B | 22 | 1 | 154 | False |
| Schwann | 27 | 20 | 1,404 | True |
| SpecSMC | 27 | 15 | 5,943 | False |
| Steroidogenic Stromal | 25 | 5 | 605 | False |
| Stromal 1 | 28 | 28 | 371,998 | True |
| T Cell | 28 | 18 | 1,587 | False |
| VenEnd | 28 | 24 | 4,497 | True |
| VenSMC | 26 | 14 | 2,895 | False |
| cDC | 26 | 2 | 237 | False |

All eligible primary/sensitivity designs were full rank with residual degrees of freedom.

### Primary P1 transport, exact age per 10 years

| Cell type | Donors | Effect (95% CI); nominal p; BH q |
|---|---:|---|
| ArtEnd | 28 | 0.06447 [-0.05218, 0.1811]; p=0.2626; q=0.8364; matched-null p=0.004975; q=0.02239 |
| ArtSMC | 28 | 0.01239 [-0.08262, 0.1074]; p=0.7884; q=0.8364; matched-null p=0.4229; q=0.5565 |
| CapEnd | 28 | 0.01654 [-0.1095, 0.1426]; p=0.7872; q=0.8364; matched-null p=0.5572; q=0.602 |
| LympEnd | 26 | 0.04098 [-0.07553, 0.1575]; p=0.4695; q=0.8364; matched-null p=0.08955; q=0.1881 |
| Macrophage | 28 | 0.03355 [-0.06571, 0.1328]; p=0.4889; q=0.8364; matched-null p=0.07463; q=0.1881 |
| Pericyte | 28 | -0.01203 [-0.132, 0.1079]; p=0.8364; q=0.8364; matched-null p=0.602; q=0.602 |
| Schwann | 20 | 0.02297 [-0.111, 0.1569]; p=0.7153; q=0.8364; matched-null p=0.4328; q=0.5565 |
| Stromal 1 | 28 | 0.1032 [-0.01054, 0.2169]; p=0.07298; q=0.6568; matched-null p=0.004975; q=0.02239 |
| VenEnd | 24 | 0.03732 [-0.1273, 0.2019]; p=0.6372; q=0.8364; matched-null p=0.1045; q=0.1881 |

Conventional P1 BH hits: 0. Expression-matched-null P1 BH hits: ArtEnd (q=0.02239), Stromal 1 (q=0.02239). These are distinct calibrations; a matched-null hit is not relabeled as a conventional age-associated P1 finding.

### Primary CLR composition, exact age per 10 years

| Cell type | Donors | Effect (95% CI); nominal p; BH q |
|---|---:|---|
| ArtEnd | 28 | -0.065 [-0.3546, 0.2246]; p=0.6447; q=0.8289 |
| ArtSMC | 28 | 0.04899 [-0.3524, 0.4504]; p=0.8017; q=0.8569 |
| CapEnd | 28 | -0.04064 [-0.5046, 0.4234]; p=0.8569; q=0.8569 |
| LympEnd | 26 | -0.2418 [-0.5319, 0.0483]; p=0.09694; q=0.6652 |
| Macrophage | 28 | 0.06806 [-0.1318, 0.2679]; p=0.4857; q=0.8189 |
| Pericyte | 28 | -0.1434 [-0.3421, 0.0553]; p=0.1478; q=0.6652 |
| Schwann | 20 | -0.133 [-0.4633, 0.1973]; p=0.3975; q=0.8189 |
| Stromal 1 | 28 | 0.1252 [-0.2998, 0.5501]; p=0.5459; q=0.8189 |
| VenEnd | 24 | -0.2243 [-0.752, 0.3035]; p=0.381; q=0.8189 |

### Donor-pseudobulk genes and prespecified sensitivity

| Model | Tested genes across cell types | Nominal p<0.05 | BH q<0.05 |
|---|---:|---:|---:|
| primary_depth_plex_clinical | 119,819 | 7,404 | 1 |
| technical_depth_plex_sensitivity | 119,819 | 6,645 | 0 |
| unadjusted_sensitivity | 119,819 | 7,957 | 0 |
| secondary_primary_plus_parity | 119,819 | 8,994 | 0 |

The primary family has 1 BH hits. The 20 smallest nominal p-values among them are shown for orientation; the complete non-cherry-picked table is `results/followup-ovary-primary-gene-fdr-hits.csv`.

| Cell type | Gene | Effect/decade (95% CI) | p | BH q |
|---|---|---:|---:|---:|
| ArtSMC | GABRB2 | 0.5388 [0.387, 0.6905] | 3.77e-07 | 0.04517 |

| Sensitivity model | Primary hits with same direction | Same-direction fraction | Also q<0.05 in sensitivity |
|---|---:|---:|---:|
| technical_depth_plex_sensitivity | 1/1 | 1 | 0 |
| unadjusted_sensitivity | 1/1 | 1 | 0 |
| secondary_primary_plus_parity | 1/1 | 1 | 0 |

Hellinger displacement is **untestable** for every coarse cell type because this postmenopausal cohort has no 20–40-year reference. There was no young-versus-old dichotomization and no borrowing of Jin's reference. Parity appears only in the declared secondary adjustment; age and parity are distinct estimands.

## Untestable arms and interpretation limits

`results/followup-ovary-untestable-arms.csv` contains all 188 blocked/untestable model rows with `outcomes_examined=False`. No failed arm was converted to a zero effect or non-significant result. Non-significance does not demonstrate equivalence, particularly for Jin's eight donors. Cell-level inference, causal language, cross-cohort pooling, and premenopausal extrapolation were prohibited.

## Reproduction and artifacts

```bash
python scripts/prepare_followup_ovary.py
python scripts/fetch_followup_ovary.py --cohort GSE202601
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_jin.py build-cache
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_jin.py analyze
python scripts/fetch_followup_ovary.py --cohort GSE329417
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/prepare_followup_ovary_postmenopausal.py decompress
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/prepare_followup_ovary_postmenopausal.py inspect
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_postmenopausal.py build-cache
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/ovary --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_ovary_postmenopausal.py analyze
python scripts/run_followup_ovary_jin.py finalize
python scripts/run_followup_ovary_postmenopausal.py finalize
python scripts/report_followup_ovary.py
```

Core code: `src/followup_ovary.py`, `scripts/fetch_followup_ovary.py`, `scripts/prepare_followup_ovary.py`, `scripts/prepare_followup_ovary_postmenopausal.py`, `scripts/run_followup_ovary_jin.py`, `scripts/run_followup_ovary_postmenopausal.py`, and this generator. Plans, download hashes, decompressed hash, and count-validation manifests are under `data/manifest/` and `results/followup-ovary-*`. Aggregate tables include endpoint, gene-family, sensitivity, untestable-arm, and resource summaries. Tests are in `tests/test_followup_ovary.py`.

Verification: **149 passed, 1 skipped**. `git diff --check` passed. No commit or push was made. No heart file was modified.
