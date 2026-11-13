# Follow-up dataset checklist

Status meanings: **RUN** = suitable now after metadata checks; **LIMITED** = only the stated restricted analysis; **BLOCKED** = do not use until the stated issue is resolved; **SKIP** = not appropriate for Agescape. Every RUN still requires raw-count and donor/sample verification before modeling.

| # | Dataset | Status | Exact action |
|---:|---|---|---|
| 1 | Tabula Sapiens 2.0, GSE306755 | RUN, Step 2 | Within-person organ interactions and skin after donor/organ crosswalk. 10x primary; Smart-seq separate. No proteomic organ-age claim. |
| 2 | Tabula Sapiens 1.0 | SKIP | Superseded by 2.0 and already represented. No new analysis. |
| 3 | GTEx v8 bulk | LIMITED, Step 10 | Bulk tissue gene-level age-bin models only. No cell fractions, per-cell drift, or single-cell claims. Exact ages require Reza's dbGaP authorization. |
| 4 | GTEx single-nucleus pilot | SKIP | Only 16 donors and ages unverified. No analysis until ages are verified. |
| 5 | HCATA integration, Bartz 2025 | LIMITED | Discovery index only. Do not treat as a dataset or quote the unverified 92M-cell figure. |
| 6 | Sole-Boldo skin, GSE130973 | LIMITED | Five male donors only. Young/old or descriptive replication only, no continuous age slope. |
| 7 | Zou eyelid skin, HRA000395 | BLOCKED | Controlled GSA terms and grouped ages. No download or modeling. |
| 8 | Dermal fibroblasts, GSE307377 | SKIP | Cultured cells and bulk RNA-seq. Not suitable for single-cell aging claims. |
| 9 | PsychAD lifespan DLPFC, syn52396927 | DONE | Already used. Do not duplicate; retain as primary brain cohort. |
| 10 | PsychAD full cross-disorder, syn60084804 | DONE/LIMITED | Already used through the neurotypical subset. Disease-mixed donors require control filtering. |
| 11 | SEA-AD, Gabitto 2024 | RUN, Step 9 | Interneuron donor-pseudobulk replication with diagnosis/pathology filters. No primary young-reference Hellinger without author-approved new rule. |
| 12 | Mathys 2019 ROSMAP, syn18485175 | SKIP | Controlled, AD-mixed, ages unverified. No analysis. |
| 13 | Anderson DLPFC multiome, GSE214637 | LIMITED | Fifteen-donor chromatin/RNA coupling pilot only. No aging slope or universal-aging claim. |
| 14 | Wang superior frontal gyrus, 2024 | BLOCKED | Accession, ages, and access unverified. No download. |
| 15 | Filippov PBMC integration | RUN, Step 5 | Only GSE157007, GSE213516, and GSE214546. Cohort-specific young/old or age models and replication. GSA/Synapse arms excluded. |
| 17 | Supercentenarian T cells, OUKA | SKIP | T-cell-only extreme-age groups. Not suitable for general aging analysis. |
| 18 | Kedlian muscle, E-MTAB-13874 | RUN, Step 4a | Seventeen-donor group contrast with exact-age verification. Cross-species muscle replication, not a powered continuous-age study. |
| 19 | Perez muscle, GSE167186 | RUN, Step 4b | Separate small snRNA group contrast. Verify exact ages and released filtered count. Bulk male-only caveat. |
| 20 | Reconstructive-surgery muscle, GSE143704 | BLOCKED | Publication and metadata unverified. No analysis. |
| 21 | Jia lifespan heart atlas | BLOCKED, Step 8 | Do not use until accession is found in the paper's data-availability statement. |
| 22 | Kanemaru heart, 2023 | BLOCKED | Donor count/age conflicts and accession unresolved. No analysis. |
| 23 | Gao left-ventricle integration, 2026 | BLOCKED | Accession, component cohorts, and ages unverified. No analysis. |
| 24 | Rahman healthy liver, 2026 | BLOCKED | Accession, counts, and ages unverified. No analysis. |
| 25 | SenNet liver multiome | BLOCKED | Repository ID and terms unresolved. No download. |
| 26 | Lake kidney atlas | BLOCKED pending access, Step 3 | Verify full accession and ages first. Reza must handle KPMP/NIDDK authorization. If obtained, run signature, CLR, and eligible Hellinger separately from existing kidney cohorts. |
| 27 | Sikkema HLCA | DONE | Already used as integrated lung cohort. No duplicate download or analysis. |
| 28 | Xu lung parenchyma subset | DONE/OVERLAP | Subset of HLCA already represented. Only a distinct clock analysis would justify new work. |
| 29 | De Man lung aging, GSE281219 + GSE136831 controls | RUN FIRST, Step 1 | Exact published entropy reproduction plus direct Shannon entropy, then eligible P1/Hellinger/composition/off-program. Baylor and IPF controls separately for inference. |
| 30 | Raredon lung 2019, GSE133747 | LIMITED | Nine donors only. Optional targeted replication, not age-slope estimation. |
| 31 | Travaglini lung 2020 | SKIP | Three donors. No aging inference. |
| 32 | Zhou adipose, HRA001279 | BLOCKED | Donor counts, ages, and GSA terms unresolved. No analysis. |
| 33 | Adipose-derived stem cells, GSE276163/GSE276164 | SKIP | Six cultured bulk samples. Not suitable for single-cell aging. |
| 34 | Burclaff gut, 2022 | SKIP | Three donors and accession unresolved. No gut-aging inference. |
| 35 | Enge pancreas, GSE81547 | SKIP | Eight donors, maximum age 54, developmental span. Not an aging cohort. |
| 36 | Adelman HSC-enriched, GSE104406 | SKIP | Bulk sorted-population RNA. No single-cell analysis. |
| 37 | HSC expansion, GSE302126 | BLOCKED | Exact ages and donor metadata unresolved. No analysis. |
| 38 | Yan retina, GSE148077 | BLOCKED | Seven donors and ages unresolved from primary record. No analysis. |
| 39 | Yang retinal aging, 2026 | BLOCKED | Accession unresolved and diabetic retinopathy confounding. No analysis. |
| 40 | Butler RPE, GSE159435 | SKIP | Bulk RNA-seq and small N. Not suitable for single-cell aging. |
| 41 | HRCA integration, Li 2026 | BLOCKED | Accession and ages unresolved. Use only as a literature index for now. |
| 42 | Jin ovary multiome, GSE202601 | RUN, Step 7a | Eight-donor young/aged group contrast and exploratory RNA-chromatin coupling. No powered continuous-age claim. |
| 43 | Postmenopausal ovary, GSE329417 | RUN, Step 7b | Verify demultiplexed metadata. Within-postmenopausal age association only if exact ages support it; otherwise descriptive QC only. |
| 44 | Cui testis, GSE254315 | RUN, Step 6 | Resolve 25 samples versus 35 donors first. Pool-level groups only unless numeric pool ages are defensible. Germline-composition caveat mandatory. |

Entry 16 in the source catalog is a pointer, not a dataset, so it is intentionally omitted from this 43-row checklist.

## Execution order

1. De Man lung
2. Tabula Sapiens 2.0
3. Lake kidney if Reza obtains access, otherwise record blocked and continue
4. Kedlian then Perez muscle
5. Three open PBMC GEO cohorts
6. Cui testis after mapping
7. Jin then postmenopausal ovary
8. Jia heart only if accession is verified
9. SEA-AD interneuron extension
10. GTEx bulk last
