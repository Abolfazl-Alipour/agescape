# Kedlian et al. 2024 Skeletal Muscle (`E-MTAB-13874`): Documented Download & Analysis Blocker

## 1. Status and Decision

- **Status**: `BLOCKED` (sustained server throughput ~0.25 MB/s over HTTPS on the required 4.03 GB raw-count `.h5ad`; single bounded resume check terminated per overnight budget; continuous-age analysis also metadata-blocked because only 5-year age brackets are released publicly).
- **Evidence files**:
  - [Resume attempt telemetry JSON](../results/followup-kedlian-resume-attempt.json)
  - [Prior guarded download telemetry JSON](../results/followup-muscle-kedlian-download-resource.json)
  - [Alternate processed-object audit](../reports/followup-kedlian-alternate-inputs.md)
  - [Primary muscle metadata audit](../reports/followup-muscle-metadata.md)

---

## 2. Single Bounded Resume Attempt (`2026-11-15`)

Per instructions allowing **one bounded resume attempt only**, `scripts/resource_guard.py download` was invoked once against the existing partial file (`154,009,600` bytes, `3.83%` from `2` prior attempts) with HTTP 206 overlap verification (`--verify-resume-overlap`) and hard guardrails (`--max-time 60`, `--speed-limit 1048576` [1.0 MB/s], `--speed-time 25`, `--min-free-gb 30`, `--max-file-gb 5`, `--min-size-bytes 4025846616`):

| Metric | Value | Source |
|---|---:|---|
| Target asset | `SKM_human_raw_cells2nuclei_2023-06-22.h5ad` | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Expected size | `4,025,846,616` bytes (`4.03 GB`) | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Bytes before resume attempt | `154,009,600` (`3.83%`) | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Bytes after resume attempt | `212,157,824` (`5.27%`) | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Bytes transferred in attempt | `+58,148,224` bytes (`+58.15 MB`) | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Elapsed duration (including curl retries) | `229.41 s` | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Effective average throughput | `0.2535 MB/s` | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Exit code and curl message | `rc = 28` (`Operation too slow. Less than 1048576 bytes/sec transferred the last 25 seconds`) | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |
| Completed | `false` | [Resume attempt JSON](../results/followup-kedlian-resume-attempt.json) |

At `0.2535 MB/s`, fetching the remaining `3,813,688,792` bytes would take ~4.2 hours of continuous transfer even without connection drops. Per instructions, **no further automated downloads were attempted**.

---

## 3. Why Smaller Processed Atlas Objects Cannot Substitute

As verified by schema byte-range checks in [the alternate-inputs audit](../reports/followup-kedlian-alternate-inputs.md) and cataloged in [the muscle metadata report](../reports/followup-muscle-metadata.md), the smaller CellGeni processed `.h5ad` downloads (`SKM_human_pp_cells2nuclei_2023-06-22.h5ad`, `SKM_MuSC_human_2023-06-22.h5ad`, `SKM_myonuclei2myofiber_human_2023-06-22.h5ad`, `SKM_myonuclei_human_2024-04-14.h5ad`, `SKM_immune_human_2023-06-22.h5ad`, `SKM_fibroblasts_Schwann_human_2023-06-22.h5ad`, and `SKM_Endothelium_SMC_human_2023-06-22.h5ad`) hold `float32` non-integer processed values rather than verified integer raw counts. Analyzing non-integer normalized values would violate the repository raw-count integer gate ([scripts/run_followup_kedlian.py](../scripts/run_followup_kedlian.py) requires the `4,025,846,616`-byte raw object).

---

## 4. Cohort Structure and Additional Metadata Limits (`E-MTAB-13874`)

From primary BioStudies SDRF / Supplementary Table 1 verification in [the muscle metadata report](../reports/followup-muscle-metadata.md):
- **Donors & libraries**: `17` unique intercostal-muscle organ donors (`10` male, `7` female; `8` young across `17-20` to `35-40` age bands; `9` old across `50-55` to `70-75` age bands), `43` libraries total (`21` scRNA covering `12` donors [5 young, 7 old]; `22` snRNA covering `13` donors [6 young, 7 old]).
- **No exact chronological ages publicly released**: SDRF and Supplementary Table 1 release only 5-year age brackets (`17-20`, `20-25`, `25-30`, `35-40`, `50-55`, `55-60`, `60-65`, `65-70`, `70-75`), blocking continuous age-per-decade modeling; only a discrete young-vs-old categorical contrast would be valid.
- **Perez overlap caveat**: Kedlian et al. explicitly reuse Perez et al. (`GSE167186`, already analyzed in [the muscle analysis report](../reports/followup-muscle-analysis.md)), and anonymized public IDs do not permit verifying cross-study donor disjointness.

---

## 5. Manual User Unblock Action (Optional)

To unblock a discrete young-vs-old re-analysis of `E-MTAB-13874`, download the raw `.h5ad` (`4,025,846,616` bytes) via browser or external downloader from:
- URL: `https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`
- Target path (relative to repo root): `../agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`

Once placed with exact byte size `4,025,846,616`, verify integer counts and run `scripts/run_followup_kedlian.py`.
