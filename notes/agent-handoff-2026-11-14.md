# Agescape agent handoff, 2026-11-14

## Purpose and current state

This handoff transfers the follow-up program after completion of the manuscript review copy and all planned follow-up steps except GTEx expression analysis, the blocked Kedlian expression arm, and the final did/did-not report.

- Repository: `/home/reza/projects/agescape-revision`
- GitHub: `https://github.com/Abolfazl-Alipour/agescape`
- Branch: `master`
- Last pushed commit before this handoff: `d07fa18` (`Verify GTEx bulk age-bin analysis design`)
- External data/scratch: `/home/reza/projects/agescape-revision-data`
- Python: `/home/reza/miniconda3/envs/agescape/bin/python`
- Test command: `/home/reza/miniconda3/envs/agescape/bin/python -m pytest -q`
- Last full test result: **155 passed, 1 skipped**, with one existing AnnData `ImplicitModificationWarning` in the synthetic muscle test.
- Last disk check: about **191 GB free** on the 598 GB root filesystem. RAM was about 49 GB. No heavy background jobs were running.
- Worktree was clean before creation of this handoff. Fetch and verify it again before starting.

The complete execution plan is `notes/followup-execution-plan.md`; the 43-dataset suitability review is `notes/followup-dataset-checklist.md`. Follow those documents and the reports cited below rather than starting a new analysis program.

## Immediate next actions

1. **Run GTEx bulk expression analysis, the last planned dataset analysis.** The metadata gate and frozen design are complete in `reports/followup-gtex-metadata.md` and `data/manifest/gtex_v8.json`. No expression matrix has been downloaded and no GTEx outcome has been examined.
2. **Commit and push GTEx with clear effect sizes, uncertainty, multiplicity, sensitivity outcomes, and untestable arms.** If the full frozen model cannot be implemented faithfully, report the blocker explicitly; do not silently substitute a smaller analysis or a different estimand.
3. **Write the complete final did/did-not report**, recommended as `reports/followup-final-summary.md`, and give Reza the same summary in the final response. Include all completed effect sizes, robustness failures, blocked datasets, unattempted datasets, contradictions with the manuscript, test output, and ending commit hash.
4. **Leave Kedlian explicitly blocked unless Reza supplies the raw file.** The exact browser-download requirement is below. Do not mark Step 4 complete without it.
5. **Do not edit the manuscript from this follow-up work without Reza's review.** In particular, SEA-AD narrows the external support for the earlier interneuron result and must be shown to Reza before any manuscript change.

## GTEx: exact starting point

Read `data/manifest/gtex_v8.json` before writing code. Its design was frozen before expression access.

- Public count file: `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz`
- URL: `https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz`
- Expected bytes: **917,091,994**
- TPM is not the primary input.
- Tissues and public donor counts: Whole Blood 755; Muscle - Skeletal 803; Lung 578; Kidney - Cortex 85; Ovary 180; Testis 361.
- Public ages are six ordered decade bins, 20–29 through 70–79. Exact ages are protected and must not be requested or imputed.
- Primary: separate tissue-specific negative-binomial count models with median-of-ratio size-factor offsets and empirical-Bayes dispersion sharing. Model ordinal age-bin order plus sex where estimable, RIN, ischemic time, exonic rate, collection site, isolation protocol, and the offset exactly as specified in the manifest.
- Gene filter: raw count at least 10 in at least `ceil(0.10 * analysis_samples)` samples per tissue.
- Primary multiplicity: **one BH family across every retained tissue-by-gene age Wald test in all six tissues**. Do not replace this with within-tissue correction after seeing outcomes.
- Sensitivities are already frozen: categorical age with pooled 60–79 (kidney ineligible under the ten-per-level rule), excluding 70–79, autolysis/Hardy adjustment, nonnegative blood ischemic time, and 76-bp-only samples. Read the exact rules in the manifest.
- Kidney has one 20–29 and four 70–79 donors; ovary/testis have four/nine 70–79 donors. Report sparse-bin and leverage diagnostics.
- No cell fractions, per-cell drift, Hellinger, entropy, single-cell composition, exact-age slope, or deconvolution-as-measured-composition claim is permitted. Testis bulk cannot separate germ-cell abundance from cell-intrinsic regulation.
- Verify the GCT sample IDs against the 17,382 public RNA-seq IDs, count integrity, gene identifiers/build, model rank, dispersion, and fit before reporting age effects.

Public metadata are already under `/home/reza/projects/agescape-revision-data/gtex-metadata/`. Do not request dbGaP access. Public GTEx downloads do not need browser assistance.

## Resource safety: mandatory

A disk-full scare and an earlier machine shutdown changed the operating procedure. Do not revert to unguarded downloads or whole-atlas in-memory reads.

- Use `scripts/resource_guard.py` for every expression download and heavy job.
- Maintain a **30 GB free-space reserve**, a documented per-file cap, and an **8 GB maximum analysis memory** unless a smaller cap suffices.
- The guard takes a nonblocking exclusive heavy-job lock, uses `systemd-run --user --scope`, sets `MemoryMax`, and sets `MemorySwapMax=0`.
- Run only one heavy job or download at a time. Small metadata verification can proceed separately, but never run two expression jobs together.
- Use backed/direct HDF5 or streamed MatrixMarket reads and small donor batches. Do not initialize a large AnnData object merely to read metadata.
- The 45 GB TS2 all-cells download was abandoned and must not be resumed. TS2 per-organ objects were used instead.
- Check `df -h /`, `free -h`, and active processes before and after heavy work. Do not delete source data or results to clear space without Reza's instruction.

Examples:

```bash
python scripts/resource_guard.py download \
  --url <PUBLIC_URL> --output ../agescape-revision-data/<file> \
  --min-free-gb 30 --max-file-gb <CAP>

python scripts/resource_guard.py run \
  --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- \
  /home/reza/miniconda3/envs/agescape/bin/python <script> <args>
```

`tests/test_resource_guard.py` covers the guard. The resource commit is `58f16ba`.

## Completed follow-up results

### De Man lung, GSE281219/GSE136831 controls

Commit: `e5ae814`. Report: `reports/followup-gse281219.md`.

- Only the 28 IPF controls were analyzed. The 32 Baylor raw matrices lack released Baylor annotations/integrated object, so a full 60-donor reproduction is blocked.
- No primary entropy result had q < 0.05.
- Published-KL gCap slope was -0.002005/year, 95% CI -0.004181 to +0.000170, p=0.0687, q=0.2747. It was negative in all 20 thinning seeds but significant in only 6; this is sampling-sensitive, not a confirmed entropy increase.
- No full-depth or 1,500-UMI Shannon association was significant.

### Tabula Sapiens 2.0

Commits: `983b37e`, `da6b297`, `a335764`, `1c527c0`. Reports: `reports/followup-tabula-sapiens-2.md`.

- Supplementary metadata reconcile 1,136,218 cells, 24 donors, and 28 tissues. Only 22 donors have curated ages; TSP30 has a 26-versus-321 cross-sheet conflict resolved in favor of the curated 26 and recorded explicitly.
- No method-specific organ pair passed both six age-verified shared donors and age SD >=12. The planned within-person multi-organ claim is therefore untestable.
- Skin has four donors and is blocked. The local skin H5AD is old four-donor TS1 data, not a new TS2 cohort.
- Only 10x blood NK cells passed the metadata candidate screen. Seven donors/2,960 cells entered exploratory expression analysis.
- No gene passed BH among 32,830 tested genes; minimum q=0.8876988605.
- P1 effect was +0.0317244765 per 10 years, 95% CI -0.1338723124 to +0.1973212655, conventional p=0.4964. The matched-null empirical p/q was 0.009950, but that is signature specificity relative to random gene sets, not a significant direct age association.
- Hellinger and off-program expression were untestable with one young donor. Residual df was only 2; no equivalence or multi-organ claim is allowed.

### Lake kidney

Commits: `bb90af4`, `9e8c340`. Reports: `reports/followup-lake-kidney.md`, `reports/followup-lake-kidney-analysis.md`.

- GSE183279 is the full processed SuperSeries; GSE183276 is only scCv3. Public processed matrices do not require KPMP/NIDDK authorization.
- Restricted public analyses used 18 scCv3 exact-age adults and 9 primary snCv3 exact-age adults; the 14-year-old donor entered only the declared all-exact-age sensitivity.
- Nine scCv3 and 13 snCv3 cell types were eligible. Assays were analyzed separately.
- No conventional P1 result survived assay-wise BH: minimum q=0.504 scCv3 and 0.904 snCv3.
- Matched-null P1 flagged three cell types per assay, but their conventional donor-model CIs included zero; these are not confirmed replication.
- No CLR result survived BH: minimum q=0.212 scCv3 and 0.954 snCv3.
- No scCv3 Hellinger result survived BH: minimum q=0.301. snCv3 Hellinger was untestable with one eligible young-reference donor.
- Existing SCR/SN objects overlap Lake releases; they are not independent replication cohorts.

### Muscle: Perez complete, Kedlian blocked

Commits: `6763430`, `4cc91e2`, `665fbc2`. Reports: `reports/followup-muscle-metadata.md`, `reports/followup-muscle-analysis.md`, `reports/followup-kedlian-alternate-inputs.md`.

- Perez: 17 donors, 6 young/11 old, 134,486 nuclei; all counts verified. No released barcode-level annotation was found.
- Whole-nucleus donor-pseudobulk models produced 14 old-versus-young and 8 continuous-age BH hits among 20,343 genes. These are **composition-confounded**, not cell-intrinsic results.
- Perez P1, cell-type models, CLR, and Hellinger remain blocked. The separate 72-man bulk arm was not used as single-cell evidence.
- Kedlian explicitly reused Perez; do not pool them or call them independent replications.
- Kedlian expression was not analyzed. Its 4,025,846,616-byte raw object repeatedly downloaded at about 0.1 MB/s with disconnects; only 154,009,600 bytes were acquired. No expression values were read.
- All seven smaller official processed alternatives lack raw/count layers and have noninteger X samples. They are not equivalent substitutes.

Kedlian requires Reza's browser download:

- URL: `https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`
- Destination: `/home/reza/projects/agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`
- Existing partial file is unusable until completed. Do not delete it if planning to resume.
- Frozen pipeline: `scripts/run_followup_kedlian.py`; plan: `data/manifest/muscle_analysis_plan.json`. Run build-cache, analyze, finalize, then `scripts/report_followup_muscle.py` after the complete hash/size gate passes.

### PBMC

Commits: `5e8bc6f`, `8cd58bd`. Reports: `reports/followup-pbmc-metadata.md`, `reports/followup-pbmc-analysis.md`.

- Only GSE213516 (Zhu) and GSE214546 (Thomson) were modeled. GSE157007 was excluded from inference with two young donors and no barcode-label map.
- Strict primary is young <=40 versus old >60. Atlas inclusion of age 60 and exact-age models are separate sensitivities. Strict N=10 Zhu and 23 Thomson.
- CD8-naive CLR decreased in both: Zhu -1.6056, q=0.0150; Thomson -1.4238, q=0.0000647.
- Non-classical-monocyte CLR increased in both: Zhu +0.5402, q=0.0890; Thomson +0.9213, q=0.000908.
- Thomson batch-adjusted targeted composition effects remained directionally consistent, but both had exploratory all-label batch q=0.139044. Do not call the batch sensitivity FDR-significant.
- Targeted CD8-naive genes replicated incompletely: CLIC1, S100A4, CCL5 passed Zhu's six-test BH; no targeted gene passed Thomson's. PRF1 and GZMA were negative in Thomson.
- Broad P1 was cohort-dependent and not batch-robust. Thomson strict survivors: primary 11/13, depth 8/13, CMV 11/13, pool 9/13, but batch 0/13 (minimum batch q=0.109444). Zhu had no strict-primary conventional P1 survivor.
- Strict depth-adjusted Hellinger produced three negative hits in Zhu and none in Thomson. There is no universal positive displacement result.
- MAIT was untestable in both released annotation universes.
- Zhu assembly provenance remains unresolved; its gene results are symbol-level only.

### Cui testis

Commit: `fbdd31d`. Report: `reports/followup-cui-testis.md`.

- The 25 GEO libraries are 25 separate donors; the other 10 paper-atlas donors are external-only. N=25, not 35.
- Donor8 is 40 in the primary figure and 39 in GEO; the primary/sensitivity choice is documented.
- All 75 count files passed integrity checks. Raw matrices contain 148,740 barcodes and 805,486,323 UMIs.
- Paper-threshold QC retains 95,521 barcodes, 3,609 fewer than GEO's 99,130-cell summary; this discrepancy remains unresolved.
- No public donor-and-barcode cell-type annotation crosswalk was found. P1, cell-type pseudobulk, CLR, and Hellinger were **not run**. No whole-testis substitute was used.
- Germ-cell-stage composition can mimic broad displacement; no general aging claim is permitted. Incomplete-spermatogenesis-adjusted models were untestable.

### Ovary

Commits: `0610c34`, `bf46ba2`. Reports: `reports/followup-ovary-metadata.md`, `reports/followup-ovary-analysis.md`.

Jin GSE202601:

- Eight donors, four young/four reproductively-aged; all 42,568 annotations joined to raw matrices.
- Six cell types passed base floors; LEC/TC were untestable.
- No primary P1, CLR, or Hellinger result passed BH.
- Three primary gene hits passed the model-wide BH family: SMC CCDC86 -1.552, GC TMEM229B -3.677, IC AC005702.1 -2.004; each q=0.04474. Two of the three also passed the separate unadjusted sensitivity family; do not treat that as independent replication.
- Ovary5 is 47 in snRNA metadata versus 49 in GEO/snATAC. No continuous-age model was fitted.
- RNA and ATAC used different nuclei. Same-nucleus coupling is impossible; donor-matched linkage remains blocked without validated peak-gene links.

Postmenopausal GSE329417:

- 28 exact-age donors aged 50–84, 64 samples, 439,011 nuclei; donors are the inferential units. SN_0086 technical replicates were combined in donor aggregation.
- Only CellBender-corrected `layers/raw_counts` were used. Original uncorrected Cell Ranger matrices are unavailable, so an uncorrected sensitivity was impossible.
- Nine cell types passed floors. No conventional P1 or CLR result passed BH. Matched-null ArtEnd/Stromal 1 results are not conventional age associations.
- ArtSMC GABRB2 was the sole primary gene hit: +0.5388 log1p-CPM/decade, 95% CI 0.3870–0.6905, q=0.04517. Direction stayed positive in sensitivities, but no sensitivity model passed BH. Do not call it robust established aging biology.
- Hellinger was untestable with no 20–40 reference. The dataset remains a bioRxiv preprint.

### Jia heart

Commit: `2bb297e`. Report: `reports/followup-jia-heart.md`.

- Accession resolved from the primary data-availability statement: Figshare `10.6084/m9.figshare.31939227`, version 1, 132 public files totaling 3,431,641,059 bytes.
- Full atlas is 29 donors/54 LV-RV tissues: 24 Jia donors plus five external HCA donors.
- Adult exact continuous-age analysis is blocked: only six internal postnatal snRNA donors have exact ages; AG6 has five-year bands.
- AG4/AG5 versus AG6 is completely confounded by internal/external source. No clean adult group claim is possible.
- Jia barcode-level annotations were not found. HS246/HS247 sex mapping conflicts.
- No expression was downloaded. Do not substitute another heart dataset or manufacture a descriptive result from source-confounded groups.

### SEA-AD interneuron extension

Commit: `5ca2735`. Report: `reports/followup-sea-ad-interneurons.md`.

- New donor-pseudobulk analysis only; no repeat of existing P1/signature/composition work and no Hellinger/reference relaxation.
- 15 no-dementia consensus controls with Not AD/Low ADNC, exact ages 72–102, 8 female/7 male. Narrow age SD=9.35 fails the general 12-year power-floor preference and limits inference.
- No gene passed SEA-AD subtype-wide FDR among 1,923 VIP, 1,831 pvalb, or 1,890 sst genes.
- Frozen both-FDR set replication: VIP 29/71 same direction, block p=0.8921; pvalb 19/64, p=0.9620; sst 86/240, p=1.0000. All three-subtype BH q=1.0; zero same-direction SEA-AD-FDR genes.
- Primary like-for-like OLS genome-wide correlations: VIP r=0.0054, 95% block-bootstrap CI -0.0777 to +0.1003; pvalb r=-0.1775, CI -0.3078 to -0.0440; sst r=0.0390, CI -0.0607 to +0.1179.
- Results were estimator-, assay-, and donor-sensitive. Leave-one-out was unstable. VIP discovery-only uncorrected p=0.0240 is not robust replication.
- Mouse-derived P1 coverage was only 4 VIP, 3 pvalb, and 6 sst genes, so no mouse-P1 verdict is possible.
- **Manuscript implication:** SEA-AD does not independently support the earlier DLPFC interneuron program. Describe it as cohort-specific with limited external generalizability when Reza decides how to revise; it does not invalidate the larger original MSSM/HBCC result and does not prove no aging.

### GTEx

Commit: `d07fa18`. Metadata only; expression not started. See the immediate next actions above.

## Earlier core revision facts relevant to the final report

The manuscript review copy is complete at `report/plos-one/final-plos-one.tex` and `.pdf`, titled **“No uniform transcriptional response detected in mouse and human aging”**. At completion it had 20 pages, a 267-word abstract, four figures, and 50 references. Follow-up analyses have not been incorporated and the manuscript must not be edited without Reza's review.

Previously verified core results:

- MSSM matched transport null: no survivor among the original three global transport survivors.
- Mouse Scrublet-on within-tissue concordance p=0.0440; cross-tissue p=0.1528. Human p=0.7003; fly p=0.001; worm p=0.003; zebrafish p=0.3606; naked mole rat exact p=0.70.
- Original interneuron replication: VIP r=0.7217/76 both-FDR genes; pvalb r=0.6839/66; sst r=0.7614/249; block p=0.001 each. Other brain cell types also change, so no exclusivity claim. SEA-AD now limits its external generalizability.
- Scrublet-on versus off was immaterial across all 16 mouse tissues.

Do not edit `report/manuscript-aging-mechanisms.md`. Do not submit to a journal, contact reviewers/repositories for access, or make anything public beyond the GitHub pushes Reza has already authorized.

## Analyses not attempted or explicitly excluded

- No unverified liver, SenNet, gut, retinal, HRCA, superior-frontal-gyrus, or GSA-Human data were used.
- Tiny C-grade skin/adipose/other datasets were not used to manufacture general human conclusions.
- No controlled raw sequencing, exact GTEx ages, FASTQ, new experiments, longitudinal claims, or causal claims.
- No same-nucleus Jin RNA-ATAC coupling, postmenopausal Hellinger, clean Jia adult aging model, Cui cell-resolved model, or Perez cell-type model was possible from released inputs.
- No TS2 all-cells or multi-organ analysis, no Baylor De Man annotation analysis, and no completed Kedlian expression analysis.
- No final follow-up manuscript integration, web-app update, or overall final did/did-not report has been done yet.

## Git and verification procedure

1. Start with:

```bash
cd /home/reza/projects/agescape-revision
git fetch origin master
git status --short --branch
git rev-parse HEAD origin/master
df -h /
free -h
```

2. If `origin/master` advanced, review and merge it before new work. Do not force-push or overwrite another agent's work.
3. Run one analysis step at a time. Commit and push each metadata gate, result, and report separately with clear messages.
4. Run focused tests and the full suite on the delivered work. Also run `git diff --check` and manifest/path security checks.
5. Every number must trace to a committed aggregate result file. Do not report downloaded metadata as biology or non-significance as equivalence.
6. Do not store credentials in the repo, remote URL, logs, reports, or commits. Reza supplied a PAT earlier; if available from the authorized prior session, use it only in a one-time push URL with this form and never print it:

```bash
git push https://x-access-token:<PAT>@github.com/Abolfazl-Alipour/agescape.git HEAD:master
```

If the credential is unavailable, ask Reza rather than disabling authentication or inventing a workaround.

## Final response Reza is still owed

After GTEx and any feasible blocked-arm resolution, report:

- What was completed, with effect sizes and confidence intervals.
- What survived robustness checks and what died or was untestable.
- Every blocked dataset and exact reason, especially Kedlian.
- What was not attempted.
- Every contradiction or narrowing of manuscript claims, especially SEA-AD versus the earlier interneuron replication and PBMC batch sensitivity.
- Full test output, disk/resource status, and ending commit hash.

Do not present this handoff or the metadata gates as completed biology.
