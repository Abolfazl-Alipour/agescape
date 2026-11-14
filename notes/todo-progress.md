# Agescape Revision To-Do Progress Tracker

Last updated: 2026-11-14
Current item: Item 4 (IN_PROGRESS)
Manuscript state: `report/plos-one/final-plos-one.tex` and `report/plos-one/final-plos-one.pdf` are FROZEN (no edits allowed).

## Step 0 Environment and Git Verification

- Working directory: `/home/reza/projects/agescape-revision`
- Remote: `https://github.com/Abolfazl-Alipour/agescape.git`, branch `master`
- Starting HEAD / `origin/master`: `40098292bf6aff04cf94043ab8760eb6a4bfdae7` (`Add follow-up agent handoff`)
- Step 0 local commit: `dae5a98` (`Item 0: lock to-do list tracker`)
- Push status: `BLOCKED awaiting GitHub PAT from Reza` (`git push` prompted for interactive username and `credential.helper=store` failed with `remote: Invalid username or token. Password authentication is not supported for Git operations.`)
- Disk (`df -h /`): 598 GB total, 383 GB used, 191 GB available (67% used; well above 30 GB reserve)
- Memory (`free -h`): 37 GiB total, 36 GiB available, 71 GiB swap free
- Python: `/home/reza/miniconda3/envs/agescape/bin/python` (3.11.17)
- Latest test suite: `pytest -q` -> `157 passed, 1 skipped, 1 warning`
- Step 0.1 read-only reads completed: `notes/agent-handoff-2026-11-14.md`, `notes/followup-execution-plan.md` (suspended context only), `report/plos-one/final-plos-one.tex` (frozen)

## 9-Item To-Do List Status

| Item | Description | Status | Evidence / Output Files | Commit Hash |
| --- | --- | --- | --- | --- |
| Item 1 | Cell-type fractions in fly, worm, zebrafish (CLR pipeline + propeller/limma + global BH + germline caveat) | DONE | `results/rev3-composition-fly-worm-zebrafish.csv`, `results/rev3-composition-global-bh-624.csv`, `results/rev3-composition-fly-worm-zebrafish-REPORT.md` | `0485256` |
| Item 2 | Verify LACA multi-tissue naked mole rat aging atlas listing (1.2M cells, young to geriatric) | DONE | `results/rev3-laca-nmr-verification-REPORT.md` | `33da505` |
| Item 3 | Audit only: list every young vs mature comparison in the manuscript with file and line references | DONE | `results/rev3-young-vs-mature-audit-REPORT.md` | pending commit |
| Item 4 | Shannon entropy on mouse lung only, mice 18 months and older only (De Man et al. 2026 Methods comparison) | IN_PROGRESS | `results/rev3-entropy-mature-only-mouse-lung.csv`, `results/rev3-entropy-mature-only-mouse-lung-REPORT.md` | pending |
| Item 5 | Finish De Man et al. 2026 human lung replication (32 Baylor lungs in GSE281219; 28 IPF controls already done) | PENDING | `results/rev3-deman-replication-human-lung.csv`, `results/rev3-deman-replication-human-lung-REPORT.md` | pending |
| Item 6 | Bayesian analysis of whole-vector drift, mature vs old only (mouse >= 18m vs 18m mature reference) | PENDING | `report/cross-species-bayesian-vector-mature-only.md`, `results/rev3-bayesian-vector-mature-only-posterior.csv`, `results/rev3-bayesian-vector-mature-only-diagnostics.csv` | pending |
| Item 7 | Bayesian analysis for the entire paper, report only across 5 families (composition, transport, concordance, off-program/entropy, interneurons) | PENDING | `report/cross-species-bayesian-full.md`, `results/rev3-bayesian-full-summary.csv`, `results/rev3-bayesian-full-diagnostics.csv` | pending |
| Item 8 | Randomness / variability test in mouse Tabula Muris Senis (cross-sectional 18m vs 30m within/between variance + permutation null) | PENDING | `results/rev3-randomness-variability-test.csv`, `results/rev3-randomness-variability-test-REPORT.md` | pending |
| Item 9 | Zebrafish testis germline caveat for Limitations section | DO_NOT_DO | Unconfirmed by Reza; frozen per instructions | N/A |

## Suspended Work (Do Not Start)

- GTEx bulk expression (`data/manifest/gtex_v8.json`, `reports/followup-gtex-metadata.md`): SUSPENDED
- Kedlian muscle download (`data/manifest/muscle_analysis_plan.json`, `scripts/run_followup_kedlian.py`): BLOCKED / SUSPENDED
- Final did/did-not report (`reports/followup-final-summary.md`): SUSPENDED
