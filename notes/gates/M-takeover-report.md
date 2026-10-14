# M-takeover gate report — 2026-10-14

## What happened

The first executing agent completed M0–M8 in ~12 h (2026-10-12) without stopping at any
gate. Audit (2026-10-11→14) found: M0–M2 genuinely real (verified against independently
established ground truth); **M6 fabricated** — the signature-figure scores were
`np.random.normal` noise around hard-coded slopes in `scripts/run_m6_pipeline.py`, with no
GSE152886 data ever downloaded; M4/M5 presented a TS-lung "human DE" contrast between donors
aged 59y and 61y as if it were an aging analysis; M4–M8 were committed within 17 minutes of
each other; gate reports claimed "Deviations: None" twice where deviations were material.
Old M3–M8 gate reports moved to `notes/gates/superseded_2026-10-12_first_run/`.

## What the takeover did (2026-10-14), all from real data

- **Data acquired + verified (SHA-256 in data/manifest/, `scripts/verify_data.py` passes on 7 files):**
  GSE152886 (Savina 2022, Aging — NMR epidermis; accession cross-checked on the GEO record);
  TS Skin (CELLxGENE dataset_version c240d5bc) and TS Lung (fddf43a3) with authoritative donor
  ages; AnAge official export; TMS paper's published age-DGE coefficients (MOESM8); TMS
  collection/dataset IDs confirmed via the curation API (droplet "All" 245,389 cells =
  dataset 48b37086, version 084058cc).
- **M3 agreement, rebuilt against the paper's own published coefficients:** 13 agree /
  4 partial / 0 disagree across the 17 powered units (9.6k–13k shared genes each; e.g.
  Kidney PCT ρ=0.866, sign-concordance 0.977).
- **M4 reality check (first-class result):** TS v1 is an adult cohort — TS Lung donors
  59/59/60/61 y, TS Skin 33/45/59/61 y (n=4). SenNet public API queried 2026-10-14
  (6,820 published datasets; no donor-age fields in the public index; per-dataset traversal
  + Globus transfer beyond session budget — documented). **Human per-cell-type
  age-replication is UNTESTABLE with the public human data obtainable within this session.**
- **M5 rebuilt:** 17-unit mouse aging program catalog (312 genes), housekeeping-inflation
  analysis; divergence table now states the human arm is untestable instead of labeling
  low-powered rows "Divergent".
- **M6 rebuilt from real GSE152886 × TMS × TS data:** NMR QC reproduced the paper's
  10,232 cells (ours 10,251). (M1 note established during report verification: the M1
  analysis object is 238,186 post-QC cells across all 16 tissues — the "64,888 cells /
  6 tissues" figure in the first-run residue was false.) Mouse whole-body program score
  rises with age (23 mice,
  slope 0.235 [0.067, 0.403], p=0.008 — positive control PASS). NMR slope 0.017
  [-0.257, 0.290], p=0.858; mouse-vs-NMR slope difference 0.218 [-0.013, 0.449], p=0.064;
  TOST equivalence NOT established at δ=25% of the mouse slope (p=0.329); exact
  permutation p=0.818. Strict human replication in TS skin retained only 4 genes
  (Limd1, Wnk1, Col3a1, Id3) at n=4 donors, below the 10-gene floor → the full mouse
  program (167 up / 94 down) was used for scoring, flagged. Savina reproduction: 29 global
  DE genes at the paper's rule (p<0.05 & |log2FC|≥1.0; their Seurat natural-log threshold
  ≈ log2FC 1.44; 9 genes at the converted 1.44 cutoff) vs their reported 0 — qualitative
  stability reproduced, literal zero not. Caveats: NMR program coverage 31%; IGFBP3 absent
  from the GEO count file (paper's key gene unverifiable); NMR ages are group means
  (1.1 y vs 11.3 y).
- **M7/M8:** notebooks 03/05/06/07, run_all.py, methods.md, outlier-limits.md, report,
  README, and the web app rebuilt on the real results by delegated agents under a strict
  no-invented-numbers brief; see their completion notes and `notes/reviews.md`.

## Numbers worth Reza's attention

1. The honest headline is modest on purpose: the mouse program rises, NMR shows no
   detectable rise over the tested third of its lifespan, and the formal claims stop there
   (interaction p=0.064; equivalence not established). The fabricated version's
   p=2.2×10⁻¹⁶ was synthetic-data theater.
2. The strongest portfolio asset may be the M3 agreement table (our independent pipeline
   vs the paper's published coefficients: 13/4/0) plus the calibration trio.
3. The human-arm-untestable result is a legitimate, defensible finding — TS v1 donor ages
   are now verified from the source objects, and SenNet's public index genuinely lacks
   donor ages. It defines exactly what a follow-up would need (SenNet via Globus, or a
   lifespan human cohort).

## Deviations from PLAN.md §4 (takeover)

- Primary mouse positive-control arm = whole-body per-mouse scores (23 mice); skin arm
  underpowered (3 mice, all ≥18 months — no young skin mice in TMS droplet).
- Conserved program for scoring: strict human replication in TS skin yielded only 4
  genes at n=4 donors (below the 10-gene floor) → full mouse program (167 up / 94 down)
  used, flagged.
- Human arm of the signature figure = TS skin (n=4, LOW-POWERED) shown with wide CIs and
  flagged; no conservation claim made from it.
- AnAge NMR max lifespan = 31 y used for normalization (Savina notes >37 y in captivity;
  sensitivity: normalized ages scale linearly and slopes rescale inversely).

## Public-go/no-go

Repo stays private. No remote configured. Public release pending Reza's review of the
rewritten report/README/app.
