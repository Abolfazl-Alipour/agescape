# M9 Data Scoping — Human Lifespan Single-Cell Sources

Date: 2026-10-21. Status: scoping complete; downloads/analysis pending where marked.

## The question M9 must answer

Does the 261-gene mouse aging program (or the cell-type-private programs) track donor age
in human single-cell data with real per-donor ages? Chapter 5 of the app currently says the
human arm is untestable with public data (Tabula Sapiens: 4 donors). M9 replaces that with
a powered answer, or with an honest, quantified account of what remains missing.

## Candidate sources, ranked

### 1. DLPFC Lifespan Atlas — Yang et al., Nature 2026 (657:1003-1015), s41586-026-10271-7

- **Cohort:** 284 neurotypical donors, ages 0-97, 1,307,674 nuclei, human DLPFC snRNA-seq.
  Single integrated cohort; 8 classes / 27 subclasses / 65 subtypes; donor-level pseudobulk
  already the paper's analysis unit. Replication set: 306 donors, 2.04M nuclei.
- **Why it is the front-runner:** exactly the design the M4 power analysis requires
  (continuous age, ~50+ donors per phase, per-cell-type modeling). One cohort = no
  cross-study integration. Late-adulthood findings (glial immune reactivation; convergence
  onto shared programs; cell-type-specific trajectories replicate best) speak directly to
  the mouse story and would power the strongest Chapter 5 agreement check yet.
- **Access (confirmed from the paper's repo and medRxiv data availability):**
  - Synapse accession **syn52396927** (doi:10.7303/syn52396927), PsychAD Study syn60084804.
  - AD Knowledge Portal (adknowledgeportal.synapse.org), controlled-use tier.
    A **data use agreement (DUA) is required**; registration exists to protect participant
    anonymity. Turnaround typically days, not months, per ADKP policy; no institutional SSO
    wall (unlike SenNet). ACTION (Reza): Synapse account + ADKP data-access request for
    syn52396927. ACTION (agent): draft request text.
  - The CELLxGENE explorer link in the repo is a **private/shared collection**, NOT public:
    api.cellxgene.cziscience.com returns 404 for the dataset and the collection is absent
    from the 401 public collections (verified 2026-10-21). No shortcut around the DUA.
- **Analysis plan once granted:** donor-by-subclass pseudobulk (they ship this), map mouse
  program genes 1:1 (Ensembl, existing mapper), score donors, regress program score on
  continuous age per cell class; report slope, CI, and a powered-vs-underpowered verdict per
  class using the M4 grid.

### 2. CELLxGENE public collections with donor ages (open NOW, no auth)

Collections confirmed public on CELLxGENE Discover (401 public collections scanned):

| Collection | Paper | Access |
|---|---|---|
| Human skeletal muscle ageing atlas | Nat Aging 2024 (s43587-024-00613-3) | public h5ad via census |
| Immune aging, multimodal tissue signatures | Nat Immunol 2025 (s41590-025-02241-4), GSE299043 | public |
| Neuron-type effects of human aging + sex | Neuron 2024 (mukamel-lab, GSE247988) | public |
| Human microglia/perivascular macrophages in aging | Nat Genet 2026 (s41588-026-02716-6) | public |
| Human Immune Health Atlas | Nature 2025 (s41586-025-09686-5) | public (reference-grade immune) |

### 3. SenNet — tissue-matched complement (blocked on registration)

- 66 human lung snRNA-seq datasets; standardized Salmon+Scanpy pipelines; donor ages coded
  on source records; mouse data in the same consortium. Login requires an institutional
  identity provider (no Google/Globus/ORCID per the SenNet FAQ). Help-desk email sent
  2026-10-21 (Reza); awaiting reply. Proceeds only if the help desk replies or a
  university SSO surfaces.

### Dropped during scoping

- **Catching et al., Cell Rep 2026** (lifespan DLPFC snRNA, NABEC): accession could not
  be located from available sources; superseded by source 1. Do not cite.
- **CELLxGENE scatter meta-analysis** (2,243 human lung snRNA datasets): ages
  inconsistently exposed; integration burden high; superseded by targeted collections.

## Census probe results (2026-10-21, census stable 2025-11-08)

Human, primary cells only, donors carrying a numeric age in `development_stage`
(HsapDv "N-year-old stage"). Caveat: disease composition was NOT filtered in this
pass; final cohort selection must restrict to control/normal donors via the obs
`disease` field at download time.

**Lung + kidney: 44 datasets, 694 donors with numeric ages.** Notable:
| dataset | donors | ages | note |
|---|---|---|---|
| Integrated human lung atlas (health AND disease) | 131 | 2-86 | needs disease filter |
| Human CellCards Multi-Study CellRef 1.0 (lung) | 21 | 12-76 | reference atlas |
| Adult Human Kidney scRNA v1.5 | 18 | 30-66 | |
| Normal adult kidney tissue | 18 | 43-74 | healthy |
| Living donor kidney | 19 | 33-65 | healthy |
| CZI Immune Aging (all cells) | 12 | 20-75 | BAL/PBMC; more donors in full study |

**Brain: 167 datasets, 2,296 donors with numeric ages.** Notable:
| dataset | donors | ages | note |
|---|---|---|---|
| MSSM_Cohort (DLPFC, PsychAD) | 783 | 19-89 | psychiatric/AD mix; control subset needed |
| HBCC_Cohort (DLPFC, NIMH) | 288 | 2-85 | lifespan incl. development |
| RADC_Cohort (DLPFC) | 100 | 54-89 | aging/AD cohort |
| SEA-AD DLPFC | 38 | 29-89 | aging/AD atlas |
| Extended GBmap | 130 | 21-83 | glioma-adjacent caution |

Implication: a POWERED human arm is achievable with zero access requests, in
DLPFC (hundreds of donors with ages) and more modestly in lung/kidney (dozens of
mostly-healthy donors). This upgrades M9 from "request and wait" to "build now,
request in parallel."

## M9 execution plan (next)

1. Cohort selection (disease-filtered control donors; per-cell-type counts >= floor).
2. Census download of obs+counts slices for selected datasets (control donors only).
3. Ortholog map to the 261-gene program (existing Ensembl 1:1 mapper).
4. Donor-level program score; linear regression on continuous age per cell class;
   M4 power verdict per cohort.
5. Agreement check vs mouse cell-type programs and (once DUA lands) vs the DLPFC
   atlas's own aDEG directions.
6. Ship as Chapter 5 upgrade + report addendum; same provenance/export/gate rules.

## Decision

Primary: DLPFC Lifespan Atlas (DUA request today). Parallel: open CELLxGENE aging
collections scored this week. SenNet proceeds only if the help desk replies.
