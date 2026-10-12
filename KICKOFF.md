# KICKOFF — Conserved Aging Explorer, Milestone M0

You are the executing coding agent for the Conserved Aging Explorer project, working in
`/home/reza/projects/agescape` (this is your repo root and working directory).

## Read these files, in this order, before writing any code

1. `/home/reza/projects/agescape/HANDOFF.md` — your operating manual. Follow it exactly.
2. `/home/reza/projects/agescape/PLAN.md` — the pre-registered decisions. When HANDOFF and PLAN
   conflict, PLAN wins and you note the conflict in your gate report.
3. `/home/reza/projects/agescape/notes/reviews.md` — review log; context for why decisions exist.

## Non-negotiables (from PLAN.md §1 and HANDOFF.md §0)

- Never fabricate a number, citation, checksum, or result. Every value traces to an executed cell.
- Never put credentials/tokens anywhere in the repo. SenNet public API is token-free.
- Never make anything public. No remote pushes Reza hasn't approved.
- No ad-hoc package installs after the environment is locked (M0 ends with a lock).
- Stop and report if: stuck >30 min on one error; a step demands credentials; a DoD item is
  unsatisfiable; you want to change a pre-registered choice you can't robustly justify.

## Your task now: Milestone M0 only — Scaffold & verification

Execute HANDOFF.md §3 (Milestone M0) step by step. In summary:

1. Create the repo tree per PLAN.md §3, including `scripts/` and `notes/gates/`. Write the
   `.gitignore` exactly as given in HANDOFF.md §3 (note: `figures/` is NOT ignored).
2. Copy local data into `data/raw/` (do not modify the originals):
   - FROM `/home/reza/projects/aging_hypothalamus/tabula-muris-senis-droplet-official-raw-obj.h5ad`
     TO `data/raw/tms_droplet_raw.h5ad`
     (SHA-256 must be `86e45e88f9c475253792f1bbcdf5851a5be0ef35da72fced607ed261aebc099d`)
   - FROM `/home/reza/projects/aging_hypothalamus/data/tms_brain.h5ad`
     TO `data/raw/tms_facs_brain_nonmyeloid.h5ad`
     (SHA-256 must be `31ee7ffe9fd11fd128925c2dae862399f5c54f27aef9c026dbeaa077ea143339`)
3. Delete the two dead 0-byte files at
   `/home/reza/projects/aging_hypothalamus/data/23872787` and `.../23872787.1` and note the
   deletion in your gate report.
4. Write `environment.yml` per PLAN.md §4.1 (Python 3.11; channels conda-forge + bioconda,
   strict priority; packages listed there — note: doublet detection is `scrublet`). Create env
   `agescape` with conda or mamba at `/home/reza/miniconda3`. Register the kernelspec:
   `conda run -n agescape python -m ipykernel install --user --name agescape`
5. Smoke-test imports per HANDOFF.md §3, then write the lock: `conda list --explicit > environment.lock`
6. Write `scripts/verify_data.py`, `scripts/fetch_data.py`, and `data/manifest/*.json` per
   HANDOFF.md §3 (manifest schema included there).
7. Write `src/__init__.py`, `src/set_seeds.py`, `src/io.py`, `src/plotting.py`.
8. Write `notes/age-bins.md` per PLAN.md §4.2 — this is binding for all later DE work.
9. `git init` if needed; make the initial commit. No remote.

## Definition of done (M0) — every item needs evidence

- `python scripts/verify_data.py` exits 0 on both data files.
- `conda run -n agescape python -c "import scanpy, anndata, pydeseq2, src; print(scanpy.__version__)"` passes.
- `jupyter kernelspec list` shows `agescape`.
- `environment.lock` committed; `.gitignore` correct; tree matches PLAN.md §3.
- `notes/age-bins.md` exists and matches PLAN.md §4.2.
- Dead files deleted.

## When M0 is done — HALT

Do NOT start Milestone M1. Write your gate report to `notes/gates/M0-report.md` using the template
in HANDOFF.md §13, then stop and give Reza a short summary with the checklist evidence. Resume only
when Reza replies with `M0 approved`. Conditional feedback means not approved: implement the
conditions, re-check, re-report.

Machine context: 49 GB RAM, 12 cores, ~247 GB free disk, conda at `/home/reza/miniconda3`.
