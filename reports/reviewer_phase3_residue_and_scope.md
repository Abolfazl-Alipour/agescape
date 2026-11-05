# Submission residue and scope recommendations

## Residue removal list

The external draft text was not available locally, so exact locations below refer to the current repository manuscript at the Phase 3 starting state. The author should repeat this check against the actual submission document.

| Item | Location | Disposition |
|---|---|---|
| Target repository in the title block | `report/manuscript-aging-mechanisms.md:5` | Remove from the title block. Keep the repository in Data and Code Availability. |
| Interactive explorer URL in correspondence | `report/manuscript-aging-mechanisms.md:4` | Remove from correspondence unless the journal explicitly supports a project link. |
| Interactive explorer link in availability | `report/manuscript-aging-mechanisms.md:227` | Remove if the explorer is no longer part of the scientific deliverable. Keep the GitHub repository and data-access statement. |
| Standalone Source lines | No standalone Source line found in the current repository manuscript | Do not carry over any such lines from the external draft. Verify the external document directly. |
| Commit hashes | No commit hash found in the current repository manuscript | Keep commit hashes in reproducibility reports, not in manuscript prose. |
| Test counts and CI status | No test-count or CI-status statement found in the current repository manuscript | Keep these in verification reports, not in the submission narrative. |
| Author-decisions appendix | No author-decisions appendix found in the current repository manuscript | Do not include internal decision logs in the submission. |
| Draft banner or draft-status label | No draft banner found in the current repository manuscript | Remove any such banner from the external document. |
| Pipeline labels that are also biological score names | P0 and P1 occur throughout the Methods and Results | Do not remove these automatically. They define signature tiers and are scientific labels, not repository residue. |
| Unresolved Zhang citation | `report/manuscript-aging-mechanisms.md:220` | This is a reference correction, not residue. The DOI and title could not be verified in the Phase 3 metadata audit and must be checked before submission. |

## Scope recommendation

Main claims should be restricted to mouse and human analyses. These arms have donor-level inference, direct composition and expression measurements, and the most detailed technical sensitivity analyses. PsychAD must not be described as independent replication of public MSSM or HBCC results because exact identifier-confirmed donor overlap exists.

Fly, worm, and zebrafish belong in an exploratory supplement. Fly and worm use pools, age windows and assays differ, and worm effect size depends on the reference. Zebrafish fine-grained analyses use six fish, and reproductive-state changes are entangled with age. NMR belongs in descriptive supplement only: five animals, skin-only data, and coarse compartments cannot establish stability or a between-species difference. Roux worm results remain descriptive because libraries are technical replicates.

The proliferation moderator analysis belongs in supplement as a post hoc exploratory result. The methylation analysis should remain a benchmark reproduction, not validation. None of these qualifications prevents mentioning the broader results, but they should not drive the title or abstract's central claim.

## Title recommendation

Recommended title:

**Cell-type-specific compositional and transcriptional changes in mouse and human aging**

This matches the supported claims without asserting conservation, universality, fibrosis, or transcriptional-information loss. If the exploratory species are included prominently in the revised report, a broader alternative is:

**Cell-type-specific aging across single-cell atlases: shared gene effects and context-dependent compositional change**

The first title is preferable for a focused main report. The author makes the final scope and title decision.
