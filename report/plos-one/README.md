# PLOS ONE manuscript

Source of truth: `final-plos-one.tex`.

Compiled review copy: `final-plos-one.pdf`, built with Tectonic 0.15.0.

The current title is "No uniform transcriptional response detected in mouse and human aging". The abstract is 287 words, and the PDF is 33 pages with four figures and 68 cited references.

The manuscript retains the scoped conclusion that no uniform response was detected under its stringent operational construct, without settling upstream causes or excluding shared mechanisms with cell-type-specific consequences. The three additional mature-only findings are labeled secondary and non-prespecified.

The optimization uses the 2,000-draw matched-null analysis as primary, retains the 624-test composition family with core and exploratory counts and sensitivity disclosures, and narrows unsupported trajectory, mechanism, disease-progression, and theory claims. `workflow-manifest.tsv` records the ordered analysis workflow; `python3 scripts/run_plos_workflow.py` from the repository root verifies its commands and committed outputs.

Affiliation is Kavox Labs. Funding, competing interests, and ethics statements are author-confirmed in the source. No journal submission has been made.

To rebuild from this directory:

```bash
tectonic final-plos-one.tex
```

The older `report/manuscript-aging-mechanisms.md` is not part of this manuscript build.
