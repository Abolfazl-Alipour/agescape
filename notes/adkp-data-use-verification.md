# AD Knowledge Portal data-use verification

Checked on 2026-10-31.

## Sources inspected

- AD Knowledge Portal, Data Access: https://adknowledgeportal.synapse.org/DataAccess/
- DLPFC Lifespan Atlas: Synapse `syn52396927`, DOI `10.7303/syn52396927`.
- Pseudobulk input: Synapse `syn62064518`.

The Data Access page was rendered in a browser and its current text was read. It states that individual-level human data in raw, processed, or derived form must not be redistributed outside the AD Knowledge Portal. It also provides requirements for acknowledging data use and for data-availability statements in research outputs.

## Terms applied to Agescape

The M10 outputs are cohort-level regression coefficients, standard errors, confidence intervals, p values, q values, cell counts, and cohort summaries. They do not contain donor identifiers, donor-level expression values, donor-level scores, or per-donor plotting data. These aggregate statistics are the research output permitted by the portal's publication and attribution framework. No controlled individual-level file or derivative donor-level table may be committed or exported.

The portal requires this acknowledgment in publicly shared research outputs:

> The results published here are in whole or in part based on data obtained from the AD Knowledge Portal (NIH-NIA, U24AG061340).

The portal also requires the study-specific acknowledgment supplied on the relevant study page and a data-availability statement naming the AD Knowledge Portal and the study or dataset DOI. The manuscript's data-availability section must retain the acknowledgment, identify `syn52396927` and `syn62064518`, and state that only aggregate results are shared.

## Repository disposition

The controlled RDS and supplementary donor metadata remain under the gitignored `data/` directory. The donor-overlap audit writes only counts by cohort and match category. It does not write matched identifiers or donor-level rows. No downgrade of the aggregate M10 results was required by the inspected terms.
