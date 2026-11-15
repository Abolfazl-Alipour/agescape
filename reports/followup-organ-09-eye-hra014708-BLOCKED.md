# Dataset 09 Access Audit: Ophthalmology Science 2026 (HRA014708 / Figshare 28757048)

## Status
Partial block only. GSA-Human `HRA014708` raw reads are controlled access (not applied for, not bypassed), and the author-clustered `/RPE.h5ad` is not deposited on Figshare. The donor-level analysis was run from the 18 open Figshare Cell Ranger count matrices; see [`reports/followup-organ-09-eye-hra014708.md`](followup-organ-09-eye-hra014708.md).

## Access Verification Details
1. CNCB-NGDC GSA-Human `HRA014708` (`https://ngdc.cncb.ac.cn/gsa-human/browse/HRA014708`): HTTP request returned an anti-bot SecurityCheck interstitial; access is governed by GSA-Human data access approval. No bypass and no application.
2. Figshare `doi:10.6084/m9.figshare.28757048` (`https://api.figshare.com/v2/articles/28757048`): 18 count matrix triplets (293,480 raw barcodes, 285,878 QC cells, integer UMIs), `README.docx`, `Script.py`. Expression sex matched Table 1 on 17 of 18 samples (mismatch: donor D2 (sample(s) 2, 55 yr): Table 1 sex M, expression sex F (QC-cell pseudobulk XIST = 19198 UMIs, TSIX = 915 UMIs, RPS4Y1 + DDX3Y + UTY + KDM5D + EIF1AY = 0 UMIs)).
3. The Figshare description's cohort counts (13 donors, 8 male, 5 female; 7 DR, 3 DM, 8 non-DM samples) differ from Table 1 (12 donors; 8 DR, 3 DM, 7 non-DM samples). Not resolvable from deposited files; Table 1 mapping used.
