# Follow-up Blocker Report: Cultured Adipose Stem Cells (`GSE276163` / `GSE276164`)

**Status**: BLOCKED / NOT ANALYZED (`notes/followup-dataset-checklist.md` #33: `SKIP`; cultured passaged bulk RNA-seq + severe sample/sex swap and contamination)
**Accessions**: [`GSE276163`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276163) (bulk mRNA-seq of passaged ASCs, `GSM8492245` to `GSM8492259`) and [`GSE276164`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276164) (small RNA-seq of pooled sEVs, `GSM8492260` to `GSM8492264`)
**Primary publication**: *Cytotechnology* (2025), "Aging and cell expansion enhance microRNA diversity in small extracellular vesicles produced from human adipose-derived stem cells", DOI `10.1007/s10616-024-00675-6`, `PMID: 39665045`, `PMCID: PMC11631832`
**Audit table**: [`followup-organ-11-adipose-gse276163-audit.csv`](../results/followup-organ-11-adipose-gse276163-audit.csv)

## 1. Primary-Source Metadata Verification

- **Donors and ages (`PMC11631832` Table 1)**: `6` human donors total (`3` young: `Donor 1` Male `26`, `Donor 2` Male `25`, `Donor 3` Female `26`; `3` old: `Donor 4` Male `79`, `Donor 5` Female `80`, `Donor 6` Male `73`).
- **Experimental unit and culture arrest**: Primary subcutaneous adipose-derived stem cells (ASCs) were expanded in vitro to passages `P5`, `P10`, and `P15`. Old donor ASCs ceased proliferating around passages `13-14`, so only `P5` and `P10` were sequenced for the `3` old donors (`6` libraries), whereas `P5`, `P10`, and `P15` were sequenced for the `3` young donors (`9` libraries), yielding `15` repeated-measure bulk RNA-seq libraries (`GSE276163_ASCs_mRNA_readcount.txt.gz`).
- **`GSE276164` companion matrix**: Contains `5` pooled group-by-passage sEV miRNA-seq conditions (`Young P5`, `Young P10`, `Young P15`, `Old P5`, `Old P10`; `GSM8492260` to `GSM8492264`), not donor-replicated mRNA-seq (`n=1` pooled sample per age-passage condition).

## 2. Documented Blockers (Why No Modeling Was Run)

1. **Uncorrected sex-label mismatches and missing female library in `GSE276163`**: Based on Table 1 (`2F / 4M` across `6` donors) and GEO sample titles (`Sample 1..15`), there should be `5` female libraries (`Donor 3` at `P5/P10/P15` = `Sample 7, 8, 9` plus `Donor 5` at `P5/P10` = `Sample 12, 13`) and `10` male libraries. Direct inspection of sex-chromosome counts in [`audit.csv`](../results/followup-organ-11-adipose-gse276163-audit.csv) reveals only **`4` `XIST+` female libraries** and **`11` Y-expressing male libraries**:
   - `Sample 3` (`GSM8492247`, labeled `Donor 1-passage 15`, Male `26`): expresses `XIST = 2,345.48` with `0.00` counts across all 6 Y genes (`UTY, DDX3Y, KDM5D, RPS4Y1, EIF1AY, ZFY`).
   - `Sample 8` (`GSM8492252`, labeled `Donor 3-passage 10`, Female `26`): expresses `XIST = 0.00` and Y-gene sum `3,066.63` (`UTY=155, DDX3Y=625, KDM5D=243, RPS4Y1=1,830.63, EIF1AY=87, ZFY=126`).
   - `Sample 9` (`GSM8492253`, labeled `Donor 3-passage 15`, Female `26`): expresses `XIST = 2.00` and Y-gene sum `3,548.18` (`UTY=140, DDX3Y=549, KDM5D=234, RPS4Y1=2,387.18, EIF1AY=91, ZFY=147`).
2. **Sample 13 male Y-chromosome contamination**: `Sample 13` (`GSM8492257`, labeled `Donor 5-passage 10`, Female `80`) expresses `XIST = 2,799.00` alongside trace Y-chromosome reads totaling `115.60` (`RPS4Y1=78.60, DDX3Y=19.00, EIF1AY=8.00, KDM5D=5.00, UTY=4.00, ZFY=1.00`), whereas `Sample 12` (`Donor 5-passage 5`) has `0.00` across all 6 Y genes.
3. **In-vitro passaging / missing old P15 / fractional RSEM counts**: Even ignoring the sample-label corruption, `GSE276163` profiles cultured expanded ASCs (with senescence onset driven heavily by passage number and complete loss of the `P15` arm in old donors), uses RSEM expected fractional counts (`154,576` non-integer entries out of `1,025,550` cells across `68,370` genes), and spans only `n=6` donors (`3` young, `3` old). No single-cell or bulk aging claim can be supported from `GSE276163` or `GSE276164`.
