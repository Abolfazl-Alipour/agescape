# Mature interneuron GO enrichment: follow-up report

Date: 2026-11-24. Owner-ordered on the 2026-11-24 voice call (about 10:35 PT: "Okay, do it").
Script: `scripts/run_mature_interneuron_go.py`.
Results: `results/followup-mature-interneuron-go.csv` (12,414 rows) and
`results/followup-mature-interneuron-go.meta.json` (run metadata).

Terms used here: GO is the Gene Ontology, a catalog that groups genes by the
biological functions they take part in. BP is Biological Process, one branch
of that catalog (what a gene does in the cell, not where it sits or what
molecule it binds). FDR is the false discovery rate, the share of claimed
hits expected to be false; q is the FDR-adjusted p value. BH is the
Benjamini-Hochberg correction, the method used to turn raw p values into q
values within each gene list. Fisher exact is the one-sided hypergeometric
test used for each term.

## 1. Headline result: no term passes in any list

No GO Biological Process term reaches q below 0.05 in any of the 14 gene
lists tested. That covers the per-type both-FDR lists, the positive-only and
negative-only splits, the 49-gene union, and the four mature MSSM discovery
sets. The null is honest, not a failure to run: the tests ran over about
6,600 to 7,400 terms per list, and the smallest q values are 0.11 to 0.14 in
the best-powered lists and 1.0 in most small lists. Mature only: every input
comes from the mature (donors age 30 and older) interneuron results. No
full-range (young-inclusive) gene list is used anywhere.

## 2. Per-type significant terms (both-FDR lists)

There are no significant terms. The closest term per type is listed so the
null is checkable, not just asserted.

1. VIP, 4 genes, universe 1,920 mapped genes, 6,803 terms tested. Smallest p
   0.0021, smallest q 1.0. The top rows are single-gene terms driven by
   CLSTN3 (for example negative regulation of lipid storage, GO:0010888,
   1 of 1 term gene in the list, fold enrichment 480). One gene cannot
   support a function claim.
2. Parvalbumin, 8 genes, universe 1,886, 6,687 terms tested. Smallest p
   4.7e-05, smallest q 0.315. The top term is positive regulation of muscle
   cell apoptotic process (GO:0010661, genes CAMK2D and PPARG, 2 of 3 term
   genes, fold enrichment 157). It misses the bar by a factor of about 6.
3. Somatostatin, 41 genes, universe 1,931, 6,841 terms tested. Smallest p
   2.0e-05, smallest q 0.140. The top term is regulation of calcium ion
   transmembrane transport via high voltage-gated calcium channel
   (GO:1902514, genes CACNA2D1, CACNB2, CAMK2D, NOS1AP, 4 of 9 term genes,
   fold enrichment 20.9). It misses the bar by a factor of about 3.
4. GABA both-FDR is empty: 0 genes, so no test was run for GABA in this
   direction set. GABA appears only in the MSSM discovery set below.

## 3. Positive-only and negative-only splits

There are no significant terms in any split. The closest term per split:

1. Positive VIP, 3 genes: smallest p 0.0078 (RNA methylation, GO:0001510,
   gene TRMT11 alone), smallest q 1.0.
2. Positive parvalbumin, 6 genes: smallest p 2.5e-05, smallest q 0.169
   (same muscle cell apoptotic process term as the full parvalbumin list,
   genes CAMK2D and PPARG).
3. Positive somatostatin, 35 genes: smallest p 4.3e-04, smallest q 1.0
   (membrane depolarization during cardiac muscle cell action potential,
   GO:0086012, genes CACNA2D1, CACNB2, CAMK2D).
4. Negative VIP, 1 gene (CLSTN3): smallest p 5.2e-04, smallest q 0.506.
   A one-gene list cannot pass a correction over 6,803 terms.
5. Negative parvalbumin, 2 genes (CLSTN2, SYP): smallest p 0.0011,
   smallest q 1.0. The top terms are synaptic vesicle terms driven by SYP
   alone.
6. Negative somatostatin, 6 genes: smallest p 0.0031, smallest q 1.0.
   The top rows are single-gene terms (NOS1AP, ILDR2, MAN1A1).

## 4. MSSM discovery sets (where power lives)

There are no significant terms here either. These four lists are the largest
(36 to 355 genes), so this is the informative null, not a power artifact
of the 4-gene and 8-gene lists.

1. VIP, 137 mapped genes of 139 (unmapped: BRD10, ENSG00000249738),
   universe 1,877, 6,773 terms tested. Smallest p 3.8e-04, smallest q 1.0.
   The top terms are antimicrobial humoral response (GO:0019730, genes APP,
   CXCL14, VIP, 3 of 3 term genes, fold enrichment 13.7) and G
   protein-coupled receptor signaling pathway (GO:0007186, 15 of 90 term
   genes, p 0.0016, fold enrichment 2.3).
2. Parvalbumin, 127 mapped genes of 128 (unmapped: BRD10), universe 1,827,
   6,639 terms tested. Smallest p 0.0013, smallest q 1.0. The top terms are
   retina layer formation (GO:0010842, genes ATP8A2, DSCAM, PTPRM, 3 of 4)
   and axon guidance (GO:0007411, 9 of 47 term genes, p 0.0040, fold
   enrichment 2.8).
3. Somatostatin, 352 mapped genes of 355 (unmapped: ENSG00000228222,
   ENSG00000233928, ENSG00000237773), universe 1,870, 6,753 terms tested.
   Smallest p 2.0e-05, smallest q 0.112. Six terms sit at q 0.112, all
   calcium or ion transport or muscle related:
   1. intracellular calcium ion homeostasis (GO:0006874, 17 of 33, p
      2.0e-05, fold enrichment 2.7),
   2. monoatomic ion transport (GO:0006811, 42 of 128, p 6.3e-05, fold
      enrichment 1.7),
   3. monoatomic ion transmembrane transport (GO:0034220, 40 of 120, p
      6.4e-05, fold enrichment 1.8),
   4. calcium ion homeostasis (GO:0055074, 17 of 36, p 8.4e-05, fold
      enrichment 2.5),
   5. regulation of cardiac muscle cell action potential (GO:0098901, 9 of
      13, p 9.5e-05, fold enrichment 3.7),
   6. regulation of muscle system process (GO:0090257, 19 of 43, p 9.9e-05,
      fold enrichment 2.4).
   This is the closest approach to significance in the whole run, and it
   still misses q 0.05 by more than a factor of 2.
4. GABA, 36 genes, universe 1,835, 6,846 terms tested. Smallest p 0.0014,
   smallest q 1.0. The top term is adenylate cyclase-inhibiting G
   protein-coupled receptor signaling pathway (GO:0007193, genes APP,
   GABBR1, GRM8, 3 of 12).

## 5. Union of both-FDR genes (49 genes)

No significant term. Universe 2,355 mapped genes (2,442 symbols), 7,362
terms tested. Smallest p 1.9e-05, smallest q 0.142. The top term is the same
calcium channel term as the somatostatin list (GO:1902514, genes CACNA2D1,
CACNB2, CAMK2D, NOS1AP, 4 of 9, fold enrichment 21.4). 26 terms reach raw
p below 0.01, and all but one have q 1.0; the nominal leaders are calcium
transport, muscle and cardiac terms driven by the same 4 to 6 genes
(CACNA2D1, CACNB2, CAMK2D, NOS1AP, PPARG, FOXO3, PPP1R12B), plus DNA
damage tolerance (GO:0006301, genes REV3L and WDR33, p 0.0025) and positive
regulation of synapse assembly (GO:0051965, genes CLSTN2, CLSTN3, EFNA5,
p 0.0055).

## 6. Shared-term table across types

At q below 0.05 the shared table is empty: no term is significant in any
type, so no term can be shared. At nominal raw p below 0.01 (uncorrected,
shown only to check whether the types point the same way at all), the
overlap is nearly zero:

1. MSSM discovery, VIP and parvalbumin: 0 shared terms.
2. MSSM discovery, VIP and somatostatin: 1 shared term, calcium ion
   transport (GO:0006816, VIP p 0.0086, somatostatin p 0.0014).
3. MSSM discovery, parvalbumin and somatostatin: 1 shared term, retina
   layer formation (GO:0010842, parvalbumin p 0.0012, somatostatin
   p 0.0012). This is a generic annotation of the same three genes, not a
   shared aging program.
4. MSSM discovery, all three types: 0 shared terms.
5. Both-FDR, all three types: 0 shared terms.
6. Both-FDR, VIP and parvalbumin: 2 shared terms, both lipid terms driven
   by different single genes (CLSTN3 in VIP, PPARG in parvalbumin):
   negative regulation of lipid storage (GO:0010888) and negative
   regulation of lipid localization (GO:1905953).
7. Both-FDR, VIP and somatostatin: 0 shared terms.
8. Both-FDR, parvalbumin and somatostatin: 1 shared term, regulation of
   muscle system process (GO:0090257, parvalbumin p 0.0097, somatostatin
   p 0.0020), driven in both types by CAMK2D, so the overlap is one shared
   gene, not two types converging on one function.

## 7. Universe sizes and sources

Gene lists come from `results/followup-mature-interneuron-genes.csv`
(658 rows: genes FDR-significant at q below 0.05 in at least one of MSSM or
HBCC, mature unmoderated WLS, donors age 30 and older). Counts were
re-verified from the file before running and match the task brief exactly:
both-FDR rows 53 (VIP 4, parvalbumin 8, somatostatin 41), 49 unique genes;
positive-only 40 unique genes (VIP 3, parvalbumin 6, somatostatin 35);
negative-only 9 unique genes (VIP 1, parvalbumin 2, somatostatin 6);
TRMT11 is the only gene significant in all three types in both cohorts;
GABA both-FDR 0; MSSM discovery FDR sets VIP 139, parvalbumin 128,
somatostatin 355, GABA 36.

Universes come from `results/followup-mature-interneuron.csv` (45,357 rows,
tested genes per cohort and subtype, model rev2_voom_wls). Tested-gene
counts (n_genes_tested) match the brief: MSSM VIP 1,932, parvalbumin 1,873,
somatostatin 1,906, GABA 1,883; HBCC VIP 1,916, parvalbumin 1,866,
somatostatin 1,895, GABA 1,848. The 658-row FDR-only file was never used as
a universe.

| List | Genes (mapped) | Universe symbols (mapped to Entrez) | Terms tested | Min q |
| --- | --- | --- | --- | --- |
| both-FDR VIP | 4 (4) | 1,986 (1,920) | 6,803 | 1.0 |
| both-FDR parvalbumin | 8 (8) | 1,942 (1,886) | 6,687 | 0.315 |
| both-FDR somatostatin | 41 (41) | 1,975 (1,931) | 6,841 | 0.140 |
| positive VIP | 3 (3) | 1,986 (1,920) | 6,803 | 1.0 |
| positive parvalbumin | 6 (6) | 1,942 (1,886) | 6,687 | 0.169 |
| positive somatostatin | 35 (35) | 1,975 (1,931) | 6,841 | 1.0 |
| negative VIP | 1 (1) | 1,986 (1,920) | 6,803 | 0.506 |
| negative parvalbumin | 2 (2) | 1,942 (1,886) | 6,687 | 1.0 |
| negative somatostatin | 6 (6) | 1,975 (1,931) | 6,841 | 1.0 |
| union | 49 (49) | 2,442 (2,355) | 7,362 | 0.142 |
| MSSM discovery VIP | 139 (137) | 1,932 (1,877) | 6,773 | 1.0 |
| MSSM discovery parvalbumin | 128 (127) | 1,873 (1,827) | 6,639 | 1.0 |
| MSSM discovery somatostatin | 355 (352) | 1,906 (1,870) | 6,753 | 0.112 |
| MSSM discovery GABA | 36 (36) | 1,883 (1,835) | 6,846 | 1.0 |

The universe for both-FDR, positive, and negative lists is the MSSM plus
HBCC tested-gene union for that subtype. The union list uses the union of
the three subtype universes. The MSSM discovery lists use the MSSM tested
genes for that subtype only. Mapped counts are lower than symbol counts
because some tested symbols have no Entrez ID in the symbol map (mostly
Ensembl IDs such as ENSG00000228222 and the symbol BRD10); those genes
cannot carry GO annotations and are excluded from both the list and the
universe. The counts are recorded per list in the metadata JSON.

## 8. Tool and GO versions

1. First attempt: clusterProfiler via BiocManager on R 4.3.3. org.Hs.eg.db
   3.18.0 installed, but clusterProfiler failed: fgsea C++ compilation
   failed, which broke DOSE, enrichplot, and clusterProfiler. No
   clusterProfiler result was produced.
2. Fallback used, as ordered: Python GOATOOLS 1.6.9 with scipy 1.18.1,
   Python 3.12. GOATOOLS supplied the GO DAG and the gene2go association
   reader. The Fisher test and the BH correction were computed directly
   (scipy fisher_exact, hand-coded BH) so the universe is exactly the
   tested genes; GOATOOLS' own study class restricts the population to
   annotated genes, which would have changed the universe.
3. GO ontology: go-basic.obo release 2026-07-26 (41,378 terms), downloaded
   2026-11-24 from current.geneontology.org.
4. GO annotations: NCBI gene2go, downloaded 2026-11-24 (file dated Nov 24
   2026), human taxid 9606 only, 446,672 annotations, 20,362 genes,
   17,987 GO terms. Direct annotations were propagated to Biological
   Process ancestors through the DAG; the BP root GO:0008150 was excluded.
   18,130 human genes carry at least one BP annotation; 13,551 BP terms
   appear in the annotation.
5. Symbol to Entrez mapping: org.Hs.eg.db 3.18.0 (Entrez Gene 2023-Sep-11,
   GO 2023-07-27), official symbol first, alias as fallback. Note the
   vintage gap: the symbol map is from 2023 while the GO ontology and
   annotations are from 2026. Unmapped list genes: BRD10 (VIP and
   parvalbumin discovery), ENSG00000249738 (VIP discovery),
   ENSG00000228222, ENSG00000233928, ENSG00000237773 (somatostatin
   discovery). No both-FDR, positive, negative, or union gene is unmapped.
6. Method: one-sided Fisher exact (hypergeometric, greater) per term per
   list; BH correction within each list over all BP terms with at least one
   universe gene annotated (about 6,600 to 7,400 tests per list); rows
   written for every term with overlap of at least 1; q below 0.05 is
   significant. The old full-range file
   `results/reviewer-phase1-interneuron-go-enrichment.csv` was used only as
   a format reference, never as an input.
7. Two script bugs were found during the run and fixed before the final
   run: GOATOOLS namespace values are `biological_process`, not `BP`
   (first run tested 0 terms), and `get_all_parents()` returns GO ID
   strings, not term objects (second run crashed). The final run is the
   third run; its log is `/tmp/go_run3.log` on Beast.

## 9. Verification

1. The top somatostatin discovery term (GO:0006874, N 1,870, n 352, K 33,
   k 17) was recomputed by hand: Fisher p 1.99207e-05, hypergeometric
   survival function 1.99207e-05, expected 6.21, all equal to the CSV row.
2. The BH code was checked against statsmodels `multipletests`
   (fdr_bh) on a 6,753-value random vector: maximum absolute difference
   2.2e-16.
3. Annotations were spot-checked in the gene2go file: TRMT11 (Entrez
   60487) carries a direct GO:0001510 RNA methylation annotation, and
   GO:0006874 has 146 direct human annotations, consistent with 33 in the
   somatostatin universe.
4. The gene lists and universes were re-counted from the source CSVs
   before running; every count matches the brief (Section 7).

## 10. Caveats

1. BH over about 6,700 terms is strict. The somatostatin discovery calcium
   signal (q 0.112) would pass a looser bar but does not pass q 0.05, and
   it rests on 17 overlapping genes out of a 352-gene list.
2. The small lists (1 to 8 genes) have almost no power by design. Their
   nulls say nothing beyond "too few genes to test". The informative null
   is the MSSM discovery sets, and those are null too.
3. GO annotations are heavily IEA (electronic, unreviewed) in gene2go, and
   the ontology changed between the 2023 symbol map and the 2026 GO files.
   A rerun with clusterProfiler and org.Hs.eg.db end to end would be a
   useful cross-check if the fgsea build is ever fixed, but it should not
   be expected to overturn a null this broad.
4. Nominal p values below 0.01 in Section 6 are uncorrected and are shown
   only to test sharing. With about 6,700 tests per list, about 67 terms
   reach p below 0.01 by chance alone per list, so one or two shared
   nominal terms across types is consistent with chance.

## 11. Verdict: do the three types share functions?

No. No GO Biological Process term is significant in any interneuron type,
so there is no shared significant function to report. The weaker nominal
view does not rescue sharing either: across the three MSSM discovery sets,
no term reaches raw p below 0.01 in all three types, and the pairwise
overlaps are one term each (calcium ion transport for VIP and somatostatin;
retina layer formation for parvalbumin and somatostatin), which is what
chance alone predicts over about 6,700 tests. The one recurring theme is
calcium ion transport and homeostasis, strongest in somatostatin, where it
is the closest near-miss in the whole run (six terms at q 0.112, driven by
calcium channel and pump genes such as CACNA2D1, CACNB2, CAMK2D, ATP2A2,
and RYR2), and it also tops the both-FDR somatostatin list and the union
(q 0.14). That is a somatostatin-leaning hint, not a shared program: VIP
and parvalbumin show no calcium signal of their own, the calcium terms do
not survive correction anywhere, and the remaining nominal leaders point in
different directions per type (G protein-coupled receptor signaling in VIP,
axon guidance in parvalbumin). The defensible sentence for the paper is
that the mature age-related genes in the three interneuron types show no
shared GO Biological Process enrichment, with a nominal, uncorrected
calcium-transport signal confined to somatostatin.
