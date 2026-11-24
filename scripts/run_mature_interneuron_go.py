#!/usr/bin/env python3
"""Mature interneuron GO enrichment (Biological Process only).

Owner-ordered 2026-11-24 (voice call, about 10:35 PT: "Okay, do it").
MATURE ONLY: all inputs are the mature (donors age 30 and older) interneuron
results. No full-range (young-inclusive) gene list is used anywhere.

Inputs (repo-relative):
  results/followup-mature-interneuron-genes.csv  658 rows, genes FDR-significant
      (q < 0.05) in at least one of MSSM / HBCC, mature unmoderated WLS.
  results/followup-mature-interneuron.csv        45,357 rows, tested genes per
      cohort x interneuron subtype (rev2_voom_wls).

Gene lists tested (short_type in files: VIP / pvalb / sst / GABA):
  direction_set  list_name      genes
  both_FDR       VIP            both-FDR (q<0.05 in MSSM and HBCC), 4
  both_FDR       parvalbumin    both-FDR, 8
  both_FDR       somatostatin   both-FDR, 41
  positive       VIP            positive-only both-FDR, 3
  positive       parvalbumin    positive-only both-FDR, 6
  positive       somatostatin   positive-only both-FDR, 35
  negative       VIP            negative-only both-FDR, 1
  negative       parvalbumin    negative-only both-FDR, 2
  negative       somatostatin   negative-only both-FDR, 6
  union          union          union of both-FDR genes, 49 unique
  MSSM_discovery VIP            mature MSSM FDR (q_mssm<0.05), 139
  MSSM_discovery parvalbumin    mature MSSM FDR, 128
  MSSM_discovery somatostatin   mature MSSM FDR, 355
  MSSM_discovery GABA           mature MSSM FDR, 36
  (GABA both-FDR / positive / negative are empty: 0 genes, no test run.)

Universe per list: the tested genes for that subtype. For the both-FDR,
positive and negative lists the universe is the MSSM + HBCC tested-gene union
for that subtype; for the union list it is the union of the three subtype
universes; for MSSM_discovery it is the MSSM tested genes for that subtype.

Method:
  GO Biological Process only. Annotations: NCBI human gene2go (taxid 9606),
  current download, propagated to ancestors through the GO DAG (GOATOOLS
  GODag, current go-basic.obo). The GO root term GO:0008150 is excluded.
  One-sided Fisher exact (hypergeometric, greater) per term, computed with
  scipy on universe N = tested genes (unannotated universe genes count as
  non-members). Benjamini-Hochberg correction within each list over all BP
  terms with at least one universe gene annotated (n_terms_tested). Rows are
  written for every term with overlap >= 1; q < 0.05 is significant.
  Small lists (1-8 genes) have almost no power: they are run anyway and the
  null is reported honestly.

Tool note: clusterProfiler + org.Hs.eg.db were attempted first via
BiocManager. org.Hs.eg.db 3.18.0 installed, but clusterProfiler failed
(fgsea C++ compilation failed, which broke DOSE / enrichplot /
clusterProfiler). Fallback per task order: Python GOATOOLS with a current GO
OBO file and human gene2go annotations. GOATOOLS is used for the GO DAG and
the gene2go association reader; the Fisher test and BH are computed directly
so the universe is exactly the tested genes (GOATOOLS' own study class
restricts the population to annotated genes).

Outputs:
  results/followup-mature-interneuron-go.csv   (phase-1 column format plus
      direction_set; genes are ';'-separated symbols)
  <meta>.json  run metadata: versions, GO release, counts, unmapped genes.

Usage (from repo root):
  python3 scripts/run_mature_interneuron_go.py \
    --obo /tmp/go_data/go-basic.obo \
    --gene2go-human /tmp/go_data/gene2go_human.tsv \
    --symbol-map /tmp/go_data/symbol_map_symbol2entrez.tsv \
    --alias-map /tmp/go_data/symbol_map_alias2entrez.tsv \
    --out results/followup-mature-interneuron-go.csv \
    --meta results/followup-mature-interneuron-go.meta.json
"""
import argparse
import csv
import json
import sys
from collections import defaultdict

GO_ROOT_BP = "GO:0008150"
TYPE_LABEL = {"VIP": "VIP", "pvalb": "parvalbumin", "sst": "somatostatin",
              "GABA": "GABA"}
TYPES = ["VIP", "pvalb", "sst", "GABA"]


def load_symbol_maps(symbol_path, alias_path):
    """Return (symbol -> entrez, alias -> entrez)."""
    sym = {}
    with open(symbol_path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            sym.setdefault(row["SYMBOL"], row["ENTREZID"])
    ali = {}
    with open(alias_path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            ali.setdefault(row["ALIAS"], row["ENTREZID"])
    return sym, ali


def map_symbols(symbols, sym, ali):
    """Map gene symbols to Entrez IDs. Returns (entrez->symbol, unmapped)."""
    entrez_to_symbol = {}
    unmapped = []
    for s in sorted(set(symbols)):
        eid = sym.get(s)
        if eid is None:
            eid = ali.get(s)
        if eid is None:
            unmapped.append(s)
        else:
            entrez_to_symbol.setdefault(int(eid), s)
    return entrez_to_symbol, unmapped


def benjamini_hochberg(pvalues):
    """Standard BH step-up adjusted p-values, same order as input."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvalues[i] * m / (rank + 1))
        adj[i] = min(prev, 1.0)
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="results/followup-mature-interneuron-genes.csv")
    ap.add_argument("--universe", default="results/followup-mature-interneuron.csv")
    ap.add_argument("--obo", required=True)
    ap.add_argument("--gene2go-human", required=True)
    ap.add_argument("--symbol-map", required=True)
    ap.add_argument("--alias-map", required=True)
    ap.add_argument("--out", default="results/followup-mature-interneuron-go.csv")
    ap.add_argument("--meta", default="results/followup-mature-interneuron-go.meta.json")
    args = ap.parse_args()

    import goatools
    from goatools.obo_parser import GODag
    from goatools.associations import read_ncbi_gene2go
    from scipy.stats import fisher_exact

    # ---- gene lists and universes from the mature files ----
    with open(args.genes) as fh:
        grows = list(csv.DictReader(fh))
    with open(args.universe) as fh:
        urows = list(csv.DictReader(fh))

    both = [r for r in grows
            if float(r["q_bh_mssm"]) < 0.05 and float(r["q_bh_hbcc"]) < 0.05]
    pos = [r for r in both
           if float(r["slope_mssm"]) > 0 and float(r["slope_hbcc"]) > 0]
    neg = [r for r in both
           if float(r["slope_mssm"]) < 0 and float(r["slope_hbcc"]) < 0]
    mssm_fdr = [r for r in grows if float(r["q_bh_mssm"]) < 0.05]

    # universe symbols per subtype: tested genes, MSSM+HBCC union, and MSSM only
    univ_union = {t: set() for t in TYPES}
    univ_mssm = {t: set() for t in TYPES}
    for r in urows:
        t = r["short_type"]
        univ_union[t].add(r["gene"])
        if r["cohort"] == "MSSM_DLPFC":
            univ_mssm[t].add(r["gene"])

    def genes_of(rows, t):
        return sorted(set(r["gene"] for r in rows if r["short_type"] == t))

    lists = []  # (direction_set, list_name, gene symbols, universe symbols)
    for t in ["VIP", "pvalb", "sst"]:
        lists.append(("both_FDR", TYPE_LABEL[t], genes_of(both, t), univ_union[t]))
        lists.append(("positive", TYPE_LABEL[t], genes_of(pos, t), univ_union[t]))
        lists.append(("negative", TYPE_LABEL[t], genes_of(neg, t), univ_union[t]))
    lists.append(("union", "union",
                  sorted(set(r["gene"] for r in both)),
                  set().union(*[univ_union[t] for t in ["VIP", "pvalb", "sst"]])))
    for t in TYPES:
        lists.append(("MSSM_discovery", TYPE_LABEL[t],
                      genes_of(mssm_fdr, t), univ_mssm[t]))

    # ---- GO annotations: GOATOOLS DAG + human gene2go ----
    godag = GODag(args.obo)
    with open(args.obo) as fh:
        obo_version = None
        for line in fh:
            if line.startswith("data-version:"):
                obo_version = line.split(":", 1)[1].strip()
                break
    assoc = read_ncbi_gene2go(args.gene2go_human, taxids=[9606])
    # assoc: entrez int -> set of GO ids (direct annotations, all namespaces)
    # Propagate to BP ancestors only; drop the BP root.
    gene_bp = {}
    for eid, gos in assoc.items():
        terms = set()
        for go in gos:
            term = godag.get(go)
            # GOATOOLS namespace values are 'biological_process',
            # 'molecular_function', 'cellular_component' (not BP/MF/CC).
            if term is None or term.namespace != "biological_process":
                continue
            terms.add(go)
            # get_all_parents() returns GO ID strings, not term objects.
            for parent_id in term.get_all_parents():
                parent = godag.get(parent_id)
                if (parent is not None
                        and parent.namespace == "biological_process"
                        and parent.id != GO_ROOT_BP):
                    terms.add(parent.id)
        terms.discard(GO_ROOT_BP)
        if terms:
            gene_bp[eid] = terms
    bp_terms = sorted(set().union(*gene_bp.values())) if gene_bp else []
    term_name = {go: godag[go].name for go in bp_terms}

    # ---- symbol mapping ----
    sym, ali = load_symbol_maps(args.symbol_map, args.alias_map)

    # ---- enrichment per list ----
    out_rows = []
    meta_lists = {}
    for direction, label, genes, univ_syms in lists:
        list_map, list_unmapped = map_symbols(genes, sym, ali)
        univ_map, univ_unmapped = map_symbols(univ_syms, sym, ali)
        list_eids = set(list_map)
        univ_eids = set(univ_map)
        n_list = len(list_eids)
        n_univ = len(univ_eids)
        # term -> universe genes annotated (propagated BP)
        term_univ_genes = defaultdict(set)
        for eid in univ_eids:
            for go in gene_bp.get(eid, ()):
                term_univ_genes[go].add(eid)
        tested = sorted(term_univ_genes)  # BP terms with >=1 universe gene
        pvals, stats = [], []
        for go in tested:
            k_genes = term_univ_genes[go] & list_eids
            k = len(k_genes)
            K = len(term_univ_genes[go])
            table = [[k, n_list - k],
                     [K - k, n_univ - K - n_list + k]]
            p = fisher_exact(table, alternative="greater")[1] if n_list else 1.0
            expected = n_list * K / n_univ if n_univ else 0.0
            fold = (k / n_list) / (K / n_univ) if n_list and K and n_univ else 0.0
            pvals.append(p)
            stats.append((go, k, K, expected, fold, k_genes))
        qs = benjamini_hochberg(pvals) if pvals else []
        n_sig = 0
        for (go, k, K, expected, fold, k_genes), p, q in zip(stats, pvals, qs):
            if k == 0:
                continue
            if q < 0.05:
                n_sig += 1
            out_rows.append({
                "analysis": "mature_interneuron_go",
                "list_name": label,
                "direction_set": direction,
                "go_id": go,
                "term": "%s (%s)" % (term_name[go], go),
                "n_list": n_list,
                "n_universe": n_univ,
                "n_term": K,
                "overlap": k,
                "expected": expected,
                "fold_enrichment": fold,
                "pvalue": p,
                "genes": ";".join(sorted(list_map[e] for e in k_genes)),
                "q_bh": q,
            })
        meta_lists["%s:%s" % (direction, label)] = {
            "n_list_genes": len(genes),
            "n_list_mapped_entrez": n_list,
            "n_universe_genes": len(univ_syms),
            "n_universe_mapped_entrez": n_univ,
            "n_universe_with_bp_annotation": len(
                [e for e in univ_eids if e in gene_bp]),
            "n_terms_tested": len(tested),
            "n_significant_q_lt_0.05": n_sig,
            "unmapped_list_genes": list_unmapped,
            "unmapped_universe_genes": sorted(univ_unmapped),
        }
        print("%s %s: n_list=%d n_univ=%d terms_tested=%d sig_q05=%d "
              "unmapped_list=%s" % (direction, label, n_list, n_univ,
                                    len(tested), n_sig, list_unmapped),
              flush=True)

    out_rows.sort(key=lambda r: (r["direction_set"], r["list_name"], r["q_bh"],
                                 r["pvalue"]))
    fields = ["analysis", "list_name", "direction_set", "go_id", "term",
              "n_list", "n_universe", "n_term", "overlap", "expected",
              "fold_enrichment", "pvalue", "genes", "q_bh"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    meta = {
        "tool": "Python GOATOOLS (GODag + read_ncbi_gene2go) with scipy "
                "one-sided Fisher exact and hand-coded Benjamini-Hochberg",
        "goatools_version": goatools.__version__,
        "python": sys.version.split()[0],
        "go_obo_file": args.obo,
        "go_obo_release": obo_version,
        "gene2go_file": args.gene2go_human,
        "gene2go_source": "NCBI gene2go, taxid 9606 (human), direct "
                          "annotations, NOT-qualified rows excluded by "
                          "GOATOOLS reader",
        "symbol_map_source": "org.Hs.eg.db 3.18.0 (Entrez Gene 2023-Sep-11, "
                             "GO 2023-07-27); official SYMBOL first, ALIAS "
                             "as fallback",
        "ontology": "Gene Ontology Biological Process only; root GO:0008150 "
                    "excluded; annotations propagated to BP ancestors",
        "test": "one-sided Fisher exact (hypergeometric, greater); BH within "
                "each list over all BP terms with >=1 universe gene; rows "
                "written for overlap >= 1; q < 0.05 significant",
        "inputs": {
            "genes_file": args.genes,
            "genes_rows": len(grows),
            "universe_file": args.universe,
            "universe_rows": len(urows),
            "both_fdr_rows": len(both),
            "both_fdr_unique_genes": len(set(r["gene"] for r in both)),
        },
        "n_bp_terms_in_annotation": len(bp_terms),
        "n_human_genes_with_bp_annotation": len(gene_bp),
        "lists": meta_lists,
    }
    with open(args.meta, "w") as fh:
        json.dump(meta, fh, indent=2)
    print("wrote %s (%d rows) and %s" % (args.out, len(out_rows), args.meta))


if __name__ == "__main__":
    main()
