#!/usr/bin/env python3

from pathlib import Path
from difflib import SequenceMatcher

import pandas as pd
import requests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"

ENTRIES = [
    ("Atlas", "A single-cell transcriptomic atlas characterizes ageing tissues in the mouse", "10.1038/s41586-020-2496-1"),
    ("Atlas", "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans", "10.1126/science.abl4896"),
    ("Atlas", "Aging Fly Cell Atlas identifies exhaustive aging features at cellular resolution", "10.1126/science.adg0934"),
    ("Atlas", "Aging atlas reveals cell-type-specific effects of pro-longevity strategies", "10.1038/s43587-024-00631-1"),
    ("Atlas", "Individual cell types in C. elegans age differently and activate distinct cell-protective responses", None),
    ("Atlas", "Germ cell progression through zebrafish spermatogenesis declines with age", "10.1242/dev.204319"),
    ("Atlas", "Single-cell transcriptomics reveals age-resistant maintenance of cell identities, stem cell compartments and differentiation trajectories in long-lived naked mole-rats skin", "10.18632/aging.204054"),
    ("Atlas", "Single-cell transcriptomics of the naked mole-rat reveals unexpected features of mammalian immunity", "10.1371/journal.pbio.3000528"),
    ("Atlas", "Lifespan single-cell transcriptomic atlas of the human prefrontal cortex", "10.1038/s41586-026-10271-7"),
    ("Atlas", "An integrated cell atlas of the human lung in health and disease", "10.1038/s41591-023-02327-2"),
    ("Shared programs", "Universal transcriptomic hallmarks of mammalian ageing and mortality", "10.1038/s41586-026-10542-3"),
    ("Shared programs", "Characterization of the conserved mammalian aging signature through single-cell cross-species meta-analysis", None),
    ("Shared programs", "Transcriptional activation of Jun and Fos members of the AP-1 complex is a conserved signature of immune aging that contributes to inflammaging", "10.1111/acel.13792"),
    ("Epigenetic information", "Loss of epigenetic information as a cause of mammalian aging", "10.1016/j.cell.2022.12.027"),
    ("Epigenetic information", "Reprogramming to recover youthful epigenetic information and restore vision", "10.1038/s41586-020-2975-4"),
    ("Transcriptional variability", "Aging increases cell-to-cell transcriptional variability upon immune stimulation", "10.1126/science.aah4115"),
    ("Transcriptional variability", "Single-cell analysis of human pancreas reveals transcriptional signatures of aging and somatic mutation patterns", "10.1016/j.cell.2017.09.004"),
    ("Clock", "Universal DNA methylation age across mammalian tissues", "10.1038/s43587-023-00462-6"),
    ("Clock", "DNA methylation age of human tissues and cell types", "10.1186/gb-2013-14-10-r115"),
    ("Clock", "Genome-wide methylation profiles reveal quantitative views of human aging rates", "10.1016/j.molcel.2012.10.016"),
    ("Clock", "An epigenetic biomarker of aging for lifespan and healthspan", "10.18632/aging.101414"),
    ("Clock", "DNA methylation GrimAge strongly predicts lifespan and healthspan", "10.18632/aging.101684"),
    ("Clock", "A multi-tissue full lifespan epigenetic clock for mice", "10.18632/aging.101590"),
    ("Turnover", "The distribution of cellular turnover in the human body", "10.1038/s41591-020-01182-9"),
    ("Turnover", "Evidence for cardiomyocyte renewal in humans", "10.1126/science.1164680"),
    ("Longevity", "Naked mole-rat mortality rates defy Gompertzian laws by not increasing with age", "10.7554/eLife.31157"),
    ("Senescence", "A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues", "10.1038/s41467-022-32552-1"),
    ("Methods", "propeller: testing for differences in cell type proportions in single cell data", "10.1093/bioinformatics/btac582"),
    ("Methods", "limma powers differential expression analyses for RNA-sequencing and microarray studies", "10.1093/nar/gkv007"),
    ("Methods", "Scanpy: large-scale single-cell gene expression data analysis", "10.1186/s13059-017-1382-0"),
    ("Methods", "The Gene Ontology resource: enriching a GOld mine", "10.1093/nar/gkaa1113"),
    ("Methods", "HAGR: the Human Ageing Genomic Resources", "10.1093/nar/gki017"),
    ("Methods", "Controlling the false discovery rate: a practical and powerful approach to multiple testing", "10.1111/j.2517-6161.1995.tb02031.x"),
]


def normalize(text):
    return " ".join("".join(character.lower() if character.isalnum() else " " for character in text).split())


def metadata(work):
    authors = work.get("author", [])
    author_text = "; ".join(" ".join([a.get("given", ""), a.get("family", "")]).strip() for a in authors[:6])
    if len(authors) > 6:
        author_text += "; et al."
    years = work.get("published-print", work.get("published-online", work.get("issued", {}))).get("date-parts", [[None]])
    return {
        "authors": author_text, "year": years[0][0],
        "journal": "; ".join(work.get("container-title", [])),
        "matched_title": "; ".join(work.get("title", [])),
        "doi": work.get("DOI"), "verification_url": f"https://doi.org/{work.get('DOI')}",
    }


def main():
    rows = []
    for category, title, doi in ENTRIES:
        try:
            if doi:
                response = requests.get(f"https://api.crossref.org/works/{doi}", timeout=60)
                response.raise_for_status()
                work = response.json()["message"]
            else:
                response = requests.get("https://api.crossref.org/works",
                                        params={"query.bibliographic": title, "rows": 1}, timeout=60)
                response.raise_for_status()
                work = response.json()["message"]["items"][0]
                doi = work.get("DOI")
            row = {"category": category, "requested_title": title, **metadata(work)}
            row["title_similarity"] = SequenceMatcher(None, normalize(title), normalize(row["matched_title"])).ratio()
            row["verification_status"] = "verified" if row["title_similarity"] >= 0.85 else "review required"
        except (requests.RequestException, KeyError, IndexError) as error:
            row = {"category": category, "requested_title": title, "authors": "", "year": None,
                   "journal": "", "matched_title": "", "doi": doi, "verification_url": "",
                   "title_similarity": 0.0, "verification_status": f"unresolved: {type(error).__name__}"}
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(RESULTS / "reviewer-phase3-reference-metadata.csv", index=False)
    lines = ["# DOI-verified reference list", "",
             "Each entry below was matched to Crossref metadata with title similarity at least 0.85.", ""]
    for index, row in result[result["verification_status"] == "verified"].iterrows():
        lines.append(f"{index + 1}. {row['authors']} ({int(row['year'])}). {row['matched_title']}. "
                     f"*{row['journal']}*. https://doi.org/{row['doi']}")
    unresolved = result[result["verification_status"] != "verified"]
    lines.extend(["", "## Unresolved manuscript citation", ""])
    for _, row in unresolved.iterrows():
        lines.append(f"- {row['requested_title']}: not verified. Do not reuse the current citation without correction.")
    (REPO / "reports" / "reviewer_phase3_references.md").write_text("\n".join(lines) + "\n")
    print(result[["category", "requested_title", "doi", "title_similarity", "verification_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
