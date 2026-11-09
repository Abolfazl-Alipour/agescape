#!/usr/bin/env python3
"""Check structural requirements for the PLOS ONE manuscript source and PDF."""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SOURCE = REPO / "report" / "plos-one" / "final-plos-one.tex"
PDF = SOURCE.with_suffix(".pdf")


def main():
    text = SOURCE.read_text()
    abstract = text.split(r"\section*{Abstract}", 1)[1].split(r"\section*{Introduction}", 1)[0]
    plain = re.sub(r"\\[a-zA-Z]+|[{}]", " ", abstract)
    words = re.findall(r"\b[\w'-]+\b", plain)
    assert len(words) <= 300, f"Abstract has {len(words)} words"
    assert not re.search(r"\b[qp]\s*=|\b[+-]?\d+\.\d+\b", abstract, re.I)
    assert "—" not in text and "–" not in text
    assert "Figure 5" not in text and "S1" not in text

    body, references = text.split(r"\section*{References}", 1)
    cited = {int(n) for group in re.findall(r"\[([0-9, ]+)\]", body)
             for n in group.replace(" ", "").split(",")}
    listed = {int(n) for n in re.findall(r"^([0-9]+)\.", references, re.M)}
    assert cited == listed, f"Citation mismatch: uncited={listed-cited}, missing={cited-listed}"

    figures = re.findall(r"\\includegraphics\[[^]]*\]\{([^}]+)\}", text)
    assert len(figures) == 4, f"Expected four figures, found {len(figures)}"
    for figure in figures:
        assert (SOURCE.parent / figure).is_file(), f"Missing figure {figure}"

    sections = re.findall(r"\\section\*\{([^}]+)\}", text)
    expected = ["Abstract", "Introduction", "Materials and Methods", "Results", "Discussion",
                "Conclusions", "Data Availability", "Funding", "Competing Interests",
                "Ethics Statement", "Supporting Information", "References"]
    assert sections == expected, f"Section order mismatch: {sections}"

    residue = [r"reviewer response", r"Phase [0-3]", r"first draft", r"\b[a-f0-9]{40}\b",
               r"\.csv(?::|\b)", r"still running", r"verdict pending", r"reanalysis is pending",
               r"\bnot redone\b", r"\breference.only\b", r"\bsuperseded\b",
               r"\boriginal analysis\b", r"\bthis revision\b"]
    for pattern in residue:
        assert not re.search(pattern, text, re.I), f"Process residue: {pattern}"
    assert "Affiliation: Kavox Labs" in text
    assert "The computational analyses were performed on the Kavox Labs cloud computing infrastructure (Amazon Web Services)." in text
    assert "The author is the founder of Kavox Labs." in text
    assert text.count("Pending author confirmation.") == 0
    assert "Draft statement, pending author confirmation." in text

    pages = None
    try:
        from pypdf import PdfReader
        pages = len(PdfReader(str(PDF)).pages)
    except ImportError:
        pass
    print(f"PASS: abstract {len(words)} words; four figures; {len(listed)} cited references")
    print("PASS: section order, no forbidden dashes, no listed process residue")
    print(f"PDF pages: {pages if pages is not None else 'pypdf unavailable'}")


if __name__ == "__main__":
    sys.exit(main())
