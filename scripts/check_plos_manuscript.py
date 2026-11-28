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
    # 2026-11-23: Reza dropped the 300-word abstract cap ("I don't care how
    # long the abstract gets"). The count is still printed, not enforced.
    # 2026-11-09: Reza ordered p-values and Bayesian results INTO the abstract; the no-statistics rule is retired.
    assert "—" not in text and "–" not in text
    # Three main figures (transport Figure D1 is listed as a file in Appendix D, 2026-11-28). Supporting figures (S1 Fig, S2 Fig, ...) are
    # allowed only in the Supporting Information section, where they are
    # listed as files (figures-si/...) with captions, not as \includegraphics.
    assert "Figure 5" not in text
    main_body, si_section = text.split(r"\section*{Supporting Information}", 1)
    si_section = si_section.split(r"\section*{References}", 1)[0]
    assert not re.search(r"\bS\d+ Fig\b", main_body), "SI reference outside SI section"
    assert not re.search(r"\bS\d+\b", main_body), "SI label outside SI section"

    body, references = text.split(r"\section*{References}", 1)
    cited = {int(n) for group in re.findall(r"\[([0-9, ]+)\]", body)
             for n in group.replace(" ", "").split(",")}
    listed = {int(n) for n in re.findall(r"^([0-9]+)\.", references, re.M)}
    assert cited == listed, f"Citation mismatch: uncited={listed-cited}, missing={cited-listed}"

    figures = re.findall(r"\\includegraphics\[[^]]*\]\{([^}]+)\}", main_body)
    assert len(figures) == 3, f"Expected three main figures, found {len(figures)}"
    for figure in figures:
        assert (SOURCE.parent / figure).is_file(), f"Missing figure {figure}"
    assert r"\includegraphics" not in si_section, "SI figures must be listed as files, not embedded"
    # Every SI figure file listed in the SI section must exist (PDF and PNG).
    si_files = re.findall(r"\(figures-si/([A-Za-z0-9_.\-]+\.pdf)\)", si_section)
    assert si_files, "No SI figure files listed in the Supporting Information section"
    for name in si_files:
        pdf = SOURCE.parent / "figures-si" / name
        assert pdf.is_file(), f"Missing SI figure {pdf}"
        assert pdf.with_suffix(".png").is_file(), f"Missing SI figure PNG {pdf.with_suffix('.png')}"
    si_refs = [int(n) for n in re.findall(r"\bS(\d+) Fig\b", si_section)]
    assert si_refs == list(range(1, len(si_refs) + 1)), f"SI references not sequential from S1: {si_refs}"
    assert len(si_refs) == len(si_files), "SI reference count does not match listed SI files"

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
    assert "pending author confirmation" not in text  # ethics confirmed by author 2026-11-09

    pages = None
    try:
        from pypdf import PdfReader
        pages = len(PdfReader(str(PDF)).pages)
    except ImportError:
        pass
    print(f"PASS: abstract {len(words)} words; three figures; {len(listed)} cited references")
    print(f"PASS: {len(si_files)} SI figures listed and present; section order, no forbidden dashes, no listed process residue")
    print(f"PDF pages: {pages if pages is not None else 'pypdf unavailable'}")


if __name__ == "__main__":
    sys.exit(main())
