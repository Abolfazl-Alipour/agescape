import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_no_controlled_data_or_tokens_in_git():
    # 1. Run git status / ls-files to ensure no sensitive files are tracked
    res = subprocess.run(["git", "ls-files"], cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    tracked_files = res.stdout.splitlines()

    banned_patterns = [
        ".synapseConfig",
        ".env",
        ".rds",
        "syn62064518",
        "syn52396927",
        "lifespan_pseudobulk.rds",
        "m10-dlpfc-whole-score-vs-age.png",
    ]

    for f in tracked_files:
        for b in banned_patterns:
            assert b not in f, f"Banned sensitive file pattern '{b}' tracked in git: {f}"


def test_zero_innerhtml_outerhtml_in_app():
    app_js = (REPO_ROOT / "app" / "app.js").read_text()
    assert "innerHTML" not in app_js, "innerHTML found in app/app.js!"
    assert "outerHTML" not in app_js, "outerHTML found in app/app.js!"
    assert "document.write" not in app_js, "document.write found in app/app.js!"


def test_zero_inline_styles_in_index_html():
    index_html = (REPO_ROOT / "app" / "index.html").read_text()
    import re
    inline_styles = re.findall(r'style="[^"]*"', index_html)
    assert len(inline_styles) == 0, f"Inline styles found in app/index.html: {inline_styles}"


def test_no_third_party_fonts_in_index_html():
    index_html = (REPO_ROOT / "app" / "index.html").read_text()
    assert "fonts.googleapis.com" not in index_html, "External Google Fonts link found in app/index.html!"
    assert "fonts.gstatic.com" not in index_html, "External gstatic link found in app/index.html!"
    assert (REPO_ROOT / "app" / "vendor" / "fonts" / "fonts.css").exists(), "Vendored fonts.css missing!"


def test_csp_headers_snippet_exists():
    snippet = REPO_ROOT / "docs" / "cloudflare-_headers.snippet"
    assert snippet.exists(), "CSP headers snippet missing!"
    content = snippet.read_text()
    assert "Content-Security-Policy" in content, "CSP header missing in snippet!"
    assert "default-src 'self'" in content, "default-src 'self' missing in CSP snippet!"


def test_no_local_user_paths_in_tracked_manifests():
    import json
    for p in (REPO_ROOT / "data" / "manifest").glob("*.json"):
        text = p.read_text()
        assert "/home/reza" not in text, f"Local path /home/reza found in manifest {p}!"
    for p in (REPO_ROOT / "results").glob("*.log"):
        text = p.read_text()
        assert "/home/reza" not in text, f"Local path /home/reza found in log {p}!"
