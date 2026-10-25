import urllib.request, re, os
from pathlib import Path

out_dir = Path("app/vendor/fonts")
out_dir.mkdir(parents=True, exist_ok=True)

url = "https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600&family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,500;0,6..72,600;0,6..72,700;1,6..72,400;1,6..72,600&family=Plus+Jakarta+Sans:wght@400;500;600;700;800&display=swap"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"})
css = urllib.request.urlopen(req).read().decode("utf-8")

# Extract font blocks
blocks = re.findall(r"/\* ([^*]+) \*/\s*(@font-face\s*\{[^}]+\})", css)
latin_blocks = [b for comment, b in blocks if "latin" in comment.lower() and "ext" not in comment.lower()]

new_css = []
for block in latin_blocks:
    urls = re.findall(r"url\((https://[^)]+)\)", block)
    for u in urls:
        fname = u.split("/")[-1]
        local_path = out_dir / fname
        if not local_path.exists():
            print(f"Downloading {fname}...")
            urllib.request.urlretrieve(u, local_path)
        block = block.replace(u, fname)
    new_css.append(block)

with open(out_dir / "fonts.css", "w") as f:
    f.write("\n\n".join(new_css))

woff_count = len(list(out_dir.glob("*.woff2")))
print(f"Vendored fonts successfully: {woff_count} woff2 files in {out_dir}.")
