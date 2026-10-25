#!/usr/bin/env bash
# Sync the deployable Agescape web app into a personal-website repo as a subpath.
# The app is fully self-contained and relative-pathed (verified 2026-10-19), so it
# works under any subpath, e.g. https://<your-domain>/agescape/.
#
# Usage: bash scripts/sync_to_website.sh /path/to/website-repo [subpath]
#   subpath defaults to "agescape". Then, in the website repo:
#     git add <subpath> && git commit && git push   (push = publish via Cloudflare)
set -euo pipefail

DEST_ROOT="${1:?usage: sync_to_website.sh /path/to/website-repo [subpath]}"
SUBPATH="${2:-agescape}"

[ -d "$DEST_ROOT/.git" ] || { echo "error: $DEST_ROOT is not a git repo"; exit 1; }
[ -f app/index.html ] || { echo "error: run from the agescape repo root"; exit 1; }

mkdir -p "$DEST_ROOT/$SUBPATH"
rsync -a --delete --exclude 'DATA_CONTRACT.md' --exclude 'provenance.json' \
  app/index.html app/app.js app/style.css app/favicon.svg app/vendor app/data app/figures \
  "$DEST_ROOT/$SUBPATH/"

echo "synced app/ -> $DEST_ROOT/$SUBPATH"
echo "next: cd $DEST_ROOT && git status --short $SUBPATH && git add $SUBPATH"
echo "PUBLISH HAPPENS ON PUSH — confirm with Reza before 'git push'."
