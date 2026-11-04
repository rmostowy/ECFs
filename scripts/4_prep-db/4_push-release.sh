#!/usr/bin/env bash
#
# 4_push-release.sh
#
# Creates a GitHub release of this repository (using the GitHub CLI, `gh`) with
# the database files built by 3_make-release.sh attached:
#
#   ecf-hhsuitedb.tar.gz, ecf-a3m.tar.gz, ecf-metadata.csv, SHA256SUMS.txt
#
# Code and database share ONE version series (see the README): the release is
# the git tag VERSION, so its source archive (the code) and the database files
# belong to the same release. Publishing the release makes Zenodo archive the
# source code of the tag and issue a DOI for this version.
#
# The release text is taken from release-notes/<version>.md in this repository.
# By default the release is created as a DRAFT, so it can be checked on GitHub
# and published by hand. Nothing is archived on Zenodo until it is published.
#
# The commit that gets tagged is the current commit (HEAD), which must be
# committed and pushed first, so that the notes, README, CITATION.cff and
# scripts of the release are in the archive.
#
# Usage:
#   ./4_push-release.sh VERSION [-d RELEASE_DIR] [-n NOTES_FILE] [-p] [-f]
#
#   VERSION       release tag, e.g. v1.0.0
#   -d RELEASE_DIR  directory with the release files (output of 3_make-release.sh)
#                   default: $HOME/Dropbox/Projects/Data/ECFs/v4/output/release
#   -n NOTES_FILE   release notes (Markdown)
#                   default: <repo>/release-notes/VERSION.md
#   -p              publish immediately instead of creating a draft
#   -f              proceed although the working tree is not clean or there
#                   are unpushed commits
#
# RELEASE_DIR can also be set through the environment variable OUT.
#
# Requires: gh (https://cli.github.com), logged in with `gh auth login`.

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

# The repository root is found from the location of this script, so the
# script can be started from any directory.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

OUT="${OUT:-$HOME/Dropbox/Projects/Data/ECFs/v4/output/release}"
NOTES=""
DRAFT=(--draft)
FORCE=0

[ $# -ge 1 ] || die "usage: $0 VERSION [-d RELEASE_DIR] [-n NOTES_FILE] [-p] [-f]"
case "$1" in
  -h|--help) sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
esac
VERSION="$1"; shift
[[ "$VERSION" =~ ^v[0-9]+\.[0-9]+\.[0-9]+([-.][0-9A-Za-z.]+)?$ ]] \
  || die "VERSION should look like v1.0.0 (got '$VERSION')"

while getopts "d:n:pfh" opt; do
  case "$opt" in
    d) OUT="$OPTARG" ;;
    n) NOTES="$OPTARG" ;;
    p) DRAFT=() ;;
    f) FORCE=1 ;;
    h) sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "usage: $0 VERSION [-d RELEASE_DIR] [-n NOTES_FILE] [-p] [-f]" ;;
  esac
done
[ -n "$NOTES" ] || NOTES="$REPO/release-notes/$VERSION.md"

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
command -v gh >/dev/null 2>&1 || die "GitHub CLI 'gh' not found (https://cli.github.com)"
gh auth status >/dev/null 2>&1 || die "gh is not logged in; run: gh auth login"

[ -f "$NOTES" ] || die "release notes not found: $NOTES"
if grep -n "<TODO" "$NOTES" >/dev/null; then
  grep -n "<TODO" "$NOTES" >&2
  die "release notes still contain <TODO ...> placeholders: $NOTES"
fi

FILES=(ecf-hhsuitedb.tar.gz ecf-a3m.tar.gz ecf-metadata.csv SHA256SUMS.txt)
for f in "${FILES[@]}"; do
  [ -f "$OUT/$f" ] || die "missing release file: $OUT/$f (run 3_make-release.sh first)"
  # GitHub accepts release assets of up to 2 GiB each
  [ "$(wc -c < "$OUT/$f")" -lt 2147483648 ] || die "$f is larger than GitHub's 2 GiB asset limit"
done

# The files must still match their checksums
echo "Verifying checksums ..."
( cd "$OUT"
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 -c SHA256SUMS.txt
  else
    sha256sum -c SHA256SUMS.txt
  fi
) || die "checksum verification failed; rebuild the release with 3_make-release.sh"

# gh works out the target repository from the git remote of the current directory
cd "$REPO"

# `gh release view` does not find drafts, so list all releases (drafts included)
if gh release list --limit 1000 --json tagName --jq '.[].tagName' | grep -qx -- "$VERSION"; then
  die "a release (possibly a draft) with tag $VERSION already exists on GitHub; check: gh release list"
fi
if git rev-parse -q --verify "refs/tags/$VERSION" >/dev/null; then
  die "a local git tag $VERSION already exists"
fi

# The release is tagged on the current commit, which GitHub must already have:
# uncommitted or unpushed changes would not be part of the release (nor of the
# Zenodo archive).
PROBLEMS=0
if [ -n "$(git status --porcelain)" ]; then
  echo "The working tree has uncommitted or untracked changes:" >&2
  git status --short >&2
  PROBLEMS=1
fi
if [ -n "$(git rev-list '@{u}..HEAD' 2>/dev/null)" ]; then
  echo "There are local commits that have not been pushed." >&2
  PROBLEMS=1
fi
if [ "$PROBLEMS" -eq 1 ]; then
  if [ "$FORCE" -eq 1 ]; then
    echo "WARNING: continuing because of -f" >&2
  else
    die "commit and push first (or use -f)"
  fi
fi
TARGET="$(git rev-parse HEAD)"

# ---------------------------------------------------------------------------
# Create the release
# ---------------------------------------------------------------------------
ASSETS=()
for f in "${FILES[@]}"; do ASSETS+=("$OUT/$f"); done

echo "Creating release $VERSION ${DRAFT[*]:-(published)} on commit ${TARGET:0:7} ..."
gh release create "$VERSION" \
  "${ASSETS[@]}" \
  --title "ECFs $VERSION: analysis code and ECF database" \
  --notes-file "$NOTES" \
  --target "$TARGET" \
  ${DRAFT[@]+"${DRAFT[@]}"}

echo
if [ ${#DRAFT[@]} -gt 0 ]; then
  echo "The release is a draft: check it on GitHub, then publish it."
  echo "Publishing creates the tag $VERSION and makes Zenodo archive the code (new version DOI)."
else
  echo "The release is published; Zenodo will now archive it (new version DOI)."
fi
