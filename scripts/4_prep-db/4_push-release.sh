#!/usr/bin/env bash
#
# 4_push-release.sh
#
# Uploads the release files built by 3_make-release.sh to GitHub as a release
# of this repository (using the GitHub CLI, `gh`):
#
#   ecf-hhsuitedb.tar.gz, ecf-a3m.tar.gz, ecf-metadata.csv, SHA256SUMS.txt
#
# The release text is taken from release-notes/<version>.md in this repository.
# By default the release is created as a DRAFT, so it can be checked on GitHub
# and published by hand.
#
# Usage:
#   ./4_push-release.sh VERSION [-d RELEASE_DIR] [-n NOTES_FILE] [-p]
#
#   VERSION       release tag, e.g. v1.0.0
#   -d RELEASE_DIR  directory with the release files (output of 3_make-release.sh)
#                   default: $HOME/Dropbox/Projects/Data/ECFs/v4/output/release
#   -n NOTES_FILE   release notes (Markdown)
#                   default: <repo>/release-notes/VERSION.md
#   -p              publish immediately instead of creating a draft
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

[ $# -ge 1 ] || die "usage: $0 VERSION [-d RELEASE_DIR] [-n NOTES_FILE] [-p]"
case "$1" in
  -h|--help) sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
esac
VERSION="$1"; shift
[[ "$VERSION" =~ ^v[0-9]+\.[0-9]+\.[0-9]+([-.][0-9A-Za-z.]+)?$ ]] \
  || die "VERSION should look like v1.0.0 (got '$VERSION')"

while getopts "d:n:ph" opt; do
  case "$opt" in
    d) OUT="$OPTARG" ;;
    n) NOTES="$OPTARG" ;;
    p) DRAFT=() ;;
    h) sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "usage: $0 VERSION [-d RELEASE_DIR] [-n NOTES_FILE] [-p]" ;;
  esac
done
[ -n "$NOTES" ] || NOTES="$REPO/release-notes/$VERSION.md"

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
command -v gh >/dev/null 2>&1 || die "GitHub CLI 'gh' not found (https://cli.github.com)"
gh auth status >/dev/null 2>&1 || die "gh is not logged in; run: gh auth login"

[ -f "$NOTES" ] || die "release notes not found: $NOTES"

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

if gh release view "$VERSION" >/dev/null 2>&1; then
  die "a release with tag $VERSION already exists on GitHub"
fi

# The tag is created by GitHub from the last pushed commit of the default
# branch, so uncommitted or unpushed changes would not be part of the release.
if [ -n "$(git status --porcelain)" ]; then
  echo "WARNING: the working tree has uncommitted changes" >&2
fi
if [ -n "$(git rev-list '@{u}..HEAD' 2>/dev/null)" ]; then
  echo "WARNING: there are local commits that have not been pushed" >&2
fi

# ---------------------------------------------------------------------------
# Create the release
# ---------------------------------------------------------------------------
ASSETS=()
for f in "${FILES[@]}"; do ASSETS+=("$OUT/$f"); done

echo "Creating release $VERSION ${DRAFT[*]:-(published)} ..."
gh release create "$VERSION" \
  "${ASSETS[@]}" \
  --title "ECF database $VERSION" \
  --notes-file "$NOTES" \
  ${DRAFT[@]+"${DRAFT[@]}"}
