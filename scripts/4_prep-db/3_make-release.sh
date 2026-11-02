#!/usr/bin/env bash
#
# 3_make-release.sh
#
# Packages the ECF database (built by 1_process-db.R and 2_ecf-metadata.R)
# into release files:
#
#   ecf-hhsuitedb.tar.gz   HH-suite database (hhsuitedb/)
#   ecf-a3m.tar.gz         per-ECF A3M alignments (a3m/)
#   ecf-metadata.csv       ECF metadata table
#   SHA256SUMS.txt         SHA-256 checksums of the three files above
#
# Usage:
#   ./3_make-release.sh [-d DB_DIR] [-o OUT_DIR]
#
#   -d DB_DIR   database directory (config.R: database.dir)
#               default: $HOME/Dropbox/Projects/Data/ECFs/v4/output/database
#   -o OUT_DIR  where the release files are written
#               default: release/ next to DB_DIR, i.e.
#               $HOME/Dropbox/Projects/Data/ECFs/v4/output/release
#
# DB_DIR and OUT_DIR can also be set through the environment variables
# DB and OUT.

set -euo pipefail

DB="${DB:-$HOME/Dropbox/Projects/Data/ECFs/v4/output/database}"
OUT="${OUT:-}"

while getopts "d:o:h" opt; do
  case "$opt" in
    d) DB="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    h) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Usage: $0 [-d DB_DIR] [-o OUT_DIR]" >&2; exit 1 ;;
  esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

# default output: output/release, a sibling of output/database
[ -n "$OUT" ] || OUT="$(dirname "$DB")/release"

# ---------------------------------------------------------------------------
# Checks on the input database
# ---------------------------------------------------------------------------
[ -d "$DB/hhsuitedb" ]       || die "missing directory: $DB/hhsuitedb"
[ -d "$DB/a3m" ]             || die "missing directory: $DB/a3m"
[ -f "$DB/ecf-metadata.csv" ] || die "missing file: $DB/ecf-metadata.csv"

# 1_process-db.R writes 6 files (a3m / hhm / cs219, each as ffdata + ffindex)
n_ff=$(find "$DB/hhsuitedb" -maxdepth 1 -name 'ecf_db_*.ff*' | wc -l | tr -d ' ')
[ "$n_ff" -eq 6 ] || die "expected 6 ecf_db_*.ffdata/ffindex files in $DB/hhsuitedb, found $n_ff"

n_a3m=$(find "$DB/a3m" -maxdepth 1 -name '*.a3m' | wc -l | tr -d ' ')
[ "$n_a3m" -gt 0 ] || die "no .a3m files in $DB/a3m"

# Every ECF in the database should have an entry in the metadata table
# (the table has a header line, hence the -1)
n_meta=$(( $(wc -l < "$DB/ecf-metadata.csv" | tr -d ' ') - 1 ))
[ "$n_meta" -eq "$n_a3m" ] || echo "WARNING: $n_a3m a3m files but $n_meta rows in ecf-metadata.csv" >&2

# ---------------------------------------------------------------------------
# Build the release files
# ---------------------------------------------------------------------------
mkdir -p "$OUT"
cd "$OUT"

# macOS files carry extended attributes (the '@' in ls -l); this keeps them
# (as ._* entries) out of the tarballs
export COPYFILE_DISABLE=1

rm -f ecf-hhsuitedb.tar.gz ecf-a3m.tar.gz ecf-metadata.csv SHA256SUMS.txt

echo "Packing hhsuitedb ..."
tar -czf ecf-hhsuitedb.tar.gz -C "$DB" hhsuitedb
echo "Packing a3m ($n_a3m files) ..."
tar -czf ecf-a3m.tar.gz -C "$DB" a3m
cp "$DB/ecf-metadata.csv" .

# ---------------------------------------------------------------------------
# Checksums
# ---------------------------------------------------------------------------
FILES=(ecf-a3m.tar.gz ecf-hhsuitedb.tar.gz ecf-metadata.csv)
if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "${FILES[@]}" > SHA256SUMS.txt
else
  sha256sum "${FILES[@]}" > SHA256SUMS.txt
fi

# ---------------------------------------------------------------------------
# Verify the archives and the checksums
# ---------------------------------------------------------------------------
tar -tzf ecf-hhsuitedb.tar.gz >/dev/null || die "ecf-hhsuitedb.tar.gz is not a valid archive"
tar -tzf ecf-a3m.tar.gz       >/dev/null || die "ecf-a3m.tar.gz is not a valid archive"
if tar -tzf ecf-hhsuitedb.tar.gz ecf-a3m.tar.gz 2>/dev/null | grep -q '/\._'; then
  die "archives contain macOS ._* metadata entries"
fi
if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 -c SHA256SUMS.txt
else
  sha256sum -c SHA256SUMS.txt
fi

echo
echo "Release files in $OUT:"
ls -lh "${FILES[@]}" SHA256SUMS.txt
