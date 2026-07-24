#!/usr/bin/env bash
#
# make3-pymol.sh
#
# Step 3 of the RBP structural-delineation figure: renders each focal
# RBP's AF3 model (cartoon, grey base colour) with its ECF / ECOD / Pfam
# segments coloured according to pymol-domain-colors.csv (produced at the
# end of make2-plot.R). Framing and ray-tracing settings match convert.py
# so all structure panels share the same look.
#
# Model files are expected to be named "fold_<reprseq_id>.cif", matching
# the AF3 representatives directory, e.g. fold_reprseq019028.cif.
#
# Usage:
#   ./make3-pymol.sh [options]
#
# Options:
#   -t, --table       path to pymol-domain-colors.csv
#                      (default: .../paper-ecf/part6_rbp-struct/pymol-domain-colors.csv)
#   -m, --models-dir  directory containing fold_<reprseq>.cif files
#                      (default: the Fig7 AF3 representatives directory)
#   -o, --out-dir     output directory for .png/.pse
#                      (default: .../paper-ecf/part6_rbp-struct/pymol-out)
#   -r, --resolution  ray-trace resolution in pixels, square (default: 3000)
#   -d, --dpi         PNG dpi (default: 300)
#   --tracks          comma-separated domain.type subset to paint,
#                      e.g. "ECF,ECOD" to skip Pfam (default: all present)
#   --no-ray          skip ray-tracing for a fast draft render
#   --mode            per-segment | combined | both  (default: per-segment)
#                      per-segment: one PNG per protein x segment, that one
#                        segment coloured and everything else grey - use
#                        this if segments are too small/numerous to tell
#                        apart when all painted on the same image at once.
#                      combined:    one PNG per protein, all segments
#                        painted together (the original behaviour).
#                      both:        both of the above.
#                      Every image for a given protein shares the same
#                      camera orientation, so per-segment PNGs line up.
#   --buffer          Angstrom margin around the structure for camera
#                      fitting (default: 15). One shared value for every
#                      protein - viewport is matched to --resolution
#                      before fitting, so this shouldn't need per-structure
#                      tuning even for elongated tail fibres.
#
# Requires PyMOL on PATH (e.g. `conda install -c conda-forge pymol-open-source`,
# or Schrodinger's PyMOL app with its command-line tool linked).
#
# Example:
#   ./make3-pymol.sh --resolution 4000 --tracks ECF,ECOD
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DATA_DIR="/Users/rmostowy/MGG Dropbox/Rafal Mostowy/Projects/Data/ECFs/v4/paper-ecf/part6_rbp-struct"
TABLE="$DATA_DIR/pymol-domain-colors.csv"
MODELS_DIR="/Users/rmostowy/MGG Dropbox/Rafal Mostowy/Projects/Data/ECFs/v4/AF3/Fig7/2026-06/AF3/representatives"
OUT_DIR="$DATA_DIR/pymol-out"
RESOLUTION=3000
DPI=300
TRACKS=""
NO_RAY=0
MODE="per-segment"
BUFFER=15

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--table)      TABLE="$2"; shift 2 ;;
    -m|--models-dir) MODELS_DIR="$2"; shift 2 ;;
    -o|--out-dir)    OUT_DIR="$2"; shift 2 ;;
    -r|--resolution) RESOLUTION="$2"; shift 2 ;;
    -d|--dpi)        DPI="$2"; shift 2 ;;
    --tracks)        TRACKS="$2"; shift 2 ;;
    --no-ray)        NO_RAY=1; shift ;;
    --mode)          MODE="$2"; shift 2 ;;
    --buffer)        BUFFER="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,40p' "$0" | sed 's/^#//; s/^ //'
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if ! command -v pymol >/dev/null 2>&1; then
  echo "Error: 'pymol' not found on PATH." >&2
  exit 1
fi
if [ ! -f "$TABLE" ]; then
  echo "Error: colour table not found at: $TABLE" >&2
  exit 1
fi
if [ ! -d "$MODELS_DIR" ]; then
  echo "Error: models directory not found at: $MODELS_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Embedded PyMOL renderer (written to a temp file so CSV parsing gets a
# real csv module rather than fragile shell/awk splitting - table fields
# such as Pfam descriptions can themselves contain commas).
# ---------------------------------------------------------------------------
TMP_PY="$(mktemp /tmp/make3-pymol.XXXXXX.py)"
trap 'rm -f "$TMP_PY"' EXIT

cat > "$TMP_PY" << 'PYEOF'
import argparse
import csv
import os
import sys
from collections import defaultdict

from pymol import cmd


def parse_args():
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    p = argparse.ArgumentParser()
    p.add_argument("--table", required=True)
    p.add_argument("--models-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--resolution", type=int, default=3000)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--tracks", default=None)
    p.add_argument("--no-ray", action="store_true")
    p.add_argument("--mode", choices=["per-segment", "combined", "both"], default="per-segment")
    p.add_argument("--buffer", type=float, default=15,
                    help="Angstrom margin around the structure passed to cmd.zoom(); "
                         "one shared value for every protein.")
    return p.parse_args(argv)


def load_table(path, tracks):
    by_protein = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"protein", "domain.type", "domain.id", "res.start", "res.end", "color"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit("Table is missing expected columns: %s" % ", ".join(sorted(missing)))
        for row in reader:
            if tracks and row["domain.type"] not in tracks:
                continue
            by_protein[row["protein"]].append(row)
    return by_protein


def hex_to_rgb01(hexcode):
    hexcode = (hexcode or "").strip().lstrip("#")
    if len(hexcode) != 6:
        return None
    return [int(hexcode[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def find_model(models_dir, protein):
    candidate = os.path.join(models_dir, "fold_%s.cif" % protein)
    return candidate if os.path.isfile(candidate) else None


def sanitize(text):
    return str(text).replace(" ", "").replace("/", "-").replace(os.sep, "-")


def render(protein, rows, models_dir, out_dir, resolution, dpi, do_ray, mode, buffer):
    cif = find_model(models_dir, protein)
    if cif is None:
        print("[skip] %s: no fold_%s.cif in %s" % (protein, protein, models_dir))
        return

    cmd.reinitialize()
    cmd.load(cif, protein)
    cmd.hide("everything", protein)
    cmd.show("cartoon", protein)

    # Base colour, same grey as convert.py
    cmd.set_color("custom_grey", [102 / 255.0, 102 / 255.0, 102 / 255.0])
    cmd.color("custom_grey", protein)

    # Match the viewport to the resolution we'll actually ray-trace at
    # *before* orienting/zooming. orient()/zoom() fit the camera to
    # whatever viewport is current; if that differs from the square
    # resolution used later in cmd.ray(), long/thin structures (e.g. tail
    # fibres) get cropped along whichever axis the mismatch shrinks -
    # exactly the "some structures are cut" symptom. Fixing the viewport
    # first makes the fit correct for every structure automatically, so
    # nothing needs per-structure manual retuning.
    cmd.viewport(resolution, resolution)

    # Fix the camera once per protein so every PNG for it - combined or
    # per-segment - is framed identically and can be compared/stacked.
    cmd.orient(protein)
    cmd.center(protein)
    cmd.zoom(protein, buffer=buffer)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ambient", 0.5)
    cmd.set("antialias", 2)

    def snapshot(png_out):
        if do_ray:
            cmd.ray(resolution, resolution)
        cmd.png(png_out, dpi=dpi)

    if mode in ("combined", "both"):
        n_painted = 0
        for i, row in enumerate(rows):
            rgb = hex_to_rgb01(row.get("color"))
            if rgb is None:
                print("[warn] %s: missing/bad colour for %s %s (%s-%s), left grey" %
                      (protein, row["domain.type"], row["domain.id"], row["res.start"], row["res.end"]))
                continue
            color_name = "c_%s_all_%d" % (protein, i)
            cmd.set_color(color_name, rgb)
            cmd.color(color_name, "%s and resi %s-%s" % (protein, row["res.start"], row["res.end"]))
            n_painted += 1

        png_out = os.path.join(out_dir, "%s__all-segments.png" % protein)
        snapshot(png_out)
        cmd.save(os.path.join(out_dir, "%s__all-segments.pse" % protein))
        print("[done] %s combined: %d/%d segments painted -> %s" % (protein, n_painted, len(rows), png_out))

        cmd.color("custom_grey", protein)  # reset before any per-segment renders

    if mode in ("per-segment", "both"):
        for i, row in enumerate(rows):
            cmd.color("custom_grey", protein)  # reset every other segment to grey

            rgb = hex_to_rgb01(row.get("color"))
            if rgb is None:
                print("[warn] %s: missing/bad colour for %s %s (%s-%s), skipped" %
                      (protein, row["domain.type"], row["domain.id"], row["res.start"], row["res.end"]))
                continue

            color_name = "c_%s_seg%d" % (protein, i)
            cmd.set_color(color_name, rgb)
            sel = "%s and resi %s-%s" % (protein, row["res.start"], row["res.end"])
            cmd.color(color_name, sel)

            tag = sanitize("%s_%s_%s-%s" % (row["domain.type"], row["domain.id"],
                                             row["res.start"], row["res.end"]))
            png_out = os.path.join(out_dir, "%s__%s.png" % (protein, tag))
            snapshot(png_out)
            print("[done] %s segment %s -> %s" % (protein, tag, png_out))

        cmd.color("custom_grey", protein)  # leave the object clean afterwards


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    tracks = set(t.strip() for t in args.tracks.split(",")) if args.tracks else None

    by_protein = load_table(args.table, tracks)
    if not by_protein:
        sys.exit("No rows loaded from %s (check --tracks filter)" % args.table)

    for protein, rows in sorted(by_protein.items()):
        render(protein, rows, args.models_dir, args.out_dir,
               args.resolution, args.dpi, not args.no_ray, args.mode, args.buffer)

    cmd.quit()


main()
PYEOF

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
PYMOL_ARGS=( --table "$TABLE" --models-dir "$MODELS_DIR" --out-dir "$OUT_DIR"
             --resolution "$RESOLUTION" --dpi "$DPI" --mode "$MODE" --buffer "$BUFFER" )
[ -n "$TRACKS" ] && PYMOL_ARGS+=( --tracks "$TRACKS" )
[ "$NO_RAY" -eq 1 ] && PYMOL_ARGS+=( --no-ray )

pymol -cq "$TMP_PY" -- "${PYMOL_ARGS[@]}"

echo "Done. Rendered images in: $OUT_DIR"
