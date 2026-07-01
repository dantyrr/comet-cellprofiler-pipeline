#!/bin/bash
# Tile a COMET image by sample ID.
# Usage:  bash tile.sh SAMPLE_ID [tile_size] [overlap]
#   e.g.  bash tile.sh ICV-T2_3
#         bash tile.sh IP-C1_3 4000 100
#
# Finds "mBrain Final runs_<SAMPLE_ID>.ome.tiff" anywhere under the Comet
# data folder, tiles it, and writes <output root>/<SAMPLE_ID>_tiles_4k/.
#
# This script runs locally (e.g. on your Mac) as a pre-processing step
# before uploading tiles to the HPC cluster — it is NOT part of the Slurm
# pipeline in hpc/. It depends on tiling/make_tiles.py, which is a
# separate stub in this repo (see tiling/README.md) until the real
# implementation is added here alongside this script.

set -e

ID="$1"
SIZE="${2:-4000}"
OVERLAP="${3:-100}"

if [ -z "$ID" ]; then
  echo "Usage: bash tile.sh SAMPLE_ID [tile_size] [overlap]"
  echo "  e.g. bash tile.sh ICV-T2_3"
  exit 1
fi

# --- config ---
# Defaults below are dtyrrell's local Mac setup. Override via env vars if
# the data location, script location, or output location ever change.
BASE="${COMET_BOX_DATA_ROOT:-/Users/dtyrrell/Library/CloudStorage/Box-Box/00_Tyrrell Lab/04_DATA/01_OUR_DATA/04_Comet/Nick's Comet data}"
SCRIPT="${MAKE_TILES_SCRIPT:-$(dirname "$0")/make_tiles.py}"
OUT_ROOT="${COMET_TILES_OUTPUT_ROOT:-$HOME/Desktop}"
# --------------

if [ ! -f "$SCRIPT" ]; then
  echo "ERROR: make_tiles.py not found at $SCRIPT"
  echo "It's a stub in this repo until the real implementation is added —"
  echo "see tiling/README.md. Set MAKE_TILES_SCRIPT to point elsewhere if needed."
  exit 1
fi

echo "Searching for sample '$ID' ..."
IN=$(find "$BASE" -name "mBrain Final runs_${ID}.ome.tiff" -print -quit 2>/dev/null)

if [ -z "$IN" ]; then
  echo "ERROR: could not find 'mBrain Final runs_${ID}.ome.tiff' under:"
  echo "  $BASE"
  echo "Check the sample ID (and that Box has materialized the file), or set"
  echo "COMET_BOX_DATA_ROOT to the correct data folder."
  exit 1
fi

OUTID=$(echo "$ID" | tr '-' '_')
OUT="$OUT_ROOT/${OUTID}_tiles_4k"

echo "Found:  $IN"
echo "Output: $OUT  (${SIZE}px tiles, ${OVERLAP}px overlap)"
echo "-----------------------------------------------------"

python3 "$SCRIPT" "$IN" "$OUT" --tile "$SIZE" --overlap "$OVERLAP"

echo "-----------------------------------------------------"
N=$(ls "$OUT"/*.ome.tiff 2>/dev/null | wc -l | tr -d ' ')
echo "Done. Wrote $N tiles to $OUT"
echo "Next: upload the folder '$OUT' to Cheaha via the OnDemand Files app."
