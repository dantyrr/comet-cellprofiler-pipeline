#!/bin/bash
# Tile every COMET image in Nick_Raw_Images_8-26-26 into <ID>_tiles_4k/ folders.
#
# Usage:
#   bash tile_all.sh                 # tile all 8 images
#   bash tile_all.sh ICV-C1_3 IP-T2_3   # tile only these sample IDs
#
# Tiles are 4000 px with 100 px overlap, matching the validated ICV run --
# merge_brain.py's dedup assumes the overlap recorded in tile_manifest.json.
#
# Each image is tiled and reported one at a time so you can upload and delete
# as you go; all 8 at once needs well over 100 GB of free disk.

set -euo pipefail

PROJECT="${COMET_PROJECT_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
IMGDIR="${COMET_IMAGE_DIR:?Set COMET_IMAGE_DIR to the folder holding the whole-slide OME-TIFFs}"
# Interpreter with tifffile+numpy. Override with COMET_PYTHON.
PYTHON="${COMET_PYTHON:-python3}"
SIZE=4000
OVERLAP=100

[ -x "$PYTHON" ] || { echo "ERROR: python not found at $PYTHON"; exit 1; }
[ -d "$IMGDIR" ] || { echo "ERROR: image dir not found: $IMGDIR"; exit 1; }

if [ $# -gt 0 ]; then
    IDS=("$@")
else
    IDS=()
    for f in "$IMGDIR"/mBrain\ Final\ runs_*.ome.tiff; do
        b=$(basename "$f" .ome.tiff)
        IDS+=("${b#mBrain Final runs_}")
    done
fi

echo "Project: $PROJECT"
echo "Samples: ${IDS[*]}"
echo "Tiles:   ${SIZE} px, ${OVERLAP} px overlap"
echo

for ID in "${IDS[@]}"; do
    IN="$IMGDIR/mBrain Final runs_${ID}.ome.tiff"
    OUTID=$(echo "$ID" | tr '-' '_')          # ICV-C1_3 -> ICV_C1_3
    OUT="$PROJECT/${OUTID}_tiles_4k"

    echo "=============================================================="
    echo "$ID  ->  $(basename "$OUT")"

    if [ ! -f "$IN" ]; then
        echo "  SKIP: not found: $IN"
        continue
    fi
    if [ ! -s "$IN" ]; then
        echo "  SKIP: file is 0 bytes (still copying from Box?)"
        continue
    fi
    if [ -f "$OUT/tile_manifest.json" ]; then
        N=$(ls "$OUT"/tile_r*_c*.ome.tiff 2>/dev/null | wc -l | tr -d ' ')
        echo "  SKIP: already tiled ($N tiles). Delete $OUT to redo."
        continue
    fi

    echo "  Verifying channel map before tiling..."
    if ! "$PYTHON" "$PROJECT/tools/verify_channel_map.py" "$IN"; then
        echo "  ABORT $ID: channel map does not match the pipeline. Not tiling."
        echo "  Fix the plane indices in the .cppipe before continuing."
        continue
    fi

    "$PYTHON" "$PROJECT/tiling/make_tiles.py" "$IN" "$OUT" --tile "$SIZE" --overlap "$OVERLAP"

    N=$(ls "$OUT"/tile_r*_c*.ome.tiff 2>/dev/null | wc -l | tr -d ' ')
    SZ=$(du -sh "$OUT" | cut -f1)
    echo "  Done: $N tiles, $SZ"
    echo "  Free disk: $(df -h "$PROJECT" | tail -1 | awk '{print $4}')"
    echo
done

echo "=============================================================="
echo "Tiling complete. Upload the *_tiles_4k folders to Cheaha, then run:"
echo "  bash submit_all_brains.sh ICV_C1_3 ICV_C3_3 ICV_T1_3 ICV_T2_3 \\"
echo "                            IP_C1_3 IP_C2_3 IP_T2_3 IP_T3_3"
