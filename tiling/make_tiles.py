#!/usr/bin/env python3
"""
Cut a large multi-channel OME-TIFF into a grid of tiles for CellProfiler.
Writes tile_r{row}_c{col}.ome.tiff files with the channel axis labeled T,
so the same pipeline reads them like the test crop.

Usage:
    python3 make_tiles.py INPUT.ome.tiff OUTPUT_DIR [--tile 6000] [--overlap 200]

--tile     tile size in px (square). 6000 is a safe size for 48GB RAM.
           NOTE: production runs (see tiling/tile.sh) call this with
           --tile 4000 --overlap 100, which is what the rest of the
           pipeline (hpc/, analysis/merge_brain.py) was validated against.
           The 6000/200 defaults here are just this script's own fallback
           if called directly without those flags.
--overlap  px of overlap between adjacent tiles (avoids cutting cells at seams;
           dedupe later using tile origin + object centroid).
"""
import argparse, os, json
import numpy as np
import tifffile

p = argparse.ArgumentParser()
p.add_argument("inp")
p.add_argument("outdir")
p.add_argument("--tile", type=int, default=6000)
p.add_argument("--overlap", type=int, default=200)
a = p.parse_args()

os.makedirs(a.outdir, exist_ok=True)

with tifffile.TiffFile(a.inp) as tif:
    series = tif.series[0]
    arr = series.asarray()
    axes = series.axes

print(f"Loaded shape={arr.shape} axes={axes} dtype={arr.dtype}")
yi, xi = axes.index("Y"), axes.index("X")
H, W = arr.shape[yi], arr.shape[xi]
step = a.tile - a.overlap

manifest = []
row = 0
for y0 in range(0, H, step):
    col = 0
    for x0 in range(0, W, step):
        y1, x1 = min(H, y0 + a.tile), min(W, x0 + a.tile)
        # skip tiny edge slivers
        if (y1 - y0) < a.tile * 0.3 or (x1 - x0) < a.tile * 0.3:
            col += 1; continue
        sl = [slice(None)] * arr.ndim
        sl[yi], sl[xi] = slice(y0, y1), slice(x0, x1)
        crop = arr[tuple(sl)]
        # skip near-empty tiles (no tissue) to save time — check DAPI (plane 0)
        if crop[0].mean() < 1.0:
            print(f"  skip empty r{row}_c{col} (Y{y0} X{x0})")
            col += 1; continue
        fn = os.path.join(a.outdir, f"tile_r{row}_c{col}.ome.tiff")
        tifffile.imwrite(fn, crop, photometric="minisblack",
                         metadata={"axes": "TYX"}, compression="zlib")
        manifest.append({"file": os.path.basename(fn), "row": row, "col": col,
                         "x0": x0, "y0": y0, "x1": x1, "y1": y1})
        print(f"  wrote r{row}_c{col}  Y[{y0}:{y1}] X[{x0}:{x1}]  shape={crop.shape}")
        col += 1
    row += 1

with open(os.path.join(a.outdir, "tile_manifest.json"), "w") as f:
    json.dump({"source": a.inp, "tile": a.tile, "overlap": a.overlap,
               "full_H": H, "full_W": W, "tiles": manifest}, f, indent=2)
print(f"\nWrote {len(manifest)} tiles + tile_manifest.json to {a.outdir}")
