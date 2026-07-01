# Tiling — TODO, not included in this handoff

`make_tiles.py` and `tile.sh` in this directory are **stubs**. The real
tiling script that produced the validated production tiles lives on the
researcher's local machine and wasn't available when this repo was
packaged. Ask the repo owner (dtyrrell) for it, or write your own — but
if you write your own, it **must** satisfy the output contract below,
since `hpc/run_array.sh` and `analysis/merge_brain.py` depend on it
exactly.

## Output contract

```
Input:  one large multi-channel OME-TIFF (COMET whole-slide image)
Output: a directory named <image_name>_tiles_4k/ containing:
  - tile_r{row}_c{col}.ome.tiff  (one file per tile, all channels preserved)
  - tile_manifest.json           (see schema below)
```

### `tile_manifest.json` schema

```json
{
  "source": "<absolute path to original whole-slide OME-TIFF>",
  "tile": 4000,
  "overlap": 100,
  "full_H": 0,
  "full_W": 0,
  "tiles": [
    {
      "file": "tile_r0_c0.ome.tiff",
      "row": 0, "col": 0,
      "x0": 0, "y0": 0,
      "x1": 0, "y1": 0
    }
  ]
}
```

- `tile`: tile edge length in pixels (4000 px in validated production runs).
- `overlap`: overlap between adjacent tiles, in pixels (100 px in validated
  production runs).
- `full_H` / `full_W`: full source image height/width in pixels.
- `tiles[].x0,y0,x1,y1`: tile's bounding box in **global image coordinates**
  (top-left / bottom-right corner).

## Overlap requirement — read this before reimplementing tiling

Tiles **must** overlap by a fixed pixel amount (100 px was used in
validation) so objects near tile boundaries — especially ramified glia —
aren't truncated by the cut. `analysis/merge_brain.py` relies on this fixed
overlap to deduplicate objects: it computes each tile's "inner zone" (the
non-overlapping core region) and keeps only objects whose global centroid
falls inside their own tile's inner zone.

**If you implement tiling with zero overlap, the dedup step in
`merge_brain.py` is not meaningful and must be skipped/disabled** — with
zero overlap there's no redundant region to dedup against, and the
inner-zone math will incorrectly drop objects near tile edges.
