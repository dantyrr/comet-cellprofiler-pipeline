# Tiling

`tile.sh` is the real, working wrapper script — it runs locally (on
dtyrrell's Mac), finds the source OME-TIFF for a given sample ID under the
lab's Box folder, and calls `make_tiles.py` to produce tiles + manifest.

`make_tiles.py` itself is still a **stub**. The actual tiling
implementation lives at `~/Downloads/make_tiles.py` on the researcher's
local machine and wasn't included when this repo was packaged — paste it
into `tiling/make_tiles.py` next time you're back in this repo. Until
then, `tile.sh` will fail with a clear error pointing here.

Whatever `make_tiles.py` implementation you use, it **must** satisfy the
output contract below (call signature: `make_tiles.py <input.ome.tiff>
<output_dir> --tile <size> --overlap <overlap>`), since `hpc/run_array.sh`
and `analysis/merge_brain.py` depend on it exactly.

## Running tile.sh

```bash
bash tiling/tile.sh SAMPLE_ID [tile_size] [overlap]
# e.g.
bash tiling/tile.sh ICV-T2_3
bash tiling/tile.sh IP-C1_3 4000 100
```

Defaults (all overridable via env var, see `.env.example`):

| Var | Default | Purpose |
|---|---|---|
| `COMET_BOX_DATA_ROOT` | dtyrrell's Box "Nick's Comet data" folder | where source whole-slide OME-TIFFs live |
| `MAKE_TILES_SCRIPT` | `tiling/make_tiles.py` (next to `tile.sh`) | tiling implementation to call |
| `COMET_TILES_OUTPUT_ROOT` | `$HOME/Desktop` | where `<SAMPLE_ID>_tiles_4k/` gets written |

After tiling, the output folder needs to be uploaded to the HPC cluster
(Cheaha, via the OnDemand Files app) into `$COMET_PROJECT_ROOT` before
running `hpc/run_array.sh`.

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
