#!/usr/bin/env python3
"""
Merge per-tile CellProfiler outputs into per-brain CSVs with overlap-zone dedup.

Usage:
    python3 merge_brain.py <BRAIN_NAME> [--allow-incomplete]

Example:
    python3 merge_brain.py ICV_C1_3

Refuses to merge (exit 2) if any tile in the manifest produced no output, since
merging then yields a brain that looks complete but is short. --allow-incomplete
overrides. Writes <BRAIN>_tile_coverage.csv recording which tiles contributed.

Project root is read from the COMET_PROJECT_ROOT environment variable
(see .env.example / README). You can also override it by passing it as
a second CLI argument: `python3 merge_brain.py ICV_C1_3 /path/to/project`.

Inputs:
    - <PROJECT>/<BRAIN>_tiles_4k/tile_manifest.json
    - <PROJECT>/results/<BRAIN>/tile_*/MyExpt_*.csv

Outputs:
    - <PROJECT>/results/<BRAIN>/merged/<BRAIN>_<class>.csv  (one per object class, dedup'd)
    - <PROJECT>/results/<BRAIN>/merged/<BRAIN>_summary.csv  (pre/post dedup counts)

Dedup logic:
    For each tile, compute the "inner zone" — the region of its 4000x4000 area
    that doesn't overlap with neighbors. With 100-px overlap, the inner zone
    excludes 50 px on each side that has a neighbor. Edge tiles (col=0, col=max,
    row=0, row=max) extend their inner zone to the global image boundary on the
    edgeward side. Each global pixel is then owned by exactly one tile.

    An object is kept if its global centroid (tile.x0 + Location_Center_X,
    tile.y0 + Location_Center_Y) falls within ITS tile's inner zone.
    Otherwise it's a duplicate that another tile owns instead.

    NOTE: this dedup logic assumes tiles were generated with a fixed pixel
    overlap (see tiling/README.md). If a tiling implementation uses zero
    overlap, this dedup step is not meaningful and should be skipped.
"""

import json
import sys
import glob
import os
import pandas as pd
from pathlib import Path

if len(sys.argv) not in (2, 3):
    print(f"Usage: {sys.argv[0]} <BRAIN_NAME> [PROJECT_ROOT]", file=sys.stderr)
    print(f"Example: {sys.argv[0]} ICV_C1_3", file=sys.stderr)
    sys.exit(1)

ALLOW_INCOMPLETE = "--allow-incomplete" in sys.argv
if ALLOW_INCOMPLETE:
    sys.argv.remove("--allow-incomplete")

brain = sys.argv[1]

# ----- resolve project root: CLI arg > env var > current directory -----
if len(sys.argv) == 3:
    PROJECT = Path(sys.argv[2])
else:
    PROJECT = Path(os.environ.get("COMET_PROJECT_ROOT", "."))

if not PROJECT.exists():
    print(
        f"ERROR: project root does not exist: {PROJECT}\n"
        f"Set COMET_PROJECT_ROOT env var, or pass the path as a second argument.",
        file=sys.stderr,
    )
    sys.exit(1)

tiles_dir   = PROJECT / f"{brain}_tiles_4k"
results_dir = PROJECT / "results" / brain
merged_dir  = results_dir / "merged"
manifest_fp = tiles_dir / "tile_manifest.json"

if not manifest_fp.exists():
    print(f"ERROR: manifest not found at {manifest_fp}", file=sys.stderr)
    sys.exit(1)
if not results_dir.exists():
    print(f"ERROR: results dir not found at {results_dir}", file=sys.stderr)
    sys.exit(1)

merged_dir.mkdir(parents=True, exist_ok=True)

# ----- read manifest, compute inner zones -----
manifest = json.loads(manifest_fp.read_text())
overlap = manifest["overlap"]
full_H  = manifest["full_H"]
full_W  = manifest["full_W"]
tiles   = manifest["tiles"]
half_overlap = overlap // 2

max_row = max(t["row"] for t in tiles)
max_col = max(t["col"] for t in tiles)

# index tiles by name for fast lookup, and add inner-zone bounds
tile_meta = {}
for t in tiles:
    name = t["file"].replace(".ome.tiff", "")
    inner_x_min = t["x0"] if t["col"] == 0       else t["x0"] + half_overlap
    inner_x_max = t["x1"] if t["col"] == max_col else t["x1"] - half_overlap
    inner_y_min = t["y0"] if t["row"] == 0       else t["y0"] + half_overlap
    inner_y_max = t["y1"] if t["row"] == max_row else t["y1"] - half_overlap
    tile_meta[name] = {
        "row": t["row"], "col": t["col"],
        "x0": t["x0"], "y0": t["y0"],
        "inner_x_min": inner_x_min, "inner_x_max": inner_x_max,
        "inner_y_min": inner_y_min, "inner_y_max": inner_y_max,
    }

print(f"Brain: {brain}")
print(f"  Project root: {PROJECT}")
print(f"  Manifest tiles: {len(tiles)}")
print(f"  Full image: {full_W} x {full_H} px")
print(f"  Tile overlap: {overlap} px ({half_overlap} on each side of midline)")
print(f"  Output: {merged_dir}")

# ----- completeness check: every manifest tile must have produced output -----
#
# CellProfiler writes MyExpt_Image.csv on every successful run regardless of how
# many objects were found, so its absence means that tile did not complete --
# as opposed to a per-class CSV being absent, which legitimately happens when a
# tile contains none of that object (e.g. no CD8 T cells).
#
# Without this gate the merge silently skips failed tiles and produces a
# brain that looks normal but is short. That happened in the 8-28-26 run:
# one tile was killed after a node hung, and the merged brain was ~2.3% short
# with nothing in the output to indicate it.
missing = [n for n in tile_meta if not (results_dir / n / "MyExpt_Image.csv").exists()]
if missing:
    print(f"\n*** INCOMPLETE: {len(missing)} of {len(tile_meta)} tiles produced no output ***",
          file=sys.stderr)
    for n in sorted(missing):
        print(f"      {n}", file=sys.stderr)
    print("\n  These tiles failed or were never run. Merging now would produce a",
          file=sys.stderr)
    print("  brain that looks complete but is short by those tiles' cells.",
          file=sys.stderr)
    print("  Re-run them, then merge again. Array indices are the tiles' position",
          file=sys.stderr)
    print("  in `ls <BRAIN>_tiles_4k/tile_r*_c*.ome.tiff | sort`.", file=sys.stderr)
    if not ALLOW_INCOMPLETE:
        print("\n  Refusing to merge. Pass --allow-incomplete to override.\n",
              file=sys.stderr)
        sys.exit(2)
    print("\n  --allow-incomplete given; merging anyway.\n", file=sys.stderr)
else:
    print(f"  Completeness: all {len(tile_meta)} tiles produced output")

# ----- discover which object classes exist -----
sample_tile = next(n for n in tile_meta
                   if (results_dir / n / "MyExpt_Image.csv").exists())
sample_csvs = sorted(glob.glob(str(results_dir / sample_tile / "MyExpt_*.csv")))
if not sample_csvs:
    print(f"ERROR: no MyExpt_*.csv in {results_dir / sample_tile}", file=sys.stderr)
    sys.exit(1)
classes = [os.path.basename(c).replace("MyExpt_", "").replace(".csv", "") for c in sample_csvs]
print(f"  Object classes: {classes}\n")

# ----- merge each class -----
summary_rows = []
for cls in classes:
    print(f"=== {cls} ===")

    # special-case Image.csv and Experiment.csv -- one row per tile, no dedup
    is_perimage = cls in ("Image", "Experiment")

    per_tile_dfs = []
    n_per_tile = {}
    for name, meta in tile_meta.items():
        fp = results_dir / name / f"MyExpt_{cls}.csv"
        if not fp.exists():
            continue
        try:
            df = pd.read_csv(fp)
        except pd.errors.EmptyDataError:
            n_per_tile[name] = 0
            continue
        if len(df) == 0:
            n_per_tile[name] = 0
            continue
        df.insert(0, "tile_name", name)
        df.insert(1, "tile_row",  meta["row"])
        df.insert(2, "tile_col",  meta["col"])
        if not is_perimage:
            # add global coords if Location_Center_X/Y exist
            if "Location_Center_X" in df.columns and "Location_Center_Y" in df.columns:
                df["global_x"] = meta["x0"] + df["Location_Center_X"]
                df["global_y"] = meta["y0"] + df["Location_Center_Y"]
        n_per_tile[name] = len(df)
        per_tile_dfs.append(df)

    if not per_tile_dfs:
        print(f"  No data found across any tile. Skipping.")
        summary_rows.append({"class": cls, "n_raw": 0, "n_dedup": 0, "n_dropped": 0})
        continue

    merged = pd.concat(per_tile_dfs, ignore_index=True)
    n_raw = len(merged)

    if is_perimage:
        # no dedup, just write
        out_fp = merged_dir / f"{brain}_{cls}.csv"
        merged.to_csv(out_fp, index=False)
        print(f"  rows: {n_raw} (no dedup for {cls})  -> {out_fp.name}")
        summary_rows.append({"class": cls, "n_raw": n_raw, "n_dedup": n_raw, "n_dropped": 0})
        continue

    # ---- dedup by inner-zone ownership ----
    if "global_x" not in merged.columns:
        # class without Location_Center (e.g., empty CSVs) — write unmodified
        out_fp = merged_dir / f"{brain}_{cls}.csv"
        merged.to_csv(out_fp, index=False)
        print(f"  rows: {n_raw} (no Location_Center cols, kept all)  -> {out_fp.name}")
        summary_rows.append({"class": cls, "n_raw": n_raw, "n_dedup": n_raw, "n_dropped": 0})
        continue

    # vectorized inner-zone check per row using its own tile's bounds
    keep = (
        (merged["global_x"] >= merged["tile_name"].map(lambda n: tile_meta[n]["inner_x_min"])) &
        (merged["global_x"] <  merged["tile_name"].map(lambda n: tile_meta[n]["inner_x_max"])) &
        (merged["global_y"] >= merged["tile_name"].map(lambda n: tile_meta[n]["inner_y_min"])) &
        (merged["global_y"] <  merged["tile_name"].map(lambda n: tile_meta[n]["inner_y_max"]))
    )
    dedup = merged[keep].reset_index(drop=True)
    n_dedup = len(dedup)
    n_dropped = n_raw - n_dedup

    out_fp = merged_dir / f"{brain}_{cls}.csv"
    dedup.to_csv(out_fp, index=False)
    print(f"  raw: {n_raw}  dedup: {n_dedup}  dropped: {n_dropped} ({100*n_dropped/n_raw if n_raw else 0:.1f}%)  -> {out_fp.name}")
    summary_rows.append({"class": cls, "n_raw": n_raw, "n_dedup": n_dedup, "n_dropped": n_dropped})

# ----- summary -----
# record which tiles contributed to each class, so a short brain is traceable
coverage = pd.DataFrame(
    [{"tile": n,
      "has_output": (results_dir / n / "MyExpt_Image.csv").exists()}
     for n in sorted(tile_meta)])
coverage.to_csv(merged_dir / f"{brain}_tile_coverage.csv", index=False)

summary_df = pd.DataFrame(summary_rows)
summary_fp = merged_dir / f"{brain}_summary.csv"
summary_df.to_csv(summary_fp, index=False)
print(f"\nSummary written to: {summary_fp}")
print(summary_df.to_string(index=False))
