#!/usr/bin/env python3
"""
Microglia-centred spatial neighbourhood analysis.

For every microglion, count the cells of each type within a radius, then compare
neighbourhood composition between treated and control within each delivery route.

Microglia-centred because they are the numerous population (13k-22k per brain);
T-cell-centred versions of the same question rest on 9-14 CD8 cells per brain and
are not powered.

Usage:
    python3 spatial_neighborhood.py DATA_DIR [--radius-um 50] [--px-um 0.28]
                                             [--split-by CD68] [--out FILE]

Coordinates are `global_x`/`global_y` -- whole-slide pixels, valid across tile
boundaries because merge_brain.py deduplicates on the global centroid.
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

TYPES = ["Microglia", "Astrocytes", "Neurons", "Tcells"]


def load(data_dir, ct, sample):
    f = os.path.join(data_dir, ct, sample + ".csv")
    return pd.read_csv(f) if os.path.exists(f) else None


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("data_dir")
    p.add_argument("--radius-um", type=float, default=50.0)
    p.add_argument("--px-um", type=float, default=0.28,
                   help="microns per pixel; confirm from OME PhysicalSizeX")
    p.add_argument("--split-by", default=None,
                   help="microglial marker to split high/low on (median within brain)")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    R = a.radius_um / a.px_um
    print(f"radius {a.radius_um} um = {R:.0f} px  (at {a.px_um} um/px)\n")

    man = pd.read_csv(os.path.join(a.data_dir, "samples.csv"))
    samples = man[man.CellType == "Microglia"][["File", "BrainID", "Group"]]
    samples["sample"] = samples.File.str.replace(".csv", "", regex=False)

    rows, per_cell = [], []
    for r in samples.itertuples():
        mg = load(a.data_dir, "Microglia", r.sample)
        if mg is None or len(mg) == 0:
            continue
        pts = mg[["global_x", "global_y"]].values
        rec = {"Group": r.Group, "BrainID": r.BrainID, "n_microglia": len(mg)}
        counts = {}
        for ct in TYPES:
            other = load(a.data_dir, ct, r.sample)
            if other is None or len(other) == 0:
                counts[ct] = np.zeros(len(mg)); continue
            tree = cKDTree(other[["global_x", "global_y"]].values)
            n = np.array([len(x) for x in tree.query_ball_point(pts, R)], float)
            if ct == "Microglia":
                n -= 1                      # don't count itself
            counts[ct] = n
            rec[f"mean_{ct}"] = n.mean()
        # T cells split by subtype
        tc = load(a.data_dir, "Tcells", r.sample)
        if tc is not None and len(tc):
            for sub in ("CD4", "CD8"):
                s = tc[tc.Subtype == sub]
                if len(s):
                    tree = cKDTree(s[["global_x", "global_y"]].values)
                    rec[f"mean_{sub}"] = np.mean([len(x) for x in tree.query_ball_point(pts, R)])
                else:
                    rec[f"mean_{sub}"] = 0.0
        rows.append(rec)

        d = pd.DataFrame({k: v for k, v in counts.items()})
        d["Group"], d["BrainID"] = r.Group, r.BrainID
        if a.split_by and a.split_by in mg.columns:
            d[a.split_by] = mg[a.split_by].values
        per_cell.append(d)

    t = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    cols = ["Group", "BrainID", "n_microglia"] + [c for c in t.columns if c.startswith("mean_")]
    print("Mean neighbours per microglion, by brain")
    print(t[cols].round(3).to_string(index=False))

    print(f"\nGroup means, and treated/control ratio within route "
          f"(radius {a.radius_um} um)")
    g = t.groupby("Group")[[c for c in t.columns if c.startswith("mean_")]].mean()
    fc = pd.DataFrame({"ICV T/C": g.loc["ICV-T"] / g.loc["ICV-C"],
                       "IP T/C": g.loc["IP-T"] / g.loc["IP-C"]}).round(3)
    fc["same direction"] = np.where(np.sign(fc["ICV T/C"] - 1) == np.sign(fc["IP T/C"] - 1),
                                    np.where(fc["ICV T/C"] > 1, "UP both", "DOWN both"), "")
    print(g.round(3).to_string())
    print(fc.to_string())

    pc = pd.concat(per_cell, ignore_index=True)
    if a.split_by:
        print(f"\nNeighbourhood by microglial {a.split_by} (high/low split at each brain's median)")
        pc["hi"] = pc.groupby("BrainID")[a.split_by].transform(lambda s: s > s.median())
        sub = pc.groupby(["Group", "hi"])[TYPES].mean().round(3)
        sub.index = sub.index.set_levels([f"{a.split_by}-low", f"{a.split_by}-high"], level=1)
        print(sub.to_string())

    if a.out:
        pc.to_csv(a.out, index=False)
        print(f"\nper-microglion counts -> {a.out}  ({len(pc):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
