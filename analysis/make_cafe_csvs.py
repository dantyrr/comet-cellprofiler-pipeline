#!/usr/bin/env python3
"""
Build reduced-set per-sample analysis CSVs from merged per-brain CellProfiler output.

Reproduces the "*_both_reduced_set" column selection used in the 7-2-26 ICV analysis
and extends it to Neurons and T cells.

Columns written per cell type
-----------------------------
Microglia, Astrocytes   (measurements live directly on the object)
    AreaShape_Area, AreaShape_MeanRadius, AreaShape_MajorAxisLength,
    AreaShape_Eccentricity, AreaShape_FormFactor
    ObjectSkeleton_NumberTrunks      <- joined in from <BRAIN>_<Class>Seeds.csv
    <marker> for every Intensity_MeanIntensity_<marker> on the object

Neurons, Tcells         (no measurements of their own -- they are children of Cells)
    the same 5 AreaShape columns and all markers, joined in from <BRAIN>_Cells.csv
    via Parent_Cells -> Cells.ObjectNumber, matched within a tile.
    Tcells additionally get: Subtype (CD4/CD8/unassigned), CD8_CD4_ratio.
    No skeleton columns -- these are compact objects with no meaningful skeleton.

Also written: global_x, global_y (whole-slide coordinates, for spatial analysis)
and cell_id.

Usage
-----
    python3 make_cafe_csvs.py OUTPUT_DIR MERGED_DIR [MERGED_DIR ...]

    python3 make_cafe_csvs.py "csv files for analysis_8-26-26" merged_*

Brain ID, route (ICV/IP) and condition (C/T) are parsed from the merged directory
name, and written to OUTPUT_DIR/samples.csv so downstream scripts never have to
parse filenames.
"""
import re
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd

REDUCED_SHAPE = [
    "AreaShape_Area",
    "AreaShape_MeanRadius",
    "AreaShape_MajorAxisLength",
    "AreaShape_Eccentricity",
    "AreaShape_FormFactor",
]

# cell type -> (object csv suffix, seeds csv suffix or None, parent object or None)
CELL_TYPES = {
    "Microglia":  ("Microglia",  "MicrogliaSeeds",  None),
    "Astrocytes": ("Astrocytes", "AstrocyteSeeds",  None),
    "Neurons":    ("Neurons",    None,              "Cells"),
    "Tcells":     ("Tcells",     None,              "Cells"),
}

JOIN_KEY = "tile_name"


def parse_brain(dirname):
    """merged_ICV_C1_3 / ICV_C1_3  ->  (brain_id, route, condition)."""
    stem = os.path.basename(str(dirname).rstrip("/"))
    stem = re.sub(r"^merged_", "", stem)
    m = re.match(r"(ICV|IP)[_-]([CT])(\d+)_(\d+)$", stem)
    if not m:
        raise ValueError(
            f"cannot parse route/condition from directory name {stem!r}; "
            "expected e.g. merged_ICV_C1_3 or merged_IP-T3_3"
        )
    route, cond = m.group(1), m.group(2)
    return stem, route, cond


def marker_columns(df):
    """Intensity_MeanIntensity_Iba1 -> Iba1, in a stable order."""
    out = {}
    for c in df.columns:
        if c.startswith("Intensity_MeanIntensity_"):
            out[c] = c[len("Intensity_MeanIntensity_"):]
    return dict(sorted(out.items(), key=lambda kv: kv[1]))


def trunks_column(seeds):
    """ObjectSkeleton_NumberTrunks_MicrogliaSkeleton -> the raw column name."""
    for c in seeds.columns:
        if c.startswith("ObjectSkeleton_NumberTrunks"):
            return c
    return None


# Default multiplier for the per-brain CD3 threshold (see derive_tcells).
TCELL_CD3_K = 2.0

# CD8/CD4 ratio above which a T cell is called CD8.
#
# The pipeline uses 1.0, validated on ICV data where the separation is clean:
# ICV control brains contain no cell above 0.86, and ICV treated CD8 cells start
# at 1.48. In the IP batch the ratio distribution is compressed -- IP control
# brains reach 1.00-1.52 and treated cells pile up just above 1.0 -- so 1.0 no
# longer separates signal from noise there. 2.0 sits above the IP noise band and
# well below the treated signal, and takes control CD8 counts to 0,0,0,1.
TCELL_CD8_RATIO = 2.0


def derive_tcells(parent, k):
    """
    Re-derive T cells from the Cells object using a PER-BRAIN CD3 threshold.

    The pipeline gates T cells on an absolute CD3 MeanIntensity >= 0.022. CD3
    staining intensity shifts between acquisition batches, so that fixed cutoff
    lands at very different points in each brain's CD3 distribution -- at or above
    the 99th percentile in some brains, well below it in others. Across these 8
    brains it swings T-cell counts ~78x (0.08% to 6.26% of cells) with no
    biological cause, and puts a control brain above every treated one.

    Only the CD3 criterion is replaced, with CD3 >= k * (that brain's median CD3).
    This is the same normalisation the pipeline already applies to CD8/CD4
    subtyping, for the same reason. NeuN, Area and FormFactor are unchanged from
    the validated pipeline values.
    """
    cd3 = parent["Intensity_MeanIntensity_CD3"]
    med = float(cd3.median())
    keep = (
        (cd3 >= k * med)
        & (parent["Intensity_MeanIntensity_NeuN"] <= 0.050)
        & parent["AreaShape_Area"].between(100, 950)
        & (parent["AreaShape_FormFactor"] >= 0.70)
    )
    return parent[keep].copy(), med


def build_tcells_derived(parent, k, cd8_ratio=TCELL_CD8_RATIO):
    """Reduced-set columns for T cells re-derived from Cells."""
    sub, med = derive_tcells(parent, k)
    markers = marker_columns(parent)
    keep = [c for c in REDUCED_SHAPE if c in parent.columns]
    missing = [c for c in REDUCED_SHAPE if c not in parent.columns]

    cols = [JOIN_KEY, "ObjectNumber"] + keep + list(markers)
    if "Math_CD8_CD4_ratio" in parent.columns:
        cols.append("Math_CD8_CD4_ratio")
    if "global_x" in parent.columns:
        cols += ["global_x", "global_y"]

    out = sub[cols].rename(columns={**markers,
                                    "Math_CD8_CD4_ratio": "CD8_CD4_ratio"})
    # CD8 vs CD4 by the pipeline's own ratio rule
    if "CD8_CD4_ratio" in out.columns:
        r = pd.to_numeric(out["CD8_CD4_ratio"], errors="coerce")
        out["Subtype"] = np.where(r >= cd8_ratio, "CD8", "CD4")
        out.loc[~np.isfinite(r), "Subtype"] = "unassigned"
    return out, missing, med


def read(merged_dir, brain, suffix):
    fp = Path(merged_dir) / f"{brain}_{suffix}.csv"
    if not fp.exists():
        return None
    try:
        df = pd.read_csv(fp)
    except pd.errors.EmptyDataError:
        return None
    return df if len(df) else None


def build_direct(obj, seeds):
    """Microglia / Astrocytes: measurements are already on the object."""
    markers = marker_columns(obj)
    keep = [c for c in REDUCED_SHAPE if c in obj.columns]
    missing = [c for c in REDUCED_SHAPE if c not in obj.columns]

    out = obj[[JOIN_KEY, "ObjectNumber"] + keep + list(markers)].copy()
    out = out.rename(columns=markers)

    if "global_x" in obj.columns:
        out["global_x"] = obj["global_x"].values
        out["global_y"] = obj["global_y"].values

    n_before = len(out)
    if seeds is not None:
        tcol = trunks_column(seeds)
        if tcol:
            s = seeds[[JOIN_KEY, "ObjectNumber", tcol]].rename(
                columns={tcol: "ObjectSkeleton_NumberTrunks"}
            )
            out = out.merge(s, on=[JOIN_KEY, "ObjectNumber"], how="inner")
        else:
            missing.append("ObjectSkeleton_NumberTrunks")
    else:
        missing.append("ObjectSkeleton_NumberTrunks")

    dropped = n_before - len(out)
    return out, missing, dropped


def build_child(child, parent):
    """Neurons / Tcells: pull shape + markers across from the parent Cells object."""
    if "Parent_Cells" not in child.columns:
        raise ValueError("child object has no Parent_Cells column")

    markers = marker_columns(parent)
    keep = [c for c in REDUCED_SHAPE if c in parent.columns]
    missing = [c for c in REDUCED_SHAPE if c not in parent.columns]

    pcols = [JOIN_KEY, "ObjectNumber"] + keep + list(markers)
    if "Math_CD8_CD4_ratio" in parent.columns:
        pcols.append("Math_CD8_CD4_ratio")
    p = parent[pcols].rename(columns={**markers,
                                      "ObjectNumber": "Parent_Cells",
                                      "Math_CD8_CD4_ratio": "CD8_CD4_ratio"})

    base = [JOIN_KEY, "ObjectNumber", "Parent_Cells"]
    if "global_x" in child.columns:
        base += ["global_x", "global_y"]

    n_before = len(child)
    out = child[base].merge(p, on=[JOIN_KEY, "Parent_Cells"], how="inner")
    return out, missing, n_before - len(out)


def tcell_subtype(out, merged_dir, brain):
    """Label each T cell CD4 / CD8 / unassigned from the subtype object CSVs."""
    out["Subtype"] = "unassigned"
    for label in ("CD4", "CD8"):
        sub = read(merged_dir, brain, f"{label}_Tcells")
        if sub is None or "Parent_Cells" not in sub.columns:
            continue
        keys = set(zip(sub[JOIN_KEY], sub["Parent_Cells"]))
        hit = [k in keys for k in zip(out[JOIN_KEY], out["Parent_Cells"])]
        out.loc[hit, "Subtype"] = label
    return out


def main(argv):
    argv = list(argv)
    tcell_k = TCELL_CD3_K
    if "--tcell-abs" in argv:                 # use the pipeline's absolute CD3 gate
        argv.remove("--tcell-abs")
        tcell_k = None
    cd8_ratio = TCELL_CD8_RATIO
    if "--cd8-ratio" in argv:
        i = argv.index("--cd8-ratio")
        cd8_ratio = float(argv[i + 1])
        del argv[i:i + 2]
    if "--tcell-k" in argv:
        i = argv.index("--tcell-k")
        tcell_k = float(argv[i + 1])
        del argv[i:i + 2]

    if len(argv) < 3:
        print(__doc__)
        return 1

    outroot = Path(argv[1])
    merged_dirs = argv[2:]
    print("T cell gate: " + ("pipeline absolute CD3 >= 0.022"
                             if not tcell_k else
                             f"per-brain CD3 >= {tcell_k} x brain median"))
    print(f"CD8 call: CD8_CD4_ratio >= {cd8_ratio}")
    outroot.mkdir(parents=True, exist_ok=True)

    manifest = []
    for md in merged_dirs:
        if not os.path.isdir(md):
            print(f"SKIP {md}: not a directory")
            continue
        brain, route, cond = parse_brain(md)
        group = f"{route}-{cond}"
        print(f"\n=== {brain}   route={route} condition={cond} group={group} ===")

        for cname, (suffix, seeds_suffix, parent_name) in CELL_TYPES.items():
            gate_note = ""

            if cname == "Tcells" and tcell_k:
                # re-derive from Cells with a per-brain CD3 threshold
                parent = read(md, brain, "Cells")
                if parent is None:
                    print(f"  {cname:11s} Cells missing -- skipped")
                    continue
                out, missing, med = build_tcells_derived(parent, tcell_k, cd8_ratio)
                dropped = 0
                gate_note = f"  [CD3 >= {tcell_k}x{med:.4f} = {tcell_k*med:.4f}]"
            else:
                obj = read(md, brain, suffix)
                if obj is None:
                    print(f"  {cname:11s} no data ({brain}_{suffix}.csv missing or empty) -- skipped")
                    continue

                if parent_name is None:
                    seeds = read(md, brain, seeds_suffix) if seeds_suffix else None
                    out, missing, dropped = build_direct(obj, seeds)
                else:
                    parent = read(md, brain, parent_name)
                    if parent is None:
                        print(f"  {cname:11s} parent {parent_name} missing -- skipped")
                        continue
                    out, missing, dropped = build_child(obj, parent)
                    if cname == "Tcells":
                        out = tcell_subtype(out, md, brain)

            out = out.drop(columns=[JOIN_KEY, "ObjectNumber", "Parent_Cells"],
                           errors="ignore")
            out.insert(0, "cell_id", [f"{brain}_{cname}_{i}" for i in range(len(out))])

            cdir = outroot / cname
            cdir.mkdir(exist_ok=True)
            fp = cdir / f"{group}_{brain}.csv"
            out.to_csv(fp, index=False)

            note = ""
            if dropped:
                note += f"  ({dropped} dropped on join)"
            if missing:
                note += f"  MISSING: {', '.join(missing)}"
            print(f"  {cname:11s} {len(out):7,d} cells x {len(out.columns)} cols -> {fp.name}{note}{gate_note}")

            manifest.append({
                "CellType": cname, "File": fp.name, "BrainID": brain,
                "Route": route, "Condition": cond, "Group": group,
                "nCells": len(out),
            })

    if not manifest:
        print("\nNothing written.")
        return 1

    mf = pd.DataFrame(manifest)
    mfp = outroot / "samples.csv"
    mf.to_csv(mfp, index=False)
    print(f"\nWrote manifest: {mfp}")
    print(mf.pivot_table(index="BrainID", columns="CellType", values="nCells",
                         aggfunc="sum", fill_value=0).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
