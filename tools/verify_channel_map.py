#!/usr/bin/env python3
"""
Verify that a COMET OME-TIFF's embedded channel names match the marker -> T-plane
map hardcoded in the CellProfiler pipeline's NamesAndTypes module.

The pipeline selects channels by plane INDEX, not by name. If an acquisition batch
ordered its cycles differently, the pipeline will silently measure the wrong marker.
Run this on every image before tiling -- especially across cohorts.

Usage:
    python3 verify_channel_map.py IMAGE.ome.tiff [IMAGE2.ome.tiff ...]
    python3 verify_channel_map.py Nick_Raw_Images_8-26-26/*.ome.tiff

Exit code 0 if every image matches, 1 if any mismatch.
"""
import re
import sys
import os
import tifffile

# marker -> 0-based T-plane index, read from COMET_28ch NamesAndTypes rules
EXPECTED = {
    "DAPI": 0, "CD206": 5, "CD4": 6, "NeuN": 7, "CD8": 8, "CD3": 10,
    "GFAP": 14, "APOE": 15, "Iba1": 17, "Arg1": 19, "CD44": 20,
    "CD74": 21, "CD11b": 23, "CD68": 24, "MAA": 26, "SNA": 27,
}

# OME channel names are not always spelled exactly as the pipeline names them.
# Map pipeline name -> regex that should match the OME <Channel Name="...">.
ALIASES = {
    "Arg1": r"arg(inase)?\s*-?1",
    "MAA":  r"\bma[al]\b|mal\s*i",
    "SNA":  r"\bsna\b|snl|ebl",
    "CD11b": r"cd\s*11\s*b",
    "CD206": r"cd\s*206",
    "Iba1": r"iba\s*-?1",
}


def channel_names(path):
    with tifffile.TiffFile(path) as tif:
        xml = tif.ome_metadata or ""
        shape, axes = tif.series[0].shape, tif.series[0].axes
    names = re.findall(r'<Channel[^>]*?Name="([^"]*)"', xml)
    return names, shape, axes


def matches(pipeline_name, ome_name):
    ome = ome_name.strip().lower()
    pat = ALIASES.get(pipeline_name)
    if pat:
        return re.search(pat, ome) is not None
    return re.search(r"\b" + re.escape(pipeline_name.lower()) + r"\b", ome) is not None


def check(path):
    try:
        names, shape, axes = channel_names(path)
    except Exception as e:
        print(f"  ERROR reading: {e}")
        return False

    print(f"  shape={shape} axes={axes} channels_in_ome={len(names)}")
    if not names:
        print("  ERROR: no <Channel Name=...> entries in OME-XML -- cannot verify.")
        print("  NOTE: make_tiles.py does not copy channel names into the tiles it")
        print("        writes, so tile_r*_c*.ome.tiff can never be verified this way.")
        print("        Run this on the ORIGINAL whole-slide OME-TIFF instead.")
        print("        (Tiles do preserve plane ORDER, so verifying the original")
        print("         is sufficient -- the tiles inherit whatever it says.)")
        return False

    ok = True
    for marker, idx in sorted(EXPECTED.items(), key=lambda kv: kv[1]):
        if idx >= len(names):
            print(f"  MISMATCH  plane {idx:2d}  expected {marker:6s}  -- plane does not exist")
            ok = False
            continue
        actual = names[idx]
        if matches(marker, actual):
            print(f"  ok        plane {idx:2d}  {marker:6s} <- {actual!r}")
        else:
            print(f"  MISMATCH  plane {idx:2d}  expected {marker:6s}  found {actual!r}")
            ok = False

    if not ok:
        print("\n  Full OME channel order for this image:")
        for i, n in enumerate(names):
            print(f"    {i:2d}  {n}")
    return ok


def main(paths):
    if not paths:
        print(__doc__)
        return 1
    all_ok = True
    for p in paths:
        print(f"\n=== {os.path.basename(p)} ===")
        if os.path.getsize(p) == 0:
            print("  SKIP: file is 0 bytes (still copying?)")
            all_ok = False
            continue
        all_ok &= check(p)
    print("\n" + "=" * 60)
    print("ALL IMAGES MATCH THE PIPELINE CHANNEL MAP" if all_ok
          else "PROBLEMS FOUND -- do not tile until resolved (see above)")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
