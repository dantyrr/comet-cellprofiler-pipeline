# COMET CellProfiler Pipeline

CellProfiler segmentation and phenotyping of 28-channel COMET multiplex
immunofluorescence mouse brain images — microglia, astrocytes, neurons, and
CD4/CD8 T cells — run at scale on a Slurm HPC cluster by tiling each
whole-slide image and processing tiles as a job array, then merging results
per brain with overlap-zone deduplication.

This repo is private and is set up for personal reproducibility (so future-me
can rebuild the environment and rerun the pipeline), not for external
distribution — see "Known limitations" below for what's missing if that
changes later.

## Prerequisites

- Singularity or Apptainer (container runtime)
- Slurm (job scheduler) — scripts target Slurm's `sbatch`/`--array`
- conda or mamba (for `merge_brain.py`'s Python environment)
- CellProfiler **4.2.8** specifically — see `container/README.md` for why the
  version is pinned

## Quickstart

1. **Build or pull the container:**

   ```bash
   cd container
   singularity build cellprofiler.sif cellprofiler.def
   # or: singularity pull cellprofiler.sif docker://cellprofiler/cellprofiler:4.2.8
   ```

2. **Tile your images:**

   ```bash
   bash tiling/tile.sh SAMPLE_ID
   ```

   `tiling/tile.sh` is real and working — it finds the source whole-slide
   OME-TIFF for a sample ID under Box and calls `tiling/make_tiles.py`.
   **`make_tiles.py` itself is still a stub** (see `tiling/README.md`) —
   paste the real implementation in from `~/Downloads/make_tiles.py`
   before this step will actually run. Whatever implementation you use,
   it must produce output matching the contract documented in
   `tiling/README.md` (tile naming, 100px overlap, `tile_manifest.json`
   schema).

3. **Set up your environment:**

   ```bash
   conda env create -f environment.yml
   conda activate comet-pipeline
   cp .env.example .env   # then edit .env with your actual paths
   export COMET_PROJECT_ROOT=/path/to/your/project
   export SLURM_PARTITION=your_partition_name   # optional, defaults to "express"
   ```

   Your project root should contain the built `cellprofiler.sif`, the
   `<image>_tiles_4k/` directories from tiling, and this repo's
   `pipeline/`, `hpc/`, and `analysis/` folders (or symlinks to them).

4. **Run one test tile** to confirm the container and pipeline work before
   committing to a full array:

   ```bash
   IMAGE=your_brain_name sbatch --export=ALL,IMAGE=$IMAGE --array=1-1 hpc/run_array.sh
   ```

5. **Submit the full run** across one or more brains:

   ```bash
   ./hpc/submit_all_brains.sh brain1 brain2 brain3
   ```

   This chains each brain's segmentation array → merge job → next brain's
   array, capping cluster concurrency and isolating per-brain failures.

## Scientific rationale

See [`docs/channel_map_and_thresholds.md`](docs/channel_map_and_thresholds.md)
for the full channel map, threshold values, and the validation evidence
behind the ratio-based CD4/CD8 T cell subtyping. Do not change any threshold,
filter, or gate value without re-validating against real data — the current
values were empirically validated across multiple test images including
positive/negative controls.

## Known limitations

- **Validated at 28-channel panel only.** Other COMET panels or acquisition
  batches require re-verifying the channel map from the image's OME-XML —
  do not assume the T-plane offsets transfer:

  ```python
  import tifffile
  with tifffile.TiffFile(path) as tif:
      xml = tif.ome_metadata
  # then parse <Channel Name="..."> entries in order
  ```

- **CellProfiler version pinned to 4.2.8.** Other versions are not
  guaranteed to reproduce identical segmentation results (adaptive
  thresholding behavior has changed across CP versions historically).

- **`tiling/make_tiles.py` not included** — only a stub, plus the
  documented output contract (`tiling/README.md`). `tiling/tile.sh` (the
  wrapper that finds source images and calls it) is real and included.
  Any tiling implementation you use must preserve the manifest schema and
  the fixed-pixel tile overlap, or `analysis/merge_brain.py`'s dedup step
  won't work correctly.

## Repo layout

```
pipeline/   the CellProfiler .cppipe pipeline (do not edit thresholds without re-validating)
container/  Singularity build recipe for CellProfiler 4.2.8
tiling/     tile.sh (real) + make_tiles.py (stub) + the output contract
hpc/        Slurm array + multi-brain submission scripts
analysis/   merge_brain.py — per-brain CSV merge with overlap dedup
docs/       full scientific rationale (channel map, thresholds, validation)
```

## Citation / acknowledgment

TODO — add lab name and grant number here before this repo is ever made
public.
