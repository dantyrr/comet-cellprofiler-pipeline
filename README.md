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

   Finds the source whole-slide OME-TIFF for a sample ID under Box and
   calls `tiling/make_tiles.py`, which cuts it into 4000px tiles (100px
   overlap) + `tile_manifest.json`. Both scripts are real/working — see
   `tiling/README.md` for env var overrides and the output contract.

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

## Documentation

| Doc | What it covers |
|---|---|
| [`docs/RUNBOOK.md`](docs/RUNBOOK.md) | **End-to-end operational guide** — every command from raw OME-TIFF to clustered output: environments, channel-map verification, tiling, upload, SSH, Slurm submission, a full troubleshooting playbook (held jobs, `DependencyNeverSatisfied`, hung nodes, timeouts), download, and local analysis. Includes an incident log from the 8-brain run. |
| [`docs/results_8-26-26.md`](docs/results_8-26-26.md) | Results of the 8-brain ICV + IP analysis, and the T-cell gating changes with their evidence and validation. |
| [`docs/channel_map_and_thresholds.md`](docs/channel_map_and_thresholds.md) | Channel map and per-object threshold rationale. |

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

- **Tiling (`tiling/tile.sh` + `tiling/make_tiles.py`) is included and
  working**, tuned to dtyrrell's Mac/Box setup by default. If you ever
  swap in a different tiling implementation, it must preserve the
  manifest schema and the fixed-pixel tile overlap (see
  `tiling/README.md`), or `analysis/merge_brain.py`'s dedup step won't
  work correctly.

## Repo layout

```
pipeline/   the CellProfiler .cppipe pipeline (do not edit thresholds without re-validating)
container/  Singularity build recipe for CellProfiler 4.2.8
tiling/     tile.sh + make_tiles.py (both real) + the output contract
hpc/        Slurm array + multi-brain submission scripts
analysis/   merge_brain.py — per-brain CSV merge with overlap dedup
docs/       runbook, results writeup, channel map + threshold rationale
tools/      verify_channel_map.py -- check OME channel names vs pipeline plane indices
analysis/   merge_brain.py, make_cafe_csvs.py (reduced-set export), cluster_cells.py
```

## Citation / acknowledgment

TODO — add lab name and grant number here before this repo is ever made
public.
