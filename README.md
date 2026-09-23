# COMET CellProfiler Pipeline

CellProfiler segmentation and phenotyping of 28-channel COMET multiplex
immunofluorescence mouse brain images — microglia, astrocytes, neurons and
CD4/CD8 T cells — run at scale on a Slurm cluster by tiling each whole-slide
image, processing tiles as a job array, and merging per brain with overlap-zone
deduplication. Downstream: single-cell clustering and spatial neighbourhood
analysis.

Private repo, set up for reproducibility rather than distribution.

---

## Status — read this first

**Current pipeline:** [`pipeline/COMET_28ch_seeded.cppipe`](pipeline/COMET_28ch_seeded.cppipe)
**Current results:** [`docs/results_8-28-26.md`](docs/results_8-28-26.md)

The original astrocyte segmentation was **wrong** and has been replaced. It
thresholded GFAP directly, and the objects it produced were GFAP-*depleted*
relative to surrounding tissue (0.81× the tissue mean) — it was segmenting
neuropil, not cells. GFAP fine processes tile essentially the whole neuropil, so
no threshold on GFAP alone can isolate astrocyte cell bodies.

Both glial populations are now segmented by **growing each nucleus into
marker-positive territory**, which guarantees one nucleus per object. Astrocytes
went from 0.81× to ~5× GFAP enrichment. Microglia counts fell from an implausible
20.8% of all cells to ~5.7%.

Consequences:

- [`docs/results_8-26-26.md`](docs/results_8-26-26.md) is **superseded** and
  carries a warning banner. Its headline astrocyte finding does not survive.
- [`pipeline/COMET_track1plus2_28ch.cppipe`](pipeline/COMET_track1plus2_28ch.cppipe)
  is kept for provenance. Use the seeded pipeline for new work.
- Neuron segmentation was **not** affected and is unchanged.

## If you are new to this project

1. **[`docs/results_8-28-26.md`](docs/results_8-28-26.md)** — what was found, what
   was wrong and fixed, and explicitly what does *not* hold up. Start here.
2. **[`docs/RUNBOOK.md`](docs/RUNBOOK.md)** — every command from raw image to
   clustered output, plus a troubleshooting playbook and an incident log from a
   real 533-tile run.
3. **[`docs/channel_map_and_thresholds.md`](docs/channel_map_and_thresholds.md)** —
   why each threshold is what it is.
4. **[`docs/spatial_analysis.md`](docs/spatial_analysis.md)** — the spatial method
   and four pitfalls that produce false results.

The most useful thing to understand early is the **n = 2 problem**. There are two
brains per group, so most treated-vs-control comparisons have within-group spread
exceeding the between-group difference, and cluster-frequency comparisons in
particular are not stable — three clustering configurations gave three different
answers for the same population. The one robust finding is a *within-brain*
comparison, where each brain is its own control. Section 4 of the results doc
lays out which claims survive and which don't.

## What each cell type is for

The four populations are not all being asked the same question, and the analysis
differs accordingly.

| population | n (8 brains) | what we want from it |
|---|---|---|
| **Microglia** | 132,271 | The richest population. Activation state from the marker panel, morphology (size, ramification, skeleton), and **spatial relationships** — which cells they sit near, and how that changes with activation. Most conclusions rest on these. |
| **Astrocytes** | 35,480 | Same intent as microglia — state, morphology, spatial context — but the population is GFAP-defined, so it is a subset (reactive/fibrous astrocytes) rather than all astrocytes. |
| **Neurons** | 460,410 | Mostly a **reference population**: a denominator for composition, and context for the glial neighbourhood analysis. The panel carries little neuronal biology beyond NeuN, so we are not phenotyping neuronal states. |
| **T cells** | 1,162 | **Counting and phenotyping only** — how many are there, and are they CD4 or CD8. Far too few (9–14 CD8 per brain) for morphology or clustering claims. They also appear as a *neighbour* in the microglial spatial analysis, where the question is asked from the T-cell side so that each cell contributes a fraction rather than a count. |

This is why, for example, the morphology work (skeletons, trunk counts, size) is
applied to microglia and astrocytes but not T cells, and why the T-cell section of
the results is about gating and counts rather than states.

## What changed from the original pipeline, and why

| change | reason |
|---|---|
| **Marker panel on `Cells` extended 4 → 15 markers** | `Neurons` and `Tcells` are children of `Cells` and previously carried *no* expression or morphology at all — 12 columns of position and parentage. Without this they could not be phenotyped. No threshold or gate was touched; object counts are bit-identical. |
| **Astrocytes and microglia re-segmented by nucleus seeding** | Thresholding GFAP produced objects that were GFAP-*depleted* (0.81× tissue). Seeding also guarantees one nucleus per object, which removed the merged process networks (~10% of objects, some spanning several cells with 39 primary processes). |
| **T-cell gate made per-brain relative** (`CD3 ≥ 1.5 × brain median`) **and a coreceptor required** (CD4 or CD8 > 2× median) | The absolute `CD3 ≥ 0.022` gate swung counts 78-fold across brains for non-biological reasons — a control brain outranked every treated one, and one brain called 6.26% of all cells T cells. Requiring a coreceptor let the CD3 cut drop while *increasing* specificity. |
| **CD8 ratio cutoff set per cohort** (ICV 1.0, IP 1.5) | The ratio distribution is compressed in the IP batch: ICV controls contain no cell above 0.86, but IP controls reach 1.52. |
| **Clustering can run on a feature subset** (`--cluster-features`) | Morphology features dominated the embedding and created a cluster defined purely by object size — which turned out to be merged multi-cell objects. Clustering on expression while *carrying* morphology for plotting separates the two roles. |
| **`run_array.sh`: `express`/2h → `short`/12h, plus a hung-node watchdog** | One tile exceeded 2h and was killed; runtime tracks object density, not file size, so slow tiles can't be predicted. Separately, nodes wedged four times — the job now times itself out, excludes that node, and requeues. |
| **`submit_all_brains.sh` chains on the array, not the merge** | The merge uses `afterok`, so one failed tile left it permanently pending, and because it never terminated the next brain's `afterany` never fired. This deadlocked the run twice. |

## Prerequisites

- Singularity or Apptainer
- Slurm (scripts target `sbatch` / `--array`)
- conda or mamba
- CellProfiler **4.2.8** specifically — see [`container/README.md`](container/README.md)

## Quickstart

```bash
# 1. container
cd container && singularity build cellprofiler.sif cellprofiler.def
#    or: singularity pull cellprofiler.sif docker://cellprofiler/cellprofiler:4.2.8

# 2. ALWAYS verify the channel map first -- the pipeline selects channels by
#    plane INDEX, so a different batch silently measures the wrong marker
python3 tools/verify_channel_map.py /path/to/*.ome.tiff

# 3. tile (4000 px, 100 px overlap -- the overlap is required by the merge step)
bash tiling/tile_all.sh

# 4. environment
conda env create -f environment.yml && conda activate comet-pipeline
export COMET_PROJECT_ROOT=/path/to/your/project

# 5. ONE tile first, and check the numbers before committing to a full run
sbatch --array=1-1 --export=ALL,IMAGE=your_brain hpc/run_array.sh
#    astrocyte median GFAP should be ~0.10 and microglia Iba1 ~0.04.
#    If astrocyte GFAP comes back near 0.024, the seeding is not working.

# 6. full run
bash hpc/submit_all_brains.sh brain1 brain2 brain3
```

Then locally:

```bash
python3 analysis/make_cafe_csvs.py OUTDIR run_*/merged_* --extra-morphology
python3 analysis/cluster_cells.py OUTDIR Microglia
python3 analysis/spatial_neighborhood.py OUTDIR --radius-um 25 --split-by CD68
```

## Documentation

| Doc | Covers |
|---|---|
| [`docs/results_8-28-26.md`](docs/results_8-28-26.md) | **Current results.** The segmentation fix, the spatial finding, marker-level treated-vs-control, and what does not hold up. |
| [`docs/RUNBOOK.md`](docs/RUNBOOK.md) | End-to-end operations: environments, channel verification, tiling, upload, SSH, Slurm, troubleshooting, download, analysis. Incident log included. |
| [`docs/spatial_analysis.md`](docs/spatial_analysis.md) | Spatial method, and the four pitfalls (density, n=2, radius, circularity). |
| [`docs/channel_map_and_thresholds.md`](docs/channel_map_and_thresholds.md) | Channel map, threshold rationale, CD4/CD8 validation evidence. |
| [`docs/results_8-26-26.md`](docs/results_8-26-26.md) | ⚠️ Superseded. Kept to record what changed and why. |

## Known limitations

- **n = 2 brains per group.** Treated-vs-control results are descriptive. No
  significance testing — per-cell tests would be pseudoreplication and 2-vs-2
  brain-level tests have no power.
- **GFAP under-counts astrocytes.** It labels fibrous and reactive astrocytes
  strongly; protoplasmic gray-matter astrocytes are GFAP-dim and fall below any
  usable threshold. These are "GFAP+ astrocytes", not all astrocytes.
- **Validated on the 28-channel panel only.** Re-verify the channel map from
  OME-XML for any other panel or batch — `tools/verify_channel_map.py` does this,
  but note that **tiles do not retain channel names**, so it must run against an
  original whole-slide file.
- **Acquisition batches differ.** The IP cohort has compressed dynamic range
  relative to ICV (CD3 peak signal differs ~7×). Absolute intensity gates do not
  transfer between batches — this is why T-cell gating is now per-brain relative.
- **CellProfiler pinned to 4.2.8.** Adaptive thresholding behaviour has changed
  across versions.
- **No anatomical annotation.** Spatial coordinates are raw pixels with no region
  labels.
- **Pixel size assumed 0.28 µm/px**, inferred rather than read from OME
  `PhysicalSizeX`. Confirm before reporting distances.

## Where to pick up

Open threads, roughly in order of value:

- **Confirm the pixel size** from OME `PhysicalSizeX`. Cheap, and every stated
  distance depends on it.
- **Continuous rather than median-split** versions of the spatial test, to show
  the T-cell association is graded with CD74/CD68 level.
- **A permutation null** for colocalisation, to put a direct p-value on it.
- **Anatomical ROIs** (QuPath/ImageJ → point-in-polygon on `global_x`/`global_y`).
  This is what would exclude a regional confound in the spatial result.
- **More animals.** Most treated-vs-control questions here are limited by n=2, not
  by the analysis.

## Repo layout

```
pipeline/   COMET_28ch_seeded.cppipe        <- current, nucleus-seeded
            COMET_track1plus2_28ch.cppipe   <- original, kept for provenance
container/  Singularity recipe for CellProfiler 4.2.8
tiling/     tile.sh, tile_all.sh, make_tiles.py + the manifest contract
tools/      verify_channel_map.py           channel names vs plane indices
hpc/        run_array.sh (with hung-node watchdog), submit_all_brains.sh,
            release_held.sh
analysis/   merge_brain.py                  per-brain merge + overlap dedup
            make_cafe_csvs.py               reduced-set per-sample CSVs
            cluster_cells.py                scale/PCA/Harmony/UMAP/Leiden
            spatial_neighborhood.py         neighbourhood composition
docs/       results, runbook, spatial method, channel map
```

## Citation / acknowledgment

TODO — add lab name and grant number before this repo is ever made public.
