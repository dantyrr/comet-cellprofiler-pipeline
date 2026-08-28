# RUNBOOK — end-to-end, as actually run

Complete operational record of the 8-brain COMET run (ICV + IP cohorts,
2026-08-26/27): every command, in order, from raw whole-slide images through
to clustered single-cell output — plus every failure hit along the way and how
it was recovered.

Results and biological interpretation: [`results_8-26-26.md`](results_8-26-26.md).
Scientific rationale for thresholds: [`channel_map_and_thresholds.md`](channel_map_and_thresholds.md).

Paths below are the ones actually used (`dtyrrell` on UAB Cheaha). Substitute
your own; the scripts themselves read `COMET_PROJECT_ROOT` and friends.

---

## Contents

1. [Environments](#1-environments)
2. [Stage 1 — Local: verify channel map, then tile](#2-stage-1--local-verify-channel-map-then-tile)
3. [Stage 2 — Upload to Cheaha](#3-stage-2--upload-to-cheaha)
4. [Stage 3 — Connect and pre-flight](#4-stage-3--connect-and-pre-flight)
5. [Stage 4 — Smoke test before committing hours](#5-stage-4--smoke-test-before-committing-hours)
6. [Stage 5 — Launch the full run](#6-stage-5--launch-the-full-run)
7. [Stage 6 — Monitoring and the troubleshooting playbook](#7-stage-6--monitoring-and-the-troubleshooting-playbook)
8. [Stage 7 — Verify completion](#8-stage-7--verify-completion)
9. [Stage 8 — Download](#9-stage-8--download)
10. [Stage 9 — Local analysis](#10-stage-9--local-analysis)
11. [Appendix A — Incident log from the 8-26-26 run](#appendix-a--incident-log-from-the-8-26-26-run)
12. [Appendix B — One-liner reference card](#appendix-b--one-liner-reference-card)

---

## 1. Environments

Two conda environments. They are kept separate because the tiling step needs
only `tifffile`/`numpy`, while the analysis step pulls in the whole scverse
stack.

| Purpose | Env | Key versions used |
|---|---|---|
| Tiling, channel verification | `ome-tiff` | tifffile 2024.12.12, numpy 2.2.6 |
| Clustering | `scanpy` | scanpy 1.11.4, anndata 0.12.2, pandas 2.3.3, numpy 2.2.6, python-igraph 0.11.9, harmonypy, seaborn 0.13.2 |
| `make_cafe_csvs.py` | either | pandas + numpy only |

```bash
conda env create -f environment.yml     # or, minimally:
conda create -n ome-tiff python=3.11 tifffile numpy
conda create -n scanpy   python=3.11 scanpy anndata harmonypy python-igraph leidenalg seaborn scikit-learn
```

Resolve the interpreter paths once and reuse them:

```bash
export COMET_PYTHON_TILE=/Users/dtyrrell/miniconda3/envs/ome-tiff/bin/python
export COMET_PYTHON_SC=/Users/dtyrrell/miniconda3/envs/scanpy/bin/python
```

> **seaborn ≥ 0.13 breaking change.** `barplot(..., errwidth=)` was removed in
> 0.13. `cluster_cells.py` uses `err_kws={"linewidth": ...}`. Older analysis
> scripts written against seaborn 0.12 will crash here.

On Cheaha, only Singularity and the stock `python3` (with pandas) are needed —
`merge_brain.py` imports nothing else.

---

## 2. Stage 1 — Local: verify channel map, then tile

### 2.1 Verify the channel map — do this before anything else

The pipeline selects channels by **plane index**, not by name. A different
acquisition batch that ordered its cycles differently will be silently measured
as the wrong marker, and every downstream number will look plausible and be wrong.

```bash
$COMET_PYTHON_TILE tools/verify_channel_map.py /path/to/images/*.ome.tiff
```

Exits non-zero and prints the full OME channel order on any mismatch.

> **Tiles cannot be verified this way.** `make_tiles.py` writes tiles with
> `metadata={"axes": "TYX"}` and does **not** copy `<Channel Name=...>` entries,
> so tiles carry plane order but no marker names. Always verify against an
> original whole-slide OME-TIFF. Tiles inherit plane order from it, so verifying
> the original is sufficient.

**If the originals are unavailable** (as happened here — tiles were already on
the cluster), verify indirectly by marker biology on one dense tile per cohort.
Load the tile, restrict to tissue pixels, and check that markers which should
co-localise do:

```python
import tifffile, numpy as np
a = tifffile.imread("tile_r1_c5.ome.tiff")[:, ::2, ::2].astype(np.float32)
tissue = a[0] > 8                      # plane 0 = DAPI
def r(i, j): return np.corrcoef(a[i][tissue], a[j][tissue])[0, 1]
r(17, 23)   # Iba1 ~ CD11b : expect HIGH (both myeloid)
r(17, 24)   # Iba1 ~ CD68  : expect HIGH
r(14, 7)    # GFAP ~ NeuN  : expect ~0 (distinct lineages)
r(17, 1)    # Iba1 ~ secondary-only control : expect ~0
```

Observed for ICV_T2_3 / IP_T3_3: Iba1~CD11b **+0.636 / +0.570**, Iba1~CD68
**+0.387 / +0.333**, GFAP~NeuN **+0.071 / −0.049**, Iba1~control **+0.016 /
+0.076**; per-plane coverage rank agreement **ρ = 0.844**. The myeloid triad
holds, lineage-distinct pairs sit at zero, controls are dark → plane order is
the same in both cohorts.

> **Pick a dense tile.** Sparse edge tiles and bright autofluorescent tissue
> fragments both give useless answers — in a fragment every channel lights up
> together and all correlations approach 1.

### 2.2 Tile

4000 px tiles, 100 px overlap. The overlap is **not optional**: `merge_brain.py`
deduplicates using each tile's non-overlapping inner zone and assumes it.

```bash
export COMET_IMAGE_DIR=/path/to/Nick_Raw_Images_8-26-26
export COMET_PYTHON="$COMET_PYTHON_TILE"
bash tiling/tile_all.sh                       # all images
bash tiling/tile_all.sh ICV-C1_3 IP-T2_3      # or named samples
```

`tile_all.sh` re-verifies the channel map per image, refuses to tile on
mismatch, and skips images already tiled. Under the hood:

```bash
python3 tiling/make_tiles.py INPUT.ome.tiff OUTDIR --tile 4000 --overlap 100
```

Each brain yields ~53–78 tiles (533 across the 8). Empty tiles are skipped
automatically (`crop[0].mean() < 1.0`).

**Disk:** the 8 source images are ~80 GB; tiles add ~100–200 GB because all 28
planes are preserved per tile. Tile a few brains, upload, delete locally, repeat.

---

## 3. Stage 2 — Upload to Cheaha

Two options. `scp` from your own terminal:

```bash
scp -r ICV_C1_3_tiles_4k dtyrrell@cheaha.rc.uab.edu:/data/user/dtyrrell/COMET_6-29-26/
```

…or the OnDemand web Files app (https://rc.uab.edu), which is easier for large
folders and survives flaky connections.

Push the pipeline and scripts too — **and note that changing the pipeline means
changing `run_array.sh` as well**, since it names the `.cppipe`. Uploading only
one of the two silently runs the old pipeline:

```bash
scp pipeline/COMET_track1plus2_28ch.cppipe hpc/run_array.sh hpc/submit_all_brains.sh analysis/merge_brain.py \
    dtyrrell@cheaha.rc.uab.edu:/data/user/dtyrrell/COMET_6-29-26/
```

> A `client_global_hostkeys_prove_confirm: server gave bad signature for ECDSA
> key 1: error in libcrypto` warning on macOS is **benign** — it concerns
> caching extra host keys for future rotation, not this session's
> authentication. A `REMOTE HOST IDENTIFICATION HAS CHANGED` message is an
> entirely different thing and must never be ignored.

---

## 4. Stage 3 — Connect and pre-flight

```bash
ssh dtyrrell@cheaha.rc.uab.edu          # password + Duo push
```

Verify the upload arrived intact — compare against locally computed checksums:

```bash
md5sum COMET_track1plus2_28ch.cppipe run_array.sh
grep -n "PIPELINE=" run_array.sh        # must name the pipeline you just uploaded
```

Confirm tiles and manifests. A missing `tile_manifest.json` will not surface
until the merge at the very end of a multi-hour run:

```bash
cd /data/user/dtyrrell/COMET_6-29-26 && mkdir -p logs
for d in *_tiles_4k; do echo "$d: $(ls $d/tile_r*_c*.ome.tiff 2>/dev/null | wc -l) tiles"; done
for d in *_tiles_4k; do [ -f "$d/tile_manifest.json" ] && echo "$d: manifest OK" || echo "$d: *** MANIFEST MISSING ***"; done
```

Reference counts for this dataset:

| brain | tiles | | brain | tiles |
|---|---|---|---|---|
| ICV_C1_3 | 53 | | IP_C1_3 | 73 |
| ICV_C3_3 | 66 | | IP_C2_3 | 64 |
| ICV_T1_3 | 64 | | IP_T2_3 | 78 |
| ICV_T2_3 | 58 | | IP_T3_3 | 77 |

`mkdir -p logs` matters — the array writes `logs/cp_%A_%a.out` and fails
immediately if the directory is absent.

---

## 5. Stage 4 — Smoke test before committing hours

Run three tiles from each cohort before launching ~500:

```bash
sbatch --array=1-3 --job-name=smoke_IP  --export=ALL,IMAGE=IP_T3_3  run_array.sh
sbatch --array=1-3 --job-name=smoke_ICV --export=ALL,IMAGE=ICV_T2_3 run_array.sh
squeue -u $USER
```

**Confirm the pipeline you think is running is the one running.** Count marker
columns per tile — do *not* `sort -u` across tiles, because the union of old
{CD3,CD4,CD8,NeuN} and the new 15 is still 15, which hides a stale pipeline:

```bash
for B in ICV_T2_3 IP_T3_3; do echo "=== $B ==="; new=0; old=0
  for f in results/$B/tile_r*_c*/MyExpt_Cells.csv; do
    n=$(head -1 "$f" | grep -o 'MeanIntensity_[A-Za-z0-9]*' | sort -u | wc -l)
    [ "$n" -ge 15 ] && new=$((new+1)) || old=$((old+1)); done
  echo "  15-marker tiles: $new"; echo "  4-marker tiles:  $old"; done
```

Check for errors:

```bash
grep -ilE "error|traceback|killed" logs/cp_*.err | head
```

---

## 6. Stage 5 — Launch the full run

```bash
cd /data/user/dtyrrell/COMET_6-29-26
bash submit_all_brains.sh ICV_C1_3 ICV_C3_3 ICV_T1_3 ICV_T2_3 \
                          IP_C1_3 IP_C2_3 IP_T2_3 IP_T3_3
```

16 submissions: 8 segmentation arrays (20 tiles concurrent each) plus 8
dependent merges. Brains run sequentially. Expect **6–12 hours** for 533 tiles.

Safe to log out — Slurm is server-side and does not care about your SSH session.

> **Concurrency and the memory QOS.** Brains are chained deliberately. All six
> arrays running at once (which happened accidentally, see Appendix A) hits
> `QOSMaxMemoryPerUser` at 20 × 32 GB per array and starves everything.

---

## 7. Stage 6 — Monitoring and the troubleshooting playbook

Health check — run this on every check-in. It releases held tasks and reports
whether the chain is alive:

```bash
squeue -u $USER -t PD -h -o "%i %r" | grep -i held | awk '{print $1}' | while read j; do scontrol release "$j"; done
squeue -u $USER -t PD -h -o "%i %r" | grep -iE "held|NeverSatisfied" || echo "chain healthy"
squeue -u $USER -h | wc -l | xargs echo "jobs left:"
```

### 7.1 `launch failed requeued held`

A node failed to start the task; Slurm requeued it and **held** it. Held jobs
never run until released, and while held the array never terminates.

```bash
scontrol release <JOBID>
# or release everything held:
squeue -u $USER -t PD -h -o "%i %r" | grep -i held | awk '{print $1}' | while read j; do scontrol release "$j"; done
```

Occurred twice here — once as a single task, once as a block of 13 (a node
dropping out mid-array).

### 7.2 `DependencyNeverSatisfied` on a merge job

The array finished with at least one task not in `COMPLETED`, so the merge's
`afterok` can never be met. **This also happens when a task was requeued and
then completed successfully** — the requeue transition alone breaks `afterok`.

Recover manually — verify tiles first, because `merge_brain.py` skips missing
files silently and would produce a quietly undercounted brain:

```bash
B=IP_C1_3; N=73
new=0; tot=0
for f in results/$B/tile_r*_c*/MyExpt_Cells.csv; do [ -f "$f" ] || continue; tot=$((tot+1))
  n=$(head -1 "$f" | grep -o 'MeanIntensity_[A-Za-z0-9]*' | sort -u | wc -l)
  [ "$n" -ge 15 ] && new=$((new+1)); done
echo "$B: $new / $tot (want $N / $N)"
# only if it matches:
scancel <MERGE_JOBID> && python3 merge_brain.py $B
```

With the `PREV_JOB=$SEG_JOB` fix this costs only that brain's merge. With the
old `PREV_JOB=$MERGE_JOB` it deadlocks every remaining brain.

### 7.3 Timeouts and hung nodes

A tile that exceeds the wall clock is killed with `TIMEOUT`.

```bash
sacct -j <ARRAYID> -X --format=JobID%18,State%20,Elapsed,ExitCode --noheader | grep -v COMPLETED
sacct -j <ARRAYID> --format=JobID%18,State%18,Elapsed,ReqMem,MaxRSS --units=G | grep -iE "fail|timeout|oom"
```

**Distinguishing a hung node from a slow tile is the hard part.** What does
*not* work: log length or output-file count. CellProfiler writes all CSVs at the
end, so a healthy task at 10 minutes and a wedged one both show 3 log lines and
only `file_list.txt`. Runtime is also unrelated to file size or tile position —
the tile that timed out here was 134 MB, not among the largest, and sat at row 0.

What *does* work:

- **Elapsed time against siblings still running.** Not against finished ones —
  a fast neighbour proves nothing.
- **`MaxRSS`.** A hung task stays near 0.2 GB; real segmentation climbs into GB.
  The 2-hour timeout showed `MaxRSS 0.19G` — it was stalled, not computing.
- **The decisive test:** requeue it. The tile that burned 2 hours on `c0142`
  completed in **3:45** on another node.

Be patient before acting — a tile here ran past 24 minutes while its siblings
took 8–13, and finished fine. Requeue only at a large multiple of sibling
runtime, because **requeuing costs that brain its automatic merge**.

```bash
# requeue, keeping it off a suspect node (three steps avoids a scheduling race)
scontrol requeuehold <JOBID>
scontrol update jobid=<JOBID> ExcNodeList=c0142
scontrol release <JOBID>
```

Pre-emptively exclude a known-bad node from a still-pending array:

```bash
scontrol update jobid=<ARRAYID> ExcNodeList=c0142
```

> A single bad node (`c0142`) caused four separate stalls in this run. If one
> node keeps wedging, report it to Research Computing — it will silently eat
> time on every future submission, not just this one.

### 7.4 Never let two runs write the same output

`run_array.sh` does `rm -f "$OUTDIR"/MyExpt_*.csv` before running. Two jobs on
the same tile interleave and corrupt it. Before resubmitting anything, confirm
nothing old is still queued:

```bash
squeue -u $USER -o "%.18i %.9P %.14j %.8T %.10M %.22R"
scancel <OLD_JOBIDS>
```

> Cancelling pending jobs can *release* others: `scancel -u $USER -t PENDING`
> cancelled a dead merge, which satisfied the `afterany` it was blocking and
> **started** the array we were trying to cancel. Always re-check `squeue` after
> a bulk cancel.

---

## 8. Stage 7 — Verify completion

Check marker counts, not just file counts — a stale merged folder from an older
run has the right number of CSVs and the wrong pipeline behind them:

```bash
cd /data/user/dtyrrell/COMET_6-29-26
for B in ICV_C1_3 ICV_C3_3 ICV_T1_3 ICV_T2_3 IP_C1_3 IP_C2_3 IP_T2_3 IP_T3_3; do
  M=results/$B/merged
  if [ -f "$M/${B}_Cells.csv" ]; then
    n=$(ls $M/*.csv | wc -l)
    m=$(head -1 $M/${B}_Cells.csv | grep -o 'MeanIntensity_[A-Za-z0-9]*' | sort -u | wc -l)
    c=$(($(wc -l < $M/${B}_Cells.csv) - 1))
    printf "%-10s %2d CSVs  %2d markers  %8d cells\n" "$B" "$n" "$m" "$c"
  else echo "$B: *** NO MERGED OUTPUT ***"; fi
done
sacct -u $USER --starttime YYYY-MM-DD -X --format=JobID%14,JobName%14,State%12 \
      --state=TIMEOUT,FAILED,NODE_FAIL,OUT_OF_MEMORY --noheader | head
```

Expected for this dataset — **13 CSVs for controls with zero CD8 T cells**
(CellProfiler writes no `CD8_Tcells.csv` at all), 14 otherwise:

| brain | CSVs | markers | cells |
|---|---|---|---|
| ICV_C1_3 | 13 | 15 | 110,104 |
| ICV_C3_3 | 13 | 15 | 132,751 |
| ICV_T1_3 | 14 | 15 | 129,508 |
| ICV_T2_3 | 14 | 15 | 101,812 |
| IP_C1_3 | 14 | 15 | 140,817 |
| IP_C2_3 | 14 | 15 | 119,095 |
| IP_T2_3 | 14 | 15 | 142,521 |
| IP_T3_3 | 14 | 15 | 161,451 |

Dedup drop rates should land at 3–7% per class.

---

## 9. Stage 8 — Download

**Run locally.** Transfers only `merged/*.csv`, skipping per-tile output
entirely; resumable — re-run the same command if it drops.

```bash
mkdir -p ~/Projects/ICV_experiment_Nick/run_8-26-26
rsync -avh --prune-empty-dirs --include='*/' --include='merged/*.csv' --exclude='*' \
  dtyrrell@cheaha.rc.uab.edu:/data/user/dtyrrell/COMET_6-29-26/results/ \
  ~/Projects/ICV_experiment_Nick/run_8-26-26/
```

127 files, ~11.7 GB. Download into a **new directory** — don't overwrite a
previous run's `merged_*` folders, which are your only backup.

Reshape `<BRAIN>/merged/` into the `merged_<BRAIN>/` layout the analysis expects:

```bash
cd ~/Projects/ICV_experiment_Nick/run_8-26-26
for d in */; do b=${d%/}; [ -d "$b/merged" ] && mv "$b/merged" "merged_$b" && rmdir "$b"; done
```

Verify against the counts the cluster reported, so a truncated transfer can't
slip through:

```bash
for B in ICV_C1_3 ICV_C3_3 ICV_T1_3 ICV_T2_3 IP_C1_3 IP_C2_3 IP_T2_3 IP_T3_3; do
  M=merged_$B
  printf "%-10s %8d cells  %2d markers\n" "$B" \
    "$(( $(wc -l < $M/${B}_Cells.csv) - 1 ))" \
    "$(head -1 $M/${B}_Cells.csv | grep -o 'MeanIntensity_[A-Za-z0-9]*' | sort -u | wc -l)"
done
```

---

## 10. Stage 9 — Local analysis

### 10.1 Build reduced-set per-sample CSVs

```bash
python3 analysis/make_cafe_csvs.py "csv files for analysis_8-26-26" run_8-26-26/merged_*
```

Writes `<OUT>/<CellType>/<Group>_<BrainID>.csv` for Astrocytes, Microglia,
Neurons and Tcells, plus `<OUT>/samples.csv` — a manifest carrying BrainID,
Route, Condition and Group, so downstream code never parses filenames.

Joins performed, none of them obvious:

- **Microglia / Astrocytes** — `ObjectSkeleton_NumberTrunks` comes from
  `<BRAIN>_MicrogliaSeeds.csv` / `<BRAIN>_AstrocyteSeeds.csv`, joined on
  `(tile_name, ObjectNumber)`. **Skeleton measurements live on the Seeds
  objects, not the parent.** Drops ~0.1% where numbering doesn't correspond.
- **Neurons / T cells** — shape and markers come from `<BRAIN>_Cells.csv` via
  `Parent_Cells → ObjectNumber`, matched within a tile. These objects carry no
  measurements of their own.
- **T cells** — re-gated per brain (see [Appendix B](#appendix-b--one-liner-reference-card)
  and `results_8-26-26.md`), with `Subtype` from the CD8/CD4 ratio.

`global_x` / `global_y` ride along on every cell for spatial analysis.

Gate overrides:

```bash
--tcell-k 2.5      # stricter CD3 (default 2.0 x brain median)
--cd8-ratio 3.0    # stricter CD8 call (default 2.0)
--tcell-abs        # revert to the pipeline's absolute CD3 >= 0.022
```

### 10.2 Cluster

```bash
D="csv files for analysis_8-26-26"
$COMET_PYTHON_SC analysis/cluster_cells.py "$D" Microglia
$COMET_PYTHON_SC analysis/cluster_cells.py "$D" Astrocytes
$COMET_PYTHON_SC analysis/cluster_cells.py "$D" Neurons
$COMET_PYTHON_SC analysis/cluster_cells.py "$D" Tcells --resolution 0.4 --neighbors 15
```

scale → PCA → Harmony (batch = brain) → UMAP → Leiden. Options:
`--resolution --neighbors --min-dist --metric --n-pcs --no-harmony --groups --outdir`.

Reproduce the older ICV-only comparison:

```bash
$COMET_PYTHON_SC analysis/cluster_cells.py "$D" Microglia --groups ICV-C,ICV-T
```

Observed scale and runtime (14-core Mac, Harmony is the bottleneck):

| cell type | cells | features | Harmony | clusters |
|---|---|---|---|---|
| Microglia | 218,135 | 17 | 385 s | 10 |
| Astrocytes | 120,352 | 17 | 210 s | 9 |
| Neurons | 460,410 | 20 | ~600 s | 10 |
| T cells | 652 | 20 | 0.2 s | 9 (5 at res 0.2) |

Neurons takes ~30 min end to end and peaks near 4–5 GB RSS. Run it in the
background and leave it alone.

Outputs per cell type in `<D>/<CellType>/output/`: UMAPs (clusters, group,
sample, split-by-group, T-cell subtype), marker dotplot, cluster-frequency
barplot, count / frequency / median-expression tables, and
`adata_raw.h5ad` + `adata_harmony.h5ad` + `adata_clustered.h5ad`.

> For 652 T cells, resolution 0.4 gives 9 clusters (~72 cells each) — finer than
> that sample size supports. A second pass at `--resolution 0.2 --outdir .../output_res0.2`
> gives 5, which is more defensible.

---

## Appendix A — Incident log from the 8-26-26 run

Everything that went wrong, in order. Nothing here was a pipeline bug; all of it
was cluster behaviour or operator sequencing.

| # | What happened | Cause | Fix |
|---|---|---|---|
| 1 | 1 task `launch failed requeued held` | node failed to start it | `scontrol release` |
| 2 | `39831123_2` `TIMEOUT` at 02:00:07, `MaxRSS 0.19G` | node `c0142` wedged, not a hard tile | moved to `short`/12h; tile then ran in **3:45** elsewhere |
| 3 | ICV_C3_3 merge `DependencyNeverSatisfied`, all 6 remaining brains stalled | merge used `afterok`; a permanently-pending merge never fires the `afterany` the next brain waits on | changed chaining to `PREV_JOB=$SEG_JOB`; manual merge |
| 4 | 6 old arrays suddenly ran concurrently, new jobs stuck on `QOSMaxMemoryPerUser` | `scancel -t PENDING` cancelled the dead merge, releasing the arrays it blocked | cancelled old arrays by ID before the new ones could collide |
| 5 | IP_C1_3 merge `DependencyNeverSatisfied` though **all 73 tasks `COMPLETED`** | a requeue breaks `afterok` even on eventual success | manual merge |
| 6 | IP_C2_3 task stuck 78 min (siblings 2–6 min) on `c0142` | hung node | requeue → completed quickly |
| 7 | 13 IP_T3_3 tasks held at once | node dropped out mid-array | bulk `scontrol release` |
| 8 | IP_T3_3 task at 3:06 with array otherwise done, on `c0142` | hung node | `requeuehold` + `ExcNodeList=c0142` + `release` |
| 9 | IP_T2_3 task ran 24 min vs siblings' 8–13 | genuinely slow tile | **left alone; finished fine** |

Lessons that changed the code in this repo:

1. **Chain brains on the array, not the merge** (`hpc/submit_all_brains.sh`).
2. **2 hours is not enough headroom** (`hpc/run_array.sh` → `short`/12h).
3. **Requeue is not free** — it forfeits that brain's automatic merge.
4. **Verify with marker counts, never file counts** — stale merged folders look complete.
5. **Patience beats requeuing** — incident 9 would have cost a merge for nothing.

## Appendix B — One-liner reference card

```bash
# health check (release held, report stalls, count remaining)
squeue -u $USER -t PD -h -o "%i %r" | grep -i held | awk '{print $1}' | while read j; do scontrol release "$j"; done
squeue -u $USER -t PD -h -o "%i %r" | grep -iE "held|NeverSatisfied" || echo "chain healthy"

# what failed, since a date
sacct -u $USER --starttime 2026-08-26 -X --format=JobID%14,JobName%14,State%12,Elapsed \
      --state=TIMEOUT,FAILED,NODE_FAIL,OUT_OF_MEMORY --noheader

# is a running task working or hung?  (near-zero MaxRSS => hung)
sstat -j <JOBID> --format=JobID%16,AveCPU,MaxRSS,MaxVMSize --noheader

# which tile is array index N?
ls <BRAIN>_tiles_4k/tile_r*_c*.ome.tiff | sort | sed -n Np

# block until a job finishes
while squeue -u $USER -h -j <JOBID> | grep -q .; do sleep 120; done; echo DONE

# per-brain progress: how many tiles are on the current pipeline?
B=IP_T3_3; new=0; tot=0
for f in results/$B/tile_r*_c*/MyExpt_Cells.csv; do [ -f "$f" ] || continue; tot=$((tot+1))
  n=$(head -1 "$f" | grep -o 'MeanIntensity_[A-Za-z0-9]*' | sort -u | wc -l)
  [ "$n" -ge 15 ] && new=$((new+1)); done; echo "$B: $new / $tot"
```

### T-cell gating, in one place

| | pipeline | used here | why |
|---|---|---|---|
| T cell | `CD3 >= 0.022` absolute | `CD3 >= 2.0 x brain median` | absolute gate swung counts **78x** across brains; a control outranked every treated brain |
| CD8 call | `CD8_CD4_ratio >= 1.0` | `>= 2.0` | 1.0 separates cleanly in ICV but not IP, where the ratio distribution is compressed |
| NeuN / Area / FormFactor | `<=0.050` / `100–950` / `>=0.70` | **unchanged** | validated, not batch-sensitive |

Validation against the CD8-knockout design: 3 of 4 control brains call **exactly
zero** CD8; treated brains call 4–9; 652 T cells total (625 CD4 / 27 CD8).

Full evidence in [`results_8-26-26.md`](results_8-26-26.md).
