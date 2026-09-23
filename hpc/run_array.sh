#!/bin/bash
#SBATCH --job-name=cp_array
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --requeue
#SBATCH --output=logs/cp_%A_%a.out
#SBATCH --error=logs/cp_%A_%a.err

# --------------------------------------------------------------------------
# Self-healing against flaky nodes.
#
# Two distinct failure modes were seen in the 8-26-26 run:
#
#  1. "launch failed requeued held" -- Slurm requeues the task itself but parks
#     it HELD, where it sits forever. --requeue above permits the retry; the
#     hold still needs `scontrol release` (see release_held.sh).
#
#  2. A node accepts the job and then wedges. Slurm sees a healthy RUNNING job
#     and does nothing; one tile burned 3h06 this way, another 2h. The tell is
#     MaxRSS stuck near 0.2G while real segmentation climbs into GB. Slurm
#     cannot detect this, so the job times ITSELF out and requeues.
#
# Healthy tiles took 2-30 min. CP_TIMEOUT_MIN=60 is ~2x the slowest legitimate
# tile observed, so it fires on hangs and not on slow work. On the final
# attempt the timeout is dropped entirely, so a genuinely slow tile still
# completes under the 12h wall clock rather than looping forever.
# --------------------------------------------------------------------------
MAX_RESTARTS=${CP_MAX_RESTARTS:-3}
TIMEOUT_MIN=${CP_TIMEOUT_MIN:-60}
RESTARTS=${SLURM_RESTART_COUNT:-0}

module load Singularity/3.5.2-GCC-5.4.0-2.26

PROJECT="${COMET_PROJECT_ROOT:?Set COMET_PROJECT_ROOT to your project directory}"
PIPELINE=$PROJECT/pipeline/COMET_28ch_seeded.cppipe
SIF=$PROJECT/cellprofiler.sif

if [[ -z "$IMAGE" ]]; then
    echo "ERROR: IMAGE variable not set. Submit with --export=ALL,IMAGE=<brain_name>" >&2
    exit 1
fi

TILES_DIR=$PROJECT/${IMAGE}_tiles_4k
RESULTS_DIR=$PROJECT/results/${IMAGE}

TILE=$(ls $TILES_DIR/tile_r*_c*.ome.tiff | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [[ -z "$TILE" ]]; then
    echo "ERROR: No tile for array index $SLURM_ARRAY_TASK_ID in $TILES_DIR" >&2
    exit 1
fi
TILE_NAME=$(basename "$TILE" .ome.tiff)

OUTDIR=$RESULTS_DIR/$TILE_NAME
mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/MyExpt_*.csv

FILE_LIST="$OUTDIR/file_list.txt"
echo "$TILE" > "$FILE_LIST"

echo "=== $IMAGE / $TILE_NAME ==="
echo "Node:     $(hostname)"          # recorded so repeat offenders are findable
echo "Attempt:  $((RESTARTS+1)) of $((MAX_RESTARTS+1))"
echo "Input:    $TILE"
echo "Output:   $OUTDIR"
echo "Started:  $(date)"
T_START=$(date +%s)

# NOTE: timeout execs singularity directly. Do not wrap this in `bash -c`,
# which starts a fresh shell that does not inherit non-exported variables --
# $SIF/$PIPELINE/$FILE_LIST/$OUTDIR would all arrive empty.
if [[ "$RESTARTS" -ge "$MAX_RESTARTS" ]]; then
    echo "Final attempt: running WITHOUT the internal timeout (12h wall clock applies)."
    singularity exec --bind "${COMET_BIND:-$PROJECT}" "$SIF" \
        cellprofiler -c -r -p "$PIPELINE" --file-list="$FILE_LIST" -o "$OUTDIR/"
    RC=$?
else
    timeout -k 30s "${TIMEOUT_MIN}m" \
        singularity exec --bind "${COMET_BIND:-$PROJECT}" "$SIF" \
        cellprofiler -c -r -p "$PIPELINE" --file-list="$FILE_LIST" -o "$OUTDIR/"
    RC=$?
fi

ELAPSED=$(( $(date +%s) - T_START ))

# A hung tile can come back as 124 (timeout gave up), or 137/143 when
# --kill-after escalates to KILL and timeout propagates the signal status
# instead. Relying on the exit code alone missed a real hang in the 8-28-26
# run -- tile_r7_c3 ran 01:00:53 against a 60m limit, exited 137, and was not
# retried. Elapsed time is the reliable signal, so gate on that too.
if [[ $RC -eq 124 || $RC -eq 137 || $RC -eq 143 ]] \
   || { [[ $RC -ne 0 ]] && [[ $ELAPSED -ge $((TIMEOUT_MIN*60 - 60)) ]]; }; then
    echo "=== NO PROGRESS: ran ${ELAPSED}s (limit ${TIMEOUT_MIN}m), rc=$RC on $(hostname) -- hung node ==="
    echo "=== requeueing (attempt $((RESTARTS+1)) of $((MAX_RESTARTS+1))) ==="
    scontrol update jobid="$SLURM_JOB_ID" ExcNodeList="$(hostname)" 2>/dev/null \
        && echo "excluded $(hostname) for the retry"
    scontrol requeue "$SLURM_JOB_ID"
    exit 1
fi

if [[ $RC -ne 0 ]]; then
    echo "=== CellProfiler exited $RC on $(hostname) ==="
    exit $RC
fi

echo "Finished: $(date)"
echo "=== Done with $IMAGE / $TILE_NAME ==="
