#!/bin/bash
# Submit array job + dependent merge job for each remaining brain.
# Each brain's array runs in parallel (within its 20-task cap).
# Brains submit sequentially: brain N+1 starts after brain N's array completes.
# This caps total cluster concurrency at ~20 and provides clean per-brain failure isolation.

set -e

# --- config (see .env.example / README for details) ---
PROJECT="${COMET_PROJECT_ROOT:?Set COMET_PROJECT_ROOT env var to your project directory}"
PARTITION="${SLURM_PARTITION:-express}"
# --------------------------------------------------------

cd "$PROJECT"
mkdir -p logs

PREV_JOB=""
for IMAGE in "$@"; do
    DIR=${IMAGE}_tiles_4k
    if [[ ! -d $DIR ]]; then
        echo "SKIP $IMAGE: missing $DIR"
        continue
    fi
    N=$(ls $DIR/tile_r*_c*.ome.tiff 2>/dev/null | wc -l)
    if [[ $N -eq 0 ]]; then
        echo "SKIP $IMAGE: no tiles in $DIR"
        continue
    fi
    mkdir -p results/${IMAGE}

    # Dependency: each brain's array waits for the previous brain's ARRAY, not
    # its merge. Chaining on the merge deadlocks: the merge uses afterok, so a
    # single failed/requeued tile leaves it permanently PENDING with
    # DependencyNeverSatisfied -- and since it never terminates, the afterany
    # that the next brain waits on never fires either, stalling every remaining
    # brain. Hanging the chain off the array (afterany, fires on any outcome)
    # means one bad tile costs that brain its automatic merge and nothing more.
    DEP_FLAG=""
    [[ -n "$PREV_JOB" ]] && DEP_FLAG="--dependency=afterany:$PREV_JOB"

    # Submit segmentation array
    SEG_JOB=$(sbatch --parsable $DEP_FLAG \
                     --array=1-${N}%20 \
                     --job-name=cp_${IMAGE} \
                     --partition=$PARTITION \
                     --export=ALL,IMAGE=${IMAGE} \
                     hpc/run_array.sh)
    echo "Submitted $IMAGE segmentation as job $SEG_JOB ($N tiles, after ${PREV_JOB:-now})"

    # Submit dependent merge job
    MERGE_JOB=$(sbatch --parsable \
                       --dependency=afterok:$SEG_JOB \
                       --job-name=merge_${IMAGE} \
                       --partition=$PARTITION \
                       --time=00:30:00 \
                       --mem=8G \
                       --cpus-per-task=1 \
                       --output=logs/merge_%j.out \
                       --error=logs/merge_%j.err \
                       --wrap="python3 analysis/merge_brain.py ${IMAGE}")
    echo "  -> merge as job $MERGE_JOB (runs after $SEG_JOB succeeds)"

    PREV_JOB=$SEG_JOB
done
echo ""
echo "All submitted. Watch with: squeue -u \$USER"
echo ""
echo "If a brain's merge shows DependencyNeverSatisfied, its array had a failed"
echo "or requeued tile. Recover with (see docs/RUNBOOK.md):"
echo "  scancel <MERGE_JOBID> && python3 analysis/merge_brain.py <BRAIN>"
