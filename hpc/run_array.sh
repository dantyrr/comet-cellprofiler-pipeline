#!/bin/bash
#SBATCH --job-name=cp_array
#SBATCH --partition=express
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/cp_%A_%a.out
#SBATCH --error=logs/cp_%A_%a.err

# NOTE: this module name/version is specific to UAB Cheaha. Adjust or
# remove for other clusters (e.g. `module load apptainer` elsewhere).
module load Singularity/3.5.2-GCC-5.4.0-2.26

# --- config (see .env.example / README for details) ---
PROJECT="${COMET_PROJECT_ROOT:?Set COMET_PROJECT_ROOT env var to your project directory}"
PARTITION="${SLURM_PARTITION:-express}"
# --------------------------------------------------------

PIPELINE=$PROJECT/pipeline/COMET_track1plus2_28ch.cppipe
SIF=$PROJECT/cellprofiler.sif

# IMAGE must be passed in via --export=ALL,IMAGE=<brain_name>
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

echo "=== Processing $IMAGE / $TILE_NAME ==="
echo "Input:  $TILE"
echo "Output: $OUTDIR"

# NOTE: --partition is set at #SBATCH time above (Slurm doesn't support
# overriding partition/mem/time via env var at submit time for this script
# when run standalone with `sbatch run_array.sh`). When submitted via
# submit_all_brains.sh, that script passes --partition=$PARTITION explicitly
# on the sbatch command line, which overrides the #SBATCH default here.
# If running this script directly with `sbatch run_array.sh`, edit the
# #SBATCH lines above to match your cluster (partition name, mem, time).

singularity exec --bind "$PROJECT" "$SIF" \
    cellprofiler -c -r \
    -p "$PIPELINE" \
    --file-list="$FILE_LIST" \
    -o "$OUTDIR/"

echo "=== Done with $IMAGE / $TILE_NAME ==="
