#!/bin/bash
# Release any tasks Slurm parked in a held state after a launch failure, and
# report anything stuck on an unsatisfiable dependency.
#
# Run periodically while a long array is in flight:
#   watch -n 300 bash release_held.sh
# or from a login node every few hours.
N=0
for J in $(squeue -u "$USER" -t PD -h -o "%i %r" | grep -i held | awk '{print $1}'); do
    scontrol release "$J" && echo "released $J" && N=$((N+1))
done
[[ $N -eq 0 ]] && echo "nothing held"
STUCK=$(squeue -u "$USER" -t PD -h -o "%i %r" | grep -i NeverSatisfied)
if [[ -n "$STUCK" ]]; then
    echo "WARNING: dependency can never be satisfied -- these merges need running by hand:"
    echo "$STUCK"
fi
echo "jobs remaining: $(squeue -u "$USER" -h | wc -l)"
