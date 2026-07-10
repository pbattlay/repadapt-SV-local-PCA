#!/bin/bash
#SBATCH --job-name=nuke_s2s3
#SBATCH --output=errout/%x.out
#SBATCH --error=errout/%x.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:10:00

# =============================================================
# Purge all script2 + script3 outputs for one species so the next
# run starts clean. Files from different runs carry run-specific names
# (chr:start-end), so without this they accumulate rather than overwrite.
#
# RUN THIS ONCE, BEFORE LAUNCHING THE SCRIPT 2 ARRAY — never from inside
# script2 itself. script2 is a SLURM array (one task per chromosome running
# in parallel); a purge inside it would have each task delete the others'
# output mid-run. Chain it as a dependency so the array cannot start until
# the purge succeeds:
#
#   clean=$(sbatch --parsable \
#       -D "${base_dir}/${raid}" \
#       "${base_dir}/scripts/nuke_script2_and_3_results.sh" \
#       "${base_dir}/config_${raid}.sh")
#   sbatch --dependency=afterok:${clean} \
#       --array="1-$(wc -l < ${base_dir}/${raid}/ref/chrs.txt)" \
#       -D "${base_dir}/${raid}" \
#       "${base_dir}/scripts/script2.sh" "${base_dir}/config_${raid}.sh"
#
# Or just run it by hand before the array — same effect.
# =============================================================

set -euo pipefail

source "$1"

# Refuse to run if the species identity is not fully resolved (guards against
# an unset variable turning a relative purge into something catastrophic).
[[ -n "${base_dir:-}" && -n "${raid:-}" ]] || {
    echo "ERROR: base_dir and/or raid unset — refusing to purge." >&2
    exit 1
}

# Operate from inside the verified species directory; set -e aborts here if it
# does not exist, so the find calls below can only ever touch this species.
cd "${base_dir}/${raid}"

echo "[$(date)] Purging script2/3 outputs for ${raid} in $(pwd)"

# Delete files only (any depth); keep the directory tree intact so script2 can
# write straight into mds/plots, local_pca/raw, breakpoints/plots, etc. without
# needing to recreate them.
for d in mds local_pca breakpoints final; do
    if [[ -d "${d}" ]]; then
        n=$(find "${d}" -mindepth 1 -type f | wc -l)
        find "${d}" -mindepth 1 -type f -delete
        echo "  ${d}/: removed ${n} file(s)"
    else
        echo "  ${d}/: not present (skipped)"
    fi
done

echo "[$(date)] Purge complete — ${raid} ready for a fresh script2 run"
