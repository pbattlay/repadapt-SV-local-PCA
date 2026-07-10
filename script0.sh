#!/bin/bash
#SBATCH --job-name=script0_map
#SBATCH --output=script0_map_%A_%a.out
#SBATCH --error=script0_map_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

# =====================================================================
# script0 — mappability mask generator (per-reference, DECOUPLED)
#
# Builds a per-base, uniquely-mappable keep-list (genmap) for EVERY FASTA in a
# directory, writing <prefix>.mappable.bed NEXT TO each FASTA. The mask is
# reference-level, so species sharing a reference share one .bed — just point
# each species config's MAPPABILITY_MASK at the right file.
#
#   Usage (recommended — one array task per FASTA; count comes from script0 itself
#   so it can never drift from the manifest the worker indexes):
#     N=$(bash script0.sh --count /path/to/fasta_dir)
#     sbatch --array="1-${N}" script0.sh /path/to/fasta_dir
#
#   Or as a single serial job over every FASTA in the directory:
#     sbatch script0.sh /path/to/fasta_dir        # loops internally
#     bash   script0.sh /path/to/fasta_dir        # interactive
#
# Output is GENOME-WIDE (all contigs in the FASTA); there is no chrs.txt here
# (this is reference-level, not per-species). The old script0 per-window retention
# printout is dropped — masking simply removes SNPs in script1a, and windows are
# then gated only by MIN_SNPS_PER_WINDOW in script2 (sparse masked regions fall
# out there). No separate window-level mappability step.
#
# Coordinates are 0-based BED (genmap/bedtools native; consumed by plink2
# `--extract bed0`). This is the ONE pipeline input that is 0-based, unlike the
# 1-based truth/manual .tsv inputs.
# =====================================================================

# ---- genmap parameters (script0-local; no per-species config) -------------
GENMAP_BIN="/home/pbattlay/.conda/envs/genmap/bin/genmap"  # genmap (conda env, NOT a module)
GENMAP_K=100                # k-mer length (~ Illumina WGS read length)
GENMAP_E=2                  # mismatches allowed when scoring k-mer uniqueness
MAPPABILITY_THRESHOLD=1     # per-base cutoff (1 = uniquely mappable). Lower toward 0.5 to
                            # keep 2-copy regions; check retention (script1a) BEFORE going below 1.
FASTA_GLOBS=("*.fasta" "*.fa" "*.fasta.gz" "*.fa.gz")
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---- inputs ----------------------------------------------------------------
# Usage:  script0.sh [--count] [fasta_dir]
#   fasta_dir defaults to the working directory, so `sbatch --chdir=DIR script0.sh`
#   needs no positional arg. `--count` prints the number of FASTAs (for sizing the
#   array) from the SAME manifest, then exits — dependency-free, safe to call at submit.
COUNT_ONLY=0
if [[ "${1:-}" == "--count" ]]; then COUNT_ONLY=1; shift; fi
FASTA_DIR="${1:-.}"
[[ -d "${FASTA_DIR}" ]] || { echo "ERROR: not a directory: ${FASTA_DIR}" >&2; exit 1; }

# ---- build a stable, sorted FASTA manifest (single source of truth) -------
mapfile -t FASTAS < <(
  { for g in "${FASTA_GLOBS[@]}"; do
      find "${FASTA_DIR}" -maxdepth 1 -type f -name "$g" 2>/dev/null
    done; } | sort -u
)
N=${#FASTAS[@]}

if [[ "${COUNT_ONLY}" -eq 1 ]]; then echo "${N}"; exit 0; fi
[[ "${N}" -gt 0 ]] || { echo "ERROR: no FASTA files (${FASTA_GLOBS[*]}) in ${FASTA_DIR}" >&2; exit 1; }
echo "[$(date)] script0 — ${N} FASTA(s) found in ${FASTA_DIR}"

# bedtools is a module here (sort + merge); genmap is a conda binary. Adjust to your cluster.
module load bedtools
[[ -x "${GENMAP_BIN}" ]]        || { echo "ERROR: genmap not found/executable at ${GENMAP_BIN}" >&2; exit 1; }
command -v bedtools >/dev/null  || { echo "ERROR: bedtools not on PATH (fix the module load)." >&2; exit 1; }

# ---- decide which FASTA(s) this invocation handles ------------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
  (( idx >= 0 && idx < N )) || { echo "ERROR: array index ${SLURM_ARRAY_TASK_ID} out of range 1..${N}" >&2; exit 1; }
  TARGETS=( "${FASTAS[$idx]}" )
else
  echo "  No SLURM_ARRAY_TASK_ID set — processing all ${N} serially."
  echo "  To parallelise instead:  sbatch --array=1-${N} $0 ${FASTA_DIR}"
  TARGETS=( "${FASTAS[@]}" )
fi

# ---- per-FASTA worker ------------------------------------------------------
process_one() {
  local fasta="$1" dir prefix fa work="" bg mask idxdir n bp alt
  dir="$(dirname "${fasta}")"
  prefix="$(basename "${fasta}")"; prefix="${prefix%.gz}"; prefix="${prefix%.*}"
  mask="${dir}/${prefix}.mappable.bed"
  bg="${dir}/${prefix}.mappability.bedgraph"
  echo "[$(date)]   ${prefix}"

  # genmap needs a decompressed FASTA
  if [[ "${fasta}" == *.gz ]]; then
    work="$(mktemp --suffix=.fasta)"; gunzip -c "${fasta}" > "${work}"; fa="${work}"
  else
    fa="${fasta}"
  fi

  # STEP 1+2: genmap mappability -> per-base keep-list (cached bedgraph -> skip)
  if [[ ! -s "${bg}" ]]; then
    idxdir="$(mktemp -d)"; rmdir "${idxdir}"          # genmap errors if the index dir exists
    echo "[$(date)]     genmap index"
    "${GENMAP_BIN}" index -F "${fa}" -I "${idxdir}" > /dev/null
    echo "[$(date)]     genmap map (K=${GENMAP_K}, E=${GENMAP_E})"
    "${GENMAP_BIN}" map -K "${GENMAP_K}" -E "${GENMAP_E}" \
        -I "${idxdir}" -O "${dir}/${prefix}.mappability" -bg --threads "${THREADS}"
    rm -rf "${idxdir}"
    if [[ ! -s "${bg}" ]]; then                       # some builds name the file differently
      alt="$(ls "${dir}/${prefix}.mappability"*.bedgraph "${dir}/${prefix}.mappability"*.bg 2>/dev/null | head -n1 || true)"
      [[ -n "${alt:-}" ]] && bg="${alt}"
    fi
  else
    echo "[$(date)]     bedgraph exists — skipping genmap"
  fi
  [[ -s "${bg}" ]] || { echo "ERROR: genmap produced no bedgraph for ${prefix}" >&2; [[ -n "${work}" ]] && rm -f "${work}"; return 1; }

  # Per-base keep-list (genome-wide); 0-based BED, touching intervals merged.
  awk -v t="${MAPPABILITY_THRESHOLD}" 'BEGIN{OFS="\t"} $4>=t {print $1,$2,$3}' "${bg}" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${mask}"

  [[ -n "${work}" ]] && rm -f "${work}"
  n=$(wc -l < "${mask}"); bp=$(awk '{s+=$3-$2} END{print s+0}' "${mask}")
  echo "[$(date)]     -> ${mask}  (${n} intervals, ${bp} bp mappable)"
  [[ "${n}" -gt 0 ]] || echo "WARNING: empty mask for ${prefix} — threshold too strict or no callable sequence." >&2
}

for f in "${TARGETS[@]}"; do process_one "${f}"; done
echo "[$(date)] script0 complete"
