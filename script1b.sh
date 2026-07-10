#!/bin/bash
#SBATCH --job-name=script1b
#SBATCH --output=errout/%x.out
#SBATCH --error=errout/%x.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00

set -euo pipefail

source "$1"
echo "[$(date)] Script 1b — ${raid}: sample exclusion"

module load plink/2.0-alpha

bed_prefix="vcf/${raid}.filtered"
if [[ ! -s "${bed_prefix}.bed" ]]; then
    echo "ERROR: ${bed_prefix}.bed not found — run script1a first." >&2
    exit 1
fi

# Sample exclusion drops individuals only. Site-level missingness (--geno) is
# NOT recomputed, so the site set stays as determined on the full sample in 1a.
# (MAF *is* recomputed later, in script1c, because allele frequencies change
# once samples are removed.)
#
# The user-facing REMOVE_SAMPLES file is a plain list of IIDs (one per line,
# matching vcf/plots/pca_preimpute.tsv). plink2 --remove wants FID<TAB>IID and
# all FIDs are 0 (--const-fid 0), so we build a temporary 2-column file.
if [[ -n "${REMOVE_SAMPLES:-}" && -s "${REMOVE_SAMPLES}" ]]; then
    n_removed=$(grep -cve '^[[:space:]]*$' "${REMOVE_SAMPLES}")
    echo "[$(date)]   Excluding ${n_removed} sample(s) listed in ${REMOVE_SAMPLES}"

    remove_tmp="vcf/.${raid}.remove_fidiid.tmp"
    awk 'NF { print 0"\t"$1 }' "${REMOVE_SAMPLES}" > "${remove_tmp}"

    plink2 \
        --bfile "${bed_prefix}" \
        --remove "${remove_tmp}" \
        --threads 2 \
        --allow-extra-chr \
        --make-bed \
        --out "vcf/${raid}.kept"

    rm -f "${remove_tmp}"

    N_IND=$(wc -l < "vcf/${raid}.kept.fam")
    echo "[$(date)] Script 1b complete — ${raid}"
    echo "  ${N_IND} individuals retained -> vcf/${raid}.kept.{bed,bim,fam}"
else
    echo "[$(date)]   No samples to exclude (REMOVE_SAMPLES unset or empty)"
    echo "[$(date)] Script 1b complete — ${raid}"
    echo "  No exclusion applied; script1c will read vcf/${raid}.filtered directly"
fi

echo "  --> Run script1c (array job) to MAF-refilter, split, and impute."
