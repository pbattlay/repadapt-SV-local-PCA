#!/bin/bash
#SBATCH --job-name=script1c
#SBATCH --output=errout/%x-%a.out
#SBATCH --error=errout/%x-%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00

set -euo pipefail

source "$1"
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ref/chrs.txt)
echo "[$(date)] Script 1c — ${chr}: MAF refilter + export + imputation + BED"

module load plink/2.0-alpha java

# Pick the analysis BED: script1b only writes ${raid}.kept when samples were
# excluded; otherwise read the full filtered set written by script1a.
if [[ -n "${REMOVE_SAMPLES:-}" && -s "${REMOVE_SAMPLES}" ]]; then
    bed_in="vcf/${raid}.kept"
else
    bed_in="vcf/${raid}.filtered"
fi
if [[ ! -s "${bed_in}.bed" ]]; then
    echo "ERROR: ${bed_in}.bed not found — run script1a (and script1b if excluding samples) first." >&2
    exit 1
fi

# =============================================================
# STEP 1: EXTRACT CHROMOSOME + RE-APPLY MAF -> per-chr VCF
# =============================================================
# --maf is recomputed on the retained sample set (frequencies change after
# --remove). --geno is intentionally NOT re-applied: per-site missingness was
# fixed on the full sample in script1a.
echo "[$(date)]   Exporting ${chr} (MAF >= ${MAF})"
plink2 \
    --bfile "${bed_in}" \
    --chr "${chr}" \
    --threads 4 \
    --allow-extra-chr \
    --maf "${MAF}" \
    --export vcf bgz id-paste=iid \
    --out "vcf/${chr}.filtered"

# =============================================================
# STEP 2: BEAGLE IMPUTATION
# =============================================================
echo "[$(date)]   Imputing with BEAGLE"
java \
    -server \
    -Xms${BEAGLE_MEM}g -Xmx${BEAGLE_MEM}g \
    -XX:+UseParallelGC -XX:ParallelGCThreads=${BEAGLE_THREADS} \
    -jar "${BEAGLE_JAR}" \
    out="vcf/${chr}_imputed" \
    gt="vcf/${chr}.filtered.vcf.gz" \
    nthreads=${BEAGLE_THREADS} \
    window=${BEAGLE_WINDOW}
echo "  Imputation complete"

# =============================================================
# STEP 3: PLINK BED CONVERSION
# =============================================================
echo "[$(date)]   Converting to PLINK BED"
plink2 \
    --vcf "vcf/${chr}_imputed.vcf.gz" \
    --threads 4 \
    --set-missing-var-ids @:# \
    --const-fid 0 \
    --allow-extra-chr \
    --make-bed \
    --out "vcf/${chr}"

N_IND=$(wc -l < "vcf/${chr}.fam")
N_SNP=$(wc -l < "vcf/${chr}.bim")
echo "  BED: ${N_IND} individuals, ${N_SNP} SNPs"

# =============================================================
# CLEANUP intermediates
# =============================================================
rm -f \
    "vcf/${chr}.filtered.vcf.gz" \
    "vcf/${chr}_imputed.vcf.gz" \
    "vcf/${chr}_imputed.log"
echo "[$(date)]   Intermediates removed"

echo "[$(date)] Script 1c complete — ${chr}"
