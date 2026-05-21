#!/bin/bash
#SBATCH --job-name=script1
#SBATCH --output=errout/%x-%a.out
#SBATCH --error=errout/%x-%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00

set -euo pipefail

source "$1"
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ref/chrs.txt)
echo "[$(date)] Script 1 — ${chr}"

module load plink/2.0-alpha bcftools java

# =============================================================
# STEP 1: VCF FILTERING
# =============================================================
echo "[$(date)]   Filtering VCF"
plink2 \
    --vcf "${vcf}" \
    --chr "${chr}" \
    --threads 4 \
    --set-missing-var-ids @:# \
    --const-fid 0 \
    --allow-extra-chr \
    --var-min-qual ${MIN_QUAL} \
    --vcf-min-gq   ${MIN_GQ} \
    --vcf-min-dp   ${MIN_DP} \
    --geno          ${GENO} \
    --max-alleles 2 \
    --maf           ${MAF} \
    --export vcf bgz id-paste=iid \
    --out "vcf/${chr}"

N_SITES=$(zcat "vcf/${chr}.vcf.gz" | grep -vc "^#")
echo "  SNPs after filtering: ${N_SITES}"

# =============================================================
# STEP 2: SORT
# =============================================================
echo "[$(date)]   Sorting"
bcftools sort "vcf/${chr}.vcf.gz" -Oz -o "vcf/${chr}.sorted.vcf.gz"

# =============================================================
# STEP 3: BEAGLE IMPUTATION
# =============================================================
echo "[$(date)]   Imputing with BEAGLE"
java \
    -server \
    -Xms${BEAGLE_MEM}g -Xmx${BEAGLE_MEM}g \
    -XX:+UseParallelGC -XX:ParallelGCThreads=${BEAGLE_THREADS} \
    -jar "${BEAGLE_JAR}" \
    out="vcf/${chr}_imputed" \
    gt="vcf/${chr}.sorted.vcf.gz" \
    nthreads=${BEAGLE_THREADS} \
    window=${BEAGLE_WINDOW}
echo "  Imputation complete"

# =============================================================
# STEP 4: PLINK BED CONVERSION
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
    "vcf/${chr}.vcf.gz" \
    "vcf/${chr}.vcf.gz.tbi" \
    "vcf/${chr}.sorted.vcf.gz" \
    "vcf/${chr}_imputed.vcf.gz" \
    "vcf/${chr}_imputed.log"
echo "[$(date)]   Intermediates removed"

echo "[$(date)] Script 1 complete — ${chr}"
