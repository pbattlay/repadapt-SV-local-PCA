### REPADAPT SV — LOCAL PCA - CONFIG
account="om62"
base_dir="/home/pbattlay/om62_scratch2/rasvlpca_auto_2"
raid="rawg0053"
vcf="/home/pbattlay/om62_scratch2/aa_2025/vcf/aa_0625_DP3.vcf.gz"

# Reference
min_chr_length=10000000    # minimum scaffold length to include (bp); filters out small scaffolds

# VCF filtering
MIN_QUAL=30    # minimum site QUAL score
MIN_GQ=20      # minimum sample genotype quality
MIN_DP=5       # minimum sample read depth
GENO=0.3       # maximum missingness per site (plink2 --geno)
MAF=0.01       # minimum minor allele frequency

# BEAGLE
BEAGLE_JAR="${base_dir}/scripts/beagle.jar"    # path to beagle.jar
BEAGLE_MEM=68                                  # JVM heap memory (GB)
BEAGLE_THREADS=4                               # imputation threads
BEAGLE_WINDOW=20.0                             # sliding window size for imputation (cM)

# Windowing
WINDOW_SIZE=10000          # local PCA window size (bp); scale MIN_SNPS_PER_WINDOW proportionally
MIN_SNPS_PER_WINDOW=20     # minimum SNPs required for a window to be included (~0.002 × WINDOW_SIZE)

# MDS scan
K_MDS=40           # number of MDS axes to compute and scan
NPC_DIST=2         # number of PCs used to compute pairwise window distances
Z_THRESH=5         # MAD Z-score threshold for flagging outlier windows
WINDOW_GAP=20      # maximum gap (in windows) allowed within a candidate run; scales inversely with WINDOW_SIZE
MIN_WIN=4          # minimum number of outlier windows to retain a candidate
N_PERM=1000        # number of permutations for candidate significance test
PERM_P=0.01        # permutation p-value threshold

# Local PCA clustering
MIN_CLUSTER_N=5              # minimum samples per cluster for a valid three-class solution
MIN_SEARCH_N=3               # minimum samples in a cluster during gap-statistic search
MIN_GAP_SCALED_CUTOFF=0.05   # minimum scaled gap statistic for the best clustering axis
GAP_PRODUCT_SCALED_CUTOFF=0.003  # minimum product of gap statistics across both cut points
HET_P_CUTOFF=0.05            # p-value threshold for heterozygote excess test in middle cluster
HET_MIN_GROUP_N=5            # minimum group size required to run heterozygote test
PREFER_UNROTATED=TRUE        # prefer unrotated PC axes over 45-degree rotated axes when both pass
UNROTATED_MARGIN=0.95        # unrotated axis must score at least this fraction of rotated to be preferred
DIP_P_THRESH=0.05            # minimum dip test p-value (non-unimodality) for the clustering axis
R2_K3_MIN=0.85               # minimum R² for three-class genotype clustering quality
COR_H_AXIS_MAX=0.5           # maximum correlation between heterozygosity and the clustering axis
MIN_INV_FREQ=0.03            # minimum frequency of the rarer homozygote class

# LD breakpoint refinement
LD_BIN_SIZE=50000      # bin size for LD profile (bp); controls breakpoint resolution
COLLAPSE_BUFFER=0.5    # proportional buffer for second-collapse overlap check (fraction of haploblock length);
                       # allows nearby fragments of the same SV to merge across TE gaps
LD_MIN_SNR=10          # minimum sig_r2/bg_r2 to attempt LD extension; below this, MDS candidate
                       # coordinates are used — prevents over-extension on high-background LD chromosomes
