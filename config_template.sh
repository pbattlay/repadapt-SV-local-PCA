### REPADAPT SV — LOCAL PCA - CONFIG
account="om62"
base_dir="/home/pbattlay/om62_scratch2/REPADAPT_LPCA"
raid="rawg0053"
vcf="/home/pbattlay/om62_scratch2/REPADAPT_LPCA/VCF/rawg0053_Ambrosia_artemisiifolia_Battlay.vcf.gz"

# Reference
min_chr_length=10000000    # minimum scaffold length to include (bp); filters out small scaffolds

# Optional inputs (presence = active; a script ignores any file it can't find).
# Defaults live under inputs/; to disable one, remove/empty the file (not the variable).
#   mappability_mask : 0-based BED keep-list of uniquely-mappable intervals, built by script0
#                      (per-reference, written next to the FASTA — point this at it). Removes
#                      collapsed-paralog SNPs the quality filters miss. ABSENT = no masking +
#                      a loud warning in script1a (you then keep paralog/balanced-het artefacts).
#   remove           : sample IIDs to exclude (one per line, no header), applied in script1b.
#   truth            : chr start end (1-based, no header) — marked on script3 overview plots ONLY.
#   manual_candidates: chr start end (1-based, no header) — extra regions injected into the
#                      script2 candidate set and tested for clustering like MDS candidates.
MAPPABILITY_MASK="${base_dir}/FASTA/rawg0053_Ambrosia_artemisiifolia_Battlay.mappable.bed"
REMOVE_SAMPLES="${base_dir}/inputs/${raid}_remove.tsv"
TRUTH_SET="${base_dir}/inputs/${raid}_truth.tsv"
MANUAL_CANDIDATES="${base_dir}/inputs/${raid}_manual_candidates.tsv"

# VCF filtering
MIN_QUAL=30         # minimum site QUAL score
MIN_GQ=20           # minimum sample genotype quality
MIN_DP=3            # minimum sample read depth
GENO=0.3            # maximum missingness per site (plink2 --geno)
MIND=0.9            # maximum missingness per individual (plink2 --mind); drops failed samples
MAF=0.01            # minimum minor allele frequency

# BEAGLE
BEAGLE_JAR="${base_dir}/scripts/beagle.jar"    # path to beagle.jar
BEAGLE_MEM=68                                  # JVM heap memory (GB)
BEAGLE_THREADS=4                               # imputation threads
BEAGLE_WINDOW=20.0                             # sliding window size for imputation (cM)

# Windowing
WINDOW_SIZE=10000          # local PCA window size (bp). SNP count per window varies with local
                           # density; the bp window gives a location-stable detectable-size floor
                           # (smallest SV ~= MIN_WIN x WINDOW_SIZE).
MIN_SNPS_PER_WINDOW=20     # minimum SNPs for a window to be included. NOTE: windows near this floor
                           # are rank-deficient when below the sample size (N individuals); set well
                           # above 20 (toward ~N) if your windows are SNP-poor.

# MDS scan
K_MDS=40           # number of MDS axes to compute and scan
NPC_DIST=2         # number of PCs used to compute pairwise window distances
Z_THRESH=5         # MAD Z-score threshold for flagging outlier windows
WINDOW_GAP=20      # maximum gap (in windows) allowed within a candidate run; scales inversely with WINDOW_SIZE
MIN_WIN=4          # minimum number of outlier windows to retain a candidate (smallest detectable
                   # SV ~= MIN_WIN x WINDOW_SIZE)
N_PERM=1000        # number of permutations for candidate significance test
PERM_P=0.01        # permutation p-value threshold

# Local PCA clustering
MIN_CLUSTER_N=5              # minimum samples per cluster for a valid three-class solution
MIN_SEARCH_N=3               # minimum samples in a cluster during gap-statistic search
MIN_GAP_SCALED_CUTOFF=0.07   # minimum scaled gap statistic for the best clustering axis
GAP_PRODUCT_SCALED_CUTOFF=0.007  # minimum product of gap statistics across both cut points
HET_P_CUTOFF=0.05            # p-value threshold for heterozygote excess test in middle cluster
HET_MIN_GROUP_N=5            # minimum group size required to run heterozygote test
DIP_P_THRESH=0.05            # minimum dip test p-value (non-unimodality) for the clustering axis
R2_K3_MIN=0.85               # minimum R² for three-class genotype clustering quality
COR_H_AXIS_MAX=0.75           # maximum correlation between heterozygosity and the clustering axis
MIN_INV_FREQ=0.03            # minimum frequency of the rarer homozygote class

# LD breakpoint refinement
# slated for removal – currently switched off with LD_MIN_SNR=10000
LD_BIN_SIZE=50000      # bin size for LD profile (bp); controls breakpoint resolution
LD_MIN_SNR=10000       # minimum sig_r2/bg_r2 to attempt LD extension; below this, MDS candidate
                       # coordinates are used — prevents over-extension on high-background LD chromosomes.
                       # Lowered 10 -> 5: extension-only, so this sharpens breakpoints with no new
                       # candidates and no background cost. Raise back toward 10 if edges over-extend.
LD_EXT_FRAC=0.25       # extension stop threshold as a fraction of the way from background to signal LD:
                       # walk out while smoothed r2 >= bg + LD_EXT_FRAC*(sig - bg). LOWER = reach further
                       # (risks walking into the LD halo past the true breakpoints; inflates SV size).

# Candidate collapse — first pass (Section C): merges spatially OVERLAPPING candidates
# whose inversion genotypes are correlated. The required overlap is its own safety, so a
# moderate correlation gate is enough here. This pass is the within-region workhorse.
R_THRESH=0.8           # minimum |genotype correlation| to merge two OVERLAPPING candidates (Section C only).
                       # Leave at 0.8 — its safety comes from the overlap co-requirement, not the threshold.

# Candidate collapse — long-range pass (Section D2): merges candidates that need NOT overlap,
# using a physical distance reach (COLLAPSE_BUFFER, bp) and a stricter correlation gate
# (R_THRESH_LONG). Each accepted merge is validated by re-clustering the union and is REVERTED
# to its constituent fragments if the union does not re-cluster as one inversion (script2 Section D2).
COLLAPSE_BUFFER=10000000    # max edge-to-edge gap (bp) for a long-range merge; Inf = chromosome-wide reach,
                       # no spatial limit (R_THRESH_LONG + the re-cluster guard carry all the safety).
                       # PHYSICAL distance, not a proportion — small fragments of a large block can now
                       # reach the block body regardless of their own size. Cap (e.g. 60000000) only if
                       # you want a hard spatial bound; real inversions span up to tens of Mb.
R_THRESH_LONG=0.95      # minimum |genotype correlation| for a LONG-RANGE merge (Section D2). Higher than
                       # R_THRESH because there is no overlap to lean on — the correlation must carry all
                       # the discrimination. 0.9 sits in the gap between within-inversion r (>=0.95) and
                       # distinct/supergene r (<=0.76) on sunflower. Set HIGHER per-species where background
                       # LD is elevated (e.g. capeweed reaches r~0.91 on some within-scaffold distinct pairs).
