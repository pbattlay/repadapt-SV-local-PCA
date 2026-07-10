#!/bin/bash
#SBATCH --job-name=script1a
#SBATCH --output=errout/%x.out
#SBATCH --error=errout/%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00

set -euo pipefail

source "$1"
echo "[$(date)] Script 1a — ${raid}: filter whole VCF + diagnostic plots"

module load plink/2.0-alpha R

mkdir -p vcf/plots

# Comma-separated list of analysis scaffolds (>= min_chr_length) for --chr
chr_list=$(paste -sd, ref/chrs.txt)

# Mappability mask (OPTIONAL; MAPPABILITY_MASK from config, built per-reference by
# script0). When present we remove paralog/repeat SNPs HERE — upstream of --mind/--geno
# and of BEAGLE imputation in script1c — so per-individual and per-site missingness and
# allele frequencies are never computed on artefactual sites, and imputation never smears
# them across the region. ABSENT = the pipeline still runs, but those collapsed-paralog
# SNPs (balanced-het artefacts that mimic inversion heterozygotes) stay in.
# (plink2 'bed0' = 0-based BED, the genmap/bedtools convention. If your plink2 build
#  rejects 'bed0', restrict the VCF with a bcftools --regions-file pre-pass instead.)
mask_args=()
if [[ -n "${MAPPABILITY_MASK:-}" && -s "${MAPPABILITY_MASK}" ]]; then
    # Restrict the (genome-wide) mask to the analysis chromosomes BEFORE plink reads it.
    # plink rejects any chromosome in an --extract file that isn't in the loaded dataset:
    # it aborts with "Invalid chromosome code" rather than ignoring it (--allow-extra-chr
    # only governs the main dataset). The genome-wide mask carries every FASTA contig —
    # including scaffolds below min_chr_length and chromosomes the VCF doesn't declare — so
    # plink dies on the first one outside its --chr set. Filtering to ref/chrs.txt (exactly
    # the chromosomes loaded here) drops those and shrinks the file. If nothing survives,
    # the mask's names don't match the VCF at all (e.g. an un-renamed 'lcl|...' FASTA).
    mask_use="vcf/mappable.analysis.bed"
    awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' ref/chrs.txt "${MAPPABILITY_MASK}" > "${mask_use}"
    n_mask=$(wc -l < "${mask_use}")
    if [[ "${n_mask}" -gt 0 ]]; then
        echo "[$(date)]   Mappability mask: ${MAPPABILITY_MASK} (${n_mask} intervals on analysis scaffolds)"
        mask_args=(--extract bed0 "${mask_use}")
    else
        {
          echo ""
          echo "  ############################################################"
          echo "  ## WARNING: mappability mask has NO intervals on the analysis"
          echo "  ## chromosomes (ref/chrs.txt). The chromosome names in the mask"
          echo "  ## don't match the VCF — e.g. FASTA 'lcl|Chr1' vs renamed VCF"
          echo "  ## 'lcl_Chr1'. Running WITHOUT masking. Fix: rename the FASTA"
          echo "  ## headers to match the VCF and regenerate the mask with script0"
          echo "  ## (or rename column 1 of the .bed)."
          echo "  ##   mask = ${MAPPABILITY_MASK}"
          echo "  ############################################################"
          echo ""
        } >&2
    fi
else
    {
      echo ""
      echo "  ############################################################"
      echo "  ## WARNING: no mappability mask."
      echo "  ##   MAPPABILITY_MASK = ${MAPPABILITY_MASK:-<unset>}"
      echo "  ## Running WITHOUT mappability masking — collapsed-paralog SNPs"
      echo "  ## (high-QUAL/GQ/DP balanced-het artefacts that mimic inversion"
      echo "  ## heterozygotes) will NOT be removed and can create false"
      echo "  ## local-PCA structure. Build one with script0 and point"
      echo "  ## MAPPABILITY_MASK at it."
      echo "  ############################################################"
      echo ""
    } >&2
fi

# =============================================================
# STEP 1: VCF FILTERING -> genome-wide BED
# =============================================================
# Filters all analysis scaffolds in one pass and writes a single plink1 BED.
# This BED is the source of truth for: the plots below, AND script1b (which
# subsets samples + splits per chromosome from it).
# NOTE: sample exclusion (REMOVE_SAMPLES) is NOT applied here — it happens in
# script1b after PCA review, so site missingness is computed on the full set.
echo "[$(date)]   Filtering VCF -> BED"
plink2 \
    --vcf "${vcf}" \
    --chr "${chr_list}" \
    --threads 4 \
    --set-missing-var-ids @:# \
    --const-fid 0 \
    --allow-extra-chr \
    ${mask_args[@]+"${mask_args[@]}"} \
    --var-min-qual ${MIN_QUAL} \
    --vcf-min-gq   ${MIN_GQ} \
    --vcf-min-dp   ${MIN_DP} \
    --geno          ${GENO} \
    --mind          ${MIND} \
    --max-alleles 2 \
    --maf           ${MAF} \
    --make-bed \
    --out "vcf/${raid}.filtered"

N_IND=$(wc -l < "vcf/${raid}.filtered.fam")
N_SNP=$(wc -l < "vcf/${raid}.filtered.bim")
echo "  After filtering: ${N_IND} individuals, ${N_SNP} SNPs"

# =============================================================
# STEP 2: DIAGNOSTIC PLOTS (R, reading the filtered BED)
# =============================================================
echo "[$(date)]   Generating plots"
export raid min_chr_length WINDOW_SIZE MIN_SNPS_PER_WINDOW

Rscript --vanilla - << 'REOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(BEDMatrix)
})

raid                <- Sys.getenv("raid")
min_chr_length      <- as.integer(Sys.getenv("min_chr_length"))
window_size         <- as.integer(Sys.getenv("WINDOW_SIZE"))
min_snps_per_window <- as.integer(Sys.getenv("MIN_SNPS_PER_WINDOW"))

theme_set(theme_minimal(base_size = 11))

bed_prefix <- file.path("vcf", paste0(raid, ".filtered"))
bim_path   <- paste0(bed_prefix, ".bim")

# ── Plot 1: Scaffold size distribution ──────────────────────────────────────
contigs       <- fread("ref/contigs.tsv", header = FALSE,
                       col.names = c("chr", "length"))
total_bp      <- sum(contigs$length)
retained_bp   <- sum(contigs$length[contigs$length >= min_chr_length])
prop_retained <- retained_bp / total_bp

p1 <- ggplot(contigs, aes(x = length / 1e6)) +
  geom_histogram(bins = 60, fill = "#4393C3", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = min_chr_length / 1e6,
             color = "#D6604D", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("%.1f%% of genome retained on scaffolds \u2265%.0f Mb",
                           prop_retained * 100, min_chr_length / 1e6),
           size = 3.5) +
  labs(x = "Scaffold length (Mb)", y = "Number of scaffolds",
       title = sprintf("%s — scaffold size distribution", raid))

ggsave("vcf/plots/scaffold_sizes.pdf", p1, width = 7, height = 5)
message("[plot] scaffold_sizes.pdf")

# ── Plot 2: SNPs per window (post-filtering, pre-imputation) ─────────────────
if (file.exists(bim_path)) {
  bim <- fread(bim_path, header = FALSE,
               col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
  bim[, window_id := floor((pos - 1L) / window_size)]
  win_counts <- bim[, .N, by = .(chr, window_id)]
  prop_pass  <- mean(win_counts$N >= min_snps_per_window)

  p2 <- ggplot(win_counts, aes(x = N)) +
    geom_histogram(bins = 60, fill = "#4393C3", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = min_snps_per_window,
               color = "#D6604D", linetype = "dashed", linewidth = 0.7) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("%.1f%% of windows pass (\u2265%d SNPs)",
                             prop_pass * 100, min_snps_per_window),
             size = 3.5) +
    labs(x = "SNPs per window", y = "Number of windows",
         title = sprintf("%s — post-filtering SNPs per %s kb window",
                         raid, format(window_size / 1000, big.mark = ",")))

  ggsave("vcf/plots/snps_per_window.pdf", p2, width = 7, height = 5)
  message("[plot] snps_per_window.pdf")
} else {
  message("[plot] No .bim found — skipping snps_per_window.pdf")
}

# ── Plot 3: PCA (100k random SNPs, pre-imputation) ───────────────────────────
if (file.exists(paste0(bed_prefix, ".bed"))) {
  bm          <- BEDMatrix(bed_prefix, simple_names = TRUE)
  samples     <- rownames(bm)
  n_snp_total <- ncol(bm)
  n_target    <- min(100000L, n_snp_total)

  set.seed(42)
  cols     <- sort(sample.int(n_snp_total, n_target))
  geno_mat <- as.matrix(bm[, cols, drop = FALSE])

  # Vectorised mean-imputation of missing genotypes (substantial pre-imputation)
  col_means <- colMeans(geno_mat, na.rm = TRUE)
  col_means[is.nan(col_means)] <- 0
  na_idx <- which(is.na(geno_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) geno_mat[na_idx] <- col_means[na_idx[, "col"]]

  # Drop any zero-variance columns before PCA (no NAs remain after imputation)
  col_var <- colMeans(geno_mat^2) - colMeans(geno_mat)^2
  keep    <- which(col_var > 0)
  if (length(keep) < ncol(geno_mat)) geno_mat <- geno_mat[, keep, drop = FALSE]
  n_snps_used <- ncol(geno_mat)

  pca_res <- prcomp(geno_mat, center = TRUE, scale. = FALSE)
  pve     <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

  pc_dt <- data.table(sample = samples,
                      PC1    = pca_res$x[, 1],
                      PC2    = pca_res$x[, 2])
  fwrite(pc_dt, "vcf/plots/pca_preimpute.tsv", sep = "\t")
  message("[table] pca_preimpute.tsv")

  p3 <- ggplot(pc_dt, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.6, size = 1.1, color = "grey35") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("n = %d", length(samples)), size = 3.5) +
    labs(x = sprintf("PC1 (%.1f%%)", pve[1]),
         y = sprintf("PC2 (%.1f%%)", pve[2]),
         title = sprintf("%s — PCA, pre-imputation (%s SNPs)", raid,
                         format(n_snps_used, big.mark = ","))) +
    coord_fixed()

  ggsave("vcf/plots/pca_preimpute.pdf", p3, width = 6, height = 6)
  message(sprintf("[plot] pca_preimpute.pdf (%d SNPs, %d samples)",
                  n_snps_used, length(samples)))
} else {
  message("[plot] No BED found — skipping PCA")
}
REOF

echo "[$(date)] Script 1a complete — ${raid}"
echo "  Filtered BED: vcf/${raid}.filtered.{bed,bim,fam}"
echo "  Plots:        vcf/plots/{scaffold_sizes,snps_per_window,pca_preimpute}.pdf"
echo "  PCA table:    vcf/plots/pca_preimpute.tsv"
echo "  --> Review the PCA, then list any samples to drop in your REMOVE_SAMPLES file before running script1b."
