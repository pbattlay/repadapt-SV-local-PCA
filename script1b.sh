#!/bin/bash
#SBATCH --job-name=script1b
#SBATCH --output=errout/%x.out
#SBATCH --error=errout/%x.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00

set -euo pipefail

source "$1"
export min_chr_length WINDOW_SIZE MIN_SNPS_PER_WINDOW raid

echo "[$(date)] Script 1b — ${raid}: diagnostic plots"

module load R

mkdir -p vcf/plots

# =============================================================
# STEP 1 — R plots (PCA done in R via BEDMatrix)
# =============================================================
echo "[$(date)]   Generating plots"

Rscript --vanilla - << 'REOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(BEDMatrix)
})

min_chr_length      <- as.integer(Sys.getenv("min_chr_length"))
window_size         <- as.integer(Sys.getenv("WINDOW_SIZE"))
min_snps_per_window <- as.integer(Sys.getenv("MIN_SNPS_PER_WINDOW"))
raid                <- Sys.getenv("raid")

theme_set(theme_minimal(base_size = 11))

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

# ── Plot 2: PCA via BEDMatrix + prcomp ──────────────────────────────────────
bim_files <- list.files("vcf", pattern = "\\.bim$", full.names = TRUE)
bed_files <- sub("\\.bim$", ".bed", bim_files)

if (length(bim_files) > 0) {
  # Count SNPs per chromosome for proportional sampling
  snp_counts <- sapply(bim_files, function(f)
    as.integer(system(sprintf("wc -l < '%s'", f), intern = TRUE)))
  total_snps <- sum(snp_counts)
  n_target   <- min(100000L, total_snps)
  n_per_chr  <- pmax(1L, round(snp_counts / total_snps * n_target))

  # Load sampled columns from each chromosome BED file
  set.seed(42)
  geno_list <- mapply(function(bed_f, n_snps, n_samp) {
    if (!file.exists(bed_f) || n_snps == 0L) return(NULL)
    bm   <- BEDMatrix(bed_f)
    cols <- sort(sample.int(ncol(bm), min(n_samp, ncol(bm))))
    as.matrix(bm[, cols, drop = FALSE])
  }, bed_files, snp_counts, n_per_chr, SIMPLIFY = FALSE)

  geno_mat    <- do.call(cbind, Filter(Negate(is.null), geno_list))
  n_samples   <- nrow(geno_mat)
  n_snps_used <- ncol(geno_mat)

  # Mean-impute any residual NAs (rare post-BEAGLE)
  na_cols <- which(colSums(is.na(geno_mat)) > 0L)
  for (j in na_cols) {
    m <- mean(geno_mat[, j], na.rm = TRUE)
    geno_mat[is.na(geno_mat[, j]), j] <- if (is.na(m)) 0 else m
  }

  pca_res <- prcomp(geno_mat, center = TRUE, scale. = FALSE)
  pve     <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
  pc_dt   <- data.table(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2])

  p2 <- ggplot(pc_dt, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.5, size = 0.9, color = "grey35") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("n = %d", n_samples), size = 3.5) +
    labs(x = sprintf("PC1 (%.1f%%)", pve[1]),
         y = sprintf("PC2 (%.1f%%)", pve[2]),
         title = sprintf("%s — PCA (%s SNPs)", raid,
                         format(n_snps_used, big.mark = ","))) +
    coord_fixed()

  ggsave("vcf/plots/pca.pdf", p2, width = 6, height = 6)
  message(sprintf("[plot] pca.pdf (%d SNPs, %d samples)", n_snps_used, n_samples))
} else {
  message("[plot] No BED files found — skipping pca.pdf")
}

# ── Plot 3: SNPs per window ───────────────────────────────────────────────────
bim_files_w <- list.files("vcf", pattern = "\\.bim$", full.names = TRUE)

if (length(bim_files_w) > 0) {
  bim_all <- rbindlist(lapply(bim_files_w, function(f) {
    dt <- fread(f, header = FALSE,
                col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
    dt[, .(chr, pos)]
  }))
  bim_all[, window_id := floor((pos - 1L) / window_size)]
  win_counts <- bim_all[, .N, by = .(chr, window_id)]
  prop_pass  <- mean(win_counts$N >= min_snps_per_window)

  p4 <- ggplot(win_counts, aes(x = N)) +
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

  ggsave("vcf/plots/snps_per_window.pdf", p4, width = 7, height = 5)
  message("[plot] snps_per_window.pdf")
} else {
  message("[plot] No .bim files found — skipping snps_per_window.pdf")
}

REOF

# =============================================================
# CLEANUP
# =============================================================
echo "[$(date)] Script 1b complete — ${raid}"
echo "  Plots written to vcf/plots/"
