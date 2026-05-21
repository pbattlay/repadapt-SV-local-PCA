#!/bin/bash
#SBATCH --job-name=script3
#SBATCH --output=errout/%x.out
#SBATCH --error=errout/%x.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00

set -euo pipefail

source "$1"
echo "[$(date)] Script 3 — ${raid}: merge, heatmap, finalise plots"

module load R

mkdir -p "final/plots"

Rscript --vanilla - "${raid}" << 'REOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

raid <- commandArgs(trailingOnly = TRUE)[1]
chrs <- readLines("ref/chrs.txt")

# =============================================================
# SECTION 1 — Merge breakpoints
# =============================================================
message("[merge] reading breakpoints")
bp_list <- lapply(chrs, function(chr) {
  f <- sprintf("breakpoints/%s.breakpoints.tsv", chr)
  if (file.exists(f)) fread(f) else NULL
})
merged <- rbindlist(Filter(Negate(is.null), bp_list), fill = TRUE)
if (nrow(merged) == 0L) stop("No breakpoint results found in breakpoints/")

# Sort by chromosome (chrs.txt order) then position
chr_levels  <- chrs[chrs %in% unique(merged$chr)]
merged[, chr := factor(chr, levels = chr_levels)]
setorder(merged, chr, merged_start)
merged[, cand_idx := .I]   # 1-based index in sorted order

out_bp <- sprintf("final/%s_merged.tsv", raid)
fwrite(merged, out_bp, sep = "\t", quote = FALSE)
message(sprintf("[merge] %d candidates → %s", nrow(merged), out_bp))

# =============================================================
# SECTION 2 — Merge genotypes
# =============================================================
message("[merge] reading genotypes")
geno_list <- lapply(chrs, function(chr) {
  f <- sprintf("breakpoints/%s.genotypes.tsv", chr)
  if (file.exists(f)) fread(f) else NULL
})
geno_list <- Filter(Negate(is.null), geno_list)
if (length(geno_list) == 0L) stop("No genotype files found in breakpoints/")

merged_geno <- Reduce(function(a, b) merge(a, b, by = "sample", all = TRUE),
                      geno_list)

# Reorder candidate columns to match sorted breakpoints order
geno_ids    <- intersect(merged$region_id, colnames(merged_geno))
merged_geno <- merged_geno[, c("sample", geno_ids), with = FALSE]

out_geno <- sprintf("final/%s_merged_genotypes.tsv", raid)
fwrite(merged_geno, out_geno, sep = "\t", quote = FALSE)
message(sprintf("[merge] %d samples × %d candidates → %s",
                nrow(merged_geno), length(geno_ids), out_geno))

# =============================================================
# SECTION 3 — Candidate × candidate LD heatmap
# =============================================================
n_cands <- nrow(merged)
if (n_cands < 2L) {
  message("[heatmap] fewer than 2 candidates — skipping")
} else {
  message(sprintf("[heatmap] computing %d × %d r² matrix", n_cands, n_cands))

  # Genotype matrix (samples × candidates, sorted column order)
  gmat    <- as.matrix(merged_geno[, geno_ids, with = FALSE])
  cor_mat <- cor(gmat, use = "pairwise.complete.obs")
  r2_mat  <- cor_mat ^ 2

  # Long format with positional indices
  r2_dt   <- as.data.table(as.table(r2_mat))
  setnames(r2_dt, c("rid_x", "rid_y", "r2"))
  idx_dt  <- data.table(region_id = merged$region_id, idx = merged$cand_idx)
  r2_dt   <- merge(r2_dt, idx_dt, by.x = "rid_x", by.y = "region_id")
  setnames(r2_dt, "idx", "x_idx")
  r2_dt   <- merge(r2_dt, idx_dt, by.x = "rid_y", by.y = "region_id")
  setnames(r2_dt, "idx", "y_idx")

  # Output LD table — upper triangle only (excluding diagonal), long format
  # Annotated with chr and position for downstream use
  pos_dt <- merged[, .(region_id, chr = as.character(chr),
                       bp_start = merged_start, bp_end = merged_end)]
  ld_out <- r2_dt[x_idx < y_idx][
    , .(rid_x, rid_y, r2 = round(r2, 6))]
  ld_out <- merge(ld_out, pos_dt, by.x = "rid_x", by.y = "region_id")
  setnames(ld_out, c("chr", "bp_start", "bp_end"), c("chr_x", "start_x", "end_x"))
  ld_out <- merge(ld_out, pos_dt, by.x = "rid_y", by.y = "region_id")
  setnames(ld_out, c("chr", "bp_start", "bp_end"), c("chr_y", "start_y", "end_y"))
  setcolorder(ld_out, c("rid_x", "chr_x", "start_x", "end_x",
                         "rid_y", "chr_y", "start_y", "end_y", "r2"))
  setorder(ld_out, chr_x, start_x, chr_y, start_y)
  out_ld <- sprintf("final/%s_LD_heatmap.tsv", raid)
  fwrite(ld_out, out_ld, sep = "\t", quote = FALSE)
  message(sprintf("[heatmap] %d candidate pairs → %s", nrow(ld_out), out_ld))

  # Chromosome annotation: start, end, midpoint index per chr
  chr_ann <- merged[, .(
    start_idx = min(cand_idx),
    end_idx   = max(cand_idx),
    mid_idx   = mean(cand_idx)
  ), keyby = chr]
  chr_ann[, chr_num  := .I]
  chr_ann[, fill_col := fifelse(chr_num %% 2L == 1L, "grey25", "grey65")]

  # Boundary lines between chromosomes
  chr_bounds <- chr_ann$end_idx[-nrow(chr_ann)] + 0.5

  # Chromosome annotation strip bounds (below and left of heatmap)
  s_lo <- -0.45
  s_hi <- -0.07

  ht <- ggplot(r2_dt, aes(x = x_idx, y = y_idx, fill = r2)) +

    # Chromosome annotation strips — annotate() bypasses aes() so no scale conflict
    annotate("rect",
             xmin  = chr_ann$start_idx - 0.5, xmax = chr_ann$end_idx + 0.5,
             ymin  = s_lo, ymax = s_hi,
             fill  = chr_ann$fill_col, color = NA) +
    annotate("rect",
             xmin  = s_lo, xmax = s_hi,
             ymin  = chr_ann$start_idx - 0.5, ymax = chr_ann$end_idx + 0.5,
             fill  = chr_ann$fill_col, color = NA) +

    # r² tiles
    geom_tile(color = "white", linewidth = 0.15) +

    # Chromosome boundary lines
    geom_vline(xintercept = chr_bounds, color = "grey55", linewidth = 0.25) +
    geom_hline(yintercept = chr_bounds, color = "grey55", linewidth = 0.25) +

    # Sequential colour scale (0–1 for r²)
    scale_fill_gradient(
      low      = "white", high = "#D6604D",
      limits   = c(0, 1), name = expression(r^2),
      breaks   = c(0, 0.25, 0.5, 0.75, 1),
      na.value = "grey85"
    ) +

    # Axes: chromosome labels at group midpoints
    scale_x_continuous(
      breaks = chr_ann$mid_idx, labels = as.character(chr_ann$chr),
      limits = c(s_lo - 0.02, n_cands + 0.5), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = chr_ann$mid_idx, labels = as.character(chr_ann$chr),
      limits = c(s_lo - 0.02, n_cands + 0.5), expand = c(0, 0)
    ) +

    coord_fixed() +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 8),
      axis.title      = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.margin     = margin(5, 10, 5, 5)
    )

  pdf_dim <- max(6, n_cands * 0.22 + 3)
  out_hm  <- sprintf("final/%s_LD_heatmap.pdf", raid)
  ggsave(out_hm, ht, width = pdf_dim, height = pdf_dim, limitsize = FALSE)
  message(sprintf("[heatmap] → %s (%.1f × %.1f in)", out_hm, pdf_dim, pdf_dim))
}

# =============================================================
# SECTION 4 — Copy plots to final/plots/
# =============================================================
message("[plots] copying to final/plots/")
n_ok  <- 0L
n_mis <- 0L
for (i in seq_len(nrow(merged))) {
  chr <- as.character(merged$chr[i])
  rid <- merged$region_id[i]
  pairs <- list(
    c(sprintf("breakpoints/plots/%s_%s.png",  chr, rid),
      sprintf("final/plots/%s_%s_bp.png",     chr, rid)),
    c(sprintf("local_pca/plots/%s_cluster.png", rid),
      sprintf("final/plots/%s_%s_lp.png",     chr, rid))
  )
  for (pr in pairs) {
    if (file.exists(pr[1])) {
      file.copy(pr[1], pr[2], overwrite = TRUE)
      n_ok <- n_ok + 1L
    } else {
      warning(sprintf("plot not found: %s", pr[1]))
      n_mis <- n_mis + 1L
    }
  }
}
message(sprintf("[plots] %d copied, %d not found", n_ok, n_mis))

REOF

echo "[$(date)] Script 3 complete — ${raid}"
