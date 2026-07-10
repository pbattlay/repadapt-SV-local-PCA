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

export TRUTH_SET="${TRUTH_SET:-}"

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

  # Canvas grows with candidate count; scale fonts + colourbar with it (no cap, so
  # per-tile resolution is preserved) — else labels/legend vanish at fit-to-page.
  pdf_dim <- max(6, n_cands * 0.22 + 3)
  sf      <- pdf_dim / 6

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
      na.value = "grey85",
      guide    = guide_colourbar(barwidth  = unit(0.35 * sf, "in"),
                                 barheight = unit(2.5  * sf, "in"))
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
    theme_minimal(base_size = 9 * sf) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8 * sf),
      axis.text.y     = element_text(size = 8 * sf),
      axis.title      = element_blank(),
      panel.grid      = element_blank(),
      legend.title    = element_text(size = 11 * sf),
      legend.text     = element_text(size = 8 * sf),
      legend.position = "right",
      plot.margin     = margin(5, 10, 5, 5)
    )

  out_hm  <- sprintf("final/%s_LD_heatmap.pdf", raid)
  ggsave(out_hm, ht, width = pdf_dim, height = pdf_dim, limitsize = FALSE)
  message(sprintf("[heatmap] → %s (%.1f × %.1f in)", out_hm, pdf_dim, pdf_dim))
}

# =============================================================
# SECTION 3b — Genome overview (candidates per chromosome)
# =============================================================
# One chromosome per row; per candidate, vertically + horizontally nested bars:
# merged (green) over LD/bp (red) over MDS (blue). Backbones span every analysis
# scaffold (contigs.tsv). Optional truth track (gold) from TRUTH_SET — chr,start,end,
# 1-based, no header. Runs for any candidate count (>=1).
if (!file.exists("ref/contigs.tsv")) {
  message("[overview] ref/contigs.tsv not found — skipping overview")
} else {
  contigs <- fread("ref/contigs.tsv", header = FALSE, col.names = c("chr", "length"))
  ov_chrs <- chrs[chrs %in% contigs$chr]                     # analysis scaffolds, chrs.txt order
  back    <- contigs[chr %in% ov_chrs]
  back[, chr := factor(chr, levels = ov_chrs)]
  setorder(back, chr)
  back[, cy := length(ov_chrs) - as.integer(chr) + 1L]       # first chr at top
  cy_map  <- setNames(back$cy, as.character(back$chr))

  mv <- copy(merged)
  mv[, cy := cy_map[as.character(chr)]]

  # Uniform-thickness lanes stacked per chromosome (like the MDS-scan plots) rather
  # than centered/nested — order top->bottom follows pipeline: truth, MDS, LD ext, merged.
  bh <- 0.06                                       # uniform bar half-height
  lanes <- rbindlist(list(
    mv[, .(cy, xmin = start,        xmax = end,        center = cy + 0.09, type = "MDS")],
    mv[, .(cy, xmin = bp_start,     xmax = bp_end,     center = cy - 0.09, type = "LD extension")],
    mv[, .(cy, xmin = merged_start, xmax = merged_end, center = cy - 0.27, type = "merged")]
  ), use.names = TRUE)

  # Optional truth track: neutral, not coloured — a black bar above plus a
  # translucent grey region shade (as in the parameter-sweep plots).
  tshade     <- NULL
  truth_path <- Sys.getenv("TRUTH_SET")
  if (nzchar(truth_path) && file.exists(truth_path) && file.info(truth_path)$size > 0) {
    truth_dt <- tryCatch(
      fread(truth_path, header = FALSE,
            colClasses = c("character", "integer", "integer"),
            col.names = c("chr", "start", "end")),
      error = function(e) { message("[overview] could not read TRUTH_SET (", e$message, ")"); NULL })
    if (!is.null(truth_dt)) {
      truth_dt <- truth_dt[chr %in% ov_chrs]
      if (nrow(truth_dt)) {
        truth_dt[, cy := cy_map[as.character(chr)]]
        lanes  <- rbindlist(list(lanes,
          truth_dt[, .(cy, xmin = start, xmax = end, center = cy + 0.27, type = "truth")]),
          use.names = TRUE)
        tshade <- truth_dt[, .(cy, xmin = start, xmax = end)]      # grey region highlight
        message(sprintf("[overview] %d truth region(s) overlaid", nrow(truth_dt)))
      }
    }
  }

  # Colours standardised with the MDS-scan plots; truth neutral (black).
  lane_levels <- c("truth", "MDS", "LD extension", "merged")
  pal         <- c(truth = "black", MDS = "steelblue",
                   `LD extension` = "firebrick", merged = "forestgreen")
  lanes[, type := factor(type, levels = lane_levels)]
  present     <- lane_levels[lane_levels %in% as.character(lanes$type)]

  ov <- ggplot() +
    geom_segment(data = back,
                 aes(x = 0, xend = length / 1e6, y = cy, yend = cy),
                 color = "grey80", linewidth = 0.6)
  if (!is.null(tshade))
    ov <- ov + geom_rect(data = tshade,
                         aes(xmin = xmin / 1e6, xmax = xmax / 1e6,
                             ymin = cy - 0.36, ymax = cy + 0.36),
                         fill = "grey70", alpha = 0.4)
  ov <- ov +
    geom_rect(data = lanes,
              aes(xmin = xmin / 1e6, xmax = xmax / 1e6,
                  ymin = center - bh, ymax = center + bh, fill = type)) +
    scale_fill_manual(values = pal, limits = lane_levels, breaks = present, name = NULL) +
    scale_y_continuous(breaks = back$cy, labels = as.character(back$chr),
                       expand = c(0, 0.6)) +
    scale_x_continuous(expand = c(0.02, 0)) +
    labs(x = "Position (Mb)", y = NULL,
         title    = sprintf("%s — candidates per chromosome", raid),
         subtitle = if (!is.null(tshade)) "truth = black bar + grey region" else NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "top")

  ov_h   <- max(4, length(ov_chrs) * 0.45 + 1.5)
  out_ov <- sprintf("final/%s_overview.pdf", raid)
  ggsave(out_ov, ov, width = 12, height = ov_h, limitsize = FALSE)
  message(sprintf("[overview] %d chromosomes → %s", length(ov_chrs), out_ov))
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
