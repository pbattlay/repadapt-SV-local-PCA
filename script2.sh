#!/bin/bash
#SBATCH --job-name=script2
#SBATCH --output=errout/%x-%a.out
#SBATCH --error=errout/%x-%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00

set -euo pipefail

source "$1"
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ref/chrs.txt)
echo "[$(date)] Script 2 — ${chr}"

module load R/4.4.0-mkl plink/2.0-alpha

# Export config for R
export CHR="${chr}"
export WINDOW_SIZE MIN_SNPS_PER_WINDOW
export K_MDS NPC_DIST Z_THRESH WINDOW_GAP MIN_WIN N_PERM PERM_P
export MIN_CLUSTER_N MIN_SEARCH_N MIN_GAP_SCALED_CUTOFF GAP_PRODUCT_SCALED_CUTOFF
export HET_P_CUTOFF HET_MIN_GROUP_N
export DIP_P_THRESH R2_K3_MIN COR_H_AXIS_MAX MIN_INV_FREQ
export LD_BIN_SIZE COLLAPSE_BUFFER LD_MIN_SNR R_THRESH R_THRESH_LONG LD_EXT_FRAC
export MANUAL_CANDIDATES="${MANUAL_CANDIDATES:-}"

Rscript --vanilla - << 'RSCRIPT'

suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
  library(lostruct)
  library(ggplot2)
  library(diptest)
})

# =============================================================
# CONFIG
# =============================================================
chr                    <- Sys.getenv("CHR")
window_size            <- as.integer(Sys.getenv("WINDOW_SIZE"))
min_snps               <- as.integer(Sys.getenv("MIN_SNPS_PER_WINDOW"))
k_mds                  <- as.integer(Sys.getenv("K_MDS"))
npc_dist               <- as.integer(Sys.getenv("NPC_DIST"))
z_thresh               <- as.numeric(Sys.getenv("Z_THRESH"))
window_gap             <- as.integer(Sys.getenv("WINDOW_GAP"))
min_win                <- as.integer(Sys.getenv("MIN_WIN"))
n_perm                 <- as.integer(Sys.getenv("N_PERM"))
perm_p                 <- as.numeric(Sys.getenv("PERM_P"))
min_cluster_n          <- as.integer(Sys.getenv("MIN_CLUSTER_N"))
min_search_n           <- as.integer(Sys.getenv("MIN_SEARCH_N"))
min_gap_scaled_cutoff  <- as.numeric(Sys.getenv("MIN_GAP_SCALED_CUTOFF"))
gap_prod_scaled_cutoff <- as.numeric(Sys.getenv("GAP_PRODUCT_SCALED_CUTOFF"))
het_p_cutoff           <- as.numeric(Sys.getenv("HET_P_CUTOFF"))
het_min_group_n        <- as.integer(Sys.getenv("HET_MIN_GROUP_N"))
dip_p_thresh           <- as.numeric(Sys.getenv("DIP_P_THRESH"))
r2_k3_min              <- as.numeric(Sys.getenv("R2_K3_MIN"))
cor_h_axis_max         <- as.numeric(Sys.getenv("COR_H_AXIS_MAX"))
min_inv_freq           <- as.numeric(Sys.getenv("MIN_INV_FREQ"))
ld_bin_size            <- as.integer(Sys.getenv("LD_BIN_SIZE"))
ld_min_snr             <- as.numeric(Sys.getenv("LD_MIN_SNR", unset = "10"))
collapse_buffer        <- as.numeric(Sys.getenv("COLLAPSE_BUFFER", unset = "Inf")) # physical bp gap; Inf = chromosome-wide
r_thresh               <- as.numeric(Sys.getenv("R_THRESH", unset = "0.8"))        # Section C (overlap) collapse only
r_thresh_long          <- as.numeric(Sys.getenv("R_THRESH_LONG", unset = "0.9"))   # Section D2 long-range collapse
ld_ext_frac            <- as.numeric(Sys.getenv("LD_EXT_FRAC", unset = "0.25"))

ts <- function() format(Sys.time(), "%H:%M:%S")
cat(sprintf("[%s] R started — %s\n", ts(), chr))

# =============================================================
# FUNCTIONS: MDS SCAN
# =============================================================

detect_candidates <- function(mds, axes, z_thresh, window_gap, min_win, n_perm, perm_p) {
  N       <- nrow(mds)
  all_idx <- seq_len(N)

  longest_run <- function(sorted_idx) {
    if (length(sorted_idx) < 1L) return(0L)
    gaps <- diff(sorted_idx) - 1L
    max(tabulate(cumsum(c(1L, gaps > window_gap))))
  }

  results <- list()

  for (ax in axes) {
    if (!ax %in% colnames(mds)) next
    vals <- mds[[ax]]
    med  <- median(vals, na.rm = TRUE)
    mad  <- median(abs(vals - med), na.rm = TRUE)
    if (!is.finite(mad) || mad == 0) next
    z <- (vals - med) / (mad * 1.4826)

    for (dir in c("pos", "neg")) {
      flag_idx <- sort(if (dir == "pos") which(z > z_thresh) else which(z < -z_thresh))
      K <- length(flag_idx)
      if (K < min_win) next

      L_obs  <- longest_run(flag_idx)
      L_perm <- replicate(n_perm, longest_run(sort(sample(all_idx, K))))
      p_val  <- mean(L_perm >= L_obs)
      if (p_val >= perm_p) next

      # Per-run length floor from the same permutation null: a run must be longer
      # than the longest streak noise produces at this window_gap / K / N. This
      # tracks window_gap automatically (no need to retune min_win when the gap
      # changes); min_win remains the hard lower bound for the SV size floor.
      run_floor <- max(min_win, as.integer(ceiling(quantile(L_perm, 1 - perm_p))))

      gaps   <- diff(flag_idx) - 1L
      run_id <- cumsum(c(1L, gaps > window_gap))

      for (rid in unique(run_id)) {
        sel <- run_id == rid
        if (sum(sel) < run_floor) next
        idx   <- flag_idx[sel]
        z_run <- abs(z[idx])
        results[[length(results) + 1L]] <- data.frame(
          chr             = mds$chr[idx[1]],
          start           = min(mds$window_start[idx]),
          end             = max(mds$window_end[idx]),
          axis            = ax,
          dir             = dir,
          n_windows       = sum(sel),
          K               = K,
          run_floor       = run_floor,
          peak_z          = round(max(z_run),  2),
          mean_z          = round(mean(z_run), 2),
          sd_z            = round(sd(z_run),   2),
          max_gap_windows = if (length(idx) > 1L) max(diff(idx) - 1L) else 0L,
          perm_p          = round(p_val, 4)
        )
      }
    }
  }

  if (!length(results)) return(NULL)
  out <- rbindlist(results)
  out[order(-peak_z)]
}

# =============================================================
# FUNCTIONS: LOCAL PCA
# =============================================================

safe_wilcox_p <- function(x, y, alternative = "greater") {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  tryCatch(wilcox.test(x, y, alternative = alternative, exact = FALSE)$p.value,
           error = function(e) NA_real_)
}

make_lpca_table <- function(raw_file) {
  raw <- as.data.table(fread(raw_file))
  if (ncol(raw) <= 6) stop("PLINK .raw file has no genotype columns")
  gt_dt  <- raw[, -(1:6), with = FALSE]
  gt_mat <- as.matrix(gt_dt)
  storage.mode(gt_mat) <- "numeric"
  H <- rowMeans(gt_mat == 1, na.rm = TRUE)
  variant_ok <- apply(gt_mat, 2, function(v) {
    v <- v[is.finite(v)]
    length(v) >= 2 && is.finite(var(v)) && var(v) > 0
  })
  gt_mat <- gt_mat[, variant_ok, drop = FALSE]
  if (ncol(gt_mat) < 2) stop("fewer than two variable genotype columns after filtering")
  col_means  <- colMeans(gt_mat, na.rm = TRUE)
  missing_idx <- which(!is.finite(gt_mat), arr.ind = TRUE)
  if (nrow(missing_idx) > 0) gt_mat[missing_idx] <- col_means[missing_idx[, 2]]
  covmat <- cov(t(gt_mat))
  if (any(!is.finite(covmat))) stop("covariance matrix contains non-finite values")
  eig <- eigen(covmat, symmetric = TRUE)
  pcs <- as.data.table(eig$vectors[, 1:2, drop = FALSE])
  setnames(pcs, c("PC1", "PC2"))
  pcs[, PC1_rot45 := (PC1 - PC2) / sqrt(2)]
  pcs[, PC2_rot45 := (PC1 + PC2) / sqrt(2)]
  cbind(data.table(FID = raw[[1]], IID = raw[[2]], H = H), pcs)
}

empty_axis_score <- function(axis, y, n_total) {
  total_range <- suppressWarnings(diff(range(y, na.rm = TRUE)))
  if (!is.finite(total_range)) total_range <- NA_real_
  data.table(
    axis = axis, valid = FALSE,
    cut1 = NA_real_, cut2 = NA_real_,
    n_L = NA_integer_, n_MID = NA_integer_, n_R = NA_integer_,
    gap_L_MID = NA_real_, gap_MID_R = NA_real_,
    min_gap = NA_real_, gap_product = NA_real_,
    total_range = total_range,
    min_gap_scaled = NA_real_, gap_product_scaled = NA_real_,
    H_L = NA_real_, H_MID = NA_real_, H_R = NA_real_,
    R2_k3 = NA_real_, dip_stat = NA_real_, dip_p = NA_real_,
    cor_H_axis_L = NA_real_, cor_H_axis_R = NA_real_,
    p_inv = NA_real_, maf = NA_real_, F_IS = NA_real_,
    gap_pass = FALSE, quality_pass = FALSE, overall_pass = FALSE
  )
}

score_axis_by_gaps <- function(axis, df,
                               min_cluster_n         = 5,
                               min_search_n          = 3,
                               min_gap_scaled_cutoff = 0.05,
                               gap_product_scaled_cutoff = 0.003,
                               dip_p_thresh          = 0.05,
                               r2_k3_min             = 0.85,
                               cor_h_axis_max        = 0.5,
                               min_inv_freq          = 0.03) {
  y  <- df[[axis]]; H <- df$H
  ok <- is.finite(y); n <- sum(ok)
  if (n < 3 * min_search_n) return(empty_axis_score(axis, y, n))
  vals <- sort(y[ok]); gaps <- diff(vals); total_range <- diff(range(vals))
  if (!is.finite(total_range) || total_range <= 0 || length(gaps) < 2)
    return(empty_axis_score(axis, y, n))
  gap_pos    <- seq_along(gaps)
  candidates <- CJ(gap1 = gap_pos, gap2 = gap_pos)[gap1 < gap2]
  candidates[, `:=`(n_L = gap1, n_MID = gap2 - gap1, n_R = n - gap2)]
  candidates <- candidates[n_L >= min_search_n & n_MID >= min_search_n & n_R >= min_search_n]
  if (!nrow(candidates)) return(empty_axis_score(axis, y, n))
  candidates[, `:=`(gap_L_MID = gaps[gap1], gap_MID_R = gaps[gap2])]
  candidates[, `:=`(min_gap = pmin(gap_L_MID, gap_MID_R),
                    gap_product = gap_L_MID * gap_MID_R)]
  candidates[, `:=`(min_gap_scaled     = min_gap / total_range,
                    gap_product_scaled = gap_product / total_range^2,
                    gap_sum            = gap_L_MID + gap_MID_R)]
  setorder(candidates, -min_gap_scaled, -gap_product_scaled, -gap_sum)
  best <- candidates[1]
  cut1    <- mean(c(vals[best$gap1], vals[best$gap1 + 1]))
  cut2    <- mean(c(vals[best$gap2], vals[best$gap2 + 1]))
  cluster <- rep(NA_character_, length(y))
  cluster[ok & y <= cut1]            <- "L"
  cluster[ok & y > cut1 & y <= cut2] <- "MID"
  cluster[ok & y > cut2]             <- "R"
  n_L   <- sum(cluster == "L",   na.rm = TRUE)
  n_MID <- sum(cluster == "MID", na.rm = TRUE)
  n_R   <- sum(cluster == "R",   na.rm = TRUE)
  n_tot <- n_L + n_MID + n_R
  grand_mean <- mean(y[ok], na.rm = TRUE)
  total_ss   <- sum((y[ok] - grand_mean)^2, na.rm = TRUE)
  ns_3  <- c(n_L, n_MID, n_R)
  mus_3 <- c(mean(y[cluster == "L"],   na.rm = TRUE),
             mean(y[cluster == "MID"], na.rm = TRUE),
             mean(y[cluster == "R"],   na.rm = TRUE))
  R2_k3 <- if (total_ss > 0) sum(ns_3 * (mus_3 - grand_mean)^2) / total_ss else NA_real_
  safe_cor <- function(hv, av) {
    ok2 <- is.finite(hv) & is.finite(av)
    if (sum(ok2) < 4L) return(NA_real_)
    cor(hv[ok2], av[ok2], method = "spearman")
  }
  cor_H_L   <- safe_cor(H[cluster == "L"], y[cluster == "L"])
  cor_H_R   <- safe_cor(H[cluster == "R"], y[cluster == "R"])
  max_cor_H <- suppressWarnings(max(abs(c(cor_H_L, cor_H_R)), na.rm = TRUE))
  if (!is.finite(max_cor_H)) max_cor_H <- NA_real_
  dip_res  <- tryCatch(dip.test(y[ok], simulate.p.value = FALSE), error = function(e) NULL)
  dip_stat <- if (!is.null(dip_res)) as.numeric(dip_res$statistic) else NA_real_
  dip_p    <- if (!is.null(dip_res)) dip_res$p.value               else NA_real_
  n_min_hom <- min(n_L, n_R)
  p_inv     <- (2 * n_min_hom + n_MID) / (2 * n_tot)
  exp_het   <- 2 * p_inv * (1 - p_inv)
  F_IS      <- if (exp_het > 0) 1 - (n_MID / n_tot) / exp_het else NA_real_
  gap_pass     <- n_L >= min_cluster_n & n_MID >= min_cluster_n & n_R >= min_cluster_n &
                  isTRUE(best$min_gap_scaled >= min_gap_scaled_cutoff) &
                  isTRUE(best$gap_product_scaled >= gap_product_scaled_cutoff)
  quality_pass <- isTRUE(!is.na(dip_p)     && dip_p     <  dip_p_thresh)  &
                  isTRUE(!is.na(R2_k3)     && R2_k3     >= r2_k3_min)     &
                  isTRUE(!is.na(max_cor_H) && max_cor_H <  cor_h_axis_max) &
                  isTRUE(!is.na(p_inv)     && p_inv     >= min_inv_freq)
  data.table(
    axis = axis, valid = TRUE, cut1 = cut1, cut2 = cut2,
    n_L = n_L, n_MID = n_MID, n_R = n_R,
    gap_L_MID = best$gap_L_MID, gap_MID_R = best$gap_MID_R,
    min_gap = best$min_gap, gap_product = best$gap_product,
    total_range = total_range,
    min_gap_scaled = best$min_gap_scaled, gap_product_scaled = best$gap_product_scaled,
    H_L = mean(H[cluster == "L"],   na.rm = TRUE),
    H_MID = mean(H[cluster == "MID"], na.rm = TRUE),
    H_R = mean(H[cluster == "R"],   na.rm = TRUE),
    R2_k3 = round(R2_k3, 4), dip_stat = round(dip_stat, 6), dip_p = round(dip_p, 6),
    cor_H_axis_L = round(cor_H_L, 4), cor_H_axis_R = round(cor_H_R, 4),
    p_inv = round(p_inv, 4), maf = round(min(p_inv, 1 - p_inv), 4),
    F_IS = round(F_IS, 4),
    gap_pass = gap_pass, quality_pass = quality_pass,
    overall_pass = gap_pass & quality_pass
  )
}

cluster_local_pca_gaps <- function(df,
                                   min_cluster_n         = 5,
                                   min_search_n          = 3,
                                   min_gap_scaled_cutoff = 0.05,
                                   gap_product_scaled_cutoff = 0.003,
                                   dip_p_thresh          = 0.05,
                                   r2_k3_min             = 0.85,
                                   cor_h_axis_max        = 0.5,
                                   min_inv_freq          = 0.03) {
  df <- as.data.table(copy(df))
  df[, PC1_rot45 := (PC1 - PC2) / sqrt(2)]
  df[, PC2_rot45 := (PC1 + PC2) / sqrt(2)]
  axes <- c("PC1", "PC2", "PC1_rot45", "PC2_rot45")
  axis_scores <- rbindlist(lapply(axes, score_axis_by_gaps,
    df = df,
    min_cluster_n         = min_cluster_n,
    min_search_n          = min_search_n,
    min_gap_scaled_cutoff = min_gap_scaled_cutoff,
    gap_product_scaled_cutoff = gap_product_scaled_cutoff,
    dip_p_thresh          = dip_p_thresh,
    r2_k3_min             = r2_k3_min,
    cor_h_axis_max        = cor_h_axis_max,
    min_inv_freq          = min_inv_freq
  ), fill = TRUE)

  passing <- axis_scores[overall_pass == TRUE]
  if (!nrow(passing)) {
    fallback <- axis_scores[gap_pass == TRUE]
    if (!nrow(fallback)) fallback <- axis_scores[valid == TRUE]
    if (!nrow(fallback)) stop("no valid 3-cluster split found on any axis")
    setorder(fallback, -min_gap_scaled, -gap_product_scaled)
    best <- fallback[1]
  } else {
    # Highest gap-score among all four axes (PC1, PC2 and their 45-deg rotations);
    # no preference for unrotated axes.
    setorder(passing, -min_gap_scaled, -gap_product_scaled); best <- passing[1]
  }

  ax      <- best$axis
  y       <- df[[ax]]
  cluster <- rep(NA_character_, nrow(df))
  cluster[is.finite(y) & y <= best$cut1]            <- "L"
  cluster[is.finite(y) & y > best$cut1 & y <= best$cut2] <- "MID"
  cluster[is.finite(y) & y > best$cut2]             <- "R"
  df[, cluster := cluster]

  list(
    selected_axis  = ax,
    selected_score = best,
    axis_scores    = axis_scores,
    bounds         = list(L_start = -Inf, L_end = best$cut1,
                          MID_start = best$cut1, MID_end = best$cut2,
                          R_start = best$cut2, R_end = Inf),
    data           = df
  )
}

classify_local_pca_cluster <- function(res) {
  s <- res$selected_score
  data.table(
    cluster_pass   = isTRUE(s$overall_pass),
    cluster_reason = if (isTRUE(s$overall_pass)) "passed all filters" else
                     if (!isTRUE(s$gap_pass))    "gap metrics below threshold" else
                     "failed quality filters"
  )
}

test_heterozygosity <- function(df, p_cutoff = 0.05, min_group_n = 5) {
  H_L   <- df[cluster == "L",   H]
  H_MID <- df[cluster == "MID", H]
  H_R   <- df[cluster == "R",   H]
  H_HOM <- df[cluster %in% c("L", "R"), H]
  n_L   <- sum(is.finite(H_L)); n_MID <- sum(is.finite(H_MID)); n_R <- sum(is.finite(H_R))
  mean_L   <- mean(H_L,   na.rm = TRUE)
  mean_MID <- mean(H_MID, na.rm = TRUE)
  mean_R   <- mean(H_R,   na.rm = TRUE)
  p_mid_gt_L    <- safe_wilcox_p(H_MID, H_L,   "greater")
  p_mid_gt_R    <- safe_wilcox_p(H_MID, H_R,   "greater")
  p_mid_gt_homs <- safe_wilcox_p(H_MID, H_HOM, "greater")
  enough_n <- n_L >= min_group_n && n_MID >= min_group_n && n_R >= min_group_n
  pass <- enough_n &&
    is.finite(p_mid_gt_L) && is.finite(p_mid_gt_R) &&
    mean_MID > mean_L && mean_MID > mean_R &&
    p_mid_gt_L <= p_cutoff && p_mid_gt_R <= p_cutoff
  data.table(
    het_pass      = pass,
    p_mid_gt_L    = round(p_mid_gt_L,    4),
    p_mid_gt_R    = round(p_mid_gt_R,    4),
    p_mid_gt_homs = round(p_mid_gt_homs, 4),
    mean_H_L      = round(mean_L,   4),
    mean_H_MID    = round(mean_MID, 4),
    mean_H_R      = round(mean_R,   4)
  )
}

make_failure_row <- function(cand_row, region_id, reason) {
  message("  FAILED: ", reason)
  cbind(cand_row, data.table(
    region_id          = region_id,
    lpca_axis          = NA_character_,
    L_end              = NA_real_,
    MID_end            = NA_real_,
    cluster_pass       = FALSE,
    het_pass           = FALSE,
    min_gap_scaled     = NA_real_,
    gap_product_scaled = NA_real_,
    n_L = NA_integer_, n_MID = NA_integer_, n_R = NA_integer_,
    mean_H_L = NA_real_, mean_H_MID = NA_real_, mean_H_R = NA_real_,
    p_mid_gt_L = NA_real_, p_mid_gt_R = NA_real_, p_mid_gt_homs = NA_real_,
    R2_k3 = NA_real_, dip_p = NA_real_,
    cor_H_axis_L = NA_real_, cor_H_axis_R = NA_real_,
    maf = NA_real_, F_IS = NA_real_
  ))
}

plot_cluster <- function(res, cc, hc, region_id, png_file) {
  df     <- res$data
  ax     <- res$selected_axis
  sc     <- res$selected_score
  bd     <- res$bounds
  label  <- paste0(
    "cluster: ", ifelse(cc$cluster_pass, "PASS", "FAIL"),
    "\nhet: ",   ifelse(hc$het_pass,     "PASS", "FAIL"),
    "\naxis: ",  ax,
    "\nmin_gap_scaled: ",     signif(sc$min_gap_scaled,     3),
    "\ngap_product_scaled: ", signif(sc$gap_product_scaled, 3),
    "\nn: L=", sc$n_L, " MID=", sc$n_MID, " R=", sc$n_R
  )
  p <- ggplot(df, aes(x = .data[[ax]], y = H, colour = cluster)) +
    geom_point(size = 1.5, alpha = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = c(bd$L_end, bd$MID_end), linetype = "dashed") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.05,
             label = label, size = 2.8) +
    theme_classic() +
    labs(title = region_id, x = ax, y = "Heterozygosity", colour = "Cluster")
  ggsave(png_file, p, width = 6.5, height = 5.2, dpi = 200)
}

# =============================================================
# FUNCTIONS: COLLAPSE
# =============================================================

assign_genotypes <- function(region_id, lpca_axis, L_end, MID_end) {
  f     <- file.path("local_pca/tables", paste0(region_id, "_lpca.tsv"))
  lpca  <- fread(f)
  scores <- lpca[[lpca_axis]]
  geno  <- ifelse(scores < L_end, 0L, ifelse(scores < MID_end, 1L, 2L))
  setNames(as.integer(geno), lpca$IID)
}

collapse_candidates <- function(cands, geno_mat, r_thresh = 0.8) {
  # Representative = WIDEST candidate (physical span), source-agnostic so a manual
  # region's full extent wins over a narrower correlated MDS fragment; ties prefer
  # the MDS candidate (real axis/stats). This is the only place a non-representative's
  # coordinates are dropped (D2 unions extents), so span here is necessary + sufficient.
  src_rank <- if ("candidate_source" %in% names(cands)) (cands$candidate_source != "mds") else logical(nrow(cands))
  ord      <- order(cands$end - cands$start, src_rank, decreasing = c(TRUE, FALSE), method = "radix")
  cands    <- cands[ord]; geno_mat <- geno_mat[ord, , drop = FALSE]
  group    <- rep(NA_integer_, nrow(cands)); gid <- 1L
  for (i in seq_len(nrow(cands))) {
    if (!is.na(group[i])) next
    group[i] <- gid
    for (j in seq_len(nrow(cands))) {
      if (i == j || !is.na(group[j])) next
      if (cands$end[i] < cands$start[j] || cands$start[i] > cands$end[j]) next
      r <- cor(geno_mat[i, ], geno_mat[j, ], use = "pairwise.complete.obs")
      if (!is.na(r) && abs(r) >= r_thresh) group[j] <- gid
    }
    gid <- gid + 1L
  }
  rbindlist(lapply(unique(group), function(g) {
    sub  <- cands[group == g]; best <- sub[1]
    best[, supporting_axes := paste(unique(sub$axis), collapse = ",")]
    best[, n_collapsed      := nrow(sub)]
    best
  }))[order(-(end - start))]
}

collapse_candidates_bp <- function(bp_table, geno_wide, r_thresh = 0.9,
                                   max_gap = Inf) {
  # Sort by R2_k3: best-clustering representative chosen when merging
  ord   <- order(bp_table$R2_k3, decreasing = TRUE, na.last = TRUE)
  cands <- bp_table[ord]
  geno_mat <- do.call(rbind, lapply(cands$region_id, function(rid)
    as.numeric(geno_wide[[rid]])))
  rownames(geno_mat) <- cands$region_id

  n <- nrow(cands)

  # ---- Single-linkage grouping (union-find) ------------------------------------
  # Replaces the earlier seed/star grouping, which measured the gap from one seed to
  # every other candidate and assigned in a single hop. That fragmented large real
  # inversions two ways: a block wider than max_gap could not be held together by a
  # seed sitting off-centre, and a bridge fragment was consumed by whichever seed
  # reached it first, so an r=1.0 block could fracture at its midpoint.
  #
  # Here every PAIR is linked iff its edge-to-edge physical gap is <= max_gap AND
  # |genotype r| >= r_thresh; groups are the connected components of that graph.
  # Linkage is transitive, so a chain of consecutive linked fragments merges as one
  # block even when the outermost members are far apart and decorrelated-at-distance.
  # A single inversion broken into pieces — or split by an internal detection hole
  # (centromere/repeat stretch with no callable window) narrower than max_gap —
  # re-assembles. The r gate halts a chain at the inversion boundary, where the next
  # locus decorrelates. max_gap therefore only needs to bridge the largest WITHIN-
  # block hole, not span the block diameter; distinct-but-near loci are held apart by
  # the r gate, and the Section-D2 re-cluster/trim guard backstops any over-reach.
  # Candidates with NA bp coords never link, so they remain singletons (as before).
  parent <- seq_len(n)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  for (i in seq_len(n - 1L)) {
    bi <- cands[i]
    if (is.na(bi$bp_start) || is.na(bi$bp_end)) next
    for (j in seq(i + 1L, n)) {
      bj <- cands[j]
      if (is.na(bj$bp_start) || is.na(bj$bp_end)) next
      # edge-to-edge gap (0 when the intervals overlap); physical bp, not proportional
      gap <- max(0, bj$bp_start - bi$bp_end, bi$bp_start - bj$bp_end)
      if (gap > max_gap) next
      r <- cor(geno_mat[i, ], geno_mat[j, ], use = "pairwise.complete.obs")
      if (!is.na(r) && abs(r) >= r_thresh) {
        ri <- find(i); rj <- find(j)
        if (ri != rj) parent[rj] <- ri   # cands is R2_k3-sorted; lower index keeps the root
      }
    }
  }
  group <- vapply(seq_len(n), find, integer(1))

  # Per-fragment annotation: merged_start/end = union of group
  merge_map <- rbindlist(lapply(unique(group), function(g) {
    sub <- cands[group == g]
    data.table(region_id    = sub$region_id,
               merged_start = min(sub$bp_start, na.rm = TRUE),
               merged_end   = max(sub$bp_end,   na.rm = TRUE),
               n_merged     = nrow(sub))
  }))

  # Collapsed: one row per group, best R2_k3, union bp coords
  collapsed2 <- rbindlist(lapply(unique(group), function(g) {
    sub  <- cands[group == g]; best <- copy(sub[1])
    best[, merged_start   := min(sub$bp_start, na.rm = TRUE)]
    best[, merged_end     := max(sub$bp_end,   na.rm = TRUE)]
    best[, n_merged       := nrow(sub)]
    best[, merged_regions := paste(sub$region_id, collapse = ",")]
    all_axes <- unique(unlist(strsplit(paste(sub$supporting_axes, collapse = ","), ",")))
    best[, supporting_axes := paste(all_axes, collapse = ",")]
    best
  }))[order(-R2_k3)]

  list(collapsed = collapsed2, merge_map = merge_map)
}

# =============================================================
# FUNCTIONS: LD BREAKPOINT REFINEMENT
# =============================================================

compute_ld_chr <- function(geno_c, sd_g, samp, inv_geno, bim, bin_size) {
  ig      <- inv_geno[samp]
  ok_samp <- is.finite(ig)
  ig_c    <- ig[ok_samp] - mean(ig[ok_samp])
  sd_ig   <- sqrt(sum(ig_c^2))
  if (sd_ig == 0) return(NULL)
  # Vectorised r² via matrix-vector multiply (MKL-accelerated)
  num <- as.numeric(ig_c %*% geno_c[ok_samp, ])
  r2  <- ifelse(sd_g > 1e-10, (num / (sd_ig * sd_g))^2, NA_real_)
  bins    <- seq(1L, max(bim$pos, na.rm = TRUE) + bin_size, by = bin_size)
  bin_idx <- findInterval(bim$pos, bins)
  bin_mid <- (bins[-length(bins)] + bins[-1]) / 2
  r2_by_bin <- tapply(r2, bin_idx, mean, na.rm = TRUE)
  valid   <- as.integer(names(r2_by_bin))
  data.table(pos = bin_mid[valid], mean_r2 = round(as.numeric(r2_by_bin), 4))
}

find_breakpoints <- function(profile, cand_start, cand_end, ld_min_snr = 10,
                             ld_ext_frac = 0.25) {
  ord <- order(profile$pos)
  pos <- profile$pos[ord]; r2v <- profile$mean_r2[ord]; n <- length(r2v)
  if (sum(!is.na(r2v)) < 3) return(list(bp_start = cand_start, bp_end = cand_end,
                                          bg_r2 = NA_real_, sig_r2 = NA_real_))

  # Smooth with NA pass-through
  smoothed <- as.numeric(stats::filter(r2v, rep(1/3, 3), sides = 2))
  smoothed[is.na(smoothed)] <- r2v[is.na(smoothed)]

  # Background: 10th percentile of full scan — robust when inversion extends
  # through much of the buffer zone (large inversions contaminate flanking estimate)
  bg_r2  <- as.numeric(quantile(r2v, 0.10, na.rm = TRUE))
  inside <- pos >= cand_start & pos <= cand_end
  sig_r2 <- if (any(inside & !is.na(r2v))) median(r2v[inside], na.rm = TRUE) else
             max(r2v, na.rm = TRUE)

  # SNR guard: if the LD signal within the MDS candidate is not meaningfully
  # above chromosome background, the genotype assignments are too noisy to
  # trust for extension — fall back to MDS candidate coordinates
  if (!is.na(sig_r2) && !is.na(bg_r2) && bg_r2 > 0 &&
      sig_r2 / bg_r2 < ld_min_snr) {
    return(list(bp_start = cand_start, bp_end = cand_end,
                bg_r2 = round(bg_r2, 4), sig_r2 = round(sig_r2, 4)))
  }

  thresh <- bg_r2 + (sig_r2 - bg_r2) * ld_ext_frac

  # Anchor scan at the peak smoothed r2 within the MDS candidate — more robust
  # than the geometric centre when the candidate covers only part of an inversion
  inside_idx <- which(inside & !is.na(smoothed))
  ci <- if (length(inside_idx) > 0L) inside_idx[which.max(smoothed[inside_idx])] else
        which.min(abs(pos - (cand_start + cand_end) / 2L))

  # Scan left — NA bins are skipped unconditionally (assumed to be TE/repeat
  # regions within the SV). Only a non-NA bin below threshold stops the scan.
  # Track last_sig: the outermost above-threshold bin seen so far.
  li <- ci; last_sig_l <- ci
  while (li > 1L) {
    li_new <- li - 1L
    val    <- smoothed[li_new]
    if (is.na(val))    { li <- li_new                              # skip NA
    } else if (val >= thresh) { li <- li_new; last_sig_l <- li_new # extend
    } else             { break }                                    # stop
  }
  li <- last_sig_l

  # Scan right
  ri <- ci; last_sig_r <- ci
  while (ri < n) {
    ri_new <- ri + 1L
    val    <- smoothed[ri_new]
    if (is.na(val))    { ri <- ri_new
    } else if (val >= thresh) { ri <- ri_new; last_sig_r <- ri_new
    } else             { break }
  }
  ri <- last_sig_r

  # Extension-only: LD refinement can extend beyond the MDS candidate but never
  # shrink it. The MDS candidate defines a lower bound on the SV extent.
  bp_start <- min(pos[li], cand_start)
  bp_end   <- max(pos[ri], cand_end)

  list(bp_start = bp_start, bp_end = bp_end,
       bg_r2 = round(bg_r2, 4), sig_r2 = round(sig_r2, 4))
}

plot_mds_with_breakpoints <- function(mds_dt, chr, frags, merged_start, merged_end,
                                      n_merged, axis, dir, out_file) {
  # frags: data.table with one row per collapsed fragment (start,end = MDS extent;
  # bp_start,bp_end = LD-refined extent). All fragments are drawn (blue = MDS,
  # red = LD refined); the green bar is the merged union when n_merged > 1.
  ax_vals <- mds_dt[[axis]]
  y_range <- range(ax_vals, na.rm = TRUE); y_span <- diff(y_range)
  bar_y   <- y_range[2] + y_span * 0.07
  bar_gap <- y_span * 0.045

  was_merged <- isTRUE(n_merged > 1L) && !is.na(merged_start) && !is.na(merged_end)

  title_str <- sprintf("%s:%d-%d [%s %s]  (%d fragment%s)",
                       chr, as.integer(merged_start), as.integer(merged_end),
                       axis, dir, nrow(frags), if (nrow(frags) > 1L) "s" else "")
  sub_str <- if (was_merged)
    "blue = MDS candidates  |  red = LD refined  |  green = merged union"
  else
    "blue = MDS candidate  |  red = LD refined"

  p <- ggplot(mds_dt, aes(x = window_mid, y = .data[[axis]])) +
    geom_point(size = 0.3, colour = "grey40", alpha = 0.6)

  # Every fragment's MDS (blue) and LD-refined (red) bars
  for (k in seq_len(nrow(frags))) {
    p <- p + annotate("segment", x = frags$start[k], xend = frags$end[k],
                      y = bar_y, yend = bar_y, colour = "steelblue", linewidth = 2.5)
    if (!is.na(frags$bp_start[k]) && !is.na(frags$bp_end[k]))
      p <- p + annotate("segment", x = frags$bp_start[k], xend = frags$bp_end[k],
                        y = bar_y - bar_gap, yend = bar_y - bar_gap,
                        colour = "firebrick", linewidth = 2.5)
  }
  if (was_merged)
    p <- p + annotate("segment", x = merged_start, xend = merged_end,
                      y = bar_y - 2 * bar_gap, yend = bar_y - 2 * bar_gap,
                      colour = "forestgreen", linewidth = 2.5)

  p <- p +
    coord_cartesian(ylim = c(y_range[1] - y_span * 0.05, bar_y + y_span * 0.05)) +
    theme_classic(base_size = 10) +
    labs(title = title_str, subtitle = sub_str, x = "Genomic position", y = axis)
  ggsave(out_file, p, width = 10, height = 2.8, dpi = 150)
}

# Best MDS axis to display a candidate on: among the candidate's supporting axes,
# the one whose windows inside the region are most outlying (highest mean |Z|).
best_axis_for_region <- function(mds_dt, r_start, r_end, cand_axes) {
  cand_axes <- unique(cand_axes[cand_axes %in% names(mds_dt)])
  if (!length(cand_axes)) return(NA_character_)
  inwin <- mds_dt$window_mid >= r_start & mds_dt$window_mid <= r_end
  if (!any(inwin)) return(cand_axes[1])
  scores <- vapply(cand_axes, function(ax) {
    v   <- mds_dt[[ax]]
    med <- median(v, na.rm = TRUE); mad <- median(abs(v - med), na.rm = TRUE)
    if (!is.finite(mad) || mad == 0) return(NA_real_)
    z <- (v - med) / (mad * 1.4826)
    mean(abs(z[inwin]), na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(scores))) return(cand_axes[1])
  cand_axes[which.max(scores)]
}

# Re-genotype a candidate on its full collapsed region [m_start, m_end] using a
# fresh local PCA of the whole region. Returns NULL if the 3-cluster solution
# does not pass (caller then keeps the fragment-based genotype).
regenotype_region <- function(chr, m_start, m_end, sample_ids) {
  if (is.na(m_start) || is.na(m_end)) return(NULL)
  rid        <- paste(chr, as.integer(m_start), as.integer(m_end), sep = "_")
  raw_prefix <- file.path("breakpoints/regeno", rid)
  cmd <- sprintf(
    "plink2 --bfile vcf/%s --chr %s --from-bp %d --to-bp %d --allow-extra-chr --export A --threads 4 --out %s 2>/dev/null",
    chr, chr, as.integer(m_start), as.integer(m_end), raw_prefix)
  if (system(cmd) != 0) return(NULL)
  raw_file <- paste0(raw_prefix, ".raw")
  if (!file.exists(raw_file) || file.info(raw_file)$size == 0) return(NULL)

  tryCatch({
    lpca <- make_lpca_table(raw_file)
    r <- cluster_local_pca_gaps(lpca,
      min_cluster_n = min_cluster_n, min_search_n = min_search_n,
      min_gap_scaled_cutoff = min_gap_scaled_cutoff,
      gap_product_scaled_cutoff = gap_prod_scaled_cutoff,
      dip_p_thresh = dip_p_thresh, r2_k3_min = r2_k3_min,
      cor_h_axis_max = cor_h_axis_max, min_inv_freq = min_inv_freq)
    cc <- classify_local_pca_cluster(r)
    hc <- test_heterozygosity(r$data, p_cutoff = het_p_cutoff, min_group_n = het_min_group_n)
    if (!isTRUE(cc$cluster_pass) || !isTRUE(hc$het_pass)) return(NULL)  # failed -> fallback
    sc <- r$selected_score; bd <- r$bounds; ax <- r$selected_axis
    scores <- lpca[[ax]]
    geno   <- ifelse(scores < bd$L_end, 0L, ifelse(scores < bd$MID_end, 1L, 2L))
    list(
      geno      = setNames(as.integer(geno), lpca$IID),
      lpca_axis = ax, L_end = bd$L_end, MID_end = bd$MID_end,
      stats = data.table(
        cluster_pass = cc$cluster_pass, het_pass = hc$het_pass,
        min_gap_scaled = sc$min_gap_scaled, gap_product_scaled = sc$gap_product_scaled,
        n_L = sc$n_L, n_MID = sc$n_MID, n_R = sc$n_R,
        mean_H_L = hc$mean_H_L, mean_H_MID = hc$mean_H_MID, mean_H_R = hc$mean_H_R,
        p_mid_gt_L = hc$p_mid_gt_L, p_mid_gt_R = hc$p_mid_gt_R, p_mid_gt_homs = hc$p_mid_gt_homs,
        R2_k3 = sc$R2_k3, dip_p = sc$dip_p,
        cor_H_axis_L = sc$cor_H_axis_L, cor_H_axis_R = sc$cor_H_axis_R,
        maf = sc$maf, F_IS = sc$F_IS))
  }, error = function(e) { message("    re-genotype failed: ", e$message); NULL })
}

# Reactive end-trim for a failed merge. When a multi-fragment union does NOT
# re-cluster, iteratively drop the WORST END fragment and re-genotype the shrunken
# span. Only END members matter: regenotype_region() re-clusters the contiguous span
# [min bp_start, max bp_end] over ALL SNPs in it, so dropping an interior member
# changes neither the span nor the SNP set (a no-op) — only the leftmost/rightmost
# fragment moves a boundary. "Worst" = the end fragment with the lowest mean |r| to
# the rest of the current core (the foreign-locus signature). Returns status "kept"
# with the surviving core + the dropped ends, or "reverted" if no trimmed core
# re-clusters. geno_tab = Section-C genotypes (the vectors that formed the group).
trim_regenotype <- function(chr, member_ids, bp_tab, geno_tab, samp_ids) {
  meta <- rbindlist(lapply(member_ids, function(rid) {
    r <- bp_tab[region_id == rid]
    if (!nrow(r)) return(NULL)
    data.table(region_id = rid, bp_start = r$bp_start[1], bp_end = r$bp_end[1])
  }))
  if (!nrow(meta)) return(list(status = "reverted", dropped = character(0)))
  setorder(meta, bp_start)
  geno_of <- function(rid) as.numeric(geno_tab[[rid]])

  dropped <- character(0)
  repeat {
    m_start <- min(meta$bp_start, na.rm = TRUE)
    m_end   <- max(meta$bp_end,   na.rm = TRUE)
    rg      <- regenotype_region(chr, m_start, m_end, samp_ids)
    if (!is.null(rg))
      return(list(status = "kept", members = meta$region_id,
                  m_start = m_start, m_end = m_end, rg = rg, dropped = dropped))
    if (nrow(meta) <= 2L) break          # a 2-member core already failed — nothing to salvage

    # mean |r| of each END fragment to the rest of the current core; drop the lower
    ends <- c(meta$region_id[1L], meta$region_id[nrow(meta)])
    core_cor <- vapply(ends, function(e) {
      others <- setdiff(meta$region_id, e)
      ge <- geno_of(e)
      v  <- vapply(others, function(o)
        abs(suppressWarnings(cor(ge, geno_of(o), use = "pairwise.complete.obs"))),
        numeric(1))
      mean(v, na.rm = TRUE)
    }, numeric(1))
    core_cor[!is.finite(core_cor)] <- -Inf       # a non-correlating end is the most foreign
    worst   <- ends[which.min(core_cor)]
    dropped <- c(dropped, worst)
    meta    <- meta[region_id != worst]
  }
  list(status = "reverted", dropped = dropped)
}

# =============================================================
# SECTION A: MDS SCAN
# =============================================================
cat(sprintf("[%s] === Section A: MDS scan ===\n", ts()))

# Chromosome length
contigs <- fread("ref/contigs.tsv", col.names = c("chr_name", "length"))
chr_len <- contigs[chr_name == chr, length]
if (!length(chr_len) || is.na(chr_len)) stop(sprintf("chr %s not in ref/contigs.tsv", chr))

# Load BED
bed_prefix <- file.path("vcf", chr)
bg  <- BEDMatrix(bed_prefix, simple_names = TRUE)
bim <- fread(paste0(bed_prefix, ".bim"), col.names = c("chr_col","snp","cm","pos","a1","a2"))
N_IND  <- nrow(bg); N_SNPS <- nrow(bim)
cat(sprintf("  %d individuals, %d SNPs\n", N_IND, N_SNPS))

cat(sprintf("  Loading genotype matrix...\n"))
gt_full <- as.matrix(bg)

# Windows: fixed physical size (window_size bp) along the chromosome. SNP count
# per window varies with local density; windows with fewer than min_snps SNPs are
# excluded (sparse windows give unstable local PCA). bp windows give a location-
# stable size floor (min detectable SV ~= MIN_WIN x window_size).
win_start  <- seq(1L, chr_len, by = window_size)
win_end    <- pmin(win_start + window_size - 1L, chr_len)
N_WIN      <- length(win_start)
snp_counts <- tabulate(findInterval(bim$pos, win_start), nbins = N_WIN)

# Eigenstuff: per-window covariance PCA
ES_LEN  <- 1L + 2L + 2L * N_IND
es_mat  <- matrix(NA_real_, nrow = N_WIN, ncol = ES_LEN)

for (i in seq_len(N_WIN)) {
  if (snp_counts[i] < min_snps) next
  snp_idx <- which(bim$pos >= win_start[i] & bim$pos <= win_end[i])
  if (length(snp_idx) < min_snps) next
  gt <- gt_full[, snp_idx, drop = FALSE]
  col_means   <- colMeans(gt, na.rm = TRUE)
  missing_idx <- which(!is.finite(gt), arr.ind = TRUE)
  if (nrow(missing_idx)) gt[missing_idx] <- col_means[missing_idx[, 2]]
  col_vars <- apply(gt, 2, var)
  if (any(is.na(col_vars)) || all(col_vars < 1e-10)) next
  covmat <- cov(t(gt))
  pca    <- eigen(covmat, symmetric = TRUE)
  es_mat[i, ] <- c(sum(covmat^2), pca$values[1], pca$values[2],
                   pca$vectors[, 1], pca$vectors[, 2])
}

valid_idx <- which(!is.na(es_mat[, 1])); N_VALID <- length(valid_idx)
cat(sprintf("  Valid windows: %d / %d\n", N_VALID, N_WIN))

if (N_VALID < 3) {
  cat(sprintf("  WARNING: fewer than 3 valid windows on %s — stopping\n", chr))
  quit(status = 0)
}

es_valid <- es_mat[valid_idx, ]
rm(gt_full, es_mat); gc()

# MDS
pcdist  <- pc_dist(es_valid, npc = npc_dist)
k       <- min(k_mds, N_VALID - 1L)
mds_fit <- cmdscale(pcdist, k = k)

# MDS output table
mds_out <- data.table(
  chr          = chr,
  window_start = win_start[valid_idx],
  window_mid   = (win_start[valid_idx] + win_end[valid_idx]) %/% 2L,
  window_end   = win_end[valid_idx],
  n_snps       = snp_counts[valid_idx]
)
for (ax in seq_len(k))           mds_out[[paste0("MDS", ax)]] <- mds_fit[, ax]
for (ax in seq(k + 1L, k_mds))  mds_out[[paste0("MDS", ax)]] <- NA_real_

fwrite(mds_out, sprintf("mds/%s.mds", chr), sep = "\t", quote = FALSE)
cat(sprintf("  Written mds/%s.mds\n", chr))
rm(es_valid, pcdist, mds_fit); gc()

# Candidate detection
axes  <- paste0("MDS", seq_len(k))
cands <- detect_candidates(mds_out, axes, z_thresh, window_gap, min_win, n_perm, perm_p)
if (!is.null(cands)) cands[, candidate_source := "mds"]

# ---- Optional manual candidates -------------------------------------------
# chr,start,end (1-based, no header) injected here, BEFORE the empty-candidate
# early-exit, so they are tested even on chromosomes the MDS scan flags as empty.
# axis="manual"/dir=NA (Section B re-derives its own local-PCA axis); detection
# stats are NA. They still face the Section B cluster/het gates — manual injection
# bypasses DETECTION only, not validation.
manual_path <- Sys.getenv("MANUAL_CANDIDATES")
if (nzchar(manual_path) && file.exists(manual_path) && file.info(manual_path)$size > 0) {
  cur_chr <- chr
  mc <- tryCatch(
    fread(manual_path, header = FALSE,
          colClasses = c("character", "integer", "integer"),
          col.names = c("chr", "start", "end")),
    error = function(e) { message("  WARNING: could not read MANUAL_CANDIDATES (",
                                  e$message, ") — skipping"); NULL })
  if (!is.null(mc)) {
    mc <- mc[chr == cur_chr]
    if (nrow(mc)) {
      mc[, `:=`(axis = "manual", dir = NA_character_,
                n_windows = NA_integer_, K = NA_integer_, run_floor = NA_integer_,
                peak_z = NA_real_, mean_z = NA_real_, sd_z = NA_real_,
                max_gap_windows = NA_integer_, perm_p = NA_real_,
                candidate_source = "manual")]
      cands <- if (is.null(cands)) mc else rbindlist(list(cands, mc), fill = TRUE)
      cat(sprintf("  Injected %d manual candidate(s) on %s\n", nrow(mc), cur_chr))
    }
  }
}

# Write candidates (empty or populated)
empty_cands <- data.table(chr=character(), start=integer(), end=integer(),
                           axis=character(), dir=character(),
                           n_windows=integer(), K=integer(), run_floor=integer(),
                           peak_z=numeric(), mean_z=numeric(), sd_z=numeric(),
                           max_gap_windows=integer(), perm_p=numeric(),
                           candidate_source=character())
if (is.null(cands) || !nrow(cands)) {
  cat(sprintf("  No candidates on %s\n", chr))
  fwrite(empty_cands, sprintf("mds/candidates/%s.candidates.tsv", chr), sep = "\t")
  quit(status = 0)
}

# Stable, UNIQUE region_id per candidate ROW (coordinate-based; suffixed _2,_3,... on
# any coordinate collision — e.g. a manual region coinciding with an MDS one — so
# region_id is never a duplicated key in Sections B-D2 rownames/colnames/lists).
# region_ids[] is row-aligned to cands and drives both the export and per-row loops.
region_ids <- paste(cands$chr, cands$start, cands$end, sep = "_")
if (anyDuplicated(region_ids))
  region_ids <- ave(region_ids, region_ids, FUN = function(x)
    if (length(x) == 1L) x else paste0(x, c("", paste0("_", 2:length(x)))))

fwrite(cands, sprintf("mds/candidates/%s.candidates.tsv", chr), sep = "\t", quote = FALSE)
cat(sprintf("  Candidate rows: %d (mds=%d, manual=%d)\n", nrow(cands),
            sum(cands$candidate_source == "mds"), sum(cands$candidate_source == "manual")))

# MDS plots per candidate
for (i in seq_len(nrow(cands))) {
  cand    <- cands[i]
  ax_col  <- cand$axis
  if (!ax_col %in% names(mds_out)) ax_col <- "MDS1"   # manual candidates: show on MDS1
  out_png <- sprintf("mds/plots/%s_%d_%d_%s_%s.png",
                     chr, cand$start, cand$end, cand$axis, cand$dir)
  tryCatch({
    ax_vals <- mds_out[[ax_col]]
    y_range <- range(ax_vals, na.rm = TRUE); y_span <- diff(y_range)
    bar_y   <- y_range[2] + y_span * 0.07
    p <- ggplot(mds_out, aes(x = window_mid, y = .data[[ax_col]])) +
      geom_point(size = 0.3, colour = "grey40", alpha = 0.6) +
      annotate("segment", x = cand$start, xend = cand$end,
               y = bar_y, yend = bar_y, colour = "steelblue", linewidth = 2.5) +
      coord_cartesian(ylim = c(y_range[1] - y_span * 0.05, bar_y + y_span * 0.05)) +
      theme_classic(base_size = 10) +
      labs(title = sprintf("%s:%d-%d [%s %s]", chr, cand$start, cand$end, cand$axis, cand$dir),
           x = "Genomic position", y = ax_col)
    ggsave(out_png, p, width = 10, height = 2.8, dpi = 150)
  }, error = function(e) message("  Plot failed: ", out_png, " — ", e$message))
}
cat(sprintf("[%s] Section A complete\n", ts()))

# =============================================================
# SECTION B: LOCAL PCA
# =============================================================
cat(sprintf("[%s] === Section B: Local PCA ===\n", ts()))

# PLINK export — one .raw per region_id (region_ids[] is unique and row-aligned to
# cands; duplicate-coordinate candidates export separately to distinct files).
for (j in seq_len(nrow(cands))) {
  ur <- cands[j]; rid <- region_ids[j]
  cmd <- sprintf(
    "plink2 --bfile vcf/%s --chr %s --from-bp %d --to-bp %d --allow-extra-chr --export A --threads 4 --out local_pca/raw/%s 2>/dev/null",
    ur$chr, ur$chr, ur$start, ur$end, rid)
  if (system(cmd) != 0) message("  WARNING: plink2 export failed for ", rid)
}

# Local PCA + clustering for each candidate row
summary_rows <- vector("list", nrow(cands))
for (i in seq_len(nrow(cands))) {
  cand_row  <- cands[i]
  region_id <- region_ids[i]
  raw_file  <- sprintf("local_pca/raw/%s.raw", region_id)
  lpca_file <- sprintf("local_pca/tables/%s_lpca.tsv", region_id)
  ax_file   <- sprintf("local_pca/axis_scores/%s_axis_scores.tsv", region_id)
  png_file  <- sprintf("local_pca/plots/%s_cluster.png", region_id)
  cat(sprintf("  [%d/%d] %s\n", i, nrow(cands), region_id))

  summary_rows[[i]] <- tryCatch({
    if (!file.exists(raw_file) || file.info(raw_file)$size == 0)
      stop("missing or empty .raw file")
    lpca <- make_lpca_table(raw_file)
    fwrite(copy(lpca)[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 8) else v)],
           lpca_file, sep = "\t")
    res <- cluster_local_pca_gaps(lpca,
      min_cluster_n = min_cluster_n, min_search_n = min_search_n,
      min_gap_scaled_cutoff = min_gap_scaled_cutoff,
      gap_product_scaled_cutoff = gap_prod_scaled_cutoff,
      dip_p_thresh = dip_p_thresh, r2_k3_min = r2_k3_min,
      cor_h_axis_max = cor_h_axis_max, min_inv_freq = min_inv_freq)
    cc <- classify_local_pca_cluster(res)
    hc <- test_heterozygosity(res$data, p_cutoff = het_p_cutoff, min_group_n = het_min_group_n)
    fwrite(res$axis_scores, ax_file, sep = "\t")
    tryCatch(plot_cluster(res, cc, hc, region_id, png_file),
             error = function(e) message("  Cluster plot failed: ", e$message))
    sc <- res$selected_score; bd <- res$bounds
    cbind(cand_row, data.table(
      region_id          = region_id,
      lpca_axis          = res$selected_axis,
      L_end              = bd$L_end,
      MID_end            = bd$MID_end,
      cluster_pass       = cc$cluster_pass,
      het_pass           = hc$het_pass,
      min_gap_scaled     = sc$min_gap_scaled,
      gap_product_scaled = sc$gap_product_scaled,
      n_L = sc$n_L, n_MID = sc$n_MID, n_R = sc$n_R,
      mean_H_L = hc$mean_H_L, mean_H_MID = hc$mean_H_MID, mean_H_R = hc$mean_H_R,
      p_mid_gt_L = hc$p_mid_gt_L, p_mid_gt_R = hc$p_mid_gt_R, p_mid_gt_homs = hc$p_mid_gt_homs,
      R2_k3 = sc$R2_k3, dip_p = sc$dip_p,
      cor_H_axis_L = sc$cor_H_axis_L, cor_H_axis_R = sc$cor_H_axis_R,
      maf = sc$maf, F_IS = sc$F_IS
    ))
  }, error = function(e) make_failure_row(cand_row, region_id, e$message))
}

lpca_summary <- rbindlist(summary_rows, fill = TRUE)
fwrite(lpca_summary, sprintf("local_pca/%s.candidates.local_pca.tsv", chr),
       sep = "\t", quote = FALSE)
cat(sprintf("  Written local_pca/%s.candidates.local_pca.tsv\n", chr))
cat(sprintf("[%s] Section B complete\n", ts()))

# =============================================================
# SECTION C: COLLAPSE
# =============================================================
cat(sprintf("[%s] === Section C: Collapse ===\n", ts()))

passing <- lpca_summary[cluster_pass == TRUE & het_pass == TRUE]
cat(sprintf("  Passing: %d / %d\n", nrow(passing), nrow(lpca_summary)))

if (!nrow(passing)) {
  cat(sprintf("  No passing candidates on %s — stopping\n", chr))
  quit(status = 0)
}

if (nrow(passing) == 1) {
  passing[, c("supporting_axes", "n_collapsed") := .(axis, 1L)]
  collapsed <- passing
} else {
  geno_list  <- lapply(seq_len(nrow(passing)), function(i)
    assign_genotypes(passing$region_id[i], passing$lpca_axis[i],
                     passing$L_end[i], passing$MID_end[i]))
  sample_ids <- names(geno_list[[1]])
  geno_mat   <- do.call(rbind, lapply(geno_list, `[`, sample_ids))
  rownames(geno_mat) <- passing$region_id
  collapsed  <- collapse_candidates(passing, geno_mat, r_thresh = r_thresh)
}

cat(sprintf("  Collapsed: %d candidates\n", nrow(collapsed)))
fwrite(collapsed, sprintf("local_pca/collapsed/%s.collapsed.tsv", chr),
       sep = "\t", quote = FALSE)

# Genotype matrix
sample_ids <- names(assign_genotypes(collapsed$region_id[1], collapsed$lpca_axis[1],
                                     collapsed$L_end[1], collapsed$MID_end[1]))
geno_cols  <- lapply(seq_len(nrow(collapsed)), function(i)
  assign_genotypes(collapsed$region_id[i], collapsed$lpca_axis[i],
                   collapsed$L_end[i], collapsed$MID_end[i])[sample_ids])
geno_wide  <- cbind(data.table(sample = sample_ids), as.data.table(do.call(cbind, geno_cols)))
setnames(geno_wide, c("sample", collapsed$region_id))
fwrite(geno_wide, sprintf("local_pca/collapsed/%s.genotypes.tsv", chr),
       sep = "\t", quote = FALSE)
cat(sprintf("  Genotype matrix: %d samples x %d candidates\n",
            nrow(geno_wide), nrow(collapsed)))
cat(sprintf("[%s] Section C complete\n", ts()))

# =============================================================
# SECTION D: LD BREAKPOINT REFINEMENT (chromosome-wide, no plots yet)
# =============================================================
cat(sprintf("[%s] === Section D: LD breakpoint refinement ===\n", ts()))

geno_wide_full <- fread(sprintf("local_pca/collapsed/%s.genotypes.tsv", chr))

# Load chromosome BED and precompute centred genotype matrix once
# (freed after MDS so reload here; centred form used for all candidates)
ld_bg   <- BEDMatrix(file.path("vcf", chr), simple_names = TRUE)
ld_bim  <- fread(paste0(file.path("vcf", chr), ".bim"),
                 col.names = c("chr_col","snp","cm","pos","a1","a2"))
ld_samp <- rownames(ld_bg)

cat(sprintf("  Loading chromosome genotypes for LD (%d SNPs)...\n", nrow(ld_bim)))
geno_chr <- as.matrix(ld_bg)
col_means_ld <- colMeans(geno_chr, na.rm = TRUE)
na_idx <- which(!is.finite(geno_chr), arr.ind = TRUE)
if (nrow(na_idx) > 0) geno_chr[na_idx] <- col_means_ld[na_idx[, 2]]
geno_c <- sweep(geno_chr, 2, col_means_ld)   # column-centred
sd_g   <- sqrt(colSums(geno_c^2))
rm(geno_chr); gc()

bp_rows <- vector("list", nrow(collapsed))

for (i in seq_len(nrow(collapsed))) {
  cand      <- collapsed[i]
  region_id <- cand$region_id
  cat(sprintf("  [%d/%d] %s\n", i, nrow(collapsed), region_id))

  bp_rows[[i]] <- tryCatch({
    inv_geno <- setNames(as.numeric(geno_wide_full[[region_id]]),
                         geno_wide_full$sample)
    profile  <- compute_ld_chr(geno_c, sd_g, ld_samp, inv_geno, ld_bim, ld_bin_size)

    if (is.null(profile) || !nrow(profile)) {
      message("    No LD profile computed")
      return(cbind(cand, data.table(bp_start = NA_integer_, bp_end = NA_integer_,
                                    bg_r2 = NA_real_, sig_r2 = NA_real_)))
    }

    bps <- find_breakpoints(profile, cand$start, cand$end,
                            ld_min_snr = ld_min_snr, ld_ext_frac = ld_ext_frac)
    fwrite(profile, sprintf("breakpoints/%s_%s_ld_profile.tsv", chr, region_id),
           sep = "\t", quote = FALSE)
    cbind(cand, data.table(bp_start = bps$bp_start, bp_end = bps$bp_end,
                            bg_r2 = bps$bg_r2, sig_r2 = bps$sig_r2))
  }, error = function(e) {
    message("  LD scan failed: ", e$message)
    cbind(cand, data.table(bp_start = NA_integer_, bp_end = NA_integer_,
                            bg_r2 = NA_real_, sig_r2 = NA_real_))
  })
}

bp_table <- rbindlist(bp_rows, fill = TRUE)
cat(sprintf("[%s] Section D complete\n", ts()))

# =============================================================
# SECTION D2: SECOND COLLAPSE + PLOTTING
# =============================================================
cat(sprintf("[%s] === Section D2: Second collapse + plotting ===\n", ts()))

if (nrow(bp_table) > 1) {
  d2      <- collapse_candidates_bp(bp_table, geno_wide_full,
                                    r_thresh = r_thresh_long, max_gap = collapse_buffer)
  final_cands  <- d2$collapsed
  bp_annotated <- merge(bp_table, d2$merge_map, by = "region_id", all.x = TRUE)
} else {
  bp_table[, c("merged_start","merged_end","n_merged","merged_regions") :=
             .(bp_start, bp_end, 1L, region_id)]
  final_cands  <- copy(bp_table)
  bp_annotated <- copy(bp_table)
}

cat(sprintf("  After second collapse: %d candidates (from %d)\n",
            nrow(final_cands), nrow(bp_table)))

# -------------------------------------------------------------
# Re-genotype each collapsed survivor on its full region [merged_start, merged_end]
# and VALIDATE the merge:
#   PASS                 -> keep the merge as one candidate; adopt the whole-region
#                           genotype + stats (genotype_source = "merged").
#   FAIL & n_merged > 1  -> REVERT the long-range merge: emit each constituent fragment
#                           as its own candidate, restoring the fragment's coordinates
#                           and Section-C genotype (genotype_source = "fragment_reverted").
#                           A long-range join survives only if the union re-clusters as
#                           one inversion; otherwise it is undone, COORDINATES INCLUDED.
#   FAIL & n_merged == 1 -> singleton whose LD-refined region did not re-cluster; keep it
#                           with its fragment genotype (genotype_source = "fragment").
# A reverted merge expands one row into several, so final_cands is rebuilt below rather
# than mutated in place.
# -------------------------------------------------------------
dir.create("breakpoints/regeno", showWarnings = FALSE, recursive = TRUE)
samp_ids <- geno_wide_full$sample

out_rows <- list(); out_geno <- list()
n_pass <- 0L; n_revert <- 0L; n_trim <- 0L
for (i in seq_len(nrow(final_cands))) {
  fc <- final_cands[i]
  rg <- regenotype_region(chr, fc$merged_start, fc$merged_end, samp_ids)

  if (!is.null(rg)) {
    # Union re-clusters -> keep the (possibly merged) candidate, adopt whole-region call
    row <- copy(fc)
    row[, genotype_source := "merged"]
    row[, lpca_axis := rg$lpca_axis]
    if ("L_end"   %in% names(row)) row[, L_end   := rg$L_end]
    if ("MID_end" %in% names(row)) row[, MID_end := rg$MID_end]
    for (cn in names(rg$stats))
      if (cn %in% names(row)) row[, (cn) := rg$stats[[cn]]]
    out_rows[[length(out_rows) + 1L]] <- row
    out_geno[[length(out_geno) + 1L]] <- as.integer(rg$geno[samp_ids])
    n_pass <- n_pass + 1L

  } else if (isTRUE(fc$n_merged > 1L)) {
    # Union fails -> REACTIVE END-TRIM before reverting. Peel the worst end fragment
    # (least correlated with the core) and re-genotype the shrunken span; keep the
    # surviving core as a merge and emit the peeled end(s) as singletons. Only if no
    # trimmed core re-clusters do we fall back to the full revert.
    member_ids <- strsplit(as.character(fc$merged_regions), ",")[[1]]
    tr <- trim_regenotype(chr, member_ids, bp_table, geno_wide_full, samp_ids)

    if (isTRUE(tr$status == "kept") && length(tr$members) >= 2L &&
        length(tr$dropped) > 0L) {
      # Trimmed core re-clusters -> keep it (whole-region call on the trimmed span);
      # label it with the best-clustering surviving fragment (a survivor, so its id
      # never collides with a dropped piece's row).
      rg     <- tr$rg
      surv   <- bp_table[region_id %in% tr$members]
      rep_id <- surv[order(-R2_k3)][1]$region_id
      row <- copy(fc)
      row[, `:=`(region_id       = rep_id,
                 merged_start    = tr$m_start,
                 merged_end      = tr$m_end,
                 merged_regions  = paste(tr$members, collapse = ","),
                 n_merged        = length(tr$members),
                 genotype_source = "merged_trimmed",
                 lpca_axis       = rg$lpca_axis)]
      if ("L_end"   %in% names(row)) row[, L_end   := rg$L_end]
      if ("MID_end" %in% names(row)) row[, MID_end := rg$MID_end]
      for (cn in names(rg$stats))
        if (cn %in% names(row)) row[, (cn) := rg$stats[[cn]]]
      out_rows[[length(out_rows) + 1L]] <- row
      out_geno[[length(out_geno) + 1L]] <- as.integer(rg$geno[samp_ids])
      n_pass <- n_pass + 1L

      # peeled end fragment(s) -> singletons (own coords + Section-C genotype)
      for (frid in tr$dropped) {
        frow <- bp_table[region_id == frid]
        if (!nrow(frow)) next
        frow <- copy(frow[1])
        frow[, `:=`(merged_start    = bp_start,
                    merged_end      = bp_end,
                    n_merged        = 1L,
                    merged_regions  = region_id,
                    genotype_source = "fragment_trimmed")]
        out_rows[[length(out_rows) + 1L]] <- frow
        out_geno[[length(out_geno) + 1L]] <- as.integer(geno_wide_full[[frid]])
      }
      n_trim <- n_trim + length(tr$dropped)

    } else {
      # No trimmed core re-clusters -> REVERT: restore each constituent fragment
      for (frid in strsplit(as.character(fc$merged_regions), ",")[[1]]) {
        frow <- bp_table[region_id == frid]
        if (!nrow(frow)) next
        frow <- copy(frow[1])
        frow[, `:=`(merged_start    = bp_start,
                    merged_end      = bp_end,
                    n_merged        = 1L,
                    merged_regions  = region_id,
                    genotype_source = "fragment_reverted")]
        out_rows[[length(out_rows) + 1L]] <- frow
        out_geno[[length(out_geno) + 1L]] <- as.integer(geno_wide_full[[frid]])
      }
      n_revert <- n_revert + 1L
    }

  } else {
    # Singleton whose LD-refined region did not re-cluster -> keep fragment genotype
    row <- copy(fc)
    row[, genotype_source := "fragment"]
    out_rows[[length(out_rows) + 1L]] <- row
    out_geno[[length(out_geno) + 1L]] <- as.integer(geno_wide_full[[fc$region_id]])
  }
}
final_cands     <- rbindlist(out_rows, fill = TRUE)
geno_cols_final <- out_geno
cat(sprintf("  Region re-genotyping: %d merge(s)/region(s) kept, %d end-fragment(s) trimmed, %d merge(s) reverted -> %d final candidates\n",
            n_pass, n_trim, n_revert, nrow(final_cands)))

# Write final breakpoints (now includes genotype_source + whole-region stats)
fwrite(final_cands, sprintf("breakpoints/%s.breakpoints.tsv", chr),
       sep = "\t", quote = FALSE)

# Write final genotype matrix (whole-region genotype where re-genotyping passed)
geno_final <- cbind(data.table(sample = samp_ids),
                    as.data.table(do.call(cbind, geno_cols_final)))
setnames(geno_final, c("sample", final_cands$region_id))
fwrite(geno_final, sprintf("breakpoints/%s.genotypes.tsv", chr),
       sep = "\t", quote = FALSE)

# -------------------------------------------------------------
# Plots: one per final candidate, on its best supporting axis (highest
# mean |Z| over the region), showing every collapsed fragment's MDS
# (blue) and LD-refined (red) bars plus the merged union (green).
# -------------------------------------------------------------
n_plots <- 0L
for (i in seq_len(nrow(final_cands))) {
  fc       <- final_cands[i]
  frag_ids <- strsplit(as.character(fc$merged_regions), ",")[[1]]
  frags    <- bp_annotated[region_id %in% frag_ids, .(start, end, bp_start, bp_end)]
  if (!nrow(frags))
    frags <- bp_annotated[region_id == fc$region_id, .(start, end, bp_start, bp_end)]
  sup_axes <- strsplit(as.character(fc$supporting_axes), ",")[[1]]
  best_ax  <- best_axis_for_region(mds_out, fc$merged_start, fc$merged_end, sup_axes)
  if (is.na(best_ax) || !best_ax %in% names(mds_out)) best_ax <- "MDS1"   # manual-only: MDS1
  out_png  <- sprintf("breakpoints/plots/%s_%s.png", chr, fc$region_id)
  tryCatch({
    plot_mds_with_breakpoints(
      mds_out, chr, frags, fc$merged_start, fc$merged_end,
      n_merged = if (!is.null(fc$n_merged) && !is.na(fc$n_merged)) fc$n_merged else 1L,
      axis = best_ax, dir = fc$dir, out_file = out_png)
    n_plots <- n_plots + 1L
  }, error = function(e) message("  Plot failed: ", e$message))
}
cat(sprintf("  Plots written: %d\n", n_plots))
cat(sprintf("[%s] Script 2 complete — %s\n", ts(), chr))

RSCRIPT

echo "[$(date)] Script 2 complete — ${chr}"
