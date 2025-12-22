# =========================================================================
# SETUP: Load Packages and Define Paths
# =========================================================================
graphics.off()

# Load necessary libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(vegan)         # For betadisper
library(ade4)          # For Cailliez correction

# Define paths
base_path <- "/Users/mpwilliams/Research/BakrAwa/analyses/Scripts"
qpwave_path <- file.path(base_path, "qpWave/rds/qpWave_regionalpairwise_PROCESSED.rds")
f3_path <- file.path(base_path, "f3stats/outgroupF3/rds/oF3_OldAfrica__qpWaveRegionalInds.rds")
output_dir <- file.path(base_path, "F3_Distribution_Results")
plot_dir <- file.path(output_dir, "Within_Region_Analysis")

# Create output directories for each correction method
plot_dir_cailliez <- file.path(plot_dir, "Cailliez_Correction")
plot_dir_sqrt <- file.path(plot_dir, "Sqrt_Transformation")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(plot_dir_cailliez)) dir.create(plot_dir_cailliez, recursive = TRUE)
if (!dir.exists(plot_dir_sqrt)) dir.create(plot_dir_sqrt, recursive = TRUE)

# =========================================================================
# STEP 1: Data Loading and Preprocessing - USE 1-est DIRECTLY
# =========================================================================

cat("Loading and processing data...\n")

combined_data <- readRDS(file = qpwave_path)
oF3_data <- data.table(readRDS(file = f3_path))

# Clean F3 Data: Remove NANs and self-comparisons
oF3_data <- oF3_data[pop2 != "NAN" & pop3 != "NAN" & pop2 != pop3]
oF3_data = oF3_data[est < 0.3]

# Use 1-est directly (this is the 1 - F3 value, where higher = less affinity)
if ("1-est" %in% names(oF3_data)) {
  setnames(oF3_data, "1-est", "one_minus_f3")
  cat("Renamed '1-est' to 'one_minus_f3'\n")
} else if ("est" %in% names(oF3_data)) {
  oF3_data[, one_minus_f3 := 1 - est]
  cat("Created 'one_minus_f3' from '1 - est'\n")
} else {
  stop("Neither '1-est' nor 'est' found in data!")
}

cat("1-F3 column summary:\n")
print(summary(oF3_data$one_minus_f3))

# Clean Combined Data
combined_data <- combined_data[!is.na(Pop1_sample_ID)]
setDT(combined_data)

# -------------------------------------------------------------------------
# STEP 2: Robust ID Extraction (For BOTH Individuals)
# -------------------------------------------------------------------------

oF3_data[, id_1 := sub(".*_&_", "", pop2)]
oF3_data[, id_2 := sub(".*_&_", "", pop3)]

# Prepare Metadata Reference
meta_ref <- unique(combined_data[, .(sample_ID, Date_mean_BP, PCA_RegionalGroup)])
meta_ref <- meta_ref[!is.na(sample_ID) & sample_ID != ""]

# Define Meta-Region Map
meta_map <- list(
  C_N_NW_SW_Anatolia = c("cAnatolia_ChL", "swAnatolia_BA", "cAnatolia_IA", "nwAnatolia_ChL", 
                         "swAnatolia_IA", "cAnatolia_Hist", "nAnatolia_Hist", "nAnatolia_BA", 
                         "nAnatolia_ChL", "cAnatolia_BA"),
  S_Levant = c("sLevant_BA", "sLevant_IA"),
  N_Levant_S_E_Anatolia = c("nLevant_ChL", "nLevant_BA", "eAnatolia_BA", "eAnatolia_ChL",
                            "sAnatolia_BA", "sAnatolia_IA"),
  SE_Anatolia_N_Mesopotamia = c("seAnatolia_BA", "nMesopotamia_BA", "seAnatolia_Hist", "seAnatolia_ChL"),
  Caucasus = c("sCaucasus_IA", "sCaucasus_ChL", "sCaucasus_BA", "sCaucasus_Hist"),
  Zagros = c("nwZagros_IA", "cZagros_ChL", "nwZagros_BA")
)

rev_map <- data.table(
  PCA_RegionalGroup = unlist(meta_map),
  meta_region = rep(names(meta_map), lengths(meta_map))
)

meta_ref <- merge(meta_ref, rev_map, by = "PCA_RegionalGroup", all.x = FALSE)

meta_ref[, Time_Period := fcase(
  Date_mean_BP > 4200 & Date_mean_BP < 7700, "TP1_4200_7700BP",
  Date_mean_BP >= 3400 & Date_mean_BP <= 4200, "TP2_3400_4200BP",
  Date_mean_BP < 3400, "TP3_LT_3400BP",
  default = NA
)]
meta_ref <- meta_ref[!is.na(Time_Period)]

# -------------------------------------------------------------------------
# STEP 3: Merge Metadata to F3 Data (Double Merge)
# -------------------------------------------------------------------------

plot_data <- merge(oF3_data, meta_ref, by.x = "id_1", by.y = "sample_ID", all.x = FALSE)
setnames(plot_data, c("meta_region", "Time_Period", "PCA_RegionalGroup"), 
         c("Region_1", "TP_1", "Group_1"))

plot_data <- merge(plot_data, meta_ref, by.x = "id_2", by.y = "sample_ID", all.x = FALSE)
setnames(plot_data, c("meta_region", "Time_Period", "PCA_RegionalGroup"), 
         c("Region_2", "TP_2", "Group_2"))

plot_data <- plot_data[TP_1 == TP_2]

cat("Plot data dimensions after merging:", nrow(plot_data), "rows\n")

# -------------------------------------------------------------------------
# STEP 4: Define Consistent Colors for Regions
# -------------------------------------------------------------------------
unique_regions_all <- unique(meta_ref$meta_region)
region_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(unique_regions_all)), unique_regions_all)

# =========================================================================
# FUNCTION: Create Separate Region Legend
# =========================================================================

create_region_legend <- function(region_colors, output_file, title_text = "Region Legend") {
  
  unique_regs <- names(region_colors)
  n_regions <- length(unique_regs)
  
  # PDF version
  pdf(output_file, width = 4, height = max(3, n_regions * 0.5 + 1))
  par(mar = c(1, 1, 2, 1))
  plot(NA, xlim = c(0, 10), ylim = c(0, n_regions + 1), 
       xlab = "", ylab = "", axes = FALSE, main = title_text)
  
  for (k in seq_along(unique_regs)) {
    y_pos <- n_regions - k + 1
    rect(0.5, y_pos - 0.4, 1.5, y_pos + 0.4, 
         col = region_colors[unique_regs[k]], border = "black", lwd = 0.5)
    text(2, y_pos, unique_regs[k], adj = 0, cex = 1)
  }
  dev.off()
  
  # PNG version
  png_file <- sub("\\.pdf$", ".png", output_file)
  png(png_file, width = 400, height = max(300, n_regions * 50 + 100), res = 100)
  par(mar = c(1, 1, 2, 1))
  plot(NA, xlim = c(0, 10), ylim = c(0, n_regions + 1), 
       xlab = "", ylab = "", axes = FALSE, main = title_text)
  
  for (k in seq_along(unique_regs)) {
    y_pos <- n_regions - k + 1
    rect(0.5, y_pos - 0.4, 1.5, y_pos + 0.4, 
         col = region_colors[unique_regs[k]], border = "black", lwd = 0.5)
    text(2, y_pos, unique_regs[k], adj = 0, cex = 1)
  }
  dev.off()
}

# =========================================================================
# FUNCTION: Create Figure 2-style Heatmap with 45° Rotated Dendrogram
# =========================================================================

create_fig2_style_heatmap <- function(mat, row_regions, region_colors, title_text, output_file) {
  
  # 1. Hierarchical clustering with Ward's method
  dist_mat <- as.dist(mat)
  hc <- hclust(dist_mat, method = "ward.D2")
  
  # Reorder matrix by dendrogram
  ord <- hc$order
  mat_ordered <- mat[ord, ord]
  row_regions_ordered <- row_regions[ord]
  sample_names <- rownames(mat_ordered)
  n <- nrow(mat_ordered)
  
  # Get value range for color scale
  upper_tri_vals <- mat_ordered[upper.tri(mat_ordered, diag = FALSE)]
  val_range <- range(upper_tri_vals, na.rm = TRUE)
  
  # Color function
  n_colors <- 100
  col_func <- colorRampPalette(c("darkred", "red", "orange", "yellow", "lightyellow"))(n_colors)
  
  val_to_col <- function(v) {
    if (is.na(v)) return("gray50")
    idx <- round((v - val_range[1]) / diff(val_range) * (n_colors - 1)) + 1
    idx <- max(1, min(n_colors, idx))
    return(col_func[idx])
  }
  
  # 2. Extract dendrogram segments and apply 45° rotation
  get_dendrogram_segments <- function(hc) {
    n_leaves <- length(hc$order)
    merge <- hc$merge
    height <- hc$height
    order <- hc$order
    
    leaf_pos <- numeric(n_leaves)
    leaf_pos[order] <- 1:n_leaves
    
    segments_list <- list()
    seg_idx <- 1
    node_x <- numeric(n_leaves - 1)
    node_y <- numeric(n_leaves - 1)
    
    get_x <- function(node) { if (node < 0) return(leaf_pos[-node]) else return(node_x[node]) }
    get_y <- function(node) { if (node < 0) return(0) else return(node_y[node]) }
    
    for (i in 1:(n_leaves - 1)) {
      child1 <- merge[i, 1]; child2 <- merge[i, 2]
      x1 <- get_x(child1); x2 <- get_x(child2)
      y1 <- get_y(child1); y2 <- get_y(child2)
      
      node_x[i] <- (x1 + x2) / 2
      node_y[i] <- height[i]
      
      segments_list[[seg_idx]] <- c(x1, y1, x1, height[i]); seg_idx <- seg_idx + 1
      segments_list[[seg_idx]] <- c(x2, y2, x2, height[i]); seg_idx <- seg_idx + 1
      segments_list[[seg_idx]] <- c(x1, height[i], x2, height[i]); seg_idx <- seg_idx + 1
    }
    do.call(rbind, segments_list)
  }
  
  dend_segs <- get_dendrogram_segments(hc)
  max_height <- max(hc$height)
  scale_factor <- 0.7
  height_scale <- (n * scale_factor) / max_height
  
  transform_to_rotated <- function(x, y) {
    y_scaled <- y * height_scale
    x_new <- x - y_scaled
    y_new <- x + y_scaled
    return(c(x_new, y_new))
  }
  
  dend_segs_transformed <- matrix(0, nrow = nrow(dend_segs), ncol = 4)
  for (i in 1:nrow(dend_segs)) {
    start <- transform_to_rotated(dend_segs[i, 1], dend_segs[i, 2])
    end <- transform_to_rotated(dend_segs[i, 3], dend_segs[i, 4])
    dend_segs_transformed[i, ] <- c(start[1], start[2], end[1], end[2])
  }
  
  tip_positions <- matrix(0, nrow = n, ncol = 2)
  for (i in 1:n) {
    transformed <- transform_to_rotated(i, 0)
    tip_positions[i, ] <- transformed
  }
  tip_colors <- region_colors[row_regions_ordered]
  
  # 3. Create the plot
  dend_extent <- n * scale_factor
  x_min <- 1 - dend_extent - 0.5; x_max <- n + 0.5
  y_min <- 0.5; y_max <- n + dend_extent + 0.5
  
  # PDF Output
  pdf(output_file, width = 20, height = 40)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 10))
  
  # Legend Panel
  par(mar = c(8, 2, 8, 1))
  legend_vals <- seq(val_range[1], val_range[2], length.out = n_colors)
  image(1, legend_vals, t(as.matrix(legend_vals)), col = col_func, axes = FALSE, xlab = "", ylab = "")
  axis(2, at = pretty(legend_vals, n = 5), las = 1, cex.axis = 0.9)
  mtext(expression(1 - F[3]), side = 2, line = 3, cex = 1)
  mtext("Low affinity", side = 3, line = 0.5, cex = 0.8, font = 2)
  mtext("High affinity", side = 1, line = 0.5, cex = 0.8, font = 2)
  box()
  
  # Main Plot Panel
  par(mar = c(10, 2, 4, 6))
  plot(NA, xlim = c(x_min, x_max), ylim = c(y_min, y_max), xlab = "", ylab = "", axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
  
  # Heatmap
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = "gray50", border = NA)
      } else if (j > i) {
        col <- val_to_col(mat_ordered[i, j])
        rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = col, border = NA)
      }
    }
  }
  
  # Dendrogram
  for (i in 1:nrow(dend_segs_transformed)) {
    segments(dend_segs_transformed[i, 1], dend_segs_transformed[i, 2],
             dend_segs_transformed[i, 3], dend_segs_transformed[i, 4], lwd = 0.6, col = "black")
  }
  points(tip_positions[, 1], tip_positions[, 2], pch = 21, bg = tip_colors, col = "black", cex = 1.2, lwd = 0.5)
  
  # Labels
  text(x = 1:n, y = 0.2, labels = sample_names, srt = 90, adj = 1, xpd = TRUE, cex = 0.4)
  text(x = n + 0.7, y = 1:n, labels = sample_names, adj = 0, xpd = TRUE, cex = 0.4)
  
  # Region Bar
  for (i in 1:n) {
    rect(n + 4, i - 0.45, n + 4.6, i + 0.45, col = region_colors[row_regions_ordered[i]], border = NA, xpd = TRUE)
  }
  title(main = title_text, line = 2, cex.main = 1.3, font.main = 2)
  dev.off()
  
  # PNG Output (Simplified for brevity, same logic as PDF)
  png_file <- sub("\\.pdf$", ".png", output_file)
  png(png_file, width = 4000, height = 4500, res = 150)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 10))
  par(mar = c(8, 2, 8, 1))
  image(1, legend_vals, t(as.matrix(legend_vals)), col = col_func, axes = FALSE, xlab = "", ylab = "")
  axis(2, at = pretty(legend_vals, n = 5), las = 1, cex.axis = 0.9)
  box()
  par(mar = c(10, 2, 4, 6))
  plot(NA, xlim = c(x_min, x_max), ylim = c(y_min, y_max), xlab = "", ylab = "", axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = "gray50", border = NA)
      else if (j > i) rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = val_to_col(mat_ordered[i, j]), border = NA)
    }
  }
  for (i in 1:nrow(dend_segs_transformed)) segments(dend_segs_transformed[i, 1], dend_segs_transformed[i, 2], dend_segs_transformed[i, 3], dend_segs_transformed[i, 4], lwd = 0.6, col = "black")
  points(tip_positions[, 1], tip_positions[, 2], pch = 21, bg = tip_colors, col = "black", cex = 1.2, lwd = 0.5)
  text(x = 1:n, y = 0.2, labels = sample_names, srt = 90, adj = 1, xpd = TRUE, cex = 0.4)
  text(x = n + 0.7, y = 1:n, labels = sample_names, adj = 0, xpd = TRUE, cex = 0.4)
  for (i in 1:n) rect(n + 4, i - 0.45, n + 4.6, i + 0.45, col = region_colors[row_regions_ordered[i]], border = NA, xpd = TRUE)
  title(main = title_text, line = 2, cex.main = 1.3, font.main = 2)
  dev.off()
  
  # 4. Create SEPARATE legend file
  legend_file <- sub("\\.pdf$", "_Legend.pdf", output_file)
  unique_regs <- unique(row_regions_ordered)
  subset_colors <- region_colors[unique_regs]
  create_region_legend(subset_colors, legend_file, title_text = paste0("Region Legend\n", basename(output_file)))
}

# =========================================================================
# FUNCTION: Run betadisper analysis with specified correction method
# =========================================================================

run_betadisper_analysis <- function(mat, groups, tp, correction_method, output_dir, region_colors) {
  
  cat(sprintf("\n--- Running betadisper analysis for %s (%s correction) ---\n", tp, correction_method))
  
  # Create distance matrix based on correction method
  if (correction_method == "cailliez") {
    dist_mat_raw <- as.dist(mat)
    dist_mat <- cailliez(dist_mat_raw)
    cat(sprintf("  Applied Cailliez correction (constant c2 = %.6f)\n", attr(dist_mat, "c2")))
  } else if (correction_method == "sqrt") {
    mat_sqrt <- sqrt(mat)
    dist_mat <- as.dist(mat_sqrt)
    cat("  Applied square root transformation to distances\n")
  } else {
    stop("Unknown correction method. Use 'cailliez' or 'sqrt'.")
  }
  
  # Run betadisper
  bd <- betadisper(dist_mat, groups, type = "median")
  
  # ANOVA test
  anova_result <- anova(bd)
  anova_p <- anova_result$`Pr(>F)`[1]
  anova_f <- anova_result$`F value`[1]
  
  # Permutation test
  perm_result <- permutest(bd, permutations = 999)
  perm_p <- perm_result$tab$`Pr(>F)`[1]
  perm_f <- perm_result$tab$F[1]
  
  # Tukey HSD pairwise comparisons
  tukey_result <- TukeyHSD(bd)
  tukey_df <- as.data.frame(tukey_result$group)
  tukey_df$comparison <- rownames(tukey_df)
  tukey_df$p.signif <- ifelse(tukey_df$`p adj` < 0.001, "***",
                              ifelse(tukey_df$`p adj` < 0.01, "**",
                                     ifelse(tukey_df$`p adj` < 0.05, "*", "ns")))
  tukey_df$Time_Period <- tp
  tukey_df$Correction_Method <- correction_method
  
  # Extract group dispersions
  dispersion_summary <- data.frame(
    Region = names(bd$group.distances),
    Mean_Distance_to_Centroid = tapply(bd$distances, bd$group, mean),
    Median_Distance_to_Centroid = tapply(bd$distances, bd$group, median),
    SD_Distance = tapply(bd$distances, bd$group, sd),
    N = as.numeric(table(bd$group)),
    Time_Period = tp,
    Correction_Method = correction_method
  )
  dispersion_summary <- dispersion_summary[order(dispersion_summary$Mean_Distance_to_Centroid), ]
  dispersion_summary$Rank_Homogeneity <- 1:nrow(dispersion_summary)
  
  # Create plots
  # 1. Base R boxplot
  pdf(file.path(output_dir, paste0("Betadisper_Boxplot_", correction_method, "_", tp, ".pdf")), width = 12, height = 9)
  par(mar = c(12, 5, 4, 2))
  boxplot(bd, main = paste0("Distance to Centroid by Region\n", tp, " (", correction_method, " correction)"),
          xlab = "", ylab = "Distance to centroid", col = region_colors[levels(bd$group)],
          las = 2, cex.axis = 0.8, cex.lab = 1.1, cex.main = 1.2)
  mtext("Region", side = 1, line = 10, cex = 1.1)
  dev.off()
  
  # 2. PCoA plot
  pdf(file.path(output_dir, paste0("Betadisper_PCoA_", correction_method, "_", tp, ".pdf")), width = 10, height = 10)
  plot(bd, main = paste0("PCoA of Group Dispersions\n", tp, " (", correction_method, " correction)"),
       col = region_colors[levels(bd$group)], hull = FALSE, ellipse = TRUE, conf = 0.68)
  dev.off()
  
  # 3. ggplot boxplot
  dist_df <- data.frame(Distance = bd$distances, Region = bd$group)
  p_box <- ggplot(dist_df, aes(x = reorder(Region, Distance, FUN = median), y = Distance, fill = Region)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    scale_fill_manual(values = region_colors) +
    labs(title = paste0("Within-Region Genetic Dispersion: ", tp),
         subtitle = paste0(correction_method, " correction | Permutation p = ", signif(perm_p, 3)),
         x = "Region (ordered by median distance)", y = "Distance to Regional Centroid") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(output_dir, paste0("Betadisper_ggplot_", correction_method, "_", tp, ".pdf")), plot = p_box, width = 10, height = 8)
  
  return(list(betadisper = bd, anova_p = anova_p, anova_f = anova_f, perm_p = perm_p, perm_f = perm_f,
              tukey = tukey_df, dispersion = dispersion_summary, permutation = perm_result))
}

# =========================================================================
# STEP 5: Analysis Loop - ALL TESTS USE 1-F3 (one_minus_f3)
# =========================================================================

within_region_dt <- plot_data[Region_1 == Region_2]
time_periods <- unique(within_region_dt$TP_1)

# Storage for results
all_stats_results <- list()
all_betadisper_results_cailliez <- list(); all_tukey_results_cailliez <- list(); all_dispersion_summaries_cailliez <- list()
all_betadisper_results_sqrt <- list(); all_tukey_results_sqrt <- list(); all_dispersion_summaries_sqrt <- list()

cat("\n=== Starting Analysis Loop ===\n")

# Create master region legend
master_legend_file <- file.path(plot_dir, "Region_Legend_Master.pdf")
create_region_legend(region_colors, master_legend_file, title_text = "Region Legend")

for (tp in time_periods) {
  cat(sprintf("\n========== Processing %s ==========\n", tp))
  
  subset_dt <- within_region_dt[TP_1 == tp]
  unique_regs_in_tp <- unique(subset_dt$Region_1)
  
  if (nrow(subset_dt) < 10 || length(unique_regs_in_tp) < 2) {
    cat("Skipping: Insufficient data.\n"); next
  }
  
  # 5.1 Statistical Testing (Kruskal-Wallis & Wilcoxon)
  kw_test <- kruskal.test(one_minus_f3 ~ Region_1, data = subset_dt)
  kw_p <- signif(kw_test$p.value, 3)
  
  pairs <- combn(unique_regs_in_tp, 2, simplify = FALSE)
  
  stats_list <- lapply(pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    vals1 <- subset_dt[Region_1 == g1]$one_minus_f3
    vals2 <- subset_dt[Region_1 == g2]$one_minus_f3
    w_test <- wilcox.test(vals1, vals2)
    med1 <- median(vals1, na.rm = TRUE)
    med2 <- median(vals2, na.rm = TRUE)
    winner <- ifelse(med1 < med2, g1, g2)
    if (w_test$p.value > 0.05) winner <- "No Significant Difference"
    
    data.frame(
      Time_Period = tp, Group1 = g1, Group2 = g2,
      p = w_test$p.value,
      p.signif = ifelse(w_test$p.value < 0.001, "***", 
                        ifelse(w_test$p.value < 0.01, "**", 
                               ifelse(w_test$p.value < 0.05, "*", "ns"))),
      Median_Group1_1minusF3 = med1, 
      Median_Group2_1minusF3 = med2,
      More_Homogenous_Region = winner,
      Kruskal_Wallis_p = kw_p
    )
  })
  
  period_stats <- do.call(rbind, stats_list)
  all_stats_results[[tp]] <- period_stats
  
  # 5.2 Violin/Butterfly Plots
  p <- ggplot(subset_dt, aes(x = Region_1, y = one_minus_f3, fill = Region_1)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.5) +
    stat_compare_means(comparisons = pairs, method = "wilcox.test", 
                       label = "p.signif", hide.ns = TRUE, vjust = 0.5) +
    labs(title = paste0("Within-Region Genetic Affinity: ", tp),
         subtitle = paste0("Global Kruskal-Wallis p = ", kw_p, "\n(Lower 1-F3 = Higher Genetic Affinity)"),
         x = "Region", y = expression(1 - F[3])) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_fill_manual(values = region_colors)
  
  ggsave(filename = file.path(plot_dir, paste0("Violin_WithinRegion_1minusF3_", tp, ".pdf")), plot = p, width = 12, height = 8)
  
  # 5.3 Heatmap & Betadisper
  heatmap_dt <- plot_data[TP_1 == tp, .(id_1, id_2, one_minus_f3, Region_1, Region_2)]
  d1 <- heatmap_dt[, .(id_A = id_1, id_B = id_2, one_minus_f3)]
  d2 <- heatmap_dt[, .(id_A = id_2, id_B = id_1, one_minus_f3)]
  sym_dt <- unique(rbind(d1, d2))
  sym_dt <- sym_dt[, .(one_minus_f3 = mean(one_minus_f3, na.rm = TRUE)), by = .(id_A, id_B)]
  
  mat <- pivot_wider(sym_dt, names_from = id_B, values_from = one_minus_f3, values_fill = NA) %>%
    column_to_rownames("id_A") %>% as.matrix()
  
  mat <- mat[rowSums(is.na(mat)) < (ncol(mat)/3), colSums(is.na(mat)) < (nrow(mat)/3)]
  common_ids <- intersect(rownames(mat), colnames(mat))
  mat <- mat[common_ids, common_ids]
  diag(mat) <- NA
  
  mat_for_clust <- mat
  mat_for_clust[is.na(mat_for_clust)] <- median(mat_for_clust, na.rm = TRUE)
  
  if (nrow(mat) > 5 && ncol(mat) > 5) {
    row_ids <- rownames(mat)
    row_regions <- meta_ref$meta_region[match(row_ids, meta_ref$sample_ID)]
    
    output_file <- file.path(plot_dir, paste0("Fig2_Heatmap_1minusF3_", tp, ".pdf"))
    create_fig2_style_heatmap(mat_for_clust, row_regions, region_colors, paste0("1-F3(X,Y; Mbuti.DG) - ", tp), output_file)
    
    # Betadisper
    if (!isSymmetric(mat_for_clust, tol = 1e-10)) mat_for_clust <- (mat_for_clust + t(mat_for_clust)) / 2
    groups <- factor(row_regions)
    
    # Cailliez
    tryCatch({
      res_c <- run_betadisper_analysis(mat_for_clust, groups, tp, "cailliez", plot_dir_cailliez, region_colors)
      all_betadisper_results_cailliez[[tp]] <- res_c
      all_tukey_results_cailliez[[tp]] <- res_c$tukey
      all_dispersion_summaries_cailliez[[tp]] <- res_c$dispersion
    }, error = function(e) cat(sprintf("  ERROR in Cailliez analysis for %s: %s\n", tp, e$message)))
    
    # Sqrt
    tryCatch({
      res_s <- run_betadisper_analysis(mat_for_clust, groups, tp, "sqrt", plot_dir_sqrt, region_colors)
      all_betadisper_results_sqrt[[tp]] <- res_s
      all_tukey_results_sqrt[[tp]] <- res_s$tukey
      all_dispersion_summaries_sqrt[[tp]] <- res_s$dispersion
    }, error = function(e) cat(sprintf("  ERROR in Sqrt analysis for %s: %s\n", tp, e$message)))
  }
}

# =========================================================================
# STEP 6: Save Intermediate Results Tables
# =========================================================================

cat("\n=== Saving Intermediate Results ===\n")

# Wilcoxon
final_stats_df <- do.call(rbind, all_stats_results)
final_stats_df$p.adj <- p.adjust(final_stats_df$p, method = "BH")
write.table(final_stats_df, file = file.path(plot_dir, "Supplementary_Data_8.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Cailliez
if (length(all_betadisper_results_cailliez) > 0) {
  global_tests_cailliez <- data.frame(
    Time_Period = names(all_betadisper_results_cailliez),
    Correction_Method = "Cailliez",
    Permutation_F = sapply(all_betadisper_results_cailliez, function(x) x$perm_f),
    Permutation_p = sapply(all_betadisper_results_cailliez, function(x) x$perm_p),
    ANOVA_F = sapply(all_betadisper_results_cailliez, function(x) x$anova_f),
    ANOVA_p = sapply(all_betadisper_results_cailliez, function(x) x$anova_p)
  )
  write.table(global_tests_cailliez, file = file.path(plot_dir_cailliez, "Betadisper_Global_Tests_Cailliez.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (length(all_tukey_results_cailliez) > 0) {
    tukey_c <- do.call(rbind, all_tukey_results_cailliez)
    write.table(tukey_c, file = file.path(plot_dir_cailliez, "Betadisper_TukeyHSD_Cailliez.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (length(all_dispersion_summaries_cailliez) > 0) {
    disp_c <- do.call(rbind, all_dispersion_summaries_cailliez)
    write.table(disp_c, file = file.path(plot_dir_cailliez, "Betadisper_Dispersion_Summary_Cailliez.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# Sqrt
if (length(all_betadisper_results_sqrt) > 0) {
  global_tests_sqrt <- data.frame(
    Time_Period = names(all_betadisper_results_sqrt),
    Correction_Method = "Sqrt",
    Permutation_F = sapply(all_betadisper_results_sqrt, function(x) x$perm_f),
    Permutation_p = sapply(all_betadisper_results_sqrt, function(x) x$perm_p),
    ANOVA_F = sapply(all_betadisper_results_sqrt, function(x) x$anova_f),
    ANOVA_p = sapply(all_betadisper_results_sqrt, function(x) x$anova_p)
  )
  write.table(global_tests_sqrt, file = file.path(plot_dir_sqrt, "Betadisper_Global_Tests_Sqrt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (length(all_tukey_results_sqrt) > 0) {
    tukey_s <- do.call(rbind, all_tukey_results_sqrt)
    write.table(tukey_s, file = file.path(plot_dir_sqrt, "Betadisper_TukeyHSD_Sqrt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (length(all_dispersion_summaries_sqrt) > 0) {
    disp_s <- do.call(rbind, all_dispersion_summaries_sqrt)
    write.table(disp_s, file = file.path(plot_dir_sqrt, "Betadisper_Dispersion_Summary_Sqrt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# =========================================================================
# STEP 7: Create Final Collated Master Table
# =========================================================================

cat("\n=== Creating Final Collated Master Table ===\n")

# 1. Start with the base Wilcoxon stats
final_collated <- as.data.frame(final_stats_df)

# Helper function to create a sorted key for merging
create_merge_key <- function(g1, g2) {
  apply(cbind(as.character(g1), as.character(g2)), 1, function(x) paste(sort(x), collapse = "_"))
}
final_collated$MergeKey <- create_merge_key(final_collated$Group1, final_collated$Group2)

# 2. Merge Cailliez Results
if (exists("global_tests_cailliez")) {
  global_c_sub <- global_tests_cailliez[, c("Time_Period", "ANOVA_p", "Permutation_p")]
  names(global_c_sub) <- c("Time_Period", "Cailliez_Global_ANOVA_p", "Cailliez_Global_Perm_p")
  final_collated <- merge(final_collated, global_c_sub, by = "Time_Period", all.x = TRUE)
} else {
  final_collated$Cailliez_Global_ANOVA_p <- NA; final_collated$Cailliez_Global_Perm_p <- NA
}

if (length(all_tukey_results_cailliez) > 0) {
  tukey_c <- do.call(rbind, all_tukey_results_cailliez)
  tukey_c$MergeKey <- sapply(tukey_c$comparison, function(x) {
    parts <- unlist(strsplit(x, "-")); paste(sort(parts), collapse = "_")
  })
  tukey_c_sub <- tukey_c[, c("Time_Period", "MergeKey", "p adj")]
  names(tukey_c_sub)[3] <- "Cailliez_TukeyHSD_p_adj"
  final_collated <- merge(final_collated, tukey_c_sub, by = c("Time_Period", "MergeKey"), all.x = TRUE)
} else {
  final_collated$Cailliez_TukeyHSD_p_adj <- NA
}

# 3. Merge Sqrt Results
if (exists("global_tests_sqrt")) {
  global_s_sub <- global_tests_sqrt[, c("Time_Period", "ANOVA_p", "Permutation_p")]
  names(global_s_sub) <- c("Time_Period", "Sqrt_Global_ANOVA_p", "Sqrt_Global_Perm_p")
  final_collated <- merge(final_collated, global_s_sub, by = "Time_Period", all.x = TRUE)
} else {
  final_collated$Sqrt_Global_ANOVA_p <- NA; final_collated$Sqrt_Global_Perm_p <- NA
}

if (length(all_tukey_results_sqrt) > 0) {
  tukey_s <- do.call(rbind, all_tukey_results_sqrt)
  tukey_s$MergeKey <- sapply(tukey_s$comparison, function(x) {
    parts <- unlist(strsplit(x, "-")); paste(sort(parts), collapse = "_")
  })
  tukey_s_sub <- tukey_s[, c("Time_Period", "MergeKey", "p adj")]
  names(tukey_s_sub)[3] <- "Sqrt_TukeyHSD_p_adj"
  final_collated <- merge(final_collated, tukey_s_sub, by = c("Time_Period", "MergeKey"), all.x = TRUE)
} else {
  final_collated$Sqrt_TukeyHSD_p_adj <- NA
}

# 4. Final Formatting and Export
desired_order <- c(
  "Time_Period", "Group1", "Group2",
  "Median_Group1_1minusF3", "Median_Group2_1minusF3", "More_Homogenous_Region",
  "p", "p.adj", "p.signif",
  "Kruskal_Wallis_p",
  "Cailliez_Global_ANOVA_p", "Cailliez_Global_Perm_p", "Cailliez_TukeyHSD_p_adj",
  "Sqrt_Global_ANOVA_p", "Sqrt_Global_Perm_p", "Sqrt_TukeyHSD_p_adj"
)

cols_to_keep <- intersect(desired_order, names(final_collated))
final_collated <- final_collated[, cols_to_keep]
final_collated <- final_collated[order(final_collated$Time_Period, final_collated$p), ]

output_file_final <- file.path(plot_dir, "final_of3_group_stats.txt")
write.table(final_collated, file = output_file_final, sep = "\t", quote = FALSE, row.names = FALSE)

output_file_final <- file.path(plot_dir, "Supplementary_Data_8.txt")
write.table(final_collated, file = output_file_final, sep = "\t", quote = FALSE, row.names = FALSE)


cat(sprintf("Successfully saved final collated stats to:\n  %s\n", output_file_final))
cat("Analysis Complete!\n")
