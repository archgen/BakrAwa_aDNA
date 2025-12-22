# ==============================================================================
# Setup and Configuration
# ==============================================================================

# Load necessary libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(maps)
library(grid)
library(gridExtra)

# Define Directories
# Note: Ensure these paths exist or create them before running
work_dir      <- "/Users/mpwilliams/Research/BakrAwa/analyses/Scripts/qpWave/"
out_fig_dir   <- file.path(work_dir, "Figures")
out_table_dir <- file.path(work_dir, "Tables")

# Set working directory
setwd(work_dir)

# ==============================================================================
# Data Loading and Metadata Processing
# ==============================================================================

# Load main data
# Base path defined for reference, though specific file is loaded below
base_path <- "./rds/qpWave_regionalpairwise_Diverse15Right_1240K_public_v50.0_BakrAwa_pmd3_Cayonu_SouthernArc__indsLabels_Transtigris"
combined_data <- readRDS(file = "./rds/qpWave_regionalpairwise_PROCESSED.rds")

# Load metadata map
map_file_path <- "/Users/mpwilliams/Research/BakrAwa/analyses/Scripts/smartPCA/Figures/pca_group_shape_color_map.txt"
map_file <- fread(map_file_path)

# Check for missing groups (Diagnostic)
missing_groups <- unique(combined_data[!(combined_data$PCA_RegionalGroup %in% map_file$PCA_Group)]$PCA_RegionalGroup)
# print(missing_groups) # Uncomment to view

# Define new rows for missing regions
new_rows <- data.table(
  PCA_Group = c(
    "cAnatolia_ChL", "cAnatolia_BA", "cAnatolia_IA", "cAnatolia_Hist",
    "swAnatolia_BA", "swAnatolia_IA",
    "nwAnatolia_ChL",
    "nAnatolia_Hist", "nAnatolia_BA", "nAnatolia_ChL"
  ),
  Region = c(
    "cAnatolia", "cAnatolia", "cAnatolia", "cAnatolia",
    "swAnatolia", "swAnatolia",
    "nwAnatolia",
    "nAnatolia", "nAnatolia", "nAnatolia"
  ),
  shape = c(
    17, 18, 19, 20,
    11, 12,
    5,
    6, 7, 8
  ),
  colour = c(
    "#C71585", "#C71585", "#C71585", "#C71585",
    "#FF4500", "#FF4500",
    "#DC143C",
    "#B22222", "#B22222", "#B22222"
  )
)

# Update map file
map_file <- rbind(map_file, new_rows)

# Verify update (Diagnostic)
remaining_missing <- unique(combined_data[!(combined_data$PCA_RegionalGroup %in% map_file$PCA_Group)]$PCA_RegionalGroup)
# print(remaining_missing) # Uncomment to view

# ==============================================================================
# Statistical Methodology: Non-parametric Permutation Test
# ==============================================================================
#' **Hypothesis Testing for Clade Overrepresentation**
#' 
#' This analysis employs a non-parametric permutation test to determine if individuals 
#' from Bakr Awa form a qpWave clade with individuals from surrounding regions more 
#' frequently than expected by chance.
#'
#' **Test Statistic (T):** 
#' The count of "significant" pairings (where Benjamini-Hochberg corrected p-value > 0.05) 
#' between Bakr Awa individuals and a specific meta-region. This measures enrichment.
#'
#' **Null Hypothesis (H0):** 
#' There is no association between meta-region and the probability of a pairing being 
#' significant. Significant labels are randomly distributed across meta-regions, 
#' proportional to region size.
#'
#' **Alternative Hypothesis (H1):** 
#' There is an overrepresentation (enrichment) of significant pairings in the meta-region.
#'
#' **Method:** 
#' We simulate the null distribution by randomly shuffling (permuting) meta-region labels 
#' while preserving the total count of significant pairs and region sizes. This accounts 
#' for unequal group sizes and potential ancestry correlations.

# ==============================================================================
# Definitions and Helper Functions
# ==============================================================================

# Define meta-region mapping
meta_map <- list(
  C_N_NW_SW_Anatolia = c(
    "cAnatolia_ChL", "cAnatolia_BA", "swAnatolia_BA", "cAnatolia_IA", 
    "nwAnatolia_ChL", "swAnatolia_IA", "cAnatolia_Hist", "nAnatolia_Hist", 
    "nAnatolia_BA", "nAnatolia_ChL"
  ),
  S_Levant = c("sLevant_BA", "sLevant_IA"),
  N_Levant_S_E_Anatolia = c(
    "nLevant_ChL", "nLevant_BA", "eAnatolia_BA", "eAnatolia_ChL",
    "sAnatolia_BA", "sAnatolia_IA"
  ),
  SE_Anatolia_N_Mesopotamia = c(
    "seAnatolia_BA", "nMesopotamia_BA", "seAnatolia_Hist", "seAnatolia_ChL"
  ),
  Caucasus = c("sCaucasus_IA", "sCaucasus_ChL", "sCaucasus_BA", "sCaucasus_Hist"),
  Zagros = c("nwZagros_IA", "cZagros_ChL", "nwZagros_BA")
)

# Create lookup table
rev_map <- data.table(
  PCA_RegionalGroup = unlist(meta_map),
  meta_region = rep(names(meta_map), lengths(meta_map))
)

# Define all meta-regions
all_meta_regions <- names(meta_map)

#' Compute Permutation P-value
#' 
#' Calculates a one-sided p-value for enrichment of significant pairs in a target region.
#' @param dt Data.table containing 'sig' (boolean) and 'meta_region' columns.
#' @param target_region String name of the region to test.
#' @param B Integer, number of permutations.
#' @param seed Integer, random seed for reproducibility.
perm_test_enrichment <- function(dt, target_region, B = 10000, seed = 1) {
  if (nrow(dt) == 0) return(NA_real_)
  
  set.seed(seed)
  
  # Observed k for the target region
  region_dt <- dt[meta_region == target_region]
  if (nrow(region_dt) == 0) return(NA_real_)
  
  k_observed <- sum(region_dt$sig)
  
  # Permutation test
  meta_labels <- dt$meta_region
  perm_k <- replicate(B, {
    meta_perm <- sample(meta_labels)
    sum(dt$sig & meta_perm == target_region)
  })
  
  # One-sided p-value (upper tail)
  p_value <- (sum(perm_k >= k_observed) + 1) / (B + 1)
  return(p_value)
}

#' Compute Enrichment Statistics
#' 
#' Generates a summary table of pair counts, percentages, and permutation p-values.
#' @param dt Data.table containing analysis data.
#' @param q Numeric, significance threshold for p_corr.
#' @param B Integer, number of permutations.
#' @param seed Integer, random seed.
#' @param include_pval Boolean, whether to run the permutation test.
compute_enrichment_stats <- function(dt, q = 0.05, B = 10000, seed = 1, include_pval = TRUE) {
  # Copy data to avoid modifying original by reference
  dt <- copy(dt)
  
  # Merge meta_region info
  setkey(rev_map, PCA_RegionalGroup)
  setkey(dt, PCA_RegionalGroup)
  dt <- dt[rev_map, nomatch = 0]
  dt[is.na(meta_region), meta_region := "Other"]
  
  # Define significance (p_corr > q)
  dt[, sig := p_corr > q]
  
  # Compute stats per Pop1_sample_ID and meta_region
  stats_dt <- dt[, .(
    total_pairs = .N,
    sig_pairs = sum(sig),
    perc_sig = ifelse(.N > 0, (sum(sig) / .N), NA_real_)
  ), by = .(Pop1_sample_ID, meta_region)]
  
  # Pivot to wide format
  wide_stats <- dcast(
    stats_dt, 
    Pop1_sample_ID ~ meta_region, 
    value.var = c("total_pairs", "sig_pairs", "perc_sig"), 
    fill = 0
  )
  
  # Calculate permutation p-values
  if (include_pval) {
    pval_dt <- dt[, {
      if (sum(sig) == 0) {
        pvals <- rep(NA_real_, length(all_meta_regions))
      } else {
        pvals <- sapply(all_meta_regions, function(region) {
          perm_test_enrichment(.SD, target_region = region, B = B, seed = seed)
        })
      }
      data.table(meta_region = all_meta_regions, pval = pvals)
    }, by = Pop1_sample_ID]
    
    pval_wide <- dcast(pval_dt, Pop1_sample_ID ~ meta_region, value.var = "pval")
    wide_stats <- merge(wide_stats, pval_wide, by = "Pop1_sample_ID", all = TRUE)
  }
  
  # Clean up percentages where total_pairs is 0
  for (region in c(all_meta_regions, "Other")) {
    total_col <- paste0("total_pairs_", region)
    perc_col <- paste0("perc_sig_", region)
    if (total_col %in% names(wide_stats)) {
      wide_stats[get(total_col) == 0, (perc_col) := NA_real_]
    }
  }
  
  setorder(wide_stats, Pop1_sample_ID)
  return(wide_stats)
}

# ==============================================================================
# Analysis 1: Chalcolithic - Early Bronze Age (ChL-EBA)
# ==============================================================================

# Filter data: 4200 < Date < 7700 BP, SNPs > 30k, Transtigris
combined_data_bp4k <- combined_data[
  (Date_mean_BP > 4200 & Date_mean_BP < 7700) & 
    Nr_SNPs > 30000 & 
    grepl("Transtigris_", Pop1)
]

# Adjust p-values
combined_data_bp4k[, p_corr := p.adjust(p.value, method = "BH"), by = Pop1_sample_ID]

# Compute stats
summary_table_ChLEBA <- compute_enrichment_stats(
  combined_data_bp4k, q = 0.05, B = 10000, seed = 1, include_pval = TRUE
)

# Preview specific columns
print(summary_table_ChLEBA[, c(1, (ncol(summary_table_ChLEBA) - 5):ncol(summary_table_ChLEBA)), with = FALSE])

# Save results
fwrite(
  summary_table_ChLEBA, 
  file = file.path(out_table_dir, "BakrAwa_enrichment_stats_ChL-EBA.txt"), 
  sep = "\t", quote = FALSE, na = "NA"
)

# ==============================================================================
# Analysis 2: Early Bronze Age - Late Bronze Age (EBA-LBA)
# ==============================================================================

# Filter data: 3400 <= Date <= 4200 BP
combined_data_bp3k4k <- combined_data[
  (Date_mean_BP >= 3400 & Date_mean_BP <= 4200) & 
    Nr_SNPs > 30000 & 
    grepl("Transtigris_", Pop1)
]

# Adjust p-values
combined_data_bp3k4k[, p_corr := p.adjust(p.value, method = "BH"), by = Pop1_sample_ID]

# Compute stats
summary_table_EBALBA <- compute_enrichment_stats(
  combined_data_bp3k4k, q = 0.05, B = 10000, seed = 1, include_pval = TRUE
)

# Preview specific columns
print(summary_table_EBALBA[, c(1, (ncol(summary_table_EBALBA) - 5):ncol(summary_table_EBALBA)), with = FALSE])

# Save results
fwrite(
  summary_table_EBALBA, 
  file = file.path(out_table_dir, "BakrAwa_enrichment_stats_EBA-LBA.txt"), 
  sep = "\t", quote = FALSE, na = "NA"
)

# ==============================================================================
# Analysis 3: Iron Age - Historical (IA-Hist)
# ==============================================================================

# Filter data: Date < 3400 BP
combined_data_bpL3k <- combined_data[
  (Date_mean_BP < 3400) & 
    Nr_SNPs > 30000 & 
    grepl("Transtigris_", Pop1)
]

# Adjust p-values
combined_data_bpL3k[, p_corr := p.adjust(p.value, method = "BH"), by = Pop1_sample_ID]

# Compute stats
summary_table_IAHist <- compute_enrichment_stats(
  combined_data_bpL3k, q = 0.05, B = 10000, seed = 1, include_pval = TRUE
)

# Preview specific columns
print(summary_table_IAHist[, c(1, (ncol(summary_table_IAHist) - 5):ncol(summary_table_IAHist)), with = FALSE])

# Save results
fwrite(
  summary_table_IAHist, 
  file = file.path(out_table_dir, "BakrAwa_enrichment_stats_IA-Hist.txt"), 
  sep = "\t", quote = FALSE, na = "NA"
)
