#!/usr/bin/env Rscript
#conda_lib <- "/storage/work/mkw5910/.conda/envs/admixtools-env/lib/R/library"
#.libPaths(c(conda_lib, .libPaths()))

# 01_worker.R
# Usage: Rscript 01_worker.R <Individual_ID>


# =========================================================================
# SETUP: Load Packages and Define Paths
# =========================================================================
# Clear any stuck graphics devices
graphics.off()

suppressPackageStartupMessages({
  library(MCMCglmm)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(coda)
  library(ggmcmc)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(parallel)
  library(pROC)
  library(bayesplot)
})

# =========================================================================
# 1. PARSE COMMAND LINE ARGUMENTS
# =========================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No individual ID provided. Usage: Rscript 01_worker.R <Individual_ID>")
}

# The specific individual to process in this run
unique_inds <- args[1] 


# =========================================================================
# STEP 1: Data Loading and Preprocessing
# =========================================================================
# Define paths

base_path <- "/storage/group/cdh5313/default/mkw5910/aDNA/analyses/qpwave/bakrawa/"
qpwave_path <- file.path(base_path, "/rds/qpWave_regionalpairwise_PROCESSED.rds")
f3_path <- file.path(base_path, "/rds/oF3_OldAfrica__qpWaveRegionalInds.rds")

output_dir <- file.path(base_path, "/out/BayesianGLM_Results")
plot_dir <- file.path(output_dir, "/out/Diagnostic_Plots")

# Create output directories
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# Load Data
combined_data <- readRDS(file = qpwave_path)
oF3_OldAfrica__qpWaveRegionalInds <- data.table(readRDS(file = f3_path))

# Clean F3 Data
oF3_OldAfrica__qpWaveRegionalInds <- oF3_OldAfrica__qpWaveRegionalInds[
  pop2 != "NAN" & pop3 != "NAN" & pop2 != pop3
]

# Clean Combined Data
combined_data <- combined_data[!is.na(Pop1_sample_ID)]
setDT(combined_data)
combined_data[, p_corr := p.adjust(p.value, method = "BH"), by = Pop1_sample_ID]
combined_data[, sig := p_corr > 0.05]

# Calculate F3 similarity
setDT(oF3_OldAfrica__qpWaveRegionalInds)
oF3_OldAfrica__qpWaveRegionalInds[, f3_sim := 1 - est]

# -------------------------------------------------------------------------
# STEP 2: Standardize Labels (Global Kinship Matrix Prep)
# -------------------------------------------------------------------------

# Create mapping: Pop2 <-> sample_ID
id_map <- unique(combined_data[, .(Pop2, sample_ID)])

# Update F3 labels to match Combined Data
for (i in 1:nrow(id_map)) {
  target_label <- id_map$Pop2[i]
  target_id    <- id_map$sample_ID[i]

  if (is.na(target_id) || target_id == "") next

  # Update pop2 and pop3 in F3 table
  oF3_OldAfrica__qpWaveRegionalInds[grepl(target_id, pop2, fixed = TRUE), pop2 := target_label]
  oF3_OldAfrica__qpWaveRegionalInds[grepl(target_id, pop3, fixed = TRUE), pop3 := target_label]
}

# Reshape to Matrix (Long to Wide)
f3_wide <- dcast(oF3_OldAfrica__qpWaveRegionalInds, pop2 ~ pop3, value.var = "f3_sim")
f3_matrix <- as.matrix(f3_wide[, -1])
rownames(f3_matrix) <- f3_wide$pop2

# Handle Missing Diagonal Values & Symmetry
max_sim <- max(f3_matrix, na.rm = TRUE)
if (any(is.na(f3_matrix))) {
  f3_matrix[is.na(f3_matrix)] <- max_sim
}
f3_matrix[lower.tri(f3_matrix)] <- t(f3_matrix)[lower.tri(f3_matrix)]

# Ensure Positive Semi-Definiteness (Global Matrix)
f3_psd <- Matrix::nearPD(f3_matrix, corr = FALSE, keepDiag = TRUE)$mat
f3_psd <- as.matrix(f3_psd)

# -------------------------------------------------------------------------
# STEP 3: Metadata Mapping
# -------------------------------------------------------------------------
meta_map <- list(
  C_N_NW_SW_Anatolia = c("cAnatolia_ChL", "cAnatolia_BA", "swAnatolia_BA", "cAnatolia_IA", "nwAnatolia_ChL", "swAnatolia_IA",
                         "cAnatolia_Hist", "nAnatolia_Hist", "nAnatolia_BA", "nAnatolia_ChL"),
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
setkey(rev_map, PCA_RegionalGroup)
setkey(combined_data, PCA_RegionalGroup)

combined_data_MD <- merge(combined_data, rev_map, by = "PCA_RegionalGroup", all.x = TRUE)
combined_data_MD <- combined_data_MD[!is.na(meta_region)]

# =========================================================================
# STEP 4: Analysis Loop (Individuals x Time Periods)
# =========================================================================
time_periods <- c("TP1_4200_7700BP", "TP2_3400_4200BP", "TP3_LT_3400BP")

# Initialize Storage
main_results_list <- list()
supp_results_list <- list(
  "TP1_4200_7700BP" = list(),
  "TP2_3400_4200BP" = list(),
  "TP3_LT_3400BP"   = list()
)


#- MCMC Simulation Parameters
#-- Number of cores to use
n_cores = 20
#-- Number of chains to run
n_chains = 4


# #- Run A
# nitt_val <- 200000
# burnin_val <- 20000
# thin_val <-40

#- Run B
nitt_val <- 5000000
burnin_val <- 500000
thin_val <-2500

# #- Quick Test
#nitt_val <- 600
#burnin_val <- 100
#thin_val <-40


# MCMC Parameters (Prior)
#- Original Prior [COMMENT OUT IF USING DYNAMIC K ADJUSTED PRIOR]
# prior <- list(
#   R = list(V = 1, fix = 1),
#   G = list(G1 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 1000))
# )


cat("Starting Analysis Loop with Parallel Chains...\n")

for (ind in unique_inds) {

  # Initialize list to hold plots for THIS individual across ALL time periods
  ind_plots <- list()

  for (tp in time_periods) {

    # 4.1 Subset Data
    subset_dt <- switch(tp,
                        "TP1_4200_7700BP" = combined_data_MD[Pop1_sample_ID == ind & (Date_mean_BP > 4200 & Date_mean_BP < 7700) & Nr_SNPs > 30000],
                        "TP2_3400_4200BP" = combined_data_MD[Pop1_sample_ID == ind & (Date_mean_BP >= 3400 & Date_mean_BP <= 4200) & Nr_SNPs > 30000],
                        "TP3_LT_3400BP"   = combined_data_MD[Pop1_sample_ID == ind & (Date_mean_BP < 3400) & Nr_SNPs > 30000]
    )

    # 4.2 Check Data Sufficiency (Initial)
    if (nrow(subset_dt) < 5) {
      cat(sprintf("Skipping %s - %s: Insufficient data (N=%d)\n", ind, tp, nrow(subset_dt)))
      next
    }

    # 4.3 Check for Single Significant Region
    sig_regions <- unique(subset_dt[sig == TRUE, meta_region])
    if (length(sig_regions) == 1) {
      cat(sprintf("Skipping %s - %s: Only one region (%s) has significant connections.\n",
                  ind, tp, sig_regions))
      next
    }

    cat(sprintf("Processing %s - %s (N=%d)...\n", ind, tp, nrow(subset_dt)))

    # 4.4 Prepare Data Variables
    subset_dt$meta_region <- as.factor(subset_dt$meta_region)
    subset_dt$sig_num <- as.integer(subset_dt$sig)

    # ---------------------------------------------------------------------
    # 4.5 ALIGN KINSHIP MATRIX AND DATA (CORRECTED ORDER)
    # ---------------------------------------------------------------------

    # 1. CLEANUP FIRST: Remove NAs before doing anything else
    subset_dt <- subset_dt[!is.na(sig_num) & !is.na(meta_region) & !is.na(Pop2)]

    # 2. Check Data Sufficiency after cleanup
    if (nrow(subset_dt) < 5) {
      cat(sprintf("Skipping %s - %s: Insufficient data after NA cleanup (N=%d)\n", ind, tp, nrow(subset_dt)))
      next
    }

    # 3. Identify IDs present in BOTH the CLEANED data and the global matrix
    individuals_in_dt <- unique(subset_dt$Pop2)
    available_in_matrix <- rownames(f3_psd)

    # 4. Create a SORTED list of valid IDs
    valid_ids <- sort(intersect(individuals_in_dt, available_in_matrix))

    # 5. Check if we have enough valid IDs left
    if (length(valid_ids) < 3) {
      cat(sprintf("Skipping %s - %s: Too few valid IDs matching kinship matrix (N=%d)\n", ind, tp, length(valid_ids)))
      next
    }

    # 6. Subset the Data to only include valid IDs
    subset_dt <- subset_dt[Pop2 %in% valid_ids]

    # 7. Subset and Order the Kinship Matrix
    f3_sub <- f3_psd[valid_ids, valid_ids]

    # 8. Assign Factor Levels to Data
    subset_dt$id <- factor(subset_dt$Pop2, levels = valid_ids)

    # 9. Convert to data.frame and set row names (Fixes ROC error)
    subset_df <- as.data.frame(subset_dt)
    rownames(subset_df) <- 1:nrow(subset_df)

    # 10. Invert Matrix
    eigvals <- eigen(f3_sub, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigvals) < 1e-8) {
      f3_sub <- f3_sub + diag(1e-6, nrow(f3_sub))
    }
    f3inv_sparse <- as(solve(f3_sub), "dgCMatrix")

    # ---------------------------------------------------------------------
    # 4.6 Run MCMCglmm (MULTIPLE CHAINS)
    # ---------------------------------------------------------------------
    # Model Formula
    fixed_formula <- sig_num ~ meta_region

    subset_df$meta_region <- droplevels(subset_df$meta_region)


    K <- ncol(model.matrix(fixed_formula, data = subset_df))

    # 2. Define the Prior using this dynamic K
    prior <- list(
      # B: Fixed Effects Prior (The Fix for Infinite Odds Ratios)
      # Uses the dynamic K to create the correct sized matrix.
      # V = 5 penalizes the infinite values (10^60) back to reality.
      B = list(mu = rep(0, K), V = diag(K) * 5),

      # R: Residual Variance (Fixed for binary/categorical)
      R = list(V = 1, fix = 1),

      # G: Random Effects (Genetic Matrix)
      # Parameter expanded to help convergence
      G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))
    )

    # Run N chains in parallel
    chains <- mclapply(1:n_chains, function(i) {
      tryCatch({
        MCMCglmm(
          fixed = fixed_formula,
          random = ~ id,
          family = "categorical",
          data = subset_df,
          ginverse = list(id = f3inv_sparse),
          prior = prior,
          nitt = nitt_val, burnin = burnin_val, thin = thin_val,
          pr = TRUE, # <--- CRITICAL: Saves Random Effects for Confounding Check
          verbose = FALSE
        )
      }, error = function(e) return(NULL))
    }, mc.cores = n_cores)

    # Check if chains failed
    if (any(sapply(chains, is.null))) {
      cat(sprintf("One or more chains failed for %s - %s. Skipping.\n", ind, tp))
      next
    }

    # ---------------------------------------------------------------------
    # 4.7 Process Chains & Diagnostics
    # ---------------------------------------------------------------------

    # Extract Sol (Fixed Effects) and VCV (Random Effects) from both chains
    sol_list <- lapply(chains, function(m) m$Sol)
    vcv_list <- lapply(chains, function(m) m$VCV)

    # Create mcmc.list objects for Coda diagnostics
    sol_mcmc_list <- as.mcmc.list(lapply(sol_list, as.mcmc))
    vcv_mcmc_list <- as.mcmc.list(lapply(vcv_list, as.mcmc))

    # Combine chains for Inference (Boosting Effective Sample Size)
    combined_Sol <- do.call(rbind, sol_list)
    combined_VCV <- do.call(rbind, vcv_list)

    # Convert combined to mcmc object for summary functions
    combined_Sol_mcmc <- as.mcmc(combined_Sol)
    combined_VCV_mcmc <- as.mcmc(combined_VCV)

    # Sanitize file names
    ind_clean <- gsub("[^[:alnum:]_\\-]", "", ind)
    tp_clean  <- gsub("[^[:alnum:]_\\-]", "", tp)

    # =====================================================================
    # NEW: CHECK CORRELATION BETWEEN FIXED AND RANDOM EFFECTS
    # =====================================================================

    # 1. Identify columns in Sol that correspond to Random Effects (id)
    all_col_names <- colnames(combined_Sol)
    rand_eff_cols <- grep("^id\\.", all_col_names)

    if (length(rand_eff_cols) > 0) {
      # 2. Extract Posterior Means for Random Effects
      rand_eff_matrix <- combined_Sol[, rand_eff_cols]
      rand_eff_means <- colMeans(rand_eff_matrix)

      # 3. Clean up names to match Data ID
      rand_eff_ids <- gsub("^id\\.", "", names(rand_eff_means))

      # 4. Create a Dataframe for Plotting
      re_df <- data.frame(Pop2 = rand_eff_ids, Random_Effect_Mean = rand_eff_means)

      # Merge with the subset_dt to get the Meta_Region for each individual
      re_df$Meta_Region <- subset_dt$meta_region[match(re_df$Pop2, subset_dt$Pop2)]

      # 5. Generate Boxplot: Random Effects vs Fixed Effects
      confound_plot <- ggplot(re_df, aes(x = Meta_Region, y = Random_Effect_Mean, fill = Meta_Region)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        labs(title = sprintf("Confounding Check: %s | %s", ind, tp),
             subtitle = "Random Effect Estimates (Genetics) by Region",
             y = "Posterior Mean of Random Effect (id)",
             x = "Fixed Effect (Meta Region)")

      ggsave(filename = file.path(plot_dir, sprintf("Confounding_%s_%s.pdf", ind_clean, tp_clean)),
             plot = confound_plot, width = 8, height = 6)

      ind_plots[[length(ind_plots) + 1]] <- confound_plot
    }

    # =====================================================================
    # STANDARD DIAGNOSTICS (FILTERED FOR FIXED EFFECTS)
    # =====================================================================

    # Filter Sol to only keep Fixed Effects for the Trace/Density plots
    # Otherwise ggs() will try to plot hundreds of random effects
    fixed_eff_cols <- grep("^meta_region|^\\(Intercept\\)", all_col_names)
    fixed_Sol_mcmc_list <- sol_mcmc_list[, fixed_eff_cols, drop=FALSE]

    # --- Gelman-Rubin Diagnostic (Convergence Check) ---
    gelman_val <- tryCatch({
      gelman.diag(fixed_Sol_mcmc_list, multivariate = FALSE)$psrf[,1]
    }, error = function(e) return(NA))

    # 1. Trace Plot (Fixed Effects Only)
    ggs_sol <- ggs(fixed_Sol_mcmc_list)
    trace_plot <- ggs_traceplot(ggs_sol) +
      ggtitle(sprintf("Trace (Fixed Effects): %s | %s", ind, tp)) + theme_bw()
    ggsave(filename = file.path(plot_dir, sprintf("Trace_%s_%s.pdf", ind_clean, tp_clean)), plot = trace_plot, width = 8, height = 6)

    # 2. Density Plot (Fixed Effects Only)
    density_plot <- ggs_density(ggs_sol) +
      ggtitle(sprintf("Density (Fixed Effects): %s | %s", ind, tp)) + theme_bw()
    ggsave(filename = file.path(plot_dir, sprintf("Density_%s_%s.pdf", ind_clean, tp_clean)), plot = density_plot, width = 8, height = 6)

    # 3. Effective Sample Size (ESS) - Calculated on COMBINED chains (Fixed Effects)
    # We use the filtered list to avoid calculating ESS for every random effect
    combined_Fixed_Sol_mcmc <- combined_Sol_mcmc[, fixed_eff_cols, drop=FALSE]
    ess_vals <- effectiveSize(combined_Fixed_Sol_mcmc)
    ess_df <- data.frame(Parameter = names(ess_vals), ESS = as.numeric(ess_vals))
    write.table(ess_df, file = file.path(plot_dir, sprintf("ESS_%s_%s.txt", ind_clean, tp_clean)), sep = "\t", row.names = FALSE, quote = FALSE)

    # 4. Gelman-Rubin Output
    gelman_df <- data.frame(Parameter = names(gelman_val), PSRF = as.numeric(gelman_val))
    write.table(gelman_df, file = file.path(plot_dir, sprintf("Gelman_%s_%s.txt", ind_clean, tp_clean)), sep = "\t", row.names = FALSE, quote = FALSE)

    # 5. Heidelberger (on Combined Chain - Fixed Effects)
    heidel_vals <- heidel.diag(combined_Fixed_Sol_mcmc)
    heidel_df <- data.frame(Parameter = rownames(heidel_vals), unclass(heidel_vals))
    write.table(heidel_df, file = file.path(plot_dir, sprintf("Heidel_%s_%s.txt", ind_clean, tp_clean)), sep = "\t", row.names = FALSE, quote = FALSE)

    # ---------------------------------------------------------------------
    # 4.8 Posterior Predictive Check & ROC (Using Combined Estimates)
    # ---------------------------------------------------------------------

    # Use the first model object for X matrix structure
    model_template <- chains[[1]]
    Xmat <- model_template$X

    # Use COMBINED posterior samples for predictions
    # Note: We need the FULL betas (including random effects) for prediction if we were predicting specific individuals,
    # but typically PPC for GLMMs focuses on the fixed effects structure or marginal predictions.
    # However, MCMCglmm's X matrix usually only includes Fixed Effects.
    # The Random Effects are in Z.
    # Standard PPC often uses just X * Beta (Fixed) for the linear predictor if Z is not easily accessible,
    # OR we can reconstruct the full linear predictor.
    # Given the complexity, we will stick to Fixed Effects prediction for the PPC/ROC as a conservative estimate.

    # Extract Fixed Effect Betas only
    betas_fixed <- as.matrix(combined_Sol[, fixed_eff_cols])
    n_draws <- nrow(betas_fixed)

    # PPC
    yrep <- matrix(NA, nrow = n_draws, ncol = nrow(Xmat))
    for (i in 1:n_draws) {
      lin_pred <- as.vector(Xmat %*% betas_fixed[i, ])
      p <- plogis(lin_pred)
      yrep[i, ] <- rbinom(n = nrow(Xmat), size = 1, prob = p)
    }

    used_rows <- rownames(Xmat)

    # --- FIX START ---
    # MCMCglmm appends .1 (or similar) to row names in the design matrix (X).
    # We must strip this suffix to match the original dataframe row names.
    used_rows_clean <- gsub("\\.[0-9]+$", "", used_rows)
    observed_y <- subset_df[used_rows_clean, "sig_num"]
    # --- FIX END ---

    # FIX: Ensure numeric and remove any accidental NAs
    observed_y <- as.numeric(observed_y)
    valid_idx <- !is.na(observed_y)

    if(sum(valid_idx) > 0 && requireNamespace("bayesplot", quietly = TRUE)){
      y_clean <- observed_y[valid_idx]
      yrep_clean <- yrep[, valid_idx, drop=FALSE]

      tryCatch({
        # Added bins=30 to silence the warning
        ppc_plot <- bayesplot::ppc_stat(y = y_clean, yrep = yrep_clean, stat = "mean", binwidth = 0.01) +
          ggtitle(sprintf("PPC (Mean): %s | %s", ind, tp))
        ggsave(filename = file.path(plot_dir, sprintf("PPC_Mean_%s_%s.pdf", ind_clean, tp_clean)), plot = ppc_plot, width = 6, height = 4)
        ind_plots[[length(ind_plots) + 1]] <- ppc_plot
      }, error = function(e) {
        cat(sprintf("PPC Plot Failed for %s: %s\n", ind, e$message))
      })

      # =====================================================================
      # ADDED: PPC Distribution Check (Bars for Binary Data)
      # =====================================================================
      tryCatch({
        cat(sprintf("Attempting PPC Distribution Plot for %s - %s...\n", ind, tp))

        # Use a random subset of 50 draws
        n_samples_dens <- min(50, nrow(yrep_clean))
        samp_idx <- sample(nrow(yrep_clean), n_samples_dens)

        # IMPORTANT: Use drop=FALSE to maintain matrix structure
        yrep_sub <- yrep_clean[samp_idx, , drop=FALSE]

        # For Binary Data (0/1), ppc_bars is the correct visualization, not density
        # ppc_dens_overlay is for continuous data.
        ppc_dist_plot <- bayesplot::ppc_bars(y = y_clean, yrep = yrep_sub) +
          ggtitle(sprintf("PPC (Distribution): %s | %s", ind, tp))

        ggsave(filename = file.path(plot_dir, sprintf("PPC_Distribution_%s_%s.pdf", ind_clean, tp_clean)),
               plot = ppc_dist_plot, width = 6, height = 4)

        ind_plots[[length(ind_plots) + 1]] <- ppc_dist_plot
        cat("PPC Distribution Plot Saved.\n")

      }, error = function(e) {
        cat(sprintf("PPC Distribution Plot Failed for %s: %s\n", ind, e$message))
      })
      # =====================================================================
    }

    # ROC/AUC
    auc_val <- NA
    pred_matrix <- as.matrix(Xmat %*% t(betas_fixed))

    # FIX: Use rowMeans instead of colMeans to get average prediction per observation
    p_mean <- rowMeans(plogis(pred_matrix))

    if (nrow(Xmat) != length(observed_y)) {
      cat(sprintf("Warning: Row mismatch for %s. Skipping ROC.\n", ind))
    } else {
      if (requireNamespace("pROC", quietly = TRUE)) {
        if (length(unique(observed_y)) == 2) {
          tryCatch({
            roc_obj <- pROC::roc(response = observed_y, predictor = p_mean, quiet = TRUE)
            auc_val <- as.numeric(pROC::auc(roc_obj))
            png(filename = file.path(plot_dir, sprintf("ROC_%s_%s.png", ind_clean, tp_clean)), width = 6*72, height = 4*72, res = 72)
            plot(roc_obj, main = sprintf("ROC (AUC=%.2f): %s | %s", auc_val, ind, tp))
            dev.off()
          }, error = function(e) cat(sprintf("ROC Error: %s\n", e$message)))
        } else {
          cat(sprintf("Skipping ROC for %s: Invariant response.\n", ind))
        }
      }
    }

    # ---------------------------------------------------------------------
    # 4.9 Diagnostic Summary
    # ---------------------------------------------------------------------
    summary_file <- file.path(plot_dir, sprintf("DiagnosticSummary_%s_%s.txt", ind_clean, tp_clean))
    auc_print <- if(is.na(auc_val)) "NA" else round(auc_val, 3)
    cat(
      "ESS (Combined):\n", capture.output(print(ess_df)), "\n\n",
      "Gelman-Rubin PSRF:\n", capture.output(print(gelman_df)), "\n\n",
      "Heidelberger:\n", capture.output(print(heidel_df)), "\n\n",
      "AUC:\n", auc_print, "\n",
      file = summary_file
    )

    # =====================================================================
    # 4.10 Extract Results (MANUAL CALCULATION FROM COMBINED CHAINS)
    # =====================================================================

    # Calculate Statistics on Combined Posterior (Fixed Effects Only)
    combined_Fixed_Sol <- combined_Sol[, fixed_eff_cols, drop=FALSE]
    combined_Fixed_Sol_mcmc <- as.mcmc(combined_Fixed_Sol)

    post_means <- colMeans(combined_Fixed_Sol)
    post_CI <- HPDinterval(combined_Fixed_Sol_mcmc)

    # Calculate pMCMC manually (2 * min(p>0, p<0))
    p_mcmc <- apply(combined_Fixed_Sol, 2, function(x) {
      2 * min(mean(x > 0), mean(x < 0))
    })

    # Genetic Variance (Combined)
    Ve <- pi^2 / 3
    Vg_samples <- combined_VCV[, 1]
    prop_gen_var_samples <- Vg_samples / (Vg_samples + Ve)

    prop_gen_var_mean <- mean(prop_gen_var_samples)
    Vg_mean <- mean(Vg_samples)
    Vg_CI <- HPDinterval(as.mcmc(Vg_samples))

    # Meta Region Names
    meta_regions_names <- colnames(combined_Fixed_Sol)
    meta_regions_names <- gsub("meta_region", "", meta_regions_names)

    df_main <- data.frame(
      Individual_ID = ind,
      Time_Period = tp,
      Meta_Region = meta_regions_names,
      Log_Odds = round(post_means, 3),
      Log_Odds_95CI = paste0("(", round(post_CI[,1], 3), ", ", round(post_CI[,2], 3), ")"),
      Odds_Ratio = round(exp(post_means), 3),
      pMCMC = signif(p_mcmc, 2),
      Prop_Gen_Var = round(prop_gen_var_mean, 3)
    )
    main_results_list[[length(main_results_list) + 1]] <- df_main

    # Supplementary Data
    # Average DIC from chains
    dic_val <- mean(sapply(chains, function(m) m$DIC))
    n_samples <- nrow(Xmat)
    n_regions <- length(unique(subset_dt$meta_region))

    col_param <- paste0(ind, "_Parameter")
    col_value <- paste0(ind, "_Value")

    df_supp <- data.frame(
      Parameter = c("DIC_Avg", "Sample_Size", "N_Meta_Regions", "Gen_Var_Mean",
                    "Gen_Var_95CI", "Prop_Gen_Var", "Gelman_Max_PSRF"),
      Value = c(round(dic_val, 3), n_samples, n_regions, round(Vg_mean, 3),
                paste0("(", round(Vg_CI[1], 3), ", ", round(Vg_CI[2], 3), ")"),
                round(prop_gen_var_mean, 3),
                ifelse(all(is.na(gelman_val)), "NA", round(max(gelman_val, na.rm=TRUE), 3)))
    )
    colnames(df_supp) <- c(col_param, col_value)
    supp_results_list[[tp]][[ind]] <- df_supp

  } # End Time Period Loop

  # Save Dashboard
  if (length(ind_plots) > 0) {
    ind_clean <- gsub("[^[:alnum:]_\\-]", "", ind)
    n_rows <- ceiling(length(ind_plots) / 4)
    h_in <- 4.5 * n_rows
    w_in <- 25
    dashboard_filename <- file.path(plot_dir, sprintf("Dashboard_%s_All_TPs.pdf", ind_clean))
    g <- arrangeGrob(grobs = ind_plots, ncol = 4, top = textGrob(paste("Diagnostic Dashboard:", ind), gp = gpar(fontsize = 15, fontface = "bold")))
    tryCatch({ ggsave(filename = dashboard_filename, plot = g, width = w_in, height = h_in, limitsize = FALSE) }, error = function(e) NULL)
  }
} # End Individual Loop

# =========================================================================
# STEP 5: Save Outputs
# =========================================================================

if (length(main_results_list) > 0) {
  final_main_table <- do.call(rbind, main_results_list)
  write.table(final_main_table, file = file.path(output_dir, sprintf("Main_Publication_Table_%s.txt", ind_clean)), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Main table saved.\n")

  # Gold Standard Extraction
  gold_dt <- as.data.table(final_main_table)
  gold_dt[, c("CI_Lower", "CI_Upper") := tstrsplit(gsub("[()]", "", Log_Odds_95CI), ",")]
  gold_dt[, CI_Lower := as.numeric(CI_Lower)]
  gold_dt[, CI_Upper := as.numeric(CI_Upper)]

  gold_standard_combinations <- gold_dt[pMCMC < 0.05 & Prop_Gen_Var > 0.9 & (CI_Lower > 0 | CI_Upper < 0)]

  if (nrow(gold_standard_combinations) > 0) {
    write.table(gold_standard_combinations[, .(Individual_ID, Time_Period, Meta_Region, Log_Odds, Log_Odds_95CI, pMCMC, Prop_Gen_Var)],
                file = file.path(output_dir, "Gold_Standard_Combinations.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Gold Standard table saved.\n")
  }
}

for (tp in names(supp_results_list)) {
  tp_list <- supp_results_list[[tp]]
  if (length(tp_list) > 0) {
    final_supp_table <- do.call(cbind, tp_list)
    write.table(final_supp_table, file = file.path(output_dir, sprintf("Supplementary_Table__%s_%s.txt", ind_clean, tp)), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

cat("Analysis Complete. Diagnostic plots saved in:", plot_dir, "\n")
