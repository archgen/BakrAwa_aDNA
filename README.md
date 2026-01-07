# BakrAwa_aDNA
Custom scripts to perform downstream analyses of qpWave and outgroup f3-statistic tests. 

## Analysis Scripts

This repository contains the R code used to generate the statistical results and figures.

### File Descriptions

*   **`qpWave_Permutation_Enrichment_Test.R`**:
    *   **Purpose:** Tests for significant clade enrichment between the target individual and individuals from the surrounding meta-regions.
    *   **Input:** Processed qpWave RDS files and metadata maps.
    *   **Output:** Enrichment statistics tables for ChL-EBA, EBA-LBA, and IA-Hist periods.

*   **`Bayesian_GLM_PopStruct_Correction.R`**:
    *   **Purpose:** Runs a Bayesian GLM (MCMCglmm) to predict significant connections to the surrounding meta-regions while controlling for genetic structure.
    *   **Usage:** Designed to be run via command line (e.g., `Rscript Bayesian_GLM_Kinship_Correction.R <Individual_ID>`).
    *   **Output:** Diagnostic plots (Trace, PPC, ROC) and summary tables of Log-Odds and genetic variance.

*   **`Within_Region_F3_Diversity_Analysis.R`**:
    *   **Purpose:** Assesses genetic homogeneity within meta-regions.
    *   **Method:** Uses `1 - F3` as a dissimilarity metric and performs `betadisper` analysis (multivariate homogeneity of groups dispersions).
    *   **Output:** Violin plots, heatmaps, and statistical test results (ANOVA/Permutation).
