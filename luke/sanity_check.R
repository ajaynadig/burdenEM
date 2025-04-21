# sanity_check.R
# Script to load gene-level data, print head, and calculate summary stats.

# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)

# --- Get Annotation and Phenotype from Command Line ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No annotation provided. Usage: Rscript sanity_check.R <annotation> [phenotype]", call. = FALSE)
}

annotation <- args[1]
phenotype <- ifelse(length(args) > 1, args[2], "50_NA") # Default phenotype
per_allele <- ifelse(length(args) > 2, args[3], FALSE) # Default per_allele

message(paste("Running sanity check for:", annotation, "annotation and", phenotype, "phenotype..."))
# ---------------------------------------------------

# Define the path to the gene-level data file using parsed arguments
gene_data_path <- paste0("burdenEM_example_output/gene_level_data_genebass_", phenotype, "_", annotation, ".rds")

# Check if the file exists before attempting to load
if (!file.exists(gene_data_path)) {
  stop(paste("Error: Gene data file not found at", gene_data_path,
             "\nPlease ensure main.R has been run successfully for this annotation and phenotype."))
}

# Load the RDS file
gene_data <- readRDS(gene_data_path)

# Calculate and print means
mean_burden_score <- mean(gene_data$burden_score, na.rm = TRUE)
mean_burden_score_ld <- mean(gene_data$burden_score_ld, na.rm = TRUE)
mean_burden_score_no_ld <- mean(gene_data$burden_score_no_ld, na.rm = TRUE)

# message("\nMean Burden Scores:")
# message(paste("  Mean burden_score:", round(mean_burden_score, 4)))
# message(paste("  Mean burden_score_ld:", round(mean_burden_score_ld, 4)))
# message(paste("  Mean burden_score_no_ld:", round(mean_burden_score_no_ld, 4)))

# Calculate and print mean/min LD correction factor if column exists
# if ("ld_correction_factor" %in% names(gene_data)) {
#   mean_ld_factor <- mean(gene_data$ld_correction_factor, na.rm = TRUE)
#   min_ld_factor <- min(gene_data$ld_correction_factor, na.rm = TRUE)
  # message("\nLD Correction Factor Stats:")
  # message(paste("  Mean ld_correction_factor (burden_score_ld / burden_score_no_ld):", round(mean_ld_factor, 4)))
  # message(paste("  Min ld_correction_factor:", round(min_ld_factor, 4)))

# Calculate and print additional means
mean_gene_intercept <- mean(gene_data$gene_intercept, na.rm = TRUE)
mean_sq_effect_estimate <- mean(gene_data$effect_estimate^2, na.rm = TRUE)
mean_sq_effect_se <- mean(gene_data$effect_se^2, na.rm = TRUE)

# multiplier = ifelse(per_allele, gene_data$burden_score, 1)
implied_h2 <- sum(gene_data$effect_estimate^2/gene_data$effect_se^2 - 1) * mean_sq_effect_se
implied_h2_se <- sqrt(3 * nrow(gene_data)) * mean_sq_effect_se

message("\nAdditional Mean Values:")
# message(paste("  Mean gene_intercept:", signif(mean_gene_intercept, 4)))
message(paste("  Mean squared effect_estimate:", signif(mean_sq_effect_estimate, 4)))
message(paste("  Mean squared effect_se:", signif(mean_sq_effect_se, 4)))
message(paste("  Implied H2:", signif(implied_h2, 4)))
message(paste("  Implied H2 SE:", signif(implied_h2_se, 4)))

message("\nSanity check script finished.")
