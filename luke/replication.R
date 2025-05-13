#!/usr/bin/env Rscript

# --- Libraries ---
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))

perform_replication <- function(annotation, ukbb_pheno, aou_pheno, aou_ancestry) {

    ukbb_model_file <- paste0("burdenEM_example_output/burdenEM_fit_genebass_", ukbb_pheno, "_", annotation, ".rds")
    if(file.exists(ukbb_model_file)){
        ukbb_model <- readRDS(ukbb_model_file)
        ukbb_model$df <- ukbb_model$df %>% dplyr::filter(!is.na(gene_id))
    } else {
        stop(paste("UKBB model file not found:", ukbb_model_file))
    }

    aou_model_file <- paste0("burdenEM_example_output/burdenEM_fit_aou_", aou_ancestry, "_", aou_pheno, "_", annotation, ".rds")
    if(file.exists(aou_model_file)){
        aou_model <- readRDS(aou_model_file)
        aou_model$df <- aou_model$df %>% dplyr::filter(!is.na(gene))
    } else {
        stop(paste("AoU gene data file not found:", aou_model_file))
    }

    # --- Compute posterior-mean effect sizes for UKBB gene data ---
    source("R/model.R") 
    ukbb_model$df$posterior_mean <- posterior_expectation2(ukbb_model)

    # --- 5. Merge AoU Gene Data with UKBB Gene Data ---
    if (!("gene_id" %in% names(ukbb_model$df))) stop("Column 'gene_id' not found in UKBB gene data.")
    if (!("gene" %in% names(aou_model$df))) stop("Column 'gene' (expected to contain Ensembl IDs) not found in AoU gene data.")
    
    merged_gene_data <- inner_join(ukbb_model$df, aou_model$df, by = c("gene_id" = "gene"), suffix = c("", ".aou"))

    # --- 6. Adjust UKBB posterior means, bin, and summarize with AoU effects ---
    if (!("posterior_mean" %in% names(merged_gene_data))) {
        stop("Column 'posterior_mean' (from UKBB) not found in merged_gene_data.")
    }
    if (!("effect_estimate.aou" %in% names(merged_gene_data))) {
        stop("Column 'effect_estimate.aou' not found in merged_gene_data. Check merge suffixes.")
    }

    breaks <- ukbb_model$component_endpoints
    breaks[1] <- -Inf
    breaks[length(breaks)] <- Inf

    merged_gene_data <- merged_gene_data %>%
        mutate(
            burden_score_multiplier = sqrt(burden_score.aou / burden_score),
            burden_score_multiplier = ifelse(is.finite(burden_score_multiplier), burden_score_multiplier, 1),
            adjusted_posterior_mean = posterior_mean * burden_score_multiplier,
            bin = cut(adjusted_posterior_mean, breaks = breaks)
        )

    combined_adjusted_bin_stats <- merged_gene_data %>%
        group_by(bin) %>%
        summarise(
            count = n(), 
            expected = mean(adjusted_posterior_mean, na.rm = TRUE),
            mean_aou = mean(effect_estimate.aou, na.rm = TRUE),
            mean_ukbb = mean(effect_estimate, na.rm = TRUE),
            .groups = 'drop'
        )

    print(combined_adjusted_bin_stats)
}

# --- Script Execution ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript replication.R <annotation> <ukbb_pheno> <aou_pheno> <aou_ancestry>", call. = FALSE)
}

perform_replication(args[1], args[2], args[3], args[4])
