#!/usr/bin/env Rscript

# --- Libraries ---
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
source("R/model.R") 

# Helper function for binomial confidence intervals using the Clopper-Pearson exact method
binomial_ci <- function(x, n, conf_level = 0.95) {
    alpha <- 1 - conf_level
    lower <- qbeta(alpha/2, x, n - x + 1)
    upper <- qbeta(1 - alpha/2, x + 1, n - x)
    data.frame(estimate = x/n, lower = lower, upper = upper)
}

compute_posterior_mean_replication <- function(ukbb_model, aou_model, convert_to_per_allele = FALSE) {

    if(convert_to_per_allele){
        effect_size_fn <- function(x, row) x / sqrt(row$burden_score)
        effect_sq_fn <- function(x, row) x^2 / row$burden_score
    } else {
        effect_size_fn <- function(x) x
        effect_sq_fn <- function(x) x^2
    }
    
    # --- Compute posterior-mean effect sizes for UKBB gene data ---
    ukbb_model$df$posterior_mean <- posterior_expectation2(ukbb_model, function_to_integrate = effect_size_fn)
    ukbb_model$df$posterior_mean_sq <- posterior_expectation2(ukbb_model, function_to_integrate = effect_sq_fn)
    ukbb_model$df$posterior_variance <- ukbb_model$df$posterior_mean_sq - (ukbb_model$df$posterior_mean^2)

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

    # Calculate adjusted posterior mean and variance first
    # merged_gene_data <- merged_gene_data %>% 
    #     mutate(
    #         burden_score_multiplier = sqrt(burden_score.aou / burden_score),
    #         burden_score_multiplier = ifelse(is.finite(burden_score_multiplier), burden_score_multiplier, 1),
    #         adjusted_posterior_mean = posterior_mean * burden_score_multiplier,
    #         adjusted_posterior_variance = posterior_variance * (burden_score_multiplier^2)
    #     )

    # --- Define breaks based on quantiles of adjusted_posterior_mean ---
    if (nrow(merged_gene_data) > 0 && sum(!is.na(merged_gene_data$posterior_mean)) > 0) {
        quantile_probs <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
        apm_values <- merged_gene_data$posterior_mean[!is.na(merged_gene_data$posterior_mean)]
        if (length(apm_values) > 0) {
            quantile_values <- quantile(apm_values, probs = quantile_probs, na.rm = FALSE, type = 7)
            breaks <- unique(c(-Inf, quantile_values, Inf))
            breaks <- sort(breaks)
        } else { # All values were NA after filtering
            warning("All posterior_mean values are NA. Using a single bin.")
            breaks <- c(-Inf, Inf)
        }
    } else {
        warning("No valid data for posterior_mean to compute quantile breaks. Using a single bin.")
        breaks <- c(-Inf, Inf)
    }

    # Add bin column using the new dynamic breaks
    merged_gene_data <- merged_gene_data %>%
        mutate(
            bin = cut(posterior_mean, breaks = breaks, include.lowest = TRUE, right = TRUE)
        )

    mean_effect_fn <- function(beta_persd, burden_score) sum(beta_persd * sqrt(burden_score)) / sum(burden_score)
    mean_effect_se_fn <- function(beta_persd_se, burden_score) sqrt(sum(beta_persd_se^2 * burden_score)) / sum(burden_score)

    combined_bin_stats <- merged_gene_data %>%
        group_by(bin) %>%
        summarise(
            count = n(),             
            expected = mean(posterior_mean),
            expected_se = if_else(count > 0,
                                  sqrt(mean(posterior_variance) / count),
                                  NA_real_),
            mean_aou = if_else(count > 0, 
                                mean_effect_fn(effect_estimate.aou, burden_score.aou), 
                                NA_real_),
            mean_aou_se = if_else(count > 0,
                                  mean_effect_se_fn(effect_se.aou, burden_score.aou),
                                  NA_real_),
            mean_ukbb = if_else(count > 0, 
                                mean_effect_fn(effect_estimate, burden_score), 
                                NA_real_),
            mean_ukbb_se = if_else(count > 0,
                                  mean_effect_se_fn(effect_se, burden_score),
                                  NA_real_),
            .groups = 'drop'
        )

    print(combined_bin_stats)
    return(combined_bin_stats)
}

#' Compute expected and observed replication rates between UKBB and AoU studies
#' 
#' @param ukbb_model Fitted burdenEM model from UK Biobank data
#' @param aou_model Fitted burdenEM model from All of Us data
#' @param pval_threshold Two-tailed p-value threshold for filtering significant genes in UKBB (NULL for 0.05/num_genes)
#' @param replication_alpha One-tailed alpha level for replication (default: 0.01)
#' @return List with expected and observed replication rates, and number of genes analyzed
compute_pvalue_replication <- function(ukbb_model, aou_model, pval_threshold = NULL, replication_alpha = 0.01) {
    # Filter UKBB genes by p-value threshold if specified
    if (is.null(pval_threshold)) {
        pval_threshold <- 0.05 / nrow(ukbb_model$df)
    }
    ukbb_model$df <- ukbb_model$df %>%
        mutate(p_value = 2 * pnorm(-abs(effect_estimate / effect_se))) %>%
        filter(p_value < pval_threshold)
    
    # Calculate variance ratio (ratio of mean variances)
    var_ratio <- mean(aou_model$df$effect_se^2) / mean(ukbb_model$df$effect_se^2)
    
    # Define 1-tailed power function for replication
    power_function <- function(beta, row) {
        z_critical <- qnorm(replication_alpha, lower.tail = FALSE)
        ncp <- abs(beta) / (row$effect_se * sqrt(var_ratio))
        pnorm(z_critical, mean = ncp, lower.tail = FALSE)
    }
    
    # Early return if no significant genes
    if (nrow(ukbb_model$df) == 0) {
        return(list(
            expected_rate = NA_real_,
            observed_rate = NA_real_,
            n_genes = 0,
            merged_data = data.frame()
        ))
    }
    
    # Calculate expected replication rate for significant genes
    ukbb_model$df$expected_replication <- posterior_expectation2(
        ukbb_model,
        function_to_integrate = power_function
    )
    
    # Merge with AoU data and calculate observed replications
    merged_data <- inner_join(
        ukbb_model$df %>% select(gene_id, effect_estimate, expected_replication, effect_se),
        aou_model$df %>% select(gene, effect_estimate.aou = effect_estimate, effect_se.aou = effect_se),
        by = c("gene_id" = "gene")
    ) %>%
    mutate(
        z_score_aou = effect_estimate.aou / effect_se.aou,
        p_value_aou = pnorm(-abs(z_score_aou)),  # One-tailed p-value
        same_direction = sign(effect_estimate) == sign(effect_estimate.aou),
        replicated = same_direction & (p_value_aou < replication_alpha)
    )
    
    # Return results
    list(
        expected_rate = mean(merged_data$expected_replication, na.rm = TRUE),
        observed_rate = mean(merged_data$replicated, na.rm = TRUE),
        n_genes = nrow(merged_data),
        merged_data = merged_data
    )
}

perform_replication <- function(annotation, ukbb_pheno, aou_pheno, aou_ancestry, convert_to_per_allele = FALSE) {

    ukbb_model_file <- paste0("fitted_models/burdenEM_fit_genebass_", ukbb_pheno, "_", annotation, ".rds")
    if(file.exists(ukbb_model_file)){
        ukbb_model <- readRDS(ukbb_model_file)
        ukbb_model$df <- ukbb_model$df %>% dplyr::filter(!is.na(gene_id))
    } else {
        stop(paste("UKBB model file not found:", ukbb_model_file))
    }

    aou_model_file <- paste0("fitted_models/burdenEM_fit_aou_", aou_ancestry, "_", aou_pheno, "_", annotation, ".rds")
    if(file.exists(aou_model_file)){
        aou_model <- readRDS(aou_model_file)
        aou_model$df <- aou_model$df %>% dplyr::filter(!is.na(gene))
    } else {
        stop(paste("AoU gene data file not found:", aou_model_file))
    }

    mean_summary <- compute_posterior_mean_replication(ukbb_model, aou_model, convert_to_per_allele = convert_to_per_allele)
    pval_summary <- compute_pvalue_replication(ukbb_model, aou_model)

    return(list(mean_summary = mean_summary, pval_summary = pval_summary))

}

# --- Script Execution (only if run directly, not sourced) ---
if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 4) {
      stop("Usage: Rscript replication.R <annotation> <ukbb_pheno> <aou_pheno> <aou_ancestry>", call. = FALSE)
    }
    
    perform_replication(args[1], args[2], args[3], args[4])
}
