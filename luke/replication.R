#!/usr/bin/env Rscript

# --- Libraries ---
# Install future and furrr if not already installed
if (!requireNamespace("future", quietly = TRUE)) {
  message("Installing 'future' package...")
  install.packages("future", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("furrr", quietly = TRUE)) {
  message("Installing 'furrr' package...")
  install.packages("furrr", repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
source("R/model.R")

# --- Utility Functions ---
# Helper function for binomial confidence intervals using the Clopper-Pearson exact method
binomial_ci <- function(x, n, conf_level = 0.95) {
    alpha <- 1 - conf_level
    lower <- qbeta(alpha/2, x, n - x + 1)
    upper <- qbeta(1 - alpha/2, x + 1, n - x)
    data.frame(estimate = x/n, lower = lower, upper = upper)
}

#' Compute expected and observed replication rates between a primary and a replication model
#'
#' @param primary_model_obj Fitted burdenEM model from the primary study
#' @param replication_model_obj Fitted burdenEM model from the replication study
#' @param pval_threshold Two-tailed p-value threshold for filtering significant genes in the primary model (NULL for 0.05/num_genes)
#' @param replication_alpha One-tailed alpha level for replication (default: 0.01)
#' @return List with expected and observed replication rates, and number of genes analyzed
compute_replication_metrics_for_pair <- function(primary_model_obj, replication_model_obj, pval_threshold = NULL, replication_alpha = 0.01, verbose_param = FALSE) {
    # Ensure df exists and is not empty in both models

    if (is.null(primary_model_obj$df) || nrow(primary_model_obj$df) == 0 || 
        is.null(replication_model_obj$df) || nrow(replication_model_obj$df) == 0) {
        if (verbose_param) message("  One or both models lack a 'df' component or 'df' is empty.")
        return(list(expected_rate = NA_real_, observed_rate = NA_real_, n_genes_tested_replication = 0L, n_significant_in_primary = 0L, merged_data = data.frame()))
    }

    # Ensure $df components are data.frames for reliable dplyr operation
    primary_model_obj$df <- as.data.frame(primary_model_obj$df)
    replication_model_obj$df <- as.data.frame(replication_model_obj$df)

    # Filter primary model genes by p-value threshold if specified
    if (is.null(pval_threshold)) {
        pval_threshold <- 0.05 / nrow(primary_model_obj$df)
    }
    # Calculate p-values for primary model
    primary_significant_genes_df <- primary_model_obj$df %>%
        rowwise() %>%
        mutate(
            one_tailed_p_val = primary_model_obj$pval_function(cur_data_all()),
        ) %>%
        ungroup() %>%
        filter(2 * pmin(one_tailed_p_val, 1 - one_tailed_p_val) < pval_threshold)

    if (verbose_param) message("  Found ", nrow(primary_significant_genes_df), " significant genes in primary model.")
    # Early return if no significant genes in primary study
    if (nrow(primary_significant_genes_df) == 0) {
        return(list(expected_rate = NA_real_, observed_rate = NA_real_, n_genes_tested_replication = 0L, n_significant_in_primary = 0L, merged_data = data.frame()))
    }


    # Calculate samplesize ratio for power calculations
    primary_rel_samplesize <- primary_model_obj$df %>%
        rowwise() %>%
        mutate(samplesize_value = primary_model_obj$rel_samplesize_function(cur_data_all())) %>%
        ungroup() %>%
        summarise(mean_samplesize = mean(samplesize_value)) %>%
        pull(mean_samplesize)
    
    replication_rel_samplesize <- replication_model_obj$df %>%
        rowwise() %>%
        mutate(samplesize_value = replication_model_obj$rel_samplesize_function(cur_data_all())) %>%
        ungroup() %>%
        summarise(mean_samplesize = mean(samplesize_value)) %>%
        pull(mean_samplesize)
    
    samplesize_ratio <- replication_rel_samplesize / primary_rel_samplesize    
    power_function <- replication_model_obj$get_power_function(replication_alpha, samplesize_ratio)

    if (verbose_param) message(sprintf("  Sample size ratio between replication and primary: %.3f", samplesize_ratio))

    # Calculate expected replication rate for significant genes
    temp_primary_model_for_posterior_exp <- primary_model_obj # Copy structure
    temp_primary_model_for_posterior_exp$df <- primary_significant_genes_df 
    primary_significant_genes_df <- primary_significant_genes_df %>% 
        mutate(expected_replication = posterior_expectation2(
            model = temp_primary_model_for_posterior_exp,
            function_to_integrate = power_function
        ))

    # Merge with replication data and calculate observed replications
    merged_data <- inner_join(
        primary_significant_genes_df %>% select(gene, expected_replication, one_tailed_p_val),
        replication_model_obj$df %>%
        rowwise() %>%
        mutate(replication_one_tailed_p_val = replication_model_obj$pval_function(cur_data_all())) %>%
        ungroup() %>%
        select(gene, replication_one_tailed_p_val),
        by = "gene"
    ) %>%
    mutate(
        replicated = ifelse(
            one_tailed_p_val < 0.5,
            replication_one_tailed_p_val < replication_alpha,
            replication_one_tailed_p_val > 1 - replication_alpha
        )
    )

    if (verbose_param) message("  Replicated ", sum(merged_data$replicated), " out of ", nrow(primary_significant_genes_df), " significant genes in primary model.")
    return(
        list(
            expected_rate = mean(merged_data$expected_replication, na.rm = TRUE),
            observed_rate = mean(merged_data$replicated, na.rm = TRUE),
            n_genes_tested_replication = nrow(merged_data),
            n_significant_in_primary = nrow(primary_significant_genes_df)
        )
    )
}

# --- Helper function to process a single primary-replication study pair ---
process_replication_pair <- function(primary_study_row, replication_study_row, annotation_param, verbose_param, pval_threshold_param = NULL) {
    if (verbose_param) {
        message(sprintf("Processing replication pair: Primary ('%s' - %s) vs Replication ('%s' - %s) for annotation '%s'", 
                        primary_study_row$identifier, primary_study_row$dataset, 
                        replication_study_row$identifier, replication_study_row$dataset, annotation_param))
    }

    primary_model_path <- stringr::str_replace(primary_study_row$model_filename, "<ANNOTATION>", annotation_param)
    replication_model_path <- stringr::str_replace(replication_study_row$model_filename, "<ANNOTATION>", annotation_param)

    # Initialize results list
    results <- list(
        primary_abbreviation = primary_study_row$abbreviation,
        primary_dataset = primary_study_row$dataset,
        replication_abbreviation = replication_study_row$abbreviation,
        replication_dataset = replication_study_row$dataset,
        annotation = annotation_param,
        expected_replication_rate = NA_real_,
        observed_replication_rate = NA_real_,
        n_genes_tested_replication = NA_integer_,
        n_significant_in_primary = NA_integer_,
        error_message = NA_character_
    )

    primary_model_obj <- tryCatch({
        model <- readRDS(primary_model_path)
        model
    }, error = function(e) {
        if (verbose_param) message(sprintf("  Error loading primary model '%s': %s", primary_model_path, e$message))
        return(NULL)
    })

    replication_model_obj <- tryCatch({
        model <- readRDS(replication_model_path)
        model
    }, error = function(e) {
        if (verbose_param) message(sprintf("  Error loading replication model '%s': %s", replication_model_path, e$message))
        return(NULL)
    })

    # Check if models loaded successfully before proceeding
    if (is.null(primary_model_obj)) {
        results$error_message <- sprintf("Failed to load primary model: %s", primary_model_path)
    } else if (is.null(replication_model_obj)) {
        results$error_message <- sprintf("Failed to load replication model: %s (Primary: %s OK)", replication_model_path, primary_model_path)
    }

    if (!is.null(results$error_message) && !is.na(results$error_message)) {
        # Return early if there was a loading error
        if (verbose_param) { 
            message(sprintf("  Model loading failed. Error: %s", results$error_message))
        }
        return(as.data.frame(results))
    }

    # Compute replication metrics, wrapped in tryCatch
    tryCatch({
        replication_calc_results <- compute_replication_metrics_for_pair(
            primary_model_obj,
            replication_model_obj,
            pval_threshold = pval_threshold_param,
            verbose_param = verbose_param
        )
        results$expected_replication_rate <- replication_calc_results$expected_rate
        results$observed_replication_rate <- replication_calc_results$observed_rate
        results$n_genes_tested_replication <- replication_calc_results$n_genes_tested_replication
        results$n_significant_in_primary <- replication_calc_results$n_significant_in_primary
        # If successful, error_message remains NA (its initial value)
    }, error = function(e) {
        error_msg_text <- paste("Error during compute_replication_metrics_for_pair:", conditionMessage(e))
        results$error_message <<- error_msg_text # Use <<- to assign to 'results' in the parent environment
        if (verbose_param) {
            message(error_msg_text)
        }
    })

    return(as.data.frame(results))
}

# --- Main function to calculate replication metrics across study pairs ---
calculate_replication_metrics <- function(studies_df, annotation_param, primary_dataset_identifier, 
                                        pval_threshold_param = NULL, verbose_param = FALSE, run_sequentially_param = FALSE) {
    if (!is.data.frame(studies_df) || nrow(studies_df) == 0) {
        if (verbose_param) message("Input studies_df is empty or not a dataframe. Returning empty results.")
        # Define expected columns for an empty result
        expected_cols <- c("primary_abbreviation", "primary_dataset", 
                           "replication_abbreviation", "replication_dataset", 
                           "annotation", "expected_replication_rate", "observed_replication_rate", 
                           "n_genes_tested_replication", "n_significant_in_primary", "error_message")
        return(stats::setNames(data.frame(matrix(ncol = length(expected_cols), nrow = 0)), expected_cols))
    }


    pair_list <- list()
    unique_identifiers <- unique(studies_df$identifier)

    for (id in unique_identifiers) {
        # Using base R subsetting to avoid dplyr::filter weirdness with model_filename
        trait_studies <- studies_df[studies_df$identifier == id, , drop = FALSE]
        primary_study_candidates <- trait_studies %>% filter(dataset == primary_dataset_identifier)
        
        if (nrow(primary_study_candidates) == 1) {
            primary_study_row <- primary_study_candidates[1, , drop = FALSE]
            replication_studies_for_trait <- trait_studies %>% filter(dataset != primary_dataset_identifier)
            
            if (nrow(replication_studies_for_trait) > 0) {
                for (j in 1:nrow(replication_studies_for_trait)) {
                    replication_study_row <- replication_studies_for_trait[j, , drop = FALSE]
                    pair_list[[length(pair_list) + 1]] <- list(
                        primary = primary_study_row, 
                        replication = replication_study_row
                    )
                }
            } else {
                if (verbose_param) message(sprintf("Info: Trait '%s' has a primary study in '%s' but no replication studies. Skipping.", id, primary_dataset_identifier))
            }
        } else if (nrow(primary_study_candidates) > 1) {
            if (verbose_param) message(sprintf("Warning: Trait '%s' has %d studies matching primary_dataset_identifier '%s'. Skipping this trait for replication.", 
                                               id, nrow(primary_study_candidates), primary_dataset_identifier))
        } else {
            if (verbose_param) message(sprintf("Info: Trait '%s' has no study matching primary_dataset_identifier '%s'. Skipping this trait for replication.", 
                                               id, primary_dataset_identifier))
        }
    }

    if (length(pair_list) == 0) {
        if (verbose_param) message("No valid primary-replication pairs found to process.")
        expected_cols <- c("primary_identifier", "primary_dataset", "primary_trait_type", "primary_model_path",
                           "replication_identifier", "replication_dataset", "replication_trait_type", "replication_model_path",
                           "annotation", "expected_replication_rate", "observed_replication_rate", 
                           "n_genes_tested_replication", "n_significant_in_primary", "error_message")
        return(stats::setNames(data.frame(matrix(ncol = length(expected_cols), nrow = 0)), expected_cols))
    }

    if (!run_sequentially_param) {
        if (verbose_param) {
            message(sprintf("Setting up parallel execution for replication calculation using up to %d workers for %d pairs.", 
                            future::availableCores(), length(pair_list)))
        }
        future::plan(future::multisession)
        all_results_df <- furrr::future_map_dfr(
            .x = pair_list,
            .f = ~process_replication_pair(.x$primary, .x$replication, annotation_param, verbose_param, pval_threshold_param = pval_threshold_param),
            .progress = verbose_param
        )
    } else {
        if (verbose_param) {
            message(sprintf("Running replication calculation in sequential mode for %d pairs.", length(pair_list)))
        }
        all_results_df <- purrr::map_dfr(
            .x = pair_list,
            .f = ~process_replication_pair(.x$primary, .x$replication, annotation_param, verbose_param, pval_threshold_param = pval_threshold_param)
        )
    }
    
    return(all_results_df)
}

# Old perform_replication and script execution block are removed as they are superseded by CLI structure.
