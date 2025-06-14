#!/usr/bin/env Rscript

# --- Libraries ---
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(future))
source("R/model.R") 

# Helper function for binomial confidence intervals using the Clopper-Pearson exact method
binomial_ci <- function(x, n, conf_level = 0.95) {
    alpha <- 1 - conf_level
    lower <- qbeta(alpha/2, x, n - x + 1)
    upper <- qbeta(1 - alpha/2, x + 1, n - x)
    data.frame(estimate = x/n, lower = lower, upper = upper)
}

compute_effect_replication_metrics_for_pair <- function(primary_model, replication_model, verbose = FALSE, convert_to_per_allele = TRUE) {

    if(convert_to_per_allele){
        effect_size_fn <- function(x, row) x / sqrt(row$burden_score)
        effect_sq_fn <- function(x, row) x^2 / row$burden_score
        mean_effect_fn <- function(beta_persd, burden_score) sum(beta_persd * sqrt(burden_score)) / sum(burden_score)
        mean_effect_se_fn <- function(beta_persd_se, burden_score) sqrt(sum(beta_persd_se^2 * burden_score)) / sum(burden_score)
    } else {
        effect_size_fn <- function(x) x
        effect_sq_fn <- function(x) x^2
        mean_effect_fn <- function(beta_persd, burden_score) mean(beta_persd)
        mean_effect_se_fn <- function(beta_persd_se, burden_score) sqrt(mean(beta_persd_se^2))
    }
    
    # --- Compute posterior-mean effect sizes for UKBB gene data ---
    primary_model$df$posterior_mean <- posterior_expectation2(primary_model, function_to_integrate = effect_size_fn)
    primary_model$df$posterior_mean_sq <- posterior_expectation2(primary_model, function_to_integrate = effect_sq_fn)
    primary_model$df$posterior_variance <- primary_model$df$posterior_mean_sq - (primary_model$df$posterior_mean^2)

    # --- 5. Merge AoU Gene Data with UKBB Gene Data ---
    if (!("gene" %in% names(primary_model$df))) stop("Column 'gene' not found in primary model gene data.")
    if (!("gene" %in% names(replication_model$df))) stop("Column 'gene' not found in replication model gene data.")
    
    primary_df_for_join <- primary_model$df %>%
      dplyr::select(gene, posterior_mean, posterior_variance, effect_estimate, burden_score, effect_se)

    replication_df_for_join <- replication_model$df %>%
      dplyr::select(gene, effect_estimate, burden_score, effect_se)
    
    merged_gene_data <- inner_join(primary_df_for_join, replication_df_for_join, by = c("gene" = "gene"), suffix = c(".primary", ".replication"))

    # --- 6. Adjust UKBB posterior means, bin, and summarize with AoU effects ---
    if (!("posterior_mean" %in% names(merged_gene_data))) {
        stop("Column 'posterior_mean' (from UKBB) not found in merged_gene_data.")
    }
    if (!("effect_estimate.replication" %in% names(merged_gene_data))) {
        stop("Column 'effect_estimate.replication' not found in merged_gene_data. Check merge suffixes.")
    }


    # --- Define breaks based on quantiles of posterior_mean ---
    if (sum(!is.na(merged_gene_data$posterior_mean)) > 0) {
        quantile_probs <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
        apm_values <- merged_gene_data$posterior_mean[!is.na(merged_gene_data$posterior_mean)]
        if (length(apm_values) > 0) {
            quantile_values <- quantile(apm_values, probs = quantile_probs, na.rm = FALSE, type = 7)
            breaks <- unique(c(-Inf, quantile_values, Inf))
            breaks <- sort(breaks)
        } else {
            if (verbose_param) message("  Warning: All posterior_mean values are NA after merging. Using a single bin.")
            breaks <- c(-Inf, Inf)
        }
    } else {
        if (verbose_param) message("  Warning: No valid data for posterior_mean to compute quantile breaks after merging. Using a single bin.")
        breaks <- c(-Inf, Inf)
    }

    merged_gene_data <- merged_gene_data %>%
        mutate(bin = cut(posterior_mean, breaks = breaks, include.lowest = TRUE, right = TRUE))

    combined_bin_stats <- merged_gene_data %>%
        group_by(bin) %>%
        summarise(
            count = n(),             
            expected_posterior_mean = mean(posterior_mean),
            expected_posterior_mean_se = if_else(count > 0,
                                  sqrt(mean(posterior_variance) / count),
                                  NA_real_),
            observed_effect_replication = if_else(count > 0, 
                                mean_effect_fn(effect_estimate.replication, burden_score.replication), 
                                NA_real_),
            observed_effect_replication_se = if_else(count > 0,
                                  mean_effect_se_fn(effect_se.replication, burden_score.replication),
                                  NA_real_),
            mean_effect_primary = if_else(count > 0, 
                                mean_effect_fn(effect_estimate.primary, burden_score.primary), 
                                NA_real_),
            mean_effect_primary_se = if_else(count > 0,
                                  mean_effect_se_fn(effect_se.primary, burden_score.primary),
                                  NA_real_),
            .groups = 'drop'
        )

    return(combined_bin_stats)
}

# --- Helper function to process a single primary-replication study pair ---
process_effect_replication_pair <- function(primary_study_row, replication_study_row, annotation_param, verbose_param) {
    if (verbose_param) {
        message(sprintf("Processing effect replication pair: Primary ('%s' - %s) vs Replication ('%s' - %s) for annotation '%s'", 
                        primary_study_row$identifier, primary_study_row$dataset, 
                        replication_study_row$identifier, replication_study_row$dataset, annotation_param))
    }

    primary_model_path <- stringr::str_replace(primary_study_row$model_filename, "<ANNOTATION>", annotation_param)
    replication_model_path <- stringr::str_replace(replication_study_row$model_filename, "<ANNOTATION>", annotation_param)

    primary_model_obj <- NULL
    replication_model_obj <- NULL
    error_message_val <- NA_character_

    tryCatch({
        primary_model_obj <- readRDS(primary_model_path)
    }, error = function(e) {
        error_message_val <<- sprintf("Error loading primary model '%s': %s", primary_model_path, e$message)
        if (verbose_param) message(paste("  ", error_message_val))
    })

    if (!is.null(error_message_val) && !is.na(error_message_val)) {
      # If primary failed, create an empty tibble with expected columns for this pair
      return(tibble(
          abbreviation = primary_study_row$abbreviation,
          primary_dataset = primary_study_row$dataset,
          replication_dataset = replication_study_row$dataset,
          annotation = annotation_param,
          bin = NA_character_,
          count = NA_integer_,
          expected_posterior_mean = NA_real_,
          expected_posterior_mean_se = NA_real_,
          observed_effect_replication = NA_real_,
          observed_effect_replication_se = NA_real_,
          mean_effect_primary = NA_real_,
          mean_effect_primary_se = NA_real_
      ))
    }

    tryCatch({
        replication_model_obj <- readRDS(replication_model_path)
    }, error = function(e) {
        error_message_val <<- sprintf("Error loading replication model '%s': %s", replication_model_path, e$message)
        if (verbose_param) message(paste("  ", error_message_val))
    })
    
    # Prepare results tibble shell, to be filled or returned with error
    results_shell <- tibble(
        abbreviation = primary_study_row$abbreviation,
        primary_dataset = primary_study_row$dataset,
        replication_dataset = replication_study_row$dataset,
        annotation = annotation_param,
        bin = NA_character_,
        count = NA_integer_,
        expected_posterior_mean = NA_real_,
        expected_posterior_mean_se = NA_real_,
        observed_effect_replication = NA_real_,
        observed_effect_replication_se = NA_real_,
        mean_effect_primary = NA_real_,
        mean_effect_primary_se = NA_real_
    )

    if (!is.null(primary_model_obj) && !is.null(replication_model_obj)) {
        pair_metrics <- compute_effect_replication_metrics_for_pair(primary_model_obj, replication_model_obj, verbose=verbose_param)
        
        if (nrow(pair_metrics) > 0) {
            results <- pair_metrics %>%
                mutate(
                    abbreviation = primary_study_row$abbreviation,
                    primary_dataset = primary_study_row$dataset,
                    replication_dataset = replication_study_row$dataset,
                    annotation = annotation_param
                ) %>%
                select(abbreviation, primary_dataset, replication_dataset, annotation, 
                       bin, count, expected_posterior_mean, expected_posterior_mean_se, 
                       observed_effect_replication, observed_effect_replication_se, 
                       mean_effect_primary, mean_effect_primary_se)
            return(results)
        } else {
            if(verbose_param && is.na(error_message_val)) { # if no prior loading error was already printed
                 message("  Info: Computation for pair resulted in no data (e.g. no common genes, missing columns, or all NA bins).")
            }
            return(results_shell)
        }
    } else {
        return(results_shell)
    }
}

# --- Main Exported Function ---
calculate_effect_replication_metrics <- function(studies_df, annotation_param, primary_dataset_identifier, 
                                               verbose_param = FALSE, run_sequentially_param = FALSE) {
    if (!is.data.frame(studies_df) || nrow(studies_df) == 0) {
        if (verbose_param) message("Input studies_df is empty or not a dataframe. Returning empty results.")
        # Define expected columns for an empty result, matching process_effect_replication_pair output
        return(tibble(
            primary_abbreviation = character(), primary_dataset = character(), 
            replication_abbreviation = character(), replication_dataset = character(), 
            annotation = character(), bin = character(), count = integer(), 
            expected_posterior_mean = double(), expected_posterior_mean_se = double(), 
            observed_effect_replication = double(), observed_effect_replication_se = double(), 
            mean_effect_primary = double(), mean_effect_primary_se = double(), 
            error_message = character()
        ))
    }

    # Filter for primary studies
    primary_studies <- studies_df %>%
        filter(dataset == primary_dataset_identifier)

    if (nrow(primary_studies) == 0) {
        if (verbose_param) message(sprintf("No primary studies found for dataset_identifier '%s'.", primary_dataset_identifier))
        return(tibble()) # Match empty structure above
    }

    pair_list <- list()
    unique_identifiers <- unique(studies_df$identifier) # Use all unique traits from original studies_df

    for (id in unique_identifiers) {
        current_trait_primary_studies <- primary_studies %>%
            filter(identifier == id)

        if (nrow(current_trait_primary_studies) == 1) {
            primary_study_row <- current_trait_primary_studies[1, , drop = FALSE]
            
            replication_study_candidates <- studies_df %>%
                filter(identifier == id, dataset != primary_dataset_identifier)

            if (nrow(replication_study_candidates) > 0) {
                for (j in 1:nrow(replication_study_candidates)) {
                    replication_study_row <- replication_study_candidates[j, , drop = FALSE]
                    pair_list[[length(pair_list) + 1]] <- list(
                        primary = primary_study_row, 
                        replication = replication_study_row
                    )
                }
            } else {
                if (verbose_param) message(sprintf("Info: Trait '%s' has a primary study in '%s' but no replication studies. Skipping for effect replication.", id, primary_dataset_identifier))
            }
        } else if (nrow(current_trait_primary_studies) > 1) {
            if (verbose_param) message(sprintf("Warning: Trait '%s' has %d studies matching primary_dataset_identifier '%s'. Skipping this trait for effect replication.", 
                                               id, nrow(current_trait_primary_studies), primary_dataset_identifier))
        } else {
            # This case should not be reached if primary_studies is filtered correctly and id comes from studies_df
            # but as a safeguard:
            if (verbose_param) message(sprintf("Info: Trait '%s' has no study matching primary_dataset_identifier '%s' (unexpected). Skipping.", 
                                               id, primary_dataset_identifier))
        }
    }

    if (length(pair_list) == 0) {
        if (verbose_param) message("No valid primary-replication pairs found to process for effect replication.")
        return(tibble()) # Match empty structure above
    }

    if (!run_sequentially_param && future::availableCores() > 1) {
        if (verbose_param) {
            message(sprintf("Setting up parallel execution for effect replication calculation using up to %d workers for %d pairs.", 
                            future::availableCores(), length(pair_list)))
        }
        future::plan(future::multisession, workers = min(future::availableCores(), length(pair_list), 16)) # Cap workers
        all_results_df <- furrr::future_map_dfr(
            .x = pair_list,
            .f = ~process_effect_replication_pair(.x$primary, .x$replication, annotation_param, verbose_param),
            .progress = verbose_param,
            .options = furrr_options()
        )
    } else {
        if (verbose_param) {
            message(sprintf("Running effect replication calculation in sequential mode for %d pairs.", length(pair_list)))
        }
        all_results_df <- purrr::map_dfr(
            .x = pair_list,
            .f = ~process_effect_replication_pair(.x$primary, .x$replication, annotation_param, verbose_param)
        )
    }
    
    return(all_results_df)
}
