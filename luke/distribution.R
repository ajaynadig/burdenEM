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
library(future)
library(furrr)
library(dplyr)
library(readr) 
suppressPackageStartupMessages(library(stringr))
library(tibble) # For as_tibble, if needed, or for general tibble operations

# --- Source custom distribution functions ---
source("luke/distribution_functions.R")
source("R/estimate_heritability.R") # Source heritability functions for h2 context

# --- Helper function to process a single study for distribution metrics ---
process_single_study_for_distribution <- function(study_row, current_annotation, current_verbose, include_hp_for_1_gene = TRUE) {
  trait_name_from_df <- study_row$identifier # Assuming 'identifier' is the preferred trait name column
  dataset_from_df <- study_row$dataset
  model_filename_pattern <- study_row$model_filename
  actual_model_filename <- stringr::str_replace(model_filename_pattern, "<ANNOTATION>", current_annotation)

  if (current_verbose) {
    message(sprintf("Processing distribution for study: %s, dataset: %s, annotation: %s", 
                    trait_name_from_df, dataset_from_df, current_annotation))
    message(sprintf("  Attempting to load model from: %s", actual_model_filename))
  }

  if (!file.exists(actual_model_filename)) {
    if (current_verbose) message(sprintf("  Model file not found: %s. Skipping.", actual_model_filename))
    return(NULL)
  }
  
  model <- tryCatch({
    readRDS(actual_model_filename)
  }, error = function(e) {
    if (current_verbose) message(sprintf("  Error loading model '%s': %s", actual_model_filename, e$message))
    return(NULL)
  })

  if (is.null(model)) return(NULL)

  if (current_verbose) {
    message(sprintf("  Model loaded for %s. Structure check: model$df exists? %s, model$df$features exists? %s", 
                    trait_name_from_df, !is.null(model$df), !is.null(model$df$features)))
    if (!is.null(model$df) && !is.null(model$df$features)) {
        message(sprintf("  Number of features (genes) in model: %d", length(model$df$features)))
    } else {
        message("  Model df or features is NULL.")
    }
  }

  # --- TryCatch block for all metric calculations for this study ---
  processed_data <- tryCatch({
    if (current_verbose) message(sprintf("  Calculating distribution metrics for %s", trait_name_from_df))

    gene_qf <- get_gene_qf_betasq(model, right_tail=TRUE)
    variance_cdf <- get_variance_cdf_betasq(model, right_tail=TRUE)

    needed_genes_calculator <- get_needed_genes_fn(model)
    test_heritability_proportions <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999)
    
    plot_data_points_list <- list()
    num_total_genes <- length(model$df$features) 

    proportions_to_plot <- test_heritability_proportions
    if (include_hp_for_1_gene) {
      hp_for_1_gene <- NA
      if (num_total_genes > 0) {
        beta_sq_for_top_gene_fraction <- gene_qf(1/num_total_genes) 
        if (!is.na(beta_sq_for_top_gene_fraction)){
            hp_for_1_gene <- variance_cdf(beta_sq_for_top_gene_fraction)
        }
      }

      if (!is.na(hp_for_1_gene)) {
        proportions_to_plot <- sort(unique(c(hp_for_1_gene, proportions_to_plot)))
      }
    }

    for (hp_g in proportions_to_plot) {
      if(is.na(hp_g)) next
      fraction_genes_result <- needed_genes_calculator(hp_g)
      if(is.na(fraction_genes_result)) next

      current_num_genes = fraction_genes_result * num_total_genes
      
      if (current_num_genes >= (1 - 1e-9)) {
        plot_data_points_list[[length(plot_data_points_list) + 1]] <- 
          data.frame(H2_Percent = hp_g * 100, 
                     Num_Genes = current_num_genes)
      }
    }

    if (length(plot_data_points_list) > 0) {
      plot_df_h2_vs_genes <- dplyr::bind_rows(plot_data_points_list)
      plot_df_h2_vs_genes$Num_Genes <- pmax(1, plot_df_h2_vs_genes$Num_Genes)
      
      study_info_df <- as.data.frame(study_row)[, c("abbreviation", "dataset"), drop = FALSE]
      study_info_replicated_df <- study_info_df[rep(1, nrow(plot_df_h2_vs_genes)), , drop = FALSE]
      final_study_data_df <- dplyr::bind_cols(study_info_replicated_df, plot_df_h2_vs_genes)

      final_study_data_df$Num_Genes_SE <- NA_real_
      final_study_data_df$xmin_genes <- NA_real_
      final_study_data_df$xmax_genes <- NA_real_

      # SE calculation block remains commented out
#      if (!is.null(model$information) && !is.null(model$null_index) && nrow(final_study_data_df) > 0) {
#        # ... SE code ... 
#      }
      if (current_verbose) message(sprintf("  Finished processing distribution for %s. Found %d data points.", trait_name_from_df, nrow(final_study_data_df)))
      final_study_data_df # Return the successfully processed data
    } else {
      if (current_verbose) message(sprintf("  No data points to process for distribution for %s.", trait_name_from_df))
      empty_study_df <- as.data.frame(study_row)[, c("abbreviation", "dataset"), drop = FALSE]
      empty_study_df$H2_Percent <- NA_real_
      empty_study_df$Num_Genes <- NA_real_
      empty_study_df$Num_Genes_SE <- NA_real_
      empty_study_df$xmin_genes <- NA_real_
      empty_study_df$xmax_genes <- NA_real_
      empty_study_df[0, , drop=FALSE] # Return empty df with correct columns
    }
  }, error = function(e) {
    if (current_verbose) message(sprintf("  ERROR during metric calculation for study %s (model: %s): %s", 
                                      trait_name_from_df, actual_model_filename, e$message))
    return(NULL) # Return NULL if any error occurs during metric calculation for this study
  })
  
  return(processed_data)
}

# --- Main function to be called by CLI for calculating distribution metrics ---
calculate_distribution_metrics_for_studies <- function(studies_df, annotation, verbose = FALSE, run_sequentially_param = FALSE, include_hp_for_1_gene = TRUE) {
  if (!is.data.frame(studies_df) || nrow(studies_df) == 0) {
    if (verbose) message("Input studies_df is empty or not a dataframe. Returning empty results for distribution.")
    # Define expected columns for an empty result
    # Base columns from studies_df + specific metric columns
    expected_cols_df <- studies_df[0, c("abbreviation", "dataset"), drop = FALSE]
    metric_cols <- c("H2_Percent", "Num_Genes", "Num_Genes_SE", "xmin_genes", "xmax_genes")
    for(mc in metric_cols) expected_cols_df[[mc]] <- numeric(0)
    return(expected_cols_df)
  }

  if (!run_sequentially_param) {
    if (verbose) {
      message(sprintf("Setting up parallel execution for distribution calculation using up to %d workers.", future::availableCores()))
    }
    future::plan(future::multisession)
    all_results_df <- furrr::future_map_dfr(
      .x = 1:nrow(studies_df),
      .f = ~process_single_study_for_distribution(studies_df[.x, , drop = FALSE], annotation, verbose, include_hp_for_1_gene),
      .progress = verbose,
      .options = furrr::furrr_options(seed = TRUE) # For reproducibility if any RNG is used internally
    )
  } else {
    if (verbose) {
      message("Running distribution calculation in sequential mode.")
    }
    all_results_df <- purrr::map_dfr(
      .x = 1:nrow(studies_df),
      .f = ~process_single_study_for_distribution(studies_df[.x, , drop = FALSE], annotation, verbose, include_hp_for_1_gene)
    )
  }
  
  if (verbose) message(sprintf("Distribution metrics calculation complete. Total rows: %d", nrow(all_results_df)))
  return(all_results_df)
}
