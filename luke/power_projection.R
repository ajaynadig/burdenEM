#!/usr/bin/env Rscript

if (!requireNamespace("future", quietly = TRUE)) {
  install.packages("future", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("furrr", quietly = TRUE)) {
  install.packages("furrr", repos = "http://cran.us.r-project.org")
}
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(future)
library(furrr)

source("R/power.R")

# Sample size multipliers: 0.5x to ~1Mx current size
SAMPLE_SIZE_FACTORS <- 2^(-1:20)

is_binary_power_projection_model <- function(model) {
  if (is.null(model) || is.null(model$df)) {
    return(FALSE)
  }

  df_names <- names(model$df)
  ("CAC_cases" %in% df_names) && ("expected_CAC_cases" %in% df_names)
}

process_single_study_for_power_projection <- function(study_row, current_annotation, current_verbose) {
  trait_name_from_df <- if ("identifier" %in% names(study_row)) study_row$identifier else study_row$description
  dataset_from_df <- study_row$dataset
  model_filename_pattern <- study_row$model_filename

  if (current_verbose) {
    message(sprintf("Processing Power Projection for: %s (Dataset: %s, Annotation: %s)",
                    trait_name_from_df, dataset_from_df, current_annotation))
  }

  if (length(model_filename_pattern) != 1 || !is.character(model_filename_pattern) || is.na(model_filename_pattern)) {
    if (current_verbose) {
      message(sprintf("  Skipping %s (%s): model_filename is invalid, NA, or missing.", trait_name_from_df, dataset_from_df))
    }
    return(NULL)
  }

  model_file <- stringr::str_replace(model_filename_pattern, "<ANNOTATION>", current_annotation)

  if (!file.exists(model_file)) {
    if (current_verbose) {
      message(sprintf("  Model file not found: %s", model_file))
    }
    return(NULL)
  }

  model <- tryCatch({
    readRDS(model_file)
  }, error = function(e) {
    if (current_verbose) {
      message(sprintf("  Error loading %s for %s (%s): %s",
                      model_file, trait_name_from_df, dataset_from_df, e$message))
    }
    return(NULL)
  })

  if (is.null(model)) {
    return(NULL)
  }

  if (is_binary_power_projection_model(model)) {
    if (current_verbose) {
      message(sprintf("  Skipping binary trait model: %s", model_file))
    }
    return(NULL)
  }

  tryCatch({
    # Extract component endpoints and compute aggregated mixture weights
    components <- model$component_endpoints
    features <- do.call(rbind, model$df$features)
    mixing_weights <- features %*% model$delta  # n_genes x n_components
    aggregated_weights <- colMeans(mixing_weights)
    aggregated_weights[aggregated_weights < 0] <- 0
    aggregated_weights <- aggregated_weights / sum(aggregated_weights)

    # Get number of genes and effective sample size
    n_genes <- nrow(features)
    rel_samplesizes <- sapply(1:n_genes, function(i) model$rel_samplesize_function(model$df[i, ]))
    mean_nn <- mean(rel_samplesizes)

    # Exome-wide significance threshold
    alpha <- 0.05 / n_genes

    # Run NTPR across sample size factors
    result <- NTPR(
      bb = components,
      pp = aggregated_weights,
      nn = mean_nn * SAMPLE_SIZE_FACTORS,
      alpha = alpha
    )

    # Build output data frame (one row per sample size factor)
    data.frame(
      abbreviation = study_row$abbreviation,
      dataset = dataset_from_df,
      sample_size_factor = SAMPLE_SIZE_FACTORS,
      effective_sample_size = mean_nn * SAMPLE_SIZE_FACTORS,
      total_power = as.vector(result$total_power),
      total_ntpr = as.vector(result$total_ntpr),
      total_h2gwas = as.vector(result$total_h2gwas),
      total_power_pos = as.vector(result$total_power_pos),
      total_power_neg = as.vector(result$total_power_neg),
      n_genes = n_genes,
      mean_effective_n = mean_nn,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    if (current_verbose) {
      message(sprintf("  Error processing %s for %s (%s): %s",
                      model_file, trait_name_from_df, dataset_from_df, e$message))
    }
    return(NULL)
  })
}

calculate_power_projection_for_studies <- function(studies_df, annotation_param, verbose_param = FALSE, run_sequentially_param = FALSE) {

  if (!is.data.frame(studies_df) || nrow(studies_df) == 0) {
    if (verbose_param) {
      message("Input studies_df is empty or not a dataframe. Returning empty results.")
    }
    return(tibble::tibble())
  }

  trait_col_name <- if ("identifier" %in% names(studies_df)) "identifier" else "description"
  required_cols <- c(trait_col_name, "dataset", "model_filename")
  missing_cols <- setdiff(required_cols, names(studies_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in studies_df:", paste(missing_cols, collapse = ", ")))
  }

  if (length(annotation_param) != 1 || !is.character(annotation_param) || is.na(annotation_param)) {
    stop("annotation_param must be a single, non-NA string.")
  }

  if (!run_sequentially_param) {
    if (verbose_param) {
      message(sprintf("Setting up parallel execution for power projection using up to %d workers.", future::availableCores()))
    }
    future::plan(future::multisession)

    results_list <- furrr::future_map(
      .x = 1:nrow(studies_df),
      .f = ~process_single_study_for_power_projection(studies_df[.x, ], annotation_param, verbose_param),
      .progress = verbose_param
    )
  } else {
    if (verbose_param) {
      message("Running power projection calculation in sequential mode.")
    }
    results_list <- purrr::map(
      .x = 1:nrow(studies_df),
      .f = ~process_single_study_for_power_projection(studies_df[.x, ], annotation_param, verbose_param)
    )
  }

  valid_results_list <- purrr::compact(results_list)

  if (length(valid_results_list) == 0) {
    if (verbose_param) {
      message("No studies could be processed successfully. Returning empty results.")
    }
    return(tibble::tibble())
  }

  results_df <- dplyr::bind_rows(valid_results_list)
  return(results_df)
}
