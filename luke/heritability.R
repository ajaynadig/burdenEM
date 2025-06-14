#!/usr/bin/env Rscript

# Load necessary libraries
# Install future and furrr if not already installed
if (!requireNamespace("future", quietly = TRUE)) {
  message("Installing 'future' package...")
  install.packages("future", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("furrr", quietly = TRUE)) {
  message("Installing 'furrr' package...")
  install.packages("furrr", repos = "http://cran.us.r-project.org")
}
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(future)
library(furrr)
# readr could be removed if studies_df is always pre-loaded
# library(readr) 

source("R/estimate_heritability.R")

# Helper function to process a single study row for heritability estimation
process_single_study_for_heritability <- function(study_row, current_annotation, current_verbose) {
  trait_name_from_df <- if ("identifier" %in% names(study_row)) study_row$identifier else study_row$description
  dataset_from_df <- study_row$dataset
  model_filename_pattern <- study_row$model_filename

  if (current_verbose) {
    message(sprintf("Processing Heritability for: %s (Dataset: %s, Annotation: %s)", 
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

  tryCatch({
    model <- readRDS(model_file)
    h2_output <- estimate_heritability_components(model, verbose = current_verbose)
    
    result <- data.frame(
      abbreviation = study_row$abbreviation,
      dataset = dataset_from_df,
      total_h2 = h2_output$total_h2$mean,
      total_h2_se = h2_output$total_h2$se,
      positive_h2 = h2_output$positive_h2$mean,
      positive_h2_se = h2_output$positive_h2$se,
      negative_h2 = h2_output$negative_h2$mean,
      negative_h2_se = h2_output$negative_h2$se,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(h2_output$annot_h2) && !is.null(h2_output$annot_h2$mean) && length(h2_output$annot_h2$mean) > 0) {
      feature_names <- names(h2_output$annot_h2$mean)
      if (is.null(feature_names) && length(h2_output$annot_h2$mean) > 0) {
        feature_names <- paste0("feature", seq_along(h2_output$annot_h2$mean))
      }
      
      for (i in seq_along(feature_names)) {
        feat <- feature_names[i]
        # Ensure column names are syntactically valid
        col_name_h2 <- make.names(paste0(feat, "_h2"))
        col_name_se <- make.names(paste0(feat, "_h2_se"))
        result[[col_name_h2]] <- h2_output$annot_h2$mean[i]
        result[[col_name_se]] <- h2_output$annot_h2$se[i]
      }
    }
    
    return(result)
  }, error = function(e) {
    if (current_verbose) {
      message(sprintf("  Error processing %s for %s (%s): %s", 
                      model_file, trait_name_from_df, dataset_from_df, e$message))
    }
    return(NULL)
  })
}

#' Calculate heritability for a list of studies
#'
#' @param studies_df A dataframe where each row represents a study. Expected columns:
#'   `identifier` (or `description`) for trait name, `dataset`, and `model_filename` 
#'   (a pattern string with <ANNOTATION> placeholder).
#' @param annotation_param The specific annotation to use (e.g., "pLoF", "synonymous").
#' @param verbose_param Boolean, if TRUE, print progress messages. Defaults to FALSE.
#' @return A dataframe consolidating heritability results for all processed studies.
#'   Each row corresponds to a study, with columns for various h2 estimates.
#'   Returns an empty tibble if no studies can be processed or if input is invalid.
calculate_heritability_for_studies <- function(studies_df, annotation_param, verbose_param = FALSE, run_sequentially_param = FALSE) {
  
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
  
  if (length(annotation_param) != 1 || !is.character(annotation_param) || is.na(annotation_param)){
      stop("annotation_param must be a single, non-NA string.")
  }

  if (!run_sequentially_param) {
    if (verbose_param) {
      message(sprintf("Setting up parallel execution for heritability calculation using up to %d workers.", future::availableCores()))
    }
    future::plan(future::multisession) 
    
    results_list <- furrr::future_map(
      .x = 1:nrow(studies_df), 
      .f = ~process_single_study_for_heritability(studies_df[.x, ], annotation_param, verbose_param),
      .progress = verbose_param 
    )
  } else {
    if (verbose_param) {
      message("Running heritability calculation in sequential mode.")
    }
    results_list <- purrr::map(
      .x = 1:nrow(studies_df), 
      .f = ~process_single_study_for_heritability(studies_df[.x, ], annotation_param, verbose_param)
    )
  }
  
  valid_results_list <- purrr::compact(results_list)
  
  if (length(valid_results_list) == 0) {
    if (verbose_param) {
      message("No studies could be processed successfully or yielded valid results. Returning empty results.")
    }
    return(tibble::tibble())
  }
  
  results_df <- dplyr::bind_rows(valid_results_list)
  
  return(results_df)
}

