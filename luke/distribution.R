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

calculate_estimated_distribution <- function(study_row, current_annotation, current_verbose, gene_proportions, gene_counts=NULL) {
  model_filename_pattern <- study_row$model_filename
  actual_model_filename <- stringr::str_replace(model_filename_pattern, "<ANNOTATION>", current_annotation)
  model <- readRDS(actual_model_filename)

  if (is.null(gene_counts)) {
    gene_counts <- nrow(model$df) * gene_proportions
  }

  cumulative_variance_fn <- get_cumulative_variance_fn(model)
  
  return(data.frame(
    needed_genes = gene_counts,
    variance = cumulative_variance_fn(gene_proportions)
  ))
}

calculate_estimated_distribution_for_studies <- function(studies_df, annotation, gene_proportions, gene_counts = NULL, verbose = FALSE, run_sequentially_param = FALSE) {
  if (!run_sequentially_param) {
    if (verbose) {
      message(sprintf("Setting up parallel execution for distribution calculation using up to %d workers.", future::availableCores()))
    }
    future::plan(future::multisession)
    all_results_df <- furrr::future_map_dfr(
      .x = 1:nrow(studies_df),
      .f = ~{
        res <- calculate_estimated_distribution(studies_df[.x, , drop = FALSE], annotation, verbose, gene_proportions, gene_counts)
        res$abbreviation <- studies_df$abbreviation[.x]
        res$dataset <- studies_df$dataset[.x]
        return(res)
      },
      .progress = verbose,
      .options = furrr::furrr_options()
    )
  } else {
    if (verbose) {
      message("Running distribution calculation in sequential mode.")
    }
    all_results_df <- purrr::map_dfr(
      .x = 1:nrow(studies_df),
      .f = ~{
        res <- calculate_estimated_distribution(studies_df[.x, , drop = FALSE], annotation, verbose, gene_proportions, gene_counts)
        res$abbreviation <- studies_df$abbreviation[.x]
        res$dataset <- studies_df$dataset[.x]
        return(res)
      }
    )
  }
  
  all_results_df <- dplyr::select(all_results_df, abbreviation, dataset, dplyr::everything())
  if (verbose) message(sprintf("Distribution metrics calculation complete. Total rows: %d", nrow(all_results_df)))
  return(all_results_df)
}