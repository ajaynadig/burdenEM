# sim_calibration.R: Calculate calibration of posterior-mean effects in simulations

# This script defines two main exported functions:
# 1. calculate_calibration_metrics(genes_df, model, ...)
#    – Returns a tibble summarising mean true vs. estimated effects in bins of
#      the posterior-mean estimate.
# 2. calculate_calibration_metrics_for_studies(studies_df, annotation_param, ...)
#    – Convenience wrapper that iterates over the rows of a studies.tsv data
#      frame, loads the appropriate genes file and fitted model for each row,
#      and binds the per-study calibration summaries.
#
# The helper is intended for use by simulations/sim_cli.R under the new
# "calibration" sub-command.

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
  library(tibble)
})

# We need posterior_expectation2 which lives in R/model.R.  This file is one
# directory up from simulations/, so we source it here (only if the symbol is
# missing to avoid duplicate sourcing when called repeatedly).
if (!exists("posterior_expectation2", mode = "function")) {
  source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "R", "model.R"))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b  # simple helper

# -----------------------------------------------------------------------------
# Core per-study computation ----------------------------------------------------
# -----------------------------------------------------------------------------

calculate_calibration_metrics <- function(genes_df,
                                          model,
                                          true_effect_col = "effectSize",
                                          per_allele_effects = FALSE,
                                          verbose = FALSE) {
  # Basic sanity check ---------------------------------------------------------
  if (!("gene" %in% names(genes_df))) stop("Column 'gene' not found in genes_df.")
  if (!("gene" %in% names(model$df))) stop("Column 'gene' not found in model$df.")
  if (!(true_effect_col %in% names(genes_df))) {
    stop(sprintf("Column '%s' not found in genes_df.", true_effect_col))
  }

  # Choose effect-size transformation -----------------------------------------
  effect_size_fn <- if (per_allele_effects) {
    function(x, row) model$to_per_allele_effects(row, x)
  } else {
    function(x, row) model$to_per_allele_effects(row, x) * sqrt(row$burden_score)
  }
  true_effect_size_fn <- if (per_allele_effects) {
    function(x, burden_score) x
  } else {
    function(x, burden_score) x * sqrt(burden_score)
  }

  # Compute posterior-mean effects per gene -----------------------------------
  model$df$posterior_mean <- posterior_expectation2(model, function_to_integrate = effect_size_fn)

  # Combine model estimates with truth via a merge ---------------------------
  model_df_for_join <- model$df %>% 
    select(gene, posterior_mean, burden_score) %>% 
    mutate(gene = as.character(gene))
  genes_df_for_join <- genes_df %>% 
    select(gene, !!sym(true_effect_col)) %>% 
    mutate(gene = as.character(gene))

  combined_df <- inner_join(model_df_for_join, genes_df_for_join, by = "gene") %>%
    rename(true_effect = !!sym(true_effect_col))

  # Define bins on the estimated effect ---------------------------------------
  quantile_probs <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  est_vals <- combined_df$posterior_mean[!is.na(combined_df$posterior_mean)]

  if (length(est_vals) > 0) {
    q_vals  <- quantile(est_vals, probs = quantile_probs, na.rm = FALSE, type = 7)
    breaks  <- unique(sort(c(-Inf, q_vals, Inf)))
  } else {
    breaks  <- c(-Inf, Inf)
  }

  combined_df <- combined_df %>%
    mutate(bin_factor = cut(posterior_mean, breaks = breaks, include.lowest = TRUE, right = TRUE)) %>%
    mutate(bin = paste0("bin", as.integer(bin_factor)))

  # Summarise per bin ----------------------------------------------------------
  summary_df <- combined_df %>%
    group_by(bin) %>%
    summarise(
      n_genes              = n(),
      mean_estimated_effect = mean(posterior_mean),
      se_estimated_effect   = sd(posterior_mean) / sqrt(n_genes),
      mean_true_effect      = mean(true_effect_size_fn(true_effect, burden_score)),
      se_true_effect        = sd(true_effect_size_fn(true_effect, burden_score)) / sqrt(n_genes),
      .groups = "drop"
    )

  return(summary_df)
}

# -----------------------------------------------------------------------------
# Wrapper over a studies.tsv data frame ---------------------------------------
# -----------------------------------------------------------------------------

calculate_calibration_metrics_for_studies <- function(studies_df,
                                                      annotation_param,
                                                      verbose_param = FALSE,
                                                      run_sequentially_param = TRUE,
                                                      per_allele_effects = FALSE) {
  # Load helper libraries for parallel-safe mapping
  if (!requireNamespace("furrr", quietly = TRUE)) install.packages("furrr")
  if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
  library(furrr)
  library(future)

  process_one_study <- function(study_row) {
    if (verbose_param) {
      message(sprintf("Processing calibration for study %s (dataset %s)",
                      study_row$identifier, study_row$dataset))
    }

    # ----- Locate model & genes files ---------------------------------------
    model_path <- stringr::str_replace(study_row$model_filename, "<ANNOTATION>", annotation_param)

    if (!file.exists(model_path)) {
      warning(sprintf("Model file not found for study %s: %s", study_row$identifier, model_path))
      return(NULL)
    }

    # genes file path derived from sumstats pattern as elsewhere -------------
    genes_filename <- gsub("sumstats", "genes", study_row$sumstats_filename_pattern, ignore.case = TRUE)
    if (!startsWith(genes_filename, "/") && !startsWith(genes_filename, "~")) {
      genes_filename <- file.path(dirname(studies_df$sumstats_filename_pattern[1]), genes_filename)
    }
    if (!file.exists(genes_filename)) {
      warning(sprintf("Genes file not found for study %s: %s", study_row$identifier, genes_filename))
      return(NULL)
    }

    # ----- Load data ---------------------------------------------------------
    model_obj  <- readRDS(model_path)
    genes_df   <- data.table::fread(genes_filename)

    # ----- Compute calibration metrics --------------------------------------
    calib_tbl  <- calculate_calibration_metrics(genes_df, model_obj, verbose = verbose_param, per_allele_effects = per_allele_effects)

    # Attach study identifier -------------------------------------------------
    calib_tbl$abbreviation <- study_row$abbreviation
    calib_tbl <- calib_tbl %>% dplyr::select(abbreviation, dplyr::everything())
    return(calib_tbl)
  }

  # Iterate over studies, optionally in parallel -----------------------------
  if (!run_sequentially_param && future::availableCores() > 1) {
    if (verbose_param) message("Running calibration in parallel mode")
    future::plan(future::multisession, workers = min(future::availableCores(), nrow(studies_df), 16))
    res <- furrr::future_map_dfr(1:nrow(studies_df), ~process_one_study(studies_df[.x, ]))
  } else {
    res <- purrr::map_dfr(1:nrow(studies_df), ~process_one_study(studies_df[.x, ]))
  }
  return(res)
}

meta_analyze_calibration <- function(calibration_df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for meta_analyze_calibration.")
  }

  # Summarize across replicates, grouping by abbreviation and bin
  summary_df <- calibration_df %>%
    dplyr::group_by(abbreviation, bin) %>% 
    dplyr::summarise(
      n_replicates = n(),
      post_mean.mean = mean(mean_estimated_effect, na.rm = FALSE),
      post_mean.sd = sd(mean_estimated_effect, na.rm = FALSE),
      true_mean.mean = mean(mean_true_effect, na.rm = FALSE),
      true_mean.sd = sd(mean_true_effect, na.rm = FALSE),
      .groups = "drop"
    )

  # Ensure original order of abbreviations and bins is preserved
  if (nrow(calibration_df) > 0) {
    ordered_keys <- calibration_df %>% 
      select(abbreviation, bin) %>% 
      distinct()
    
    summary_df <- summary_df %>%
      mutate(abbreviation = factor(abbreviation, levels = unique(ordered_keys$abbreviation)),
             bin = factor(bin, levels = unique(ordered_keys$bin))) %>%
      arrange(abbreviation, bin)
  }
  
  return(summary_df)
}

# End of sim_calibration.R
