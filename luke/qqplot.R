#!/usr/bin/env Rscript

# --- Libraries ---
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(future))

# Source model functions
source("R/model.R")

#' Compute QQ plot data for observed effect sizes vs expected from the model
#'
#' For each gene, computes:
#' - Z-score from observed effect estimate (continuous) or case counts (binary)
#' - Posterior mean effect size from the fitted model
#' - Expected quantiles for both
#'
#' @param model Fitted burdenEM model object
#' @param verbose Logical for verbose output
#' @return Data frame with QQ plot data
compute_qqplot_data_for_model <- function(model, verbose = FALSE) {
    if (is.null(model) || is.null(model$df)) {
        if (verbose) message("  Model or model$df is NULL, returning empty data frame.")
        return(tibble())
    }
    
    # Convert to regular data frame if needed
    df <- as.data.frame(model$df)

    # Detect trait type based on available columns
    is_binary <- ("CAC_cases" %in% names(df)) && ("expected_CAC_cases" %in% names(df))
    is_continuous <- ("effect_estimate" %in% names(df)) || ("gamma_perSD" %in% names(df))

    if (is_binary) {
        stop("compute_qqplot_data_for_model() only supports quantitative traits; binary trait models are not supported.")
    }

    # Define functions for posterior expectation
    effect_size_fn <- function(x, row) x
    effect_sq_fn <- function(x, row) x^2

    # Compute posterior mean effect sizes for each gene
    model$df$posterior_mean <- posterior_expectation2(model, function_to_integrate = effect_size_fn)
    model$df$posterior_mean_sq <- posterior_expectation2(model, function_to_integrate = effect_sq_fn)
    model$df$posterior_variance <- model$df$posterior_mean_sq - (model$df$posterior_mean^2)
    model$df$posterior_sd <- sqrt(pmax(model$df$posterior_variance, 0))
    df <- as.data.frame(model$df)

    if (is_continuous) {
        # Continuous trait: compute Z-scores
        effect_col <- if ("effect_estimate" %in% names(df)) "effect_estimate" else "gamma_perSD"
        df$effect_estimate <- df[[effect_col]]
        
        # Check if we have mean_n and mean_variant_intercept for proper scaling (like R/qqplot.R)
        has_scaling_cols <- ("mean_n" %in% names(df)) && ("mean_variant_intercept" %in% names(df))
        
        # Also check for effect_se as fallback
        se_col <- if ("effect_se" %in% names(df)) "effect_se" else if ("gamma_perSD_se" %in% names(df)) "gamma_perSD_se" else NULL
        
        if (has_scaling_cols) {
            # Use R/qqplot.R formula: Z = gamma_perSD * sqrt(mean_n / mean_variant_intercept)
            scaling_factor <- sqrt(df$mean_n / df$mean_variant_intercept)
            df$Z_observed <- df$effect_estimate * scaling_factor
            df$Z_posterior_mean <- df$posterior_mean * scaling_factor
            df$effect_se <- 1 / scaling_factor
            if (verbose) message("  Detected continuous trait model (using mean_n/mean_variant_intercept)")
        } else if (!is.null(se_col)) {
            # Fallback: use effect_estimate / effect_se
            df$effect_se <- df[[se_col]]
            df$Z_observed <- df$effect_estimate / df$effect_se
            df$Z_posterior_mean <- df$posterior_mean / df$effect_se
            if (verbose) message("  Detected continuous trait model (using effect_se)")
        } else {
            if (verbose) message("  Required columns for continuous trait not found (need mean_n/mean_variant_intercept or effect_se).")
            return(tibble())
        }
        
        trait_type <- "continuous"
        
    } else {
        if (verbose) message("  Could not detect trait type (missing required columns).")
        return(tibble())
    }
    
    n_genes <- nrow(df)
    
    # Ensure burden_score column exists
    if (!("burden_score" %in% names(df))) {
        df$burden_score <- NA_real_
    }
    
    # For proper QQ plots, EACH series should be sorted by itself
    # This ensures both observed and posterior lines are smooth
    
    # QQ data for observed Z-scores (sorted by observed)
    qq_observed <- df %>%
        arrange(Z_observed) %>%
        mutate(
            rank = row_number(),
            expected_quantile = qnorm((rank - 0.5) / n_genes),
            trait_type = trait_type,
            Z_value = Z_observed,
            type = "observed",
            value_type = "observed_Z"
        ) %>%
        select(gene, expected_quantile, rank, trait_type, effect_estimate, effect_se,
               burden_score, posterior_mean, posterior_sd, Z_value, type, value_type)
    
    # QQ data for posterior mean Z-scores (sorted by posterior)
    qq_posterior <- df %>%
        arrange(Z_posterior_mean) %>%
        mutate(
            rank = row_number(),
            expected_quantile = qnorm((rank - 0.5) / n_genes),
            trait_type = trait_type,
            Z_value = Z_posterior_mean,
            type = "posterior_mean",
            value_type = "posterior_mean_Z"
        ) %>%
        select(gene, expected_quantile, rank, trait_type, effect_estimate, effect_se,
               burden_score, posterior_mean, posterior_sd, Z_value, type, value_type)
    
    # Combine both - now they share the same expected_quantile for each gene
    qq_combined <- bind_rows(qq_observed, qq_posterior) %>%
        rename(Z_observed = Z_value)
    
    return(qq_combined)
}


#' Generate expected QQ quantiles from the fitted model via simulation
#'
#' Simulates effect sizes from the prior distribution and adds noise
#' to generate expected Z-scores for QQ comparison.
#'
#' @param model Fitted burdenEM model
#' @param n_genes Number of genes (for quantile computation)
#' @param n_draws Number of Monte Carlo draws
#' @return Data frame with expected quantiles
compute_expected_quantiles_from_model <- function(model, n_genes, n_draws = 1e6) {
    # Get average sample size and variance intercept from model data
    mean_n <- mean(model$df$mean_n)
    mean_intercept <- mean(model$df$mean_variant_intercept)
    
    # Get component weights (average across genes)
    features <- do.call(rbind, model$df$features)
    gene_component_contributions <- features %*% model$delta
    weights <- colMeans(gene_component_contributions)
    
    # Generate draws from the effect size distribution
    draws_per_component <- as.integer(floor(n_draws * weights))
    
    draws <- unlist(lapply(seq_along(model$component_endpoints), function(i) {
        if (draws_per_component[i] > 0) {
            # Draw uniformly from [0, endpoint] for each component
            model$component_endpoints[i] * runif(draws_per_component[i])
        } else {
            numeric(0)
        }
    }))
    
    # Add measurement noise: observed Z = true_effect * sqrt(n/intercept) + noise
    n_draws_actual <- length(draws)
    noise_sd <- 1  # Standard normal noise after scaling
    Z_simulated <- draws * sqrt(mean_n / mean_intercept) + rnorm(n_draws_actual, sd = noise_sd)
    
    # Compute expected quantiles at the same ranks as the data
    probs <- ((1:n_genes) - 0.5) / n_genes
    expected_Z <- quantile(Z_simulated, probs = probs)
    
    expected_df <- tibble(
        rank = 1:n_genes,
        expected_quantile = qnorm(probs),
        expected_Z_from_model = as.numeric(expected_Z),
        type = "model_expected"
    )
    
    return(expected_df)
}


#' Process a single study for QQ plot computation
#'
#' @param study_row A single row from studies_df
#' @param annotation The annotation to use
#' @param verbose Logical for verbose output
#' @return Data frame with QQ plot data for this study
process_qqplot_for_study <- function(study_row, annotation, verbose = FALSE) {
    if (verbose) {
        message(sprintf("Processing QQ plot for study: %s (%s) with annotation '%s'",
                        study_row$abbreviation, study_row$dataset, annotation))
    }
    
    model_path <- stringr::str_replace(study_row$model_filename, "<ANNOTATION>", annotation)
    
    model <- NULL
    error_msg <- NA_character_
    
    tryCatch({
        model <- readRDS(model_path)
    }, error = function(e) {
        error_msg <<- sprintf("Error loading model '%s': %s", model_path, e$message)
        if (verbose) message(paste("  ", error_msg))
    })
    
    if (is.null(model)) {
        return(tibble(
            abbreviation = study_row$abbreviation,
            dataset = study_row$dataset,
            annotation = annotation,
            gene = NA_character_,
            Z_observed = NA_real_,
            expected_quantile = NA_real_,
            rank = NA_integer_,
            type = NA_character_,
            trait_type = NA_character_,
            value_type = NA_character_,
            effect_estimate = NA_real_,
            effect_se = NA_real_,
            burden_score = NA_real_,
            posterior_mean = NA_real_,
            posterior_sd = NA_real_,
            error_message = error_msg
        ))
    }
    
    # Compute QQ data
    qq_data <- compute_qqplot_data_for_model(model, verbose = verbose)
    
    if (nrow(qq_data) == 0) {
        return(tibble(
            abbreviation = study_row$abbreviation,
            dataset = study_row$dataset,
            annotation = annotation,
            gene = NA_character_,
            Z_observed = NA_real_,
            expected_quantile = NA_real_,
            rank = NA_integer_,
            type = NA_character_,
            trait_type = NA_character_,
            value_type = NA_character_,
            effect_estimate = NA_real_,
            effect_se = NA_real_,
            burden_score = NA_real_,
            posterior_mean = NA_real_,
            posterior_sd = NA_real_,
            error_message = "No QQ data computed"
        ))
    }
    
    # Add study identifiers
    qq_data <- qq_data %>%
        mutate(
            abbreviation = study_row$abbreviation,
            dataset = study_row$dataset,
            annotation = annotation,
            error_message = NA_character_
        ) %>%
        select(abbreviation, dataset, annotation, everything())
    
    return(qq_data)
}


#' Main function to calculate QQ plot data for all studies
#'
#' @param studies_df Data frame with study information
#' @param annotation The annotation to use
#' @param verbose Logical for verbose output
#' @param run_sequentially_param If TRUE, run sequentially instead of in parallel
#' @return Data frame with QQ plot data for all studies
calculate_qqplot_for_studies <- function(studies_df, annotation, verbose = FALSE, run_sequentially_param = FALSE) {
    if (!is.data.frame(studies_df) || nrow(studies_df) == 0) {
        if (verbose) message("Input studies_df is empty or not a dataframe. Returning empty results.")
        return(tibble())
    }
    
    if (!run_sequentially_param) {
        if (verbose) {
            message(sprintf("Setting up parallel execution for QQ plot calculation using up to %d workers.",
                            future::availableCores()))
        }
        future::plan(future::multisession)
        
        all_results <- furrr::future_map_dfr(
            .x = 1:nrow(studies_df),
            .f = ~process_qqplot_for_study(studies_df[.x, , drop = FALSE], annotation, verbose),
            .progress = verbose,
            .options = furrr::furrr_options()
        )
    } else {
        if (verbose) {
            message("Running QQ plot calculation in sequential mode.")
        }
        all_results <- purrr::map_dfr(
            .x = 1:nrow(studies_df),
            .f = ~process_qqplot_for_study(studies_df[.x, , drop = FALSE], annotation, verbose)
        )
    }
    
    if (verbose) {
        message(sprintf("QQ plot calculation complete. Total rows: %d", nrow(all_results)))
    }
    
    return(all_results)
}
