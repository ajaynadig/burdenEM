library(glmmTMB)
library(dplyr)

estimate_overdispersion <- function(variant_data, intercept_frequency_bin_edges, verbose = FALSE) {
  # Ensure required columns are present
  required_cols <- c("AC_cases", "expected_count", "AF", "gene")
  if (!all(required_cols %in% names(variant_data))) {
    stop("Input variant_data must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Ensure 'gene' is treated as a factor
  if (!is.factor(variant_data$gene)) {
      variant_data$gene <- as.factor(variant_data$gene)
  }

  # Define frequency bins based on provided edges
  variant_data <- variant_data %>%
    mutate(
      freq_bin = cut(AF, breaks = intercept_frequency_bin_edges, include.lowest = TRUE, right = FALSE)
    )

  # Add a small epsilon for stability if expected_count can be zero
#   epsilon <- 1e-10

#   # Fit the negative binomial model with a single dispersion parameter
#   model <- glmmTMB(
#     formula = AC_cases ~ offset(log(expected_count + epsilon)),
#     family = nbinom2,
#     dispformula = ~ 1, # Use a single dispersion parameter for all bins
#     data = variant_data
#   )
#   k_value <- sigma(model)
  
#   # Calculate overdispersion alpha = 1/k
#   overdispersion <- ifelse(k_value > 0, 1 / k_value, NA_real_)
  
#   # Check if overdispersion calculation resulted in NA
#   if (is.na(overdispersion)) {
#       stop("Overdispersion estimation failed or resulted in NA.")
#   }

  variant_data$overdispersion <- NA_real_
  unique_freq_bins <- levels(variant_data$freq_bin)
  for (current_bin_name in unique_freq_bins) {
    in_bin <- variant_data$freq_bin == current_bin_name
    beta <- mean(variant_data$AC_cases[in_bin] * variant_data$expected_count[in_bin]) / 
            mean(variant_data$expected_count[in_bin]^2)
    mu <- beta * variant_data$expected_count[in_bin]
    residual <- (variant_data$AC_cases[in_bin] - mu)^2 - mu
    overdispersion <- mean(residual * mu^2) / mean(mu^4)
    
    if (verbose) {
        cat("Freq bin:", as.character(current_bin_name), "\n")
        cat("  Calibration slope (beta) :", beta, "\n")
        cat("  Moment-based overdispersion (alpha) before clamping at 0:", overdispersion, "\n")
    }
    overdispersion <- max(overdispersion, 0)
    variant_data$overdispersion[in_bin] <- overdispersion
  }
    
    return(variant_data)
}

#' Aggregate Variant Data to Gene-Level Lists
#'
#' Groups variant-level data by gene and functional category, and summarizes
#' key columns (AC_cases, expected_count, overdispersion) into list-columns,
#' where each element of the list is a vector of values for the variants in that gene.
#' It also calculates n_variants and brings forward the first N and prevalence.
#'
#' @param variant_data A data frame containing variant-level information.
#'   Must include 'gene', 'functional_category', 'AC_cases', 'expected_count',
#'   'overdispersion', 'N', and 'prevalence'.
#' @param verbose Logical, if TRUE, prints progress messages.
#' @return A data frame aggregated by gene and functional_category, with
#'   list-columns for 'AC_cases', 'expected_count', 'overdispersion',
#'   and columns for 'n_variants', 'N', 'prevalence'.
aggregate_variants_to_gene_lists <- function(variant_data, verbose = FALSE) {
  if (verbose) message("Aggregating variant data to gene-level lists...")

  required_cols <- c("gene", "functional_category", "AC_cases", "expected_count", "overdispersion", "N", "prevalence")
  missing_cols <- setdiff(required_cols, names(variant_data))
  if (length(missing_cols) > 0) {
    stop(paste("variant_data is missing one or more required columns for list aggregation:",
               paste(missing_cols, collapse = ", ")))
  }

  gene_level_data <- variant_data %>%
    dplyr::group_by(gene, functional_category) %>%
    dplyr::summarise(
      AC_cases = list(AC_cases),
      expected_count = list(expected_count),
      overdispersion = list(overdispersion),
      n_variants = dplyr::n(),
      N = dplyr::first(N),
      prevalence = dplyr::first(prevalence), # Assumes prevalence is constant
      burden_score = sum(variant_variance),
      .groups = 'drop'
    ) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::mutate(
      CAC_cases = sapply(AC_cases, function(x) sum(unlist(x))),
      expected_CAC_cases = sapply(expected_count, function(x) sum(unlist(x)))
    )

  if (verbose) message(paste("Aggregated to", nrow(gene_level_data), "gene-category rows with list-columns."))
  return(gene_level_data)
}

#' Calculate Gene-Level Likelihoods for Negative Binomial Model
#'
#' Calculates the sum of negative binomial log-likelihoods for all variants in a gene
#' across a grid of effect sizes (beta_vec).
#' This function is designed to be used as the `likelihood_fn` in `fit_burdenem_model`.
#'
#' @param gene_data_row A single row from the gene-level data frame. Expected to have
#'   list-columns: `AC_cases` (list of allele counts in cases for variants in the gene),
#'   `expected_count` (list of expected counts for variants in the gene under H0),
#'   `overdispersion` (list of overdispersion parameters for variants in the gene).
#' @param beta_vec A numeric vector of effect sizes (log rate ratios) to evaluate.
#' @return A numeric vector of total actual likelihoods (sum of variant likelihoods, exponentiated)
#'   for the gene, one for each value in `beta_vec`.
gene_likelihood <- function(gene_data_row, beta_vec) {
  # Extract lists from the gene_data_row
  ac_cases <- gene_data_row$AC_cases[[1]]
  ac_expected <- gene_data_row$expected_count[[1]]
  overdispersion <- gene_data_row$overdispersion[[1]]

  # Calculate sum of log-likelihoods for each beta
  total_sum_log_likelihood_per_beta <- sapply(beta_vec, function(current_beta) {
    mu <- ac_expected * exp(current_beta)
    theta <- ifelse(overdispersion <= 0 | !is.finite(overdispersion), Inf, 1 / overdispersion)

    log_likelihoods_for_variants <- stats::dnbinom(
      x = ac_cases,
      mu = mu,
      size = theta,
      log = TRUE
    )
    # Return the sum of log-likelihoods for this beta
    sum(log_likelihoods_for_variants, na.rm = TRUE) 
  })

  # Normalize to prevent underflow/overflow before exponentiating
  max_log_likelihood <- max(total_sum_log_likelihood_per_beta, na.rm = TRUE)

  # If max_log_likelihood is -Inf (all log-likelihoods were -Inf or NA),
  # then all exponentiated likelihoods will be 0.
  # If max_log_likelihood is NA (e.g. all inputs to max were NA), result will be NA vector.
  if (is.na(max_log_likelihood) || !is.finite(max_log_likelihood)) {
    # If max is -Inf, exp(-Inf - (-Inf)) is NaN, should be 0.
    # If max is NA, exp(NA - NA) is NA.
    if (identical(max_log_likelihood, -Inf)) {
        return(rep(0, length(total_sum_log_likelihood_per_beta)))
    } else {
        # This will correctly produce a vector of NAs if max_log_likelihood was NA
        return(exp(total_sum_log_likelihood_per_beta - max_log_likelihood)) 
    }
  }

  normalized_log_likelihoods <- total_sum_log_likelihood_per_beta - max_log_likelihood
  final_likelihoods <- exp(normalized_log_likelihoods)
  
  return(final_likelihoods)
}

