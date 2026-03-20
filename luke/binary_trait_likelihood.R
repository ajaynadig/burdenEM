library(dplyr)
library(VGAM)

estimate_overdispersion_poisson <- function(variant_data, intercept_frequency_bin_edges, verbose = FALSE) {
  # Ensure required columns are present
  required_cols <- c("AC_cases", "expected_count", "AF", "gene", "dataset")
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

  # Overdispersion is estimated per (freq_bin, dataset). Different cohorts have
  # different sequencing platforms, sample ascertainment, and population structure
  # residuals, so pooling across cohorts biases the estimate. This mirrors how
  # continuous-trait intercepts are computed per-dataset in calculate_variant_intercept().
  # For single-cohort studies this reduces to the same computation as before.
  variant_data$overdispersion <- NA_real_
  unique_freq_bins <- levels(variant_data$freq_bin)
  for (current_bin_name in unique_freq_bins) {
    in_bin <- which(variant_data$freq_bin == current_bin_name)
    for (ds in unique(variant_data$dataset[in_bin])) {
      idx <- in_bin[variant_data$dataset[in_bin] == ds]
      if (length(idx) == 0) next
      beta_cal <- mean(variant_data$AC_cases[idx] * variant_data$expected_count[idx]) /
              mean(variant_data$expected_count[idx]^2)
      mu <- beta_cal * variant_data$expected_count[idx]
      residual <- (variant_data$AC_cases[idx] - mu)^2 - mu
      overdispersion <- mean(residual * mu^2) / mean(mu^4)

      if (verbose) {
          cat("Freq bin:", as.character(current_bin_name), "[", ds, "]\n")
          cat("  Calibration slope (beta) :", beta_cal, "\n")
          cat("  Moment-based overdispersion (alpha) before clamping at 0:", overdispersion, "\n")
      }
      overdispersion <- max(overdispersion, 0)
      variant_data$overdispersion[idx] <- overdispersion
    }
  }

    return(variant_data)
}

estimate_overdispersion_binomial <- function(variant_data, intercept_frequency_bin_edges, verbose = FALSE) {
  # Ensure required columns are present
  required_cols <- c("AC_cases", "expected_count", "AF", "gene", "prevalence", "dataset")
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

  # Per-dataset estimation (see note in estimate_overdispersion_poisson).
  variant_data$overdispersion <- NA_real_
  unique_freq_bins <- levels(variant_data$freq_bin)
  for (current_bin_name in unique_freq_bins) {
    in_bin <- which(variant_data$freq_bin == current_bin_name)
    for (ds in unique(variant_data$dataset[in_bin])) {
      idx <- in_bin[variant_data$dataset[in_bin] == ds]
      if (length(idx) == 0) next
      n <- variant_data$expected_count[idx] / variant_data$prevalence[idx]
      p <- variant_data$prevalence[idx]
      npq <- n * p * (1 - p)
      residual <- (variant_data$AC_cases[idx] - n*p)^2 - npq
      rho <- mean(residual * npq * (n-1)) / mean((npq * (n-1))^2)
      if(is.na(rho)) rho <- 0

      if (verbose) {
          cat("Freq bin:", as.character(current_bin_name), "[", ds, "]\n")
          cat("  Between-class correlation (rho) before clamping at [0,1):", rho, "\n")
      }
      rho <- min(1-1e-9, max(rho, 0))
      variant_data$overdispersion[idx] <- rho
    }
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
      burden_score = if("variant_variance" %in% names(variant_data)) sum(variant_variance) else NA_real_,
      # Pair each variant's variance with its own cohort's N and P² before summing (must precede list() conversions)
      h2_numer_scale = if("variant_variance" %in% names(variant_data)) sum(variant_variance * N * prevalence^2) else NA_real_,
      N = list(N),
      prevalence = list(prevalence),
      .groups = 'drop'
    ) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::mutate(
      CAC_cases = sapply(AC_cases, function(x) sum(unlist(x))),
      expected_CAC_cases = sapply(expected_count, function(x) sum(unlist(x))),
      # Pooled prevalence for rate_ratio cap
      mean_prevalence = mapply(function(ec, prev) sum(ec) / sum(ec / prev), expected_count, prevalence),
      h2_denom = mapply(function(n, prev) {
        pop <- unique(data.frame(N = n, prevalence = prev))
        sum(pop$N * pop$prevalence * (1 - pop$prevalence))
      }, N, prevalence)
    )

  if (verbose) message(paste("Aggregated to", nrow(gene_level_data), "gene-category rows with list-columns."))
  return(gene_level_data)
}

# Exponentiate a vector of log-likelihoods up to multiplication by a constant without causing underflow
stable_exp <- function(x){
  if(length(x)==0) return(x)
  xmax <- max(x, na.rm=TRUE)
  if (!is.infinite(xmax)) {
    return(exp(x - xmax))
  }
  if(xmax==Inf) stop("Argument to stable_exp contains Inf values")
  # All values are -Inf or NA: gene has zero likelihood under all grid effects.
  stop("stable_exp: all log-likelihood values are -Inf or NA. This should not happen with properly preprocessed data.")
}


poisson_log_likelihood <- function(gene_data_row, beta_vec) {
  # Extract lists from the gene_data_row
  ac_cases <- gene_data_row[["AC_cases"]][[1]]
  ac_expected <- gene_data_row[["expected_count"]][[1]]

  # Calculate sum of log-likelihoods for each beta
  total_sum_log_likelihood_per_beta <- sapply(beta_vec, function(current_beta) {
    mu <- ac_expected * exp(current_beta)

    log_likelihoods_for_variants <- stats::dpois(
      x = ac_cases,
      lambda = mu,
      log = TRUE
    )
    # Return the sum of log-likelihoods for this beta
    sum(log_likelihoods_for_variants)
  })

  return(total_sum_log_likelihood_per_beta)
}

negative_binomial_log_likelihood <- function(gene_data_row, beta_vec) {
  # Extract lists from the gene_data_row
  ac_cases <- gene_data_row[["AC_cases"]][[1]]
  ac_expected <- gene_data_row[["expected_count"]][[1]]
  overdispersion <- gene_data_row[["overdispersion"]][[1]]

  # Calculate sum of log-likelihoods for each beta
  total_sum_log_likelihood_per_beta <- sapply(beta_vec, function(current_beta) {
    mu <- ac_expected * exp(current_beta)
    theta <- ifelse(overdispersion <= 0 | !is.finite(overdispersion), Inf, 1 / overdispersion)

    log_likelihoods_for_variants <- stats::dnbinom(
      x = ac_cases,
      mu = pmin(mu, theta),
      size = theta,
      log = TRUE
    )
    # Return the sum of log-likelihoods for this beta
    sum(log_likelihoods_for_variants)
  })

  return(total_sum_log_likelihood_per_beta)
}

binomial_log_likelihood <- function(gene_data_row, beta_vec) {
  # Extract lists from the gene_data_row
  cac_cases <- gene_data_row[["AC_cases"]][[1]]
  cac_total <- round(gene_data_row[["expected_count"]][[1]] / gene_data_row[["prevalence"]][[1]])
  prevalence <- gene_data_row[["prevalence"]][[1]]

  # Calculate sum of log-likelihoods for each beta
  total_sum_log_likelihood_per_beta <- sapply(beta_vec, function(current_beta) {

    p <- prevalence * exp(current_beta)
    log_likelihoods_for_variants <- stats::dbinom(
      x = cac_cases,
      size = cac_total,
      prob = pmin(p, 1),
      log = TRUE
    )
    # Return the sum of log-likelihoods for this beta
    sum(log_likelihoods_for_variants)
  })

  return(total_sum_log_likelihood_per_beta)
}

beta_binomial_log_likelihood <- function(gene_data_row, beta_vec) {
  # Extract lists from the gene_data_row
  cac_cases <- gene_data_row[["AC_cases"]][[1]]
  cac_total <- round(gene_data_row[["expected_count"]][[1]] / gene_data_row[["prevalence"]][[1]])
  prevalence <- gene_data_row[["prevalence"]][[1]]
  overdispersion <- gene_data_row[["overdispersion"]][[1]]

  # Calculate sum of log-likelihoods for each beta
  total_sum_log_likelihood_per_beta <- sapply(beta_vec, function(current_beta) {

    p <- prevalence * exp(current_beta)
    log_likelihoods_for_variants <- dbetabinom(
      x = cac_cases,
      size = cac_total,
      prob = pmin(p, .9999), # function gives NaN for p=1
      rho = overdispersion,
      log = TRUE
    )
    # Return the sum of log-likelihoods for this beta
    sum(log_likelihoods_for_variants)
  })

  return(total_sum_log_likelihood_per_beta)
}


#' Get a gene-level likelihood function for a binary trait
#'
#' @param model_type Character string, the model type to use. Valid options are 
#' 'poisson', 'negative_binomial', 'binomial', 'beta_binomial'.
#' @param effect_size_function A function that computes effect sizes in the appropriate units from a vector of effects in possibly other units and a row of the gene dataframe
#' @return A function that operates on (1) a row of gene-level dataframe model$df and (2) a vector of effect sizes
#' @export
get_gene_likelihood <- function(model_type, to_per_allele_effects){
  if (model_type == "pois"){
    f = poisson_log_likelihood
  }else if (model_type == "nbinom"){
    f = negative_binomial_log_likelihood
  }else if (model_type == "binom"){
    f = binomial_log_likelihood
  }else if (model_type == "betabinom"){
    f = beta_binomial_log_likelihood
  }else{
    stop("Valid model options are 'pois', 'nbinom', 'binom', 'betabinom'")
  }

  result <- function(gene_data_row, beta_vec){
    beta_per_allele <- to_per_allele_effects(gene_data_row, beta_vec)
    stable_exp(f(gene_data_row, beta_per_allele))
  }
    
}
