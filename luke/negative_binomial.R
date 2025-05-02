library(glmmTMB)
library(dplyr)

estimate_overdispersion <- function(variant_data, intercept_frequency_bin_edges) {
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
  epsilon <- 1e-10

  # Fit the negative binomial model with a single dispersion parameter
  # Include 'gene' as a random intercept effect
  model <- glmmTMB(
    formula = AC_cases ~ offset(log(expected_count + epsilon)),
    family = nbinom2,
    dispformula = ~ 1, # Use a single dispersion parameter for all bins
    data = variant_data
  )

  # Extract the estimated dispersion (k, the size parameter)
  # sigma() returns the size parameter 'k' for nbinom2
  # For a single dispersion parameter (~1), sigma returns a single value
  k_value <- sigma(model)
  
  # Calculate overdispersion alpha = 1/k
  # Handle potential division by zero or issues if k is not positive
  overdispersion <- ifelse(k_value > 0, 1 / k_value, NA_real_)
  
  # Check if overdispersion calculation resulted in NA
  if (is.na(overdispersion)) {
      stop("Overdispersion estimation failed or resulted in NA.")
  }
  
  # Since we have only one overdispersion value, assign it to all variants
  variant_data$overdispersion <- overdispersion
  
  return(variant_data)
}


#' Calculate Negative Binomial Log-Likelihoods Across a Grid and Add to DataFrame
#'
#' Calculates the negative binomial log-likelihood for each variant across a grid
#' of beta values (log odds ratio) and adds these likelihoods as columns to the
#' input variant dataframe.
#'
#' @param variant_df A data frame containing variant-level data. Must include
#'   columns: 'AC_cases', 'expected_count', 'overdispersion'.
#' @param grid A numeric vector representing the grid of beta values (log ORs)
#'   over which to calculate likelihoods.
#'
#' @return The input `variant_df` with additional columns corresponding to the
#'   log-likelihood calculated for each beta value in the `grid`. Column names
#'   will be structured like "LL_betaValue", e.g., "LL_-1.5".
#' @importFrom stats dnbinom
#' @export
add_nbinom_likelihood <- function(variant_df, grid) {

    # --- Input Checks ---
    required_cols <- c("AC_cases", "expected_count", "overdispersion")
    missing_cols <- setdiff(required_cols, names(variant_df))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns in variant_df:", paste(missing_cols, collapse=", ")))
    }
    if (!is.numeric(grid)) {
        stop("'grid' must be a numeric vector.")
    }
    # Check for non-finite values in necessary input columns
    # Allowing NaNs as per user preference, but Inf/-Inf in overdispersion/expected_count can be problematic
    if(any(!is.finite(variant_df$overdispersion[!is.na(variant_df$overdispersion)]))) {
        warning("Non-finite values detected in 'overdispersion'. This might lead to issues in likelihood calculation.")
    }
     if(any(!is.finite(variant_df$expected_count[!is.na(variant_df$expected_count)]))) {
        warning("Non-finite values detected in 'expected_count'. This might lead to issues in likelihood calculation.")
    }


    # --- Helper function to calculate LL for one beta value (vectorized over variants) ---
    log_likelihood_variant_beta <- function(beta_val, expected_count_vec, AC_cases_vec, overdispersion_vec) {
        # Calculate theta (size parameter for dnbinom) from overdispersion
        # Handle overdispersion <= 0 (theta -> Inf, approaches Poisson)
        theta_vec <- ifelse(overdispersion_vec <= 0 | !is.finite(overdispersion_vec), Inf, 1 / overdispersion_vec)

        # Calculate mu (mean parameter for dnbinom)
        mu_vec <- expected_count_vec * exp(beta_val)

        # Calculate log-likelihood using dnbinom (vectorized)
        # Ensure inputs to dnbinom are valid
        # Negative mu can occur if expected_count is negative, which shouldn't happen but check.
        # size (theta) must be positive. Our ifelse handles the zero case -> Inf.
        if(any(expected_count_vec < 0, na.rm=TRUE)) {
            warning("Negative values detected in 'expected_count'. This will likely cause errors in dnbinom.")
        }

        ll <- stats::dnbinom(AC_cases_vec, mu = mu_vec, size = theta_vec, log = TRUE)

        # dnbinom returns -Inf for log=TRUE if probability is 0.
        # It can return NaN if inputs are invalid (e.g., negative size, negative mu depending on implementation).
        # No specific replacement needed here, let R handle standard numeric results.
        return(ll)
    }

    # --- Apply across the grid ---
    # sapply iterates through each beta value in 'grid'
    # For each beta, it calls log_likelihood_variant_beta, passing the *vectors* from variant_df
    # The result is a matrix: rows = variants, cols = grid points
    message(paste("Calculating likelihoods for", nrow(variant_df), "variants across", length(grid), "grid points..."))
    likelihood_matrix <- sapply(grid, function(beta_val) {
        log_likelihood_variant_beta(
            beta_val = beta_val,
            expected_count_vec = variant_df$expected_count,
            AC_cases_vec = variant_df$AC_cases,
            overdispersion_vec = variant_df$overdispersion
        )
    }) # Result is variants x grid_points matrix

    # --- Combine with input dataframe ---
    # Create informative column names for the likelihood matrix
    colnames(likelihood_matrix) <- paste0("LL_", format(grid, digits=3, scientific=FALSE)) # e.g., LL_-1.500, LL_0.000, LL_1.500

    # Check row count match before cbind
    if (nrow(variant_df) != nrow(likelihood_matrix)) {
        stop("Internal error: Row count mismatch between variant data and calculated likelihood matrix.")
    }

    # Combine the original dataframe with the new likelihood columns
    variant_df_with_ll <- cbind(variant_df, likelihood_matrix)

    message("Likelihood columns added to the dataframe.")
    return(variant_df_with_ll)
}
