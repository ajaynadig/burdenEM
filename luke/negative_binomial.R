library(glmmTMB)
library(dplyr)

# Define log-likelihood function
log_likelihood_fn <- function(expected_count, AC_case, beta, alpha) {
  theta <- 1 / alpha
  mu <- expected_count * exp(beta)
  sum(dnbinom(AC_case, mu = mu, size = theta, log = TRUE))
}

# Define grid_likelihood function
grid_likelihood <- function(variant_df, frequency_bin_edges, num_cases, grid) {
  # Ensure required columns are present
  required_cols <- c("AF", "gene", "AC_case")
  missing_cols <- setdiff(required_cols, colnames(variant_df))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from 'variant_df': ", paste(missing_cols, collapse = ", "))
  }

  # Create frequency bins
  variant_df <- variant_df %>%
    mutate(freq_bin = cut(AF, breaks = frequency_bin_edges, include.lowest = TRUE, right = FALSE))

  # Check for NA bins
  if (any(is.na(variant_df$freq_bin))) {
    stop("Some allele frequencies fall outside the specified 'frequency_bin_edges'. Please adjust the bin edges.")
  }

  # Group by gene and compute likelihoods
  likelihood_matrix <- variant_df %>%
    group_by(gene) %>%
    summarize(log_likelihood = sum(log_likelihood_fn(expected_count, AC_case, beta = grid, alpha))) %>%
    ungroup() %>%
    arrange(gene) %>%
    pull(log_likelihood) %>%
    matrix(ncol = length(grid), byrow = TRUE)

  colnames(likelihood_matrix) <- paste0("beta_", grid)

  return(likelihood_matrix)
}
