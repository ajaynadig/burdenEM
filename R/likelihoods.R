#Functions to compute likelihoods for burdenEM inference


poisson_uniform_likelihood <- function(genetic_data,
                                       component_endpoints,
                                       grid_size) {
  no_cpts = length(component_endpoints)
  no_genes = nrow(genetic_data)

  likelihood <- matrix(NA,
                       nrow = no_genes,
                       ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)

  for (kk in 1:no_cpts) {
    # Expand case_count into a matrix by replicating the vector along the columns
    case_count_matrix <- replicate(grid_size,genetic_data$case_count)

    # Calculate the rate for each mu value in the grid and each case count
    rate <- genetic_data$expected_count * t(replicate(no_genes,exp(mu_grid * component_endpoints[kk])))

    # Compute the likelihood for each case count and rate combination
    likelihoods <- dpois(case_count_matrix, rate)

    # Calculate the average likelihood for each test across the grid of mu values
    likelihood[, kk] <- rowMeans(likelihoods)
  }
  return(likelihood)
}

