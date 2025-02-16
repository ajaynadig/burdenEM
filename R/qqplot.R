compute_estimate_expected_rvas <- function(input_data,
                                           model) {
  if(is.null(input_data$annotation)){
    input_data$annotation <- 1
    print('Plotting all data points together as a single group...')
  }

  no_genes <- nrow(input_data)
  no_draws <- 1e7
  no_samples <- unique(input_data$mean_n)[1]
  mean_intercept_nn <- unique(input_data$mean_variant_intercept)[1]
  weights_agg <- colSums(model$delta)/5

  # Generate draws in a vectorized way
  draws <- unlist(sapply(1:length(model$component_endpoints),
                         function(i){
                           model$component_endpoints[i] * runif(as.integer(floor(no_draws * weights_agg[i])))
                           }))

  # Add noise
  no_draws <- length(draws)
  draws <- draws + rnorm(no_draws) * sqrt(mean_intercept_nn / no_samples)


  # Compute estimated quantiles
  est <- quantile(draws * sqrt(no_samples), ((1:no_genes) / no_genes) - 0.5 / no_genes)

  est_data <- data.frame(observed=est) %>%
    dplyr::arrange(observed) %>%
    dplyr::mutate(index = row_number(),
                  n_genes = n()) %>%
    dplyr::mutate(expected = qnorm(index/n_genes))
  return(est_data)
}

compute_true_expected_rvas <- function(input_data){
  if(is.null(input_data$annotation)){
    input_data$annotation <- 1
    print('Plotting all data points together as a single group...')
  }
  qq_data <- input_data %>%
    dplyr::group_by(annotation) %>%
    dplyr::mutate(
      Z = gamma_perSD*sqrt(mean_n/mean_variant_intercept),
      n_genes = n()) %>%
    dplyr::group_by(annotation) %>%
    dplyr::arrange(Z) %>%
    dplyr::mutate(index = row_number()) %>%
    dplyr::mutate(observed = Z,
                  expected = qnorm(index/n_genes))
  return(qq_data)
}
