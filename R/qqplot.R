compute_estimate_expected_rvas <- function(input_data,
                                           model,
                                           mean_intercept_nn=1) {

  no_genes <- nrows(input_data)
  no_samples <- mean(input_data$N)

  # Set initial parameters
  samples <- numeric(0)  # Initialize as an empty numeric vector

  # Sample sizes for each weight
  sample_sizes <- floor(no_samples * model$delta)
  sample_sizes <- as.integer(sample_sizes)  # Ensure integer values

  if (all(sample_sizes == 0)) {
    stop("All sample sizes are zero.")
  }

  # Generate samples in a vectorized way
  samples <- unlist(mapply(function(endpoint, size) {
    if (size > 0) {
      endpoint * runif(size)
    } else {
      numeric(0)
    }
  }, model$component_endpoints, sample_sizes))

  # Add noise
  samples <- samples + rnorm(no_samples, sd = sqrt(mean_intercept_nn))

  # Compute estimated quantiles
  qq_data <- quantile(samples * sqrt(mean_intercept_nn), ((1:no_genes) / no_genes) - 0.5 / no_genes)

  return(qq_data)
}

compute_true_expected_rvas <- function(input_data){
  if(is.null(input_data$annotation)){
    input_data$annotation <- 1
  }
  qq_data <- input_data %>%
    dplyr::group_by(annotation) %>%
    dplyr::mutate(SE = sd(effect_estimate),
           Z = effect_estimate/as.numeric(effect_se)) %>%
    dplyr::group_by(annotation) %>% # TO CHECK
    dplyr::arrange(Z) %>%
    dplyr::mutate(index = row_number(),
                  n_genes = n()) %>%
    dplyr::mutate(observed = Z/sqrt(1) ,
                  expected = qnorm(index/n_genes))
  return(qq_data)
}
