#Functions to manipulate burdenEM models

choose_component_endpoints_trio = function(component_endpoints,
                                 no_cpts,
                                 prevalence) {
  if (!is.null(component_endpoints)) {
    return(component_endpoints)
  } else {
    component_endpoints = seq(0,log(1/prevalence),length.out = no_cpts)
  }
}

choose_component_endpoints_rvas <- function(component_endpoints,
                                           no_cpts,
                                           input_data) {
  if (is.null(component_endpoints)) {
    lower_bound_5 = unlist(quantile(input_data$effect_estimate, probs = c(0.05)))
    upper_bound_95 = unlist(quantile(input_data$effect_estimate, probs = c(0.95)))
    component_endpoints = c(min(input_data$effect_estimate),
                            seq(lower_bound_5, upper_bound_95, length.out = (no_cpts-2)),
                            max(input_data$effect_estimate))
  }
  return(component_endpoints)
}


initialize_model <- function(likelihood_function,
                             genetic_data,
                             component_endpoints,
                             features,
                             grid_size){

  conditional_likelihood = likelihood_function(genetic_data=genetic_data,
                                               component_endpoints=component_endpoints,
                                               grid_size=grid_size)

  no_cpts = length(component_endpoints)

  if (is.null(features)) {
    features <- matrix(1, nrow = nrow(genetic_data), ncol = 1)
    rownames(features) = rownames(genetic_data)
  }

  delta_init = matrix(1, nrow = ncol(features), ncol = no_cpts)

  model = list(component_endpoints = component_endpoints,
               delta = delta_init,
               conditional_likelihood = conditional_likelihood,
               features = features,
               grid_size = grid_size)
}

posterior_expectation <- function(model,
                                  genetic_data,
                                  function_to_integrate,
                                  grid_size) {
  weights <- model$features %*% model$delta
  posteriors = weights * model$conditional_likelihood

  posteriors <- posteriors / rowSums(posteriors)

  no_tests = nrow(model$conditional_likelihood)
  no_cpts = length(model$component_endpoints)

  conditional_posterior_expectations <- matrix(NA, nrow = no_tests, ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)

  for (kk in 1:no_cpts) {
    #Expand case_count into a matrix by replicating the vector along the columns
    case_count_matrix <- replicate(grid_size,genetic_data$case_count)

    # Calculate the rate for each mu value in the grid and each case count
    rate <- genetic_data$expected_count * t(replicate(no_tests,exp(mu_grid * model$component_endpoints[kk])))

    # Compute the likelihood for each case count and rate combination: likelihood within each components
    likelihoods <- dpois(case_count_matrix, rate)

    function_vals = t(replicate(no_tests,
                                  function_to_integrate(mu_grid * model$component_endpoints[kk])))

    conditional_posterior_expectations[,kk] <- rowMeans(function_vals * likelihoods)/rowMeans(likelihoods)
  }

  posterior_expectations = rowSums(posteriors * conditional_posterior_expectations)
}


posterior_expectation_rvas <- function(model,
                                       genetic_data,
                                       function_to_integrate,
                                       grid_size) {
  weights <- model$features %*% model$delta
  posteriors = weights * model$conditional_likelihood
  posteriors_rowsums <- rowSums(posteriors)
  posteriors_rowsums <- if_else(posteriors_rowsums == 0, 1e-300, posteriors_rowsums)

  posteriors <- posteriors / posteriors_rowsums

  no_tests = nrow(model$conditional_likelihood)
  no_cpts = length(model$component_endpoints)

  conditional_posterior_expectations <- matrix(NA, nrow = no_tests, ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)

  trait_type = unique(genetic_data$trait_type)
  for(gg in 1:no_tests){
    for (kk in 1:no_cpts) {
      tmp_gene_data <- genetic_data[gg,]
      tmp_cond_lik <- likelihood_per_combo_rvas()
    }
  }

  for (kk in 1:no_cpts) {
    #Expand case_count into a matrix by replicating the vector along the columns
    case_count_matrix <- replicate(grid_size,genetic_data$case_count)

    # Calculate the rate for each mu value in the grid and each case count
    rate <- genetic_data$expected_count * t(replicate(no_tests,exp(mu_grid * model$component_endpoints[kk])))

    # Compute the likelihood for each case count and rate combination: likelihood within each components
    likelihoods <- dpois(case_count_matrix, rate)

    function_vals = t(replicate(no_tests,
                                function_to_integrate(mu_grid * model$component_endpoints[kk])))

    conditional_posterior_expectations[,kk] <- rowMeans(function_vals * likelihoods)/rowMeans(likelihoods)
  }

  posterior_expectations = rowSums(posteriors * conditional_posterior_expectations)
}
