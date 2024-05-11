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

initialize_model <- function(likelihood_function, genetic_data, component_endpoints,features, grid_size){

  conditional_likelihood = likelihood_function(genetic_data, component_endpoints,grid_size)

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
                                  function_to_integrate) {
  weights <- model$features %*% model$delta
  posteriors = weights * model$conditional_likelihood

  posteriors <- posteriors / rowSums(posteriors)

  no_tests = nrow(model$conditional_likelihood)
  no_cpts = length(model$component_endpoints)

  conditional_posterior_expectations <- matrix(NA, nrow = no_tests, ncol = no_cpts)


  for (kk in 1:no_cpts) {

  function_vals = t(replicate(no_tests,
                                  function_to_integrate(mu_grid * mixture_params[kk])))
  conditional_posterior_expectations[,kk] <- rowMeans(function_vals * model$conditional_likelihood)/rowMeans(model$conditional_likelihood)
  }

  posterior_expectations = rowSums(posteriors * conditional_posterior_expectations)
}
