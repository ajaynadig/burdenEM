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
    lower_bound = unlist(quantile(input_data$effect_estimate, probs = c(0.001)))
    upper_bound = unlist(quantile(input_data$effect_estimate, probs = c(0.999)))
    bound = max(abs(min(input_data$effect_estimate)), abs(max(input_data$effect_estimate)))
    component_endpoints = c(-bound,
                            seq(lower_bound, upper_bound, length.out = (no_cpts-2)),
                            bound)
  }
  return(component_endpoints)
}


initialize_model <- function(likelihood_function,
                             genetic_data,
                             component_endpoints,
                             features,
                             grid_size){

  conditional_likelihood = likelihood_function(genetic_data,
                                               component_endpoints,
                                               grid_size)

  no_cpts = length(component_endpoints)

  if (is.null(features)) {
    features <- matrix(1, nrow = nrow(genetic_data), ncol = 1)
    rownames(features) = rownames(genetic_data)
  }

  delta_init = matrix(1/no_cpts, nrow = ncol(features), ncol = no_cpts)

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
    case_count_matrix <- replicate(length(mu_grid),genetic_data$case_count)

    # Calculate the rate for each mu value in the grid and each case count
    rate <- genetic_data$expected_count * t(replicate(no_tests,exp(mu_grid * model$component_endpoints[kk])))

    # Compute the likelihood for each case count and rate combination
    likelihoods <- dpois(case_count_matrix, rate)

    function_vals = t(replicate(no_tests,
                                function_to_integrate(mu_grid * model$component_endpoints[kk])))

    conditional_posterior_expectations[,kk] <- rowMeans(function_vals * likelihoods)/rowMeans(likelihoods)
  }

  posterior_expectations = rowSums(posteriors * conditional_posterior_expectations)
}


posterior_expectation_rvas <- function(model,
                                  genetic_data,
                                  function_to_integrate=function(x){return(x)},
                                  grid_size=100) {
  weights <- model$features %*% model$delta
  posteriors = weights * model$conditional_likelihood

  posteriors <- posteriors / rowSums(posteriors)

  no_tests = nrow(model$conditional_likelihood)
  no_cpts = length(model$component_endpoints)

  conditional_posterior_expectations <- matrix(NA, nrow = no_tests, ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)

  for (kk in 1:no_cpts) {
    mu <- mu_grid * component_endpoints[kk]
    likelihoods <- sapply(mu, function(x){dnorm(x = genetic_data$effect_estimate - x,
                                                mean = 0,
                                                sd = genetic_data$effect_se)})
    function_vals = t(replicate(no_tests,
                                function_to_integrate(mu_grid * model$component_endpoints[kk])))
    conditional_posterior_expectations[,kk] <- rowMeans(function_vals * likelihoods)/rowMeans(likelihoods)
  }
  posterior_expectations = rowSums(posteriors * conditional_posterior_expectations)
  return(posterior_expectations)
}

effective_penetrance_func <- function(model,
                                      genetic_data,
                                      prevalence) {
  peneff_numerator =  mean(posterior_expectation(model,
                                                 genetic_data,
                                                 function(x) {
                                                   (exp(x)-1)*exp(x)
                                                 },
                                                 grid_size = 10))

  peneff_denominator = mean(posterior_expectation(model,
                                                  genetic_data,
                                                  function(x) {
                                                    (exp(x)-1)
                                                  },
                                                  grid_size = 10))

  peneff = prevalence * (peneff_numerator/peneff_denominator)
  return(peneff)
}
