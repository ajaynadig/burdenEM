# Core Expectation-Maximization scripts
EM_fit <- function(model,
                   num_iter) {
  #Pre-invert X_T %*% X to save time
  OLS_denom = solve(t(model$features) %*% model$features)
  OLS_denom_t_features <- OLS_denom %*% t(model$features)

  for (iter in 1:num_iter) {
    weights <- model$features %*% model$delta
    posteriors <- weights * model$conditional_likelihood
    posteriors <- posteriors / rowSums(posteriors)
    model$delta <- OLS_denom_t_features %*% posteriors
  }

  return(model)
}

bootstrap_EM <- function(model,
                         n_boot,
                         num_iter) {
  cat("...bootstrap EM")

  bootstrap_samples <- sapply(1:n_boot,
                              function(dummy) {
                                sample(1:nrow(model$conditional_likelihood),replace = TRUE)
                              })

  bootstrap_delta <- lapply(1:n_boot,
                                function(iter) {
                                  if (iter %% 20 == 0) {
                                    cat(paste0("...",iter))
                                  }

                                  model_boot = model
                                  model_boot$conditional_likelihood = model_boot$conditional_likelihood[bootstrap_samples[,iter],]
                                  model_boot$features = model_boot$features[bootstrap_samples[,iter],]
                                  model_boot$delta = matrix(1, nrow = nrow(model_boot$delta), ncol = ncol(model_boot$delta))


                                  boot_output <- EM_fit(model_boot,
                                                        num_iter)

                                  return(boot_output$coefs)
                                })

  return(bootstrap_delta)
}

null_EM_trio <- function(genetic_data,
                         model,
                         num_iter,
                         n_null) {

  cat("...null EM")

  null_coefs <- sapply(1:n_null,
                       function(dummy) {
                         if (dummy %% 20 == 0) {
                           cat(paste0("...",dummy))
                         }

                         genetic_data_null = genetic_data

                         genetic_data_null$case_count = rpois(nrow(genetic_data),
                                                              genetic_data$expected_count)

                         model_null = initialize_model(likelihood_function = poisson_uniform_likelihood,
                                                       genetic_data = genetic_data_null,
                                                       component_endpoints = mixture_params,
                                                       features = features,
                                                       grid_size = grid_size)


                         model_null = EM_fit(model_null,
                                             num_iter)

                         return(model_null$delta)

                       })


  return(null_coefs)
}
