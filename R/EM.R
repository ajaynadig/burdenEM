# Core Expectation-Maximization scripts
EM_fit <- function(model,
                   max_iter,
                   tol = NULL,
                   return_likelihood = TRUE) {
  #Pre-invert X_T %*% X to save time
  OLS_denom = solve(t(model$features) %*% model$features)
  OLS_denom_t_features <- OLS_denom %*% t(model$features)
  ll = rep(NA, max_iter)
  ll_change = tol + 1
  #first iteration
  weights <- model$features %*% model$delta
  posteriors <- weights * model$conditional_likelihood
  ll[1] = sum(log(rowSums(posteriors)))
  posteriors <- posteriors / rowSums(posteriors)
  model$delta <- OLS_denom_t_features %*% posteriors

  iter_count = 1


  while (TRUE){#(ll_change > tol) {

    if (iter_count >= max_iter) {break}

    weights <- model$features %*% model$delta
    posteriors <- weights * model$conditional_likelihood

    ll[iter_count + 1] = sum(log(rowSums(posteriors)))
    ll_change = abs(ll[iter_count+1]-ll[iter_count])/abs(ll[iter_count])


    posteriors <- posteriors / rowSums(posteriors)
    model$delta <- OLS_denom_t_features %*% posteriors

    iter_count = iter_count + 1
  }

  if (return_likelihood) {
    model$ll = ll
  }

  return(model)
}

bootstrap_EM <- function(model,
                         n_boot,
                         max_iter,
                         bootstrap_samples = NULL,
                         bootstrap_seeds = NULL) {
  if (is.null(bootstrap_samples)) {
    if (is.null(bootstrap_seeds)) {
      bootstrap_seeds <- 1:n_boot
    }
    bootstrap_samples <- sapply(bootstrap_seeds,
                                function(dummy) {
                                  set.seed(dummy)
                                  sample(1:nrow(model$conditional_likelihood),replace = TRUE)
                                })
    cat("...bootstrap EM")

  } else {
    cat("...bootstrap EM with user-specified samples")

  }

  bootstrap_delta <- lapply(1:n_boot,
                            function(iter) {

                              model_boot = model
                              model_boot$conditional_likelihood = model_boot$conditional_likelihood[bootstrap_samples[,iter],]
                              model_boot$features = model_boot$features[bootstrap_samples[,iter],]


                              boot_output <- EM_fit(model_boot,
                                                    max_iter)

                              if (iter %% 20 == 0) {
                                cat(paste0("...",iter,"...("))
                                cat(paste0(sum(!is.na(boot_output$ll))," iters)..."))
                              }

                              return(boot_output$delta)
                            })

  return(list(bootstrap_delta = bootstrap_delta,
              bootstrap_samples = bootstrap_samples))
}

null_EM_trio <- function(genetic_data,
                         model,
                         max_iter,
                         n_null,
                         grid_size) {

  cat("...null EM")

  null_coefs <- lapply(1:n_null,
                       function(dummy) {
                         if (dummy %% 20 == 0) {
                           cat(paste0("...",dummy))
                         }

                         genetic_data_null = genetic_data

                         genetic_data_null$case_count = rpois(nrow(genetic_data),
                                                              genetic_data$expected_count)

                         model_null = initialize_model(likelihood_function = poisson_uniform_likelihood,
                                                       genetic_data = genetic_data_null,
                                                       component_endpoints = model$component_endpoints,
                                                       features = model$features,
                                                       grid_size = grid_size)


                         model_null = EM_fit(model_null,
                                             max_iter)

                         return(model_null$delta)

                       })


  return(null_coefs)
}

null_EM_rvas <- function(genetic_data,
                         model,
                         max_iter,
                         n_null,
                         grid_size) {

  cat("...null EM")

  null_coefs <- lapply(1:n_null,
                       function(dummy) {
                         if (dummy %% 20 == 0) {
                           cat(paste0("...",dummy))
                         }

                         genetic_data_null = genetic_data

                         genetic_data_null$effect_estimate  = rnorm(nrow(genetic_data), mean = genetic_data$effect_estimate,
                                                              sd=genetic_data$effect_se/sqrt(genetic_data$N))

                         model_null = initialize_model(likelihood_function = normal_uniform_likelihood,
                                                       genetic_data = genetic_data_null,
                                                       component_endpoints = model$component_endpoints,
                                                       features = model$features,
                                                       grid_size = grid_size)


                         model_null = EM_fit(model_null,
                                             max_iter = max_iter)

                         return(model_null$delta)
                       })


  return(null_coefs)
}
