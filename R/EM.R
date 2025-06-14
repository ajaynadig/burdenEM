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

# --- New function: EM for grid-based models ---
EM_fit_grid <- function(model, max_iter = 1000) {
  print(paste("Running EM fit for grid-based model with max_iter = ", max_iter))
  features <- do.call(rbind, model$df$features)
  cdl <- model$df$likelihood %*% t(model$components)
  XtX_inv <- solve(t(features) %*% features)
  delta <- model$delta
  ll_old <- -Inf
  ll_new <- NA
  for (i in seq_len(max_iter)) {
    weights <- features %*% delta
    post <- weights * cdl
    ll_new <- sum(log(rowSums(post)))
    # if (!is.na(ll_old) && abs((ll_new - ll_old) / ll_new) < tol) break
    post <- post / rowSums(post)
    delta <- XtX_inv %*% t(features) %*% post
    ll_old <- ll_new
  }
  model$delta <- delta
  model$ll <- ll_new
  return(model)
}

bootstrap_EM <- function(model,
                         n_boot,
                         max_iter,
                         bootstrap_samples = NULL,
                         bootstrap_seeds = NULL,
                         tol = 0) {
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
                              model_boot$conditional_likelihood = model_boot$conditional_likelihood[bootstrap_samples[,iter], , drop = FALSE]
                              model_boot$features = model_boot$features[bootstrap_samples[,iter], , drop = FALSE]

                              boot_output <- EM_fit(model = model_boot, max_iter = max_iter, tol)

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
                         grid_size,
                         tol) {

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


                         model_null = EM_fit(model_null,max_iter = max_iter, tol)

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

#' Calculate Information Matrix for Delta Parameters
#'
#' Computes the observed Fisher information matrix for each feature stratum
#' as the negative Hessian of the log-likelihood with respect to the delta parameters.
#' The inverse of this matrix provides an estimate of the covariance matrix for delta.
#'
#' The delta parameters represent the mixture weights for components within each
#' feature stratum. This information matrix (and its inverse, the covariance) is
#' used for calculating standard errors of heritability estimates. The calculation
#' is based on the method described by Louis (1982) for observed information
#' from incomplete data.
#'
#' @param model A BurdenEM model object, which must have been at least initialized
#'   (e.g., via `initialize_grid_model`) and ideally fitted (though fitting primarily
#'   updates `model$delta`, which is an input here).
#'   Key model components used: `model$null_index`, `model$df$features`,
#'   `model$df$likelihood`, `model$components`, `model$delta`.
#'
#' @return A list of matrices. Each matrix in the list corresponds to a feature
#'   stratum and represents the information matrix for the delta parameters in that
#'   stratum (excluding the null component, which is used as a reference).
#'   The dimensions of each matrix are (n_components - 1) x (n_components - 1).
#' @export
information_matrices <- function(model) {
  remove_index <- model$null_index
  features <- do.call(rbind, model$df$features)
  cdl <- model$df$likelihood %*% t(model$components)
  weights <- features %*% model$delta
  gene_likelihood <- rowSums(weights * cdl)
  cdl_normalized <- cdl / gene_likelihood
  covariance <- list()
  for (i in 1:ncol(features)) {
    feature_cdl <- features[,i] * cdl_normalized
    feature_cdl <- feature_cdl[,-remove_index] - feature_cdl[,remove_index]
    covariance[[i]] <- t(feature_cdl) %*% feature_cdl
  }
  return(covariance)
}
