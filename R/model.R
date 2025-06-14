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

#' Generate Grid Points and Component Weight Distributions
#'
#' Creates a fine grid of effect sizes (\eqn{\beta}) spanning intervals defined by
#' component endpoints. It also calculates the weight distribution for each
#' component across this grid. Each component's distribution represents a
#' uniform density between 0 and its corresponding endpoint (or between two
#' endpoints if defining intervals).
#'
#' @param component_endpoints Numeric vector. Sorted endpoints defining the boundaries
#'   of the effect size components. For example, `c(0, 0.1, 0.5)` defines
#'   three components: [0, 0], [0, 0.1], and [0, 0.5].
#' @param grid_points_per_component Positive integer. The number of grid points
#'   to generate *within* each interval defined by adjacent endpoints. The total
#'   number of grid points will be `(length(component_endpoints) - 1) * grid_points_per_component`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{grid}{Numeric vector. The combined grid points across all intervals.}
#'     \item{components}{Matrix (n_components x n_grid_points). Each row represents
#'       a component, and the values are the probability weights for that
#'       component at each grid point.}
#'   }
#' @export
#' @examples
#' endpoints <- c(0, 0.1, 0.5)
#' grid_info <- make_grid(endpoints, grid_points_per_component = 5)
#' print(grid_info$grid)
#' print(grid_info$components)
make_grid <- function(component_endpoints, grid_points_per_component) {
  if (!is.numeric(component_endpoints)) stop("component_endpoints must be numeric.")
  if (any(diff(component_endpoints) < 0)) stop("component_endpoints must be sorted in increasing order.")
  if (min(component_endpoints^2) != 0) stop("component_endpoints must include 0.")
  if (!is.numeric(grid_points_per_component) || length(grid_points_per_component) != 1 ||
      grid_points_per_component < 1 || grid_points_per_component %% 1 != 0) {
    stop("grid_points_per_component must be a positive integer.")
  }
  # Build sequences for each interval
  grid_list <- lapply(seq_len(length(component_endpoints) - 1), function(i) {
    seq(component_endpoints[i], component_endpoints[i + 1], length.out = grid_points_per_component + 1)
  })
  # Combine and drop the last point of each interval to avoid duplication
  grid <- unlist(lapply(grid_list, function(s) s[-length(s)]))
  # Compute weights via adjacent differences
  diffs <- c(diff(grid), 0)
  # Build component matrix: rows=components, cols=grid points
  component_matrix <- t(sapply(seq_along(component_endpoints), function(k) {
    endpoint_val <- component_endpoints[k]
    lb <- min(0, endpoint_val); ub <- max(0, endpoint_val)
    w <- ifelse(grid >= lb & grid <= ub, diffs, 0)
    if (sum(w) > 0) w / sum(w) else w
  }))
  return(list(grid = grid, components = component_matrix))
}

#' Initialize Grid-Based BurdenEM Model
#'
#' Sets up the initial structure for the BurdenEM model using a grid-based
#' approximation for the effect size distribution. It generates the grid,
#' computes gene-specific likelihoods across the grid, prepares the input
#' data frame, and initializes mixture weights (delta).
#'
#' @param gene_data Data frame. Must contain gene-level summary statistics,
#'   including a `burden_score` or `mutation_rate` column
#'   and optionally a `features` column (list column where each element is a
#'   numeric vector of features for that gene). If `features` is missing, a
#'   column of 1s is added.
#' @param likelihood_fn Function. A function that calculates the likelihood of
#'   observing the gene's data given a vector of possible effect sizes (the grid).
#'   It should take two arguments: a single row of `gene_data` and the `grid_effects`
#'   vector, and return a vector of likelihoods (one per grid point). See
#'   `likelihoods.R` for examples.
#' @param component_endpoints Numeric vector. Sorted endpoints defining the
#'   component boundaries, passed to `make_grid`. Should include 0.
#' @param h2_function Function. A function that takes as input an effect size
#'   and a row of gene data, returning the heritability explained
#'   by that gene if it has that effect size.
#' @param grid_points_per_component Positive integer. Number of grid points per
#'   component interval, passed to `make_grid`. Defaults to 10.
#'
#' @return A list object representing the BurdenEM model, containing:
#'   \describe{
#'     \item{df}{Data frame. Modified `gene_data` including the computed likelihood matrix
#'       as a list column `likelihood` (genes x grid points) and `gene_score`.}
#'     \item{grid_effects}{Numeric vector. The grid points generated by `make_grid`.}
#'     \item{components}{Matrix. The component weight distributions from `make_grid`.}
#'     \item{delta}{Matrix (n_features x n_components). Initial mixture weights,
#'       typically uniform.}
#'     \item{h2_function}{Function. The provided heritability calculation function.}
#'     \item{component_endpoints}{Numeric vector. The provided component endpoints.}
#'     \item{null_index}{Integer. The index of the component corresponding to zero effect (endpoint = 0).}
#'     \item{drop_columns}{List of column names to drop from the output data frame.}
#'   }
#' @export
initialize_grid_model <- function(gene_data,
                                  likelihood_fn,
                                  component_endpoints,
                                  h2_function = function(beta, score) score * beta^2,
                                  grid_points_per_component = 10,
                                  drop_columns = c()) {
  # 1. Generate grid and component distributions
  grid_obj <- make_grid(component_endpoints, grid_points_per_component)
  grid_effects <- grid_obj$grid
  components <- grid_obj$components
  n_genes <- nrow(gene_data)

  likelihood <- t(sapply(seq_len(n_genes), function(i) {
    likelihood_fn(gene_data[i, ], grid_effects)
  }))

  # 3. Ensure 'features' column exists
  if (!"features" %in% colnames(gene_data)) {
    gene_data$features <- rep(1, n_genes)
    n_features <- 1
  }
  else
  {
    n_features <- length(gene_data$features[[1]])
  }

  n_components <- nrow(components)
  delta <- matrix(1/n_components, nrow = n_features, ncol = n_components)

  # 4. Build output DF with likelihood
  df <- gene_data %>%
    dplyr::select(-dplyr::any_of(drop_columns))
  df$likelihood <- likelihood

  # Index of null component
  null_index <- which(component_endpoints == 0)[1]

  # 5. Assemble model object
  model <- list(
    df = df,
    grid_effects = grid_effects,
    components = components,
    delta = delta,
    h2_function = h2_function,
    component_endpoints = component_endpoints,
    null_index = null_index
  )
  return(model)
}

initialize_model <- function(likelihood_function,
                             genetic_data,
                             component_endpoints,
                             features,
                             grid_size){

  conditional_likelihood <- likelihood_function(genetic_data,
                                                component_endpoints,
                                                grid_size)

  no_cpts <- length(component_endpoints)

  if (is.null(features)) {
    features <- matrix(1, nrow = nrow(genetic_data), ncol = 1)
    rownames(features) <- rownames(genetic_data)
  }

  delta_init <- matrix(1/no_cpts, nrow = ncol(features), ncol = no_cpts)

  model <- list(component_endpoints = component_endpoints,
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

#' Compute posterior probability over grid effects for each gene
#' @param model BurdenEM grid model as created by initialize_grid_model
#' @return Matrix of size n_genes x n_grid with row-sums = 1
grid_posterior <- function(model){
  Fmat <- do.call(rbind, model$df$features)       # genes x features
  weights_grid <- Fmat %*% model$delta %*% model$components  # genes x grid
  unnorm_post <- weights_grid * model$df$likelihood          # element-wise
  post <- unnorm_post / rowSums(unnorm_post)
  return(post)
}

#' Calculate Complete Data Likelihoods
#'
#' Computes the likelihood conditional upon assignment of each gene to each component.
#'
#' @param model A BurdenEM model object.
#'
#' @return Matrix (n_genes x n_components). Likelihood of each gene under each component.
#' @export
complete_data_likelihoods <- function(model) {
  return (model$df$likelihood %*% t(model$components))
}

#' Calculate Mixture Weights (Prior Probabilities per Component)
#'
#' Computes the prior probability for each gene belonging to each component,
#' based on the gene's features and the current mixture parameters (`delta`).
#' `P(Component_k | Features_i, delta) = Features_i %*% delta[,k]`
#'
#' @param model A BurdenEM model object.
#'
#' @return Matrix (n_genes x n_components). Prior probability (mixture weight)
#'   for each gene and component combination. Rows *do not* necessarily sum to 1.
#' @export
mixture_weights <- function(model) {
  features <- do.call(rbind, model$df$features)
  return (features %*% model$delta)
}

#' Calculate Posterior Assignment Probabilities (Responsibilities)
#'
#' Computes the posterior probability (responsibility) of each gene belonging
#' to each component, given the data.
#'
#' @param model A BurdenEM model object.
#'
#' @return Matrix (n_genes x n_components). Entry (i,j) is the probability
#'   of gene i belonging to component j.
#' @export
responsibilities <- function(model) {
  weights <- mixture_weights(model)
  cdl <- complete_data_likelihoods(model)
  assignments <- cdl * weights
  return (assignments / rowSums(assignments))
}

#' Calculate Conditional Expectation of a Function Given Component Assignment
#'
#' Computes the expected value of `function_to_integrate(beta)` for each gene,
#' conditional on that gene belonging to each component `k`.
#' The function to integrate can optionally depend on gene-specific data.
#'
#' @param model A BurdenEM model object.
#' @param function_to_integrate A function to be integrated.
#'   It can take one or two primary arguments:
#'   1. `f(beta)`: `beta` is a numeric vector of effect sizes (`model$grid_effects`).
#'      The function must return a numeric vector of the same length as `beta`.
#'   2. `f(beta, gene_data_row)`: `beta` is as above. `gene_data_row` is a single
#'      row (a data.frame or list) from `model$df` corresponding to the current gene.
#'      The function must return a numeric vector of the same length as `beta`.
#'   In both cases, the function is expected to be vectorized with respect to its first argument (`beta`).
#'   Defaults to `function(beta) beta`.
#'
#' @return Matrix (n_genes x n_components). Entry (i,j) is the expected value of
#'   `function_to_integrate(beta_values, model$df[i,])` (or just `function_to_integrate(beta_values)`
#'   if it's a single-argument function) if gene i belongs to component j.
#' @export
complete_data_means <- function(model,  function_to_integrate = function(beta) beta) {

  num_formal_args <- length(formals(function_to_integrate))

  denominator <- model$df$likelihood %*% t(model$components) # (genes x components)

  if (num_formal_args > 1) {
    # function_to_integrate depends on gene_data_row
    F_mat_list <- lapply(1:nrow(model$df), function(gene_idx) {
      function_to_integrate(model$grid_effects, model$df[gene_idx, , drop = FALSE])
    })
    F_mat <- do.call(rbind, F_mat_list)

    # Dimensionality check for F_mat
    if (!is.matrix(F_mat) || nrow(F_mat) != nrow(model$df) || ncol(F_mat) != length(model$grid_effects)) {
        stop(paste("The result of applying function_to_integrate(beta_vector, gene_data_row) per gene",
                   "does not produce a matrix of expected dimensions (n_genes x n_grid_points).",
                   "Please ensure it returns a vector of length n_grid_points."))
    }
    numerator <- (model$df$likelihood * F_mat) %*% t(model$components)
  } else {
    # function_to_integrate only depends on beta (grid_effects)
    fvals_vec <- function_to_integrate(model$grid_effects) # (grid_vector)
    numerator <- (model$df$likelihood %*% diag(fvals_vec)) %*% t(model$components)
  }

  result <- numerator / denominator
  result[denominator == 0] <- 0

  return(result)
}

#' Posterior expectation of a function of beta for each gene
#'
#' @param model Fitted burdenEM model.
#' @param function_to_integrate Function to integrate over the effect size distribution, which can have
#'   one or two arguments: an effect size vector, and optionally a row of the data frame `df`.
#' @return Numeric vector of length = number of genes.
posterior_expectation2 <- function(model, function_to_integrate = function(x) x) {
    weights <- mixture_weights(model)
    cd_means <- complete_data_means(model, function_to_integrate)
    cd_likelihoods <- complete_data_likelihoods(model)
    return(rowSums(cd_means * cd_likelihoods * weights) / rowSums(cd_likelihoods * weights))
}

#' Compute mean and standard error of a posterior mean summed across genes
#'
#' This function calculates the posterior expectation of a scalar function
#' of the model effect size, summed across all genes,
#' and provides the standard error of this quantity.
#'
#' @param model Fitted burdenEM model
#' @param function_to_integrate Function to integrate over the effect size distribution, which can have
#'   one or two arguments: a single row of `gene_data` and the `grid_effects` vector.
#' @return A list containing two elements:
#'   \item{mean}{The mean of the integrated function summed across genes within each stratum}
#'   \item{se}{The standard error of the integral within each stratum}
#'
#' @export
posterior_expectation_with_se <- function(model,
                                          function_to_integrate = function(x) x) {
    weights <- mixture_weights(model)
    cd_means <- complete_data_means(model, function_to_integrate)

    cd_likelihoods <- complete_data_likelihoods(model)
    likelihoods <- rowSums(cd_likelihoods * weights)

    means <- rowSums(cd_means * cd_likelihoods * weights) / likelihoods
    # entry i,k: d/dw_k E(f(beta_i)) = d/dw_k E(cd_means[i,z_i])
    derivatives <- (cd_likelihoods / likelihoods) * (cd_means - means)

    num_strata <- length(model$information)
    features <- do.call(rbind, model$df$features)
    mean <- rep(0, num_strata)
    se <- rep(0, num_strata)
    for (i in 1:num_strata) {
        mean[i] <- t(features[,i]) %*% means
        derivatives_i <- t(features[,i]) %*% derivatives
        derivatives_i <- derivatives_i[-model$null_index] - derivatives_i[model$null_index]

        information <- model$information[[i]]
        offset <- 1e-6 * mean(diag(information)) * diag(nrow(information))
        information <- information + offset
        se[i] <- sqrt(derivatives_i %*% solve(information, derivatives_i))
}

return(list(mean = mean,
            se = se))
}

posterior_gene_samples <- function(model,
                                   num_samples = 1) {

  weights <- model$features %*% model$delta
  posteriors = weights * model$conditional_likelihood

  posteriors <- posteriors / rowSums(posteriors)

  gene_samples = sapply(1:nrow(posteriors),
                        function(x) {
                          endpoint_samples = sample(model$component_endpoints, size = num_samples, prob =posteriors[x,], replace = TRUE)
                          loggamma_samples = runif(num_samples,min = pmin(0,endpoint_samples), max = pmax(0,endpoint_samples))
                          return(loggamma_samples)
                        })

  return(t(gene_samples))

}
