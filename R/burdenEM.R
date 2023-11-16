burdenEM <- function(mixture_params, model_type = 'uniform', effect_estimate = NULL,
                     effect_se = NULL, case_count = NULL, case_rate = NULL, p_value = NULL,
                     features = NULL, num_iter = 100) {

  # Validate inputs
  stopifnot(is.vector(mixture_params))

  no_cpts <- length(mixture_params)

  # Handling effect estimates and standard errors
  if (!is.null(effect_estimate)) {
    stopifnot(is.null(case_count))  # Ensure that case_count is not provided
    if (is.null(effect_se)) {
      stopifnot(!is.null(p_value))  # p_value must be specified if effect_se is not
      z_score <- qnorm(1 - p_value / 2) * sign(effect_estimate)
      effect_se <- abs(effect_estimate / z_score)
    }
    stopifnot(length(effect_se) == length(effect_estimate))
    data_type <- 'normal'
    no_tests <- length(effect_estimate)
  } else {
    stopifnot(!is.null(case_count))
    stopifnot(model_type == 'uniform')  # Only uniform is supported for Poisson
    if (is.null(case_rate)) {
      stopifnot(!is.null(p_value))
      # Implement poisson_rate function based on the MATLAB code
      case_rate <- poisson_rate(case_count, p_value)
    }
    stopifnot(length(case_rate) == length(case_count))
    data_type <- 'poisson'
    no_tests <- length(case_count)
  }

  # Compute likelihood
  likelihood <- matrix(0, nrow = no_tests, ncol = no_cpts)
  mu_grid <- seq(0.05, 1, by = 0.1) #THE FINENESS OF THIS GRID SHOULD PROBABLY BE A USER SPECIFIED VALUE
  if (model_type == 'uniform' && data_type == 'normal') { #CODE FOR NORMAL DATA LIKELIHOOD PROBABLY BROKEN
    for (kk in 1:no_cpts) {
      mu <- mu_grid * mixture_params[kk]
      sigma <- effect_se
      likelihood[, kk] <- rowMeans(dnorm(effect_estimate, mean = mu, sd = sigma))
    } #SHOULD THERE BE LIKELIHOOD CALCULATION FOR NORMAL DATA, NORMAL MODEL HERE?
  } else if (model_type == 'uniform' && data_type == 'poisson') {
    for (kk in 1:no_cpts) {
      # Expand case_count into a matrix by replicating the vector along the columns
      case_count_matrix <- replicate(length(mu_grid),case_count)

      # Calculate the rate for each mu value in the grid and each case count
      rate <- case_rate * t(replicate(length(case_rate),exp(mu_grid * mixture_params[kk])))

      # Compute the likelihood for each case count and rate combination
      likelihoods <- dpois(case_count_matrix, rate)

      # Calculate the average likelihood for each test across the grid of mu values
      likelihood[, kk] <- rowMeans(likelihoods)
    }
  } else {
    stop('Options for model_type are uniform or normal')
  }

  if (is.null(features)) {
    features <- matrix(1, nrow = no_tests, ncol = 1)
  }

  stopifnot(all(rowSums(features) == 1) & all(features >= 0))


  # EM algorithm
  coefs <- matrix(1, nrow = ncol(features), ncol = no_cpts)
  #Pre-invert X_T %*% X to save time
  OLS_denom = solve(t(features) %*% features)
  for (rep in 1:num_iter) {
    weights <- features %*% coefs
    posteriors <- weights * likelihood
    posteriors <- posteriors / replicate(ncol(posteriors),rowSums(posteriors))
    coefs <- OLS_denom %*% t(features) %*% posteriors
    #print(coefs)
  }

  return(list(coefs = coefs, mixture_params = mixture_params))
  #Other stuff that we probably should output...
  # Some sort of convergence indicator i.e. likelihood
  # Some QC metrics to help check that inferred distribution is plausible

}

poisson_rate <- function(x, p, steps = 100) {
  lambda <- x
  factx <- factorial(x)

  for (step in 1:steps) {
    cdf <- ppois(x, lambda, lower.tail = FALSE)
    lambda <- lambda - (cdf - (1 - p)) / (exp(-lambda) * lambda^x / factx)
  }

  return(lambda)
}

