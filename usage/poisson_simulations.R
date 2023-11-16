nn <- 10000
normal_model <- FALSE

poisson_simulation <- function(n_gene, n_covariate, n_mixcomponent) {

}




# Generating features matrix
features <- matrix(runif(nn), nrow = nn, ncol = 2)
features[, 2] <- 1 - features[, 1]

model_params <- c(-4, 1, 0, 1, 4)
coefs <- matrix(c(.01, 0, .04, 0, .9, 1, .04, 0, .01, 0), nrow = 2, byrow = TRUE)

# Mixture component assignments
weights <- features %*% coefs
weights <- weights / replicate(ncol(weights),rowSums(weights))

cpt_param <- numeric(nn)
for (ii in 1:nn) {
  cpt_param[ii] <- sample(model_params, 1, prob = weights[ii, ], replace = TRUE)
}

# Effect sizes
beta <- runif(nn)
beta <- beta * cpt_param

# Noise in the effect estimates
if (normal_model) {
  noise_var <- runif(nn)
  beta_hat <- beta + sqrt(noise_var) * rnorm(nn)

  # Estimation
  system.time({
    coefs_est <- burdenEM(model_params, effect_estimate = beta_hat, effect_se = sqrt(noise_var),
                          features = features, model_type = 'uniform')
  })
} else {
  case_rate <- 10 * runif(nn)
  lambda <- case_rate * exp(beta)
  case_count <- rpois(nn, lambda)

  # Estimation
  timediff = system.time({
    coefs_est <- burdenEM(model_params, case_count = case_count, case_rate = case_rate,
                          features = features, model_type = 'uniform')
  })
}

print(coefs_est)
print(timediff)
