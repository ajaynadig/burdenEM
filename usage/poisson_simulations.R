
#Building out liability threshold model unit conversion
#Adding N + mutation rate to burdenEM script

#Model misspecification (\delta)
#Mutation rate misspecification
#Varying N
#Simulate gammas using real expected

nn <- 10000
normal_model <- FALSE

poisson_simulation <- function(n_gene, model_params) {
  features <- matrix(runif(nn), nrow = nn, ncol = 2)
  features[, 2] <- 1 - features[, 1]

  coefs = rbind(dirmult::rdirichlet(n = 1, rep(1, length(model_params))),
                dirmult::rdirichlet(n = 1, rep(1, length(model_params))))

  weights <- features %*% coefs
  weights <- weights / replicate(ncol(weights),rowSums(weights))

  cpt_param <- numeric(n_gene)
  for (ii in 1:n_gene) {
    cpt_param[ii] <- sample(model_params, 1, prob = weights[ii, ], replace = TRUE)
  }

  # Effect sizes
  beta <- runif(n_gene)
  beta <- beta * cpt_param

  case_rate <- 10 * runif(nn)
  lambda <- case_rate * exp(beta)
  case_count <- rpois(nn, lambda)

  # Estimation
  coefs_est <- burdenEM(model_params, case_count = case_count, case_rate = case_rate,
                          features = features, model_type = 'uniform')

  return(list(true_coef = coefs,
              est_coef = coefs_est$coefs,
              mixture_params = coefs_est$mixture_params))

}



simulation_output <- lapply(1:100,
                            function(x) {
                              print(x)
                              output = poisson_simulation(10000,model_params = c(-4, -1, 0, 1, 4))
                            })

onecoef_results <- sapply(1:100,
                          function(x) {
                            return(c(simulation_output[[x]]$true_coef[1,1],
                                     simulation_output[[x]]$est_coef[1,1]))
                          })

ggplot(mapping = aes(x = onecoef_results[1,],
                     y = onecoef_results[2,]))+
  geom_abline(linetype = "dashed")+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "True Coefficient Value", y = "Estimated Coefficient Value")

