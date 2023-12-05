packages = c('dplyr', 'ggplot2')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
}

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
  if (model_type == 'normal'){
    stopifnot(data_type == 'normal')
    likelihood=sapply(model_params, function(x){dnorm(effect_estimate, mean=0, sd=sqrt(x + effect_se^2))})
    # print(head(likelihood))
  }else if(model_type == 'uniform'){
    if (data_type == 'normal') { #CODE FOR NORMAL DATA LIKELIHOOD PROBABLY BROKEN
      for (kk in 1:no_cpts) {
        mu <- mu_grid * mixture_params[kk]
        sigma <- effect_se
        likelihood[, kk] <- rowMeans(dnorm(effect_estimate, mean = mu, sd = sigma))
      } #SHOULD THERE BE LIKELIHOOD CALCULATION FOR NORMAL DATA, NORMAL MODEL HERE?
    } else if (data_type == 'poisson') {
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

simulation_results <- function(model_params, n_gene, n_features=2, data_type = 'poisson') {
  features <- matrix(runif(n_gene), nrow = n_gene, ncol = n_features)
  features[, 2] <- 1 - features[, 1] #TODO: generalize for more than 2 features (?)
  
  coefs = rbind(dirmult::rdirichlet(n = 1, rep(1, length(model_params))),
                dirmult::rdirichlet(n = 1, rep(1, length(model_params))))
  
  weights <- features %*% coefs
  weights <- weights / replicate(ncol(weights),rowSums(weights))
  
  cpt_param <- numeric(n_gene)
  for (ii in 1:n_gene) {
    cpt_param[ii] <- sample(model_params, 1, prob = weights[ii, ], replace = TRUE)
  }
  
  if(data_type == 'poisson'){
    # Effect sizes
    beta <- runif(n_gene)
    beta <- beta * cpt_param
    
    case_rate <- 10 * runif(n_gene)
    lambda <- case_rate * exp(beta)
    case_count <- rpois(n_gene, lambda)
    # Estimation
    coefs_est <- burdenEM(model_params, 
                          case_count = case_count, 
                          case_rate = case_rate,
                          features = features, 
                          model_type = 'uniform')
    
  }else if(data_type == 'normal'){
    beta <- rnorm(n_gene, 0, sqrt(cpt_param))
    # se_hat <- matrix(runif(nG*nP, min=0, max=0.1), nrow = nP, ncol = nG, byrow = TRUE) 
    se_hat <- matrix(rep(1/nG, nG), nrow = 1, ncol = nG) 
    beta_hat <- beta + se_hat*rnorm(nG, 0, 1)
    system.time({
      coefs_est <- burdenEM(model_params, 
                            model_type = 'normal', 
                            effect_estimate = beta_hat, 
                            effect_se = se_hat,
                            features = features)
    })
  }else{
    stop("'data_type' must be from ('poisson', 'normal')")
  }
  
  
  
  
  return(list(true_coef = coefs,
              est_coef = coefs_est$coefs,
              mixture_params = coefs_est$mixture_params))
  
}


simulation_plots <- function(model_params, data_type, n_gene=10000, n_features=2, rounds=100, save = FALSE, outdir = '~/Desktop/', name='simulation_plot'){
  set.seed(1000)
  simulation_output <- lapply(1:rounds,
                              function(x) {
                                if(x%%10 ==0) print(x)
                                output = simulation_results(n_gene, model_params = model_params, data_type=data_type)
                              })
  
  results <- data.frame()
  for(i in 1:n_features){
    for(j in 1:length(model_params)){
      tmp_results <- data.frame(t(sapply(1:rounds,
                                         function(x) {
                                           return(c(simulation_output[[x]]$true_coef[i,j],
                                                    simulation_output[[x]]$est_coef[i,j]))
                                         }))) %>%
        mutate(idx = j + (i-1)*length(model_params))
      results <- rbind(results, tmp_results)
    }
  }
  
  library(ggplot2)
  library(dplyr)
  p <- results %>%
    ggplot + aes(x = X1, y = X2, color=factor(idx)) +
    geom_abline(linetype = "dashed")+
    geom_point()+
    geom_smooth(method = "lm", lty = 2, se=F)+
    labs(x = "True Coefficient Value", y = "Estimated Coefficient Value", color = 'Coef index') +
    scale_color_brewer(palette = 'Set3') +
    theme_classic()
  if(save){
    png(paste0(outdir, '/', name, '_', data_type,'.png'), height = 3, width = 5, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(p)
}

# Usage: 
n_genes <- 10000
n_features <- 2 
rounds <- 100
save <- TRUE
model_params <- c(0.100, 1.325, 2.550, 3.775, 5.000)
simulation_results(model_params = model_params, data_type = 'normal', n_genes=n_genes, n_features=n_features, rounds=rounds, save=save)
simulation_results(model_params = model_params, data_type = 'poisson', n_genes=n_genes, n_features=n_features, rounds=rounds, save=save)