#Functions to compute likelihoods for burdenEM inference


poisson_uniform_likelihood <- function(genetic_data,
                                       component_endpoints,
                                       grid_size) {
  no_cpts = length(component_endpoints)
  no_genes = nrow(genetic_data)

  likelihood <- matrix(NA,
                       nrow = no_genes,
                       ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)

  for (kk in 1:no_cpts) {
    # Expand case_count into a matrix by replicating the vector along the columns
    case_count_matrix <- replicate(grid_size,genetic_data$case_count)

    # Calculate the rate for each mu value in the grid and each case count
    rate <- genetic_data$expected_count * t(replicate(no_genes,exp(mu_grid * component_endpoints[kk])))

    # Compute the likelihood for each case count and rate combination
    likelihoods <- dpois(case_count_matrix, rate)

    # Calculate the average likelihood for each test across the grid of mu values
    likelihood[, kk] <- rowMeans(likelihoods)
  }
  return(likelihood)
}

normal_uniform_likelihood <- function(genetic_data,
                                      component_endpoints,
                                      grid_size){

  no_cpts <- length(component_endpoints)
  no_genes <- nrow(genetic_data)

  likelihood <- matrix(NA,
                       nrow = no_genes,
                       ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)

  for (kk in 1:no_cpts) {
      mu <- mu_grid * component_endpoints[kk]
      likelihood[,kk] <- rowMeans(
        sapply(mu, function(x){dnorm(x = genetic_data$effect_estimate - x, mean = 0, sd = genetic_data$effect_se)}) # matrix of n_gene * n_mu
      )
  }
  return(likelihood)
}


poisson_rate <- function(x, p, steps = 100) {
  # poisson_rate computes the rate, lambda, such that:
  #   poisscdf(x,lambda,'upper') == p
  # steps (optional, default 100): number of Newton steps to take

  # Check if 'steps' argument is provided, default to 100 if not
  if (missing(steps)) {
    steps <- 100
  }

  if(x!=0){
    lambda <- x
    lfactx <- lfactorial(x)
    for (step in 1:steps) {
      cdf <- ppois(x, lambda)
      lambda <- lambda - (cdf + p - 1) / -exp(-lambda+x*log(lambda) - lfactx)
    }
  }else{
    lambda <- -log(1-p) # Solve (Poisson_PMF = p when x=0)
  }

  if(abs(ppois(x, lambda, lower.tail = F) - p) > 1e-6){
    print(paste0('Poisson CDF:', ppois(x, lambda, lower.tail = F)))
    print(paste0('\nInput Pvalue:', p))
    stop("\nThe rate does not correspond to the input pvalue, please check")
  }
  return(lambda)
}

poisson_pdf_adjusted <- function(x, lambda){
  # Replace the factorial term in the original Poisson PDF with Stirling's
  # approximation, to avoild potential overflowing issue:
  # https://en.wikipedia.org/wiki/Stirling%27s_approximation

  if(!is.na(lambda)){
    if(lambda > 50){
      pow = x*log(lambda) - lambda - x*log(x) + x
      p = exp(pow)
    }else{
      p = dpois(x, lambda)
    }
  }else{
    p  = NA
  }

  return(p)
}

poisson_uniform_likelihood_rvas <- function(genetic_data,
                                            component_endpoints,
                                            grid_size) {
  # Observe data for gene g:
  # x: number of minor alleles in casese (AC_cases)
  # p: p-value for association using SAIGE
  # gamma: beta = log(gamma) = log(e^beta)
  #
  # Suppose that SAIGE runs Poisson model with smart choice of lambda,
  # x|H_0 ~ Pois(lambda), p = poisscdf(x,lambda,'upper')
  #
  # Then we want to solve for lambda, and we can model the data as:
  # x|gamma ~ Pois(lambda * gamma) = Pois(lambda*e^beta)

  no_cpts = length(component_endpoints)
  no_genes = nrow(genetic_data)

  likelihood <- matrix(NA,
                       nrow = no_genes,
                       ncol = no_cpts)

  mu_grid = seq(0.05,1,by = 1/grid_size)
  case_count_matrix <- replicate(grid_size,genetic_data$AC_cases)

  # Calculate the rate for each mu value in the grid and each case count
  if(is.null(genetic_data$p_value)){
    warning('Pvalues from SAIGE-GENE Burden test is preferred for Poisson likelihood computation...')
    genetic_data$p_value <- 1 - ppois(abs(genetic_data$effect_estimate/genetic_data$effect_se))
  }else{
    genetic_data$p_value <- if_else(genetic_data$effect_estimate > 0, genetic_data$p_value/2, 1-genetic_data$p_value/2) # SAIGE-GENE Burden Pvalue
  }
  genetic_data$expected_count <- mapply(poisson_rate, genetic_data$AC_cases, genetic_data$p_value)

  for (kk in 1:no_cpts) {
    print(paste0('Computing likelihood for component ', kk, '/', no_cpts,'...'))
    # Expand case_count into a matrix by replicating the vector along the columns
    rate <- genetic_data$expected_count * t(replicate(no_genes,exp(mu_grid * component_endpoints[kk])))

    # Compute the likelihood for each case count and rate combination
    likelihoods <- mapply(poisson_pdf_adjusted, case_count_matrix, rate)
    likelihoods <- matrix(likelihoods, nrow = nrow(rate), ncol = ncol(rate))

    # Calculate the average likelihood for each test across the grid of mu values
    likelihood[, kk] <- rowMeans(likelihoods, na.rm = T)
  }
  return(likelihood)
}

likelihood_function_rvas <- function(trait_type){
  if(tolower(trait_type) == 'binary'){
    print('Using Poisson/Negative Binomial likelihood function for BINARY traits...')
    likelihood_func <- poisson_uniform_likelihood_rvas
  }else{
    print('Using Normal likelihood function for COUNTINUOUS traits...')
    likelihood_func <- normal_uniform_likelihood
  }
  return(likelihood_func)
}
