compute_polygenicity_rvas <- function(genetic_data,
                                      model){


  simulate_posterior <- function(dummy,weights,uniform_a,uniform_b){
    component = sample(c(1:length(weights)),size = 1, prob = weights)
    sample = runif(1, min = uniform_a[component], max = uniform_b[component])
  }

  # Sample from the effect size distribution
  a <- if_else(model$component_endpoints>0, 0, model$component_endpoints)
  b <- if_else(model$component_endpoints<0, 0, model$component_endpoints)

  beta <- unlist(sapply(1:ncol(model$features),
         function(x) {
           mixing_weights = model$delta[x,]
           data_x <- genetic_data[model$features[,x]==1,]
           n_genes <- nrow(data_x)
           sapply(1:n_genes, simulate_posterior, mixing_weights, a, b)
         }))

  # Compute the Kurtosis of Effect size distribution (the 4th moment of normal distribution integrals)
  kurtosis <- kurtosis(beta) + 3

  # Compute Polygenicity: Me = n_gene/(3*Kurtosis) (?) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6732528/
  n_genes <- nrow(genetic_data)
  polygenicity <- 3*n_genes/kurtosis

  return(polygenicity)
}
