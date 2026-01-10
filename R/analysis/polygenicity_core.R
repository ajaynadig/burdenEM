estimate_polygenicity_rvas <- function(genetic_data,
                                       model){
  library(e1071)
  n_genes <- nrow(genetic_data)
  weights_agg <- colSums(model$delta)/5
  no_draws <- 1e6

  # Generate draws in a vectorized way
  draws <- unlist(sapply(1:length(model$component_endpoints),
                         function(i){
                           model$component_endpoints[i] * runif(as.integer(floor(no_draws * weights_agg[i])))
                         }))
  cumh2 <- cumsum(sort(draws^2, decreasing = TRUE))  # Cumulative sum of sorted squared samples
  first <- which(cumh2 >= cumh2[length(cumh2)] / 2)[1]  # Find the first index where cumulative h2 reaches half
  # Sort samples squared in descending order
  x <- sort(draws^2, decreasing = TRUE)
  print(paste0('Number of draws:', length(draws)))

  polygenicity_50 <- first / length(draws) * n_genes
  polygenicity_eff  = 3*n_genes*mean(draws^2)^2 / mean(draws^4)
  n_large_gene_2e_2 <- sum(abs(draws) > 0.02)/length(draws) * n_genes
  n_large_gene_1e_2 <- sum(abs(draws) > 0.01)/length(draws) * n_genes
  n_large_gene_1e_1 <- sum(abs(draws) > 0.1)/length(draws) * n_genes
  return(list(
    polygenicity_eff = polygenicity_eff,
    polygenicity_50 = polygenicity_50,
    n_large_gene_2e_2 = n_large_gene_2e_2,
    n_large_gene_1e_2 = n_large_gene_1e_2,
    n_large_gene_1e_1 = n_large_gene_1e_1
  ))
}
