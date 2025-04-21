#heritability estimation functions
estimate_heritability_trio <- function(model,
                                       genetic_data,
                                       prevalence,
                                       gamma_scaling_factor = 1) {
  #  print(model$delta)
  annot_h2 <- sapply(1:ncol(model$features),
                     function(x) {
                       # print(model$delta)
                       mixing_weights = model$delta[x,]
                       mixing_weights[mixing_weights < 0] <- 0

                       #sample rate ratios
                       endpoint_samples = sample(model$component_endpoints, size = 100000, prob =mixing_weights, replace = TRUE)
                       loggamma_samples = runif(100000,min = pmin(0,endpoint_samples), max = pmax(0,endpoint_samples)) * gamma_scaling_factor
                       # gamma_samples_trunc = pmin(exp(loggamma_samples),1/prevalence)
                       gamma_samples = exp(loggamma_samples)

                       #sample mutation rates
                       mu_samples <- sample(genetic_data$case_rate[model$features[,x] ==1], size = 100000, replace = TRUE)

                       heritability = sum(model$features[,x]) * (prevalence * mean((gamma_samples -1)^2 * 2*mu_samples))/(1-prevalence)
                     })



  frac_h2 = annot_h2/sum(annot_h2)
  frac_expected = sapply(1:ncol(model$features),
                         function(x) {
                           sum(genetic_data$case_rate[model$features[,x] == 1])
                         })/sum(genetic_data$case_rate)

  enrichment = frac_h2/frac_expected


  return(list(total_h2 = sum(annot_h2),
              annot_h2 = annot_h2,
              frac_h2 = frac_h2,
              frac_expected = frac_expected,
              enrichment = enrichment))
}

estimate_heritability_rvas <- function(model,
                                       genetic_data,
                                       per_allele_effects = FALSE){

  if(per_allele_effects) {
    heritability_fn <- function(weights, endpoints, multipliers) {
      heritability = t(weights %*% (endpoints^2)) %*% multipliers / 3
    }
  } else {
    heritability_fn <- function(weights, endpoints, multipliers) {
      heritability = sum((weights %*% (endpoints^2))) / 3
    }
  }
  mixing_weights <- model$features %*% model$delta

  annot_h2 <- sapply(1:ncol(model$features),
                     function(x) {
                       rows <- model$features[,x]==1
                       heritability <- heritability_fn(mixing_weights[rows,], 
                                                      model$component_endpoints, 
                                                      as.matrix(genetic_data$burden_score[rows]))
                     })
  print(annot_h2)
  ids <- which(model$component_endpoints>0)
  positive_h2 <- heritability_fn(mixing_weights[,ids], 
                                  model$component_endpoints[ids], 
                                  as.matrix(genetic_data$burden_score))
  print(positive_h2)
  ids <- which(model$component_endpoints<0)
  negative_h2 <- heritability_fn(mixing_weights[,ids], 
                                  model$component_endpoints[ids], 
                                  as.matrix(genetic_data$burden_score))
  print(negative_h2)

  total_h2 = heritability_fn(mixing_weights, 
                                  model$component_endpoints, 
                                  as.matrix(genetic_data$burden_score))
  print(total_h2)
  stopifnot(all.equal(positive_h2 + negative_h2, total_h2, tolerance=1e-6))


  frac_h2 = annot_h2/total_h2
  prop_positive_h2 = sum(positive_h2)/total_h2
  prop_negative_h2 = sum(negative_h2)/total_h2
  frac_expected = sapply(1:ncol(model$features),
                         function(x) {
                           sum(genetic_data$burden_score[model$features[,x] == 1])
                         })/sum(genetic_data$burden_score)

  enrichment = frac_h2/frac_expected

  return(
    list(total_h2 = total_h2,
         annot_h2 = annot_h2,
         frac_h2 = frac_h2,
         prop_positive_h2 = prop_positive_h2,
         prop_negative_h2 = prop_negative_h2,
         frac_expected = frac_expected,
         enrichment = enrichment)
  )

}
