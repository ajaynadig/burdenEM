#heritability estimation functions
estimate_heritability_trio <- function(model,
                                       genetic_data,
                                       prevalence) {
  annot_h2 <- sapply(1:ncol(model$features),
                     function(x) {
                       mixing_weights = model$delta[x,]
                       mixing_weights[mixing_weights < 0] <- 0

                       #sample rate ratios
                       loggamma_samples = sample(model$component_endpoints, size = 100000, prob =mixing_weights, replace = TRUE)
                       gamma_samples_trunc = pmin(exp(loggamma_samples),1/prevalence)

                       #sample mutation rates
                       mu_samples <- sample(genetic_data$case_rate[model$features[,x] ==1], size = 100000, replace = TRUE)

                       heritability = sum(model$features[,x]) * (prevalence * mean((gamma_samples_trunc -1)^2 * mu_samples))/(1-prevalence)
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
                                       genetic_data){
  per_allele_factor <- unique(genetic_data$per_allele_factor)

  annot_h2 <- sapply(1:ncol(model$features),
                     function(x) {
                       mixing_weights = model$delta[x,]
                       data_x <- genetic_data[model$features[,x]==1,]
                       no_genes <- nrow(data_x)
                       heritability = sum(mixing_weights%*% (model$component_endpoints)^2*no_genes/3*per_allele_factor)
                     })

  positive_h2 <- sapply(1:ncol(model$features),
                             function(x) {
                               ids <- which(model$component_endpoints>0)
                               mixing_weights = model$delta[x,ids]
                               data_x <- genetic_data[model$features[,x]==1,]
                               no_genes <- nrow(data_x)
                               heritability = sum(mixing_weights%*% (model$component_endpoints[ids])^2*no_genes/3*per_allele_factor)
                             })
  negative_h2 <- sapply(1:ncol(model$features),
                        function(x) {
                          ids <- which(model$component_endpoints<0)
                          mixing_weights = model$delta[x,ids]
                          data_x <- genetic_data[model$features[,x]==1,]
                          no_genes <- nrow(data_x)
                          heritability = sum(mixing_weights%*% (model$component_endpoints[ids])^2*no_genes/3*per_allele_factor)
                        })



  total_h2 = sum(annot_h2)
  frac_h2 = annot_h2/total_h2
  prop_positive_h2 = sum(positive_h2)/total_h2
  prop_negative_h2 = sum(negative_h2)/total_h2
  frac_expected = sapply(1:ncol(model$features),
                         function(x) {
                           sum(genetic_data$CAF[model$features[,x] == 1], na.rm = T)
                         })/sum(genetic_data$CAF, na.rm = T)

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
