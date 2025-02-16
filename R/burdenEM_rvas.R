source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/io.R')
source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/model.R')
source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/likelihoods.R')
source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/EM.R')
source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/estimate_heritability.R')
source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/qqplot.R')

#burdenEM_rvas
burdenEM_rvas <- function(input_data,
                          features = NULL,
                          component_endpoints = NULL,
                          no_cpts = 10,
                          grid_size = 10,
                          heritability_est = TRUE,
                          polygenicity_est = TRUE,
                          num_iter =500,
                          prevalence = NULL,
                          bootstrap = TRUE,
                          n_boot = 100,
                          null_sim = TRUE,
                          n_null = 100,
                          return_likelihood = FALSE,
                          qq_plot=FALSE,
                          estimate_posteriors = FALSE) {


  if (is.null(features)) {
    features <- matrix(1, nrow = nrow(input_data), ncol = 1)
  }
  genetic_data = process_data_rvas(input_data,
                                   features)

  component_endpoints = choose_component_endpoints_rvas(component_endpoints,
                                                        no_cpts,
                                                        genetic_data)
  cat("...initializing model")
  model = initialize_model(likelihood_function = likelihood_function_rvas(trait_type=unique(input_data$trait_type)),
                           genetic_data = genetic_data,
                           component_endpoints = component_endpoints,
                           features = features,
                           grid_size = grid_size)


  #Full data EM
  cat("...running EM in full dataset")
  model = EM_fit(model,
                 num_iter)

  #Bootstrap EM
  if (bootstrap) {
    model$bootstrap_output = bootstrap_EM(model,
                                          n_boot,
                                          num_iter)

  }

  if (null_sim) {
    model$null_delta = null_EM_rvas(genetic_data,
                                    model,
                                    num_iter,
                                    n_null,
                                    grid_size)
  }

  if(polygenicity_est){
    model$polygenicity <- estimate_polygenicity_rvas(genetic_data = genetic_data, model = model)
  }

  if (heritability_est) {
    print('... Estimate heritability in the full dataset')
    model$heritability_output = estimate_heritability_rvas(model = model,
                                                           genetic_data = genetic_data)

    #bootstrap heritability estimation
    if (bootstrap) {
      bootstrap_heritability_output <- lapply(1:n_boot,
                                              function(iter) {
                                                model_boot = model
                                                model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                                boot_heritability = estimate_heritability_rvas(model = model_boot,
                                                                                               genetic_data = genetic_data[model$bootstrap_output$bootstrap_samples[,iter],])

                                              })

      bootstrap_h2_ests = sapply(1:length(bootstrap_heritability_output), function(x) bootstrap_heritability_output[[x]]$total_h2)
      heritability_CI = quantile(bootstrap_h2_ests,c(0.025,0.975))

      bootstrap_annoth2_ests = sapply(1:length(bootstrap_heritability_output), function(x) bootstrap_heritability_output[[x]]$annot_h2)
      annot_h2_CI = sapply(1:nrow(bootstrap_annoth2_ests),
                           function(i) {
                             quantile(bootstrap_annoth2_ests[i,],c(0.025,0.975))
                           })


      bootstrap_frach2_ests = sapply(1:length(bootstrap_heritability_output), function(x) bootstrap_heritability_output[[x]]$frac_h2)
      frach2_CI =  sapply(1:nrow(bootstrap_frach2_ests),
                          function(i) {
                            quantile(bootstrap_frach2_ests[i,],c(0.025,0.975))
                          })

      bootstrap_enrich_ests = sapply(1:length(bootstrap_heritability_output), function(x) bootstrap_heritability_output[[x]]$enrichment)
      enrich_CI = sapply(1:nrow(bootstrap_enrich_ests),
                         function(i) {
                           quantile(bootstrap_enrich_ests[i,],c(0.025,0.975))
                         })
      model$heritability_output$heritability_CI = heritability_CI
      model$heritability_output$annot_h2_CI = annot_h2_CI
      model$heritability_output$frach2_CI = frach2_CI
      model$heritability_output$enrich_CI = enrich_CI

    }

    #null heritability estimates
    if (null_sim) {
      null_h2_ests <- sapply(1:n_null,
                             function(iter) {
                               model_null = model
                               model_null$delta = model$null_delta[[iter]]

                               null_heritability = estimate_heritability_rvas(model = model_null,
                                                                              genetic_data = genetic_data)

                               return(null_heritability$total_h2)

                             })

      model$heritability_output$null_h2_ests = null_h2_ests
      model$heritability_output$total_h2_p = mean(model$heritability_output$null_h2_ests > model$heritability_output$total_h2 )
    }
  }

  # Generate data for QQplots
  if(qq_plot){
      model$qq_data <- compute_estimate_expected_rvas(genetic_data, model)
  }

  #Get some posterior expectations
  if(estimate_posteriors){

  } # TODO: edit

  return(model)

}
