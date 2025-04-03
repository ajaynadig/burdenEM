#burdenEM_trio
burdenEM_trio <- function(input_data,
                          features = NULL,
                          component_endpoints = NULL,
                          no_cpts = 10,
                          grid_size = 10,
                          heritability_est = TRUE,
                          max_iter =10000,
                          max_iter_boot = 1000,
                          tol = 1e-6,
                          prevalence = NULL,
                          bootstrap = TRUE,
                          bootstrap_samples = NULL,
                          n_boot = 100,
                          null_sim = TRUE,
                          n_null = 100,
                          return_likelihood = TRUE,
                          estimate_posteriors = FALSE,
                          estimate_effective_penetrance = TRUE) {



  genetic_data = process_data_trio(input_data,
                                   features)

  component_endpoints = choose_component_endpoints_trio(component_endpoints,
                                                        no_cpts,
                                                        prevalence)
  cat("...initializing model")
  model = initialize_model(likelihood_function = poisson_uniform_likelihood,
                           genetic_data = genetic_data,
                           component_endpoints = component_endpoints,
                           features = features,
                           grid_size = grid_size)


  #Full data EM
  cat("...running EM in full dataset")
  model = EM_fit(model,
                 max_iter,
                 tol = tol,
                 return_likelihood = return_likelihood)

  #Bootstrap EM
  if (bootstrap) {
    model$bootstrap_output = bootstrap_EM(model,
                                          n_boot,
                                          max_iter_boot,
                                          bootstrap_samples)

  }

  if (null_sim) {
    model$null_delta = null_EM_trio(genetic_data,
                                    model,
                                    max_iter,
                                    n_null,
                                    grid_size)
  }

  if (heritability_est) {
    #estimate heritability in the full dataset
    model$heritability_output = estimate_heritability_trio(model = model,
                                                           genetic_data = genetic_data,
                                                           prevalence = prevalence)

    #bootstrap heritability estimation
    if (bootstrap) {
      bootstrap_heritability_output <- lapply(1:n_boot,
                                              function(iter) {



                                                model_boot = model
                                                model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                                boot_heritability = estimate_heritability_trio(model = model_boot,
                                                                                               genetic_data = genetic_data[model$bootstrap_output$bootstrap_samples[,iter],],
                                                                                               prevalence = prevalence)

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

                               null_heritability = estimate_heritability_trio(model = model_null,
                                                                              genetic_data = genetic_data,
                                                                              prevalence = prevalence)

                               return(null_heritability$total_h2)

                             })

      model$heritability_output$null_h2_ests = null_h2_ests
      model$heritability_output$total_h2_p = mean(model$heritability_output$null_h2_ests > model$heritability_output$total_h2 )
    }


  }

  #Get some posterior expectations and hypothesis tests
  if (estimate_posteriors == TRUE) {
    #first, get the naive rate ratio estimates
    RR_naive = genetic_data$case_count / genetic_data$expected_count

    #get some simple poisson p values
    RR_poisson_p = ppois(genetic_data$case_count,
                         genetic_data$expected_count,
                         lower.tail = FALSE)

    #get posterior means
    RR_posterior_means <- posterior_expectation(model,
                                                genetic_data,
                                                exp,
                                                grid_size)

    #hypothesis test
    RR_problessthanequal0 <- posterior_expectation(model,
                                                   genetic_data,
                                                   function(x) {x <= 0},
                                                   grid_size)

    posterior_gene_estimate_df <- data.frame(Case_Count = genetic_data$case_count,
                                             Expected_Count = genetic_data$expected_count,
                                             Estimate = RR_naive,
                                             Estimate_Poisson_P = RR_poisson_p,
                                             Posterior_Mean = RR_posterior_means,
                                             Posterior_ProbLessEqualZero = RR_problessthanequal0)

    rownames(posterior_gene_estimate_df) <- rownames(genetic_data)

    model$posterior_gene_estimates = posterior_gene_estimate_df

  }

  if (estimate_effective_penetrance) {
    cat("...computing effective penetrance")

    peneff = effective_penetrance_func(model,
                                       genetic_data,
                                       prevalence)
    peneff_CI = NA
    if (bootstrap) {
      cat("...bootstrap effective penetrance")
      bootstrap_peneff_ests = sapply(1:length(model$bootstrap_output$bootstrap_delta),
                                     function(iter) {
                                       model_boot = model
                                       model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[,iter],]
                                       model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_samples[,iter],]
                                       model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                       genetic_data_boot = genetic_data[model$bootstrap_output$bootstrap_samples[,iter],]

                                       effective_penetrance_func(model_boot,
                                                                 genetic_data_boot,
                                                                 prevalence)
                                     })

      peneff_CI = quantile(bootstrap_peneff_ests,c(0.025,0.975))
    }

    model$penetrance = list(effective_penetrance = peneff,
                            effective_penetrance_CI = peneff_CI)
  }

  return(model)

}
