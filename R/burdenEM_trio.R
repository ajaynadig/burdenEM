#burdenEM_trio
burdenEM_trio <- function(input_data,
                     features = NULL,
                     component_endpoints = NULL,
                     no_cpts = 10,
                     grid_size = 10,
                     heritability_est = TRUE,
                     num_iter =500,
                     prevalence = NULL,
                     bootstrap = TRUE,
                     n_boot = 100,
                     null_sim = TRUE,
                     n_null = 100,
                     return_likelihood = FALSE,
                     estimate_posteriors = FALSE) {



  genetic_data = process_data_trio(input_data,
                                   features)

  component_endpoints = choose_component_endpoints_trio(component_endpoints,
                                                        no_cpts,
                                                        prevalence)

  model = initialize_model(likelihood_function = poisson_uniform_likelihood,
                           genetic_data = genetic_data,
                           component_endpoints = mixture_params,
                           features = features,
                           grid_size = grid_size)

  #Full data EM
  model = EM_fit(model,
                 num_iter)

  #Bootstrap EM
  if (bootstrap) {
    model$bootstrap_delta = bootstrap_EM(model,
                                         n_boot,
                                         num_iter)

  }

  if (null_sim) {
    model$null_delta = null_EM_trio(genetic_data,
                                    model,
                                    num_iter,
                                    n_null)
  }



}
