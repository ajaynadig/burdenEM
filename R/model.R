#Functions to manipulate burdenEM models

choose_component_endpoints_trio = function(component_endpoints,
                                 no_cpts,
                                 prevalence) {
  if (!is.null(component_endpoints)) {
    return(component_endpoints)
  } else {
    component_endpoints = seq(0,log(1/prevalence),length.out = no_cpts)
  }
}

initialize_model <- function(likelihood_function, genetic_data, component_endpoints,features, grid_size){

  conditional_likelihood = likelihood_function(genetic_data, component_endpoints,grid_size)

  no_cpts = length(component_endpoints)

  if (is.null(features)) {
    features <- matrix(1, nrow = nrow(genetic_data), ncol = 1)
    rownames(features) = rownames(genetic_data)
  }

  delta_init = matrix(1, nrow = ncol(features), ncol = no_cpts)

  model = list(component_endpoints = component_endpoints,
               delta = delta_init,
               conditional_likelihood = conditional_likelihood,
               features = features,
               grid_size = grid_size)
}
