# burdenEM.R
# Contains core functions for fitting the BurdenEM model.

# Ensure necessary functions are available
# Assuming these are sourced in the main script or environment
# source("luke/gene_features.R")
# source("R/model.R")
# source("R/EM.R")

choose_model_endpoints <- function(num_positive_components, 
                                    num_genes, 
                                    min_relevant_h2=0.01,
                                    max_gene_h2=0.1){
    smallest_positive_endpt = sqrt(min_relevant_h2/num_genes)
    largest_positive_endpt = sqrt(max_gene_h2)
    positive_endpts = exp(seq(log(smallest_positive_endpt), log(largest_positive_endpt), length.out = num_positive_components))
    # positive_endpts = c(.1, .2, .5, 1, 2, 5)
    return (c(-rev(positive_endpts), 0, positive_endpts))
}

#' Fit BurdenEM Model (Grid-Based)
#'
#' Prepares features, initializes, and fits the grid-based BurdenEM model.
#'
#' @param gene_level_data Data frame containing gene-level summary statistics,
#'                        including 'effect_estimate', 'effect_se', and optionally
#'                        'gene_score' and a feature column.
#' @param feature_col_name Optional string: name of the column for gene features. Default: NULL.
#' @param num_feature_bins Integer: number of quantile bins for the feature column. Default: 5.
#' @param burdenem_grid_size Grid size for likelihood calculation.
#' @param num_iter Number of iterations for EM algorithm.
#' @param per_allele_effects Logical, whether input effects are per-allele. Default: FALSE.
#' @param verbose Logical, whether to print detailed messages. Default: FALSE.
#' @param likelihood_fn Function defining the row-wise likelihood. Default: Normal likelihood.
#' @param base_component_endpoints Numeric vector for initial component endpoint values before scaling.
#'                                Default: c(-1, -0.5, -0.2, -0.1, -0.05, -0.02).
#'
#' @return A list object representing the fitted BurdenEM model.
#' @export
fit_burdenem_model <- function(
    gene_level_data,
    burdenem_grid_size,
    num_iter,
    per_allele_effects = FALSE,
    verbose = FALSE,
    likelihood_function = function(row, beta_vec) {
        dnorm(row$effect_estimate, mean = beta_vec, sd = row$effect_se)
    },
    h2_function = function(beta, row) beta^2,
    num_positive_components = 10,
    drop_columns = c("effect_estimate", "effect_se"),
    optimizer = "EM"
) {

    # --- Initialize Grid-Based BurdenEM Model ---
    if(verbose) message("\n--- Initializing BurdenEM Model ---")

    # Define component endpoints based on per_allele_effects and data
    component_endpoints <- choose_model_endpoints(num_positive_components, 
                                                  num_genes = nrow(gene_level_data), 
                                                  min_relevant_h2=0.001,
                                                  max_gene_h2=0.1)
    if(per_allele_effects){
        burden_score_mean <- mean(gene_level_data$burden_score, na.rm = TRUE)
        component_endpoints <- component_endpoints / sqrt(burden_score_mean)
    }
    if(verbose) { message("Component endpoints:"); message(component_endpoints) }

    burdenem_model <- initialize_grid_model(
        gene_data = gene_level_data,
        likelihood_fn = likelihood_function,
        component_endpoints = component_endpoints,
        h2_function = h2_function,
        grid_points_per_component = burdenem_grid_size,
        drop_columns = drop_columns
    )

    # --- Fit Model using new EM ---
    if(verbose) message("\n--- Running EM Fit (grid) ---")

    burdenem_model <- fit_mixture_model(burdenem_model, optimizer = optimizer, max_iter = num_iter)

    burdenem_model$information <- information_matrices(burdenem_model)

    return(burdenem_model)
}
