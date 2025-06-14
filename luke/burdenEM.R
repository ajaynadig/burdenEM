# burdenEM.R
# Contains core functions for fitting the BurdenEM model.

# Ensure necessary functions are available
# Assuming these are sourced in the main script or environment
# source("luke/gene_features.R")
# source("R/model.R")
# source("R/EM.R")

#' Fit BurdenEM Model (Grid-Based)
#'
#' Prepares features, initializes, and fits the grid-based BurdenEM model.
#'
#' @param gene_level_data Data frame containing gene-level summary statistics,
#'                        including 'effect_estimate', 'effect_se', and optionally
#'                        'gene_score' and a feature column.
#' @param feature_col_name Optional string: name of the column for gene features. Default: NULL.
#' @param num_feature_bins Integer: number of quantile bins for the feature column. Default: 5.
#' @param burdenem_grid_size Grid size for likelihood calculation. Default: 100.
#' @param num_iter Number of iterations for EM algorithm. Default: 100000.
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
    burdenem_grid_size = 100,
    num_iter = 10000,
    per_allele_effects = FALSE,
    verbose = FALSE,
    likelihood_function = function(row, beta_vec) {
        dnorm(row$effect_estimate, mean = beta_vec, sd = row$effect_se)
    },
    h2_function = function(beta, row) beta^2,
    base_component_endpoints = c(-5, -1, -0.2,  -0.05),
    customize_components=FALSE,
    drop_columns = c("effect_estimate", "effect_se")
) {

    # --- Initialize Grid-Based BurdenEM Model ---
    if(verbose) message("\n--- Initializing BurdenEM Model ---")

    # Define component endpoints based on per_allele_effects and data
    if(customize_components){
      lower_bound = unlist(quantile(gene_level_data$gamma_per_sd, probs = c(0.001)))
      upper_bound = unlist(quantile(gene_level_data$gamma_per_sd, probs = c(0.999)))
      bound = max(abs(min(gene_level_data$gamma_per_sd)), abs(max(gene_level_data$gamma_per_sd)))
      component_endpoints1 = c(-bound,
                               seq(lower_bound, upper_bound, length.out = (burdenem_no_cpts/2)),
                               bound)
      component_endpoints2 = seq(-bound, bound, length.out = burdenem_no_cpts/2)
      component_endpoints2[which(component_endpoints2 == 0)] = 1e-300
      component_endpoints <- sort(unique(c(component_endpoints1, component_endpoints2)))
    }else{
      component_endpoints <- c(base_component_endpoints, 0, -rev(base_component_endpoints))
    }

    if(!per_allele_effects){
        burden_score_mean <- mean(gene_level_data$burden_score, na.rm = TRUE)
        component_endpoints <- component_endpoints * sqrt(burden_score_mean)
    }
    if(verbose) { message("Component endpoints:"); print(component_endpoints) }

    # Ensure required model functions are available (e.g., source R/model.R, R/EM.R)
    if (!exists("initialize_grid_model")) stop("Function 'initialize_grid_model' not found. Source R/model.R")
    if (!exists("EM_fit_grid")) stop("Function 'EM_fit_grid' not found. Source R/EM.R")
    if (!exists("information_matrices")) stop("Function 'information_matrices' not found. Source R/EM.R")

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

    burdenem_model <- EM_fit_grid(burdenem_model, max_iter = num_iter)

    burdenem_model$information <- information_matrices(burdenem_model)

    return(burdenem_model)
}
