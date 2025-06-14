# --- Libraries and Sources ---
library(dplyr)
source("luke/distribution_functions.R")
source("R/estimate_heritability.R")
cat("Files sourced.\n")

# --- Configuration --- 
model_dir <- "burdenEM_example_output/"
trait_name <- "Height"
model_file_path <- paste0(model_dir, "burdenEM_fit_genebass_50_NA_pLoF.rds")

cat("--- Testing PMF Derivatives for Trait: Height ---\n")
cat(paste("Loading model from:", model_file_path, "\n"))
model <- readRDS(model_file_path)
if(is.null(model)){
  stop(paste0("Failed to load model: ", model_file_path))
}
cat("Model loaded.\n")

# --- Parameters for numerical derivative ---
component_index_to_test <- 1 
epsilon <- 1e-12
model_plus_eps <- model
model_plus_eps$delta[, component_index_to_test] <- model_plus_eps$delta[, component_index_to_test] + epsilon
x_values <- model$component_endpoints[model$component_endpoints >= 0]^2

cat(paste0("\nTesting derivatives for component_index = ", component_index_to_test, " with epsilon = ", epsilon, "\n"))

# --- Test d_gene_PMF_d_delta_i ---
# cat("\n--- Testing dg_dw ---\n")
analytical_dG_dw <- get_dG_dw(model, component_index_to_test)
gene_cdf_plus_eps <- get_gene_cdf_betasq(model_plus_eps)
gene_cdf_original <- get_gene_cdf_betasq(model)
numerical_dG_dw <- function(x) {(gene_cdf_plus_eps(x) - gene_cdf_original(x)) / epsilon}
comparison_df <- data.frame(
  x = x_values,
  analytical = analytical_dG_dw(x_values),
  numerical = numerical_dG_dw(x_values)
)
print(comparison_df)

# --- Test dv_dw --- 
cat("\n--- Testing dV_dw ---\n")
# Analytical derivative as a function
analytical_dV_dw_fn <- get_dV_dw(model, component_index_to_test)
V_fn <- get_variance_cdf_betasq(model)
V_fn_pluseps <- get_variance_cdf_betasq(model_plus_eps)
numerical_dV_dw_fn <- function(x) (V_fn_pluseps(x) - V_fn(x)) / epsilon

x_values <- c(0.01, 0.001, 0.0001, 0.00001)
comparison_df_dv_dw <- data.frame(
  grid_effect_test_point = x_values,
  analytical = analytical_dV_dw_fn(x_values),
  numerical = numerical_dV_dw_fn(x_values)
  )

print(comparison_df_dv_dw)

# --- Test get_variance_qf_betasq is inverse of get_variance_cdf_betasq ---
cat("\n--- Testing Variance QF and CDF inverse property ---\n")
VarCDF_fn <- get_variance_cdf_betasq(model, right_tail=TRUE)
VarQF_fn <- get_variance_qf_betasq(model, right_tail=TRUE)
p_values_for_inverse_test <- c(0.01, 0.1, 0.5, 0.9, 0.99)
q_from_qf <- VarQF_fn(p_values_for_inverse_test)
p_reconstructed_from_cdf <- VarCDF_fn(q_from_qf)
inverse_check_df <- data.frame(p_original = p_values_for_inverse_test, q_from_qf = q_from_qf, p_reconstructed = p_reconstructed_from_cdf)
print(inverse_check_df, digits=5)

# --- Test get_variance_pdf_betasq --- 
cat("\n--- Testing get_variance_pdf_betasq ---\n")

pdf_fn_var <- get_variance_pdf_betasq(model)
cdf_fn_var <- get_variance_cdf_betasq(model, right_tail=FALSE) # Use right_tail=FALSE for standard CDF definition
empirical_pdf_fn_var <- function(x) {
  (cdf_fn_var(x + epsilon) - cdf_fn_var(x)) / epsilon
}

x_values <- c(0.01, 0.001, 0.0001, 0.00001)
comparison_df_var_pdf <- data.frame(
  x = x_values,
  analytical = pdf_fn_var(x_values),
  numerical = empirical_pdf_fn_var(x_values)
)
print(comparison_df_var_pdf, digits=5)

# --- Test get_dVinv_dw --- 
cat("\n--- Testing get_dVinv_dw ---\n")
analytical_dVinv_dw_fn <- get_dVinv_dw(model, component_index_to_test)
Vinv_base_fn <- get_variance_qf_betasq(model, right_tail=TRUE) # Assuming right_tail=TRUE is standard for Vinv in this context
Vinv_plus_eps_fn <- get_variance_qf_betasq(model_plus_eps, right_tail=TRUE)

numerical_dVinv_dw_fn <- function(v) {
  (Vinv_plus_eps_fn(v) - Vinv_base_fn(v)) / epsilon
}

v_test_points <- p_values_for_inverse_test

comparison_df_dVinv_dw <- data.frame(
  v = v_test_points,
  analytical = analytical_dVinv_dw_fn(v_test_points),
  numerical = numerical_dVinv_dw_fn(v_test_points)
  )
print(comparison_df_dVinv_dw, digits=5)

# --- Test get_needed_genes_fn_derivative --- 
cat("\n--- Testing get_needed_genes_fn_derivative ---\n")
analytical_dF_dw_fn <- get_needed_genes_fn_derivative(model, component_index_to_test)

# Numerical derivative of get_needed_genes_fn w.r.t. weight of component_index_to_test
# Use the same model_prime and model_prime_perturbed
F_base_fn <- get_needed_genes_fn(model)
F_plus_eps_fn <- get_needed_genes_fn(model_plus_eps)

numerical_dF_dw_fn <- function(v) {
  (F_plus_eps_fn(v) - F_base_fn(v)) / epsilon
}

analytical_dF_dw_at_v <- analytical_dF_dw_fn(v_test_points)
numerical_dF_dw_at_v <- numerical_dF_dw_fn(v_test_points)

comparison_df_dF_dw <- data.frame(
  v = v_test_points,
  analytical = analytical_dF_dw_at_v,
  numerical = numerical_dF_dw_at_v,
  abs_diff = abs(analytical_dF_dw_at_v - numerical_dF_dw_at_v),
  rel_diff = ifelse(analytical_dF_dw_at_v == 0 & numerical_dF_dw_at_v == 0, 0, 
                  ifelse(analytical_dF_dw_at_v == 0, NA, abs(analytical_dF_dw_at_v - numerical_dF_dw_at_v) / abs(analytical_dF_dw_at_v)))
)
print(comparison_df_dF_dw, digits=5)

cat("\n--- Testing get_F_v_standard_errors ---\n")
needed_genes_var_fn <- get_needed_genes_var_fn(model)
needed_genes_fn <- get_needed_genes_fn(model)
needed_genes_var_at_v <- needed_genes_var_fn(v_test_points)
needed_genes_at_v <- needed_genes_fn(v_test_points)
se <- sqrt(diag(needed_genes_var_at_v))
print(data.frame(v = v_test_points, needed_genes = needed_genes_at_v, se = se), digits=5)

cat("--------------------------------------------\n")
cat("\n--- Test Complete ---
")