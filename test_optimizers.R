source("luke/main.R")

# Test with a single study
variant_pattern <- "^aou_afr_250_pLoF_.*\\.txt\\.bgz$"
genes_file <- "data/genes/ukbb_ld_corrected_burden_scores_pLoF_0.0_0.001.tsv"

cat("\n=== Running with EM optimizer (max_iter=10000) ===\n")
result_em <- run_burdenEM_rvas(
  variant_dir = "data/var_txt",
  variant_file_pattern = variant_pattern,
  ld_corrected_scores_file = genes_file,
  output_file_prefix = "fitted_models/test_em_10k",
  annotation_to_process = "pLoF",
  num_iter = 10000,
  verbose = FALSE,
  optimizer = "EM"
)

cat("\n=== Running with mixsqp optimizer ===\n")
result_mixsqp <- run_burdenEM_rvas(
  variant_dir = "data/var_txt",
  variant_file_pattern = variant_pattern,
  ld_corrected_scores_file = genes_file,
  output_file_prefix = "fitted_models/test_mixsqp",
  annotation_to_process = "pLoF",
  num_iter = 5000,  # mixsqp has its own convergence criteria
  verbose = FALSE,
  optimizer = "mixsqp"
)

cat("\n=== Comparing Results ===\n")
cat("EM log-likelihood:", result_em$ll, "\n")
cat("mixsqp log-likelihood:", result_mixsqp$ll, "\n")
cat("LL difference:", abs(result_em$ll - result_mixsqp$ll), "\n")

cat("\nDelta (mixture proportions) max absolute difference:", 
    max(abs(result_em$delta - result_mixsqp$delta)), "\n")

# Calculate heritability for both
source("R/estimate_heritability.R")

h2_em <- estimate_heritability_components(result_em)
h2_mixsqp <- estimate_heritability_components(result_mixsqp)

cat("\n=== Heritability Comparison ===\n")
cat("EM total_h2:", h2_em$total_h2$mean, "±", h2_em$total_h2$se, "\n")
cat("mixsqp total_h2:", h2_mixsqp$total_h2$mean, "±", h2_mixsqp$total_h2$se, "\n")
cat("H2 difference:", abs(h2_em$total_h2$mean - h2_mixsqp$total_h2$mean), "\n")

cat("\nEM positive_h2:", h2_em$positive_h2$mean, "±", h2_em$positive_h2$se, "\n")
cat("mixsqp positive_h2:", h2_mixsqp$positive_h2$mean, "±", h2_mixsqp$positive_h2$se, "\n")
cat("Positive H2 difference:", abs(h2_em$positive_h2$mean - h2_mixsqp$positive_h2$mean), "\n")

if (abs(h2_em$total_h2$mean - h2_mixsqp$total_h2$mean) < 0.01) {
  cat("\n✓ PASS: H2 estimates agree within 0.01\n")
} else {
  cat("\n✗ FAIL: H2 estimates differ by more than 0.01\n")
}
