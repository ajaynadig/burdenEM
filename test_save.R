source("luke/main.R")

result <- run_burdenEM_rvas(
  variant_dir = "data/var_txt",
  variant_file_pattern = "^aou_afr_250_pLoF_.*\\.txt\\.bgz$",
  ld_corrected_scores_file = "data/genes/ukbb_ld_corrected_burden_scores_pLoF_0.0_0.001.tsv",
  output_file_prefix = "fitted_models/test_em_manual",
  annotation_to_process = "pLoF",
  num_iter = 1000,
  verbose = TRUE,
  optimizer = "EM"
)

cat("\nChecking saved file from run_burdenEM_rvas...\n")
if (file.exists("fitted_models/test_em_manual.rds")) {
  cat("File size:", file.info("fitted_models/test_em_manual.rds")$size, "\n")
  if (file.info("fitted_models/test_em_manual.rds")$size > 0) {
    test_load <- readRDS("fitted_models/test_em_manual.rds")
    cat("Successfully loaded model\n")
    cat("Model has components:", names(test_load), "\n")
  } else {
    cat("ERROR: File is empty!\n")
  }
} else {
  cat("ERROR: File was not created!\n")
}
