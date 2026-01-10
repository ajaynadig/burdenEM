# R/burdenEM_rvas_analysis/check_synonymous_variance.R

# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

# Define the path to the specific synonymous variant file
# Adjust ANCESTRY and SAMPLE_N as needed based on your setup
ANCESTRY <- "NA" # Assuming "NA" based on previous logs
SAMPLE_N <- "50" # Assuming "50" based on previous logs
ANNOTATION <- "synonymous"
BIN_LOW <- "0"
BIN_HIGH <- "1e-05"
BIN_LABEL <- "0e+00_1e-05" # Label from filename

variant_file_path <- file.path(
  "/Users/lukeoconnor/Dropbox/burdenEM_results/data", # Base results directory
  "genebass", # Data source
  "var_txt", # Subdirectory for variant files
  paste0("genebass_", SAMPLE_N, "_", ANCESTRY, "_", ANNOTATION, "_low_", BIN_LOW, "_high_", BIN_HIGH, "_", BIN_LABEL, ".txt.bgz")
)

# Check if the file exists
if (!file.exists(variant_file_path)) {
  stop(paste("Error: Variant file not found at:", variant_file_path))
}

message(paste("Loading synonymous variants from:", variant_file_path))

# Load the data
# Use show_col_types = FALSE to suppress column type messages
synonymous_data <- read_tsv(variant_file_path, col_types = cols(
    chr = col_character(),
    pos = col_double(),
    ref = col_character(),
    alt = col_character(),
    rsid = col_character(),
    gene = col_character(),
    consequence = col_character(),
    impact = col_character(),
    AF = col_double(),
    se = col_double(),
    p = col_double(),
    beta = col_double(),
    tstat = col_double(),
    n = col_double(),
    ac = col_double(),
    pass = col_logical(),
    info = col_character(),
    burden_weight = col_double(),
    variant_variance = col_double(),
    functional_category = col_character()
), show_col_types = FALSE)

message(paste("Loaded", nrow(synonymous_data), "synonymous variants."))

# Calculate beta_per_sd = beta * sqrt(variant_variance)
# Use mixed case column names (AF, beta) and ensure variant_variance exists
synonymous_data_processed <- synonymous_data %>% 
  filter(!is.na(AF) & !is.na(beta) & !is.na(variant_variance)) %>% # Filter for required columns
  mutate(
    beta_orig = as.numeric(beta),           # Use beta (lowercase)
    variant_variance = as.numeric(variant_variance), # Ensure numeric
    # Calculate scaling factor from variant_variance, handle non-positive variance
    scaling_factor = ifelse(variant_variance > 0, sqrt(variant_variance), 0),
    beta_per_sd = beta_orig * scaling_factor
  ) %>% 
  # Filter out rows where beta_per_sd calculation resulted in NA/NaN/Inf
  filter(!is.na(beta_per_sd) & is.finite(beta_per_sd))

message(paste("Processed", nrow(synonymous_data_processed), "variants after calculating beta_per_sd and filtering."))

# Calculate the variance of beta_per_sd
if (nrow(synonymous_data_processed) > 1) {
  variance_beta_per_sd <- var(synonymous_data_processed$beta_per_sd, na.rm = TRUE)
  message("\n--- Variance of Beta per SD for Synonymous Variants (AF bin [0, 1e-05)) ---")
  print(variance_beta_per_sd)
} else {
  warning("Not enough valid variants (<= 1) to calculate variance.")
  message("\n--- Variance of Beta per SD for Synonymous Variants (AF bin [0, 1e-05)) ---")
  print(NA)
}

message("\nScript finished.")
