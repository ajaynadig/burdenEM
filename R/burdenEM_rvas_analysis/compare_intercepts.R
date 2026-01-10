# Script to compare BHR intercepts with intercepts from calculate_variant_intercept

# --- Libraries ---
library(dplyr)
library(readr)
library(tidyr) # For pivot_wider
library(purrr) # Added for map
library(stringr) # Added for str_match

# --- Source Files ---
source("R/intercept.R") # Assumes intercept.R is in the R directory relative to script
source("R/io.R")       # Added to load the new function

# --- Parameters ---
pheno <- "50_NA" # Phenotype to analyze
data_name <- "genebass"     # Data source (e.g., 'genebass', 'ukb')
anc <- "EUR"             # Ancestry
input_root <- path.expand('~/Dropbox/burdenEM_results/data/') # Root for input variant files
output_root <- path.expand('~/Dropbox/burdenEM_results/')    # Root where BHR results are saved

# Annotations used in the analysis (adjust if needed based on pipeline run)
# These will be used to filter the combined variant data later
annotations_to_compare <- c("pLoF", "missense|LC")

# Frequency bins for the new intercept function
frequency_bin_edges <- c(0, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1)

# --- Debugging Flags ---
# Set to target annotation or NULL to run all
# Example: debug_annotation <- "missense|LC" 
debug_annotation <- NULL # Set to NULL to run all
# Set to target bin string (matching cut() output) or NULL to run all bins for the annotation
# Example: debug_freq_bin <- "[0,1e-05)" 
debug_freq_bin <- NULL # Set to NULL to run all

# --- End Debugging Flags ---

# --- 1. Load BHR Intercepts ---
message("--- Loading BHR Intercepts ---")

# Correct path provided by user
bhr_file_path <- '/Users/lukeoconnor/Dropbox/burdenEM_results/BHR/genebass/genebass_eur_50_NA.csv'

bhr_intercept_df <- NULL # Initialize as NULL

if (file.exists(bhr_file_path)) {
  tryCatch({
    bhr_raw_data <- read_csv(bhr_file_path, show_col_types = FALSE)
    message(paste("Loaded BHR results from:", bhr_file_path))
    
    # Check required columns
    required_bhr_cols <- c("annotation", "AF_bin", "intercept")
    if (!all(required_bhr_cols %in% names(bhr_raw_data))) {
      stop(paste("Missing required columns in BHR file:", paste(setdiff(required_bhr_cols, names(bhr_raw_data)), collapse=", ")))
    }
    
    # Process BHR data
    bhr_intercept_df <- bhr_raw_data %>% 
      # Map BHR annotations to our categories
      mutate(
        functional_category = case_when(
          grepl("pLoF", annotation, ignore.case = TRUE) ~ "pLoF",
          grepl("missense", annotation, ignore.case = TRUE) ~ "missense|LC",
          TRUE ~ NA_character_ # Ignore other annotations if any
        ),
        # Standardize AF bins based on the observed format 'low_high' in the BHR file
        frequency_bin = case_when(
          # Map exact strings from BHR AF_bin column to cut() format
          AF_bin == "0e+00_1e-05" ~ "[0,1e-05)",
          AF_bin == "1e-05_1e-04" ~ "[1e-05,0.0001)",
          AF_bin == "1e-04_1e-03" ~ "[0.0001,0.001)",
          AF_bin == "1e-03_1e-02" ~ "[0.001,0.01)",
          # Add mappings for other bins if they exist in BHR output and are needed
          # AF_bin == "1e-02_0.05" ~ "[0.01,0.05)", # Example if this format existed
          # AF_bin == "0.05_0.1"  ~ "[0.05,0.1)",  # Example if this format existed
          # AF_bin == "0.1_0.5"   ~ "[0.1,0.5)",   # Example if this format existed
          # AF_bin == "0.5_1"     ~ "[0.5,1]",     # Example if this format existed
          TRUE ~ NA_character_ # Fallback if format doesn't match known patterns
        )
      ) %>% 
      filter(!is.na(functional_category) & functional_category %in% annotations_to_compare & !is.na(frequency_bin)) %>% 
      # Group by our categories and the standardized freq bin, average if multiple BHR annotations map to one category (e.g., missense types)
      group_by(functional_category, frequency_bin) %>% 
      summarize(bhr_intercept = mean(intercept, na.rm = TRUE), .groups = 'drop')
      
    if(nrow(bhr_intercept_df) == 0){
        message("Warning: No BHR intercepts found matching the target annotations and expected AF bin formats after processing.")
        bhr_intercept_df <- NULL # Set back to NULL if empty
    } else {
        message(paste("Successfully processed", nrow(bhr_intercept_df), "BHR intercepts for comparison."))
    }

  }, error = function(e) {
    message(paste("Error loading or processing BHR results:", e$message))
    bhr_intercept_df <- NULL # Ensure it's NULL on error
  })
} else {
  message(paste("BHR results file not found at:", bhr_file_path))
  bhr_intercept_df <- NULL
}

# If BHR loading failed or yielded no results, create an empty placeholder df
if (is.null(bhr_intercept_df)) {
    message("Proceeding without BHR data for comparison.")
    # Define structure for joining later, even if empty
    bhr_intercept_df <- data.frame(functional_category = character(), 
                                   frequency_bin = character(), 
                                   bhr_intercept = numeric())
}

# --- 2. Load Variant-Level Data ---
message("\n--- Loading Variant-Level Data ---")
# Corrected path based on list_dir and full_burden_pipeline_LO.R analysis
variant_dir <- file.path(input_root, data_name, 'var_txt')

# Call the centralized function from io.R
variant_data <- load_variant_files_with_category(
  variant_dir = variant_dir,
  data_name = data_name,
  pheno = pheno,
  annotations_to_process = annotations_to_compare
)

# Check if loading was successful and data is available
if (nrow(variant_data) == 0) {
  stop("No variant data loaded successfully. Check logs and file paths.")
}

# --- 3. Prepare Variant Data for Intercept Function ---
message("\n--- Preparing Variant Data ---")
# Select columns required by calculate_variant_intercept
# Functional category is now derived from filenames and already present

# Check for required columns (including the derived 'functional_category')
required_cols <- c("gene", "AF", "beta", "variant_variance", "functional_category") # Actual names + derived category
missing_cols <- setdiff(required_cols, names(variant_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in the combined variant data:", paste(missing_cols, collapse = ", ")))
}

variant_data_prepared <- variant_data %>% 
  # Filter for relevant functional categories (though already filtered during load, double-check)
  filter(functional_category %in% annotations_to_compare) %>% 
  # Ensure variant_variance is non-negative before sqrt
  filter(variant_variance >= 0) %>% 
  mutate(
    gene = as.character(gene), # Ensure gene is character
    AF = as.numeric(AF),
    beta_orig = as.numeric(beta), # Original beta
    se_orig_sq = as.numeric(variant_variance), # Original variance
    burden_weight_orig = 1, # Set original weight to 1
    se = sqrt(variant_variance), # Original SE, needed by main function
    
    # Calculate scaling factor from variant_variance, handle non-positive variance
    scaling_factor = ifelse(variant_variance > 0, sqrt(variant_variance), 0),
    beta_per_sd = beta_orig * scaling_factor # Renamed: This is beta * sqrt(var) = beta * sd

  ) %>% 
  # Select only the columns needed for intercept calculation AND main function
  # Main function needs: gene, se
  # Intercept function needs: beta_per_sd (and gene/AF/functional_category for grouping/debugging)
  select(gene, AF, beta_per_sd, se, functional_category) %>% 
  # Remove rows with NA/NaN/Inf in critical numeric columns AFTER calculation/selection
  filter(!is.na(AF) & !is.na(beta_per_sd) & !is.na(se) &
         is.finite(AF) & is.finite(beta_per_sd) & is.finite(se))

message(paste("Prepared data has", nrow(variant_data_prepared), "variants matching target annotations."))

if(nrow(variant_data_prepared) == 0) {
    stop("No variants found matching the specified functional categories after preparation.")
}

# --- 4. Calculate New Intercepts ---
message("\n--- Calculating New Intercepts ---")
new_intercepts_list <- list()

# Determine annotations to process based on debug flag
if (!is.null(debug_annotation) && debug_annotation %in% annotations_to_compare) {
  annotations_to_process <- debug_annotation
  message(paste("--- DEBUG MODE: Processing only annotation:", debug_annotation, "---"))
} else {
  annotations_to_process <- annotations_to_compare
  if (!is.null(debug_annotation)) {
    warning(paste("Debug annotation '", debug_annotation, "' not in annotations_to_compare. Running all.", sep=""))
  }
}

# Loop through the specified annotations
for (ann in annotations_to_process) {
  message(paste("... calculating for annotation:", ann))
  
  # Filter data for the current annotation
  current_annotation_data <- variant_data_prepared %>% 
    filter(functional_category == ann)
  
  if(nrow(current_annotation_data) == 0) {
    message(paste("No prepared variant data found for annotation:", ann, ". Skipping calculation."))
    new_intercepts_list[[ann]] <- list(intercepts_by_freq = NULL, N_genes = 0, N_variants = 0)
    next
  }
  
  # Call the function from intercept.R
  new_intercept_results <- calculate_variant_intercept(current_annotation_data, frequency_bin_edges)
  
  # Store results (the function now returns a list)
  new_intercepts_list[[ann]] <- new_intercept_results
}

# --- 5. Print and Compare Results ---
message("\n--- Comparison Results ---")

# Prepare New Intercepts data in a long format for joining
new_intercept_long_list <- list()
for(ann in names(new_intercepts_list)) {
  results <- new_intercepts_list[[ann]]
  
  # Check if intercepts_by_freq is NULL or empty
  if (is.null(results$intercepts_by_freq) || length(results$intercepts_by_freq) == 0) {
    message(paste("No new intercept results found for annotation:", ann))
    next # Skip to next annotation
  }
  
  # Print details for this annotation
  message(paste("\nAnnotation:", ann))
  message(paste("  Number of genes included:", results$N_genes))
  message(paste("  Number of variants included:", results$N_variants))
  message("  Intercepts per frequency bin:")
  print(results$intercepts_by_freq)
  
  # Convert named vector to data frame, extracting frequency bin from name
  # Name format is assumed to be 'functional_category_frequency_bin'
  intercept_names <- names(results$intercepts_by_freq)
  # Calculate the start position after the prefix (ann + '_')
  prefix_length <- nchar(ann) + 1
  # Extract the substring starting after the prefix
  freq_bins_extracted <- substr(intercept_names, prefix_length + 1, nchar(intercept_names))

  df <- data.frame(frequency_bin = freq_bins_extracted,
                     new_intercept = results$intercepts_by_freq,
                     functional_category = ann,
                     row.names = NULL) # Ensure row names are reset
  new_intercept_long_list[[length(new_intercept_long_list) + 1]] <- df
}

new_intercept_df <- bind_rows(new_intercept_long_list)

# --- Join BHR and New Intercepts for Side-by-Side Comparison --- 
message("\n--- Side-by-Side Comparison (BHR vs. New Intercept by Bin) ---")

if (nrow(new_intercept_df) > 0) {
  # Ensure frequency_bin columns are character for joining
  bhr_intercept_df <- bhr_intercept_df %>% mutate(frequency_bin = as.character(frequency_bin))
  new_intercept_df <- new_intercept_df %>% mutate(frequency_bin = as.character(frequency_bin))
  
  # Perform a full join
  comparison_df <- full_join(bhr_intercept_df, new_intercept_df, by = c("functional_category", "frequency_bin"))
  
  # Order for readability - Create factor levels based on unique bins found
  unique_bins_ordered <- unique(sort(c(as.character(bhr_intercept_df$frequency_bin), as.character(new_intercept_df$frequency_bin))))
  
  comparison_df <- comparison_df %>% 
    filter(!is.na(frequency_bin)) %>% # Remove rows where bin might be NA after join
    mutate(freq_bin_factor = factor(frequency_bin, levels = unique_bins_ordered)) %>% 
    arrange(functional_category, freq_bin_factor) %>% 
    select(functional_category, frequency_bin, bhr_intercept, new_intercept)
  
  print(comparison_df)
  
} else {
  message("No new intercept data was generated, cannot create comparison table.")
}

message("\nScript finished.")
