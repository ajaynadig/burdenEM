#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

source("R/estimate_heritability.R")

# Load AoU phenotype information
aou_phenotype_info <- read_csv("tables/aou_phenotype_info.csv", show_col_types = FALSE)

# Hard-coded parameters
annotation <- "synonymous"

# Derive trait_datasets from aou_phenotype_info

# Ensure required columns are present and not NA for genebass definitions
genebass_definitions <- aou_phenotype_info %>%
  filter(
    !is.na(description) & description != "",
    !is.na(ukb_phenoname) & ukb_phenoname != ""
  ) %>%
  distinct(description, ukb_phenoname)

genebass_entries <- lapply(1:nrow(genebass_definitions), function(i) {
  row <- genebass_definitions[i, ]
  list(
    name = row$description,       # Human-readable name from aou_phenotype_info
    data_name = "genebass",
    pheno = row$ukb_phenoname,    # UKB code for file path
    is_aou = FALSE
  )
})

# Ensure required columns are present and not NA for AoU definitions
valid_aou_pheno_info <- aou_phenotype_info %>%
  filter(
    !is.na(description) & description != "",
    !is.na(phenoname) & phenoname != "",
    !is.na(ancestry) & ancestry != ""
    # Consider adding: & ancestry %in% c("eur", "afr", "amr") if only these are desired
  )

aou_entries <- lapply(1:nrow(valid_aou_pheno_info), function(i) {
  row <- valid_aou_pheno_info[i, ]
  list(
    name = row$description,       # Human-readable name from aou_phenotype_info
    data_name = paste0("aou_", row$ancestry),
    pheno = row$phenoname,        # AoU specific code for file path
    is_aou = TRUE,
    ukbb_pheno = row$ukb_phenoname # Store corresponding UKB code for reference
  )
})

# Combine genebass and AoU entries
trait_datasets <- c(genebass_entries, aou_entries)
# trait_datasets <- c(genebass_entries)

# Function to process a single trait-dataset
process_trait <- function(trait_name, data_name, pheno, is_aou = FALSE, ukbb_pheno = NULL) {
  message(sprintf("Processing: %s (%s) - %s", trait_name, data_name, pheno))
  
  # Get the model file path based on dataset type
  if (is_aou) {
    aou_ancestry <- sub("aou_", "", data_name)
    model_file <- sprintf("fitted_models/burdenEM_fit_aou_%s_%s_%s.rds", 
                         aou_ancestry, pheno, annotation)
  } else {
    model_file <- sprintf("fitted_models/burdenEM_fit_%s_%s_%s.rds", 
                         data_name, pheno, annotation)
  }
  
  if (!file.exists(model_file)) {
    message(sprintf("  Model file not found: %s", model_file))
    return(NULL)
  }
  
  tryCatch({
    model <- readRDS(model_file)
    h2_output <- estimate_heritability_components(model, verbose = FALSE)
    
    # Extract basic heritability estimates
    result <- data.frame(
      trait = trait_name,
      dataset = data_name,
      total_h2 = h2_output$total_h2$mean,
      total_h2_se = h2_output$total_h2$se,
      positive_h2 = h2_output$positive_h2$mean,
      positive_h2_se = h2_output$positive_h2$se,
      negative_h2 = h2_output$negative_h2$mean,
      negative_h2_se = h2_output$negative_h2$se,
      stringsAsFactors = FALSE
    )
    
    # Add annotation-specific heritability estimates
    if (!is.null(h2_output$annot_h2$mean) && length(h2_output$annot_h2$mean) > 0) {
      feature_names <- names(h2_output$annot_h2$mean)
      if (is.null(feature_names)) {
        feature_names <- paste0("feature", seq_along(h2_output$annot_h2$mean))
      }
      
      # Create columns for each feature's h2 and se
      for (i in seq_along(feature_names)) {
        feat <- feature_names[i]
        result[[paste0(feat, "_h2")]] <- h2_output$annot_h2$mean[i]
        result[[paste0(feat, "_h2_se")]] <- h2_output$annot_h2$se[i]
      }
    }
    
    result
  }, error = function(e) {
    message(sprintf("  Error processing %s: %s", model_file, e$message))
    NULL
  })
}

# Process all trait-dataset combinations
results <- map_dfr(trait_datasets, ~ {
  if (.x$is_aou) {
    with(.x, process_trait(name, data_name, pheno, is_aou = TRUE, ukbb_pheno = ukbb_pheno))
  } else {
    with(.x, process_trait(name, data_name, pheno, is_aou = FALSE))
  }
})

# Save results
output_file <- sprintf("tables/heritability_estimates_%s.csv", annotation)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.csv(results, file = output_file, row.names = FALSE)
message(sprintf("\nResults saved to: %s", output_file))
