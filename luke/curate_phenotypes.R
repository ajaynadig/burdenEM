#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

ukb_phenotypes <- c("100017_NA", "4080_NA",  "4079_NA",  "21001_NA",  "21002_NA",  "30040_NA",  "30060_NA",  "30070_NA",  "30080_NA",  "30100_NA",
                    "30190_NA", "30200_NA", "30220_NA", "30530_NA", "30610_NA", "30650_NA", "30670_NA", "30690_NA", "30740_NA", "30760_NA", "30860_NA",
                    "30870_NA", "48_NA", "49_NA", "50_NA", "5983_NA", "WHR_custom_NA", "30730_NA", "30620_NA", "30780_NA",
                    '20002_asthma', '20002_hypertension', '20002_diabetes', '130708_NA', '131306_NA', '130706_NA', '20002_hypertension', "20002_osteoarthritis",
                    "20002_osteoarthritis", '130708','20002_diabetes', '131706_NA', '131306_NA', '20002_asthma')
aou_phenotypes <- c('3001420', 'blood-pressure-systolic-mean', 'blood-pressure-diastolic-mean', 'BMI', 'weight', '3023599', '3009744', '3019897', '3024929', '3043111',
                    '3011948', '3008342', '3013869', '3019550', '3035995', '3013721', '3013682', '3027114', '3004501', '3007070', '3020630',
                    '3022192', 'waist-circumference-mean', 'hip-circumference-mean', 'height', 'heart-rate-mean', 'WHR', "3026910", "3006923", "3028288",
                    '495', '401', '250', '250.2', '411', '250.1', 'CV_401', 'MS_708', '740', 'EM_202.2', 'EM_202', 'EM_202.1', 'CV_404', 'RE_475')

# Define input and output file paths
input_csv_path <- path.expand("~/Dropbox/burdenEM_results/allxall_ALL_phenotype_info_250k_new.csv")
output_tsv_path <- "data/shared.studies.tsv"

# Ensure output directory exists
dir.create(dirname(output_tsv_path), showWarnings = FALSE, recursive = TRUE)

# Read the input phenotype CSV
phenotype_info_raw <- read_csv(input_csv_path, col_types = cols(.default = "c"))

# Filter for specified aou_phenotypes
phenotype_info_filtered <- phenotype_info_raw %>%
  dplyr::filter(phenoname %in% aou_phenotypes)

# Create a mapping from AoU phenotype codes to UKB phenotype codes
# Ensure aou_phenotypes and ukb_phenotypes have the same length and correspond to each other
if (length(aou_phenotypes) != length(ukb_phenotypes)) {
  stop("aou_phenotypes and ukb_phenotypes vectors must have the same length for mapping.")
}
pheno_map <- data.frame(aou_code = aou_phenotypes, ukb_code = ukb_phenotypes, stringsAsFactors = FALSE)

# Join filtered data with the mapping to get ukb_code
phenotype_info_mapped <- phenotype_info_filtered %>%
  left_join(pheno_map, by = c("phenoname" = "aou_code"))

# Define target datasets
target_datasets <- c('genebass', 'aou_afr', 'aou_amr', 'aou_eas', 'aou_sas', 'aou_eur')

# Process the phenotype information
studies_data <- phenotype_info_mapped %>%
  # Create rows for each target dataset
  crossing(dataset = target_datasets) %>%
  mutate(
    identifier = description, # Original 'description' from CSV (e.g., 'Height')
    trait_type = trait_type, # Original 'trait_type' from CSV
    # dataset column is already from crossing
    abbreviation = description, # Original 'description' from CSV, used as abbreviation
    description = description_more, # Original 'description_more' from CSV
    sumstats_filename_pattern = case_when(
      dataset == 'genebass' ~ sprintf("data/var_txt/^genebass_%s_<ANNOTATION>_.*\\.txt\\.bgz$", ukb_code),
      startsWith(dataset, 'aou_') ~ sprintf("data/var_txt/^%s_%s_<ANNOTATION>_.*\\.txt\\.bgz$", dataset, phenoname),
      TRUE ~ NA_character_
    ),
    model_filename = case_when(
      dataset == 'genebass' ~ sprintf("fitted_models/burdenEM_fit_genebass_%s_<ANNOTATION>.rds", ukb_code),
      startsWith(dataset, 'aou_') ~ sprintf("fitted_models/burdenEM_fit_%s_%s_<ANNOTATION>.rds", dataset, phenoname),
      TRUE ~ NA_character_ # Should not happen
    )
  ) %>%
  # Select and order final columns for studies.tsv
  select(
    identifier,
    trait_type,
    dataset,
    description,
    abbreviation,
    sumstats_filename_pattern,
    model_filename
  ) %>%
  distinct() # Ensure unique rows if any input phenotype might lead to identical entries after processing

# Write the output to studies.tsv
write_tsv(studies_data, output_tsv_path)

message(sprintf("Successfully generated %s with %d rows.", output_tsv_path, nrow(studies_data)))
