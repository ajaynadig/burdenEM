#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

genebass_info <- read_delim('data/phenotype/genebass_pheno_results.txt.bgz')
phecode_info <- read_csv('data/phenotype/phecodeX_info.csv')
aou_info <- read_tsv('data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv')  %>%
  filter(ancestry == 'META' & pheno_sex == 'both') %>%
  merge(., phecode_info %>% select(phenoname=phecode), by = 'phenoname', all.x=T) %>%
  select(phenoname, category, disease_category, n_cases_aou = n_cases, n_controls_aou = n_controls) %>%
  distinct()

full_mapped <- read_tsv('data/phenotype/aou_ukb_matched_all_phenotype_v8.csv') %>%
  merge(., genebass_info %>% mutate(n_controls_ukb = if_else(is.na(n_controls), 0, n_controls)) %>% select(ukb_phenocode=phenocode, n_cases_ukb = n_cases, n_controls_ukb), by = 'ukb_phenocode' ) %>%
  merge(., aou_info, by = 'phenoname' )  %>%
  mutate(prevalence_aou = n_cases_aou/(n_cases_aou + n_controls_aou),
         prevalence_ukb = n_cases_ukb/(n_cases_ukb + n_controls_ukb),
         N_aou = n_cases_aou + n_controls_aou,
         N_ukb = n_cases_ukb + n_controls_ukb,
         ) %>%
  mutate(selected = N_aou > 200000 & N_ukb > 200000 & prevalence_aou > 0.01 & prevalence_ukb > 0.01 | (phenoname %in% c('3028288', 'EM_202.2'))) %>%
  mutate(selected_for_display = (selected & N_aou > 210000& prevalence_aou > 0.1 & prevalence_ukb > 0.1) &
           (phenoname %in% c('weight', 'height', 'BMI', 'blood-pressure-diastolic-mean', 'WHR') | category != 'physical_measurement') | (phenoname %in% c('3028288', 'EM_202.2')) )
table(full_mapped$selected)
table(full_mapped$selected_for_display)
write_tsv(full_mapped, 'data/all_mapped_phenotypes.tsv')
sub_mapped <- full_mapped %>%
  filter(selected)
subsub_mapped  <- full_mapped %>%
  filter(selected_for_display)
table(subsub_mapped$category)

ukb_phenotypes <- c(sub_mapped$ukb_phenocode)
aou_phenotypes <- c(sub_mapped$phenoname)
# ukb_phenotypes <- c(subsub_mapped$ukb_phenocode)
# aou_phenotypes <- c(subsub_mapped$phenoname)


# Define input and output file paths
input_csv_path <- path.expand('data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv')
output_tsv_path <- "data/shared.full.selected.studies.tsv"

# Ensure output directory exists
dir.create(dirname(output_tsv_path), showWarnings = FALSE, recursive = TRUE)

# Read the input phenotype CSV
phenotype_info_raw <- read_tsv(input_csv_path)

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
# target_datasets <- c('genebass', 'aou_afr', 'aou_amr', 'aou_eas', 'aou_sas', 'aou_eur')
target_datasets <- c('genebass', 'aou_afr', 'aou_amr', 'aou_eur', 'aou_meta', 'biobank_meta')

# Process the phenotype information
studies_data <- phenotype_info_mapped %>%
  # Create rows for each target dataset
  crossing(dataset = target_datasets) %>%
  mutate(
    identifier = description, # Original 'description' from CSV (e.g., 'Height')
    trait_type = trait_type, # Original 'trait_type' from CSV
    # dataset column is already from crossing
    abbreviation = description, # Original 'description' from CSV, used as abbreviation
    # description = description_more, # Original 'description_more' from CSV
    sumstats_filename_pattern = case_when(
      dataset == 'genebass' ~ sprintf("data/var_txt/^genebass_%s_<ANNOTATION>_.*\\.txt\\.bgz$", ukb_code),
      startsWith(dataset, 'aou_') ~ sprintf("data/var_txt/^%s_%s_<ANNOTATION>_.*\\.txt\\.bgz$", dataset, phenoname),
      startsWith(dataset, 'biobank_') ~ sprintf("data/var_txt/^%s_%s_<ANNOTATION>_.*\\.txt\\.bgz$", dataset, phenoname),
      TRUE ~ NA_character_
    ),
    model_filename = case_when(
      dataset == 'genebass' ~ sprintf("fitted_models/burdenEM_fit_genebass_%s_<ANNOTATION>.rds", ukb_code),
      startsWith(dataset, 'aou_') ~ sprintf("fitted_models/burdenEM_fit_%s_%s_<ANNOTATION>.rds", dataset, phenoname),
      startsWith(dataset, 'biobank_') ~ sprintf("fitted_models/burdenEM_fit_%s_%s_<ANNOTATION>.rds", dataset, phenoname),
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
