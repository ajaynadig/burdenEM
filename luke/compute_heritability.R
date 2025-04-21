#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

source("R/estimate_heritability.R")

# --- Get Annotation from Command Line ---
args <- commandArgs(trailingOnly = TRUE)
annotation <- args[1]
phenotype <- ifelse(length(args) > 1, args[2], "50_NA")
per_allele_effects <- ifelse(length(args) > 2, args[3], FALSE)
message(paste("Calculating heritability for", annotation, "annotation and", phenotype, "phenotype..."))

gene_file <- paste0("burdenEM_example_output/gene_level_data_genebass_",phenotype,"_", annotation, ".rds")
gene_data <- readRDS(gene_file)

model_file <- paste0("burdenEM_example_output/burdenEM_fit_genebass_",phenotype,"_", annotation, ".rds")
model <- readRDS(model_file)

h2_output <- estimate_heritability_rvas(
  model = model, 
  genetic_data = gene_data, 
  per_allele_effects = per_allele_effects
  )
message("burdenEM total heritability estimate: ", round(h2_output$total_h2, 6))

# (mean chisq - 1) * appropriate constant
message("\nMethod-of-moments estimator:")
multiplier <- if (per_allele_effects) {
  mean(gene_data$effect_se^2 * gene_data$burden_score)
} else {
  mean(gene_data$effect_se^2)
}
implied_h2_se <- sqrt(3 * nrow(gene_data)) * multiplier
implied_h2 <- sum(gene_data$effect_estimate^2/gene_data$effect_se^2 - 1) * multiplier
message(paste("  h2:", round(implied_h2, 6)))
message(paste("  h2 SE:", round(implied_h2_se, 7)))


