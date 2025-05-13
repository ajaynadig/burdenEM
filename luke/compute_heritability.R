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
annotation <- ifelse(length(args) > 1, args[1], "pLoF")
phenotype <- ifelse(length(args) > 1, args[2], "50_NA")
per_allele_effects <- ifelse(length(args) > 2, args[3], FALSE)
data_name <- ifelse(length(args) > 3, args[4], "genebass")

message(paste("Calculating heritability for", annotation, "annotation and", phenotype, "phenotype, data_name:", data_name, "..."))

model_file <- paste0("burdenEM_example_output/burdenEM_fit_", data_name, "_",phenotype,"_", annotation, ".rds")
print(model_file)
model <- readRDS(model_file)

h2_output <- estimate_heritability_components(model, verbose = FALSE)
message("burdenEM total heritability estimate: ", round(h2_output$total_h2$mean, 6), " (SE: ", round(h2_output$total_h2$se, 6), ")")
message("burdenEM positive-effect heritability estimate: ", round(h2_output$positive_h2$mean, 6), " (SE: ", round(h2_output$positive_h2$se, 6), ")")
message("burdenEM negative-effect heritability estimate: ", round(h2_output$negative_h2$mean, 6), " (SE: ", round(h2_output$negative_h2$se, 6), ")")

alt_h2 <- posterior_expectation2(model, model$h2_function)
argmax_index <- which.max(alt_h2)
# map(model$df[argmax_index, ], ~ print(.x))

# (mean chisq - 1) * appropriate constant
message("\nMethod-of-moments estimator:")
multiplier <- if (per_allele_effects) {
  mean(model$df$effect_se^2 * model$df$burden_score)
} else {
  mean(model$df$effect_se^2)
}
implied_h2_se <- sqrt(3 * nrow(model$df)) * multiplier
implied_h2 <- sum(model$df$effect_estimate^2/model$df$effect_se^2 - 1) * multiplier
message(paste("  h2:", round(implied_h2, 6)))
message(paste("  h2 SE:", round(implied_h2_se, 7)))
