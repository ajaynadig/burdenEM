#!/usr/bin/env Rscript

# sim_cli.R: Calculate true heritability or polygenicity for simulated data

# Library imports
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
library(argparse)
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
if (!requireNamespace("tools", quietly = TRUE)) install.packages("tools")
library(tools)

# Source functions for calculating true heritability and polygenicity
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(script_dir) || script_dir == ".") {
  script_dir <- getwd()
}

source(file.path(script_dir, "sim_heritability.R"))
source(file.path(script_dir, "..", "luke", "heritability.R"))

# Main function to orchestrate the process
main <- function(args) {
  command <- args$command
  studies_filepath <- args$studies_file

  if (!file.exists(studies_filepath)) {
    stop(paste("Studies file not found:", studies_filepath))
  }

  studies_df <- fread(studies_filepath)
  
  if (!("sumstats_filename_pattern" %in% names(studies_df))) {
    stop("'sumstats_filename_pattern' column not found in studies file.")
  }

  # --- Output file setup ---
  studies_file_basename <- basename(studies_filepath)
  studies_file_prefix <- file_path_sans_ext(studies_file_basename)
  studies_file_prefix <- sub("\\.studies$", "", studies_file_prefix, ignore.case = TRUE)
  
  output_dir <- file.path(dirname(studies_filepath), "tables")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- Command-specific logic ---
  if (command == "polygenicity") {
    # --- 1. Calculate TRUE polygenicity ---
    cat("--- Calculating TRUE polygenicity ---\n")
    source(file.path(script_dir, "sim_polygenicity.R"))
    true_results_list <- list()
    for (i in 1:nrow(studies_df)) {
      study_row <- studies_df[i, ]
      sumstats_pattern <- study_row$sumstats_filename_pattern
      genes_filename <- gsub("sumstats", "genes", sumstats_pattern, ignore.case = TRUE)
      if (!startsWith(genes_filename, "/") && !startsWith(genes_filename, "~")){
          genes_filename <- file.path(dirname(studies_filepath), genes_filename)
      }
      if (!file.exists(genes_filename)) {
        warning(paste("Genes file not found for study row", i, ":", genes_filename, "- skipping."))
        next
      }
      genes_df_current <- fread(genes_filename)
      cat(sprintf("Processing study %d/%d: %s, using genes file: %s\n", i, nrow(studies_df), study_row$identifier, genes_filename))
      current_result <- calculate_true_polygenicity(genes_df_current)
      if (!is.null(current_result)) {
          current_result$abbreviation <- study_row$abbreviation
          true_results_list[[length(true_results_list) + 1]] <- current_result
      }
    }

    if (length(true_results_list) == 0) {
      cat("No TRUE polygenicity results generated. Exiting.\n")
      return()
    }
    true_polygenicity_df <- bind_rows(true_results_list)
    if ("abbreviation" %in% names(true_polygenicity_df)) {
      true_polygenicity_df <- true_polygenicity_df %>% dplyr::select(abbreviation, dplyr::everything())
    }
    true_output_filename <- sprintf("%s.true.polygenicity.tsv", studies_file_prefix)
    true_output_filepath <- file.path(output_dir, true_output_filename)
    fwrite(true_polygenicity_df, true_output_filepath, sep = "\t")
    cat(paste("TRUE polygenicity results saved to:", true_output_filepath, "\n"))

    # --- 2. Calculate ESTIMATED polygenicity ---
    source(file.path(script_dir, "..", "luke", "polygenicity.R")) # For calculate_polygenicity_metrics
    cat("\n--- Calculating ESTIMATED polygenicity ---\n")
    estimated_polygenicity_df <- calculate_polygenicity_metrics(studies_df, annotation = args$annotation, verbose = TRUE, run_sequentially_param = args$no_parallel)

    if (nrow(estimated_polygenicity_df) > 0) {
      if ("abbreviation" %in% names(estimated_polygenicity_df)) {
        estimated_polygenicity_df <- estimated_polygenicity_df %>% dplyr::select(abbreviation, dplyr::everything())
      }
      estimated_output_filename <- sprintf("%s.estimated.polygenicity.tsv", studies_file_prefix)
      estimated_output_filepath <- file.path(output_dir, estimated_output_filename)
      fwrite(estimated_polygenicity_df, estimated_output_filepath, sep = "\t")
      cat(paste("ESTIMATED polygenicity results saved to:", estimated_output_filepath, "\n"))
    } else {
      cat("No ESTIMATED polygenicity results generated to save separately.\n")
    }

    # --- 3. Meta-analyze and save --- 
    cat("\n--- Meta-analyzing polygenicity results ---\n")
    meta_polygenicity_df <- meta_analyze_polygenicity(true_polygenicity_df, estimated_polygenicity_df)
    
    meta_output_filename <- sprintf("%s.meta.polygenicity.tsv", studies_file_prefix)
    meta_output_filepath <- file.path(output_dir, meta_output_filename)
    if (is.null(meta_polygenicity_df) || nrow(meta_polygenicity_df) == 0) {
      cat("WARNING: meta_polygenicity_df is NULL or empty. File will not be written.\n")
    } else {
      fwrite(meta_polygenicity_df, meta_output_filepath, sep = "\t")
    }
    cat(paste("Meta-analyzed polygenicity results saved to:", meta_output_filepath, "\n"))

  } else if (command == "heritability") {
    # --- 1. Calculate TRUE heritability ---
    cat("--- Calculating TRUE heritability ---\n")
    true_results_list <- list()
    for (i in 1:nrow(studies_df)) {
      study_row <- studies_df[i, ]
      sumstats_pattern <- study_row$sumstats_filename_pattern
      genes_filename <- gsub("sumstats", "genes", sumstats_pattern, ignore.case = TRUE)
      if (!startsWith(genes_filename, "/") && !startsWith(genes_filename, "~")){
          genes_filename <- file.path(dirname(studies_filepath), genes_filename)
      }
      if (!file.exists(genes_filename)) {
        warning(paste("Genes file not found for study row", i, ":", genes_filename, "- skipping."))
        next
      }
      genes_df <- fread(genes_filename)
      cat(sprintf("Processing study %d/%d: %s, using genes file: %s\n", i, nrow(studies_df), study_row$identifier, genes_filename))
      current_result <- calculate_true_heritability(genes_df)
      current_result$abbreviation <- study_row$abbreviation
      true_results_list[[length(true_results_list) + 1]] <- current_result
    }
    
    if (length(true_results_list) == 0) {
      cat("No TRUE heritability results generated. Exiting.\n")
      return()
    }

    true_heritability_df <- bind_rows(true_results_list)
    if ("abbreviation" %in% names(true_heritability_df)) {
      true_heritability_df <- true_heritability_df %>% dplyr::select(abbreviation, dplyr::everything())
    }
    
    true_output_filename <- sprintf("%s.true.heritability.tsv", studies_file_prefix)
    true_output_filepath <- file.path(output_dir, true_output_filename)
    fwrite(true_heritability_df, true_output_filepath, sep = "\t")
    cat(paste("TRUE heritability results saved to:", true_output_filepath, "\n"))

    # --- 2. Calculate ESTIMATED heritability ---
    cat("\n--- Calculating ESTIMATED heritability ---\n")
    # The 'abbreviation' column is returned by this function, aliased from 'identifier'
    estimated_heritability_df <- calculate_heritability_for_studies(studies_df, annotation_param = args$annotation, verbose_param = TRUE)

    if (nrow(estimated_heritability_df) > 0) {
      if ("abbreviation" %in% names(estimated_heritability_df)) {
        estimated_heritability_df <- estimated_heritability_df %>% dplyr::select(abbreviation, dplyr::everything())
      }
      estimated_output_filename <- sprintf("%s.estimated.heritability.tsv", studies_file_prefix)
      estimated_output_filepath <- file.path(output_dir, estimated_output_filename)
      fwrite(estimated_heritability_df, estimated_output_filepath, sep = "\t")
      cat(paste("ESTIMATED heritability results saved to:", estimated_output_filepath, "\n"))
    } else {
      cat("No ESTIMATED heritability results generated to save separately.\n")
    }
    
    # --- 3. Meta-analyze and save ---
    cat("\n--- Meta-analyzing results ---\n")
    meta_results_df <- meta_analyze_heritability(true_heritability_df, estimated_heritability_df)
    
    meta_output_filename <- sprintf("%s.meta.heritability.tsv", studies_file_prefix)
    meta_output_filepath <- file.path(output_dir, meta_output_filename)
    fwrite(meta_results_df, meta_output_filepath, sep = "\t")
    cat(paste("Meta-analyzed results saved to:", meta_output_filepath, "\n"))
  } else if (command == "distribution") {
    # --- 1. Calculate TRUE distribution ---
    cat("--- Calculating TRUE distribution ---\n")
    source(file.path(script_dir, "sim_distribution.R"))
    true_results_list <- list()
    for (i in 1:nrow(studies_df)) {
      study_row <- studies_df[i, ]
      sumstats_pattern <- study_row$sumstats_filename_pattern
      genes_filename <- gsub("sumstats", "genes", sumstats_pattern, ignore.case = TRUE)
      if (!startsWith(genes_filename, "/") && !startsWith(genes_filename, "~")){
          genes_filename <- file.path(dirname(studies_filepath), genes_filename)
      }
      if (!file.exists(genes_filename)) {
        warning(paste("Genes file not found for study row", i, ":", genes_filename, "- skipping."))
        next
      }
      genes_df <- fread(genes_filename)
      cat(sprintf("Processing study %d/%d: %s, using genes file: %s\n", i, nrow(studies_df), study_row$identifier, genes_filename))
      current_result <- calculate_true_distribution(genes_df)
      current_result$abbreviation <- study_row$abbreviation
      true_results_list[[length(true_results_list) + 1]] <- current_result
    }
    
    if (length(true_results_list) == 0) {
      cat("No TRUE distribution results generated. Exiting.\n")
      return()
    }

    true_distribution_df <- bind_rows(true_results_list)
    if ("abbreviation" %in% names(true_distribution_df)) {
      true_distribution_df <- true_distribution_df %>% dplyr::select(abbreviation, dplyr::everything())
    }
    
    true_output_filename <- sprintf("%s.true.distribution.tsv", studies_file_prefix)
    true_output_filepath <- file.path(output_dir, true_output_filename)
    fwrite(true_distribution_df, true_output_filepath, sep = "\t")
    cat(paste("TRUE distribution results saved to:", true_output_filepath, "\n"))

    # --- 2. Calculate ESTIMATED distribution ---
    cat("\n--- Calculating ESTIMATED distribution ---\n")
    source(file.path(script_dir, "..", "luke", "distribution.R"))
    estimated_distribution_df <- calculate_distribution_metrics_for_studies(
      studies_df = studies_df, 
      annotation = args$annotation,
      verbose = TRUE,
      include_hp_for_1_gene = FALSE,
      run_sequentially_param = args$no_parallel
    )

    if (nrow(estimated_distribution_df) > 0) {
      if ("abbreviation" %in% names(estimated_distribution_df)) {
        estimated_distribution_df <- estimated_distribution_df %>% dplyr::select(abbreviation, dplyr::everything())
      }
      estimated_output_filename <- sprintf("%s.estimated.distribution.tsv", studies_file_prefix)
      estimated_output_filepath <- file.path(output_dir, estimated_output_filename)
      fwrite(estimated_distribution_df, estimated_output_filepath, sep = "\t")
      cat(paste("ESTIMATED distribution results saved to:", estimated_output_filepath, "\n"))
    } else {
      cat("No ESTIMATED distribution results generated to save separately.\n")
    }
    
    # --- 3. Meta-analyze and save ---
    cat("\n--- Meta-analyzing results ---\n")
    meta_results_df <- meta_analyze_distribution(true_distribution_df, estimated_distribution_df)
    
    meta_output_filename <- sprintf("%s.meta.distribution.tsv", studies_file_prefix)
    meta_output_filepath <- file.path(output_dir, meta_output_filename)
    fwrite(meta_results_df, meta_output_filepath, sep = "\t")
    cat(paste("Meta-analyzed results saved to:", meta_output_filepath, "\n"))
  } else if (command == "calibration") {
    # --- Calculate calibration metrics ---
    cat("--- Calculating calibration metrics ---\n")
    source(file.path(script_dir, "sim_calibration.R"))

    calibration_df <- calculate_calibration_metrics_for_studies(
      studies_df = studies_df,
      annotation_param = args$annotation,
      verbose_param = TRUE,
      run_sequentially_param = args$no_parallel,
      per_allele_effects = args$per_allele_effects
    )

    if (nrow(calibration_df) > 0) {
      if ("abbreviation" %in% names(calibration_df)) {
        calibration_df <- calibration_df %>% dplyr::select(abbreviation, dplyr::everything())
      }
      calib_output_filename <- sprintf("%s.calibration.tsv", studies_file_prefix)
      calib_output_filepath <- file.path(output_dir, calib_output_filename)
      fwrite(calibration_df, calib_output_filepath, sep = "\t")
      cat(paste("Calibration results saved to:", calib_output_filepath, "\n"))

      # --- Meta-analyze and save ---
      cat("\n--- Meta-analyzing calibration results ---\n")
      meta_calibration_df <- meta_analyze_calibration(calibration_df)
      meta_output_filename <- sprintf("%s.meta.calibration.tsv", studies_file_prefix)
      meta_output_filepath <- file.path(output_dir, meta_output_filename)
      fwrite(meta_calibration_df, meta_output_filepath, sep = "\t")
      cat(paste("Meta-analyzed calibration results saved to:", meta_output_filepath, "\n"))
    } else {
      cat("No calibration results generated.\n")
    }
  }
}

# Argument parsing
parser <- ArgumentParser(description = "Calculate true heritability or polygenicity from simulated gene data.")
subparsers <- parser$add_subparsers(title = "commands", dest = "command")

parser_heritability <- subparsers$add_parser("heritability", help = "Calculate true and estimated heritability, and meta-analyze.")
parser_heritability$add_argument("studies_file", type = "character", help = "Path to the studies TSV file.")
parser_heritability$add_argument("-a", "--annotation", type = "character", required = TRUE, help = "Annotation to use for estimated heritability calculation (e.g., pLoF).")

parser_polygenicity <- subparsers$add_parser("polygenicity", help = "Calculate true polygenicity metrics.")
parser_polygenicity$add_argument("studies_file", type = "character", help = "Path to the studies TSV file.")
parser_polygenicity$add_argument("-a", "--annotation", type = "character", required = TRUE, help = "Annotation to use for estimated polygenicity calculation (e.g., pLoF).")
parser_polygenicity$add_argument("--no_parallel", action="store_true", default=FALSE, help="Disable parallel execution across studies (runs sequentially).")

parser_distribution <- subparsers$add_parser("distribution", help = "Calculate true and estimated distribution, and meta-analyze.")
parser_distribution$add_argument("studies_file", type = "character", help = "Path to the studies TSV file.")
parser_distribution$add_argument("-a", "--annotation", type = "character", required = TRUE, help = "Annotation to use for estimated distribution calculation (e.g., pLoF).")
parser_distribution$add_argument("--no_parallel", action="store_true", default=FALSE, help="Disable parallel execution across studies (runs sequentially).")

# Calibration subcommand
parser_calibration <- subparsers$add_parser("calibration", help = "Calculate calibration of posterior-mean effects.")
parser_calibration$add_argument("studies_file", type = "character", help = "Path to the studies TSV file.")
parser_calibration$add_argument("-a", "--annotation", type = "character", required = TRUE, help = "Annotation used when loading fitted models (e.g., pLoF).")
parser_calibration$add_argument("--no_parallel", action="store_true", default=FALSE, help="Disable parallel execution across studies (runs sequentially).")
parser_calibration$add_argument("--per_allele_effects", action="store_true", default=FALSE, help="Convert effect sizes to per-allele (e.g., beta -> beta/variant_variance).")

args <- parser$parse_args()

# Run main function
main(args)
