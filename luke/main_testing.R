# main.R
# Runs burdenEM inputting variant-level sumstats and gene-level burden scores/gene features; saves gene-level sumstats and a model object
# For help:
# Rscript main.R -h

# Paths used to set default values in command line arguments
TOP_LEVEL_DATA_DIR = "/Users/lukeoconnor/Dropbox/burdenEM_results/"
OUTPUT_DIR = "./burdenEM_results"

# --- Libraries ---
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(optparse)

# --- Source Required Functions ---
# Assuming R/* scripts are in the R/ directory relative to the project root
source("R/io.R")
source("luke/intercept.R")
source("luke/variant_to_gene.R")
source("luke/genome_wide_burden.R")
source("R/model.R")
source("R/likelihoods.R")
source("R/EM.R") # Source for EM_fit, bootstrap_EM, null_EM_rvas
source("luke/gene_features.R") # Source for get_bins
source("luke/negative_binomial.R") # Source for estimate_overdispersion

#' Run BurdenEM RVAS Workflow
#'
#' Loads variant data, processes it to gene-level, initializes the BurdenEM model,
#' fits the model using EM, and optionally performs bootstrapping and null simulations.
#'
#' @param variant_dir Directory containing variant files.
#' @param data_name Dataset name prefix (e.g., 'genebass').
#' @param pheno Phenotype name (e.g., '50_NA').
#' @param annotation_to_process Character string for the single functional category to load.
#' @param intercept_frequency_bin_edges Numeric vector defining AF bin edges, only for intercept calculation.
#' @param feature_col_name Optional string: name of the column in the LD-corrected scores file to use for gene features (e.g., "lof_oe"). Default: NULL.
#' @param num_feature_bins Integer: number of quantile bins for the feature column. Default: 5.
#' @param burdenem_no_cpts Number of components for BurdenEM.
#' @param burdenem_grid_size Grid size for likelihood calculation.
#' @param num_iter Number of iterations for EM algorithm.
#' @param bootstrap Logical, whether to perform bootstrapping.
#' @param n_boot Number of bootstrap samples if bootstrap=TRUE.
#' @param output_dir Directory to save the fitted model.
#' @param verbose Logical, whether to print detailed messages.
#' @param frequency_range Numeric vector c(min, max) for AF range filter [min,max). Default: c(0, 0.001).
#' @param per_allele_effects Logical, whether to calculate per-allele effect sizes instead of per-gene. Default: FALSE.
#' @param correct_for_ld Logical, whether to apply LD correction to burden scores and gamma. Default: FALSE.
#' @param correct_genomewide Logical, whether to apply genome-wide burden correction. Default: FALSE.
#' @param estimate_overdispersion Logical, whether to estimate overdispersion. Default: TRUE.
#' @return A list object representing the fitted BurdenEM model.
#' @export
run_burdenEM_rvas <- function(
    variant_dir = file.path(TOP_LEVEL_DATA_DIR, "genebass", "var_txt"),
    ld_corrected_scores_file = file.path(TOP_LEVEL_DATA_DIR, "utility", "ukbb_ld_corrected_burden_scores_<ANNOTATION>.tsv"),
    data_name = "genebass",
    pheno = "50_NA",
    annotation_to_process, 
    intercept_frequency_bin_edges = c(0, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3),
    frequency_range = c(0, 0.001),
    feature_col_name = "oe_lof",
    num_feature_bins = 5,
    burdenem_no_cpts = 11,
    burdenem_grid_size = 101,
    num_iter = 100000,
    bootstrap = TRUE,
    n_boot = 100,
    output_dir = OUTPUT_DIR,
    verbose = FALSE,
    per_allele_effects = FALSE,
    correct_for_ld = FALSE,
    correct_genomewide = FALSE,
    skip_overdispersion = TRUE
) {

    if(verbose) message("--- Starting BurdenEM RVAS Workflow ---")


    # --- 1. Load Variant-Level Data ---
    if(verbose) message("\n--- Loading Variant-Level Data ---")

    # Load data for the single specified annotation
    variant_data <- load_variant_files_with_category(
        variant_dir = variant_dir,
        data_name = data_name,
        pheno = pheno,
        annotations_to_process = c(annotation_to_process), # Corrected: Pass single annotation as vector
        frequency_range = frequency_range, # Pass the frequency range
        prevalence = pheno_prevalence # Pass the loaded prevalence
    )
    if(verbose) message(paste("Successfully loaded", nrow(variant_data), "variants."))
    trait_is_binary = !is.null(variant_data$AC_cases)

    # Ensure this script only runs for binary traits currently
    if (!trait_is_binary) {
        stop("This script currently only supports binary/categorical traits.")
    }

    # TODO: Implement variant LD pruning based on correct_for_ld flag

    # --- 2. Load LD-Corrected Burden Scores ---
    if(verbose) message("\n--- Loading LD-Corrected Burden Scores ---")
    
    # Construct filename dynamically based on the annotation being processed
    dynamic_ld_scores_file <- file.path(TOP_LEVEL_DATA_DIR, "data", "utility", 
        sprintf("ukbb_ld_corrected_burden_scores.tsv_%s.tsv", annotation_to_process))
    
    ld_corrected_scores_df <- data.table::fread(dynamic_ld_scores_file)%>% 
        dplyr::rename(burden_score_ld = burden_score)
    if(verbose) message(paste("Loaded", nrow(ld_corrected_scores_df), 
        "LD-corrected scores from:", dynamic_ld_scores_file))

    # --- 3. Process Variants to Gene-Level (including intercept calculation, aggregation, estimate/SE calc, filtering) ---
    # Estimate overdispersion first for binary traits
    if (!skip_overdispersion) {
        if (verbose) message("\n---Estimating overdispersion parameter ---")
        variant_data <- estimate_overdispersion(
            variant_data = variant_data,
            intercept_frequency_bin_edges = intercept_frequency_bin_edges
        )
        if(verbose) message("Overdispersion estimated:")
        print(variant_data$overdispersion[0])
    } else {
        variant_data$overdispersion <- 0
    }

    # --- 3b. Binary Trait: Calculate and Add Variant Likelihoods ---
    # Define the grid for beta values
    # TODO: Make this grid definition more flexible / parameterizable
    if (verbose) message("\n--- Processing Variants to Gene-Level ---")
    beta_grid <- seq(-3, 3, length.out = burdenem_grid_size) # Example grid
    if(verbose) message(paste("Calculating and adding variant likelihoods on a grid of size", length(beta_grid)))

    # Call the new function which returns the dataframe with likelihoods added
    variant_data_with_ll <- add_nbinom_likelihood(
        variant_df = variant_data, # Pass the data with overdispersion
        grid = beta_grid
    )
    if(verbose) message("Variant likelihood columns added.")

    # --- 3c. Binary Trait: Aggregate to Gene Level ---
    # We now already have the combined dataframe needed by the updated function
    gene_level_results <- process_variant_to_gene_binary(
        variant_data_with_ll = variant_data_with_ll # Pass the dataframe with likelihoods
    )
    gene_level_stats <- gene_level_results$gene_level_stats
    gene_level_likelihoods <- gene_level_results$gene_level_likelihoods
    if(verbose) message("Likelihoods aggregated to gene level.")

    # Placeholder for downstream steps - need to adapt model init/fit
    # Note: 'gene_level_data' is assigned stats for now, but we'll need likelihoods for model
    gene_level_data <- gene_level_stats
    message("Binary trait processing complete up to likelihood aggregation.")
    # TODO: Adapt model initialization and fitting steps for binary traits using gene_level_likelihoods


    # --- 7. Save Gene-Level Data --- 
    gene_data_filename <- paste0("gene_level_data_", data_name, "_", pheno, "_", annotation_to_process, ".rds")
    gene_data_path <- file.path(output_dir, gene_data_filename)
    # TODO: Clarify/standardize exactly what columns are saved here
    if(verbose) message(paste("\n--- Saved Gene-Level Data to:", gene_data_path, "---"))
    saveRDS(gene_level_data, file = gene_data_path)
    if(verbose) message("Gene-level data saved successfully.")
    # TODO: Save gene-level statistics (burden scores, features, posteriors?)

    # --- 8. Initialize BurdenEM Model --- 
    if(verbose) message("\n--- Initializing BurdenEM Model ---")

    # Always initialize for binary trait now
    message("Initializing model for binary trait...")

    # Define the endpoints for the mixture components
    component_endpoints <- seq(-3, 3, length.out = burdenem_no_cpts)

    # Initialize model using the precomputed likelihoods and specified endpoints
    model <- initialize_model_from_likelihood(
      grid_loglikelihoods = gene_level_likelihoods,
      beta_grid = beta_grid,
      component_endpoints = component_endpoints,
      features = NULL # TODO: Load and pass actual gene features matrix here
    )

    if(verbose) {
        message(">>> DEBUG: Initial gene assignment based on max conditional likelihood:")
        # Find the index of the component with the maximum likelihood for each gene
        max_likelihood_component_index <- apply(model$conditional_likelihood, 1, which.max)
        # Create a table of counts for each component index
        component_counts <- table(factor(max_likelihood_component_index, levels = 1:length(model$component_endpoints)))
        # Map indices back to endpoint values for clarity
        names(component_counts) <- paste0("EndPt_", round(model$component_endpoints, 3))
        print(component_counts)
    }
    if(verbose) message("BurdenEM model initialized.")

    # --- 9. Fit Model using EM --- 
    if(verbose) message("\n--- Running EM Fit ---")
    model <- EM_fit(model = model, max_iter = num_iter, tol = 0)
    if(verbose) {
        message("EM fit complete. Final delta parameters:")
        message("Component Endpoints:")
        print(round(model$component_endpoints, 3))
        # Assign endpoint names to delta columns for clarity
        delta_to_print <- model$delta
        colnames(delta_to_print) <- paste0("EndPt_", round(model$component_endpoints, 3))
        print(round(delta_to_print, 5))
    }

    
    # --- 10. Bootstrap EM (Optional) ---
    if (bootstrap) {
        if(verbose) message(paste("\n--- Running Bootstrap EM (n_boot =", n_boot, ") ---"))
        model$bootstrap_seeds = 1:n_boot
        bootstrap_output <- bootstrap_EM(
            model = model,
            n_boot = n_boot,
            max_iter = num_iter,
            bootstrap_seeds = model$bootstrap_seeds
        )
        model$bootstrap_output = bootstrap_output$bootstrap_delta
        if(verbose) message("Bootstrap EM complete.")
    }

    # --- Save Model --- 
    # Generate dynamic filename
    output_filename_gen <- paste0("burdenEM_fit_", data_name, "_", pheno, "_", annotation_to_process, ".rds")
    output_path <- file.path(output_dir, output_filename_gen)
    if(verbose) message(paste("\n--- Saving Fitted Model to:", output_path, "---"))
    saveRDS(model, file = output_path)
    if(verbose) message("Model saved successfully.")
    return(model)
}

if (!interactive()) {
    
    message("--- Parsing Command Line Arguments ---")

    # Check/install optparse
    if (!requireNamespace("optparse", quietly = TRUE)) {
        message("Installing 'optparse' package...")
        install.packages("optparse", repos = "http://cran.us.r-project.org")
    }
    library(optparse)

    option_list = list(
        make_option(c("-a", "--annotation_to_process"), type="character", default=NULL, 
            help="Single functional annotation category to process (e.g., 'pLoF', 'missense', 'synonymous') [required]", metavar="character"),
        make_option(c("-f", "--feature_col_name"), type="character", default=NULL, 
            help="Optional: Name of the column in LD-scores file for gene features (e.g., 'oe_lof')", metavar="character"),
        make_option(c("-b", "--num_feature_bins"), type="integer", default=5, 
            help="Number of bins for feature column [default %default].", metavar="integer"),
        make_option(c("-v", "--variant_dir"), type="character", 
            default=file.path(TOP_LEVEL_DATA_DIR, "data", "genebass", "var_txt"), 
            help="Directory containing variant files [default %default].", metavar="character"),
        make_option(c("-l", "--ld_scores_file"), type="character", 
            default=file.path(TOP_LEVEL_DATA_DIR, "data","utility","ukbb_ld_corrected_burden_scores.tsv"), 
            help="Path to LD-corrected burden scores file [default %default].", metavar="character"),
        make_option(c("-o", "--output_dir"), type="character", default="./burdenEM_example_output", 
            help="Output directory for results [default %default].", metavar="character"),
        make_option(c("-p", "--phenotype"), type="character", default="50_NA", 
            help="Phenotype code (default: 50_NA -> height).", metavar="character"),
        make_option(c("-d", "--data_name"), type="character", default="genebass", 
            help="Dataset name prefix (e.g., 'genebass') [default %default].", metavar="character"),
        make_option(c("-c", "--burdenem_no_cpts"), type="integer", default=11, 
            help="Number of components for BurdenEM [default %default].", metavar="integer"),
        make_option(c("-i", "--num_iter"), type="integer", default=1000, 
            help="Number of EM iterations [default %default].", metavar="integer"),
        make_option(c("--bootstrap_n"), type="integer", default=100, 
            help="Number of bootstrap samples [default %default].", metavar="integer"),
        make_option(c("--no_bootstrap"), action="store_false", default=TRUE, dest="bootstrap",
            help="Disable bootstrapping."),
        make_option(c("--per_allele_effects"), action="store_true", default=FALSE, dest="per_allele_effects",
            help="Calculate per-allele effect sizes instead of per-gene. [default: %default]"),
        make_option(c("--correct_for_ld"), action="store_true", default=FALSE, dest="correct_for_ld",
            help="Apply LD correction to burden scores and gamma."),
        make_option(c("--correct_genomewide"), action="store_true", default=FALSE, 
                    help="Apply genome-wide burden correction after variant-to-gene processing [default %default]"),
        make_option(c("-q", "--quiet"), action="store_false", default=TRUE, dest="verbose",
            help="Suppress verbose messages."),
        make_option(c("--frequency_range"), type="character", default="0,0.001",
            help="Comma-separated string for AF range filter [min,max) (e.g., '0,0.001')", metavar="character"),
        make_option(c("--intercept_frequency_bin_edges"), type="character", default="0,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3", 
            help="Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3')", metavar="character"),
        make_option(c("-g", "--burdenem_grid_size"), type="integer", default=101, 
            help="Grid size for likelihood calculation [default %default].", metavar="integer"),
        make_option(c("--skip_overdispersion"), action="store_true", default=FALSE,
            help="Skip overdispersion estimation and set it to 0 [default: %default]")
    ); 
    
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$annotation_to_process)){
        print_help(opt_parser)
        stop("Annotation type (--annotation_to_process) must be specified.", call.=FALSE)
    }

    # Ensure output directory exists
    if (!dir.exists(opt$output_dir)) {
        message(paste("Creating output directory:", opt$output_dir))
        dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    message("\n--- Running BurdenEM RVAS Workflow from Command Line ---")
    
    # Prepare arguments for the main function
    run_args <- list(
        variant_dir = opt$variant_dir,
        ld_corrected_scores_file = opt$ld_scores_file,
        data_name = opt$data_name,
        pheno = opt$phenotype,
        annotation_to_process = opt$annotation_to_process,
        intercept_frequency_bin_edges = as.numeric(strsplit(opt$intercept_frequency_bin_edges, ",")[[1]]),
        frequency_range = as.numeric(strsplit(opt$frequency_range, ",")[[1]]), # Parse frequency range
        feature_col_name = opt$feature_col_name,
        num_feature_bins = opt$num_feature_bins,
        burdenem_no_cpts = opt$burdenem_no_cpts,
        burdenem_grid_size = opt$burdenem_grid_size, # Pass grid size from CLI args
        num_iter = opt$num_iter,
        bootstrap = opt$bootstrap,
        n_boot = opt$bootstrap_n,
        output_dir = opt$output_dir,
        verbose = opt$verbose,
        per_allele_effects = opt$per_allele_effects,
        correct_for_ld = opt$correct_for_ld,
        correct_genomewide = opt$correct_genomewide,
        skip_overdispersion = opt$skip_overdispersion # Control estimation via flag
    )

    # Call the main function with dynamically prepared arguments
    fitted_model <- do.call(run_burdenEM_rvas, run_args)

    if (!is.null(fitted_model)) {
        message("\n--- Workflow Finished Successfully ---")
    } else {
         message("\n--- Workflow Finished (potentially with errors or no model returned) ---")
    }
}