# main.R
# Runs burdenEM inputting variant-level sumstats and gene-level burden scores/gene features; saves gene-level sumstats and a model object
# For help:
# Rscript main.R -h

# Paths used to set default values in command line arguments
TOP_LEVEL_DATA_DIR = "~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/data"
OUTPUT_DIR = "~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/results_2025/"

# --- Libraries ---
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(optparse)

# --- Source Required Functions ---
# Assuming R/* scripts are in the R/ directory relative to the project root
setwd('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/')
source("R/io.R")
source("luke/intercept.R")
source("luke/variant_to_gene.R")
source("luke/genome_wide_burden.R")
source("R/model.R")
source("R/likelihoods.R")
source("R/EM.R") # Source for EM_fit, bootstrap_EM, null_EM_rvas
source("luke/gene_features.R") # Source for get_bins

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
#' @return A list object representing the fitted BurdenEM model.
#' @export
run_burdenEM_rvas <- function(
    variant_dir = file.path(TOP_LEVEL_DATA_DIR, "genebass", "var_txt"),
    ld_corrected_scores_file = file.path(TOP_LEVEL_DATA_DIR, "utility", "ukbb_ld_corrected_burden_scores_<ANNOTATION>.tsv"),
    data_name = "genebass",
    pheno = "50_NA",
    annotation_to_process = 'pLoF',
    intercept_frequency_bin_edges = c(0, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3),
    frequency_range = c(0, 0.001),
    feature_col_name = "lof.oe",
    num_feature_bins = 5,
    burdenem_no_cpts = 15,
    burdenem_grid_size = 100,
    num_iter = 100000,
    bootstrap = TRUE,
    n_boot = 100,
    output_dir = OUTPUT_DIR,
    verbose = FALSE,
    per_allele_effects = FALSE,
    customize_components = FALSE,
    correct_for_ld = TRUE,
    correct_genomewide = FALSE
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
        frequency_range = frequency_range # Pass the frequency range
    )
    if(verbose) message(paste("Successfully loaded", nrow(variant_data), "variants."))

    # --- 2. Load LD-Corrected Burden Scores ---
    if(correct_for_ld){
      if(verbose) message("\n--- Loading LD-Corrected Burden Scores ---")

      # Construct filename dynamically based on the annotation being processed
      dynamic_ld_scores_file <- file.path(TOP_LEVEL_DATA_DIR, "utility", "burdenscore",
                                          sprintf("ukbb_ld_corrected_burden_scores_%s_%s.tsv", pheno, annotation_to_process))

      ld_corrected_scores_df <- data.table::fread(dynamic_ld_scores_file)%>%
        dplyr::rename(burden_score_ld = burden_score)
      if(verbose) message(paste("Loaded", nrow(ld_corrected_scores_df),
                                "LD-corrected scores from:", dynamic_ld_scores_file))
    }


    # --- 3. Process Variants to Gene-Level (including intercept calculation, aggregation, estimate/SE calc, filtering) ---
    if (verbose) message("\n--- Processing Variants to Gene-Level ---")
    gene_level_data <- process_variant_to_gene(
        variant_data = variant_data,
        frequency_bin_edges = intercept_frequency_bin_edges
    )

    # --- 4. Append LD-Corrected Scores ---
    if(correct_for_ld){
      if(verbose) message("\n--- Appending LD-Corrected Burden Scores ---")
      gene_level_data <- gene_level_data %>%
        left_join(ld_corrected_scores_df, by = c("gene", "functional_category"))
      if(verbose) message(paste("LD-corrected scores joined. Resulting columns:", paste(names(gene_level_data), collapse=", ")))
    }

    # --- 4b. Apply Genome-Wide Correction (Optional) ---
    if (correct_genomewide) {
        if (verbose) message("\n--- Applying Genome-Wide Burden Correction ---")
        gene_level_data <- correct_genomewide_burden(gene_level_data)
    }

    # --- 5. Handle LD Correction and put either per-allele or per-s.d. effect sizes as effect_estimate and effect_se ---
    if(verbose) message("\n--- Preparing Data for BurdenEM Model Input ---")
    # Pass per_allele_effects and correct_for_ld from function arguments
    message(paste("Applying per-allele effects:", per_allele_effects))
    gene_level_data <- specify_gene_effect_sizes(gene_level_data,
                                               per_allele_effects = per_allele_effects,
                                               correct_for_ld = correct_for_ld)

    # --- 6. Prepare Gene Features & Add to Gene Data ---
    if (!is.null(feature_col_name)) {
        if (verbose) message(paste("\n--- Generating", num_feature_bins, "feature bins for:", feature_col_name, "---"))
        gene_level_data <- get_bins(gene_level_data, column_name = feature_col_name, num_bins = num_feature_bins, verbose = verbose)
        gene_features_matrix <- do.call(rbind, gene_level_data$features)

    } else {
        if(verbose) message("No feature column specified. Skipping feature generation.")
        gene_features_matrix <- NULL
    }

    # --- 7. Save Gene-Level Data ---
    gene_data_filename <- paste0("gene_level_data_", data_name, "_", pheno, "_", annotation_to_process, ".rds")
    gene_data_path <- file.path(output_dir, gene_data_filename)
    saveRDS(gene_level_data, file = gene_data_path)
    if(verbose) message(paste("\n--- Saved Gene-Level Data to:", gene_data_path, "---"))

    # --- 8. Initialize BurdenEM Model ---
    if(customize_components){
      lower_bound = unlist(quantile(gene_level_data$gamma_per_sd, probs = c(0.001)))
      upper_bound = unlist(quantile(gene_level_data$gamma_per_sd, probs = c(0.999)))
      bound = max(abs(min(gene_level_data$gamma_per_sd)), abs(max(gene_level_data$gamma_per_sd)))
      component_endpoints1 = c(-bound,
                               seq(lower_bound, upper_bound, length.out = (burdenem_no_cpts/2)),
                               bound)
      component_endpoints2 = seq(-bound, bound, length.out = burdenem_no_cpts/2)
      component_endpoints2[which(component_endpoints2 == 0)] = 1e-300
      burdenem_component_endpoints <- sort(unique(c(component_endpoints1, component_endpoints2)))
    }else{
      burdenem_component_endpoints <- c(-1, -0.5, -0.2, -0.1, -0.05, -0.02)
      burdenem_component_endpoints <- c(burdenem_component_endpoints, 0, -rev(burdenem_component_endpoints))
    }
    if(!per_allele_effects){
        burdenem_component_endpoints <- burdenem_component_endpoints * sqrt(mean(gene_level_data$burden_score))
    }
    # burdenem_component_endpoints <- c(-4, -2, -1, -0.5, -0.2, -0.1, -0.05)
    message(paste("Computed component endpoints:"))
    print(burdenem_component_endpoints)

    if(verbose) message("Initializing continuous-trait BurdenEM model...")
    burdenem_likelihood_func <- likelihood_function_rvas(trait_type = "continuous")

    burdenem_model <- initialize_model(
        likelihood_function = burdenem_likelihood_func,
        genetic_data = gene_level_data, # Pass the processed data
        component_endpoints = burdenem_component_endpoints,
        features = gene_features_matrix,
        grid_size = burdenem_grid_size
    )

    # --- 9. Fit Model using EM ---
    if(verbose) message("\n--- Running EM Fit ---")
    burdenem_model <- EM_fit(model = burdenem_model, max_iter = num_iter, tol = 0)
    if(verbose) {
        message("EM fit complete. Final delta parameters (first few rows/cols):")
        print(head(burdenem_model$delta))
    }

    # --- 10. Bootstrap EM (Optional) ---
    if (bootstrap) {
        if(verbose) message(paste("\n--- Running Bootstrap EM (n_boot =", n_boot, ") ---"))
        burdenem_model$bootstrap_seeds = 1:n_boot
        bootstrap_output <- bootstrap_EM(
            model = burdenem_model,
            n_boot = n_boot,
            max_iter = num_iter,
            bootstrap_seeds = burdenem_model$bootstrap_seeds
        )
        burdenem_model$bootstrap_output = bootstrap_output$bootstrap_delta
        if(verbose) message("Bootstrap EM complete.")
    }

    # --- Save Model ---
    # Generate dynamic filename
    output_filename_gen <- paste0("burdenEM_fit_", data_name, "_", pheno, "_", annotation_to_process, ".rds")
    output_path <- file.path(output_dir, output_filename_gen)
    if(verbose) message(paste("\n--- Saving Fitted Model to:", output_path, "---"))
    saveRDS(burdenem_model, file = output_path)
    if(verbose) message("Model saved successfully.")
    return(burdenem_model)
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
        make_option(c("-p", "--phenotype"), type="character", default="50_NA,30780_NA",
            help="Phenotype code split by comma (default: 50_NA -> height).", metavar="character"),
        make_option(c("-d", "--data_name"), type="character", default="genebass",
            help="Dataset name prefix (e.g., 'genebass') [default %default].", metavar="character"),
        make_option(c("-c", "--burdenem_no_cpts"), type="integer", default=10,
            help="Number of components for BurdenEM [default %default].", metavar="integer"),
        make_option(c("-i", "--num_iter"), type="integer", default=1000,
            help="Number of EM iterations [default %default].", metavar="integer"),
        make_option(c("--bootstrap_n"), type="integer", default=100,
            help="Number of bootstrap samples [default %default].", metavar="integer"),
        make_option(c("--no_bootstrap"), action="store_false", default=TRUE, dest="bootstrap",
            help="Disable bootstrapping."),
        make_option(c("--per_allele_effects"), action="store_true", default=FALSE, dest="per_allele_effects",
            help="Calculate per-allele effect sizes instead of per-gene. [default: %default]"),
        make_option(c("--customize_components"), action="store_true", default=FALSE, dest="customize_components",
                    help="Customize selection of component endpoints using gene-level data. [default: %default]"),
        make_option(c("--correct_for_ld"), action="store_true", default=FALSE, dest="correct_for_ld",
            help="Apply LD correction to burden scores and gamma."),
        make_option(c("--correct_genomewide"), action="store_true", default=FALSE,
                    help="Apply genome-wide burden correction after variant-to-gene processing [default %default]"),
        make_option(c("-q", "--quiet"), action="store_false", default=TRUE, dest="verbose",
            help="Suppress verbose messages."),
        make_option(c("--frequency_range"), type="character", default="0,0.001",
            help="Comma-separated string for AF range filter [min,max) (e.g., '0,0.001')", metavar="character"),
        make_option(c("--intercept_frequency_bin_edges"), type="character", default="0,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3",
            help="Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3')", metavar="character")
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

    phenotypes <- str_split_1(opt$phenotype, ',')

    for(pheno in phenotypes){
      # Prepare arguments for the main function
      run_args <- list(
        variant_dir = opt$variant_dir,
        ld_corrected_scores_file = opt$ld_scores_file,
        data_name = opt$data_name,
        pheno = pheno,
        annotation_to_process = opt$annotation_to_process,
        intercept_frequency_bin_edges = as.numeric(strsplit(opt$intercept_frequency_bin_edges, ",")[[1]]),
        frequency_range = as.numeric(strsplit(opt$frequency_range, ",")[[1]]), # Parse frequency range
        feature_col_name = opt$feature_col_name,
        num_feature_bins = opt$num_feature_bins,
        burdenem_no_cpts = opt$burdenem_no_cpts,
        # burdenem_grid_size = 100, # Use default from function
        num_iter = opt$num_iter,
        bootstrap = opt$bootstrap,
        n_boot = opt$bootstrap_n,
        output_dir = opt$output_dir,
        verbose = opt$verbose,
        per_allele_effects = opt$per_allele_effects,
        customize_components = opt$customize_components,
        correct_for_ld = opt$correct_for_ld,
        correct_genomewide = opt$correct_genomewide
      )

      # Call the main function with dynamically prepared arguments
      fitted_model <- do.call(run_burdenEM_rvas, run_args)
    }

    if (!is.null(fitted_model)) {
        message("\n--- Workflow Finished Successfully ---")
    } else {
         message("\n--- Workflow Finished (potentially with errors or no model returned) ---")
    }
}


# Rscript main.R \
# -a <annotation_category> \
# -f <feature_column_name> \
# -b <num_bins> \
# -v <variant_directory_path> \
# -l <ld_scores_file_path> \
# -o <output_directory_path> \
# -p <phenotype_code> \
# -d <dataset_name> \
# -c <num_components> \
# -i <num_iterations> \
# --bootstrap_n <num_bootstraps> \
# --frequency_range <min_freq,max_freq> \
# --intercept_frequency_bin_edges <edge1,edge2,...> \
# # Optional flags (include if needed):
# # --no_bootstrap
# # --per_allele_effects
# # --correct_for_ld
# # --correct_genomewide
# # --quiet
