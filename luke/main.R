# main.R
# Runs burdenEM inputting variant-level sumstats and gene-level burden scores/gene features; saves gene-level sumstats and a model object
# For help:
# Rscript main.R -h

# Paths used to set default values in command line arguments
TOP_LEVEL_DATA_DIR <- "/Users/lukeoconnor/Dropbox/burdenEM_results/"
OUTPUT_DIR <- "./burdenEM_results"

# --- Libraries ---
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(optparse)
options(vroom.show_problems = FALSE)

# --- Source Required Functions ---
source("R/io.R")
source("luke/intercept.R")
source("luke/variant_to_gene.R")
source("R/model.R")
source("R/likelihoods.R")
source("R/EM.R") # Source for EM_fit, bootstrap_EM, null_EM_rvas
source("luke/gene_features.R") # Source for get_bins
source("luke/burdenEM.R") # Source for fit_burdenem_model
source("luke/negative_binomial.R") # Source for estimate_overdispersion

#' Process Variant Data to Gene-Level for Poisson Model (Binary Traits)
#'
#' Aggregates variant-level data to gene-level, calculating
#' combined allele count in cases (CAC_cases) and combined allele frequency (CAF).
#' It also brings N (total sample size) and prevalence to the gene level.
#'
#' @param variant_data A data frame from `load_variant_files_with_category`,
#'   expected to contain `gene`, `functional_category`, `AC_cases`, `AF` (allele frequency),
#'   `N` (sample size), and `prevalence`.
#' @param verbose Logical, if TRUE, prints progress messages.
#' @return A data frame aggregated by `gene` and `functional_category`, with
#'   columns `gene`, `functional_category`, `CAC_cases`, `CAF`, `N`, `prevalence`, and `n_variants`.
process_variant_to_gene_poisson <- function(variant_data, verbose = FALSE) {
  # Input validation
  required_cols <- c("gene", "functional_category", "AC_cases", "AF", "N", "prevalence")
  missing_cols <- setdiff(required_cols, names(variant_data))
  if (length(missing_cols) > 0) {
    stop(paste("variant_data is missing one or more required columns for Poisson processing:",
               paste(missing_cols, collapse = ", ")))
  }

  if (verbose) message("Aggregating variants to gene-level for Poisson model...")

  gene_level_data <- variant_data %>%
    dplyr::group_by(gene, functional_category) %>%
    dplyr::summarise(
      CAC_cases = sum(AC_cases),
      CAF = sum(AF),
      N = dplyr::first(N), # Assumes N is constant for the group
      prevalence = dplyr::first(prevalence), # Assumes prevalence is constant and in the data
      n_variants = dplyr::n(),
      burden_score = sum(2 * AF * (1 - AF)),
      CAC_expected = 2 * prevalence * N * CAF, 
      .groups = 'drop'
    ) %>%
    dplyr::filter(!is.na(gene))

  if (verbose) message(paste("Aggregated to", nrow(gene_level_data), "gene-category rows for Poisson model."))

  return(gene_level_data)
}

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
    variant_dir = file.path(TOP_LEVEL_DATA_DIR, "data", "genebass", "var_txt"),
    ld_corrected_scores_file = NULL,
    data_name = "genebass",
    pheno = "50_NA",
    annotation_to_process,
    intercept_frequency_bin_edges = c(0, 1e-5,1e-4, 1e-3),
    frequency_range = c(0, 0.001),
    feature_col_name = "lof.oe",
    num_feature_bins = 5,
    burdenem_no_cpts = 10, # TODO
    burdenem_grid_size = 100,
    num_iter = 1000,
    bootstrap = TRUE,
    n_boot = 100,
    output_dir = OUTPUT_DIR,
    verbose = FALSE,
    per_allele_effects = FALSE,
    correct_for_ld = FALSE,
    correct_genomewide = FALSE
) {

    # --- 1. Load Variant-Level Data ---
    if(verbose) message("\n--- Loading Variant-Level Data ---")
    variant_data <- load_variant_files_with_category(
        variant_dir = variant_dir,
        data_name = data_name,
        pheno = pheno,
        annotations_to_process = c(annotation_to_process),
        frequency_range = frequency_range
    )
    if(verbose) message(paste("Successfully loaded", nrow(variant_data), "variants."))

    # --- Detect Trait Type --- 
    required_binary_cols <- c("AC_cases", "N")
    is_binary <- all(required_binary_cols %in% names(variant_data))
    trait_type <- if(is_binary) "binary" else "continuous"
    message(paste("Detected trait type:", trait_type))

    # --- 2. Load LD-Corrected Burden Scores ---
    if(verbose) message("\n--- Loading LD-Corrected Burden Scores ---")
    lower_fmt <- format(frequency_range[1], nsmall=1, scientific=FALSE)
    upper_fmt <- format(frequency_range[2], nsmall=1, scientific=FALSE)
    dynamic_ld_scores_file <- ld_corrected_scores_file %>%
        gsub("<ANNOTATION>", annotation_to_process, ., fixed = TRUE) %>%
        gsub("<LOWER>", lower_fmt, ., fixed = TRUE) %>%
        gsub("<UPPER>", upper_fmt, ., fixed = TRUE)
    ld_corrected_scores_df <- data.table::fread(dynamic_ld_scores_file) %>%
        dplyr::rename(burden_score_ld = burden_score)
    if(verbose) message(paste("Loaded", nrow(ld_corrected_scores_df), 
        "LD-corrected scores from:", dynamic_ld_scores_file))

    # --- 3. Process Variants to Gene-Level --- 
    if (verbose) message(paste("\n--- Processing Variants to Gene-Level (", trait_type, ") ---"))
    if (is_binary) {
        gene_level_data <- process_variant_to_gene_poisson(
            variant_data = variant_data
        )
    } else {
        gene_level_data <- process_variant_to_gene(
            variant_data = variant_data, 
            frequency_bin_edges = intercept_frequency_bin_edges
        ) %>% 
            dplyr::filter(!is.na(gene))
    }

    if (nrow(gene_level_data) == 0) {
        stop("No gene-level data remained after processing and filtering. Check input data and filters.")
    }
    message(paste("Aggregated to", nrow(gene_level_data), "gene-category rows."))

    # --- 4. Append LD-Corrected Scores ---
    if(verbose) message("\n--- Appending LD-Corrected Burden Scores ---")
    merge_col <- if (grepl("^ENSG", gene_level_data$gene[1])) "gene_id" else "gene"
    gene_level_data <- gene_level_data %>%
        left_join(ld_corrected_scores_df, by = c("gene" = merge_col, "functional_category" = "functional_category"))

    # --- 5. Specify Effect Sizes (Continuous Only) & Prepare for Model --- 
    if(trait_type == "continuous"){
        if(verbose) message("\n--- Specifying Gene Effect Sizes (Continuous) ---")
        # Pass per_allele_effects and correct_for_ld from function arguments
        message(paste("Applying per-allele effects:", per_allele_effects))
        gene_level_data <- specify_gene_effect_sizes(gene_level_data, 
                                                   per_allele_effects = per_allele_effects, 
                                                   correct_for_ld = correct_for_ld)
    } else {
        if(verbose) message("\n--- Skipping Effect Size Specification (Binary) ---")
        # For binary, the relevant columns (CAC_cases, CAF, N, prevalence) are already present
    }

    # --- 6. Prepare Gene Features & Add to Gene Data --- 
    # This step uses the output 'gene_level_data' from step 3/4, independent of trait type
    # It adds the 'features' column based on 'feature_col_name'
    if (!is.null(feature_col_name)) {
        if (verbose) message(paste("\n--- Generating", num_feature_bins, "feature bins for:", feature_col_name, "---"))
        if (!(feature_col_name %in% names(gene_level_data))){
             stop(paste("Specified feature_col_name '", feature_col_name, "' not found in gene_level_data after aggregation.", sep=""))
        }
        # Ensure get_bins is available (e.g., source("luke/gene_features.R"))
        if (!exists("get_bins")) stop("Function 'get_bins' not found. Source luke/gene_features.R")
        gene_level_data <- get_bins(gene_level_data, column_name = feature_col_name, num_bins = num_feature_bins, verbose = verbose)
    } else {
        # Assign a single feature bin if no feature column is specified
        gene_level_data$features <- replicate(nrow(gene_level_data), c(1), simplify = FALSE)
         if (verbose) message("\n--- No feature column specified, assigning single feature bin ---")
    }
    
    # --- 7. Define Likelihood Function ---
    if(trait_type == "continuous") {
        message("Using Normal likelihood for continuous trait.")
        likelihood_function <- function(row, beta_vec) {
            dnorm(row$effect_estimate, mean = beta_vec, sd = row$effect_se)
        }
    } else { # Binary
        message("Using Poisson likelihood for binary trait.")
        likelihood_function <- function(row, beta_vec) {
            lambda <- row$CAC_expected * exp(beta_vec)
            dpois(x = row$CAC_cases, lambda = lambda)
        }
    }

    # Define h2 function based on per_allele_effects
    if (trait_type == "continuous") {
    if (per_allele_effects){
        if (!("burden_score" %in% names(gene_level_data))){
             stop("burden_score column required in gene_level_data when per_allele_effects is TRUE")
        }
        h2_function <- function(beta, row) row$burden_score * beta^2
    }
    else {
       h2_function <- function(beta, row) beta^2
    }
    }
    else {
        h2_function <- function(beta, row) {
        (row$prevalence * (exp(beta)-1)^2 * 2*row$CAF)/(1-row$prevalence)}
    }
    
    # --- 8. Fit BurdenEM Model --- 
    if (verbose) message("\n--- Fitting BurdenEM Model ---")
    burdenem_model <- fit_burdenem_model(
        gene_level_data = gene_level_data,
        burdenem_grid_size = burdenem_grid_size,
        num_iter = num_iter,
        per_allele_effects = per_allele_effects, # Keep for continuous endpoint scaling
        verbose = verbose,
        likelihood_function = likelihood_function,
        h2_function = h2_function
    )

    print(burdenem_model$component_endpoints)
    print(burdenem_model$delta)

    # --- 9. Save Model --- 
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
        make_option(c("-p","--pheno"), type="character", default=NULL, 
            help="Phenotype name (e.g., '50_NA', '20002_diabetes') [required]", metavar="character"),
        make_option(c("-a", "--annotation_to_process"), type="character", default=NULL, 
            help="Single functional annotation category to process (e.g., 'pLoF', 'missense', 'synonymous') [required]", metavar="character"),
        make_option(c("-f", "--feature_col_name"), type="character", default=NULL, 
            help="Optional: Name of the column in LD-scores file for gene features (e.g., 'oe_lof')", metavar="character"),
        make_option(c("-b", "--num_feature_bins"), type="integer", default=5, 
            help="Number of bins for feature column [default %default].", metavar="integer"),
        make_option(c("-v", "--variant_dir"), type="character", 
            default=file.path(TOP_LEVEL_DATA_DIR, "data", "genebass", "var_txt"), 
            help="Directory containing variant files [default %default]", metavar="character"),
        make_option(c("-l", "--ld_scores_file"), type="character", 
            default=file.path(TOP_LEVEL_DATA_DIR, "data", "utility", "ukbb_ld_corrected_burden_scores_<ANNOTATION>_<LOWER>_<UPPER>.tsv"),
            help="Path to LD-corrected burden scores file template (use <ANNOTATION> placeholder) [default %default].", metavar="character"),
        make_option(c("-o", "--output_dir"), type="character", default="./burdenEM_example_output", 
            help="Output directory for results [default %default].", metavar="character"),
        make_option(c("-d", "--data_name"), type="character", default="genebass", 
            help="Dataset name prefix (e.g., 'genebass') [default %default].", metavar="character"),
        make_option(c("-c", "--burdenem_no_cpts"), type="integer", default=10, 
            help="Number of components for BurdenEM [default %default].", metavar="integer"),
        make_option(c("-i", "--num_iter"), type="integer", default=1000, 
            help="Number of EM iterations [default %default].", metavar="integer"),
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
            help="Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3')", metavar="character")
    ); 
    
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$pheno)){
        print_help(opt_parser)
        stop("Phenotype name (--pheno) must be specified.", call.=FALSE)
    }
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
        pheno = opt$pheno,
        annotation_to_process = opt$annotation_to_process,
        intercept_frequency_bin_edges = as.numeric(strsplit(opt$intercept_frequency_bin_edges, ",")[[1]]),
        frequency_range = as.numeric(strsplit(opt$frequency_range, ",")[[1]]),
        feature_col_name = opt$feature_col_name,
        num_feature_bins = opt$num_feature_bins,
        burdenem_no_cpts = opt$burdenem_no_cpts, # Note: burdenem_no_cpts not used directly now
        burdenem_grid_size = 100,
        num_iter = opt$num_iter,
        bootstrap = opt$bootstrap,
        n_boot = opt$bootstrap_n,
        output_dir = opt$output_dir,
        verbose = opt$verbose,
        per_allele_effects = opt$per_allele_effects,
        correct_for_ld = opt$correct_for_ld,
        correct_genomewide = opt$correct_genomewide
    )

    # Call the main function with dynamically prepared arguments
    fitted_model <- do.call(run_burdenEM_rvas, run_args)

    if (!is.null(fitted_model)) {
        message("\n--- Workflow Finished Successfully ---")
    } else {
         message("\n--- Workflow Finished (potentially with errors or no model returned) ---")
    }
}