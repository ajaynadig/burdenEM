# main_cli.R
# Command-line interface for running burdenEM
# This file contains only the command-line interface code that was previously in main.R

# Source the main functions without running the CLI
source("luke/main.R")

# Only run the CLI code if this script is executed directly
message("Parsing Command Line Arguments")

# Check/install optparse
if (!requireNamespace("optparse", quietly = TRUE)) {
    message("Installing 'optparse' package...")
    install.packages("optparse", repos = "http://cran.us.r-project.org")
}
library(optparse)

# Check/install future
if (!requireNamespace("future", quietly = TRUE)) {
    message("Installing 'future' package...")
    install.packages("future", repos = "http://cran.us.r-project.org")
}
library(future)

# Check/install furrr
if (!requireNamespace("furrr", quietly = TRUE)) {
    message("Installing 'furrr' package...")
    install.packages("furrr", repos = "http://cran.us.r-project.org")
}
library(furrr)

# Check/install purrr (for sequential map, though often a dependency of furrr/dplyr)
if (!requireNamespace("purrr", quietly = TRUE)) {
    message("Installing 'purrr' package...")
    install.packages("purrr", repos = "http://cran.us.r-project.org")
}
library(purrr)

library(readr) # For reading the studies TSV
library(stringr) # For string manipulations

option_list = list(
    make_option(c("-a", "--annotation"), type="character", default=NULL, 
        help="Single functional annotation category to process (e.g., 'pLoF', 'missense', 'synonymous') [required]", metavar="character"),
    make_option(c("-f", "--feature_col_name"), type="character", default=NULL, 
        help="Optional: Name of the column in the genes file for gene features (e.g., 'oe_lof')", metavar="character"),
    make_option(c("-b", "--num_feature_bins"), type="integer", default=5, 
        help="Number of bins for feature column [default %default].", metavar="integer"),
    make_option(c("-g", "--genes_file"), type="character", default=NULL, 
        help="Path to the genes file template (use <ANNOTATION>, <LOWER>, and <UPPER> placeholders) [required].", metavar="character"),
    make_option(c("-c", "--num_positive_components"), type="integer", default=10, 
        help="Number of positive components for BurdenEM [default %default].", metavar="integer"),
    make_option(c("-i", "--num_iter"), type="integer", default=10000, 
        help="Number of EM iterations [default %default].", metavar="integer"),
    make_option(c("--per_allele_effects"), action="store_true", default=FALSE, dest="per_allele_effects",
        help="Calculate per-allele effect sizes instead of per-gene. [default: %default]"),
    make_option(c("--correct_for_ld"), action="store_true", default=FALSE, dest="correct_for_ld",
        help="Apply LD correction to burden scores and gamma."),
    make_option(c("--frequency_range"), type="character", default="0,0.001", help="Comma-separated min,max for allele frequency range (e.g., '0,0.001'). Default: %default"),
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print extra details during execution."),
    make_option(c("--no_parallel"), action="store_true", default=FALSE, help="Disable parallel execution across studies (runs sequentially)."),
    make_option(c("-n", "--name"), type="character", default=NULL, 
        help="Name for the output directory. Overrides the directory in the studies file.", metavar="character"),
    make_option(c("--intercept_frequency_bin_edges"), type="character", default="0,1e-5,1e-4,1e-3",
        help="Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3')", metavar="character")
)

opt_parser = OptionParser(option_list=option_list, usage = "%prog [options] studies_file")
opt_parsed = parse_args(opt_parser, positional_arguments = TRUE)
opt = opt_parsed$options

# Check required arguments
if (length(opt_parsed$args) == 0 || is.null(opt_parsed$args[1])) {
    print_help(opt_parser)
    stop("A studies_file positional argument must be provided.", call.=FALSE)
}
studies_file_path <- opt_parsed$args[1]

if (is.null(opt$annotation) || is.null(opt$genes_file)) {
    print_help(opt_parser)
    stop("Annotation (--annotation), and genes file (--genes_file) must be specified.", call.=FALSE)
}

# Parse frequency ranges from CLI - these will be applied to all studies
frequency_range_cli <- as.numeric(strsplit(opt$frequency_range, ",")[[1]])
intercept_frequency_bin_edges_cli <- as.numeric(strsplit(opt$intercept_frequency_bin_edges, ",")[[1]])

# Read the studies file
studies_df <- readr::read_tsv(studies_file_path, show_col_types = FALSE)

# Worker function to process a single study
process_study_cli <- function(study_row, opt_config, freq_range_cli, icept_freq_bins_cli) {
    current_study <- study_row

    # --- 1. Sumstats path handling ---
    sumstats_pattern_template <- current_study$sumstats_filename_pattern
    full_sumstats_pattern_with_anno <- stringr::str_replace(sumstats_pattern_template, "<ANNOTATION>", opt_config$annotation)
    
    derived_variant_dir <- dirname(full_sumstats_pattern_with_anno)
    derived_variant_file_pattern <- basename(full_sumstats_pattern_with_anno)

    # --- 2. Genes file handling ---
    genes_file_template <- opt_config$genes_file
    processed_genes_file <- stringr::str_replace(genes_file_template, "<DATASET>", current_study$dataset)

    # --- 3. Model output path handling ---
    model_output_path_template <- current_study$model_filename
    model_output_path_with_anno <- stringr::str_replace(model_output_path_template, "<ANNOTATION>", opt_config$annotation)

    if(opt_config$verbose) {
        message(paste("Initial model_output_path_template:", model_output_path_template))
        message(paste("Initial model_output_path_with_anno:", model_output_path_with_anno))
    }

    if (!is.null(opt_config$name)) {
        studies_dir <- dirname(studies_file_path)
        output_dir <- file.path(studies_dir, opt_config$name)
        model_output_dir <- file.path(output_dir, "models")
        model_filename_base <- paste0(current_study$identifier, ".", opt_config$annotation, ".rds")
        model_output_path_with_anno <- file.path(model_output_dir, model_filename_base)
        if(opt_config$verbose) {
            message(paste("(--name) studies_dir:", studies_dir))
            message(paste("(--name) model_output_dir:", model_output_dir))
            message(paste("(--name) final model_output_path_with_anno:", model_output_path_with_anno))
        }
    }
    
    model_output_dir <- dirname(model_output_path_with_anno)
    if (!dir.exists(model_output_dir)) {
        dir.create(model_output_dir, recursive = TRUE)
        if(opt_config$verbose) message(paste("Created model output directory:", model_output_dir, "for study:", current_study$identifier))
    }
    output_file_prefix_for_run <- tools::file_path_sans_ext(model_output_path_with_anno)
    if(opt_config$verbose) {
        message(paste("Final output_file_prefix_for_run:", output_file_prefix_for_run))
    }

    # --- 4. Prepare arguments for run_burdenEM_rvas for the current study ---
    run_args <- list(
        variant_dir = derived_variant_dir,
        variant_file_pattern = derived_variant_file_pattern,
        ld_corrected_scores_file = processed_genes_file,
        output_file_prefix = output_file_prefix_for_run,
        annotation_to_process = opt_config$annotation,
        feature_col_name = opt_config$feature_col_name,
        num_feature_bins = opt_config$num_feature_bins,
        num_positive_components = opt_config$num_positive_components,
        num_iter = opt_config$num_iter,
        per_allele_effects = opt_config$per_allele_effects,
        correct_for_ld = opt_config$correct_for_ld,
        frequency_range = freq_range_cli, 
        intercept_frequency_bin_edges = icept_freq_bins_cli,
        verbose = opt_config$verbose
    )

    run_args <- run_args[!sapply(run_args, is.null)]

    main_message <- paste0("Processing: ", current_study$identifier, " (", current_study$dataset, ") [", opt_config$annotation, "]")
    if(opt_config$verbose) {
        main_message <- paste0("Running BurdenEM for: ", current_study$identifier, " (", current_study$dataset, ") with annotation: ", opt_config$annotation, 
                           ". Output prefix: ", output_file_prefix_for_run)
    }
    message(main_message)
    do.call(run_burdenEM_rvas, run_args)
    # tryCatch({
    #     do.call(run_burdenEM_rvas, run_args)
    #     if(opt_config$verbose || !opt$no_parallel) { # Also print completion if parallel to confirm it finished
    #          message(paste0("Completed: ", current_study$identifier, " (", current_study$dataset, ")"))
    #     }
    # }, error = function(e) {
    #     message(paste0("ERROR processing study ", current_study$identifier, " (", current_study$dataset, "): ", e$message))
    # })
    return(NULL) 
}

# Main loop to process each study defined in the TSV file
message(paste("Processing", nrow(studies_df), "studies from:", studies_file_path))

if (!opt$no_parallel && nrow(studies_df) > 1) {
    # Determine number of workers for future_map
    # future::availableCores() gives logical cores, future::availableCores(logical=FALSE) gives physical
    # Using availableCores() (logical) is often a good default for multisession
    num_workers <- future::availableCores()
    if (is.na(num_workers) || num_workers <= 1) {
        message("Only one core detected or available. Running sequentially.")
        # Fallback to sequential execution
        purrr::map(1:nrow(studies_df), 
                   ~process_study_cli(studies_df[.x, ], opt, frequency_range_cli, intercept_frequency_bin_edges_cli))
    } else {
        message(paste("Setting up parallel execution using up to", num_workers, "workers for", nrow(studies_df), "studies."))
        future::plan(future::multisession, workers = num_workers) # Explicitly set workers
        
        # furrr::future_map iterates over elements, so we pass indices or rows
        # Passing indices and subsetting inside .f is common
        results_list <- furrr::future_map(
            .x = 1:nrow(studies_df), 
            .f = ~process_study_cli(studies_df[.x, ], opt, frequency_range_cli, intercept_frequency_bin_edges_cli),
            .progress = opt$verbose, # Show progress bar if verbose
            .options = furrr_options(seed = TRUE) # Ensure reproducibility if any RNG is used in worker
        )
        message("Parallel processing finished.")
        # Reset future plan to default (sequential) after use, if desired, or leave as is for potential subsequent parallel operations
        # future::plan(future::sequential) 
    }
} else {
    if (nrow(studies_df) <= 1 && !opt$no_parallel) message("Only one study to process, running sequentially.")
    else message("Running in sequential mode (--no_parallel specified or only one study/core).")
    
    purrr::map(1:nrow(studies_df), 
               ~process_study_cli(studies_df[.x, ], opt, frequency_range_cli, intercept_frequency_bin_edges_cli))
}

message("All studies processed.")
