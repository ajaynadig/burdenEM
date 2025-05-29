# main_cli.R
# Command-line interface for running burdenEM
# This file contains only the command-line interface code that was previously in main.R

# Source the main functions without running the CLI
source("luke/main.R")

# Only run the CLI code if this script is executed directly
if (!interactive() && identical(commandArgs(trailingOnly = FALSE)[1], "--file=luke/main_cli.R")) {
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
        make_option(c("--intercept_frequency_bin_edges"), type="character", default="0,1e-5,1e-4,1e-3",
            help="Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3')", metavar="character")
    )

    opt_parser = OptionParser(option_list=option_list)
    opt = parse_args(opt_parser)

    # Check required arguments
    if (is.null(opt$pheno) || is.null(opt$annotation_to_process)) {
        print_help(opt_parser)
        stop("Phenotype (--pheno) and annotation (--annotation_to_process) must be specified.", call.=FALSE)
    }

    # Parse frequency ranges
    frequency_range <- as.numeric(strsplit(opt$frequency_range, ",")[[1]])
    intercept_frequency_bin_edges <- as.numeric(strsplit(opt$intercept_frequency_bin_edges, ",")[[1]])

    # Create output directory if it doesn't exist
    if (!dir.exists(opt$output_dir)) {
        dir.create(opt$output_dir, recursive = TRUE)
    }

    # Prepare arguments for run_burdenEM_rvas
    run_args <- list(
        variant_dir = opt$variant_dir,
        ld_corrected_scores_file = opt$ld_scores_file,
        data_name = opt$data_name,
        pheno = opt$pheno,
        annotation_to_process = opt$annotation_to_process,
        feature_col_name = opt$feature_col_name,
        num_feature_bins = opt$num_feature_bins,
        burdenem_no_cpts = opt$burdenem_no_cpts,
        num_iter = opt$num_iter,
        per_allele_effects = opt$per_allele_effects,
        correct_for_ld = opt$correct_for_ld,
        correct_genomewide = opt$correct_genomewide,
        verbose = opt$verbose,
        frequency_range = frequency_range,
        intercept_frequency_bin_edges = intercept_frequency_bin_edges,
        output_dir = opt$output_dir
    )

    # Remove NULL values from the arguments
    run_args <- run_args[!sapply(run_args, is.null)]

    # Run the main function
    do.call(run_burdenEM_rvas, run_args)
}
