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
library(dplyr) # For data manipulation
library(data.table) # For fread

option_list = list(
    make_option(c("-a", "--annotation"), type="character", default=NULL,
        help="Single functional annotation category to process (e.g., 'pLoF', 'missense', 'synonymous') [required]", metavar="character"),
    make_option(c("-f", "--feature_col_name"), type="character", default=NULL,
        help="Optional: Name of the column in the genes file for gene features (e.g., 'oe_lof')", metavar="character"),
    make_option(c("-b", "--num_feature_bins"), type="integer", default=5,
        help="Number of bins for feature column [default %default].", metavar="integer"),
    make_option(c("-g", "--genes_file"), type="character", default=NULL,
        help="Path to the genes file template (optionally use <ANNOTATION>, <LOWER>, <UPPER>, and <DATASET> placeholders) [required].", metavar="character"),
    make_option(c("-c", "--num_positive_components"), type="integer", default=10,
        help="Number of positive components for BurdenEM [default %default].", metavar="integer"),
    make_option(c("--burdenem_grid_size"), type="integer", default=4,
        help="Grid size for BurdenEM [default %default].", metavar="integer"),
    make_option(c("-i", "--num_iter"), type="integer", default=5000,
        help="Number of EM iterations [default %default].", metavar="integer"),
    make_option(c("--per_allele_effects"), action="store_true", default=FALSE, dest="per_allele_effects",
        help="Calculate per-allele effect sizes instead of per-gene. [default: %default]"),
    make_option(c("-m", "--binary_trait_model_type"), type="character", default="betabinom", dest="binary_trait_model_type",
        help="Model type for binary traits (betabinom, binom, nbinom, or pois) [default: %default]"),
    make_option(c("--correct_for_ld"), action="store_true", default=FALSE, dest="correct_for_ld",
        help="Apply LD correction to burden scores and gamma."),
    make_option(c("--frequency_range"), type="character", default="0,0.001", help="Comma-separated min,max for allele frequency range (e.g., '0,0.001'). Default: %default"),
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print extra details during execution."),
    make_option(c("--no_parallel"), action="store_true", default=FALSE, help="Disable parallel execution across studies (runs sequentially)."),
    make_option(c("-n", "--name"), type="character", default=NULL,
        help="Name for the output directory. Overrides the directory in the studies file.", metavar="character"),
    make_option(c("--skip_existing"), action="store_true", default=FALSE, help="Skip studies whose output .rds already exists (after applying --name and <ANNOTATION>)."),
    make_option(c("--intercept_frequency_bin_edges"), type="character", default="0,1e-5,1e-4,1e-3",
        help="Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3')", metavar="character"),
    make_option(c("--optimizer"), type="character", default="EM",
        help="Optimization method: 'EM' or 'mixsqp' [default: %default]", metavar="character"),
    make_option(c("--removed_variants_file"), type="character", default=NULL, dest="removed_variants_file",
        help="Path template to removed variants file for LD pruning (binary traits only). Supports <DATASET> placeholder. Accepts local paths or gs:// URIs (auto-downloaded). File should have 'locus' or 'SNP' column.", metavar="character"),
    make_option(c("--gcs_var_prefix"), type="character", default=NULL, dest="gcs_var_prefix",
        help="GCS prefix for variant files (e.g., 'gs://aou_wlu/v8_analysis/BurdenEM'). When set, variant files are downloaded from GCS. Subdirectory structure: <prefix>/<subdir>/var_files/. Default: read from local paths.", metavar="character"),
    make_option(c("--num_workers"), type="integer", default=NULL, dest="num_workers",
        help="Number of parallel workers (caps memory usage). Default: min(4, availableCores()). Ignored if --no_parallel.", metavar="integer"),
    make_option(c("--no_chunk_by_identifier"), action="store_true", default=FALSE, dest="no_chunk_by_identifier",
        help="Disable phenotype-grouped processing. By default, studies are processed one phenotype (identifier) at a time so the GCS cache is reused across the 6 datasets sharing component files, then cleaned up before the next phenotype.")
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

# Sort by identifier so studies sharing variant files are processed adjacently.
# This maximizes GCS cache reuse: the 6 datasets per phenotype (3 AoU components +
# genebass + 2 metas) all reference the same component files.
studies_df <- studies_df %>% dplyr::arrange(identifier, dataset)

# --- Meta dataset component definitions ---
.meta_components <- list(
    aou_meta = c("aou_afr", "aou_amr", "aou_eur"),
    biobank_meta = c("aou_afr", "aou_amr", "aou_eur", "genebass")
)

get_component_datasets <- function(dataset) {
    if (dataset %in% names(.meta_components)) return(.meta_components[[dataset]])
    return(dataset)
}

is_meta_dataset <- function(dataset) {
    dataset %in% names(.meta_components)
}

# Helper: extract ancestry/population suffix from dataset name
# e.g., "aou_afr" -> "afr", "aou_eur" -> "eur", "genebass" -> "genebass"
get_ancestry_from_dataset <- function(dataset) {
    if (dataset %in% names(.meta_components)) {
        stop(paste("Meta dataset label", dataset, "should be expanded to components before calling get_ancestry_from_dataset"))
    }
    if (grepl("genebass", dataset, ignore.case = TRUE)) return("genebass")
    parts <- strsplit(dataset, "_")[[1]]
    return(parts[length(parts)])
}

# --- ENSG → gene symbol mapping (for biobank_meta gene name harmonization) ---
# AoU variant files use ENSG IDs; genebass uses gene symbols. When constructing
# biobank_meta on-the-fly, we need to harmonize AoU gene names to symbols.
.ensg_to_symbol <- NULL
if (any(studies_df$dataset == "biobank_meta")) {
    utility_dir <- "data/utility"
    genes_file <- list.files(utility_dir, pattern = "^genebass_ld_corrected_burden_scores_pLoF_", full.names = TRUE)[1]
    if (!is.na(genes_file)) {
        mapping_df <- data.table::fread(genes_file, select = c("gene", "gene_id"))
        .ensg_to_symbol <- setNames(mapping_df$gene, mapping_df$gene_id)
        message(paste("Loaded ENSG→symbol mapping:", length(.ensg_to_symbol), "genes"))
    } else {
        warning("Cannot find genebass genes file for ENSG→symbol mapping. biobank_meta gene harmonization may fail.")
    }
}

# --- GCS helpers for variant file access ---
.gcs_temp_dir <- file.path(tempdir(), "burdenem_var_files")

# Map dataset to GCS subdirectory
gcs_subdir_for <- function(dataset_prefix) {
    if (grepl("^genebass", dataset_prefix)) return("genebass/var_files")
    if (grepl("^aou_", dataset_prefix)) return("aou/var_files")
    stop("No GCS mapping for dataset prefix: ", dataset_prefix)
}

# Find variant files on GCS matching a dataset's sumstats pattern
# Returns local paths after downloading
gcs_download_variant_files <- function(dataset, annotation, sumstats_pattern, verbose = FALSE) {
    gcs_prefix <- opt$gcs_var_prefix
    gcs_dir <- paste0(gcs_prefix, "/", gcs_subdir_for(dataset))

    # Extract phenocode from sumstats_filename_pattern
    # Pattern looks like: data/var_txt/^aou_afr_3000963_<ANNOTATION>_.*\.txt\.bgz$
    pattern_base <- basename(sumstats_pattern)
    # Remove regex anchors and escapes to get a glob-friendly pattern
    # e.g., ^aou_afr_3000963_pLoF_.*\.txt\.bgz$ → aou_afr_3000963_pLoF_*
    clean_pattern <- sub("^\\^", "", pattern_base)                  # remove leading ^
    clean_pattern <- sub("\\$$", "", clean_pattern)                  # remove trailing $
    clean_pattern <- gsub("\\.\\*", "*", clean_pattern)              # .* → *
    clean_pattern <- gsub("\\\\\\.", ".", clean_pattern)             # \\. → .
    clean_pattern <- stringr::str_replace(clean_pattern, "<ANNOTATION>", annotation)

    gcs_glob <- paste0(gcs_dir, "/", clean_pattern)
    if (verbose) message(paste("  GCS listing:", gcs_glob))

    gcs_files <- tryCatch(
        system(paste("gsutil ls", shQuote(gcs_glob)), intern = TRUE, ignore.stderr = TRUE),
        error = function(e) character(0),
        warning = function(w) character(0)
    )
    gcs_files <- gcs_files[startsWith(gcs_files, "gs://")]

    if (length(gcs_files) == 0) {
        if (verbose) message(paste("  No GCS files found for:", gcs_glob))
        return(character(0))
    }

    # Download to temp dir (with caching — already-downloaded files are skipped)
    dir.create(.gcs_temp_dir, showWarnings = FALSE, recursive = TRUE)
    # Filter to files not already cached locally
    needed <- gcs_files[!file.exists(file.path(.gcs_temp_dir, basename(gcs_files)))]
    if (length(needed) > 0) {
        if (verbose) message(paste("  Downloading", length(needed), "files (", length(gcs_files) - length(needed), "cached)"))
        # Batch download with gsutil -m for parallelism
        gcs_args <- paste(shQuote(needed), collapse = " ")
        ret <- system(paste("gsutil -m cp", gcs_args, shQuote(.gcs_temp_dir)), intern = FALSE)
        if (ret != 0) warning("Some GCS downloads may have failed")
    } else if (verbose) {
        message(paste("  All", length(gcs_files), "files cached locally"))
    }
    local_files <- file.path(.gcs_temp_dir, basename(gcs_files))
    local_files <- local_files[file.exists(local_files)]
    return(local_files)
}

# Clean up temp files for a specific dataset/phenotype to free disk space
cleanup_gcs_temp <- function() {
    if (dir.exists(.gcs_temp_dir)) {
        unlink(.gcs_temp_dir, recursive = TRUE)
        message("Cleaned up temp variant files")
    }
}

# Load variant data for a single study (component or standalone), handling GCS if needed
load_variant_data_for_study <- function(study_row, annotation, freq_range, verbose = FALSE) {
    if (!is.null(opt$gcs_var_prefix)) {
        # GCS mode: download files, then load from temp dir
        pattern_with_anno <- stringr::str_replace(study_row$sumstats_filename_pattern, "<ANNOTATION>", annotation)
        local_files <- gcs_download_variant_files(
            dataset = study_row$dataset,
            annotation = annotation,
            sumstats_pattern = study_row$sumstats_filename_pattern,
            verbose = verbose
        )
        if (length(local_files) == 0) {
            warning(paste("No variant files found on GCS for:", study_row$dataset, study_row$identifier))
            return(data.frame(gene=character(), AF=numeric(), beta=numeric(), variant_variance=numeric(), functional_category=character()))
        }
        # Load from temp dir using the downloaded files
        variant_dir <- .gcs_temp_dir
        # Build a regex that matches the downloaded filenames
        file_pattern <- paste0("^", paste(basename(local_files), collapse = "|^"))
        # Simpler: use the original pattern's basename
        file_pattern <- basename(pattern_with_anno)
    } else {
        # Local mode
        pattern_with_anno <- stringr::str_replace(study_row$sumstats_filename_pattern, "<ANNOTATION>", annotation)
        variant_dir <- dirname(pattern_with_anno)
        file_pattern <- basename(pattern_with_anno)
    }

    load_variant_files_with_category(
        variant_dir = variant_dir,
        variant_file_pattern = file_pattern,
        annotations_to_process = c(annotation),
        frequency_range = freq_range
    )
}

# --- On-the-fly meta variant construction ---
# Instead of reading pre-generated meta variant files, load component files
# and combine them in memory. Supports both local and GCS sources.
construct_meta_variant_data <- function(study_row, annotation, freq_range, verbose) {
    meta_ds <- study_row$dataset
    identifier <- study_row$identifier
    components <- .meta_components[[meta_ds]]

    message(paste("Constructing", meta_ds, "variant data on-the-fly from components:", paste(components, collapse = ", ")))

    component_data_list <- list()
    for (comp_ds in components) {
        # Find the component study row in studies_df
        comp_study <- studies_df %>%
            dplyr::filter(identifier == !!identifier, dataset == comp_ds)

        if (nrow(comp_study) == 0) {
            warning(paste("Component study not found:", comp_ds, "for identifier:", identifier))
            next
        }
        comp_study <- comp_study[1, ]  # take first if duplicates

        if (verbose) message(paste("  Loading component:", comp_ds))

        # Load variant data (handles GCS or local transparently)
        comp_data <- tryCatch({
            load_variant_data_for_study(comp_study, annotation, freq_range, verbose)
        }, error = function(e) {
            warning(paste("Failed to load component", comp_ds, ":", e$message))
            return(NULL)
        })

        if (is.null(comp_data) || nrow(comp_data) == 0) {
            warning(paste("No variant data loaded for component:", comp_ds))
            next
        }

        # Ensure dataset column is set
        if (!"dataset" %in% names(comp_data)) {
            comp_data$dataset <- comp_ds
        }

        # For biobank_meta: harmonize ENSG → gene symbol for AoU components
        if (meta_ds == "biobank_meta" && !grepl("genebass", comp_ds) && !is.null(.ensg_to_symbol)) {
            ensg_mask <- grepl("^ENSG", comp_data$gene)
            if (any(ensg_mask)) {
                mapped <- .ensg_to_symbol[comp_data$gene[ensg_mask]]
                n_mapped <- sum(!is.na(mapped))
                n_total <- sum(ensg_mask)
                comp_data$gene[ensg_mask] <- ifelse(!is.na(mapped), mapped, comp_data$gene[ensg_mask])
                if (verbose) message(paste("    Gene name mapping:", n_mapped, "of", n_total, "ENSG IDs mapped to symbols"))
            }
        }

        message(paste("  Component", comp_ds, ":", nrow(comp_data), "variants"))
        component_data_list[[comp_ds]] <- comp_data
    }

    # Don't clean up — cached files are reused across studies for the same phenotype

    if (length(component_data_list) == 0) {
        warning(paste("No component data loaded for meta dataset:", meta_ds, "identifier:", identifier))
        return(data.frame(gene=character(), AF=numeric(), beta=numeric(), variant_variance=numeric(), functional_category=character()))
    }

    combined <- data.table::rbindlist(component_data_list, use.names = TRUE, fill = TRUE) %>%
        tibble::as_tibble()
    message(paste("Combined", meta_ds, "variant data:", nrow(combined), "total variants from", length(component_data_list), "components"))

    # Apply per-cohort LD pruning for binary traits (since we bypassed the override)
    is_binary <- all(c("AC_cases", "N") %in% names(combined))
    if (is_binary && length(.removed_variants_by_dataset) > 0 && "dataset" %in% names(combined)) {
        original_count <- nrow(combined)
        if (all(c("CHR", "POS") %in% names(combined))) {
            combined <- combined %>% dplyr::mutate(.locus = paste0(CHR, ":", POS))
        } else {
            combined <- combined %>% dplyr::mutate(.locus = NA_character_)
        }
        if (!all(is.na(combined$.locus))) {
            for (comp_ds in components) {
                removed_loci <- .removed_variants_by_dataset[[comp_ds]]
                if (!is.null(removed_loci) && length(removed_loci) > 0) {
                    combined <- combined %>%
                        dplyr::filter(!(dataset == comp_ds & .locus %in% removed_loci))
                }
            }
        }
        combined <- combined %>% dplyr::select(-.locus)
        n_removed <- original_count - nrow(combined)
        pct_removed <- if (original_count > 0) 100 * n_removed / original_count else 0
        message(sprintf("LD pruning (on-the-fly): removed %d of %d variants (%.2f%%) across components for %s",
                        n_removed, original_count, pct_removed, meta_ds))
    }

    return(combined)
}

# --- Download/load removed variants for LD pruning (binary traits only) ---
.removed_variants_by_dataset <- list()
.current_dataset <- NULL
.current_trait_type <- NULL

if (!is.null(opt$removed_variants_file)) {
    # Identify all unique datasets that have binary traits
    binary_datasets <- unique(studies_df$dataset[studies_df$trait_type == "binary"])
    # Expand meta datasets to component cohorts (aou_meta -> aou_afr, aou_amr, aou_eur)
    all_component_datasets <- unique(unlist(lapply(binary_datasets, get_component_datasets)))

    for (ds in all_component_datasets) {
        # Use ancestry suffix for <DATASET> placeholder (e.g., "aou_afr" -> "afr")
        # since removed variant files are named by ancestry, not full dataset name
        ancestry <- get_ancestry_from_dataset(ds)
        resolved_path <- stringr::str_replace(opt$removed_variants_file, "<DATASET>", ancestry)

        if (startsWith(resolved_path, "gs://")) {
            # Download from GCS
            local_file <- tempfile(pattern = paste0("removed_variants_", ds, "_"), fileext = ".tsv")
            message(paste("Downloading removed variants for", ds, "from:", resolved_path))
            system_result <- system(paste("gsutil cp", resolved_path, local_file), intern = FALSE)
            if (system_result != 0) stop(paste("Failed to download removed variants for", ds, "from:", resolved_path))
        } else {
            # Local file
            local_file <- resolved_path
            if (!file.exists(local_file)) stop(paste("Removed variants file not found:", local_file))
        }

        removed_variants_df <- data.table::fread(local_file)
        if ("locus" %in% names(removed_variants_df)) {
            removed_loci <- removed_variants_df$locus
        } else if ("SNP" %in% names(removed_variants_df)) {
            removed_loci <- sub("^([^:]+:[^:]+):.*", "\\1", removed_variants_df$SNP)
        } else {
            stop(paste("Cannot identify locus column in removed variants file for", ds))
        }
        message(paste("Loaded", length(removed_loci), "variant loci to remove for", ds))
        .removed_variants_by_dataset[[ds]] <- removed_loci
    }

    # Override load_variant_files_with_category for per-cohort LD pruning (meta datasets)
    original_load_variant_files <- load_variant_files_with_category
    load_variant_files_with_category <- function(variant_dir, variant_file_pattern, annotations_to_process, frequency_range) {
        variant_data <- original_load_variant_files(
            variant_dir = variant_dir, variant_file_pattern = variant_file_pattern,
            annotations_to_process = annotations_to_process, frequency_range = frequency_range
        )

        current_ds <- .current_dataset
        # Only filter for binary traits
        if (!identical(.current_trait_type, "binary")) return(variant_data)

        if (current_ds %in% names(.meta_components)) {
            # Meta dataset: per-cohort filtering using the 'dataset' column
            if ("dataset" %in% names(variant_data)) {
                original_count <- nrow(variant_data)
                if (all(c("CHR", "POS") %in% names(variant_data))) {
                    variant_data <- variant_data %>% dplyr::mutate(.locus = paste0(CHR, ":", POS))
                } else if ("locus" %in% names(variant_data)) {
                    variant_data <- variant_data %>% dplyr::mutate(.locus = locus)
                } else {
                    warning("Cannot construct locus for meta dataset variant filtering.")
                    return(variant_data)
                }
                for (component_ds in .meta_components[[current_ds]]) {
                    removed_loci <- .removed_variants_by_dataset[[component_ds]]
                    if (!is.null(removed_loci) && length(removed_loci) > 0) {
                        variant_data <- variant_data %>%
                            dplyr::filter(!(dataset == component_ds & .locus %in% removed_loci))
                    }
                }
                variant_data <- variant_data %>% dplyr::select(-.locus)
                n_removed <- original_count - nrow(variant_data)
                pct_removed <- if (original_count > 0) 100 * n_removed / original_count else 0
                message(sprintf("LD pruning: removed %d of %d variants (%.2f%%) across component cohorts for meta dataset %s",
                                n_removed, original_count, pct_removed, current_ds))
            }
        } else {
            # Non-meta dataset: filter using this dataset's removed variants
            removed_loci <- .removed_variants_by_dataset[[current_ds]]
            if (!is.null(removed_loci) && length(removed_loci) > 0) {
                original_count <- nrow(variant_data)
                if (all(c("CHR", "POS") %in% names(variant_data))) {
                    variant_data <- variant_data %>%
                        dplyr::mutate(.locus = paste0(CHR, ":", POS)) %>%
                        dplyr::filter(!(.locus %in% removed_loci)) %>%
                        dplyr::select(-.locus)
                } else if ("locus" %in% names(variant_data)) {
                    variant_data <- variant_data %>% dplyr::filter(!(locus %in% removed_loci))
                }
                n_removed <- original_count - nrow(variant_data)
                pct_removed <- if (original_count > 0) 100 * n_removed / original_count else 0
                message(sprintf("LD pruning: removed %d of %d variants (%.2f%%) for dataset %s",
                                n_removed, original_count, pct_removed, current_ds))
            }
        }
        return(variant_data)
    }
}

# Worker function to process a single study
.last_identifier <- ""  # Track phenotype for GCS cache cleanup

process_study_cli <- function(study_row, opt_config, freq_range_cli, icept_freq_bins_cli) {
    current_study <- study_row

    # Clean up GCS cache when switching to a new phenotype
    if (!is.null(opt_config$gcs_var_prefix) && current_study$identifier != .last_identifier) {
        if (nchar(.last_identifier) > 0) cleanup_gcs_temp()
        .last_identifier <<- current_study$identifier
    }

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
        # Keep the entire relative path from the studies file (e.g., 'fitted_models/...')
        model_output_path_with_anno <- file.path(output_dir, model_output_path_with_anno)
        if(opt_config$verbose) {
            message(paste("(--name) studies_dir:", studies_dir))
            message(paste("(--name) output_dir:", output_dir))
            message(paste("(--name) final model_output_path_with_anno:", model_output_path_with_anno))
        }
    }

    # Optionally skip if output already exists
    if (isTRUE(opt_config$skip_existing) && file.exists(model_output_path_with_anno)) {
        message(paste("Skipping existing:", model_output_path_with_anno))
        return(NULL)
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

    # --- 4. On-the-fly meta construction or GCS download ---
    preloaded_variant_data <- NULL
    if (opt_config$verbose) {
        message(paste("  GCS mode:", !is.null(opt_config$gcs_var_prefix),
                      "| is_meta:", is_meta_dataset(current_study$dataset),
                      "| dataset:", current_study$dataset))
    }
    if (is_meta_dataset(current_study$dataset)) {
        # Construct meta variant data on-the-fly from component datasets
        preloaded_variant_data <- construct_meta_variant_data(
            study_row = current_study,
            annotation = opt_config$annotation,
            freq_range = freq_range_cli,
            verbose = opt_config$verbose
        )
    } else if (!is.null(opt_config$gcs_var_prefix)) {
        # Non-meta study with GCS: download variant files, then load from temp
        preloaded_variant_data <- tryCatch({
            load_variant_data_for_study(current_study, opt_config$annotation, freq_range_cli, opt_config$verbose)
        }, error = function(e) {
            warning(paste("Failed to load variant data from GCS for", current_study$dataset, ":", e$message))
            NULL
        })
        # Don't clean up here — cached files are reused by meta studies for the same phenotype
    }

    # --- 5. Prepare arguments for run_burdenEM_rvas for the current study ---
    run_args <- list(
        variant_dir = derived_variant_dir,
        variant_file_pattern = derived_variant_file_pattern,
        ld_corrected_scores_file = processed_genes_file,
        output_file_prefix = output_file_prefix_for_run,
        annotation_to_process = opt_config$annotation,
        feature_col_name = opt_config$feature_col_name,
        num_feature_bins = opt_config$num_feature_bins,
        num_positive_components = opt_config$num_positive_components,
        burdenem_grid_size = opt_config$burdenem_grid_size,
        num_iter = opt_config$num_iter,
        per_allele_effects = opt_config$per_allele_effects,
        correct_for_ld = opt_config$correct_for_ld,
        frequency_range = freq_range_cli,
        intercept_frequency_bin_edges = icept_freq_bins_cli,
        verbose = opt_config$verbose,
        binary_trait_model_type = opt_config$binary_trait_model_type,
        optimizer = opt_config$optimizer,
        variant_data = preloaded_variant_data
    )

    run_args <- run_args[!sapply(run_args, is.null)]

    # Set globals for variant filtering override (used by load_variant_files_with_category)
    .current_dataset <<- current_study$dataset
    .current_trait_type <<- current_study$trait_type

    main_message <- paste0("Processing: ", current_study$identifier, " (", current_study$dataset, ") [", opt_config$annotation, "]")
    if(opt_config$verbose) {
        main_message <- paste0("Running BurdenEM for: ", current_study$identifier, " (", current_study$dataset, ") with annotation: ", opt_config$annotation,
                           ". Output prefix: ", output_file_prefix_for_run)
    }
    message(main_message)
    tryCatch({
        do.call(run_burdenEM_rvas, run_args)
    }, error = function(e) {
        message(paste0("ERROR processing study ", current_study$identifier, " (", current_study$dataset, "): ", e$message))
    })
    # Free per-study memory before moving on (variant data, model fits, etc.)
    rm(run_args, preloaded_variant_data)
    gc(verbose = FALSE)
    return(NULL)
}

# Main loop to process each study defined in the TSV file
message(paste("Processing", nrow(studies_df), "studies from:", studies_file_path))

# Determine effective worker count (capped to limit memory)
resolve_num_workers <- function() {
    if (opt$no_parallel || nrow(studies_df) <= 1) return(1L)
    avail <- future::availableCores()
    if (is.na(avail) || avail <= 1) return(1L)
    requested <- if (!is.null(opt$num_workers)) opt$num_workers else min(4L, avail)
    max(1L, min(as.integer(requested), as.integer(avail)))
}
num_workers <- resolve_num_workers()

run_chunk <- function(row_indices, label = NULL) {
    if (length(row_indices) == 0) return(invisible(NULL))
    if (!is.null(label)) message(paste0("--- ", label, " (", length(row_indices), " studies) ---"))
    if (num_workers <= 1) {
        purrr::map(row_indices,
                   ~process_study_cli(studies_df[.x, ], opt, frequency_range_cli, intercept_frequency_bin_edges_cli))
    } else {
        future::plan(future::multisession, workers = num_workers)
        on.exit(future::plan(future::sequential), add = TRUE)
        furrr::future_map(
            .x = row_indices,
            .f = ~process_study_cli(studies_df[.x, ], opt, frequency_range_cli, intercept_frequency_bin_edges_cli),
            .progress = FALSE,
            .options = furrr_options(seed = TRUE)
        )
    }
    # Free GCS cache between chunks (each phenotype's component files are no longer needed)
    if (!is.null(opt$gcs_var_prefix)) cleanup_gcs_temp()
    gc(verbose = FALSE)
    invisible(NULL)
}

if (opt$no_chunk_by_identifier) {
    message(paste("Processing all", nrow(studies_df), "studies in one batch with", num_workers, "worker(s)."))
    run_chunk(seq_len(nrow(studies_df)))
} else {
    # Group by phenotype identifier so all 6 datasets sharing component files
    # process together — maximizes GCS cache reuse and bounds peak disk usage.
    chunks <- split(seq_len(nrow(studies_df)), studies_df$identifier)
    # Preserve sorted order
    chunks <- chunks[unique(studies_df$identifier)]
    message(paste("Processing", length(chunks), "phenotypes (", nrow(studies_df),
                  "studies total) with", num_workers, "worker(s) per phenotype."))
    chunk_idx <- 0
    for (ident in names(chunks)) {
        chunk_idx <- chunk_idx + 1
        run_chunk(chunks[[ident]],
                  label = sprintf("[%d/%d] phenotype: %s", chunk_idx, length(chunks), ident))
    }
}

message("All studies processed.")
