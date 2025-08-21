#!/usr/bin/env Rscript
library(docopt)
library(readr)
library(dplyr)
library(stringr)

# --- Helper function to generate output table filenames ---
generate_output_tablename <- function(studies_file_path, command, dataset, annotation, name = NULL, pvalue_threshold = NULL) {
  studies_file_prefix <- sub("\\.studies\\.tsv$", "", basename(studies_file_path))
  output_dir_name <- if (!is.null(name) && nzchar(name)) name else studies_file_prefix
  top_dir <- file.path(dirname(studies_file_path), output_dir_name)

  # The output filename prefix should also be consistent
  studies_file_prefix <- output_dir_name

  base_filename_stem <- paste(studies_file_prefix, command, sep = ".")
  
  # Add pvalue_threshold if provided
  if (!is.null(pvalue_threshold)) {
    base_filename_stem <- paste0(base_filename_stem, ".pval", pvalue_threshold)
  }
  
  # Add dataset, annotation, and extension
  final_filename_stem <- paste(base_filename_stem, dataset, annotation, "tsv", sep = ".")
  
  # Construct full path
  output_path <- file.path(top_dir, "tables", final_filename_stem)
  
  return(output_path)
}

doc <- 'Step 2 CLI to analyze a fitted burdenEM model

Usage:
  cli.R heritability <studies_file> [--annotation=<ann>] [--trait_type=<tt>] [--verbose] [--no_parallel] [--name=<name>]
  cli.R polygenicity <studies_file> [--annotation=<ann>] [--trait_type=<tt>] [--verbose] [--no_parallel] [--name=<name>]
  cli.R replication <studies_file> [--primary_dataset=<pd>] [--pvalue_threshold=<pval>] [--annotation=<ann>] [--trait_type=<tt>] [--verbose] [--no_parallel] [--name=<name>]
  cli.R effect_replication <studies_file> [--primary_dataset=<pd>] [--annotation=<ann>] [--trait_type=<tt>] [--verbose] [--no_parallel] [--name=<name>]
  cli.R distribution <studies_file> [--annotation=<ann>] [--trait_type=<tt>] [--verbose] [--no_parallel] [--name=<name>]
  cli.R -h | --help
  cli.R --version

Options:
  <studies_file>             Path to the .studies.tsv file.
  -a, --annotation=<ann>     Annotation to use [default: pLoF].
  -t, --trait_type=<tt>      Trait type to process [default: all].
  -v, --verbose              Print extra output [default: FALSE].
  --no_parallel              Disable parallel processing [default: FALSE].
  -n, --name=<name>          Custom name for output directory and file prefix.
  --primary_dataset=<pd>     Identifier for the primary dataset in replication pairs (e.g., \'ukbb_eur\').
  --pvalue_threshold=<pval>  Two-tailed p-value threshold for significance in the primary study.
  -h, --help                 Show this screen.
  --version                  Show version.

'

args <- docopt(doc)

# Determine the command
command <- NULL
if (args$heritability) {
  command <- "heritability"
} else if (args$polygenicity) {
  command <- "polygenicity"
} else if (args$replication) {
  command <- "replication"
} else if (args$effect_replication) {
  command <- "effect_replication"
} else if (args$distribution) {
  command <- "distribution"
} else {
  stop("Error: No valid command specified. Use --help for usage.")
}

# Global setup for studies file and initial filtering if a command is active
studies_df <- NULL # Initialize to NULL
if (!is.null(command) || !(args$help || isTRUE(args$version))) { # Proceed if it's a real command run
    studies_file_path <- args$studies_file
    if (is.null(studies_file_path)){
        # This case should ideally be caught by docopt if <studies_file> is always required by commands
        stop("Error: <studies_file> argument is missing for the command.")
    }
    if (!file.exists(studies_file_path)) {
        stop(sprintf("Studies file not found: %s. Please provide a valid path.", studies_file_path))
    }
    if (args$verbose) message(sprintf("Reading studies file: %s", studies_file_path))
    studies_df <- readr::read_tsv(studies_file_path, show_col_types = FALSE)
    original_study_count <- nrow(studies_df)
    if (args$verbose) message(sprintf("Loaded %d studies.", original_study_count))

    # If a custom name is provided, update the model path to look in that directory
    if (!is.null(args$name) && nzchar(args$name)) {
        if (!("model_filename" %in% names(studies_df))) {
            stop("'model_filename' column not found in studies file, but --name argument was provided.")
        }
        studies_dir <- dirname(studies_file_path)
        studies_df <- studies_df %>%
            mutate(
                model_filename = file.path(
                    studies_dir,
                    args$name,
                    stringr::str_replace(model_filename, "<ANNOTATION>", args$annotation)
                )
            )
    }

    # Trait Type Filtering
    if (args$trait_type != "all") {
        if ("trait_type" %in% names(studies_df)) {
            if (args$verbose) message(sprintf("Filtering studies by trait_type: '%s'", args$trait_type))
            studies_df <- studies_df %>% dplyr::filter(trait_type == args$trait_type)
            if (args$verbose) message(sprintf("%d studies remaining after trait_type filter.", nrow(studies_df)))
        } else {
            if (args$verbose) message("Warning: 'trait_type' column not found in studies file, skipping trait_type filter.")
        }
    }

    
    if (nrow(studies_df) == 0 && original_study_count > 0 && !is.null(command)) {
        message("No studies remain after applying filters. Subcommand will not be executed with data.")
        # Depending on desired behavior, could exit here: quit(status = 0)
    }
}

# Main command execution block
if (!is.null(command)) {
    if (args$verbose) {
    message(sprintf("Executing command: '%s' with %d studies.", command, if(is.data.frame(studies_df)) nrow(studies_df) else 0))
}

# Dispatch to command-specific logic
if (command == "heritability") {
    if (nrow(studies_df) == 0 && original_study_count > 0) {
        # This means filters removed all studies.
        if (args$verbose) { message("Skipping heritability calculation: all studies were filtered out.") }
        message("No heritability results to display as all studies were filtered.")
    } else {
        # This means studies_df has rows, OR studies_df was empty from the start (original_study_count == 0).
        # In both these sub-cases, we proceed to calculate.
        if (args$verbose) { message("Sourcing heritability functions from luke/heritability.R") }
        source("luke/heritability.R")
        
        if (args$verbose) {
            message(sprintf("Calling calculate_heritability_for_studies with annotation: '%s', verbose: %s, for %d studies.", 
                            args$annotation, args$verbose, nrow(studies_df)))
        }
        
        heritability_results <- calculate_heritability_for_studies(
            studies_df = studies_df,
            annotation_param = args$annotation,
            verbose_param = args$verbose,
            run_sequentially_param = args$no_parallel
        )
        
        if (nrow(heritability_results) > 0) {
            if (args$verbose) {
                message("Heritability calculation complete. Results:")
                print(heritability_results)
            }

            # Generate output filename and save the table
            output_table_path <- generate_output_tablename(
                studies_file_path = args$studies_file,
                command = command, # 'heritability'
                dataset = "all", # Consolidated table for all datasets in studies_df
                annotation = args$annotation,
                name = args$name
            )
            
            # Ensure the output directory exists
            dir.create(dirname(output_table_path), showWarnings = FALSE, recursive = TRUE)
            
            cat("Saving results to:", output_table_path, "\n")
            readr::write_tsv(heritability_results, output_table_path)
            if (args$verbose) {
                message(sprintf("Heritability results saved to: %s", output_table_path))
            }

        } else {
            if (args$verbose) {
                message("Heritability calculation yielded no results or all studies failed. No table will be saved.")
            } else {
                message("No heritability results to display or save.")
            }
        }
    }
} else if (command == "polygenicity") {
    if (args$verbose) { message(sprintf("Sourcing polygenicity functions from luke/polygenicity.R")) }
    source("luke/polygenicity.R")

    if (args$verbose) {
        message(sprintf("Calling calculate_polygenicity_metrics with annotation: '%s', verbose: %s, for %d studies.", 
                        args$annotation, args$verbose, nrow(studies_df)))
    }
    polygenicity_results <- calculate_polygenicity_metrics(
        studies_df = studies_df,
        annotation = args$annotation, # Pass annotation from CLI
        verbose = args$verbose,    # Pass verbosity from CLI
        run_sequentially_param = args$no_parallel # Pass parallel flag
    )

    if (nrow(polygenicity_results) > 0) {
        if (args$verbose) {
            message("Polygenicity calculation complete. Results:")
            print(polygenicity_results)
        }

        output_table_path <- generate_output_tablename(
            studies_file_path = args$studies_file,
            command = command, # 'polygenicity'
            dataset = "all",   # Consolidated table
            annotation = args$annotation,
            name = args$name
        )
        
        dir.create(dirname(output_table_path), showWarnings = FALSE, recursive = TRUE)
        cat("Saving results to:", output_table_path, "\n")
        readr::write_tsv(polygenicity_results, output_table_path)
        if (args$verbose) {
            message(sprintf("Polygenicity results saved to: %s", output_table_path))
        }
    } else {
        if (args$verbose) {
            message("Polygenicity calculation yielded no results or all studies failed. No table will be saved.")
        } else {
            message("No polygenicity results to display or save.")
        }
    }
} else if (command == "replication") {
    if (is.null(args$primary_dataset)) {
        stop("Error: --primary_dataset is required for the replication command.")
    }
    if (args$verbose) { message(sprintf("Sourcing replication functions from luke/replication.R")) }
    source("luke/replication.R")

    if (args$verbose) {
        message(sprintf("Calling calculate_replication_metrics with primary_dataset: '%s', annotation: '%s', verbose: %s, for %d studies.", 
                        args$primary_dataset, args$annotation, args$verbose, nrow(studies_df)))
    }
    # Convert pval_threshold to numeric if it's not NULL
    pval_threshold_numeric <- if (!is.null(args$pvalue_threshold)) as.numeric(args$pvalue_threshold) else NULL

    replication_results <- calculate_replication_metrics(
        studies_df = studies_df,
        annotation_param = args$annotation,
        primary_dataset_identifier = args$primary_dataset,
        pval_threshold_param = pval_threshold_numeric,
        verbose_param = args$verbose,
        run_sequentially_param = args$no_parallel
    )

    if (nrow(replication_results) > 0) {
        # Reformat results before printing/saving
        replication_results <- replication_results %>%
            dplyr::mutate(abbreviation = primary_abbreviation) %>%
            dplyr::select(
                abbreviation,
                primary_dataset,
                replication_dataset,
                annotation,
                expected_replication_rate,
                observed_replication_rate,
                n_genes_tested_replication,
                n_significant_in_primary
            )

        if (args$verbose) {
            message("Replication analysis complete. Results:")
            print(replication_results)
        }

        output_table_path <- generate_output_tablename(
            studies_file_path = args$studies_file,
            command = command, 
            dataset = paste0("primary-", gsub("[^A-Za-z0-9_]", "", args$primary_dataset)), # Sanitize primary_dataset for filename
            annotation = args$annotation,
            name = args$name,
            pvalue_threshold = args$pvalue_threshold
        )
        
        dir.create(dirname(output_table_path), showWarnings = FALSE, recursive = TRUE)
        cat("Saving results to:", output_table_path, "\n")
        readr::write_tsv(replication_results, output_table_path)
        message(sprintf("Replication results saved to: %s", output_table_path))
    } else {
        if (args$verbose) {
            message("Replication analysis yielded no results or no valid pairs found. No table will be saved.")
        } else {
            message("No replication results to display or save.")
        }
    }
} else if (command == "effect_replication") {
    if (identical(args$trait_type, "binary")) {
        stop("Error: The 'effect_replication' subcommand does not support binary traits. Please use --trait_type continuous or all.")
    }
    if (is.null(args$primary_dataset)) {
        stop("Error: --primary_dataset is required for the effect_replication command.")
    }
    if (args$verbose) { message(sprintf("Sourcing effect replication functions from luke/effect_replication.R")) }
    source("luke/effect_replication.R")

    if (identical(args$trait_type, "all")) {
        if (args$verbose) { message("Warning: --trait_type is 'all'. Filtering to 'continuous' traits for effect replication as binary traits are not supported.") }
        studies_df <- studies_df %>% dplyr::filter(trait_type == "continuous")
        if (args$verbose && nrow(studies_df) == 0) {
            message("  After filtering for continuous traits, no studies remain.")
        }
    }

    if (args$verbose) {
        message(sprintf("Calling calculate_effect_replication_metrics with primary_dataset: '%s', annotation: '%s', verbose: %s, for %d studies.", 
                        args$primary_dataset, args$annotation, args$verbose, if(is.data.frame(studies_df)) nrow(studies_df) else 0))
    }
    effect_replication_results <- calculate_effect_replication_metrics(
        studies_df = studies_df,
        annotation_param = args$annotation,
        primary_dataset_identifier = args$primary_dataset,
        verbose_param = args$verbose,
        run_sequentially_param = args$no_parallel
    )

    if (!is.null(effect_replication_results) && nrow(effect_replication_results) > 0) {
        # No specific reformatting for effect_replication_results currently, print as is.
        if (args$verbose) {
            message("Effect replication analysis complete. Results:")
            print(effect_replication_results)
        }

        output_table_path <- generate_output_tablename(
            studies_file_path = args$studies_file,
            command = command, 
            dataset = paste0("primary-", gsub("[^A-Za-z0-9_]", "", args$primary_dataset)), # Sanitize primary_dataset for filename
            annotation = args$annotation,
            name = args$name
        )
        
        dir.create(dirname(output_table_path), showWarnings = FALSE, recursive = TRUE)
        cat("Saving results to:", output_table_path, "\n")
        readr::write_tsv(effect_replication_results, output_table_path)
        message(sprintf("Effect replication results saved to: %s", output_table_path))
    } else {
        if (args$verbose) {
            message("Effect replication analysis yielded no results or no valid pairs found. No table will be saved.")
        } else {
            message("No effect replication results to display or save.")
        }
    }
} else if (command == "distribution") {
    source("luke/distribution.R") # Source the relevant R script
    
    # Call the main calculation function
    num_points <- 20
    gene_proportions <- exp(seq(log(1/18000), 0, length.out = num_points))
    distribution_results <- calculate_estimated_distribution_for_studies(
        studies_df = studies_df,
        annotation = args$annotation,
        gene_proportions = gene_proportions,
        verbose = args$verbose,
        run_sequentially_param = args$no_parallel
    )
    
    # Save results if any
    if (!is.null(distribution_results) && nrow(distribution_results) > 0) {
        output_table_path <- generate_output_tablename(
            studies_file_path = args$studies_file,
            command = command, # 'distribution'
            dataset = "all", # Consolidated table for all datasets in studies_df
            annotation = args$annotation,
            name = args$name
        )
        
        dir.create(dirname(output_table_path), showWarnings = FALSE, recursive = TRUE)
        cat("Saving results to:", output_table_path, "\n")
        readr::write_tsv(distribution_results, output_table_path)
        if (args$verbose) {
            message(sprintf("Distribution metrics results saved to: %s", output_table_path))
        }
    } else {
        if (args$verbose) {
            message("Distribution metrics calculation yielded no results or all studies failed. No table will be saved.")
        } else {
            message("No distribution metrics results to display or save.")
        }
    }
}
}
