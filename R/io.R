#Functions for input/output and data preprocessing

process_data_trio <- function(input_data,
                              features) {

  if (!is.null(input_data$case_rate)) {
    input_data$expected_count = 2 * input_data$N * input_data$case_rate
  }

  if (any(is.na(input_data))) {
    stop("NAs present in input data, please check")
  } else if ( !is.null(features) & any(is.na(features))) {
    stop("NAs present in features, please check")
  }

  if (!is.null(features) & !(all(rownames(input_data) == rownames(features)))) {
    stop("features rownames do not match input data rownames, please check")
  }


  return(input_data)
}

process_data_rvas <- function(input_data,
                              features){

  if(is.null(input_data$effect_estimate)){
    if(!is.null(input_data$z_score) & !is.null(input_data$AF) & !is.null(input_data$N)){
      input_data$effect_estimate <- input_data$z_score/sqrt(2*input_data$AF*(1-input_data$AF)*(input_data$N + input_data$z_score^2))
    }else{
      stop("effect estimates is missing from the data, please check")
    }
  }else{
    if (is.null(input_data$effect_se)) {
      stopifnot(!is.null(input_data$p_value))  # p_value must be specified if effect_se is not
      input_data$z_score <- qnorm(1 - input_data$p_value / 2) * sign(input_data$effect_estimate)
      input_data$effect_se <- abs(input_data$effect_estimate / input_data$z_score)
    }
    stopifnot(length(input_data$effect_se) == length(input_data$effect_estimate))
  }

  if (any(is.na(input_data))) {
    stop("NAs present in input data, please check")
  } else if ( !is.null(features) & any(is.na(features))) {
    stop("NAs present in features, please check")
  }

  if (is.null(input_data$N)){
    warning('Sample size (`N`) is missing from the data, which might be required for heritability estimation')
  }

  if (is.null(input_data$CAF)){
    warning('Combined Allele Frequency (`CAF`) is missing from the data, which might be required for power analyses')
  }

  if (is.null(input_data$trait_type)){
    warning('Trait type (`trait_type`) is missing from the data,\nwhich should be one of ("binary", "continuous") and is required for specifying likelihood model and power calculation.
            \nForcing it to be "continuous" for now')
    input_data$trait_type <- 'continuous'
  }else{
    if(tolower(unique(input_data$trait_type)) %in% BINARY_TRAIT_TYPES){
      input_data$trait_type <- 'binary'
    }else{
      input_data$trait_type <- 'continuous'
    }
  }

  if (is.null(input_data$AC_cases)){
    warning('Allele Counts in Cases (`AC_cases`) is missing from the data, which is required for binary traits')
  }

  if (!is.null(features) ) {
    if(!(all(rownames(input_data) == rownames(features)))){
      stop("features rownames do not match input data rownames, please check")
    }
    if(!all(rowSums(features) == 1) & all(features >= 0)){
      stop("features need to be all positive with rowSum of 1, please check")
    }
  }

  return(input_data)
}



#' Load Variant Files and Assign Functional Category
#'
#' Finds variant files matching a pattern, reads them, assigns a functional
#' category based on the filename and predefined mapping rules, filters by
#' allowed annotations, and combines them into a single dataframe.
#'
#' @param variant_dir Directory containing the variant files.
#' @param data_name The dataset name prefix used in filenames (e.g., 'genebass').
#' @param pheno The phenotype name.
#' @param annotations_to_process A character vector of functional annotations to include (e.g., c("pLoF", "missense|LC")).
#' @param frequency_range Optional numeric vector of length 2 specifying the min (inclusive) and max (exclusive) allele frequency (AF) range to keep. Default is NULL (no filtering).
#' @return A dataframe containing aggregated variant data across files, filtered by AF if specified, with columns appropriate for the detected trait type.
#' @importFrom dplyr bind_rows filter mutate select case_when any_of
#' @importFrom readr read_tsv cols
#' @importFrom stringr str_match
#' @importFrom purrr map
#' @export
load_variant_files_with_category <- function(variant_dir, variant_file_pattern = NULL, data_name, pheno, annotations_to_process, frequency_range = NULL) {
  BINARY_TRAIT_TYPES <- c("binary", "categorical", "icd_first_occurrence")

  message(paste("Looking for variant files in:", variant_dir))

  if (!is.null(variant_file_pattern) && nzchar(variant_file_pattern)) {
    # If a specific pattern is provided by the CLI (via studies file), use it directly
    file_pattern <- variant_file_pattern
    message(paste("Using provided file pattern:", file_pattern))
  } else {
    # Fallback: Construct the regex pattern to match specific annotations if no direct pattern is given
    if (is.null(annotations_to_process) || length(annotations_to_process) == 0) {
      stop("No annotations provided to 'annotations_to_process' and no 'variant_file_pattern' was given.")
    }
    if (is.null(data_name) || is.null(pheno)){
        stop("data_name and pheno must be provided if variant_file_pattern is not specified.")
    }
    # Create a regex OR group for the annotations, e.g., (pLoF|missense\|LC)
    annotations_regex_part <- paste0("(", paste(annotations_to_process, collapse = "|"), ")")
    file_pattern <- paste0("^", data_name, "_", pheno, "_", annotations_regex_part, "_.*\\.txt\\.bgz$")
    message(paste("Constructed file pattern:", file_pattern))
  }

  variant_files <- list.files(path = variant_dir, pattern = file_pattern, full.names = TRUE)

  if (length(variant_files) == 0) {
    warning(paste("No variant files found matching pattern:", file_pattern, "in directory:", variant_dir))
    # Return empty dataframe with expected columns (defaulting to continuous structure)
    return(data.frame(gene=character(), AF=numeric(), beta=numeric(), variant_variance=numeric(), functional_category=character()))
  } else {
    message(paste("Found", length(variant_files), "potential variant files to process."))
  }

  detected_trait_type <- NULL # Initialize trait type detection

  variant_data_list <- purrr::map(variant_files, ~{
      file_path <- .x
      filename <- basename(file_path)
      message(paste("Processing:", filename))

      # Extract annotation part from filename
      # This part needs to be robust. If variant_file_pattern is used, data_name and pheno might not be perfectly aligned for this regex.
      # For now, we assume the <ANNOTATION> part in the pattern (which was substituted by opt$annotation_to_process)
      # IS the functional_category_label we want. This is true if variant_file_pattern was like ^dataset_pheno_ACTUALANNOTATION_.*.txt.bgz
      # If variant_file_pattern is more complex, this extraction might need to be smarter or rely on `annotations_to_process` directly.
      # The `annotations_to_process` from main.R is a single string now.
      functional_category_label <- annotations_to_process[1] # This is the single annotation string (e.g., "pLoF")
      # Old extraction logic, might be useful if functional_category_label needs to be derived differently:
      # annotation_match <- stringr::str_match(filename, paste0(data_name, "_", pheno, "_([^_]+(?:_[^_]+)*?)_")) 
      # annotation_part <- if (!is.na(annotation_match[1, 2])) annotation_match[1, 2] else NA
      # functional_category_label <- annotation_part # This was the old logic if derived from filename

      # The functional_category_label is now directly taken from annotations_to_process[1]
      # No further mapping or checking based on filename parts is needed here when variant_file_pattern is used.

      tryCatch({
          # Read file, explicitly trying to read trait_type as character
          data <- readr::read_tsv(file_path, show_col_types = FALSE, col_types = readr::cols(gene = "c", phenotype_key = 'c',
                                                                                             description = 'c', CHR = 'c' , trait_type = 'c',
                                                                                             .default = "d")) # Guess others

          if (nrow(data) > 0) {

             # --- Detect trait_type from first valid file --- #
             if (is.null(detected_trait_type)) {
               if (!"trait_type" %in% names(data)) {
                 message(paste("  -> Skipping file (missing 'trait_type' column):", basename(file_path)))
                 return(NULL)
               }
               first_trait_val <- tolower(data$trait_type[1])
               if (first_trait_val %in% c('continuous', BINARY_TRAIT_TYPES, recursive=TRUE)) {
                 detected_trait_type <<- first_trait_val # Use <<- to modify outer scope variable
                 message(paste("  -> Detected trait_type as:", detected_trait_type))
               } else {
                 stop(paste("  -> Skipping file (invalid 'trait_type' value found: ", data$trait_type[1] ,"):", basename(file_path)))
               }
             }
             # --------------------------------------------- #

             # Define required columns based on detected trait_type
             base_required_cols <- c("gene", "AF", "beta", "variant_variance")
             categorical_required_cols <- c("N", "AC_cases")
             required_cols <- if (detected_trait_type %in% BINARY_TRAIT_TYPES) c(base_required_cols, categorical_required_cols) else base_required_cols

             if (!all(required_cols %in% names(data))) {
                 missing_cols <- setdiff(required_cols, names(data))
                 message(paste("  -> Skipping file (missing required columns for", detected_trait_type, "trait:", paste(missing_cols, collapse=", "), ")"))
                 return(NULL)
             }

             data <- data %>% dplyr::mutate(functional_category = functional_category_label)
             message(paste("  -> Read", nrow(data), "variants, assigned category:", functional_category_label))

             # Recalculate AF to ensure consistency with AC_cases > 0
             if (detected_trait_type %in% BINARY_TRAIT_TYPES) {
                 data <- data %>% dplyr::mutate(AF = pmax(AF, AC_cases / (2*N), na.rm = TRUE))
             }

             # --- Add AF Filtering --- #
             if (!is.null(frequency_range)) {
               original_count <- nrow(data)
               data <- data %>% dplyr::filter(AF >= frequency_range[1] & AF < frequency_range[2])
               filtered_count <- nrow(data)
               if (original_count > filtered_count) {
                   message(paste0("   Filtered out ", original_count - filtered_count, " variants based on AF range. ", filtered_count, " remaining."))
               }
             }
             # -------------------------- #

             # --- Process based on detected trait_type --- #
             if (detected_trait_type %in% BINARY_TRAIT_TYPES) {
                data <- data %>%
                    dplyr::mutate(expected_count = 2 * N * prevalence * AF) %>%
                    dplyr::select(dplyr::any_of(c("gene", "AF", "beta", "variant_variance", "expected_count", "AC_cases", "N", "functional_category", "prevalence")))
                  
                message("Expected over observed AC_cases: ", mean(data$AC_cases)/mean(data$expected_count))

             } else { # Continuous
                data <- data %>%
                    dplyr::select(dplyr::any_of(c("gene", "AF", "beta", "variant_variance", "functional_category"))) # Select final columns
             }
             # -------------------------------- #
             data
          } else {
             message("  -> File is empty.")
             NULL
          }
      }, error = function(e) {
          message(paste("  -> Error reading file:", e$message))
          NULL
      })
  })

  variant_data_list <- variant_data_list[!sapply(variant_data_list, is.null)]
  if (length(variant_data_list) == 0) {
    stop("Failed to read any relevant variant files or all were empty/skipped/missing columns.")
  }

  # Combine valid dataframes
  combined_data <- dplyr::bind_rows(variant_data_list)
  message(paste("Loaded and combined data from", length(variant_data_list), "relevant files.", nrow(combined_data), "total variants."))


  return(combined_data)
}
