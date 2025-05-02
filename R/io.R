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
    if(tolower(unique(input_data$trait_type)) %in% c('binary', 'categorical', 'qualitative')){
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
#' @param prevalence Optional numeric value specifying the prevalence of the trait. Default is NULL.
#' @return A dataframe containing aggregated variant data across files, filtered by AF if specified, with columns appropriate for the detected trait type.
#' @importFrom dplyr bind_rows filter mutate select case_when any_of
#' @importFrom readr read_tsv cols
#' @importFrom stringr str_match
#' @importFrom purrr map
#' @export
load_variant_files_with_category <- function(variant_dir, data_name, pheno, annotations_to_process, frequency_range = NULL, prevalence = NULL) {

  message(paste("Looking for variant files in:", variant_dir))
  # Construct the regex pattern to match specific annotations
  if (is.null(annotations_to_process) || length(annotations_to_process) == 0) {
    stop("No annotations provided to 'annotations_to_process'.")
  }
  # Create a regex OR group for the annotations, e.g., (pLoF|missense\|LC)
  # Note: We assume annotation names don't need further complex regex escaping for now.
  annotations_regex_part <- paste0("(", paste(annotations_to_process, collapse = "|"), ")")
  file_pattern <- paste0("^", data_name, "_", pheno, "_", annotations_regex_part, "_.*\\.txt\\.bgz$")
  message(paste("Using file pattern:", file_pattern))

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
      annotation_match <- stringr::str_match(filename, paste0(data_name, "_", pheno, "_([^_]+(?:_[^_]+)*?)_"))
      annotation_part <- if (!is.na(annotation_match[1, 2])) annotation_match[1, 2] else NA

      # Map filename part to desired functional category label
      functional_category_label <- dplyr::case_when(
          is.na(annotation_part) ~ NA_character_,
          grepl("^pLoF", annotation_part, ignore.case = TRUE) ~ "pLoF",
          grepl("^missense", annotation_part, ignore.case = TRUE) ~ "missense|LC",
          grepl("^synonymous", annotation_part, ignore.case = TRUE) ~ "synonymous",
          TRUE ~ NA_character_ # Default for non-matching/other categories
      )

      if (is.null(functional_category_label)) {
          message(paste("  -> Skipping file (annotation pattern not recognized):", annotation_part))
          return(NULL)
      }

      if (!functional_category_label %in% annotations_to_process) {
          message(paste("  -> Skipping file (annotation '", functional_category_label, "' not in target list)"))
          return(NULL)
      }

      tryCatch({
          # Read file, explicitly trying to read trait_type as character
          data <- readr::read_tsv(file_path, show_col_types = FALSE,
                                 col_types = readr::cols(gene = "c", N="d", AC_cases="d", prevalence="d", trait_type = "c", .default = "d"))

          if (nrow(data) > 0) {

             # --- Detect trait_type from first valid file --- #
             if (is.null(detected_trait_type)) {
               if (!"trait_type" %in% names(data)) {
                 message(paste("  -> Skipping file (missing 'trait_type' column):", basename(file_path)))
                 return(NULL)
               }
               first_trait_val <- tolower(data$trait_type[1])
               if (first_trait_val %in% c('continuous', 'categorical')) {
                 detected_trait_type <<- first_trait_val # Use <<- to modify outer scope variable
                 message(paste("  -> Detected trait_type as:", detected_trait_type))
               } else {
                 stop(paste("  -> Skipping file (invalid 'trait_type' value found: ", data$trait_type[1] ,"):", basename(file_path)))
               }
             }
             # --------------------------------------------- #

             # Define required columns based on detected trait_type
             base_required_cols <- c("gene", "AF", "beta", "variant_variance")
             categorical_required_cols <- c("N", "AC_cases", "prevalence")
             required_cols <- if (detected_trait_type == 'categorical') c(base_required_cols, categorical_required_cols) else base_required_cols

             if (!all(required_cols %in% names(data))) {
                 missing_cols <- setdiff(required_cols, names(data))
                 message(paste("  -> Skipping file (missing required columns for", detected_trait_type, "trait:", paste(missing_cols, collapse=", "), ")"))
                 return(NULL)
             }

             data <- data %>% dplyr::mutate(functional_category = functional_category_label)
             message(paste("  -> Read", nrow(data), "variants, assigned category:", functional_category_label))

             # Recalculate AF to ensure consistency with AC_cases > 0
             data <- data %>% dplyr::mutate(AF = pmax(AF, AC_cases / (2*N), na.rm = TRUE))

             # Coerce AC_cases to integer
             data$AC_cases <- as.integer(data$AC_cases)

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
             if (detected_trait_type == 'categorical') {
                data <- data %>%
                    dplyr::mutate(expected_count = 2* N * prevalence * AF) %>%
                    dplyr::select(dplyr::any_of(c("gene", "AF", "beta", "variant_variance", "expected_count", "AC_cases", "N_cases", "functional_category"))) 
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
