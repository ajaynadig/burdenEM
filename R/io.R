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


