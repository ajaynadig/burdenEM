heritability_enrichment_table <- function(data,output) {
  heritability_est_df <- data.frame()
  enrich_est_df <- data.frame()

  for (prevfactornumber in 1:length(data$prev_factors)) {

    prev_factor = data$prev_factors[prevfactornumber]
    print(prev_factor)
    for (outputnum in 1:length(data$loop_vars$names)) {

      prev_nonmod = data$loop_vars$prevalences[outputnum]
      prev_mod = data$loop_vars$prevalences[outputnum] * prev_factor

      variant_class = if (grepl("PTV",data$loop_vars$subsets[outputnum])) {
        "PTV"
      } else if (grepl("Mis2",data$loop_vars$subsets[outputnum])) {
        "Mis2"
      }else if (grepl("Mis1",data$loop_vars$subsets[outputnum])) {
        "Mis1"
      }else if (grepl("Mis0",data$loop_vars$subsets[outputnum])) {
        "Mis0"
      }else if (grepl("Syn",data$loop_vars$subsets[outputnum])) {
        "Syn"
      }
      role = if (grepl("Proband",data$loop_vars$subsets[outputnum])) {
        "Proband"
      } else if (grepl("Sibling",data$loop_vars$subsets[outputnum])) {
        "Sibling"
      }

      model = output[[prevfactornumber]][[outputnum]]
      #model = burdenEM_models[[outputnum]]

      if (variant_class == "PTV") {
        heritability = model$heritability_output$total_h2 * data$ptv_scale_factor
        print(heritability)
        heritability_lower = model$heritability_output$heritability_CI[1] * data$ptv_scale_factor
        heritability_upper = model$heritability_output$heritability_CI[2] * data$ptv_scale_factor

      } else {
        heritability = model$heritability_output$total_h2
        heritability_lower = model$heritability_output$heritability_CI[1]
        heritability_upper = model$heritability_output$heritability_CI[2]
      }
      heritability_p = model$heritability_output$total_h2_p
      if (is.null(heritability_p)) {heritability_p <- NA}

      peneff = model$penetrance$effective_penetrance
      peneff_lower = model$penetrance$effective_penetrance_CI[1]
      peneff_upper = model$penetrance$effective_penetrance_CI[2]


      output_h2_df <- data.frame(name = data$loop_vars$names[outputnum],
                                 variant_class = variant_class,
                                 prev_nonmod = prev_nonmod,
                                 prev_mod = prev_mod,
                                 prev_factor = prev_factor,
                                 role = role,
                                 h2 = heritability,
                                 h2_lower = heritability_lower,
                                 h2_upper = heritability_upper,
                                 h2_p = heritability_p,
                                 peneff = peneff,
                                 peneff_lower = peneff_lower,
                                 peneff_upper = peneff_upper)
      heritability_est_df = rbind(heritability_est_df,output_h2_df)

      heritability_enrich = model$heritability_output$enrichment
      heritability_lower = model$heritability_output$enrich_CI[1,]
      heritability_upper = model$heritability_output$enrich_CI[2,]

      frac_h2 = model$heritability_output$frac_h2
      frac_expected = model$heritability_output$frac_expected

      output_enrich_df <- data.frame(name = data$loop_vars$names[outputnum],
                                     variant_class = variant_class,
                                     prev_nonmod = prev_nonmod,
                                     prev_mod = prev_mod,
                                     prev_factor = prev_factor,
                                     role = role,
                                     annot = c("LOEUF1_mu1","LOEUF1_mu2",paste0("LOEUF",2:5)),
                                     h2_enrich = heritability_enrich,
                                     h2_enrich_lower = heritability_lower,
                                     h2_enrich_upper = heritability_upper,
                                     frac_h2 = frac_h2,
                                     frac_expected = frac_expected)

      enrich_est_df = rbind(enrich_est_df, output_enrich_df)
    }

  }
  return(list(heritability = heritability_est_df,
              enrichment = enrich_est_df))
}

bootstrap_function <- function(model, function_to_bootstrap, ...) {
  args <- list(...)  # Capture additional arguments as a list

  lapply(seq_along(model$bootstrap_output$bootstrap_delta), function(iter) {
    model_boot <- model
    model_boot$conditional_likelihood <- model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[, iter], ]
    model_boot$features <- model_boot$features[model$bootstrap_output$bootstrap_samples[, iter], ]
    model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iter]]

    # If genetic_data is provided in ..., subset it
    if ("genetic_data" %in% names(args)) {
      genetic_data_boot <- args$genetic_data[model$bootstrap_output$bootstrap_samples[, iter], ]
      args$genetic_data <- genetic_data_boot  # Replace genetic_data in args with the bootstrapped version
    }

    # Call the function with the modified arguments
    do.call(function_to_bootstrap, c(list(model = model_boot), args))
  })
}

get_fraccase <- function(model,
                         genetic_data,
                         RR_thresh,
                         gamma_scaling_factor = 1) {
  sum(2*genetic_data$case_rate*posterior_expectation(model,
                                                     genetic_data,
                                                     function(x) {
                                                       exp(gamma_scaling_factor*x) * (exp(x) > RR_thresh)
                                                     },
                                                     grid_size = 10))
}

get_fraccase_df <- function(data,
                            modelPTV,genetic_dataPTV,
                            modelMis2, genetic_dataMis2,
                            gamma_scaling_factor = 1,
                            boot = TRUE,
                            RR_range = c(2:20)) {
  fraccases_effsizethresh_df = data.frame()

  for (RR_thresh in RR_range) {
    print(RR_thresh)
    frac_casesgreater_pergene_PTV =  get_fraccase(modelPTV,
                                                  genetic_dataPTV,
                                                  gamma_scaling_factor,
                                                  RR_thresh)

    frac_casesgreater_pergene_Mis2 =  get_fraccase(modelMis2,
                                                   genetic_dataMis2,
                                                   gamma_scaling_factor,
                                                   RR_thresh)
    frac_casesgreater_combined = frac_casesgreater_pergene_PTV*data$ptv_scale_factor + frac_casesgreater_pergene_Mis2


    #bootstrap
    if (boot) {
      PTV_boot <- unlist(bootstrap_function(modelPTV,
                                            get_fraccase,
                                            genetic_data = genetic_dataPTV,
                                            RR_thresh = RR_thresh,
                                            gamma_scaling_factor = gamma_scaling_factor))
      Mis2_boot <- unlist(bootstrap_function(modelMis2,
                                             get_fraccase,
                                             genetic_data = genetic_dataMis2,
                                             RR_thresh = RR_thresh,
                                             gamma_scaling_factor = gamma_scaling_factor))
      frac_casesgreater_combined_boot = PTV_boot*data$ptv_scale_factor + Mis2_boot
    } else {
      PTV_boot <- rep(NA, length(modelPTV$bootstrap_output$bootstrap_delta))
      Mis2_boot <- rep(NA, length(modelPTV$bootstrap_output$bootstrap_delta))
      frac_casesgreater_combined_boot = PTV_boot*data$ptv_scale_factor + Mis2_boot
    }

    iter_df <- data.frame(RR_thresh = RR_thresh,
                          fraccases_greater_PTV = frac_casesgreater_pergene_PTV,
                          fraccases_greater_PTV_lowerCI = quantile(PTV_boot,0.025,na.rm = TRUE),
                          fraccases_greater_PTV_upperCI = quantile(PTV_boot,0.975,na.rm = TRUE),
                          fraccases_greater_Mis2 = frac_casesgreater_pergene_Mis2,
                          fraccases_greater_Mis2_lowerCI = quantile(Mis2_boot,0.025,na.rm = TRUE),
                          fraccases_greater_Mis2_upperCI = quantile(Mis2_boot,0.975,na.rm = TRUE),
                          fraccases_greater_PTVplusMis2 = frac_casesgreater_combined,
                          fraccases_greater_PTVplusMis2_lowerCI = quantile(frac_casesgreater_combined_boot, 0.025,na.rm = TRUE),
                          fraccases_greater_PTVplusMis2_upperCI = quantile(frac_casesgreater_combined_boot, 0.975,na.rm = TRUE),
                          gamma_scaling_factor = gamma_scaling_factor)

    fraccases_effsizethresh_df = rbind(fraccases_effsizethresh_df,iter_df)

  }
  return(fraccases_effsizethresh_df)
}

ngenes_for_fraccases <- function(data,
                                 modelPTV,genetic_dataPTV,
                                 modelMis2, genetic_dataMis2,
                                 gamma_scaling_factor = 1) {

  ngenes_for_fraccases <- sapply(2:20,
                                 function(RR_thresh) {
                                   frac_casesgreater_pergene_PTV =  (2*genetic_dataPTV$case_rate*posterior_expectation(modelPTV,
                                                                                                                       genetic_dataPTV,
                                                                                                                       function(x) {
                                                                                                                         exp(gamma_scaling_factor*x) * (exp(x) > RR_thresh)
                                                                                                                       },
                                                                                                                       grid_size = 10))

                                   frac_casesgreater_pergene_Mis2 =  (2*genetic_dataMis2$case_rate*posterior_expectation(modelMis2,
                                                                                                                         genetic_dataMis2,
                                                                                                                         function(x) {
                                                                                                                           exp(gamma_scaling_factor*x) * (exp(x) > RR_thresh)
                                                                                                                         },
                                                                                                                         grid_size = 10))

                                   fraccase_greater_combined = (frac_casesgreater_pergene_PTV*scale_ptv_factor) + frac_casesgreater_pergene_Mis2

                                   fraccase_greater_combined_ordered = fraccase_greater_combined[order(fraccase_greater_combined, decreasing = TRUE)]
                                   return(cumsum(fraccase_greater_combined_ordered))
                                 })

  colnames(ngenes_for_fraccases) = paste0("RR>",2:20)
  rownames(ngenes_for_fraccases) <- 1:nrow(ngenes_for_fraccases)
  return(ngenes_for_fraccases)
}

count_genes_by_effectsize <- function(model, genetic_data, genes_to_analyze = NULL, effsize_range =  c(2,5,10,20)) {

  if (is.null(genes_to_analyze)) {
    filter = rep(TRUE, nrow(genetic_data))
  } else {
    filter = rownames(genetic_data) %in% genes_to_analyze
  }
  effsizethresh_df = data.frame()

  for (RR_thresh in effsize_range) {
    num_genes =  sum(posterior_expectation(model,
                                           genetic_data,
                                           function(x) {exp(x) > RR_thresh},
                                           grid_size = 10)[filter])


    iter_df <- data.frame(RR_thresh = RR_thresh,
                          num_genes = num_genes)

    effsizethresh_df = rbind(effsizethresh_df, iter_df)

  }
  return(effsizethresh_df)
}


mutvar_per_gene <- function(model,genetic_data,prevalence,gamma_scaling_factor) {
  ((prevalence * 2 * genetic_data$case_rate)/(1 - prevalence)) * posterior_expectation(model,
                                                                                       genetic_data,
                                                                                       function(x) {
                                                                                         (exp(gamma_scaling_factor*x) -1)^2
                                                                                       },
                                                                                       grid_size = 10)
}

aggregate_mutvar_per_gene <- function(data,
                                      modelPTV,genetic_dataPTV,
                                      modelMis2, genetic_dataMis2,
                                      prevalence,
                                      gamma_scaling_factor = 1) {


  mutvar_pergene_PTV <- mutvar_per_gene(modelPTV,genetic_dataPTV,prevalence,gamma_scaling_factor)

  mutvar_pergene_Mis2 <- mutvar_per_gene(modelMis2,genetic_dataMis2,prevalence,gamma_scaling_factor)


  mutvar_pergene_combined = (mutvar_pergene_PTV*data$ptv_scale_factor) + mutvar_pergene_Mis2
  names(mutvar_pergene_combined) <- rownames(genetic_dataPTV)

  mutvar_pergene_ordered = mutvar_pergene_combined[order(mutvar_pergene_combined, decreasing = TRUE)]

  #bootstrap
  mutvar_pergene_PTV_boot <- bootstrap_function(modelPTV,
                                                mutvar_per_gene,
                                                genetic_data = genetic_dataPTV,
                                                prevalence = prevalence,
                                                gamma_scaling_factor = 1)

  mutvar_pergene_Mis2_boot <- bootstrap_function(modelMis2,
                                                 mutvar_per_gene,
                                                 genetic_data = genetic_dataMis2,
                                                 prevalence = prevalence,
                                                 gamma_scaling_factor = 1)

  mutvar_pergene_combined_boot <- sapply(1:length(mutvar_pergene_PTV_boot),
                                         function(i) {
                                           mutvar_pergene_combined_boot = (mutvar_pergene_PTV_boot[[i]]*data$ptv_scale_factor) + mutvar_pergene_Mis2_boot[[i]]
                                           names(mutvar_pergene_combined_boot) <- rownames(genetic_dataPTV)

                                           mutvar_pergene_ordered_boot = mutvar_pergene_combined_boot[order(mutvar_pergene_combined_boot, decreasing = TRUE)]
                                         })
  rownames(mutvar_pergene_combined_boot) <- 1:nrow(mutvar_pergene_combined_boot)
  return(list(full = mutvar_pergene_ordered,
              boot = mutvar_pergene_combined_boot))
}

burdenEM_power_forecasting <- function(old_genetic_data,
                                       model,
                                       N_new_case,
                                       gamma_scaling_factor = 1,
                                       effect_size_thresholds = c(0,2,5,10,20)) {

  #sample effect sizes from posteriors
  log_gamma_samples = posterior_gene_samples(model) * gamma_scaling_factor

  #Sample some new counts
  new_case_counts <- rpois(nrow(old_genetic_data), 2 * N_new_case * exp(log_gamma_samples) * old_genetic_data$case_rate)

  #combine
  combined_counts = new_case_counts + old_genetic_data$case_count

  new_genetic_data = data.frame(case_count = combined_counts,
                                case_rate = old_genetic_data$case_rate,
                                N = old_genetic_data$N +N_new_case)
  rownames(new_genetic_data) = rownames(old_genetic_data)
  new_genetic_data$expected_count = 2 * new_genetic_data$N * new_genetic_data$case_rate

  #Estimate some poisson P values
  poisson_p <- ppois(combined_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)

  #Sample some null counts
  null_counts = rpois(nrow(new_genetic_data), new_genetic_data$expected_count)

  #Get some null P values
  null_p <- ppois(null_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)


  # QQ plot
  qq <- ggplot(mapping = aes(y = -log10(poisson_p[order(poisson_p,decreasing = TRUE)]),
                             x = -log10(null_p[order(null_p,decreasing = TRUE)])))+
    geom_point()+
    geom_abline()+
    geom_hline(yintercept = -log10( 0.05/length(poisson_p)), linetype = "dashed")+
    labs(x = "Simulated Null Quantiles", y = "Simulated Effect Quantiles")+
    theme_bw()

  #How many significant genes?
  count_table = data.frame()
  for (RR_thresh in effect_size_thresholds) {

    gene_greaterthan_thresh =posterior_expectation(model,
                                                   old_genetic_data,
                                                   function(x) {exp(x) > RR_thresh},
                                                   grid_size = 10)

    iter_df <- data.frame(threshold = RR_thresh,
                          count_fdr = sum((p.adjust(poisson_p, method = "fdr") < 0.05) * gene_greaterthan_thresh),
                          count_bonferroni = sum((poisson_p < 0.05/length(poisson_p)) * gene_greaterthan_thresh))
    count_table <- rbind(count_table,iter_df)
  }

  return(list(count_table = count_table,
              qq = qq))

}

burdenEM_power_forecasting_mod <- function(old_genetic_data,
                                           model,
                                           N_new_case,
                                           gamma_scaling_factor = 1,
                                           effect_size_thresholds = c(0,2,5,10,20)) {

  #sample effect sizes from posteriors. Define these to be the "true effect sizes"
  log_gamma_samples = posterior_gene_samples(model)

  #Scale the effect sizes for the new counts
  log_gamma_samples_scale =log_gamma_samples * gamma_scaling_factor

  #Sample some new counts
  new_case_counts <- rpois(nrow(old_genetic_data), 2 * N_new_case * exp(log_gamma_samples_scale) * old_genetic_data$case_rate)

  #combine
  combined_counts = new_case_counts + old_genetic_data$case_count

  new_genetic_data = data.frame(case_count = combined_counts,
                                case_rate = old_genetic_data$case_rate,
                                N = old_genetic_data$N +N_new_case)
  rownames(new_genetic_data) = rownames(old_genetic_data)
  new_genetic_data$expected_count = 2 * new_genetic_data$N * new_genetic_data$case_rate

  #Estimate some poisson P values
  poisson_p <- ppois(combined_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)

  #Sample some null counts
  null_counts = rpois(nrow(new_genetic_data), new_genetic_data$expected_count)

  #Get some null P values
  null_p <- ppois(null_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)


  # QQ plot
  qq <- ggplot(mapping = aes(y = -log10(poisson_p[order(poisson_p,decreasing = TRUE)]),
                             x = -log10(null_p[order(null_p,decreasing = TRUE)])))+
    geom_point()+
    geom_abline()+
    geom_hline(yintercept = -log10( 0.05/length(poisson_p)), linetype = "dashed")+
    labs(x = "Simulated Null Quantiles", y = "Simulated Effect Quantiles")+
    theme_bw()

  #How many significant genes?
  count_table = data.frame()
  for (RR_thresh in effect_size_thresholds) {

    #Rather than using the posterior distribution, simply use the sampled effect size as the 'true' effect size

    gene_greaterthan_thresh = exp(log_gamma_samples) >= RR_thresh
    iter_df <- data.frame(threshold = RR_thresh,
                          count_fdr = sum((p.adjust(poisson_p, method = "fdr") < 0.05) * gene_greaterthan_thresh),
                          count_bonferroni = sum((poisson_p < 0.05/length(poisson_p)) * gene_greaterthan_thresh))
    count_table <- rbind(count_table,iter_df)
  }

  return(list(count_table = count_table,
              qq = qq))

}

make_supptable <- function(tables,names_to_keep) {
  h2_table = tables$heritability
  enrich_table = tables$enrichment

  h2_table = h2_table[h2_table$name %in% names_to_keep,]
  enrich_table = enrich_table[enrich_table$name %in% names_to_keep,]

  # Reshape enrich_table to wide format
  enrich_wide <- reshape(enrich_table,
                         idvar = c("name", "variant_class", "prev_nonmod", "prev_mod", "prev_factor", "role"),
                         timevar = "annot",
                         direction = "wide")

  # Rename columns to follow the new naming convention
  colnames(enrich_wide) <- gsub("h2_enrich\\.", "Enrichment_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("h2_enrich_lower\\.", "Enrichment_Lower95CI_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("h2_enrich_upper\\.", "Enrichment_Upper95CI_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("frac_h2\\.", "FractionMutVar_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("frac_expected\\.", "FractionExpected_", colnames(enrich_wide))

  # Merge the reshaped enrich_wide with h2_table
  h2_table_merged <- merge(h2_table, enrich_wide,
                           by = c("name", "variant_class", "prev_nonmod", "prev_mod", "prev_factor", "role"),
                           all.x = TRUE)

  supptable = h2_table_merged[,c("name",
                                 "variant_class",
                                 "prev_mod",
                                 "role",
                                 "h2",
                                 "h2_lower",
                                 "h2_upper",
                                 "peneff",
                                 "peneff_lower",
                                 "peneff_upper",
                                 names(h2_table_merged)[16:43])]

  names(supptable) = c("Model Name",
                       "Variant Class",
                       "Prevalence",
                       "Role",
                       "MutVar",
                       "MutVar_lower95CI",
                       "MutVar_upper95CI",
                       "EffectivePenetrance",
                       "EffectivePenetrance_Lower95CI",
                       "EffectivePenetrance_Upper95CI",
                       names(h2_table_merged)[16:43])


  return(supptable)


}

#define function to get count tables
get_genetic_data <- function(i, data) {
  counts = data$counts
  subsets = data$loop_vars$subsets
  mutation_rate = data$loop_vars$mutation_rate
  N_subset = data$loop_vars$N_subset
  names = data$loop_vars$names
  print(names[i])

  case_count = rep(0, nrow(rowData(counts)))

  for (assay in subsets[[i]]) {
    if (length(data$loop_vars$datasets) > 0) {
      if (sum(data$loop_vars$datasets[[i]]) > 1) {
        assay_count = rowSums(assays(counts)[[assay]][,data$loop_vars$datasets[[i]]])
      } else {
        assay_count = assays(counts)[[assay]][,data$loop_vars$datasets[[i]]]
      }
    } else {
      assay_count = rowSums(assays(counts)[[assay]])

    }
    case_count = case_count + assay_count
  }

  mu = rowData(counts)[,mutation_rate[i]]
  #  mutation_rate = rowData(counts)$mut.ptv.gnomad.v2

  posterior_mu_factor = rowData(counts)$PosteriorMuCorrectionFactor

  N = N_subset[i]

  input_df <- data.frame(case_count = case_count,
                         case_rate = mu * posterior_mu_factor,
                         N = N)
  rownames(input_df) = rownames(rowData(counts))

  filter_impossiblecount = !(input_df$case_count > 0 & input_df$case_rate == 0)


  genetic_data = process_data_trio(input_df[filter_impossiblecount,],
                                   data$features[filter_impossiblecount,])

  return(list(genetic_data = genetic_data,
              features = data$features[filter_impossiblecount,]))

}

binomial_LL_withposterior <- function(c,data,totals,N_1,N_tot,model,genetic_data) {
  posterior_terms = posterior_expectation(model,genetic_data,function(x) {exp(x * (c-1))},10)
  probability_factors = (N_1/N_tot) * posterior_terms


  LLs = dbinom(data,totals,probability_factors,log = TRUE)
  return(sum(LLs))
}

binomial_analysis <- function(model,
                              genetic_data_total,
                              genetic_data_subsample,
                              param_grid = seq(0.6,1.35,length.out  = 30)) {

  counts_total = genetic_data_total$case_count
  counts_subsample = genetic_data_subsample$case_count
  N_total = genetic_data_total$N[1]
  N_subsample = genetic_data_subsample$N[1]
  print(N_subsample)

  LL_grid_posterior = sapply(param_grid,
                             binomial_LL_withposterior,
                             counts_subsample,
                             counts_total,
                             N_subsample,
                             N_total,
                             model,
                             genetic_data_total)

  return(list(LL_output = data.frame(val = param_grid,
                                     LL = LL_grid_posterior),
              MLE = param_grid[which.max(LL_grid_posterior)]))


}

get_scaled_heritability <- function(model,
                                    genetic_data,
                                    gamma_scaling_factor,
                                    prevalence,
                                    heritability_scaling_factor = 1,
                                    n_boot = 100) {
  h2_output_mod <- estimate_heritability_trio(model,
                                              genetic_data,
                                              prevalence,
                                              gamma_scaling_factor = gamma_scaling_factor)



  bootstrap_heritability_output_mod <- lapply(1:100,
                                              function(iter,model,genetic_data) {



                                                model_boot = model
                                                model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                                boot_heritability = estimate_heritability_trio(model = model_boot,
                                                                                               genetic_data = genetic_data[model$bootstrap_output$bootstrap_samples[,iter],],
                                                                                               prevalence = prevalence,
                                                                                               gamma_scaling_factor = gamma_scaling_factor)

                                              },
                                              model,
                                              genetic_data)

  bootstrap_h2_ests_mod = sapply(1:length(bootstrap_heritability_output_mod), function(x) bootstrap_heritability_output_mod[[x]]$total_h2)
  heritability_CI_mod = quantile(bootstrap_h2_ests_mod,c(0.025,0.975))

  return(list(h2 = h2_output_mod$total_h2 * heritability_scaling_factor,
              bootstrap_h2_ests = bootstrap_h2_ests_mod* heritability_scaling_factor,
              h2_lower = heritability_CI_mod[1]* heritability_scaling_factor,
              h2_upper = heritability_CI_mod[2]* heritability_scaling_factor))

}






