source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas_analysis/utils.R')

load_burden_test_data <- function(pheno, data_name, anc='eur', threshold=2.5e-6,
                                  cutoffs=c(2.5e-6, 10**(-5:-2), 0.05)){
  if(data_name == 'aou'){
    path <- get_aou_burden_test_path(pheno=pheno, anc=anc)
  }else{
    path <- get_genebass_burden_test_path(pheno=pheno)
  }
  burden_test_sub <- read_csv(path) %>%
    dplyr::mutate(
      effect_estimate = gamma_perSD,
      effect_se = sqrt(abs(intercept_vector)),
    ) %>%
    dplyr::mutate(
      Wald_stat = effect_estimate/effect_se,
      p_Wald = 2 * (1 - pnorm(abs(Wald_stat))),
      significant_gene = p_Wald < threshold
    ) %>%
    dplyr::mutate(
      p_bins = factor(cut(p_Wald, breaks = cutoffs, labels = FALSE, include.lowest = TRUE), levels=1:length(cutoffs))
      # p_bins = case_when(
      #   p_Wald < threshold ~ '0',
      #   p_Wald >= threshold & p_Wald < 1e-5 ~ '1',
      #   p_Wald >= 1e-5 & p_Wald < 1e-4 ~ '2',
      #   p_Wald >= 1e-4 & p_Wald < 1e-3 ~ '3',
      #   p_Wald >= 1e-3 & p_Wald < 1e-2 ~ '4',
      #   p_Wald >= 1e-2 & p_Wald < 5e-2 ~ '5',
      #   p_Wald >= 5e-2 ~ '6'
      # )
    )
  if(data_name == 'aou'){
    burden_test_sub <- burden_test_sub %>%
      mutate(
        ukb_phenotype = ukb_phenotypes[phenotype_key],
        aou_phenotype = phenotype_key,
        gene_symbol = gene_name, gene_id=gene
      )
  }else{
    burden_test_sub <- burden_test_sub %>%
      mutate(
        aou_phenotype = aou_phenotypes[phenotype_key],
        ukb_phenotype = phenotype_key,
        gene_symbol = gene
      )
  }
  burden_test_sub <- burden_test_sub
  return(burden_test_sub)
}
