source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas.R')
source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
library(emojifont)
library(rtracklayer)
library(GenomicRanges)

### CONSTANTS
pop_colors['nfe'] = color_eur
pop_colors['global'] = 'black'
pop_colors['equal_rep'] = 'gray'
pop_colors['genebass'] = 'gray20'
pop_colors['mid'] = color_mde = '#EEA9B8'

ANNOTATIONS=c('pLoF', 'missense_damaging', 'missense_benign', 'synonymous')
annotations=c('pLoF', 'missense_damaging', 'missense_benign', 'synonymous')
annotation_types2 = annotations
annotation_names2 = c('pLoF', 'Missense (Damaging)', 'Missense (Benign)', 'Synonymous')
names(annotation_names2) = annotation_types2
colors2 = c(colors, 'missense-benign' = 'gray30', 'missense_damaging' = color_mis,
            'missense-notbenign'=color_mis, 'missense_benign' = 'gray30')
annotation_fill_scale2 = scale_fill_manual(name = 'Annotation', values = colors2, breaks = annotation_types2, labels = annotation_names2)
annotation_color_scale2 = scale_color_manual(name = 'Annotation', values = colors2, breaks = annotation_types2, labels = annotation_names2)

ukb_p_bins_names <- c('[0, 2.5e-6]', '(2.5e-6, 1e-5]', '(1e-5, 1e-4]', '(1e-4, 1e-3]', '(1e-3, 1e-2]', '(1e-2, 5e-2]', '(5e-2, 1]')

pops = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'))
pops = c('African/African American (AFR)', 'American Admixed/Latino (AMR)', 'East Asian (EAS)', 'European (EUR)', 'Middle Eastern (MID)',' South Asian (SAS)', 'Other (OTH)')
names(pops) <- c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth')
label_type = labeller(annotation=annotation_names2, ancestry = pops)

sim_names <- c('realistic','more_polygenic', 'more_positive_effects','small_N','large_N','strong_selection', 'strong_popstrat','less_polygenic')
sim_labels <- c('Realistic', 'Polygenicity++', 'Positive effect++', 'N--', 'N++', 'Selection++', 'Stratification++', 'Polygenicity--')
sample_sizes <- c(1e5, 1e5, 1e5, 3e4, 5e5, 1e5, 1e5, 1e5)
simulation_colors <- c("#887a7e", "#8f89bb", "#bea0b6", "#d9887e", "#7ea8bb", "#c5a77e","#718391", "#9bd3a1")

names(sample_sizes) <- sim_names
names(sim_labels) <- sim_names
names(simulation_colors) <- sim_labels

random_phenonames <- rep(c('Quantitaitve', 'Binary (p = 0.1%)', 'Binary (p = 1%)',
                    'Binary (p = 10%)', 'Binary (p = 20%)', 'Binary (p = 50%)'), each=5)
aou_random_phenos <- c(paste0('random_0.5_continuous_', 1:5), paste0('random_0.5_0.001_', 1:5),
                       paste0('random_0.5_0.01_', 1:5), paste0('random_0.5_0.1_', 1:5),
                       paste0('random_0.5_0.2_', 1:5), paste0('random_0.5_0.5_', 1:5))
names(random_phenonames) <- aou_random_phenos

ukb_phenotypes <- c("100017_NA", "4080_NA",  "4079_NA",  "21001_NA",  "21002_NA",  "30040_NA",  "30060_NA",  "30070_NA",  "30080_NA",  "30100_NA",
                    "30190_NA", "30200_NA", "30220_NA", "30530_NA", "30610_NA", "30650_NA", "30670_NA", "30690_NA", "30740_NA", "30760_NA", "30860_NA",
                    "30870_NA", "48_NA", "49_NA", "50_NA", "5983_NA", "WHR_custom_NA", "30730_NA", "30620_NA", "30780_NA",
                    '20002_asthma', '20002_hypertension', '20002_diabetes', '130708_NA', '131306_NA', '130706_NA', '20002_hypertension', "20002_osteoarthritis",
                    "20002_osteoarthritis", '130708','20002_diabetes', '131706_NA', '131306_NA', '20002_asthma')
aou_phenotypes <- c('3001420', 'blood-pressure-systolic-mean', 'blood-pressure-diastolic-mean', 'BMI', 'weight', '3023599', '3009744', '3019897', '3024929', '3043111',
                    '3011948', '3008342', '3013869', '3019550', '3035995', '3013721', '3013682', '3027114', '3004501', '3007070', '3020630',
                    '3022192', 'waist-circumference-mean', 'hip-circumference-mean', 'height', 'heart-rate-mean', 'WHR', "3026910", "3006923", "3028288",
                    '495', '401', '250', '250.2', '411', '250.1', 'CV_401', 'MS_708', '740', 'EM_202.2', 'EM_202', 'EM_202.1', 'CV_404', 'RE_475')
names(aou_phenotypes) <- ukb_phenotypes
names(ukb_phenotypes) <- aou_phenotypes

themes <- theme_classic() + theme(
  text = element_text(family = "sans"),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 13),
  strip.text = element_text(size = 13),
  # title = element_text(size = 13, face='bold'),
)

### DATA PATHS
burdenEM_ROOT <- '~/Dropbox (Partners HealthCare)/burdenEM/'
burdenEM_RESULT_ROOT <- paste0(burdenEM_ROOT, 'burdenEM_results/')
SIM_RESULTS_PATH <- paste0(burdenEM_RESULT_ROOT, 'simulation_results_grid_size_100_num_iter_1000_n_boot_10_new_cpt_with_feature.csv')
TRUE_SIM_PATH <- '~/Dropbox (Partners HealthCare)/burdenEM_shared/simulated_sumstats/true_values.csv'

Genebass_BHR_paper_PATH <- paste0(burdenEM_RESULT_ROOT, 'BHR/genebass_bhr_paper.csv')
Genebass_BHR_PATH <- paste0(burdenEM_RESULT_ROOT, 'BHR/genebass_eur_tmp_full.csv')
Genebass_BurdenEM_PATH <- paste0(burdenEM_RESULT_ROOT, 'BurdenEM/genebass_eur_n_boot_10_grid_size_100_iter_1000_full_h2_result_tmp_burdenEM_union_cpt.csv')
Genebass_weight_PATH <- paste0(burdenEM_RESULT_ROOT, 'BurdenEM/genebass_eur_n_boot_10_grid_size_100_iter_1000_weight_result_tmp_burdenEM_union_cpt.csv')
Genebass_NTPR_PATH <- paste0(burdenEM_RESULT_ROOT, 'BurdenEM/genebass_eur_iter_1000_ntpr_result_tmp_burdenEM_union_cpt.csv')
common_var_polygenicity_PATH <- paste0(burdenEM_RESULT_ROOT, 'Weissbrod_M50_polygenicity.csv')

get_aou_bhr_path <- function(anc){
  path <- paste0(burdenEM_RESULT_ROOT, 'BHR/aou_', anc,'_tmp_full.csv')
  return(path)
}

# aou_bhr <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/aou_eur_tmp_full.csv') %>%
#   rbind(., read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/aou_amr_tmp_full.csv') ) %>%
#   rbind(., read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/aou_afr_tmp_full.csv') ) %>%
#   select(1:12) %>%
#   mutate(ukb_phenotype = ukb_phenotypes[phenotype_key])  %>%
#   merge(., aou_neff, by.x = c('ancestry', 'phenotype_key') , by.y = c('ancestry', 'phenoname') ) %>%
#   mutate(bhr_z = bhr_h2/bhr_h2_se,
#          bhr_p = 2 * (1 - pnorm(abs(bhr_z))),
#          bhr_significant = bhr_p < bhr_threshold)
#
# aou_genebass_bhr <- aou_bhr %>%
#   merge(., genebass_bhr, by.x = c('annotation', 'AF_bin', 'ukb_phenotype'), by.y = c('annotation', 'AF_bin', 'phenotype_key'), suffix = c('.aou', '.genebass'))


get_aou_burden_test_path <- function(pheno, anc){
  path <- paste0(burdenEM_RESULT_ROOT ,'Burden_test/aou_', anc,'/aou_', anc,'_', pheno,'_burden_test.csv')
  return(path)
}

get_genebass_burden_test_path <- function(pheno){
  path <- paste0(burdenEM_RESULT_ROOT ,'Burden_test/genebass_eur/genebass_eur_', pheno,'_burden_test.csv')
  return(path)
}

get_aou_burdenEM_ancestry_path <- function(anc){
  path = paste0(burdenEM_RESULT_ROOT, 'BurdenEM/aou_', anc,'_n_boot_10_grid_size_100_iter_1000_full_h2_result_tmp_burdenEM_union_cpt.csv')
  return(path)
}

AoU_BurdenEM_PATH <- paste0(burdenEM_RESULT_ROOT, 'BurdenEM/aou_eur_iter_1000_full_h2_result_tmp_burdenEM.csv') # Should eventually be Meta-analysis

### UTIL DATA
# gtf_path <- '~/Dropbox (Partners HealthCare)/analysis/lung_function_pleiotropy_yixuan/paper/Homo_sapiens.GRCh38.105.gtf'
# gtf_data <- import(gtf_path, format = "gtf")
# gtf_genes <- gtf_data[gtf_data$type == "gene"]
# gtf_gene_info <- data.frame(
#   gene_id = mcols(gtf_genes)$gene_id,
#   gene_symbol = mcols(gtf_genes)$gene_name
# )
# write_csv(gtf_gene_info, paste0(burdenEM_RESULT_ROOT, 'data/gtf_gene_info.csv'))
gtf_gene_info <- read_csv(paste0(burdenEM_RESULT_ROOT, 'data/gtf_gene_info.csv'))

aou_neff <- read_delim(paste0(burdenEM_RESULT_ROOT, 'data/aou/aou_n_eff.txt.bgz'), delim = '\t')
genebass_neff <- read_delim(paste0(burdenEM_RESULT_ROOT, 'data/genebass/genebass_n_eff.txt.bgz'), delim = '\t')
genebass_nvar <- read_csv(paste0(burdenEM_RESULT_ROOT, 'Burden_test/genebass_eur_tmp_burden_test_n_var.csv'))
genebass_pheno_info <- read_delim('~/Dropbox (Partners HealthCare)/aou/aou_genebass_meta/genebass_pheno_info.txt.bgz', delim = '\t') %>%
  mutate(phenotype = paste0(phenocode, '_', coding_description),
         N = n_cases + if_else(is.na(n_controls), 0, n_controls),
         description = if_else(phenocode=='20002', coding_description, description))

### FUNCTIONS
load_burdenEM_sim_true_values <- function(true_sim_path, n_rep=100, n_h2=5, n_scenario=8){
  true_value <- read_csv(true_sim_path)
  colnames(true_value) <- paste0('true_', colnames(true_value))
  true_value <- true_value %>%
    mutate(name = rep(rep(c('realistic','small_N','large_N','strong_selection',
                            'strong_popstrat','more_polygenic','less_polygenic',
                            'more_positive_effects'), each = n_h2), n_rep),
           rep = rep(rep(1:n_rep, each=n_h2*n_scenario)),
           true_h2_rounded = round(true_burden_h2, 3))
  return(true_value)
}

output_figure <- function(figure, folder, name, height, width){
  png(paste0(burdenEM_ROOT, folder,'/figures/', name, '.png'), width=width, height=height, units = 'in', res = 300)
  print(figure)
  dev.off()
}

powerfn <- function(x, a, n) {
  pnorm(sqrt(n) * x - sqrt(a)) + pnorm(-sqrt(n) * x - sqrt(a))
}

scaleFUN <- function(x) sprintf("%.2f", x)

percent_labels <- function(x) {
  paste0(format(round(x * 100), nsmall = 0), "%")
}

scientific_10 <- function(x) {
  # Create scientific notation, remove "e", replace with "10^", and clean up the "1 x"
  formatted <- gsub("e\\+?", " %*% 10^", scales::scientific_format()(x))  # Replace "e" with "10^" and remove "+"
  formatted <- gsub("^1 \\\\%\\*% ", "", formatted)  # Remove the leading "1 %*% " for numbers like 1 x 10^
  parse(text = formatted)
}

load_genebass_bhr_paper_results <- function(bhr_threshold=0.05){
  bhr_paper <- read_csv(Genebass_BHR_paper_PATH) %>%
    mutate(phenotype_key = str_replace(str_replace(phenotype_key, 'NA', '_NA'), '20002', '20002_')) %>%
    mutate(AF_bin = case_when(
      endsWith(summary_statistic, '_group1') ~ '0e+00_1e-05',
      endsWith(summary_statistic, '_group2') ~ '1e-05_1e-04',
      endsWith(summary_statistic, '_group3') ~ '1e-04_1e-03',
    ),
    annotation = case_when(
      grepl('pLoF', summary_statistic) ~ 'pLoF',
      grepl('missense-benign', summary_statistic) ~ 'missense_benign',
      grepl('missense-notbenign', summary_statistic) ~ 'missense_damaging',
      grepl('synonymous', summary_statistic) ~ 'synonymous',
    )
    ) %>%
    mutate(bhr_z = bhr_h2/bhr_h2_se,
           bhr_p = 2 * (1 - pnorm(abs(bhr_z))),
           bhr_significant = bhr_p < bhr_threshold)
  return(bhr_paper)
}

load_genebass_bhr_results <- function(bhr_threshold = 0.05){
  genebass_bhr <- read_csv(Genebass_BHR_PATH) %>%
    dplyr::select(1:12) %>%
    merge(., genebass_pheno_info %>% dplyr::select(phenotype, phenocode, coding, trait_type, modifier, pheno_sex, description) , by.x = 'phenotype_key', by.y = 'phenotype', all.x=T) %>%
    merge(., genebass_neff, by = colnames(genebass_neff)[1:5]) %>%
    mutate(bhr_z = bhr_h2/bhr_h2_se,
           bhr_p = 2 * (1 - pnorm(abs(bhr_z))),
           bhr_significant = bhr_p < bhr_threshold)
  return(genebass_bhr)
}

load_genebass_BurdenEM_results <- function(n_eff_threshold = 200000){
  genebass_h2 <- read_csv(Genebass_BurdenEM_PATH) %>%
    merge(., genebass_pheno_info , by = 'phenotype', all.x=T) %>%
    merge(., genebass_neff, by = colnames(genebass_neff)[1:5]) %>%
    filter(trait_type == 'continuous'& n_eff > n_eff_threshold) # TODO: edit to include binary
  return(genebass_h2)
}

load_aou_burdenEM_ancestry_data <- function(){
  aou_afr_h2 <- read_csv(get_aou_burdenEM_ancestry_path('afr')) %>%
    mutate(ancestry = 'afr', phenoname = phenotype)
  aou_amr_h2 <- read_csv(get_aou_burdenEM_ancestry_path('amr')) %>%
    mutate(ancestry = 'amr', phenoname = phenotype)
  aou_eur_h2 <- read_csv(get_aou_burdenEM_ancestry_path('eur')) %>%
    mutate(ancestry = 'eur', phenoname = phenotype)

  aou_h2 <- rbind(aou_afr_h2, aou_amr_h2, aou_eur_h2) %>%
    merge(., aou_neff, by = c('ancestry', 'phenoname') ) %>%
    merge(., full_pheno_sum %>% dplyr::select(trait_type, ancestry = pop, n_cases, n_controls, description, phenoname), by = c('ancestry', 'phenoname'))
  return(aou_h2)
}

compute_aou_bhr_h2_results <- function(bhr_data, AF_bins = c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03')){
  bhr_sum <- bhr_data %>%
    filter(AF_bin %in% AF_bins) %>%
    group_by(ancestry, phenotype_key, annotation, N) %>%
    dplyr::summarize(bhr_h2 = sum(bhr_h2),
                     bhr_h2_se = sqrt(sum(bhr_h2_se^2)))
  return(bhr_sum)
}

generate_beta_samples <- function(components, weights, n_draws){
  beta_samples <- unlist(sapply(1:length(components),
                         function(i){components[i] * runif(as.integer(floor(n_draws * weights[i])))}))
  return(beta_samples)
}

n_large_effect_gene <- function(components, weights){
  library(e1071)
  n_genes = 18000

  beta_samples <- generate_beta_samples(components, weights, n_draws=1e6)
  cumh2 <- cumsum(sort(beta_samples^2, decreasing = TRUE))  # Cumulative sum of sorted squared samples
  first <- which(cumh2 >= cumh2[length(cumh2)] / 2)[1]  # Find the first index where cumulative h2 reaches half

  polygenicity_50 <- first / length(beta_samples) * n_genes
  polygenicity_eff  = 3*n_genes*mean(beta_samples^2)^2 / mean(beta_samples^4)
  n_large_gene_1e_1_pos <- sum(abs(beta_samples) > 0.1 & beta_samples > 0)/length(beta_samples) * n_genes
  n_large_gene_5e_2_pos <- sum(abs(beta_samples) > 0.05 & beta_samples > 0)/length(beta_samples) * n_genes
  n_large_gene_2e_2_pos <- sum(abs(beta_samples) > 0.02 & beta_samples > 0)/length(beta_samples) * n_genes
  n_large_gene_1e_2_pos <- sum(abs(beta_samples) > 0.01 & beta_samples > 0)/length(beta_samples) * n_genes
  n_large_gene_1e_1_neg <- sum(abs(beta_samples) > 0.1 & beta_samples < 0)/length(beta_samples) * n_genes
  n_large_gene_5e_2_neg <- sum(abs(beta_samples) > 0.05 & beta_samples < 0)/length(beta_samples) * n_genes
  n_large_gene_2e_2_neg <- sum(abs(beta_samples) > 0.02 & beta_samples < 0)/length(beta_samples) * n_genes
  n_large_gene_1e_2_neg <- sum(abs(beta_samples) > 0.01 & beta_samples < 0)/length(beta_samples) * n_genes

  result_list <- as.integer(c(polygenicity_50, polygenicity_eff,
                              n_large_gene_1e_1_pos, n_large_gene_5e_2_pos, n_large_gene_2e_2_pos, n_large_gene_1e_2_pos,
                              n_large_gene_1e_1_neg, n_large_gene_5e_2_neg, n_large_gene_2e_2_neg, n_large_gene_1e_2_neg))
  return(result_list)
}


