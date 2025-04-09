#install BHR
# devtools::install_github("ajaynadig/bhr", force = TRUE)
packages <- c('ggplot2', 'tidyverse','dplyr','ggrepel','optparse','gridExtra',
              'data.table','grDevices', 'stringr', 'purrr', 'readr', 'ggpubr',
              'bhr', 'rtracklayer', 'GenomicRanges', 'magrittr', 'stringi')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages( p,  repos = c(CRAN = "http://cran.r-project.org") )
  }
}

source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM_old/burdenEM/wlu_test/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas.R')

# gsutil -m cp -n gs://aou_wlu/250k_analysis/burdenEM/genebass/var_files/genebass_*.txt.bgz ~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/data/genebass/var_txt/.
# gsutil -m cp gs://aou_wlu/250k_analysis/burdenEM/genebass/gene_files/genebass_*.txt.bgz ~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/data/genebass/binary_gene_txt/.
# gsutil -m cp -n gs://aou_wlu/250k_analysis/burdenEM/aou/var_files/aou_meta*.txt.bgz ~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/data/aou/var_txt/.
# gsutil -m cp -n gs://aou_wlu/250k_analysis/burdenEM/aou/gene_files/aou_*.txt.bgz ~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/data/aou/binary_gene_txt/.

############## CONSTANTS ##############
BIN_NAMES <- c('singleton', 'doubleton', '0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03', '1e-03_1e-02', '1e-02_1e-01', '2_1e-4')
ANNOTATIONS <- c('pLoF', 'missense_damaging', 'missense_benign', 'synonymous')
gtf_path <- '~/Dropbox (Partners HealthCare)/analysis/lung_function_pleiotropy_yixuan/paper/Homo_sapiens.GRCh38.105.gtf'
gtf_data <- import(gtf_path, format = "gtf")
gtf_genes <- gtf_data[gtf_data$type == "gene"]
gtf_gene_info <- data.frame(
  gene_id = mcols(gtf_genes)$gene_id,
  gene_name = mcols(gtf_genes)$gene_name
)

baseline_model <- read.table('~/Dropbox (Partners HealthCare)/burdenEM/ms_baseline_oe5.txt')
T_baseline <- baseline_model %>%
  mutate(baseline_oe5_total5 = 1 - baseline_oe1_total5 - baseline_oe2_total5 - baseline_oe3_total5 - baseline_oe4_total5) %>% # Constraint bins
  merge(., gtf_gene_info, by.x = 'gene', by.y = 'gene_id') %>%
  mutate(gene_id = gene, gene_name = gene_name) %>% select(-gene)

file_path <- function(data_name, pheno, af_bin, annotation, root, ancestry){
  if(af_bin == 'singleton'){
    lower <- 1
    upper <- 2
  }else if(af_bin == 'doubleton'){
    lower <- 2
    upper <- 3
  }else{
    lower <-  format(as.numeric(if_else(startsWith(af_bin, '0'), as.character(0), str_split(af_bin, '_') %>% map_chr(., 1))), scientific = FALSE)
    upper <- format(as.numeric(str_split(af_bin, '_') %>% map_chr(., 2)), scientific = FALSE)
  }
  if(lower == '0.00001') lower = '1e-05'
  if(upper == '0.00001') upper = '1e-05'
  suffix <- paste0('_', pheno, '_', annotation, '_low_', lower, '_high_', upper, '_', af_bin, '.txt.bgz')
  print(suffix)
  if(toupper(data_name) == 'AOU'){
    af_bin <- str_replace(af_bin, 'e+00', '')
    path <- paste0(root, 'aou_', ancestry, suffix)
  }else{
    path <- paste0(root, 'genebass', suffix)
  }
  return(path)
}

read_var_txt <- function(path, input_root=''){
  data <- read_delim(paste0(input_root, path), delim = '\t',
                     col_types = cols(beta = col_double(),
                                      phenotype_key = col_character(),
                                      N = col_integer(),
                                      POS = col_integer(),
                                      AC_cases = col_double(),
                                      AF = col_double(),
                                      AF_total = col_double(),
                                      prevalence = col_double(),
                                      variant_variance = col_double()))
  if(input_root != ''){
    index <- if_else(startsWith('genebass', path), 3, 4)
    af_bin <-str_replace(stri_split_fixed(str = str_split(path, 'high_')[[1]][2], pattern = "_", n = 2)[[1]][2], '.txt.bgz', '')
    print(af_bin)
    data <- data %>% mutate(AF_bin = af_bin)
  }
  return(data)
}

get_freq_interval <- function(freq){
  interval = case_when(
    freq < 1e-5  ~ '0e+00_1e-05',
    freq < 1e-4  ~ '1e-05_1e-04',
    freq < 1e-3  ~ '1e-04_1e-03',
  )
  return(interval)
}

############## Setup Arguments ##############
option_list <- list(
  make_option(c("-a", "--per_allele"), type="logical", default=FALSE,
              help="Whether to run burdenEM with per allele effect sizes (if not, the default is running burdenEM with per SD effect sizes)", metavar="logical"),
  make_option(c("-b", "--n_boot"), type="integer", default=10,
              help="Number of bootstrap to run for BurdenEM", metavar="integer"),
  make_option(c("-c", "--genomewide_correction"), type="logical", default=FALSE,
              help="use genomewide correction in BHR", metavar="character"),
  make_option(c("-d","--data_name"), type="character", default='genebass',
              help="Data source [default= %default]", metavar="character"),
  make_option(c("-g", "--grid_size"), type="integer", default=100,
              help="Grid size to use for BurdenEM", metavar="integer"),
  make_option(c("-i","--input_root"), type="character", default = '~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/data/',
              help="input directory where variant txt files are stored [default= %default]", metavar="character"),
  make_option(c("-k", "--n_cpt"), type="integer", default=15,
              help="Number of mixture components for BurdenEM", metavar="integer"),
  make_option(c("-m", "--n_null"), type="integer", default=1,
              help="Number of null simulations to run for BurdenEM", metavar="integer"),
  make_option(c("-n", "--n_iter"), type="integer", default=1000,
              help="Number of iterations to run for EM in BurdenEM", metavar="integer"),
  make_option(c("-o", "--output_root"), type="character", default='~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/',
              help="output directory [default= %default]", metavar="character"),
  make_option(c("-p","--phenotypes"), type="character", default="",
              help="list of phenotypes split by comma [default= %default]", metavar="character"),
  make_option(c("-r","--rerun"), type="logical", default= FALSE,
              help="-whether to overwrite existing files [default= %default]", metavar="logical"),
  make_option(c("-s","--suffix"), type="character", default="",
              help="suffix to the output files [default= %default]", metavar="character"),
  make_option(c("-t","--test"), type="logical", default= FALSE,
              help="-whether to run test [default= %default]", metavar="logical"),
  make_option(c("-x","--skip_bhr"), type="logical", default= FALSE,
              help="-whether to skip running BHR [default= %default]", metavar="logical"),
  make_option(c("-y","--skip_burden_test"), type="logical", default= FALSE,
              help="-whether to skip running Burden test [default= %default]", metavar="logical")
);


############## Load Arguments ##############
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser, positional_arguments=0);

print(str(opt))
per_allele <- opt$options$per_allele
n_boot <- opt$options$n_boot
genomewide_correction <- opt$options$genomewide_correction
data_name <- opt$options$data_name
grid_size <- opt$options$grid_size
input_root <- opt$options$input_root
n_cpt <- opt$options$n_cpt
n_null <- opt$options$n_null
n_iter <- opt$options$n_iter
output_root <- opt$options$output_root
phenotypes <- opt$options$phenotypes
overwrite <- opt$options$rerun
suffix <- opt$options$suffix
test <- opt$options$test
skip_bhr <- opt$options$skip_bhr
skip_burden_test <- opt$options$skip_burden_test


############## Process Arguments ##############
ancestries <- c('eur')
if(tolower(data_name) == 'aou') ancestries <- c('amr', 'eur')
if(genomewide_correction) suffix <- paste0('_genomewide_corrected', suffix)
if(test) suffix <- paste0(suffix, '_test')
if(phenotypes != ''){
  phenotypes <- str_split(phenotypes, ',')[[1]]
}else{
  if(tolower(data_name) == 'aou'){
    aou_util_data_path <- '~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/data/'
    phenotypes <- c('blood-pressure-systolic-mean', 'blood-pressure-diastolic-mean', "BMI", "3024929",
                    "3013721", "3027114", "3004501", "3007070", "3022192", "height", "3026910", "3028288",
                    # "495", "401", "250", '250.2', '411', 'EM_239', 'CV_401', '272.1', 'MS_708',
                    # '740', 'EM_202.2', 'EM_202', '250.1', 'EM_202.1', 'CV_404', 'RE_475', '250.2', '411',
                    "random_0.5_continuous_1", "random_0.5_continuous_2", "random_0.5_continuous_3", "random_0.5_continuous_4","random_0.5_continuous_5",
                    'random_0.5_0.001_1', 'random_0.5_0.001_2', 'random_0.5_0.001_3', 'random_0.5_0.001_4', 'random_0.5_0.001_5',
                    'random_0.5_0.01_1', 'random_0.5_0.01_2', 'random_0.5_0.01_3', 'random_0.5_0.01_4', 'random_0.5_0.01_5',
                    'random_0.5_0.1_1', 'random_0.5_0.1_2', 'random_0.5_0.1_3', 'random_0.5_0.1_4', 'random_0.5_0.1_5',
                    'random_0.5_0.2_1', 'random_0.5_0.2_2', 'random_0.5_0.2_3', 'random_0.5_0.2_4', 'random_0.5_0.2_5',
                    'random_0.5_0.5_1', 'random_0.5_0.5_2', 'random_0.5_0.5_3', 'random_0.5_0.5_4', 'random_0.5_0.5_5')
    random_phenotypes <- c("random_0.5_continuous_1", "random_0.5_continuous_2", "random_0.5_continuous_3", "random_0.5_continuous_4","random_0.5_continuous_5",
                    'random_0.5_0.001_1', 'random_0.5_0.001_2', 'random_0.5_0.001_3', 'random_0.5_0.001_4', 'random_0.5_0.001_5',
                    'random_0.5_0.01_1', 'random_0.5_0.01_2', 'random_0.5_0.01_3', 'random_0.5_0.01_4', 'random_0.5_0.01_5',
                    'random_0.5_0.1_1', 'random_0.5_0.1_2', 'random_0.5_0.1_3', 'random_0.5_0.1_4', 'random_0.5_0.1_5',
                    'random_0.5_0.2_1', 'random_0.5_0.2_2', 'random_0.5_0.2_3', 'random_0.5_0.2_4', 'random_0.5_0.2_5',
                    'random_0.5_0.5_1', 'random_0.5_0.5_2', 'random_0.5_0.5_3', 'random_0.5_0.5_4', 'random_0.5_0.5_5')
    # phenotypes <- c('blood-pressure-systolic-mean', 'blood-pressure-diastolic-mean', "BMI", "3024929",
    #                 "3013721", "3027114", "3004501", "3007070", "3022192", "height", "3026910", "3028288",
    #                 "495", "401", "250", '250.2', '411', 'EM_239', 'CV_401', '272.1', 'MS_708',
    #                 '740', 'EM_202.2', 'EM_202', '250.1', 'EM_202.1', 'CV_404', 'RE_475', '250.2', '411')
    phenotypes <- read_delim(paste0(aou_util_data_path, '250k_qc_aou_phenotype_meta_info_250k.txt.bgz'), delim = '\t') %>%
      select(-lambda_gc_exome, -lambda_gc_acaf, -lambda_gc_gene_burden_001) %>%
      filter(trait_type == 'continuous' & !startsWith(phenoname, 'random')) %>%
      select(phenoname) %>%
      unique() %>%
      unlist() %>%
      unname()
    # phenotypes <- unique(c(random_phenotypes, phenotypes))
  }else{
    phenotypes <- c('50_NA', '100022_NA', 'alcohol_intake_custom_NA', '20414_NA', "100017_NA", "21001_NA","21002_NA", '30010_NA',
                    '30020_NA', "30040_NA",  "30060_NA",  "30070_NA",  "30080_NA",  "30100_NA", "30190_NA",
                    "30200_NA", "30220_NA", "30530_NA", "30610_NA", "30620_NA", '30640_NA', "30650_NA",
                    "30670_NA", '30680_NA',"30690_NA", '30730_NA', "30740_NA", '30750_NA', "30760_NA", '30770_NA',
                    '30780_NA', "30860_NA", '30870_NA', '30600_NA', '3062_NA', '3063_NA', "48_NA", "49_NA",
                    '4079_NA', '4080_NA', "5983_NA", "WHR_custom_NA",'189_NA', '20016_NA', '20022_NA', '20023_NA',
                    '20127_NA', '23105_NA', '78_NA','130706_NA', '130708_NA', '131494_NA', '131306_NA', '131286_NA')
    phenotypes <- c('50_NA', "21001_NA", '20127_NA', '30750_NA', '30780_NA')
  }
}

gene_list <- read_delim('~/Downloads/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', delim='\t') %>%
  dplyr::select(gene, cds_length) %>%
  dplyr::arrange(desc(cds_length))

long_gene_list <- gene_list[1:round(nrow(gene_list)*0.2),] %$% gene
short_gene_list <- gene_list[round(nrow(gene_list)*0.2+1):nrow(gene_list),] %$% gene

# BHR
if(!skip_bhr){
  print('-------------Running BHR-------------')
  BIN_NAMES <- c('singleton', 'doubleton', '0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03', '1e-03_1e-02', '2_1e-04')
  # BIN_NAMES <- c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03', '1e-03_1e-02', '2_1e-04')
  BIN_NAMES <- c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03', '1e-03_1e-02')

  for(anc in ancestries){
    full_results <- data.frame()
    burden_scores <- data.frame()
    for(pheno in phenotypes){
      output_path <- paste0(output_root, 'BHR/', data_name, '/', data_name,'_', anc, '_', pheno, suffix, '.csv')
      print(output_path)
      af_bin_names <- BIN_NAMES
      if(file.exists(output_path) & (!overwrite)){
        results_holder <- read_csv(output_path)
        full_results <- rbind(full_results, results_holder)
        current_bin_names <- unique(results_holder$AF_bin)
        if(length(setdiff(BIN_NAMES, current_bin_names)) == 0){
          next
        }else{
          af_bin_names <- setdiff(BIN_NAMES, current_bin_names)
          print(af_bin_names)
        }
      }
      results_holder <- data.frame()
      for (anno in ANNOTATIONS){
        for (i in 1:length(af_bin_names)) {
          # Loading data
          path <- file_path(data_name = data_name, pheno=pheno, af_bin = af_bin_names[i], annotation=anno, root=paste0(input_root, data_name, '/var_txt/'), ancestry=anc)
          if(!file.exists(path)){
            print(paste(data_name, pheno, anno, af_bin_names[i], 'variant txt file does not exist!'))
            print(path)
            next
          }
          print(path)
          if(pheno == '50_NA' & af_bin_names[i] == '2_1e-04') next
          summary_statistics <- read_var_txt(path) %>%
            filter(!is.na(variant_variance)) %>%
            mutate(phenotype_key = pheno)
          print(paste0('n variants left:', nrow(summary_statistics)))
          if(nrow(summary_statistics) == 0){
            print(paste(data_name, pheno, anno, af_bin_names[i], 'variant txt file is empty!'))
            next
          }

          gene_position <- summary_statistics %>%
            group_by(gene) %>%
            dplyr::summarize(gene_position = min(POS))

          summary_statistics$chromosome <- summary_statistics$CHR
          summary_statistics <- summary_statistics %>%
            merge(., gene_position, by = 'gene') %>%
            merge(., gtf_gene_info, by.x = 'gene', by.y = if_else(data_name == 'genebass', 'gene_name', 'gene_id'))
          if(data_name == 'genebass') summary_statistics <- summary_statistics %>% mutate(gene_name = gene, gene = gene_id)


          n_bhr_trait <- head(summary_statistics[summary_statistics$phenotype_key == pheno,"N"],1)
          print(c(pheno, n_bhr_trait))

          summary_statistics <- summary_statistics %>%
            filter(gene_name %in% short_gene_list) # 12.10.2024 long short gene analysis
          tmp_burden_score <- summary_statistics %>%
            dplyr::group_by(gene_name) %>%
            dplyr::summarize(
              burden_score = sum(variant_variance, na.rm=T)
            ) %>%
            mutate(
              annotation = anno,
              af_bin = af_bin_names[i],
              pheno = pheno
            )
          burden_scores <- rbind(burden_scores, tmp_burden_score)

          output = BHR(mode = "univariate",
                       trait1_sumstats = summary_statistics[summary_statistics$phenotype_key == pheno,],
                       annotations = list(baseline_model),
                       genomewide_correction = genomewide_correction,
                       slope_correction = 4.55087151/n_bhr_trait,
                       custom_variant_variances = T
          )

          tmp_results <- c(anc,
                           pheno,
                           anno,
                           af_bin_names[i],
                           n_bhr_trait,
                           output$mixed_model$heritabilities[1,ncol(baseline_model)],
                           output$mixed_model$heritabilities[2,ncol(baseline_model)],
                           output$mixed_model$enrichments[1,1],
                           output$mixed_model$enrichments[2,1],
                           output$mixed_model$enrichments[1,2],
                           output$mixed_model$enrichments[2,2],
                           output$mixed_model$enrichments[1,3],
                           output$mixed_model$enrichments[2,3],
                           output$mixed_model$enrichments[1,4],
                           output$mixed_model$enrichments[2,4],
                           ((1-sum(output$mixed_model$fractions[1,]))/(1-sum(output$mixed_model$fraction_burden_score))),
                           output$significant_genes$number_significant_genes,
                           output$significant_genes$fraction_burdenh2_significant,
                           output$significant_genes$fraction_burdenh2_significant_se,
                           output$qc$intercept,
                           output$qc$intercept_se,
                           output$qc$attenuation_ratio,
                           output$qc$attenuation_ratio_se,
                           output$qc$lambda_gc,
                           output$qc$lambda_gc_se,
                           output$qc$mu_genome)
          results_holder <- rbind(results_holder, tmp_results)

          print(paste0("Finished BHR estimate for ",pheno," in summary statistic group ",anno, ' ', af_bin_names[i], '(', toupper(anc), ')'))
          # if(output$significant_genes$number_significant_genes>0){
          #   print('Outputting significant gene table...')
          #   sig_gene_table <- output$significant_genes$sig_table %>%
          #     merge(., gtf_gene_info, by.x = 'gene', by.y = 'gene_id')
          #   write_csv(sig_gene_table, paste0(output_root, 'BHR/', data_name, '/sig_gene/', data_name,'_', anc, '_', af_bin_names[i], '_', anno, '_test_', pheno, suffix, '.csv'))
          # }
        }
        if(test) break
      }
      if(nrow(results_holder) == 0) next
      results_holder[,6:ncol(results_holder)] <- sapply(results_holder[,6:ncol(results_holder)], as.numeric)
      colnames(results_holder) <- c("ancestry", "phenotype_key", "annotation", "AF_bin", "N",
                                    "bhr_h2", "bhr_h2_se",
                                    "bhr_enrichment_oe1", "bhr_enrichment_oe1_se",
                                    "bhr_enrichment_oe2", "bhr_enrichment_oe2_se",
                                    "bhr_enrichment_oe3", "bhr_enrichment_oe3_se",
                                    "bhr_enrichment_oe4", "bhr_enrichment_oe4_se",
                                    "bhr_enrichment_oe5",
                                    "n_significant_genes","fraction_h2_significant_genes", "fraction_h2_significant_genes_se",
                                    "intercept", "intercept_se", "attenuation_ratio", "attenuation_ratio_se", "lambda_gc", "lambda_gc_se", "mu_genome")
      write_csv(results_holder, paste0(output_root, 'BHR/', data_name, '/', data_name,'_', anc, '_', pheno, suffix, '.csv'))
      full_results <- rbind(full_results, results_holder)
      if(test) break
    }
    write_csv(full_results, paste0(output_root, 'BHR/', data_name, '_', anc, suffix, '_tmp_full.csv'))
    write_csv(burden_scores, paste0(output_root, 'BHR/', data_name, '_', anc, suffix, '_tmp_full_burden_score.csv'))
    if(test) break
  }
}


# Burden test
if(!skip_burden_test){
  print('-------------Running Burden Test-------------')
  for(anc in ancestries){
    data_label <- if_else(data_name == 'genebass', data_name, paste0(tolower(data_name), '_', anc))
    bhr_intercept <- read_csv(paste0(output_root, 'BHR/', data_name, '_', anc, suffix, '_tmp_full.csv')) %>%
      select('ancestry', 'phenotype_key', 'annotation', 'AF_bin', 'bhr_h2', 'intercept')
    af_bin_names <-  c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03')
    if(data_name == 'aou' & anc %in% c('afr', 'amr')) af_bin_names <-  c('singleton', '2_1e-04', '1e-04_1e-03')

    for(pheno in phenotypes) {
      if(!overwrite & file.exists(paste0(output_root, 'Burden_test/', data_name, '_', anc, '/', data_name, '_', anc, '_', pheno, suffix, '_burden_test.csv'))){
        next
      }else{
        n_var <- data.frame()
        gene_results <- data.frame()
      }
      for (anno in ANNOTATIONS){
        # Load data for each MAF bin, adding just one extra column containing the BHR intercept for that trait-MAF bin pair
        files <- list.files(paste0(input_root, data_name, '/var_txt/'), pattern = paste0("^", data_label, "_", pheno, "_", anno,".*_(", paste(gsub("([+-])", "\\\\\\1", af_bin_names), collapse = "|"), ")\\.txt\\.bgz$"))
        var_data <- files %>%
          map_dfr(read_var_txt, input_root=paste0(input_root, data_name, '/var_txt/')) %>%
          mutate(ancestry = anc, phenotype_key = pheno, annotation=anno)
        if(nrow(var_data) == 0) next
        print(paste0('[Burden test] ', data_name, ' ', toupper(anc), ' ', pheno, ' ', anno))

        gene_ac_cases <- var_data %>%
          group_by(gene) %>%
          dplyr::summarize(AC_cases = as.integer(sum(AC_cases, na.rm = T)),
                           CAF = sum(AF, na.rm=T))

        var_data <- var_data %>%
          merge(., bhr_intercept, by =c("ancestry","phenotype_key", "annotation", "AF_bin"), all.x=T) %>%
          select(gene, ancestry, annotation, phenotype_key, intercept, variant_variance, beta, N, prevalence, trait_type) %>%
          filter(complete.cases(.))
        if(nrow(var_data) == 0) next

        tmp_n_var <- data.frame(ancestry = anc, phenotype = pheno, annotation = anno, N_var = nrow(var_data))
        n_var <- rbind(n_var, tmp_n_var)

        # Group rows of table by gene (as opposed to by MAF bin)
        genes <- unique(var_data$gene)
        noGenes <- length(genes)
        nn <- mean(var_data$N, na.rm=T)
        mean_intercept_nn <- mean(var_data$intercept) * nn

        tmp_gene_results <- var_data %>%
          # X'*y/n where X is the n x m genotype matrix and y is the n x 1 phenotype vector
          mutate(Xty = var_data$variant_variance * var_data$beta) %>%
          group_by(gene, annotation, ancestry, phenotype_key, N, prevalence, trait_type) %>%
          dplyr::summarize(
            XWty = sum(Xty), # (X*W)'*y/n where W is a variants by genes matrix of burden weights
            burdenScore = sum(variant_variance),
            intercept_vector = sum(intercept * variant_variance)/burdenScore # Convert variant-wise intercept into a gene-wise intercept
            ) %>%
          mutate(
            gamma_perAllele = XWty / burdenScore, # Burden effect size estimates in per-allele units
            gamma_perSD = gamma_perAllele * sqrt(burdenScore) # Burden effect size estimates in per-sd units
            ) %>%
          merge(., gene_ac_cases, by = 'gene')
        tmp_gene_results <- tmp_gene_results %>% mutate(mean_variant_intercept = mean_intercept_nn, mean_n = nn)
        gene_results <- rbind(gene_results, tmp_gene_results)
      }
      if(nrow(gene_results) == 0) next
      gene_results <- merge(gene_results, T_baseline, by.x = 'gene', by.y = if_else(data_name == 'genebass', 'gene_name', 'gene_id'))
      write_csv(gene_results, paste0(output_root, 'Burden_test/', data_name, '_', anc, '/', data_name, '_', anc, '_', pheno, suffix, '_burden_test.csv'))
      write_csv(n_var, paste0(output_root, 'Burden_test/', data_name, '_', anc, '/', data_name, '_', anc, '_', pheno, suffix, '_burden_test_n_var.csv'))
    }
   }
}


# BurdenEM
print('-------------Running BurdenEM-------------')
for(anc in ancestries){
  if(!overwrite & file.exists(paste0(output_root, 'BurdenEM/', data_name, '_', anc,'_iter_', n_iter, suffix, '_full_h2_result_tmp_burdenEM.csv'))){
    print('Reading previous BurdenEM results...')
    h2_results <- read_csv(paste0(output_root, 'BurdenEM/', data_name, '_', anc,'_iter_', n_iter, suffix, '_full_h2_result_tmp_burdenEM.csv'))
    weight_results <- read_csv(paste0(output_root, 'BurdenEM/', data_name, '_', anc,'_iter_', n_iter, suffix, '_weight_result_tmp_burdenEM.csv'))
    real_qq_results <- read_csv(paste0(output_root, 'BurdenEM/', data_name, '_', anc,'_iter_', n_iter, suffix, '_real_qq_result_tmp_burdenEM.csv'))
    est_qq_results <- read_csv(paste0(output_root, 'BurdenEM/', data_name, '_', anc,'_iter_', n_iter, suffix, '_est_qq_result_tmp_burdenEM.csv'))
    phenotypes_exist <- unique(h2_results$phenotype)
    phenos_to_run <- phenotypes[!(phenotypes %in% phenotypes_exist)]
    colnames(h2_results) <- NULL
  }else{
    h2_results <- data.frame()
    weight_results <- data.frame()
    real_qq_results <- data.frame()
    est_qq_results <- data.frame()
    phenos_to_run <- phenotypes
  }

  for(pheno in phenos_to_run){
    print(pheno)
    burden_test_path <- paste0(output_root, 'Burden_test/', data_name, '_', anc, '/', data_name, '_', anc, '_', pheno, '_burden_test.csv')
    if(!file.exists(burden_test_path)) next
    gene_data <- read_csv(burden_test_path) %>%
      filter(gene %in% short_gene_list) %>% # 12.03.2024
      select(gene, annotation, ancestry, phenotype_key, N, prevalence, trait_type, burdenScore, intercept_vector, gamma_perAllele,
             gamma_perSD, AC_cases, CAF, baseline_oe1_total5, baseline_oe2_total5, baseline_oe3_total5, baseline_oe4_total5, baseline_oe5_total5,
             mean_variant_intercept, mean_n) %>%
      mutate(trait_type = 'continuous') %>%
      filter(complete.cases(.))

    trait_type <- gene_data %>% filter(phenotype_key == pheno) %$% trait_type %>% unique() %>% unlist()
    # if(trait_type == 'continuous') next
    trait_type = 'continuous'
    if(per_allele){
      gene_data$effect_estimate <- gene_data$gamma_perAllele
      gene_data$effect_se <- sqrt(abs(gene_data$intercept_vector/gene_data$burdenScore))
      suffix <- paste0(suffix, '_per_allele')
    }else{
      gene_data$effect_estimate <- gene_data$gamma_perSD
      gene_data$effect_se <- sqrt(abs(gene_data$intercept_vector))
    }

    for(anno in ANNOTATIONS){
      sub_data <- gene_data %>%
        filter(phenotype_key == pheno & annotation == anno)
      if(nrow(sub_data) == 0) next
      print(paste0('\n---------------', toupper(anc), ':',  pheno, ' (', anno, ') ---------------'))
      input_data <- compute_true_expected_rvas(sub_data)
      input_data$trait_type <- if_else(trait_type=='continuous', 'continuous', 'binary')

      if(trait_type != 'continuous'){
        # pheno_label <- if_else(pheno == 'alcoholintake_customNA', 'alcohol_intake_custom_NA', str_replace(str_replace(pheno, 'NA', '_NA'), '20002', '20002_'))
        missense_label <- if_else(data_name == 'aou', 'missenseLC', 'missense|LC')
        anno_label <- if_else(grepl('missense', anno), missense_label, anno)
        gene_file <- read_delim(paste0(input_root, data_name, '/binary_gene_txt/', data_label,'_gene_', pheno,'_', anno_label,'.txt.bgz'), delim = '\t') %>%
          select(gene_symbol, gene_id, p_value = Pvalue_Burden) %>%
          mutate(phenotype_key = pheno)
        input_data <- input_data %>%
          merge(., gene_file, by.x = c('gene', 'phenotype_key'), by.y = c(if_else(data_name == 'genebass', 'gene_symbol', 'gene_id'), 'phenotype_key')) %>%
          filter(complete.cases(.))
        prevalence <- input_data %$%
          prevalence %>%
          unique(.)
        cpt_lower <- min(log10(prevalence), min(input_data$effect_estimate))
        cpt_upper <- max(log10(1/prevalence), max(input_data$effect_estimate))
        component_endpoints <- seq(cpt_lower, cpt_upper, length.out = n_cpt)
      }else{
        lower_bound_5 = unlist(quantile(input_data$effect_estimate, probs = c(0.001)))
        upper_bound_95 = unlist(quantile(input_data$effect_estimate, probs = c(0.999)))
        bound = max(abs(min(input_data$effect_estimate)), abs(max(input_data$effect_estimate)))
        component_endpoints1 = c(-bound,
                                 seq(lower_bound_5, upper_bound_95, length.out = (n_cpt-2)),
                                 bound)
        bound = max(abs(min(input_data$effect_estimate)), abs(max(input_data$effect_estimate)))
        component_endpoints2 = seq(-bound, bound, length.out = n_cpt)
        component_endpoints2[which(component_endpoints2 == 0)] = 1e-300
        component_endpoints <- unique(c(component_endpoints1, component_endpoints2))
      }
      input_data$per_allele_factor <- 1
      if(per_allele){
        component_endpoints <- component_endpoints / sqrt(mean(input_data$burdenScore))
        input_data$per_allele_factor <- mean(input_data$burdenScore)
      }



      features <- data.matrix(input_data %>%
                                ungroup(.) %>%
                                select(baseline_oe1_total5, baseline_oe2_total5, baseline_oe3_total5, baseline_oe4_total5, baseline_oe5_total5))
      model <- burdenEM_rvas(input_data=input_data,
                             features = features,
                             component_endpoints = component_endpoints,
                             no_cpts = length(component_endpoints),
                             grid_size = grid_size,
                             heritability_est = TRUE,
                             polygenicity_est = TRUE,
                             num_iter =n_iter,
                             prevalence = NULL,
                             bootstrap = T,
                             n_boot = n_boot,
                             null_sim = F,
                             n_null = n_null,
                             return_likelihood = FALSE,
                             qq_plot=TRUE,
                             estimate_posteriors = FALSE)
      print(model$heritability_output)
      real_qq_results <- rbind(real_qq_results, input_data %>% select(expected, observed) %>% mutate(phenotype = pheno, annotation = anno))
      est_qq_results <- rbind(est_qq_results, model$qq_data %>% select(expected, observed) %>% mutate(phenotype = pheno, annotation = anno))
      tmp_h2_results <- data.frame(t(c(model$heritability_output$total_h2,
                                       model$heritability_output$heritability_CI,
                                       model$heritability_output$prop_positive_h2,
                                       model$heritability_output$prop_negative_h2,
                                       model$heritability_output$annot_h2,
                                       model$heritability_output$annot_h2_CI[1,],
                                       model$heritability_output$annot_h2_CI[2,],
                                       model$heritability_output$frac_h2,
                                       model$heritability_output$frach2_CI[1,],
                                       model$heritability_output$frach2_CI[2,],
                                       model$heritability_output$frac_expected,
                                       model$heritability_output$enrichment,
                                       model$heritability_output$enrich_CI[1,],
                                       model$heritability_output$enrich_CI[2,],
                                       model$polygenicity$polygenicity_eff,
                                       model$polygenicity$polygenicity_50,
                                       model$polygenicity$n_large_gene_2e_2,
                                       model$polygenicity$n_large_gene_1e_2,
                                       pheno, anno
                                       )))

      tmp_weights <- cbind(t(model$delta), components = model$component_endpoints) %>%
        as.data.frame(.) %>%
        mutate(phenotype = pheno, annotation = anno)
      h2_results <- rbind(h2_results, tmp_h2_results)
      weight_results <- rbind(weight_results, tmp_weights)
    }
    # print(h2_results)
  }
  colnames(h2_results) <- c('total_h2', 'total_h2_CI_lower', 'total_h2_CI_upper',
                            'prop_positive_h2', 'prop_negative_h2',
                            paste0('h2_', 1:ncol(features)),
                            paste0('h2_CI_lower_', 1:ncol(features)),
                            paste0('h2_CI_upper_', 1:ncol(features)),
                            paste0('prop_h2_', 1:ncol(features)),
                            paste0('prop_h2_CI_lower_', 1:ncol(features)),
                            paste0('prop_h2_CI_upper_', 1:ncol(features)),
                            paste0('prop_expected_', 1:ncol(features)),
                            paste0('enrichment_', 1:ncol(features)),
                            paste0('enrichment_CI_lower_', 1:ncol(features)),
                            paste0('enrichment_CI_upper_', 1:ncol(features)),
                            'polygenicity_eff', 'polygenicity_50',
                            'n_large_gene_2e_2', 'n_large_gene_1e_2',
                            'phenotype', 'annotation')
  View(h2_results)
  write_csv(h2_results, paste0(output_root, 'BurdenEM/', data_name, '_', anc, '_n_boot_', n_boot,'_grid_size_', grid_size,'_iter_', n_iter, suffix, '_full_h2_result_tmp_burdenEM_union_cpt.csv'))
  write_csv(weight_results, paste0(output_root, 'BurdenEM/', data_name, '_', anc, '_n_boot_', n_boot,'_grid_size_', grid_size,'_iter_', n_iter, suffix, '_weight_result_tmp_burdenEM_union_cpt.csv'))
  write_csv(real_qq_results, paste0(output_root, 'BurdenEM/', data_name, '_', anc, '_n_boot_', n_boot,'_grid_size_', grid_size,'_iter_', n_iter, suffix, '_real_qq_result_tmp_burdenEM_union_cpt.csv'))
  write_csv(est_qq_results, paste0(output_root, 'BurdenEM/', data_name, '_', anc, '_n_boot_', n_boot,'_grid_size_', grid_size,'_iter_', n_iter, suffix, '_est_qq_result_tmp_burdenEM_union_cpt.csv'))
}


# Rscript full_burden_pipeline.R -x T -y T -d genebass -r T


