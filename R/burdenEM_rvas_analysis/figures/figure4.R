source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas_analysis/utils.R')

phenos_to_select <- c('Height', 'BMI', 'HbA1c', 'Neuroticism')
names(phenos_to_select) <-  c('50_NA', '21001_NA', '30750_NA', '20127_NA')
ntpr_colors <- c('#887a7e', '#8f89bb', '#bea0b6', '#c08e82')
names(ntpr_colors) <- phenos_to_select[c('50_NA', '21001_NA', '30750_NA', '20127_NA')]

sample_size_factor <- 2^(-1:20)
threshold_factor <- 2^(-4:5)

get_ntpr_data <- function(
    Burden_test_PATH_func=get_genebass_burden_test_path,
    Weight_DATA_PATH=Genebass_weight_PATH,
    NTPR_DATA_PATH=Genebass_NTPR_PATH,
    sample_size_factor = sample_size_factor,
    threshold_factor = threshold_factor,
    overwrite=F){
  if((!file.exists(NTPR_DATA_PATH)) | overwrite){
    weight_data <- read_csv(Weight_DATA_PATH)
    ntpr_result <- data.frame()
    ntpr_threshold_result <- data.frame()
    for(pheno in unique(weight_data$phenotype)){
      print(pheno)
      gene_data <- read_csv(Burden_test_PATH_func(pheno))
      for(anno in ANNOTATIONS){
        print(anno)
        tmp_gene_data <- gene_data %>% filter(annotation == anno)
        mean_n <- unique(tmp_gene_data$mean_n)
        mean_variant_intercept <- unique(tmp_gene_data$mean_variant_intercept)
        tmp_weight_data <- weight_data %>%
          filter(annotation == anno & phenotype == pheno)
        components <- unlist(tmp_weight_data$components)
        weights_agg <-unlist(rowMeans(tmp_weight_data[, 1:5]))
        result <- NTPR(bb = components, pp = weights_agg, nn = mean_n/mean_variant_intercept*sample_size_factor, alpha = 0.05/nrow(tmp_gene_data))
        threshold_result <- NTPR(bb = components, pp = weights_agg, nn = mean_n/mean_variant_intercept, alpha = 0.05/nrow(tmp_gene_data)*threshold_factor)
        tmp_ntpr_result <- data.frame(total_ntpr = c(result$total_ntpr),
                                      total_h2gwas = c(result$total_h2gwas),
                                      total_power = c(result$total_power),
                                      total_power_pos = c(result$total_power_pos),
                                      total_power_neg = c(result$total_power_neg),
                                      sample_size = c(mean_n*sample_size_factor)) %>%
          mutate(annotation = anno, phenotype = pheno, N = mean_n, n_genes = nrow(tmp_gene_data))
        ntpr_result <- rbind(ntpr_result, tmp_ntpr_result)
        tmp_ntpr_threshold_result <- data.frame(total_ntpr = c(result$total_ntpr),
                                                total_power = c(result$total_power),
                                                sample_size = c(mean_n*sample_size_factor)) %>%
          mutate(annotation = anno, phenotype = pheno, N = mean_n, n_genes = nrow(tmp_gene_data))
        ntpr_threshold_result <- rbind(ntpr_threshold_result, tmp_ntpr_threshold_result)
      }
    }
    write_csv(ntpr_result, NTPR_DATA_PATH)
    write_csv(ntpr_threshold_result, str_replace(NTPR_DATA_PATH, 'ntpr', 'ntpr_threshold'))
  }
  ntpr_result <- read_csv(NTPR_DATA_PATH)
  return(ntpr_result)
}

load_n_large_effect_gene_data <- function(weight_data_path=Genebass_weight_PATH){
  n_large_gene_path <- str_replace(weight_data_path, 'weight', 'n_large_gene')
  if(!file.exists(n_large_gene_path)){
    weight_data <- read_csv(weight_data_path)
    weight_results <- weight_results %>%
      mutate(weight = rowSums(weight_results[, 1:5])/5)
    phenotypes_to_run <- unique(weight_results$phenotype)
    n_large_gene_result <- data.frame()
    for(anno in ANNOTATIONS){
      for(pheno in phenotypes_to_run){
        sub_weight <- weight_results %>% filter(phenotype == pheno & annotation == anno)
        component_endpoints <- sub_weight$components
        weights <- sub_weight$weight
        tmp_results <- c(pheno, anno, n_large_effect_gene(components=component_endpoints, weights=weights))
        n_large_gene_result <- rbind(n_large_gene_result, tmp_results)
      }
    }
    colnames(n_large_gene_result) <- c('phenotype', 'annotation','polygenicity_50', 'polygenicity_eff',
                                       'n_large_gene_1e_1_pos', 'n_large_gene_5e_2_pos', 'n_large_gene_2e_2_pos', 'n_large_gene_1e_2_pos',
                                       'n_large_gene_1e_1_neg', 'n_large_gene_5e_2_neg', 'n_large_gene_2e_2_neg', 'n_large_gene_1e_2_neg')
    n_large_gene_result[, 3:12] <- sapply(n_large_gene_result[, 3:12], as.numeric)
    n_large_gene_result <- n_large_gene_result %>%
      merge(., genebass_pheno_info %>% dplyr::select(phenotype, phenocode, coding, trait_type, modifier, pheno_sex, description, N) , by = 'phenotype', all.x=T)
    write_csv(n_large_gene_result, n_large_gene_path)
  }
  n_large_gene_data <- read_csv(n_large_gene_path)
  return(n_large_gene_data)
}

load_long_short_gene_data <- function(){
  process_long_short_burden_score <- function(type){
    gene_burdenscore <- read_csv(paste0(burdenEM_RESULT_ROOT, 'BHR/genebass_eur_', type, '_genes_tmp_full_burden_score.csv')) %>%
      dplyr::group_by(gene_name, annotation, pheno) %>%
      dplyr::summarize(burdenScore = sum(burden_score, na.rm=T)) %>%
      mutate(label = if_else(type == 'long', 'Long Genes (top 20%)', 'Short Genes (rest 80%)'))
    return(gene_burdenscore)
  }
  gene_burden_score <- rbind(process_long_short_burden_score('long'),
                             process_long_short_burden_score('short')) %>%
    group_by(annotation, pheno, label) %>%
    dplyr::summarize(total_burden_score = sum(burdenScore))

  long_gene_bhr <- read_csv(paste0(burdenEM_RESULT_ROOT, 'BHR/genebass_eur_long_genes_tmp_full.csv'))
  short_gene_bhr <- read_csv(paste0(burdenEM_RESULT_ROOT, 'BHR/genebass_eur_short_genes_tmp_full.csv'))
  long_gene_burdenEM <- read_csv(paste0(burdenEM_RESULT_ROOT, 'BurdenEM/genebass_eur_n_boot_10_grid_size_100_iter_1000_full_h2_result_tmp_burdenEM_union_cpt_long_gene_kept.csv'))
  short_gene_burdenEM <- read_csv(paste0(burdenEM_RESULT_ROOT, 'BurdenEM/genebass_eur_n_boot_10_grid_size_100_iter_1000_full_h2_result_tmp_burdenEM_union_cpt_long_gene_removed.csv'))

  long_short_bhr <- rbind(
    long_gene_bhr %>% mutate(label = 'Long Genes (top 20%)'),
    short_gene_bhr %>% mutate(label = 'Short Genes (rest 80%)')
    ) %>%
      filter(AF_bin %in% c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03')) %>%
      group_by(ancestry, phenotype_key, annotation, N, label) %>%
      dplyr::summarize(bhr_h2 = sum(bhr_h2),
                       bhr_h2_se = sqrt(sum(bhr_h2_se^2))) %>%
      mutate(model = 'BHR') %>%
      ungroup() %>%
      dplyr::select(pheno=phenotype_key, annotation, h2 = bhr_h2, label, model) %>%
      group_by(pheno, annotation,  model) %>%
      mutate(h2_all = sum(h2),
             prop_h2 = h2/h2_all)

  long_short_burdenEM <- rbind(
    long_gene_burdenEM %>% mutate(label = 'Long Genes (top 20%)'),
    short_gene_burdenEM %>% mutate(label = 'Short Genes (rest 80%)')
  ) %>%
    mutate(model = 'BurdenEM') %>%
    dplyr::select(pheno=phenotype, annotation, h2 = total_h2, label, model) %>%
    group_by(pheno, annotation, model) %>%
    mutate(h2_all = sum(h2),
           prop_h2 = h2/h2_all) %>%
    rbind(., long_short_bhr) %>%
    merge(., gene_burden_score, by = c('annotation', 'label', 'pheno')) %>%
    mutate(h2_per_burdenscore = h2/total_burden_score)

  return(long_short_burdenEM)
}

plot_figure4a <- function(filter = T, height=3, width=5, save=F, name='figure4a'){
  ntpr_result <- get_ntpr_data()
  if(filter) ntpr_result <- ntpr_result %>% filter(annotation == 'pLoF' & phenotype %in% names(phenos_to_select))
  figure <- ntpr_result %>%
    mutate(sample_size_factor = sample_size/N,
           phenotype = if_else(phenotype %in% names(phenos_to_select), phenos_to_select[phenotype], phenotype),
           annotation = factor(annotation, levels = annotation_types2)) %>%
    ggplot +
    aes(x = sample_size, y = total_h2gwas, color = phenotype, group=phenotype) +
    geom_point(size = 1) + geom_line() +
    geom_hline(yintercept = c(0.5, 1), lty = 2) +
    geom_vline(xintercept = c(500000, 8181170405), color='#8b733d') +
    geom_text(aes(x=500000-150000, label="UK Biobank", y=0.7), color='#8b733d', angle=90, size = 3) +
    geom_text(aes(x=8181170405-5000000000, label="World\nPopulation", y=0.7), color='#8b733d', angle=90, size = 3) +
    scale_x_log10(labels = scientific_10) +  # Log scale for x-axis with limits
    scale_y_continuous(breaks = seq(0,1, 0.25), label = percent) +
    labs(x = 'Sample size', y = NULL, title =  expression(Proportion~of~h[Burden]^2), color = NULL) + themes +
    scale_color_manual(values = ntpr_colors) +
    scale_fill_manual(values = ntpr_colors) +
    theme(legend.position = 'top') +
    facet_grid(~annotation, labeller = label_type)
  if(filter){
    figure <- figure +
      theme(strip.background = element_blank(),
            strip.text = element_blank())
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure4b <- function(filter = T, height=3, width=5, save=F, name='figure4b'){
  ntpr_result <- get_ntpr_data()
  if(filter) ntpr_result <- ntpr_result %>% filter(annotation == 'pLoF' & phenotype %in% names(phenos_to_select))
  figure <- ntpr_result %>%
    mutate(sample_size_factor = sample_size/N,
           phenotype = if_else(phenotype %in% names(phenos_to_select), phenos_to_select[phenotype], phenotype),
           annotation = factor(annotation, levels = annotation_types2)) %>%
    ggplot + aes(x = sample_size, y = n_genes*total_power, color = phenotype, group=phenotype) +
    geom_point(size = 1) + geom_line() +
    scale_x_log10(labels = scientific_10) +
    scale_y_log10(label = comma) +
    geom_vline(xintercept = c(500000, 8181170405), color='#8b733d') +
    geom_text(aes(x=500000-150000, label="UK Biobank", y=1000), color='#8b733d', angle=90, size = 3) +
    geom_text(aes(x=8181170405-5000000000, label="World\nPopulation", y=1000), color='#8b733d', angle=90, size = 3) +
    labs( x = 'Sample size', y = NULL, title = 'N exome-wide significant genes', color = NULL) +
    themes +
    scale_color_manual(values = ntpr_colors) +
    scale_fill_manual(values = ntpr_colors) +
    theme(legend.position = 'top') +
    facet_grid(~annotation, labeller = label_type)
  if(filter){
    figure <- figure +
      theme(strip.background = element_blank(),
            strip.text = element_blank())
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure4c <- function(emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                              '30760_NA', '5983_NA', '30690_NA', '30070_NA'),
                          filter=T, height=3, width=5, save=F, name='figure4c'){
  n_large_gene_result <- load_n_large_effect_gene_data(weight_data_path=Genebass_weight_PATH) %>%
    filter(N > 100000 & trait_type == 'continuous')
  if(filter) n_large_gene_result <- n_large_gene_result %>% filter(annotation == 'pLoF')
  n_large_gene_result <- n_large_gene_result  %>%
    mutate(n_large_gene_1e_2 = n_large_gene_1e_2_pos + n_large_gene_1e_2_neg) %>%
    dplyr::select(1:2, 8, 12:20) %>%
    filter(n_large_gene_1e_2 > 0) %>%
    arrange(desc(n_large_gene_1e_2))
  pheno_levels <- n_large_gene_result$phenotype
  name_format <- c('n_large_gene_1e_2_pos'='Trait-increasing', 'n_large_gene_1e_2_neg'='Trait-decreasing')
  n_large_gene_long <- n_large_gene_result  %>%
    mutate(phenotype = factor(phenotype, levels = pheno_levels)) %>%
    pivot_longer(., cols = c('n_large_gene_1e_2_pos', 'n_large_gene_1e_2_neg'), ) %>%
    mutate(name = name_format[name])

  phenotypes <- unique(n_large_gene_long$phenotype)
  colors_emphasis <- rep('gray70', length(phenotypes))
  names(colors_emphasis) <- phenotypes
  colors_emphasis[emphasis_traits] <- c("#C95693", "#7B03D1", "#D9B407", "#7B8FEA",
                                        "#B20C00", "#F36E04", "#67D720","#227E18","#43BBA7")

  n_large_gene_long[, 'emphasis'] <- as.character(n_large_gene_long[, 'phenotype'])
  n_large_gene_long[!(n_large_gene_long$phenotype %in% emphasis_traits), 'emphasis'] <- ""

  figure <- n_large_gene_long %>%
    ggplot + aes(x = phenotype, y = value, group = interaction(phenotype, name),
                 alpha=name, color = phenotype, fill=phenotype, label=description) +
    geom_col() + themes + theme(axis.text.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                legend.position = 'bottom') +
    scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
    ylim(0, 23) +
    # scale_color_manual(values = ntpr_colors) +
    # scale_fill_manual(values = ntpr_colors) +
    scale_color_manual(values = colors_emphasis) +
    scale_fill_manual(values = colors_emphasis) +
    geom_text(data = n_large_gene_long %>% filter(name == 'Trait-increasing'), aes(y=n_large_gene_1e_2 +1),
                    angle=90, alpha =1, hjust =0, vjust = 0.5, size = 1.5) +
    labs(x = '30 phenotypes with at least one large effect gene', y = NULL, title = 'N genes with heritability > 0.01%', color = NULL, fill= NULL) +
    guides(color = 'none', fill='none') +
    guides(
      alpha = guide_legend(nrow = 1, keywidth = 0.6, keyheight = 0.3)
    )
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure4d <- function(filter = T, height=3, width=5, save=F, name='figure4d'){
  ntpr_result <- get_ntpr_data()
  if(filter) ntpr_result <- ntpr_result %>% filter(annotation == 'pLoF' & phenotype %in% names(phenos_to_select))
  figure <- ntpr_result %>%
    mutate(sample_size_factor = sample_size/N,
           phenotype = if_else(phenotype %in% names(phenos_to_select), phenos_to_select[phenotype], phenotype),
           annotation = factor(annotation, levels = annotation_types2)) %>%
    ggplot + aes(x = sample_size, y = 1-total_ntpr, color = phenotype, group=phenotype) +
    geom_point(size = 1) + geom_hline(yintercept = c(0.05, 0.01), lty = 2) +
    geom_line() +
    geom_vline(xintercept = c(500000, 8181170405), color='#8b733d') +
    geom_text(aes(x=500000-150000, label="UK Biobank", y=1e-4), color='#8b733d', angle=90, size = 3) +
    geom_text(aes(x=8181170405-5000000000, label="World\nPopulation", y=1e-4), color='#8b733d', angle=90, size = 3) +
    scale_x_log10(labels = scientific_10) +
    scale_y_log10(labels = scientific_10, breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.05)) +
    labs( x = 'Sample size', y = NULL, color = NULL, title='False discovery rate (FDR)') +
    themes +
    scale_color_manual(values = ntpr_colors) +
    scale_fill_manual(values = ntpr_colors) +
    theme(legend.position = 'top') +
    facet_grid(~annotation, labeller = label_type)
  if(filter){
    figure <- figure +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = 'none')
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure4e <- function(filter = T, height=3, width=7.5, save=F, name='figure4e'){
  long_short_gene <- load_long_short_gene_data()
  if(filter) long_short_gene <- long_short_gene %>% filter(annotation == 'pLoF' & pheno %in% names(phenos_to_select))

  figure <- long_short_gene %>%
    mutate(phenotype = phenos_to_select[pheno]) %>%
    ggplot + aes(x = phenotype, y = prop_h2, color = phenotype,
                 group=interaction(phenotype, annotation),
                 fill=phenotype, alpha = label) +
    geom_col(stat='jitter') + themes + geom_hline(yintercept=0.5, lty=1, lwd = 0.3) +
    scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
    scale_y_continuous(label= percent) +
    scale_color_manual(values = ntpr_colors) +
    scale_fill_manual(values = ntpr_colors) +
    labs(x =NULL, y=expression(Proportion~h^2[Burden]), color = NULL, fill=NULL) +
    facet_grid(annotation~model, labeller = label_type) +
    theme(legend.position = 'top', legend.direction = 'vertical') +
    guides(
      fill = guide_legend(nrow = 1, keywidth = 0.6, keyheight = 0.3),
      color = guide_legend(nrow = 1, keywidth = 0.6, keyheight = 0.3),
      alpha = guide_legend(nrow = 1, keywidth = 0.6, keyheight = 0.3)
    )
  if(filter){
    figure <- figure +
      theme(strip.background.y  = element_blank(),
            strip.text.y = element_blank())
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}


plot_figure4 <- function(height = 8, width = 9, save=T, name='figure4'){
  figure4a <- plot_figure4a(filter = T, height=3.5, width=5, save=T, name='figure4a')
  figure4b <- plot_figure4b(filter = T, height=3.5, width=5, save=T, name='figure4b')
  figure4c <- plot_figure4c(emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                                '30760_NA', '5983_NA', '30690_NA', '30070_NA'),
                            filter=T, height=3, width=5, save=F, name='figure4c')
  figure4d <- plot_figure4d(filter = T, height=3.5, width=5, save=T, name='figure4d')
  figure4e <- plot_figure4e(filter = T, height=3.5, width=7.5, save=F, name='figure4e')

  figure <- ggpubr::ggarrange(ggpubr::ggarrange(figure4a, figure4b,
                                                labels = c('a', 'b'),
                                                font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                                                ncol=2, nrow = 1, common.legend = T, align='h'),
                              ggpubr::ggarrange(figure4c, figure4d ,
                                                labels = c('c', 'd'),
                                                font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                                                ncol=2, nrow = 1),
                              figure4e,
                              labels = c('', '', 'e'),
                              font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                              ncol=1, align = "v", heights = c(2, 2, 2))
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}
