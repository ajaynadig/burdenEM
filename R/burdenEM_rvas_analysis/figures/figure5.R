source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas_analysis/utils.R')

load_genebass_pos_neg_data <- function(){
  genebass_h2 <- load_genebass_BurdenEM_results(n_eff_threshold = 200000)

  genebass_h2_v2 <- genebass_h2 %>%
    mutate(description = str_replace(description, '\\)', '')) %>%
    mutate(description = str_replace(description, '\\(', '')) %>%
    mutate(pos_neg_cat = case_when(
      prop_positive_h2 >= 0.7 ~ 'More Positive',
      prop_positive_h2 <= 0.3 ~ 'More Negative',
      TRUE ~ 'Balanced'
    ))

  selected_phenos <- genebass_h2_v2 %>%
    filter(pos_neg_cat != 'Balanced') %$%
    description %>%
    unique()
  selected_pheno_name <- c('Townsend deprivation index',
                           'Mean time to identify matches',
                           'Neuroticism score', 'BMI', 'RBC distribution width',
                           "Basophill percentage", 'Albumin', 'FVC', 'FEV1',
                           'Cholesterol', 'HDL Cholesterol', 'LDL direct',
                           'Waist circumference', 'Height', 'Alcohol intake',
                           'Waist-to-hip ratio')[-c(1, 9, 10, 13, 16)]
  names(selected_pheno_name) <- selected_phenos[-c(1, 9, 10, 13, 16)]

  pLoF_info <- genebass_h2_v2 %>%
    filter(annotation == 'pLoF' & description %in% names(selected_pheno_name)) %>%
    arrange(total_h2) %>%
    dplyr::select(phenotype, description, pos_neg_cat)

  genebass_h2_v2 <- genebass_h2_v2 %>%
    filter(!is.na(description)) %>%
    filter(annotation %in% c('pLoF', 'missense_damaging'))%>%
    dplyr::select(-pos_neg_cat) %>%
    merge(., pLoF_info, by = c('phenotype', 'description'), all.y = T) %>%
    dplyr::select(phenotype,description, annotation, prop_positive_h2, prop_negative_h2, total_h2, pos_neg_cat) %>%
    mutate(annotation = factor(annotation, levels = annotations))%>%
    pivot_longer(cols = c('prop_positive_h2', 'prop_negative_h2')) %>%
    mutate(name = factor(str_split(name, '_') %>% map_chr(., 2), levels = c('positive', 'negative'), labels = c('Positive effect', 'Negative effect')) ,
           description =  factor(description, levels = unique(pLoF_info$description), labels =selected_pheno_name[unique(pLoF_info$description)]) ,
           pos_neg_cat = as.character(pos_neg_cat)
    ) %>%
    arrange(total_h2)
  return(genebass_h2_v2)
}

get_width_data_archived <- function(weight_data_path=Genebass_weight_PATH, overwrite=F){
  library(boot)
  percentile_fn <- function(data, indices) {
    sample_data <- data[indices]  # Resample data
    c(quantile(sample_data, probs = 0.10),
      quantile(sample_data, probs = 0.90))
  }
  width_path <- str_replace(weight_data_path, 'weight', 'width')
  if(! file.exists(width_path) | overwrite){
    weight_results <- read_csv(weight_data_path)
    weight_results <- weight_results %>%
      mutate(weight = rowSums(weight_results[, 1:5])/5)
    phenotypes_to_run <- unique(weight_results$phenotype)
    width_result <- data.frame(matrix(nrow = 0, ncol = 4))
    for(anno in ANNOTATIONS){
      for(pheno in phenotypes_to_run){
        sub_weight <- weight_results %>% filter(phenotype == pheno & annotation == anno)
        component_endpoints <- sub_weight$components
        weights <- sub_weight$weight
        beta_samples <- generate_beta_samples(component_endpoints, weights, n_draws=1e6)
        boot_res <- boot(beta_samples, statistic = percentile_fn, R = 100)
        ci_10th <- boot.ci(boot_res, index = 1, type = "perc")
        ci_90th <- boot.ci(boot_res, index = 2, type = "perc")
        tmp_result <- c(pheno, anno, quantile(beta_samples, 0.1), quantile(beta_samples, 0.9), ci_10th$perc[4:5], ci_90th$perc[4:5])
        width_result <- rbind(width_result, tmp_result)
      }
    }
    colnames(width_result) <- c('phenotype', 'annotation', paste0('quantile', c(10, 90)), 'CI_10th_lower', 'CI_10th_upper', 'CI_90th_lower', 'CI_90th_upper')
    width_result[, 3:8] <- sapply(width_result[, 3:8], as.numeric)
    width_result <- width_result[, 1:8] %>%
      merge(., genebass_pheno_info %>% dplyr::select(phenotype, phenocode, coding, trait_type, modifier, pheno_sex, description, N) , by.x = 'phenotype', all.x=T)
    write_csv(width_result, width_path)
  }
  width_data <- read_csv(width_path)
  return(width_data)
}


get_width_data <- function(weight_data_path=Genebass_weight_PATH, overwrite=F){
  library(boot)
  percentile_fn <- function(data, indices) {
    sample_data <- data[indices]  # Resample data
    c(quantile(sample_data, probs = 0.10),
      quantile(sample_data, probs = 0.90))
  }
  width_path <- str_replace(weight_data_path, 'weight', 'width_h2')
  if(! file.exists(width_path) | overwrite){
    weight_results <- read_csv(weight_data_path)
    weight_results <- weight_results %>%
      mutate(weight = rowSums(weight_results[, 1:5])/5)
    phenotypes_to_run <- unique(weight_results$phenotype)
    width_result <- data.frame(matrix(nrow = 0, ncol = 4))
    for(anno in ANNOTATIONS){
      for(pheno in phenotypes_to_run){
        sub_weight <- weight_results %>% filter(phenotype == pheno & annotation == anno)
        component_endpoints <- sub_weight$components
        weights <- sub_weight$weight
        beta_samples <- generate_beta_samples(component_endpoints, weights, n_draws=1e7)
        beta_square <- beta_samples^2
        sorted_beta_square <- sort(beta_square)
        cumsum_beta_square <- cumsum(sorted_beta_square)
        total_beta_square <- sum(sort(beta_square))
        threshold_10 <- total_beta_square * 0.10
        threshold_90 <- total_beta_square * 0.90
        first10 <- sqrt(sorted_beta_square[which(cumsum_beta_square >= threshold_10)[1]])
        first90 <- sqrt(sorted_beta_square[which(cumsum_beta_square >= threshold_90)[1]])
        tmp_result <- c(pheno, anno, first10, first90)
        width_result <- rbind(width_result, tmp_result)
      }
    }
    colnames(width_result) <- c('phenotype', 'annotation', paste0('quantile', c(10, 90)))
    width_result[, 3:4] <- sapply(width_result[, 3:4], as.numeric)
    width_result <- width_result[, 1:4] %>%
      merge(., genebass_pheno_info %>% dplyr::select(phenotype, phenocode, coding, trait_type, modifier, pheno_sex, description, N) , by.x = 'phenotype', all.x=T)
    write_csv(width_result, width_path)
  }
  width_data <- read_csv(width_path)
  return(width_data)
}

plot_figure5a <- function(height = 8, width = 15, save=F, name='figure5a'){
  genebass_h2 <- load_genebass_pos_neg_data()
  figure1 <- genebass_h2 %>%
    filter(pos_neg_cat == 'More Negative') %>%
    ggplot + aes(x = description, y = value*total_h2, color = annotation, group=interaction(description, annotation), fill=annotation, alpha = name) +
    geom_col(stat='identity') + themes +
    scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
    scale_y_continuous(label= percent) +
    annotation_color_scale2 +
    annotation_fill_scale2 +
    labs(x =NULL, y=NULL, color = NULL) + coord_flip() +
    facet_grid(pos_neg_cat~annotation, labeller = label_type) + guides(colour = 'none', fill= 'none')
  figure2 <- genebass_h2 %>%
    filter(pos_neg_cat == 'More Positive') %>%
    ggplot + aes(x = description, y = value*total_h2, color = annotation, group=interaction(description, annotation), fill=annotation, alpha = name) +
    geom_col(stat='identity') + themes +
    scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
    scale_y_continuous(label= percent) +
    theme(strip.background.x = element_blank(),
          strip.text.x = element_blank()) +
    annotation_color_scale2 +
    annotation_fill_scale2 +
    labs(x =NULL, y=expression(h[Burden]^2), color = NULL) + coord_flip() +
    facet_grid(pos_neg_cat~annotation, labeller = label_type) + guides(colour = 'none', fill= 'none')
  figure = ggpubr::ggarrange(figure1, figure2, nrow=2, hjust = 0, common.legend = TRUE, align = 'v',
                        font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
  )
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure5b <- function(emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                              '30760_NA', '5983_NA', '30690_NA', '30070_NA'),
                          height = 3, width = 5, save=F, name='figure5b'){
  common_poly <- read_csv(common_var_polygenicity_PATH)
  rare_poly <- genebass_ht <- load_genebass_BurdenEM_results(n_eff_threshold = 100000) %>%
    dplyr::select(phenotype, annotation, polygenicity_50) %>%
    filter(annotation == 'pLoF') %>%
    merge(., common_poly, by = 'phenotype')

  phenotypes <- unique(rare_poly$phenotype)
  colors_emphasis <- rep('gray70', length(phenotypes))
  names(colors_emphasis) <- phenotypes
  colors_emphasis[emphasis_traits] <- c("#C95693", "#7B03D1", "#D9B407", "#7B8FEA",
                                        "#B20C00", "#F36E04", "#67D720","#227E18","#43BBA7")

  rare_poly[, 'emphasis'] <- rare_poly[, 'phenotype']
  rare_poly[, 'emphasis'][!(rare_poly$phenotype %in% emphasis_traits)] <- ""
  rare_poly[, 'Trait'][!(rare_poly$phenotype %in% emphasis_traits)] <- ""

  figure <- rare_poly %>%
    ggplot + aes(x = M_50, y = polygenicity_50, label=Trait, color=emphasis) +
    geom_point() + themes +
    scale_y_continuous(label= comma) +
    scale_x_continuous(label= comma) +
    labs(y = expression(N[genes]~to~explain~"50%"*~h[Burden]^2),
         x = expression(N[SNPs]~to~explain~"50%"*~h[common]^2)) +
    scale_color_manual(values = colors_emphasis) +
    geom_text_repel(max.overlaps = 200, force = 100, size = 2,point.padding = 1,seed = 10332134) +
    labs(title = 'Polygenicity: common vs. rare variants') +
    guides(color = 'none')
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}


plot_figure5c <- function(emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                              '30760_NA', '5983_NA', '30690_NA', '30070_NA'), filter=T,
                          height = 3, width = 5, save=F, name='figure5c'){
  width_data <- get_width_data(weight_data_path=Genebass_weight_PATH, overwrite=F)  %>% filter(N > 100000)
  if(filter) width_data <- width_data %>% filter(trait_type == 'continuous')%>% filter(annotation == 'pLoF')

  phenotypes <- unique(width_data$phenotype)
  colors_emphasis <- rep('gray70', length(phenotypes))
  names(colors_emphasis) <- phenotypes
  colors_emphasis[emphasis_traits] <- c("#C95693", "#7B03D1", "#D9B407", "#7B8FEA",
                                        "#B20C00", "#F36E04", "#67D720","#227E18","#43BBA7")

  width_data[, 'emphasis'] <- width_data[, 'phenotype']
  width_data[!(width_data$phenotype %in% emphasis_traits), 'emphasis']<- ""
  width_data[!(width_data$phenotype %in% emphasis_traits), 'description'] <- ""

  figure <-  width_data %>%
    ggplot + aes(x = quantile10, y=quantile90,
                 # xmin=CI_10th_lower, xmax=CI_10th_upper,
                 # ymin=CI_90th_lower, ymax=CI_90th_upper,
                 color = emphasis, label = description) +
    geom_point(size = 1) +
    # geom_errorbar(lwd = 0.3) + geom_errorbarh(lwd = 0.3) +
    geom_smooth(method = lm, lty=1, color = 'gray', lwd = 0.3, fill = 'gray90') +
    themes + labs(x = 'Tenth percentile', y = '90th percentile', title='Width of effect size') +
    scale_color_manual(values = colors_emphasis) +
    geom_text_repel(max.overlaps = 200, force = 100, size = 2,point.padding = 1,seed = 10332134) +
    scale_y_continuous(label=scientific_10) +
    scale_x_continuous(label=scientific_10) +
    guides(color = 'none')
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}


plot_figure5 <- function(height = 8, width = 9, save=T, name='figure5'){
  figure5a <- plot_figure5a(height = 5, width = 8, save=T, name='figure5a')
  figure5b <- plot_figure5b(emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                                '30760_NA', '5983_NA', '30690_NA', '30070_NA'),
                            height = 3, width = 5, save=F, name='figure5b')
  figure5c <- plot_figure5c(emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                                '30760_NA', '5983_NA', '30690_NA', '30070_NA'), filter=T,
                            height = 3, width = 5, save=T, name='figure5c_new') # TODO

  figure <- ggpubr::ggarrange(figure5a, ggpubr::ggarrange(figure5b, figure5c,
                                                          labels = c('b', 'c'),
                                                          font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                                                          ncol=2, align = "h", heights = c(1,1)),
                              labels = c('a', ''),
                              font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                              ncol=1, align = "v", heights = c(4, 2.5))

  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}
