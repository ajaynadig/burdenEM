source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas_analysis/utils.R')

fig1_colors <- c('#887a7e', '#8f89bb', '#bea0b6')
names(fig1_colors) <- c('Realistic (N = 100K)', 'More Polygenic', 'More Positive Effect')

load_burdenEM_sim_results_fig1 <- function(simulation_path){
  true_value <- load_burdenEM_sim_true_values(TRUE_SIM_PATH)
  sim_test <- read_csv(simulation_path) %>%
    mutate(true_h2_rounded = round(true_h2, 3)) %>%
    merge(., true_value, by = c('name', 'rep', 'true_h2_rounded'))

  sim_test <- sim_test %>%
    filter(name %in% c('realistic', 'more_polygenic', 'more_positive_effects')) %>%
    mutate(name = factor(name, levels = c('realistic', 'more_polygenic', 'more_positive_effects'),
                         labels = c('Realistic (N = 100K)', 'More Polygenic', 'More Positive Effect')))
  return(sim_test)
}

load_burdenEM_aou_random_pheno_qq <- function(anc, suffix=''){
  real_qq_results <- read_csv(paste0(burdenEM_ROOT, 'burdenEM_results/BurdenEM/aou_', anc,'_n_boot_10_grid_size_100_iter_1000_real_qq_result_tmp_burdenEM_random_new_cpt', suffix,'.csv'))
  est_qq_results <- read_csv(paste0(burdenEM_ROOT, 'burdenEM_results/BurdenEM/aou_', anc,'_n_boot_10_grid_size_100_iter_1000_est_qq_result_tmp_burdenEM_random_new_cpt', suffix,'.csv'))
  real_qq_results <- real_qq_results %>%
    mutate(phenotype = random_phenonames[phenotype]) %>%
    mutate(annotation = factor(annotation, levels = ANNOTATIONS, labels = annotation_names2))
  est_qq_results <- est_qq_results %>%
    mutate(phenotype = random_phenonames[phenotype])%>%
    mutate(annotation = factor(annotation, levels = ANNOTATIONS, labels = annotation_names2))
  return(list(real_qq_results, est_qq_results))
}

plot_figure1a <- function(sim_test, height = 3, width = 12, save = F, name='figure1a'){
  figure <- sim_test %>%
    ggplot + aes(x = true_h2_rounded, y = burden_h2, color = name, group = interaction(name, true_h2_rounded)) +
    geom_pointrange(
      stat = "summary",
      fun.min = min,
      fun.max = max,
      fun = median
    ) + geom_abline(slope = 1, intercept=0, color = 'grey', lwd=0.1) +
    labs(x = NULL, y = expression(Estimated~h[Burden]^2), color = NULL) +
    scale_color_manual(values = fig1_colors) +
    scale_fill_manual(values = fig1_colors) +
    scale_y_continuous(label=percent) +
    scale_x_continuous(label = c('0%', '0.5%', '1%', '1.5%', '2%'), breaks = c(0, 0.005, 0.01, 0.015, 0.02))+
    scale_alpha_discrete(name = NULL, range = c(1, 0.2), labels = c('Positive effect', 'Negative effect')) +
    themes + facet_grid(~name) + theme(plot.margin = unit(c(0.5, 0, 0, 0.3), 'cm'),)
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure1b <- function(sim_test, height = 3.5, width = 12, save=F, name='figure1b'){
  figure <- sim_test %>%
    mutate(positive_h2 = burden_h2*prop_positive_h2,
           negative_h2 = burden_h2*prop_negative_h2) %>%
    pivot_longer(cols = c('prop_positive_h2', 'prop_negative_h2'), names_to = 'alpha') %>%
    mutate(alpha = factor(str_split(alpha, '_') %>% map_chr(., 2), levels = c('positive', 'negative'), labels = c('Positive effect', 'Negative effect')) ,
    )%>%
    group_by(true_h2_rounded, name, alpha) %>%
    dplyr::summarize(value = mean(value)) %>%
    mutate(truth = if_else(name == 'More Positive Effect', 0.7, NA))%>%
    ggplot + aes(x = true_h2_rounded, y = value, color = name, group=interaction(true_h2_rounded, name), alpha = alpha, fill = name) +
    geom_col(stat='identity') +
    geom_hline(aes(yintercept = truth), color = 'grey', lwd=0.1) +
    labs(x = expression(True~h[Burden]^2), y = expression(h[Burden]^2~"(" * "%" * ")"), color=NULL, fill=NULL) +
    scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
    scale_y_continuous(label= percent) +
    scale_x_continuous(label = c('0.5%', '1%', '1.5%', '2%'), breaks = c(0.005, 0.01, 0.015, 0.02))+
    scale_color_manual(values = fig1_colors) +
    scale_fill_manual(values = fig1_colors) +
    themes + theme(
      plot.margin = unit(c(0.5, 0, 0, 0.3), 'cm'),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.y = element_text(hjust = 1),) +
    facet_grid(~name) + guides(colour = 'none', fill= 'none')
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure1c <- function(sim_test, height = 3.5, width = 12, save=F, name='figure1c'){
  figure <- sim_test %>%
    filter(true_h2_rounded > 0 ) %>%
    group_by(name, true_h2_rounded) %>%
    dplyr::summarize(
      median_true = median(true_num_genes_fifty, na.rm = T),
      median_estimated = median(polygenicity_50, na.rm = T),
      min_true = min(true_num_genes_fifty, na.rm = T),
      min_estimated = min(polygenicity_50, na.rm = T),
      max_true = max(true_num_genes_fifty, na.rm = T),
      max_estimated = max(polygenicity_50, na.rm = T)
    ) %>%
    ggplot + aes(x = median_true, y = median_estimated, color = name, group = interaction(name, true_h2_rounded), alpha=true_h2_rounded)+
    geom_point(size = 2.5)+
    geom_pointrange(aes(ymin = min_estimated,ymax = max_estimated), lwd=0.5) +
    geom_abline(slope = 1, intercept=0, lty = 1, color = 'grey', lwd=0.1) +
    labs(x = expression(Polygenicity~": "~true~N[genes]~to~explain~"50%"*~h[Burden]^2),
         y = expression(Estimated~N[genes]), color = NULL) +
    scale_color_manual(values = fig1_colors) +
    scale_fill_manual(values = fig1_colors) +
    scale_x_continuous(label=comma) +
    scale_y_continuous(label=comma) +
    scale_alpha_continuous(name = expression(True~h[Burden]^2, labels=percent_labels), range = c(0.3, 1)) +
    themes + theme(
      plot.margin = unit(c(0.5, 0, 0, 0.3), 'cm'),
      strip.background = element_blank(),
      strip.text.x = element_blank(),) +
    facet_grid(~name, scales = 'free') + guides(colour = 'none', fill= 'none')
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure1d <-  function(anc, suffix = '',height = 4, width = 8, save=F, name='figure1d', filter=T){
  data <- load_burdenEM_aou_random_pheno_qq(anc=anc, suffix=suffix)
  if(filter){
    real_qq <- data[[1]] %>%
      filter(annotation %in% c('pLoF', 'Synonymous') & phenotype != 'Binary (p = 0.1%)' & phenotype != 'Binary (p = 1%)')
    est_qq <- data[[2]] %>%
      filter(annotation %in% c('pLoF', 'Synonymous') & phenotype != 'Binary (p = 0.1%)'& phenotype != 'Binary (p = 1%)')
  }else{
    real_qq <- data[[1]]
    est_qq <- data[[2]]
  }

  figure <- real_qq %>%
    ggplot +
    aes(x = expected, y = observed) +
    geom_abline(intercept = 0, slope = 1, color='gray40', lwd = 0.3) +
    geom_point(alpha = 0.4, size = 0.5, color = '#A1B8D7') +
    geom_line(data = est_qq, aes(x = expected, y = observed), color='#D2A010', lwd = 0.3) +
    themes  + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.8), 'cm'),) +
    labs(x = 'Expected Z score', y = 'Observed Z score') +
    facet_grid(annotation~phenotype, labeller = label_type, scale= 'fixed')
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure1 <-  function(height = 9, width = 9, save=F, name='figure1'){
  sim_test <- load_burdenEM_sim_results_fig1(SIM_RESULTS_PATH)
  figure1a <- plot_figure1a(sim_test=sim_test, save=F)
  figure1b <- plot_figure1b(sim_test=sim_test, save=F)
  figure1c <- plot_figure1c(sim_test=sim_test, save=F)
  figure1d <- plot_figure1d(anc='eur', suffix = '', save=F)
  figure1top <- ggpubr::ggarrange(figure1a, figure1b, figure1c,
                                  labels = c('a', 'b', 'c'),
                                  font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                                  ncol=1, align = "v", heights = c(1,1,1))
  if(save) output_figure(figure1top, 'Manuscript', 'figure1top', height = 5.5, width = 12)

  figure <- ggpubr::ggarrange(figure1top, figure1d,
                              labels = c('', 'd'),
                              font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                              ncol=1, align = "v", heights = c(4, 2.5))
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}
