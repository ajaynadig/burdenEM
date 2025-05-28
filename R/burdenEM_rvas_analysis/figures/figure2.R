source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas_analysis/utils.R')

load_genebass_bhr_burdenEM_data <- function(genebass_burdenEM_path=Genebass_BurdenEM_PATH, filter=T){
  genebass_bhr <- load_genebass_bhr_results(bhr_threshold = 0.05)
  genebass_h2 <- read_csv(genebass_burdenEM_path) %>%
    merge(., genebass_pheno_info , by = 'phenotype', all.x=T) %>%
    merge(., genebass_neff, by = colnames(genebass_neff)[1:5])
  genebass_bhr_burdenEM <- genebass_bhr %>%
    filter(AF_bin %in% c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03')) %>%
    group_by(trait_type, phenocode, pheno_sex, coding, modifier, phenotype_key, ancestry, annotation, N, n_eff, description) %>%
    dplyr::summarize(bhr_h2 = sum(bhr_h2),
                     bhr_h2_se = sqrt(sum(bhr_h2_se^2)),
                     bhr_h2_z = bhr_h2/bhr_h2_se) %>%
    merge(.,genebass_h2 %>% select(trait_type, phenocode, pheno_sex, coding, modifier, phenotype_key=phenotype,  annotation, N, n_eff, description, total_h2, annotation, total_h2_CI_lower, total_h2_CI_upper),
          by = c('trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier', 'phenotype_key', 'N', 'n_eff', 'description','annotation'))%>%
    mutate(annotation = factor(annotation, levels = annotation_types2)) %>%
    filter(trait_type == 'continuous')
  if(filter) genebass_bhr_burdenEM <- genebass_bhr_burdenEM %>% filter(annotation %in% c('pLoF'))
  return(genebass_bhr_burdenEM)
}

load_genebass_aou_burdenEM_data <- function(aou_burdenEM_path=AoU_BurdenEM_PATH,
                                            genebass_burdenEM_path=Genebass_BurdenEM_PATH,
                                            filter=T){
  aou_h2 <- read_csv(aou_burdenEM_path) %>%
    mutate(ancestry = 'eur', phenoname = phenotype) %>%
    merge(., aou_neff, by = c('ancestry', 'phenoname') ) %>%
    merge(., full_pheno_sum %>% select(trait_type, ancestry = pop, phenoname, n_cases, n_controls, description), by = c('ancestry', 'phenoname')) %>%
    filter(trait_type == 'continuous')

  genebass_h2 <- read_csv(genebass_burdenEM_path) %>%
    merge(., genebass_pheno_info , by = 'phenotype', all.x=T) %>%
    merge(., genebass_neff, by = colnames(genebass_neff)[1:5])

  aou_genebass_merged <- aou_h2 %>%
    mutate(ukb_phenotype = ukb_phenotypes[phenotype]) %>%
    merge(., genebass_h2 %>%
            mutate(aou_phenotype = aou_phenotypes[phenotype]),
          by.x = c('phenotype', 'annotation'), by.y = c('aou_phenotype', 'annotation'), suffix = c('.aou', '.genebass')) %>%
    mutate(annotation = factor(annotation, levels = annotations))
  if(filter) aou_genebass_merged <- aou_genebass_merged %>% filter(annotation %in% c('pLoF'))
  return(aou_genebass_merged)
}

load_aou_ancestry_bhr_data <- function(bhr_threshold=0.05, filter=T){
  # main_bhr <- read_csv(get_aou_bhr_path('meta')) # TODO: waiting for aou meta-analysis
  main_bhr <- load_genebass_bhr_results(bhr_threshold = 0.05)
  aou_bhr <- compute_aou_bhr_h2_results(read_csv(get_aou_bhr_path('afr')),
                                        AF_bins = c('singleton', '2_1e-04', '1e-04_1e-03')) %>%
    rbind(., compute_aou_bhr_h2_results(read_csv(get_aou_bhr_path('amr')),
                                        AF_bins = c('singleton', '2_1e-04', '1e-04_1e-03'))) %>%
    rbind(., compute_aou_bhr_h2_results(read_csv(get_aou_bhr_path('eur')))) %>%
    mutate(ukb_phenotype = ukb_phenotypes[phenotype_key])  %>%
    merge(., full_pheno_sum %>% select(trait_type, ancestry = pop, phenoname, n_cases, n_controls, description), by.x = c('ancestry', 'phenotype_key') , by.y = c('ancestry', 'phenoname') )%>%
    merge(., aou_neff, by.x = c('ancestry', 'phenotype_key') , by.y = c('ancestry', 'phenoname') ) %>%
    mutate(bhr_z = bhr_h2/bhr_h2_se,
           bhr_p = 2 * (1 - pnorm(abs(bhr_z))),
           bhr_significant = bhr_p < bhr_threshold,
           ancestry = pops[ancestry])
  aou_main_bhr <- aou_bhr %>%
    merge(., compute_aou_bhr_h2_results(main_bhr), by.x = c('annotation', 'ukb_phenotype'), by.y = c('annotation', 'phenotype_key'), suffix = c('.aou', '.main')) %>%
    mutate(annotation = factor(annotation, levels = annotations))
  if(filter) aou_main_bhr <- aou_main_bhr %>% filter(annotation %in% c('pLoF'))
  return(aou_main_bhr)
}

load_aou_genebass_enrichment_data <- function(){
  genebass_h2 <- load_genebass_BurdenEM_results()
  aou_h2 <- load_aou_burdenEM_ancestry_data()
  enrichment <- rbind(
    aou_h2 %>% mutate(label = '1') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_1, enrichment_lower = enrichment_CI_lower_1, enrichment_upper = enrichment_CI_upper_1),
    aou_h2 %>% mutate(label = '2') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_2, enrichment_lower = enrichment_CI_lower_2, enrichment_upper = enrichment_CI_upper_2),
    aou_h2 %>% mutate(label = '3') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_3, enrichment_lower = enrichment_CI_lower_3, enrichment_upper = enrichment_CI_upper_3),
    aou_h2 %>% mutate(label = '4') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_4, enrichment_lower = enrichment_CI_lower_4, enrichment_upper = enrichment_CI_upper_4),
    aou_h2 %>% mutate(label = '5') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_5, enrichment_lower = enrichment_CI_lower_5, enrichment_upper = enrichment_CI_upper_5),
    genebass_h2 %>% mutate(label = '1', ancestry = 'genebass') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_1, enrichment_lower = enrichment_CI_lower_1, enrichment_upper = enrichment_CI_upper_5),
    genebass_h2 %>% mutate(label = '2', ancestry = 'genebass') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_2, enrichment_lower = enrichment_CI_lower_2, enrichment_upper = enrichment_CI_upper_5),
    genebass_h2 %>% mutate(label = '3', ancestry = 'genebass') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_3, enrichment_lower = enrichment_CI_lower_3, enrichment_upper = enrichment_CI_upper_5),
    genebass_h2 %>% mutate(label = '4', ancestry = 'genebass') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_4, enrichment_lower = enrichment_CI_lower_4, enrichment_upper = enrichment_CI_upper_5),
    genebass_h2 %>% mutate(label = '5', ancestry = 'genebass') %>% select(annotation, phenotype, ancestry, label, enrichment = enrichment_5, enrichment_lower = enrichment_CI_lower_5, enrichment_upper = enrichment_CI_upper_5)
  )
  return(enrichment)
}


load_aou_genebass_phenotype_lst <- function(matched_pheno){
  genebass_h2 <- load_genebass_BurdenEM_results()
  matched_pheno <- matched_pheno %>%
    filter(ukb_phenocode %in% genebass_h2$phenocode & aou_phenocode %in% aou_h2$phenoname)

  sub_h2 <- genebass_h2 %>%
    mutate(ancestry = 'genebass') %>%
    filter(phenocode %in% unique(matched_pheno$ukb_phenocode)) %>%
    select(ancestry, annotation, total_h2, n_eff, phenotype = phenotype) %>%
    rbind(., aou_h2 %>%
            filter(phenoname %in% unique(matched_pheno$aou_phenocode)) %>%
            select(ancestry, annotation, total_h2, n_eff, phenotype=phenoname)) %>%
    mutate(annotation = factor(annotation, levels = annotation_types2),
           ancestry = factor(ancestry, levels = c('genebass', 'eur', 'afr', 'amr')))
  return(unique(sub_h2$phenotype))
}


plot_figure2a <- function(genebass_burdenEM_path=Genebass_BurdenEM_PATH,
                          color_field = 'annotation',
                          emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                              '30760_NA', '5983_NA', '30690_NA', '30070_NA'),
                          filter = T, height=3, width=5, save=F, name='figure2a'){
  genebass_bhr_burdenEM <- load_genebass_bhr_burdenEM_data(
    genebass_burdenEM_path=genebass_burdenEM_path, filter=filter
  )

  phenotypes <- unique(genebass_bhr_burdenEM$phenotype_key)
  colors_emphasis <- rep('gray70', length(phenotypes))
  names(colors_emphasis) <- phenotypes
  colors_emphasis[emphasis_traits] <- c("#C95693", "#7B03D1", "#D9B407", "#7B8FEA",
                                        "#B20C00", "#F36E04", "#67D720","#227E18","#43BBA7")
  if(color_field != 'annotation'){
    genebass_bhr_burdenEM[, color_field] <- genebass_bhr_burdenEM[, 'phenotype_key']
    genebass_bhr_burdenEM[, color_field][!(genebass_bhr_burdenEM$phenotype_key %in% emphasis_traits)] <- ""
  }
  genebass_bhr_burdenEM[, 'description'][!(genebass_bhr_burdenEM$phenotype_key %in% emphasis_traits)] <- ""


  figure <- genebass_bhr_burdenEM   %>%
    ggplot + aes(x = bhr_h2, y = total_h2, color = get(color_field), label=description) +
    geom_point() +
    geom_errorbar(aes(ymin = total_h2_CI_lower,ymax = total_h2_CI_upper), lwd=0.3, alpha=0.3) +
    geom_errorbarh(aes(xmin = bhr_h2- 1.96*as.numeric(bhr_h2_se),xmax = bhr_h2 + 1.96*as.numeric(bhr_h2_se)),
                   lwd=0.3, alpha=0.3) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    scale_y_continuous(label=percent) +
    scale_x_continuous(label=percent) +
    labs(x = expression(BHR~h[Burden]^2), y = expression(BurdenEM~h[Burden]^2)) +
    facet_wrap(~annotation, scale = 'free', labeller = label_type, ncol=2)  +
    guides(color = 'none') + themes
  if(color_field == 'annotation'){
    figure <- figure + annotation_color_scale2
  }else{
    figure <- figure + scale_color_manual(values = colors_emphasis) +
      geom_text_repel(max.overlaps = 200, force = 100, size = 2,point.padding = 1,seed = 10332134) +
      theme(strip.background = element_blank(),
            strip.text = element_blank()) +
      labs(title = 'Genebass Burden Heritability: BHR vs. BurdenEM')
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}


plot_figure2b <- function(aou_burdenEM_path=AoU_BurdenEM_PATH,
                          genebass_burdenEM_path=Genebass_BurdenEM_PATH,
                          color_field = 'annotation',
                          emphasis_traits = c('height', '3028288', '3024929', '3004501',
                                              '3013721', '3007070', '3022192', '3027114', 'blood-pressure-diastolic-mean'),
                          filter = T, height=3, width=5, save=F, name='figure2b'){

  aou_genebass_merged <- load_genebass_aou_burdenEM_data(
    aou_burdenEM_path=AoU_BurdenEM_PATH,
    genebass_burdenEM_path=Genebass_BurdenEM_PATH,
    filter=filter)

  phenotypes <- unique(aou_genebass_merged$phenotype)
  colors_emphasis <- rep('gray70', length(phenotypes))
  names(colors_emphasis) <- phenotypes
  colors_emphasis[emphasis_traits] <- c("#C95693", "#7B03D1", "#D9B407", "#7B8FEA",
                                        "#B20C00", "#F36E04", "#67D720","#227E18","#43BBA7")
  if(color_field != 'annotation'){
    aou_genebass_merged[, color_field] <- aou_genebass_merged[, 'phenotype']
    aou_genebass_merged[, color_field][!(aou_genebass_merged$phenotype %in% emphasis_traits)] <- ""
  }
  aou_genebass_merged[, 'description.genebass'][!(aou_genebass_merged$phenotype %in% emphasis_traits)] <- ""

  figure <- aou_genebass_merged %>%
    ggplot + aes(x = total_h2.genebass, y = total_h2.aou, color = get(color_field), label = description.genebass) +
    labs(x = expression(Genebass~h[Burden]^2), y = expression(italic(All~by~All)~h[Burden]^2)) +
    geom_errorbar(aes(ymin = total_h2_CI_lower.aou,ymax = total_h2_CI_upper.aou), lwd=0.3, alpha=0.3) +
    geom_errorbarh(aes(xmin =  total_h2_CI_lower.genebass,xmax = total_h2_CI_upper.genebass), lwd=0.3, alpha=0.3) +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    scale_y_continuous(label=percent) +
    scale_x_continuous(label=percent) +
    geom_point() +themes +
    facet_wrap(~annotation, scale = 'free', labeller = label_type, ncol=2)+
    guides(color = 'none')

  if(color_field == 'annotation'){
    figure <- figure + annotation_color_scale2
  }else{
    figure <- figure + scale_color_manual(values = colors_emphasis) +
      geom_text_repel(max.overlaps = 200, force = 100, size = 2,point.padding = 1,seed = 10332134) +
      theme(strip.background = element_blank(),
            strip.text = element_blank()) +
      labs(title = 'BurdenEM Burden Heritability: Genebass vs. AoU')
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}


plot_figure2c <- function(color_field = 'annotation',
                          emphasis_traits = c('height', '3028288', '3024929', '3004501',
                                              '3013721', '3007070', '3022192', '3027114', 'blood-pressure-diastolic-mean'),
                          filter = T, height=3, width=5, save=F, name='figure2c'){

  aou_ancestry_bhr <- load_aou_ancestry_bhr_data(filter=filter)

  phenotypes <- unique(aou_ancestry_bhr$phenotype_key)
  colors_emphasis <- rep('gray70', length(phenotypes))
  names(colors_emphasis) <- phenotypes
  colors_emphasis[emphasis_traits] <- c("#C95693", "#7B03D1", "#D9B407", "#7B8FEA",
                                        "#B20C00", "#F36E04", "#67D720","#227E18","#43BBA7")
  if(color_field != 'annotation'){
    aou_ancestry_bhr[, color_field] <- aou_ancestry_bhr[, 'phenotype_key']
    aou_ancestry_bhr[, color_field][!(aou_ancestry_bhr$phenotype_key %in% emphasis_traits)] <- ""
  }
  aou_ancestry_bhr[, 'description'][!(aou_ancestry_bhr$phenotype_key %in% emphasis_traits)] <- ""

  figure <- aou_ancestry_bhr %>%
    filter(bhr_h2.aou < 0.2) %>%
    ggplot + aes(x = bhr_h2.main, y = bhr_h2.aou, color = get(color_field), label = description) +
    labs(x = expression(BHR~h[Burden]^2~(Genebass)), y = expression(BHR~h[Burden]^2)) +
    geom_errorbar(aes(ymin = bhr_h2.aou- 1.96*as.numeric(bhr_h2_se.aou),
                      ymax = bhr_h2.aou + 1.96*as.numeric(bhr_h2_se.aou)), lwd=0.3, alpha=0.3) +
    geom_errorbarh(aes(xmin = bhr_h2.main- 1.96*as.numeric(bhr_h2_se.main),
                       xmax = bhr_h2.main + 1.96*as.numeric(bhr_h2_se.main)), lwd=0.3, alpha=0.3) +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    scale_y_continuous(label=percent) +
    scale_x_continuous(label=percent) +
    geom_point() +themes +
    facet_grid(annotation~ancestry.aou, scale = 'free', labeller = label_type)+
    guides(color = 'none')

  if(color_field == 'annotation'){
    figure <- figure + annotation_color_scale2 +
      scale_y_continuous(label=percent, limits = c(-0.1, 0.1)) +
      scale_x_continuous(label=percent, limits = c(-0.01, 0.03))
  }else{
    figure <- figure + scale_color_manual(values = colors_emphasis) +
      geom_text_repel(max.overlaps = 200, force = 100, size = 2,point.padding = 1,seed = 10332134)  +
      theme(strip.background.y = element_blank(),
            strip.text.y = element_blank()) +
      labs(title = 'BHR Burden Heritability across AoU ancestry groups')
  }
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure2d <- function(filter = T, height=3, width=10, save=F, name='figure2d'){
  enrichment <- load_aou_genebass_enrichment_data()
  phenotypes_to_include <- load_aou_genebass_phenotype_lst(matched_pheno = matched_pheno)
  if(filter) enrichment <- enrichment %>% filter(annotation %in% ANNOTATIONS[1:2])
  figure <- enrichment  %>%
    filter(phenotype %in% phenotypes_to_include) %>%
    mutate(annotation = factor(annotation, levels = ANNOTATIONS))%>%
    ggplot + aes(x = label, y = enrichment, color= ancestry, fill = ancestry) +
    geom_boxplot(alpha = 0.2, width = .45, position='dodge') +
    scale_x_discrete(labels = annotation_names2) +
    labs(x ='Constraint Quintiles', y = expression(h[Burden]^2~Enrichment)) + themes +
    scale_color_manual(name = 'Data (Ancestry)',values = pop_colors, breaks = c('genebass', 'eur', 'afr', 'amr'), labels = c('UKB (EUR)', paste0('AoU (', toupper(c('eur', 'afr', 'amr')), ')')))+
    scale_fill_manual(name = 'Data (Ancestry)',values = pop_colors, breaks = c('genebass', 'eur', 'afr', 'amr'), labels = c('UKB (EUR)', paste0('AoU (', toupper(c('eur', 'afr', 'amr')), ')'))) +
    facet_grid(~annotation, labeller = label_type, scale = 'free')
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}

plot_figure2 <-  function(height = 6, width = 10, save=F, name='figure2'){
  figure2a <- plot_figure2a(genebass_burdenEM_path=Genebass_BurdenEM_PATH,
                            color_field = 'emphasis',
                            emphasis_traits = c('50_NA', '30780_NA', '30100_NA', '30040_NA', '100022_NA',
                                                '30760_NA', '5983_NA', '30690_NA', '30070_NA'),
                            filter = T, height=3, width=5, save=F)
  figure2b <- plot_figure2b(aou_burdenEM_path=AoU_BurdenEM_PATH,
                            genebass_burdenEM_path=Genebass_BurdenEM_PATH,
                            color_field = 'emphasis',
                            emphasis_traits = c('height', '3028288', '3024929', '3004501',
                                                '3013721', '3007070', '3022192', '3027114', 'blood-pressure-diastolic-mean'),
                            filter = T, height=3, width=5, save=F, name='figure2b')
  figure2c <- plot_figure2c(color_field = 'emphasis',
                            emphasis_traits = c('height', '3028288', '3024929', '3004501',
                                                '3013721', '3007070', '3022192', '3027114', 'blood-pressure-diastolic-mean'),
                            filter = T, height=3, width=10, save=F, name='figure2c')
  figure2d <- plot_figure2d(filter = T, height=3, width=10, save=T, name='figure2d')

  figure <- ggpubr::ggarrange(ggpubr::ggarrange(figure2a, figure2b,
                                                labels = c('a', 'b'),
                                                font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                                                ncol=2, align = "h", heights = c(1,1)), figure2c, figure2d,
                              labels = c('', 'c', 'd'),
                              font.label = list(size = 15, color = "black", face = "bold", family = 'sans'),
                              ncol=1, align = "v", heights = c(1, 1.2, 1.2))
  if(save) output_figure(figure, 'Manuscript', name, height = height, width = width)
  return(figure)
}
