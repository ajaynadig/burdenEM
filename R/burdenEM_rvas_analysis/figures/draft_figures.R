source('~/Dropbox (Partners HealthCare)/github_repo/burdenEM/R/burdenEM_rvas_analysis/utils.R')
#### BHR results: paper vs. new
genebass_bhr <- load_genebass_bhr_results(bhr_threshold = 0.05)
bhr_paper <- load_genebass_bhr_paper_results(bhr_threshold=0.05)
genebass_bhr_full <- genebass_bhr %>%
  select(phenotype_key, AF_bin, annotation, bhr_h2, bhr_h2_se, bhr_z, bhr_p, bhr_significant, description) %>%
  merge(., bhr_paper %>% select(phenotype_key, AF_bin, annotation, bhr_h2, bhr_h2_se, bhr_z, bhr_p, bhr_significant),
        by = c('phenotype_key', 'AF_bin', 'annotation'), suffix = c('_new', '_paper'))%>%
  mutate(annotation = factor(annotation, levels = annotation_types2)) %>%
  mutate(significance = case_when(
    bhr_significant_new & bhr_significant_paper ~ 'Both',
    bhr_significant_new & !bhr_significant_paper ~ 'New',
    !bhr_significant_new & bhr_significant_paper ~ 'Paper',
    !bhr_significant_new & !bhr_significant_paper ~ 'None',
  )) %>%
  mutate(AF_bin = factor(AF_bin, levels = c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03'), labels = c('AF: [0, 0.00001)', 'AF: [0.00001, 0.0001)', 'AF: [0.0001, 0.001)')))

p <- genebass_bhr_full  %>%
  ggplot + aes(x = bhr_h2_paper, y = bhr_h2_new, color = annotation) +
  geom_point(aes(pch = significance), size = 4) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_shape_manual(name=paste0('Nominal significance (0.05)'), breaks = c('Both', 'New', 'Paper', 'None'), values=c("\u25CF", "\u25D0","\u25D1", "\u25CB")) +
  geom_errorbar(aes(ymin = bhr_h2_new- 1.96*as.numeric(bhr_h2_se_new),ymax = bhr_h2_new + 1.96*as.numeric(bhr_h2_se_new)), lwd=0.1, alpha=0.1) +
  geom_errorbarh(aes(xmin = bhr_h2_paper- 1.96*as.numeric(bhr_h2_se_paper),xmax = bhr_h2_paper + 1.96*as.numeric(bhr_h2_se_paper)), lwd=0.1, alpha=0.1) +
  annotation_color_scale2 + themes + labs(x = expression(bold(BHR~h^2~(paper))), y = expression(bold(BHR~h^2~(new))))+ facet_grid(~AF_bin, scale = 'free') +
  geom_text_repel(data = genebass_bhr_full %>% filter(bhr_h2_paper/bhr_h2_new > 1.2 | bhr_h2_paper/bhr_h2_new < 0.9), aes(label=description), size = 2, max.overlaps = 20)
p
output_figure(p, 'genebass_BHR_h2_new_old_comparison', height = 4, width = 10)

p <- genebass_bhr_full  %>%
  filter(significance == 'Both') %>%
  ggplot + aes(x = bhr_h2_paper, y = bhr_h2_new, color = annotation) +
  geom_point(aes(pch = significance), size = 4) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_shape_manual(name=paste0('Nominal significance (0.05)'), breaks = c('Both', 'New', 'Paper', 'None'), values=c("\u25CF", "\u25D0","\u25D1", "\u25CB")) +
  geom_errorbar(aes(ymin = bhr_h2_new- 1.96*as.numeric(bhr_h2_se_new),ymax = bhr_h2_new + 1.96*as.numeric(bhr_h2_se_new)), lwd=0.1, alpha=0.1) +
  geom_errorbarh(aes(xmin = bhr_h2_paper- 1.96*as.numeric(bhr_h2_se_paper),xmax = bhr_h2_paper + 1.96*as.numeric(bhr_h2_se_paper)), lwd=0.1, alpha=0.1) +
  annotation_color_scale2 + themes + labs(x = expression(bold(BHR~h^2~(paper))), y = expression(bold(BHR~h^2~(new))))+ facet_grid(~AF_bin, scale = 'free') +
  geom_text_repel(data = genebass_bhr_full %>% filter(bhr_h2_paper/bhr_h2_new > 1.2 | bhr_h2_paper/bhr_h2_new < 0.9), aes(label=description), size = 2, max.overlaps = 20) +
  guides(shape = 'none')
p
output_figure(p, 'BHR_h2_new_old_comparison_filtered', height = 4, width = 10)

##### Enrichment analysis
genebass_h2 <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BurdenEM/genebass_eur_iter_1000_full_h2_result_tmp_burdenEM.csv') %>%
  merge(., genebass_pheno_info , by = 'phenotype', all.x=T) %>%
  merge(., genebass_neff, by = colnames(genebass_neff)[1:5])
aou_afr_h2 <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BurdenEM/aou_afr_iter_1000_full_h2_result_tmp_burdenEM.csv') %>%
  mutate(ancestry = 'afr', phenoname = phenotype)
aou_amr_h2 <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BurdenEM/aou_amr_iter_1000_full_h2_result_tmp_burdenEM.csv') %>%
  mutate(ancestry = 'amr', phenoname = phenotype)
aou_eur_h2 <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BurdenEM/aou_eur_iter_1000_full_h2_result_tmp_burdenEM.csv') %>%
  mutate(ancestry = 'eur', phenoname = phenotype)

aou_h2 <- rbind(aou_afr_h2, aou_amr_h2, aou_eur_h2) %>%
  merge(., aou_neff, by = c('ancestry', 'phenoname') ) %>%
  merge(., full_pheno_sum %>% select(trait_type, ancestry = pop, phenoname, n_cases, n_controls, description), by = c('ancestry', 'phenoname'))

aou_genebass_merged <- aou_h2 %>%
  mutate(ukb_phenotype = ukb_phenotypes[phenotype]) %>%
  merge(., genebass_h2 %>%
          mutate(aou_phenotype = aou_phenotypes[phenotype]),
        by.x = c('phenotype', 'annotation'), by.y = c('aou_phenotype', 'annotation'), suffix = c('.aou', '.genebass'))

enrichment <- rbind(
  aou_h2 %>% mutate(label = '1') %>% select(annotation, phenotype, ancestry, label, h2 = enrichment_1, h2_lower = enrichment_CI_lower_1, h2_upper = enrichment_CI_upper_1),
  aou_h2 %>% mutate(label = '2') %>% select(annotation, phenotype, ancestry, label, h2 = enrichment_2, h2_lower = enrichment_CI_lower_2, h2_upper = enrichment_CI_upper_2),
  aou_h2 %>% mutate(label = '3') %>% select(annotation, phenotype, ancestry, label, h2 = enrichment_3, h2_lower = enrichment_CI_lower_3, h2_upper = enrichment_CI_upper_3),
  aou_h2 %>% mutate(label = '4') %>% select(annotation, phenotype, ancestry, label, h2 = enrichment_4, h2_lower = enrichment_CI_lower_4, h2_upper = enrichment_CI_upper_4),
  aou_h2 %>% mutate(label = '5') %>% select(annotation, phenotype, ancestry, label, h2 = enrichment_5, h2_lower = enrichment_CI_lower_5, h2_upper = enrichment_CI_upper_5)
)

p <- enrichment  %>%
  mutate(annotation = factor(annotation, levels = c('pLoF', 'missense_damaging', 'missense_benign', 'synonymous')))%>%
  ggplot + aes(x = label, y = h2, color= annotation) +
  # geom_pointrange(aes(ymin = h2_lower, ymax = h2_upper), stat = "identity", position = position_dodge(width = 0.6)) +
  geom_boxplot() +
  annotation_color_scale2 +
  scale_x_discrete(labels = annotation_names2) +
  labs(x ='Constraint Quintiles', y = expression(bold(Burden~h^2~Enrichment))) + themes +
  guides(color = guide_legend(title = "Variant Type"))+
  facet_wrap(~ancestry, labeller = label_type, scale = 'free', nrow = 3)+
  theme(strip.text = element_text(face="bold"))
p
output_figure(p, 'aou_BurdenEM_h2_enrichment', height = 6, width = 12)

enrichment <- rbind(
  genebass_h2 %>% mutate(label = '1') %>% select(annotation, phenotype, label, h2 = enrichment_1, h2_lower = enrichment_CI_lower_1, h2_upper = enrichment_CI_upper_1),
  genebass_h2 %>% mutate(label = '2') %>% select(annotation, phenotype, label, h2 = enrichment_2, h2_lower = enrichment_CI_lower_2, h2_upper = enrichment_CI_upper_2),
  genebass_h2 %>% mutate(label = '3') %>% select(annotation, phenotype, label, h2 = enrichment_3, h2_lower = enrichment_CI_lower_3, h2_upper = enrichment_CI_upper_3),
  genebass_h2 %>% mutate(label = '4') %>% select(annotation, phenotype, label, h2 = enrichment_4, h2_lower = enrichment_CI_lower_4, h2_upper = enrichment_CI_upper_4),
  genebass_h2 %>% mutate(label = '5') %>% select(annotation, phenotype, label, h2 = enrichment_5, h2_lower = enrichment_CI_lower_5, h2_upper = enrichment_CI_upper_5)
)

p <- enrichment  %>%
  mutate(annotation = factor(annotation, levels = c('pLoF', 'missense_damaging', 'missense_benign', 'synonymous')))%>%
  ggplot + aes(x = label, y = h2, color= annotation) +
  # geom_pointrange(aes(ymin = h2_lower, ymax = h2_upper), stat = "identity", position = position_dodge(width = 0.6)) +
  geom_boxplot() +
  annotation_color_scale2 +
  scale_x_discrete(labels = annotation_names2) +
  labs(x ='Constraint Quintiles', y = expression(bold(h[Burden]^2~Enrichment))) + themes +
  guides(color = guide_legend(title = "Variant Type"))+
  theme(strip.text = element_text(face="bold"))
p
output_figure(p, 'genebass_BurdenEM_h2_enrichment', height = 2.5, width = 8)


##### positive negative h2 comparison pLoF and missense
genebass_neff <- read_delim('~/Desktop/burdenEM_results/data/genebass/genebass_n_eff.txt.bgz', delim = '\t')
genebass_h2 <- read_csv('~/Desktop/burdenEM_results/BurdenEM/genebass_eur_iter_1000_full_h2_result_tmp_burdenEM.csv') %>%
  merge(., genebass_pheno_info , by = 'phenotype', all.x=T) %>%
  merge(., genebass_neff, by = colnames(genebass_neff)[1:5])
genebass_h2 <- genebass_h2 %>%
  mutate(positive_h2 = total_h2*prop_positive_h2,
         negative_h2 = total_h2*prop_negative_h2) %>%
  select(annotation, description, prop_negative_h2, prop_positive_h2) %>%
  pivot_wider(names_from = 'annotation', values_from = c('prop_negative_h2', 'prop_positive_h2'))

p <- genebass_h2 %>%
  # mutate(description = if_else(grepl(description, 'Date'), str_to_sentence(str_replace(str_split(description, 'first reported \\(') %>% map_chr(., 2), ')', '')), description))
  ggplot + aes(y = prop_positive_h2_pLoF, x = prop_positive_h2_missense_damaging) +
  geom_point(stat='identity', color = '#8b733d', size = 2) + themes +
  geom_hline(yintercept = 0.5, lty=2, color = 'gray20') +
  geom_vline(xintercept = 0.5, lty=2, color = 'gray20') +
  geom_abline(intercept = 0, slope = 1, lty=2, color = 'gray20') +
  scale_y_continuous(label= percent) +
  scale_x_continuous(label= percent) +
  # scale_x_discrete(breaks = rev(unique(description))) +
  # scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels = percent(c(0, 0.5, 1, 0.5, 1, 0.5, 1))) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(hjust = 1, size = 13),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20)) +
  labs(y = expression(bold(h[Burden]^2~from~positive~effect~size~(pLoF))), x = expression(bold(h[Burden]^2~from~positive~effect~size~(Missense)))) +
  geom_label_repel(data = genebass_h2 %>% filter((prop_positive_h2_pLoF > 0.5 & prop_positive_h2_missense_damaging < 0.5) |
                                                   (prop_positive_h2_pLoF < 0.5 & prop_positive_h2_missense_damaging > 0.5) ),
                   aes(label = description), color = '#8b733d', max.overlaps = 5, size = 2) +
  theme(plot.margin = unit(c(0.5, 0, 0,0), "cm"))
p
output_figure(p, 'pLoF_missense_proportion_positive_comparison', height = 4, width = 6)

#### Long gene short gene analysis
long_gene_burdenscore <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/genebass_eur_long_genes_tmp_full_burden_score.csv') %>%
  dplyr::group_by(gene_name, annotation, pheno) %>%
  dplyr::summarize(burdenScore = sum(burden_score, na.rm=T)) %>%
  mutate(label = 'Long Genes (top 20%)')
short_gene_burdenscore <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/genebass_eur_short_genes_tmp_full_burden_score.csv') %>%
  dplyr::group_by(gene_name, annotation, pheno) %>%
  dplyr::summarize(burdenScore = sum(burden_score, na.rm=T)) %>%
  mutate(label = 'Short Genes (rest 80%)')
gene_burden_score <- rbind(long_gene_burdenscore, short_gene_burdenscore) %>%
  group_by(annotation, pheno, label) %>%
  dplyr::summarize(
    total_burden_score = sum(burdenScore)
  )

long_gene_bhr <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/genebass_eur_long_genes_tmp_full.csv')
short_gene_bhr <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BHR/genebass_eur_short_genes_tmp_full.csv')
long_short_bhr <- rbind(
  long_gene_bhr %>% mutate(label = 'Long Genes (top 20%)'),
  short_gene_bhr %>% mutate(label = 'Short Genes (rest 80%)')
)%>%
  filter(AF_bin %in% c('0e+00_1e-05', '1e-05_1e-04', '1e-04_1e-03')) %>%
  group_by(ancestry, phenotype_key, annotation, N, label) %>%
  dplyr::summarize(bhr_h2 = sum(bhr_h2),
                   bhr_h2_se = sqrt(sum(bhr_h2_se^2))) %>%
  mutate(model = 'BHR') %>%
  ungroup() %>%
  select(pheno=phenotype_key, annotation, h2 = bhr_h2, label, model) %>%
  group_by(pheno, annotation,  model) %>%
  mutate(h2_all = sum(h2),
         prop_h2 = h2/h2_all)

long_gene_burdenEM <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BurdenEM/genebass_eur_n_boot_10_grid_size_100_iter_1000_full_h2_result_tmp_burdenEM_union_cpt_long_gene_kept.csv')
short_gene_burdenEM <- read_csv('~/Dropbox (Partners HealthCare)/burdenEM/burdenEM_results/BurdenEM/genebass_eur_n_boot_10_grid_size_100_iter_1000_full_h2_result_tmp_burdenEM_union_cpt_long_gene_removed.csv')
long_short_burdenEM <- rbind(
  long_gene_burdenEM %>% mutate(label = 'Long Genes (top 20%)'),
  short_gene_burdenEM %>% mutate(label = 'Short Genes (rest 80%)')
) %>%
  mutate(model = 'BurdenEM') %>%
  select(pheno=phenotype, annotation, h2 = total_h2, label, model) %>%
  group_by(pheno, annotation, model) %>%
  mutate(h2_all = sum(h2),
         prop_h2 = h2/h2_all) %>%
  rbind(., long_short_bhr) %>%
  merge(., gene_burden_score, by = c('annotation', 'label', 'pheno')) %>%
  mutate(h2_per_burdenscore = h2/total_burden_score)

figure <- long_short_burdenEM %>%
  filter(annotation == 'pLoF') %>%
  mutate(phenotype = phenos_to_select[pheno]) %>%
  filter(!is.na(phenotype)) %>%
  ggplot + aes(x = phenotype, y = h2, color = phenotype, group=interaction(phenotype, annotation),
               fill=phenotype, alpha = label) +
  geom_col(stat='jitter') + themes +
  scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
  scale_y_continuous(label= percent) +
  scale_color_manual(values = ntpr_colors) +
  scale_fill_manual(values = ntpr_colors) +
  labs(x =NULL, y=expression(h^2[Burden]), color = NULL, fill=NULL) +
  facet_grid(~model, labeller = label_type) +
  theme(legend.position = 'top', legend.direction = 'vertical') +
  guides(
    fill = guide_legend(nrow = 2, keywidth = 1, keyheight = 0.3),
    color = guide_legend(nrow = 2, keywidth = 1, keyheight = 0.3),
    alpha = guide_legend(nrow = 2, keywidth = 1, keyheight = 0.3)
  )
output_figure(figure, 'BurdenEM_results', 'long_short_gene_analysis_pLoF_abs_h2', height = 3, width = 7.5)

figure <- long_short_burdenEM %>%
  filter(annotation == 'pLoF') %>%
  mutate(phenotype = phenos_to_select[pheno]) %>%
  filter(!is.na(phenotype)) %>%
  ggplot + aes(x = phenotype, y = h2_per_burdenscore, color = phenotype, group=interaction(phenotype, annotation),
               fill=phenotype, alpha = label) +
  geom_point(stat='identity') + themes +
  scale_alpha_discrete(name = NULL, range = c(1, 0.2)) +
  scale_y_continuous(label= comma) +
  scale_color_manual(values = ntpr_colors) +
  scale_fill_manual(values = ntpr_colors) +
  labs(x =NULL, y=NULL, color = NULL, fill=NULL) +
  facet_grid(~model, labeller = label_type) +
  theme(legend.position = 'top', legend.direction = 'vertical') +
  guides(
    fill = guide_legend(nrow = 2, keywidth = 1, keyheight = 0.3),
    color = guide_legend(nrow = 2, keywidth = 1, keyheight = 0.3),
    alpha = guide_legend(nrow = 2, keywidth = 1, keyheight = 0.3)
  )
output_figure(figure, 'BurdenEM_results', 'long_short_gene_analysis_pLoF_h2_per_burdenscore', height = 3, width = 7.5)

p <- gene_burden_score %>%
  mutate(annotation = factor(annotation)) %>%
  ggplot + aes(x = burdenScore, color = annotation) +
  labs(x = 'Burden Score', y = NULL) +
  geom_density() + themes +
  scale_x_log10(label=comma) +
  annotation_color_scale2 + annotation_fill_scale2 +
  facet_grid(~label) +
  geom_text_repel(data = gene_mean_burden_score, aes(label = scientific(mean)), x = Inf, y = 0.5, vjust = 0, hjust = 1)

output_figure(p, 'BurdenEM_results', 'long_short_gene_burdenscore', 4, 7.5)


### Simulation figures
p1 <- sim_test %>%
  ggplot + aes(x = true_h2_rounded, y = burden_h2, fill = name, color = name,  group = interaction(name, true_h2_rounded)) +
  geom_boxplot(
    alpha = 0.8
  ) + geom_abline(slope = 1, intercept=0, lty = 2, color = '#8b733d') +
  labs(x = NULL, y = expression(bold(Estimated~h[Burden]^2)), color = NULL, fill = NULL) +
  scale_y_continuous(label=percent) +
  scale_color_manual(values = distinct_colors) +
  scale_fill_manual(values = distinct_colors) +
  scale_x_continuous(label = c('0%', '0.5%', '1%', '1.5%', '2%'), breaks = c(0, 0.005, 0.01, 0.015, 0.02),expand = expansion(mult = c(0.1, 0.1)))+
  scale_alpha_discrete(name = NULL, range = c(1, 0.2), labels = c('Positive effect', 'Negative effect')) +
  themes +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.position = c(0.25, 0.75),
    strip.text = element_text(face="bold", size = 13))
p1


p2 <- sim_test  %>%
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
  geom_pointrange(aes(ymin = min_estimated,ymax = max_estimated), lwd=0.5) + geom_abline(slope = 1, intercept=0, lty = 2, color = '#8b733d') +
  labs(x = expression(bold(Polygenicity~": "~true~N[genes]~to~explain~"50%"*~h[Burden]^2)),
       y = expression(bold(Estimated~N[genes])), color = NULL) +
  scale_color_manual(values = distinct_colors) +
  scale_fill_manual(values = distinct_colors) +
  # scale_x_continuous(label = c('0%', '0.5%', '1%', '1.5%', '2%'), breaks = c(0, 0.005, 0.01, 0.015, 0.02),expand = expansion(mult = c(0.1, 0.1)))+
  scale_alpha_discrete(name = NULL, range = c(1, 0.2), labels = c('Positive effect', 'Negative effect')) +
  themes +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = 'top',
    strip.text = element_text(face="bold", size = 13))
p2


figure <- ggpubr::ggarrange(p1, p,
                            labels = NULL,
                            nrow=1, align = "v",
                            widths = c(1,2))
output_figure(figure, 'test_sim_figure.png', height = 3, width = 12)

