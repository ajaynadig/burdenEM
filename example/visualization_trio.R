source("~/Mirror/oconnor_rotation/rare_dn_h2/scripts/set_up.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/scripts/set_up_kaplanis.R")


resource_dir <- "~/Mirror/oconnor_rotation/rare_dn_h2/github/resources/"
output_dir = "~/Mirror/oconnor_rotation/rare_dn_h2/github/outputs/"
path_to_supptables <- paste0(output_dir, "SuppTables/")
path_to_figs <-paste0(output_dir, "Figures/")


load(paste0(output_dir, "models_autism_Apr25.Rdata"))
load(paste0(output_dir,"models_ddd_Apr25.Rdata"))


current_date <- format(Sys.Date(), "%b%y")


#collate the heritability results into a summary dataframe

heritability_enrich_autism = heritability_enrichment_table(autism_data,burdenEM_models_autism)
heritability_enrich_dd = heritability_enrichment_table(kaplanis_data, burdenEM_models_DDD)

SupplementaryTable2 = make_supptable(heritability_enrich_autism,
                                     heritability_enrich_autism$heritability$name[!grepl("New",heritability_enrich_autism$heritability$name) &
                                                                                    !grepl("Fu2022",heritability_enrich_autism$heritability$name)&
                                                                                    !grepl("DDID",heritability_enrich_autism$heritability$name) &
                                                                                    !grepl("ASC",heritability_enrich_autism$heritability$name) &
                                                                                    !grepl("GeneDx",heritability_enrich_autism$heritability$name) &
                                                                                    !grepl("SPARK",heritability_enrich_autism$heritability$name)])

write.table(SupplementaryTable2,paste0(path_to_supptables, "SupplementaryTable2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

SupplementaryTable7 = make_supptable(heritability_enrich_dd,
                                     heritability_enrich_dd$heritability$name)
write.table(SupplementaryTable7,paste0(path_to_supptables, "SupplementaryTable7.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#PTV: #FF5F8A
#Mis2: darkorange
#Mis1: #F9B332
#Mis0: #EFD09F
#Syn: grey50

palette_variantclass = c("PTV" = "#FF5F8A",
                         "Mis2" = "darkorange",
                         "Mis1" = "#F9B332",
                         "Mis0" = "#EFD09F",
                         "Syn" = "grey50")

Proband_Vmu_plot <- ggplot(data = heritability_enrich_autism$heritability[heritability_enrich_autism$heritability$prev_factor ==  autism_data$prev_factors[1] &
                                                                            grepl("Combined",heritability_enrich_autism$heritability$name) &
                                                                            grepl("Proband",heritability_enrich_autism$heritability$name),],
                           mapping = aes(x= factor(variant_class, levels = c("Syn","Mis0","Mis1","Mis2","PTV")),
                                         y = h2,
                                         ymin = h2_lower,
                                         ymax = h2_upper,
                                         fill = factor(variant_class, levels = c("Syn","Mis0","Mis1","Mis2","PTV"))
                           ))+
  geom_hline(yintercept = 0)+
  geom_pointrange(position = position_dodge2(width = 0.25), shape = 21, color = "black", size = 1)+
  scale_fill_manual(values = palette_variantclass)+
  theme_bhr_legend_gridlines()+
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale")+
  theme(strip.text.y = element_text(size = 15, angle = 270),
        axis.text.x = element_text(color = rev(palette_variantclass), face = "bold"))+
  ylim(0,0.04)+
  guides(fill = "none")

Sibling_Vmu_plot <- ggplot(data = heritability_enrich_autism$heritability[heritability_enrich_autism$heritability$prev_factor ==  autism_data$prev_factors[1] &
                                                                            grepl("Combined",heritability_enrich_autism$heritability$name) &
                                                                            grepl("Sibling",heritability_enrich_autism$heritability$name),],
                           mapping = aes(x= factor(variant_class, levels = c("Syn","Mis0","Mis1","Mis2","PTV")),
                                         y = h2,
                                         ymin = h2_lower,
                                         ymax = h2_upper,
                                         fill = factor(variant_class, levels = c("Syn","Mis0","Mis1","Mis2","PTV"))
                           ))+
  geom_hline(yintercept = 0)+
  geom_pointrange(position = position_dodge2(width = 0.25), shape = 21, color = "black", size = 1)+
  scale_fill_manual(values = palette_variantclass)+
  theme_bhr_legend_gridlines()+
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale")+
  theme(strip.text.y = element_text(size = 15, angle = 270),
        axis.text.x = element_text(color = rev(palette_variantclass), face = "bold"))+
  ylim(0,0.0375)+
  guides(fill = "none")+
  ggtitle("Unaffected Siblings; N = 9567 ")

ggsave(paste0(path_to_figs,"SupplementaryFigure2.pdf"),
       Sibling_Vmu_plot,device = cairo_pdf,width =4.5, height = 3,dpi = "retina")

Prevalence_MutVar_plot <- ggplot(data = heritability_enrich_autism$heritability[grepl("Combined",heritability_enrich_autism$heritability$name) &
                                                                                  grepl("Proband",heritability_enrich_autism$heritability$name),],
                                 mapping = aes(x= factor(variant_class, levels = c("Syn","Mis0","Mis1","Mis2","PTV")),
                                               y = h2,
                                               ymin = h2_lower,
                                               ymax = h2_upper,
                                               fill = factor(paste0(round(prev_mod*100,2),"%"), levels = c("1%","2.24%","2.76%"))
                                 ))+
  geom_hline(yintercept = 0)+
  geom_pointrange(position = position_dodge2(width = 0.6), shape = 21, color = "black", size = 1)+
  # scale_fill_manual(values = palette_variantclass)+
  theme_bhr_legend_gridlines()+
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale", fill = "Prevalence")+
  theme(strip.text.y = element_text(size = 15, angle = 270),
        axis.text.x = element_text(color = rev(palette_variantclass), face = "bold"))

ggsave(paste0(path_to_figs,"SupplementaryFigure3.pdf"),
       Prevalence_MutVar_plot,device = cairo_pdf,width =6, height = 4,dpi = "retina")



#Polygenicity
MutVar_pergene = aggregate_mutvar_per_gene(autism_data,
                                           burdenEM_models_autism[[1]][[5]],get_genetic_data(5,autism_data)$genetic_data,
                                           burdenEM_models_autism[[1]][[11]],get_genetic_data(11,autism_data)$genetic_data,
                                           prevalence = autism_data$baseprev * autism_data$prev_factors[1])


#Make annotation
gnomad_information <- data.frame(fread(paste0(resource_dir,"gnomad.v2.1.1.lof_metrics.by_gene.txt"),
                                       select = c("gene","gene_id", "chromosome", "start_position",	"end_position", "pLI")))

top15_genes_id = names(MutVar_pergene$full)[1:15]
top15_genes_name = gnomad_information$gene[match(top15_genes_id,gnomad_information$gene_id)]

half_index = which(cumsum(MutVar_pergene$full) > sum(MutVar_pergene$full)/2)[1]

SupplementaryTable3 = data.frame(gene_ID = top15_genes_id,
                                 gene_name = top15_genes_name,
                                 gene_MutVar = MutVar_pergene$full[1:15],
                                 gene_FractionMutVar = MutVar_pergene$full[1:15]/sum(MutVar_pergene$full),
                                 cumulativeFractionMutVar = cumsum(MutVar_pergene$full[1:15]/sum(MutVar_pergene$full)))
write.table(SupplementaryTable3,paste0(path_to_supptables, "SupplementaryTable3.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

mutvar_cumsum = cumsum(MutVar_pergene$full)

#bootstrap
bootstrap_cumsum <- sapply(1:ncol(MutVar_pergene$boot),
                           function(x) {
                             cumsum(MutVar_pergene$boot[,x])
                           })

bootstrap_cumsum_CI <- t(sapply(1:nrow(bootstrap_cumsum),
                                function(x) {
                                  quantile(bootstrap_cumsum[x,],c(0.025,0.975))
                                }))

Polygenicity_plot <- ggplot(mapping = aes(x = 1:length(mutvar_cumsum),
                                          y = mutvar_cumsum)) +
  geom_ribbon(aes(ymin = bootstrap_cumsum_CI[,1], ymax = bootstrap_cumsum_CI[,2]),
              alpha = 0.2,
              color = "grey40",
              fill = "deepskyblue3") +  # Adjust alpha for transparency
  geom_line(color = "black") +
  theme_bhr_gridlines() +
  geom_vline(xintercept = half_index, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Number of Genes",
       y = "Cumulative Mutational Variance\nPTV + Mis2") +
  annotate('text', x = half_index + 10, y = 0.005,
           label = "15 genes explain half\nof mutational variance",
           size = 5, fontface = "italic", hjust = 0)

#Bootstrap enrichment
enrichment_viz_df <- heritability_enrich_autism$enrichment[heritability_enrich_autism$enrichment$prev_factor ==  autism_data$prev_factors[1] &
                                                             grepl("Combined",heritability_enrich_autism$enrichment$name) &
                                                             grepl("Proband",heritability_enrich_autism$enrichment$name) &
                                                             heritability_enrich_autism$enrichment$variant_class %in% c("PTV"),]
enrichment_viz_df$frac_h2_lower95 <- c(burdenEM_models_autism[[1]][[5]]$heritability_output$frach2_CI[1,]
                                       #      ,
                                       #burdenEM_models_autism[[1]][[11]]$heritability_output$frach2_CI[1,]
)
enrichment_viz_df$frac_h2_upper95 <- c(burdenEM_models_autism[[1]][[5]]$heritability_output$frach2_CI[2,]
                                       #,
                                       #burdenEM_models_autism[[1]][[11]]$heritability_output$frach2_CI[2,]
)

#Enrichment plot

enrichment_viz_df$annot_reformat = gsub("_","\n",enrichment_viz_df$annot)

EnrichmentPlot <- ggplot(enrichment_viz_df[enrichment_viz_df$annot != "LOEUF5",],
                         mapping = aes(x =factor(annot_reformat, levels = annot_reformat[variant_class == "PTV"]),
                                       y = h2_enrich,
                                       ymin = h2_enrich_lower,
                                       ymax = h2_enrich_upper,
                                       fill = factor(variant_class, levels = c("PTV")),
                                       label = paste0(round(frac_h2*100),"%")))+
  geom_col(color = "black")+
  geom_errorbar(color = "black",
                # position = position_dodge(width = 0.9),
                width = 0.2)+
  geom_label(mapping = aes(y = -0.25),
             fill = "white",
             #  position = position_dodge2(width = 0.9),
             fontface = 'bold',
             size = 4.5, hjust = 0.4, label.r = unit(0,'lines'))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_fill_manual(values = palette_variantclass)+
  theme_bhr_legend_gridlines()+
  labs(x = "Annotation", y = "Enrichment\nMutational Variance", fill = "Variant Class")+
  # ylim(0,1)+
  guides(fill = "none")+
  theme(axis.text.x = element_text(size = 15))#angle = 20, hjust = 1))



#Aside: lengths of genes
GER_genes = read.table(paste0(resource_dir, "GER_IDs.txt"), header = FALSE)

NC_genes = read.table(paste0(resource_dir,"NC_IDs.txt"), header = FALSE)

length_table = data.frame(gene = gnomad_information_v2$gene_id,
                          length = gnomad_information_v2$cds_length,
                          GER = gnomad_information_v2$gene_id %in% GER_genes$V1,
                          NC= gnomad_information_v2$gene_id %in% NC_genes$V1)
length_table$group = "Neither"
length_table$group[length_table$GER] <- "GER"
length_table$group[length_table$NC] <- "NC"

summary_table = data.frame(class = c("Neither\n(N = 16562)","GER\n(N = 2414)","NC\n(N = 728)"),
                           median = c(median(length_table$length[length_table$group == "Neither"]),
                                      median(length_table$length[length_table$group == "GER"]),
                                      median(length_table$length[length_table$group == "NC"])))

SupplementaryFigure4 <- ggplot(data = summary_table,
                               mapping = aes(x = factor(class, levels = c("Neither\n(N = 16562)","GER\n(N = 2414)","NC\n(N = 728)")),
                                             y = median))+
  geom_col(width = 0.5, color = "black", fill = "white")+
  theme_bhr_gridlines()+
  labs(x = "Gene Set", y = "Coding Length\nMedian")

ggsave(paste0(path_to_figs,"SupplementaryFigure4.pdf"),
       SupplementaryFigure4,device = cairo_pdf,width =5, height = 4,dpi = "retina")

FracCase_autism <- get_fraccase_df(autism_data,
                                   burdenEM_models_autism[[1]][[5]],get_genetic_data(5,autism_data)$genetic_data,
                                   burdenEM_models_autism[[1]][[11]],get_genetic_data(11,autism_data)$genetic_data)

SupplementaryTable4 = FracCase_autism[,c(1:7)]
names(SupplementaryTable4) <- c("Rate Ratio Threshold",
                                "Fraction of Cases with PTV",
                                "Fraction of Cases with PTV Lower95CI",
                                "Fraction of Cases with PTV Upper95CI",
                                "Fraction of Cases with Mis2",
                                "Fraction of Cases with Mis2 Lower95CI",
                                "Fraction of Cases with Mis2 Upper95CI")

write.table(SupplementaryTable4,paste0(path_to_supptables, "SupplementaryTable4.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

FracCases_plot <- ggplot(data = FracCase_autism,
                         mapping = aes(x = RR_thresh,
                                       y = fraccases_greater_PTVplusMis2,
                                       ymin = fraccases_greater_PTVplusMis2_lowerCI,
                                       ymax = fraccases_greater_PTVplusMis2_upperCI))+
  geom_hline(yintercept = 0)+
  geom_ribbon(alpha = 0.2,
              color = "grey40",
              fill = "deepskyblue3")+
  geom_line()+
  #geom_point(size = 2)+
  labs(x = "RR Threshold", y = "Fraction of Cases\nHighly Penetrant PTV/Mis2")+
  theme_bhr_legend_gridlines()+
  geom_vline(xintercept = 5, linetype = "dashed")+
  ylim(0,0.105)+
  annotate("text", x = 5.5, y = 0.095, label = "ACMG\nCriterion",
           size = 5, fontface = "italic", hjust = 0)

#Effective penetrance
peneff_df <- data.frame(peneff = c(burdenEM_models_autism[[1]][[5]]$penetrance$effective_penetrance,
                                   burdenEM_models_autism[[1]][[11]]$penetrance$effective_penetrance,
                                   burdenEM_models_autism[[1]][[17]]$penetrance$effective_penetrance),
                        peneff_lower = c(burdenEM_models_autism[[1]][[5]]$penetrance$effective_penetrance_CI[1],
                                         burdenEM_models_autism[[1]][[11]]$penetrance$effective_penetrance_CI[1],
                                         burdenEM_models_autism[[1]][[17]]$penetrance$effective_penetrance_CI[1]),
                        peneff_upper = c(burdenEM_models_autism[[1]][[5]]$penetrance$effective_penetrance_CI[2],
                                         burdenEM_models_autism[[1]][[11]]$penetrance$effective_penetrance_CI[2],
                                         burdenEM_models_autism[[1]][[17]]$penetrance$effective_penetrance_CI[2]),
                        label = c("PTV","Mis2","Mis1"))

peneff_df$RR_eff = peneff_df$peneff/(autism_data$baseprev * autism_data$prev_factors[1])
peneff_df$RR_eff_lower = peneff_df$peneff_lower/(autism_data$baseprev * autism_data$prev_factors[1])
peneff_df$RR_eff_upper = peneff_df$peneff_upper/(autism_data$baseprev * autism_data$prev_factors[1])

peneff_plot <- ggplot(data = peneff_df,
                      mapping = aes(x = factor(label, levels = c("Mis1","Mis2","PTV")),
                                    y = peneff,
                                    ymin = peneff_lower,
                                    ymax = peneff_upper,
                                    fill = factor(label, levels = c("Mis1","Mis2","PTV"))))+
  geom_hline(yintercept = 0)+
  geom_pointrange(shape = 21, color = "black", size = 1)+
  labs(x = "Variant Class", y = "Effective Penetrance")+
  scale_fill_manual(values = palette_variantclass)+
  theme_bhr_legend_gridlines()+
  theme(axis.text.x = element_text(color = rev(palette_variantclass)[3:5], face = "bold"))+
  ylim(0,1)+
  guides(fill = "none")

#Heritability by dataset
dataset_h2_results = data.frame()

for (study in c("SPARK","ASC","GeneDx")) {

  if (study %in% c("SPARK","ASC")) {
    prev_factor = autism_data$prev_factors[1]

  } else {
    prev_factor = autism_data$prev_factors[2]

  }

  #PTV
  #get index
  ptv_male_index = which(grepl(study, autism_data$loop_vars$names) &
                           grepl("PTV", autism_data$loop_vars$names) &
                           grepl("Male", autism_data$loop_vars$names))
  ptv_female_index = which(grepl(study, autism_data$loop_vars$names) &
                             grepl("PTV", autism_data$loop_vars$names) &
                             grepl("Female", autism_data$loop_vars$names))

  #get scaling factor
  binomial_LL_ptv_male = binomial_analysis(model = burdenEM_models_autism[[1]][[5]],
                                           genetic_data_total = get_genetic_data(5,autism_data)$genetic_data,
                                           genetic_data_subsample = get_genetic_data(ptv_male_index,autism_data)$genetic_data)

  binomial_LL_ptv_female = binomial_analysis(model = burdenEM_models_autism[[1]][[5]],
                                             genetic_data_total = get_genetic_data(5,autism_data)$genetic_data,
                                             genetic_data_subsample = get_genetic_data(ptv_female_index,autism_data)$genetic_data)
  #get modified h2
  h2_output_ptv_male = get_scaled_heritability(burdenEM_models_autism[[1]][[5]],
                                               gamma_scaling_factor = binomial_LL_ptv_male$MLE,
                                               genetic_data = get_genetic_data(5,autism_data)$genetic_data,
                                               prevalence = autism_data$loop_vars$prevalences[ptv_male_index] * prev_factor,
                                               heritability_scaling_factor = autism_data$ptv_scale_factor)

  h2_output_ptv_female = get_scaled_heritability(burdenEM_models_autism[[1]][[5]],
                                                 gamma_scaling_factor = binomial_LL_ptv_female$MLE,
                                                 genetic_data = get_genetic_data(5,autism_data)$genetic_data,
                                                 prevalence = autism_data$loop_vars$prevalences[ptv_female_index] * prev_factor,
                                                 heritability_scaling_factor = autism_data$ptv_scale_factor)

  #Mis2
  #get modified h2
  h2_output_mis2_male = get_scaled_heritability(burdenEM_models_autism[[1]][[11]],
                                                gamma_scaling_factor = binomial_LL_ptv_male$MLE,
                                                genetic_data = get_genetic_data(11,autism_data)$genetic_data,
                                                prevalence = autism_data$loop_vars$prevalences[ptv_male_index] * prev_factor)

  h2_output_mis2_female = get_scaled_heritability(burdenEM_models_autism[[1]][[11]],
                                                  gamma_scaling_factor = binomial_LL_ptv_female$MLE,
                                                  genetic_data = get_genetic_data(11,autism_data)$genetic_data,
                                                  prevalence = autism_data$loop_vars$prevalences[ptv_female_index] * prev_factor)


  #combined

  h2_male_combined = h2_output_ptv_male$h2 + h2_output_mis2_male$h2
  h2_female_combined = h2_output_ptv_female$h2 + h2_output_mis2_female$h2

  bootstrap_h2_combined_male= sapply(1:length(h2_output_ptv_male$bootstrap_h2_ests),
                                     function(x) h2_output_ptv_male$bootstrap_h2_ests[x] + h2_output_mis2_male$bootstrap_h2_ests[x])
  h2_male_combined_CI = quantile(bootstrap_h2_combined_male,c(0.025,0.975))

  bootstrap_h2_combined_female= sapply(1:length(h2_output_ptv_female$bootstrap_h2_ests),
                                       function(x) h2_output_ptv_female$bootstrap_h2_ests[x] + h2_output_mis2_female$bootstrap_h2_ests[x])
  h2_female_combined_CI = quantile(bootstrap_h2_combined_female,c(0.025,0.975))

  #fraction of cases with large effect variant
  male_fraccase <- get_fraccase_df(autism_data,burdenEM_models_autism[[1]][[5]],
                                   get_genetic_data(5,autism_data)$genetic_data,
                                   burdenEM_models_autism[[1]][[11]],
                                   get_genetic_data(11,autism_data)$genetic_data,
                                   gamma_scaling_factor = binomial_LL_ptv_male$MLE,
                                   RR_range = c(5))

  female_fraccase <- get_fraccase_df(autism_data,burdenEM_models_autism[[1]][[5]],
                                     get_genetic_data(5,autism_data)$genetic_data,
                                     burdenEM_models_autism[[1]][[11]],
                                     get_genetic_data(11,autism_data)$genetic_data,
                                     gamma_scaling_factor = binomial_LL_ptv_female$MLE,
                                     RR_range = c(5))

  iter_df <- data.frame(study = study,
                        sex = c("Male","Female"),
                        prev = c(autism_data$loop_vars$prevalences[ptv_male_index] * prev_factor,
                                 autism_data$loop_vars$prevalences[ptv_female_index] * prev_factor),
                        gamma_scaling_factor = c(binomial_LL_ptv_male$MLE,binomial_LL_ptv_female$MLE),
                        h2_ptv = c(h2_output_ptv_male$h2,h2_output_ptv_female$h2),
                        h2_ptv_lower = c(h2_output_ptv_male$h2_lower,h2_output_ptv_female$h2_lower),
                        h2_ptv_upper = c(h2_output_ptv_male$h2_upper,h2_output_ptv_female$h2_upper),
                        h2_mis2 = c(h2_output_mis2_male$h2,h2_output_mis2_female$h2),
                        h2_mis2_lower = c(h2_output_mis2_male$h2_lower,h2_output_mis2_female$h2_lower),
                        h2_mis2_upper = c(h2_output_mis2_male$h2_upper,h2_output_mis2_female$h2_upper),
                        h2_combined = c(h2_male_combined,h2_female_combined),
                        h2_combined_lower = c(h2_male_combined_CI[1],h2_female_combined_CI[1]),
                        h2_combined_upper = c(h2_male_combined_CI[2],h2_female_combined_CI[2]),
                        fraccase_RR5 = c(male_fraccase$fraccases_greater_PTVplusMis2[male_fraccase$RR_thresh == 5],
                                         female_fraccase$fraccases_greater_PTVplusMis2[female_fraccase$RR_thresh == 5]),
                        fraccase_RR5_lower = c(male_fraccase$fraccases_greater_PTVplusMis2_lowerCI[male_fraccase$RR_thresh == 5],
                                               female_fraccase$fraccases_greater_PTVplusMis2_lowerCI[female_fraccase$RR_thresh == 5]),
                        fraccase_RR5_upper = c(male_fraccase$fraccases_greater_PTVplusMis2_upperCI[male_fraccase$RR_thresh == 5],
                                               female_fraccase$fraccases_greater_PTVplusMis2_upperCI[female_fraccase$RR_thresh == 5]))


  iter_df$study[iter_df$study == "GeneDx"] <- "GeneDX"

  dataset_h2_results = rbind(dataset_h2_results, iter_df)

}

dataset_mutvar_maindisplay = ggplot(data = dataset_h2_results,
                                    mapping = aes(x = factor(study,levels = c("SPARK","ASC","GeneDX")),
                                                  y = h2_combined,
                                                  ymin = h2_combined_lower,
                                                  ymax = h2_combined_upper,
                                                  fill = sex))+
  geom_hline(yintercept = 0)+
  geom_pointrange(position = position_dodge(width = 0.5), color = "black", shape = 21, size = 1)+
  theme_bhr_legend_gridlines()+
  scale_fill_manual(values = c("Male" = "chartreuse3", "Female" = "#9467bd"))+
  labs(x = "Study", y = "Mutational Variance\nPTV + Mis2", color = "Sex")+
  guides(fill = "none")

dataset_fraccase_maindisplay = ggplot(data = dataset_h2_results,
                                      mapping = aes(x = factor(study,levels = c("SPARK","ASC","GeneDX")),
                                                    y = fraccase_RR5,
                                                    ymin = fraccase_RR5_lower,
                                                    ymax = fraccase_RR5_upper,
                                                    fill = sex))+
  geom_hline(yintercept = 0)+
  geom_col(color = "black", position = position_dodge(width = 1), width = 0.8)+
  geom_errorbar(position = position_dodge(width = 1), width = 0.1)+
  theme_bhr_legend_gridlines()+
  scale_fill_manual(values = c("Male" = "chartreuse3", "Female" = "#9467bd"))+
  labs(x = "Study", y = "Fraction of Cases\nRR>5, PTV + Mis2", fill = "Sex")


spacer_labeled <- wrap_elements(grid::textGrob("", gp = grid::gpar(fontsize = 25)))

combined_fig1 = (spacer_labeled | Proband_Vmu_plot) / (Polygenicity_plot | EnrichmentPlot) / (FracCases_plot | peneff_plot) / (dataset_fraccase_maindisplay | dataset_mutvar_maindisplay)+
  plot_annotation(tag_levels = 'A', tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 25))


ggsave(paste0(path_to_figs,"Fig1_all.pdf"),
       combined_fig1,device = cairo_pdf,width =12, height = 15,dpi = "retina")

SupplementaryTable5 <- dataset_h2_results[,c(1:10,
                                             13:15)]
names(SupplementaryTable5) <- c("Study",
                                "Sex",
                                "Prevalence",
                                "Scaling_Factor",
                                "MutVar_PTV",
                                "MutVar_PTV_Lower95CI",
                                "MutVar_PTV_Upper95CI",
                                "MutVar_Mis2",
                                "MutVar_Mis2_Lower95CI",
                                "MutVar_Mis2_Upper95CI",
                                "Fraction of Cases with RR>5 PTV/Mis2",
                                "Fraction of Cases with RR>5 PTV/Mis2 Lower95CI",
                                "Fraction of Cases with RR>5 PTV/Mis2 Upper95CI")
write.table(SupplementaryTable5,paste0(path_to_supptables, "SupplementaryTable5.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#DDD Comparisons

#Mutational variance comparison


heritability_enrich_dd_forcomparison = heritability_enrich_dd$heritability[heritability_enrich_dd$heritability$prev_factor == 1 &
                                                                             !(grepl("Male",heritability_enrich_dd$heritability$name)) &
                                                                             !(grepl("Female",heritability_enrich_dd$heritability$name)),]
heritability_enrich_dd_forcomparison$dx = "DD"
heritability_enrich_autism_forcomparison <- heritability_enrich_autism$heritability[heritability_enrich_autism$heritability$prev_factor ==  autism_data$prev_factors[1] &
                                                                                      grepl("Combined",heritability_enrich_autism$heritability$name) &
                                                                                      grepl("Proband",heritability_enrich_autism$heritability$name),]
heritability_enrich_autism_forcomparison$dx = "Autism"

mutvar_compare_df <- rbind(heritability_enrich_dd_forcomparison,
                           heritability_enrich_autism_forcomparison)

mutvar_compare_plot <- ggplot(data = mutvar_compare_df,
                              mapping = aes(x= factor(variant_class, levels = c("Syn","Mis0","Mis1","Mis2","PTV")),
                                            y = h2,
                                            ymin = h2_lower,
                                            ymax = h2_upper,
                                            fill = dx
                              ))+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2"))+
  geom_pointrange(position = position_dodge2(width = 0.4), shape = 21, size = 1)+
  theme_bhr_legend_gridlines()+
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale")+
  theme(strip.text.y = element_text(size = 15, angle = 270),
        axis.text.x = element_text(color = rev(palette_variantclass), face = "bold"))+
  ylim(0,0.045)+
  guides(fill = "none")+
  geom_rect(mapping = aes(xmin = 0.85, xmax = 2.4, ymin = 0.030, ymax = 0.044), color = "black", fill = "white")+
  annotate("text", x = 1, y = 0.04, label = "Autism", color = "darkgreen",
           size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = 0.034, label = "DD", color = "orchid3",
           size = 6, fontface = "bold", hjust = 0)



#Frac cases large effect
FracCase_dd<- get_fraccase_df(kaplanis_data,
                              burdenEM_models_DDD[[1]][[1]],get_genetic_data(1,kaplanis_data)$genetic_data,
                              burdenEM_models_DDD[[1]][[2]],get_genetic_data(2,kaplanis_data)$genetic_data,boot = TRUE)
SupplementaryTable8 = FracCase_dd[,c(1:7)]
names(SupplementaryTable8) <- c("Rate Ratio Threshold",
                                "Fraction of Cases with PTV",
                                "Fraction of Cases with PTV Lower95CI",
                                "Fraction of Cases with PTV Upper95CI",
                                "Fraction of Cases with Mis2",
                                "Fraction of Cases with Mis2 Lower95CI",
                                "Fraction of Cases with Mis2 Upper95CI")
write.table(SupplementaryTable8,paste0(path_to_supptables, "SupplementaryTable8.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

FracCase_autism$dx = "Autism"
FracCase_dd$dx = "DD"

compare_fraccase_dx <- rbind(FracCase_autism,FracCase_dd)

write.csv(compare_fraccase_dx,
          file =paste0('~/Mirror/oconnor_rotation/rare_dn_h2/outputs/FracCase_output_',current_date),
          quote = FALSE,
          row.names = FALSE)

FracCases_compare_plot <- ggplot(data = compare_fraccase_dx,
                                 mapping = aes(x = RR_thresh,
                                               y = fraccases_greater_PTVplusMis2,
                                               ymin = fraccases_greater_PTVplusMis2_lowerCI,
                                               ymax = fraccases_greater_PTVplusMis2_upperCI,
                                               color = dx,
                                               fill = dx))+
  scale_color_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2"))+
  scale_fill_manual(values = c("Autism" = "forestgreen", "DD" = "orchid1"))+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_ribbon(alpha = 0.2,
              color = "grey40")+
  labs(x = "RR Threshold", y = "Fraction of Cases\nHighly Penetrant PTV/Mis2")+
  theme_bhr_legend_gridlines()+
  ylim(0,0.22)+
  geom_vline(xintercept = 5, linetype = "dashed")+
  guides(color = "none", fill = "none")+
  annotate("text", x = 5.5, y = 0.196, label = "ACMG\nCriterion",
           size = 5, fontface = "italic", hjust = 0)


#Effective rate ratio
peneff_df_dd <- data.frame(peneff = c(burdenEM_models_DDD[[1]][[1]]$penetrance$effective_penetrance,
                                      burdenEM_models_DDD[[1]][[2]]$penetrance$effective_penetrance,
                                      burdenEM_models_DDD[[1]][[3]]$penetrance$effective_penetrance),
                           peneff_lower = c(burdenEM_models_DDD[[1]][[1]]$penetrance$effective_penetrance_CI[1],
                                            burdenEM_models_DDD[[1]][[2]]$penetrance$effective_penetrance_CI[1],
                                            burdenEM_models_DDD[[1]][[3]]$penetrance$effective_penetrance_CI[1]),
                           peneff_upper = c(burdenEM_models_DDD[[1]][[1]]$penetrance$effective_penetrance_CI[2],
                                            burdenEM_models_DDD[[1]][[2]]$penetrance$effective_penetrance_CI[2],
                                            burdenEM_models_DDD[[1]][[3]]$penetrance$effective_penetrance_CI[2]),
                           label = c("PTV","Mis2","Mis1"))

peneff_df_dd$RR_eff = peneff_df_dd$peneff/(kaplanis_data$baseprev * kaplanis_data$prev_factors[1])
peneff_df_dd$RR_eff_lower = peneff_df_dd$peneff_lower/(kaplanis_data$baseprev * kaplanis_data$prev_factors[1])
peneff_df_dd$RR_eff_upper = peneff_df_dd$peneff_upper/(kaplanis_data$baseprev * kaplanis_data$prev_factors[1])

peneff_df$dx = "Autism"
peneff_df_dd$dx = "DD"

peneff_df_compare <- rbind(peneff_df,peneff_df_dd)

RReff_compare_plot <- ggplot(data = peneff_df_compare,
                             mapping = aes(x = factor(label, levels = c("Mis1","Mis2","PTV")),
                                           y = RR_eff,
                                           ymin = RR_eff_lower,
                                           ymax = RR_eff_upper,
                                           fill = dx))+
  scale_fill_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2"))+
  geom_hline(yintercept = 0)+
  geom_pointrange(position = position_dodge2(width = 0.4), shape = 21, size = 1)+
  labs(x = "Variant Class", y = "Effective Rate Ratio")+
  theme_bhr_legend_gridlines()+
  theme(axis.text.x = element_text(color = rev(palette_variantclass)[-c(1,2)], face = "bold"))+
  ylim(0,42) +
  guides(fill = "none")

FigDD <- mutvar_compare_plot + FracCases_compare_plot + RReff_compare_plot+
  plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave(paste0(path_to_figs,"Figure4.pdf"),
       FigDD,device = cairo_pdf,width =14, height = 4,dpi = "retina")

#Forecasting
effectsize_thresh_allgenes = count_genes_by_effectsize(burdenEM_models_autism[[1]][[5]],
                                                       get_genetic_data(5,autism_data)$genetic_data)
library(dplyr)

Fu_sigtests <- read.csv(paste0(resource_dir,"fu_files/significance_tests_fu.csv"))
ASD72_genes <- Fu_sigtests$gene_id[Fu_sigtests$ASD72][!is.na(Fu_sigtests$gene_id[Fu_sigtests$ASD72])]
ASD185_genes <-  Fu_sigtests$gene_id[Fu_sigtests$ASD185][!is.na(Fu_sigtests$gene_id[Fu_sigtests$ASD185])]

effectsize_thresh_Fu72 = count_genes_by_effectsize(burdenEM_models_autism[[1]][[5]],
                                                   get_genetic_data(5,autism_data)$genetic_data,
                                                   genes_to_analyze = ASD72_genes)

effectsize_thresh_Fu185 = count_genes_by_effectsize(burdenEM_models_autism[[1]][[5]],
                                                    get_genetic_data(5,autism_data)$genetic_data,
                                                    genes_to_analyze = ASD185_genes)

library(MetBrewer)
forecast_df_average_expand <- read.table(paste0(output_dir,"forecasting_output_avgs_Mar42025.tsv"),
                                         sep = "\t", header = TRUE)
SupplementaryTable6 = forecast_df_average_expand[forecast_df_average_expand$threshold == 0, -1]
names(SupplementaryTable6) = c("New Data Dataset",
                               "N New Cases",
                               "Number FDR Significant Genes (mean)",
                               "Number FDR Significant Genes (sd)",
                               "Number Bonferroni Significant Genes (mean)",
                               "Number Bonferroni Significant Genes (sd)")
write.table(SupplementaryTable6,paste0(path_to_supptables, "SupplementaryTable6.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


palette = MetBrewer::met.brewer("Juarez",3,"discrete")

ForecastingPlot <- ggplot(data = forecast_df_average_expand[forecast_df_average_expand$threshold == 0,],
                          mapping = aes(x = 38088 + N_new_case,
                                        y = count_bonferroni_mean,
                                        ymin = count_bonferroni_mean - count_bonferroni_sd,
                                        ymax = count_bonferroni_mean + count_bonferroni_sd,
                                        color = factor(dataset, levels = c("SPARK","ASC","GeneDx")),
                                        fill = factor(dataset, levels = c("SPARK","ASC","GeneDx"))))+
  #geom_errorbar(width = 1000)+
  #geom_point()+
  geom_ribbon(alpha = 0.2, linetype = "blank") +
  geom_line()+
  geom_line(mapping = aes(x = 38088 + N_new_case,
                          y = count_bonferroni_mean - count_bonferroni_sd,
                          color = factor(dataset, levels = c("SPARK","ASC","GeneDx"))),
            alpha = 0.3)+
  geom_line(mapping = aes(x = 38088 + N_new_case,
                          y = count_bonferroni_mean + count_bonferroni_sd,
                          color = factor(dataset, levels = c("SPARK","ASC","GeneDx"))),
            alpha = 0.3)+
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette)+
  # facet_grid(cols = vars(threshold))+
  #ylim(0,120)+
  theme_bhr_legend_gridlines()+
  labs(x = "Sample Size", y = "Number of Significant Genes\n(PTV)", fill = "New Data\nSource")+
  guides(alpha = "none",  # Remove alpha legend
         color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))

ggsave(paste0(path_to_figs,"Figure3_bottom.pdf"),
       ForecastingPlot,device = cairo_pdf,width =9, height = 4,dpi = "retina")


