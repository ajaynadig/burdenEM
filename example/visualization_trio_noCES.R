source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/burdenEM_trio.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/EM.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/estimate_heritability.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/io.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/likelihoods.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/model.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/secondary_analysis_functions.R")

source("~/Mirror/oconnor_rotation/rare_dn_h2/scripts/set_up.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/scripts/set_up_kaplanis.R")


resource_dir <- "~/Mirror/oconnor_rotation/rare_dn_h2/github/resources/"
output_dir = "~/Mirror/oconnor_rotation/rare_dn_h2/github/outputs/"
path_to_supptables <- paste0(output_dir, "SuppTables/")
path_to_figs <-paste0(output_dir, "Figures/")


load(paste0(output_dir, "noCESset2_models_autism_Apr25.Rdata"))
load(paste0(output_dir,"noCESset2_models_ddd_Apr25.Rdata"))


current_date <- format(Sys.Date(), "%b%y")


#collate the heritability results into a summary dataframe

heritability_enrich_autism = heritability_enrichment_table(autism_data,burdenEM_models_autism)
heritability_enrich_dd = heritability_enrichment_table(kaplanis_data, burdenEM_models_DDD)

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
