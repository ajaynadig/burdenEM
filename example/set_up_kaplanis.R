#Helper function from kyle
library(SummarizedExperiment)
library(ggplot2)
library(data.table)
library(patchwork)
library(readr)

resource_dir <- "~/Mirror/oconnor_rotation/rare_dn_h2/github/resources/"
setwd(resource_dir)

theme_bhr_legend_gridlines <- function(){
  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15,  color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(fill = "white"))
}

theme_bhr_gridlines <- function(){
  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15,  color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          legend.position = "None",
          legend.direction = "horizontal",
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(fill = "white"))

}

resolve_duplicate_names_prior_to_collecting_by_gene_id <- function(input_table) {
  input_table$Gene_ID[input_table$Gene == 'HIST1H4F'] = 'ENSG00000274618'
  input_table$LOEUF[input_table$Gene == 'HIST1H4F'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'HIST1H4F'] = NA

  input_table$Gene_ID[input_table$Gene == 'ATXN7'] = 'ENSG00000163635'
  input_table$LOEUF[input_table$Gene == 'ATXN7'] = 0.328
  input_table$LOEUF_bin[input_table$Gene == 'ATXN7'] = 1

  input_table$Gene_ID[input_table$Gene == 'PRSS50'] = 'ENSG00000283706'
  input_table$LOEUF[input_table$Gene == 'PRSS50'] = 0.99
  input_table$LOEUF_bin[input_table$Gene == 'PRSS50'] = 5

  input_table$Gene_ID[input_table$Gene == 'CCDC39'] = 'ENSG00000284862'
  input_table$LOEUF[input_table$Gene == 'CCDC39'] = 0.754
  input_table$LOEUF_bin[input_table$Gene == 'CCDC39'] = 4

  input_table$Gene_ID[input_table$Gene == 'HSPA14'] = 'ENSG00000187522'
  input_table$LOEUF[input_table$Gene == 'HSPA14'] = 0.395
  input_table$LOEUF_bin[input_table$Gene == 'HSPA14'] = 1

  input_table$Gene_ID[input_table$Gene == 'IGF2'] = 'ENSG00000167244'
  input_table$LOEUF[input_table$Gene == 'IGF2'] = 1.131
  input_table$LOEUF_bin[input_table$Gene == 'IGF2'] = 6

  input_table$Gene_ID[input_table$Gene == 'PDE11A'] = 'ENSG00000128655'
  input_table$LOEUF[input_table$Gene == 'PDE11A'] = 1.404
  input_table$LOEUF_bin[input_table$Gene == 'PDE11A'] = 7

  # These occur across multiple chromosomes
  temp1 = subset(input_table, Gene %in% c('RF00017', 'RF00019'))
  temp2 = subset(input_table, !Gene %in% c('RF00017', 'RF00019'))

  if(nrow(temp1) > 0) {
    temp1$Gene = paste0(temp1$Gene, '_chr', temp1$Chrom)
    input_table = rbind(temp1, temp2)
  } else {
    input_table = temp2
  }

  # RF00019 requires further parsing; here take the one that VEP uses most often
  input_table$Gene_ID[input_table$Gene == 'RF00019_chr3'] = 'ENSG00000212392'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr3'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr3'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr6'] = 'ENSG00000200314'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr6'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr6'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr7'] = 'ENSG00000201913'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr7'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr7'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr12'] = 'ENSG00000207176'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr12'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr12'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr16'] = 'ENSG00000199668'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr16'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr16'] = NA

  return(input_table)
}

#DDD experiment

kaplanis_data <- read.table("kaplanis_files/kaplanis_variants_annotated_2024-05-15.txt", header = TRUE)
sex_info <- read.table("kaplanis_files/fordist_joint_dnm_ID_sex_2019_08_30.txt", header = TRUE)

#get the PTV nonindel : PTV ratio

num_ptv_nonindel = sum((kaplanis_data$isPTV & !kaplanis_data$isIndel)[!is.na(kaplanis_data$isPTV & !kaplanis_data$isIndel)])
num_ptv = sum(kaplanis_data$isPTV[!is.na(kaplanis_data$isPTV & !kaplanis_data$isIndel)])

kaplanis_ptv_scale_factor = num_ptv/num_ptv_nonindel

mutation_rate_table = read.table("new_count_data/ASD_gene_table_w_bespoke_mutation_rates_2024-07-24.txt",
                                 header = TRUE,
                                 sep = "\t")
kaplanis_data_fixedID <- resolve_duplicate_names_prior_to_collecting_by_gene_id(kaplanis_data)
IDs <- unique(kaplanis_data_fixedID$Gene_ID)
names = kaplanis_data_fixedID$Gene[match(IDs,kaplanis_data_fixedID$Gene_ID)]

id_match = IDs %in% mutation_rate_table$Gene_ID
name_match = names %in% mutation_rate_table$Gene

missing_ids = IDs[!(id_match | name_match)]

#First, construct the ColData


studies = c("DDD",
            "GDX",
            "RUMC")

colData <- data.frame()

for (study in studies) {
  dataset = sex_info[sex_info$study == study,]
  summary_df = data.frame(study = study,
                          N_Proband = nrow(dataset),
                          N_MaleProband = sum(dataset$sex == "M"),
                          N_FemaleProband = sum(dataset$sex == "F"))

  colData = rbind(colData,summary_df)
}

rownames(colData) = colData$study

#Now, define the rowData
genes_with_dn = unique(kaplanis_data_fixedID$Gene_ID)[!(unique(kaplanis_data_fixedID$Gene_ID) %in% missing_ids)]
rowData = data.frame(gene_id = genes_with_dn,
                     gene = kaplanis_data_fixedID$Gene[match(genes_with_dn, kaplanis_data_fixedID$Gene_ID)],
                     LOEUF = kaplanis_data_fixedID$LOEUF[match(genes_with_dn, kaplanis_data_fixedID$Gene_ID)],
                     Chrom = kaplanis_data_fixedID$Chrom[match(genes_with_dn, kaplanis_data_fixedID$Gene_ID)])

#find indices of genes_with_dn in mutation_rate_table
mutation_table_indices <- sapply(genes_with_dn,
                                 function(id) {
                                   # print(id)
                                   matched_genes = unique(kaplanis_data_fixedID$Gene[!is.na(kaplanis_data_fixedID$Gene) & kaplanis_data_fixedID$Gene_ID == id])

                                   if (id %in% mutation_rate_table$Gene_ID) {
                                     return(which(mutation_rate_table$Gene_ID == id))
                                   } else if (any(matched_genes %in% mutation_rate_table$Gene)) {
                                     matched_gene = matched_genes[matched_genes %in% mutation_rate_table$Gene ]
                                     return(which(mutation_rate_table$Gene ==  matched_gene))
                                   }

                                   return(NA)
                                 }
)

rowData = cbind(rowData, mutation_rate_table[mutation_table_indices,-c(1:5)])
rownames(rowData) = rowData$gene_id


#initialize the SummarizedExperiment

kaplanis_counts = SummarizedExperiment(rowData = rowData,
                                       colData = colData)


#Now, tally the counts
tally_function_kaplanis <- function(variantdata,gene_ids, variant_class, Sex) {
  variantdata = variantdata[!is.na(variantdata$Gene_ID),]
  if (variant_class == "PTV") {
    variantdata_subset = variantdata[(variantdata$isPTV)  & !variantdata$isIndel & variantdata$Sex == Sex,]
    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_subset$Gene_ID == id)
                  }))

  } else if (grepl("Mis",variant_class)) {
    num_criteria = parse_number(variant_class)
    variantdata_missense = variantdata[variantdata$isMis & variantdata$Sex == Sex, ]
    alphamissense_filter = as.numeric(!is.na(variantdata_missense$am_pathogenicity) & variantdata_missense$am_pathogenicity >= 0.97)
    MPC2_filter = as.numeric(!is.na(variantdata_missense$MPC_v2) & variantdata_missense$MPC_v2 >= 2)

    missense_filter = (alphamissense_filter + MPC2_filter) == num_criteria
    variantdata_missense_filter = variantdata_missense[missense_filter,]

    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_missense_filter$Gene_ID == id)
                  }))
  } else if (variant_class == "Syn") {
    variantdata_subset = variantdata[variantdata$isSyn & variantdata$Sex == Sex,]
    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_subset$Gene_ID == id)
                  }))
  }
}

count_subset_names = c("PTV_Proband_M","PTV_Proband_F",
                       "Mis2_Proband_M","Mis2_Proband_F",
                       "Mis1_Proband_M","Mis1_Proband_F",
                       "Mis0_Proband_M","Mis0_Proband_F",
                       "Syn_Proband_M","Syn_Proband_F")

subset_variant_class = stringr::str_split_i(count_subset_names,"_",1)
subset_Sex = stringr::str_split_i(count_subset_names,"_",3)

for (subsetnum in 1:length(count_subset_names)) {
  print(count_subset_names[subsetnum])
  subset_count_matrix = matrix(data = NA, nrow = nrow(rowData), ncol = nrow(colData))
  rownames(subset_count_matrix) = rownames(rowData)
  colnames(subset_count_matrix) = rownames(colData)


  for (num_count in 1:length(studies)) {
    print(num_count)
    kaplanis_data_fixedID_study = kaplanis_data_fixedID[kaplanis_data_fixedID$Study == studies[num_count],]

    subset_count_matrix[,num_count] = tally_function_kaplanis(kaplanis_data_fixedID_study,
                                                              rownames(subset_count_matrix),
                                                              subset_variant_class[subsetnum],
                                                              subset_Sex[subsetnum])
  }

  assays(kaplanis_counts)[[count_subset_names[subsetnum]]] = subset_count_matrix
}

#Define the rowData and the colData

library(SummarizedExperiment)
library(readr)

#subset to genes with complete mutation rate information
kaplanis_counts = kaplanis_counts[complete.cases(mutation_rate_table[mutation_table_indices,-c(1:5)]),]

#subset to genes with available LOEUF information
kaplanis_counts = kaplanis_counts[!is.na(rowData(kaplanis_counts)$LOEUF),]

#subset to autosomal genes
kaplanis_counts = kaplanis_counts[!(rowData(kaplanis_counts)$Chrom %in% c("X","Y")),]

#make a one-hot encoded LOEUF annotation matrix
quantiles <- quantile(rowData(kaplanis_counts)$LOEUF, probs = seq(0, 1, by = 0.2))

# Use cut to map each value to its quantile range
quantile_assignments <- cut(rowData(kaplanis_counts)$LOEUF, breaks = quantiles, include.lowest = TRUE,
                            labels = c(1:5))

for (bin in 1:5) {
  rowData(kaplanis_counts)[paste0("LOEUF",bin)] =  as.numeric(quantile_assignments == bin)
}

#Get the posterior mutation rate correction factors based on gnomad counts
#First, let's get the mutation rate correction factors
gnomad_information <- read.table("~/Mirror/oconnor_rotation/rare_dn_h2/repo/burdenEM/gnomad_reference/gnomad.v4.1.constraint_metrics.tsv", header = TRUE)

gnomad_input_df <- data.frame(case_count = gnomad_information$syn.obs[gnomad_information$canonical == "true" & grepl("ENSG",gnomad_information$gene_id)],
                              expected_count = gnomad_information$syn.exp[gnomad_information$canonical == "true"& grepl("ENSG",gnomad_information$gene_id)],
                              gene_id = gnomad_information$gene_id[gnomad_information$canonical == "true"& grepl("ENSG",gnomad_information$gene_id)])
gnomad_input_df = gnomad_input_df[complete.cases(gnomad_input_df),]

mut_calibration_syn <- burdenEM_trio(gnomad_input_df,
                                     features = NULL,
                                     component_endpoints = seq(-2,2,length.out = 31),
                                     heritability_est = FALSE,
                                     null_sim = FALSE,
                                     max_iter = 1000,
                                     bootstrap = FALSE,
                                     return_likelihood = TRUE,
                                     estimate_posteriors = TRUE,
                                     estimate_effective_penetrance = FALSE)

rowData(kaplanis_counts)$PosteriorMuCorrectionFactor = mut_calibration_syn$posterior_gene_estimates$Posterior_Mean[match(rowData(kaplanis_counts)$gene_id,
                                                                                                                         gnomad_input_df$gene_id)]
ggplot(mapping = aes(x = rowData(kaplanis_counts)$PosteriorMuCorrectionFactor))+
  geom_histogram(color = "black", fill = "white")+
  theme_bhr_legend_gridlines()+
  labs(x = "Posterior Mean\nMutation Rate Correction", y = "Count")

kaplanis_counts = kaplanis_counts[!is.na(rowData(kaplanis_counts)$PosteriorMuCorrectionFactor),]

#Define the "baseline model" features
#get the class-wise mutation rates
LOEUF = as.matrix(rowData(kaplanis_counts)[,13:17])

mu_byclass = as.matrix(rowData(kaplanis_counts)[,5:9])
cumulative_mu = rowSums(mu_byclass) * rowData(kaplanis_counts)$PosteriorMuCorrectionFactor
cumulative_mu_bin <- cbind(as.numeric(cumulative_mu > quantile(cumulative_mu, 0.8)),
                           as.numeric(cumulative_mu <= quantile(cumulative_mu, 0.8)))


kaplanis_features <- cbind(LOEUF[,1] * cumulative_mu_bin[,1],
                           LOEUF[,1] * cumulative_mu_bin[,2],
                           LOEUF[,2:5])
colnames(kaplanis_features) = c("LOEUF1_mu1",
                                "LOEUF1_mu2",
                                "LOEUF2",
                                "LOEUF3",
                                "LOEUF4",
                                "LOEUF5")


#Define looping variables for main burdenEM loop
kaplanis_loop_vars = list()


kaplanis_loop_vars$names = c("Proband PTV",
                             "Proband Mis2",
                             "Proband Mis1",
                             "Proband Mis0",
                             "Proband Syn",
                             "Proband Male PTV",
                             "Proband Female PTV",
                             "Proband Male Mis2",
                             "Proband Female Mis2",
                             "Proband Male Mis1",
                             "Proband Female Mis1",
                             "Proband Male Mis0",
                             "Proband Female Mis0",
                             "Proband Male Syn",
                             "Proband Female Syn")
kaplanis_loop_vars$N_subset = c(rep(sum(kaplanis_counts$N_Proband),5),
                                rep(c(sum(kaplanis_counts$N_MaleProband),
                                      sum(kaplanis_counts$N_FemaleProband)),
                                    5))
kaplanis_loop_vars$subsets = list(c("PTV_Proband_M","PTV_Proband_F"),
                                  c("Mis2_Proband_M","Mis2_Proband_F"),
                                  c("Mis1_Proband_M","Mis1_Proband_F"),
                                  c("Mis0_Proband_M","Mis0_Proband_F"),
                                  c("Syn_Proband_M","Syn_Proband_F"),
                                  c("PTV_Proband_M"),
                                  c("PTV_Proband_F"),
                                  c("Mis2_Proband_M"),
                                  c("Mis2_Proband_F"),
                                  c("Mis1_Proband_M"),
                                  c("Mis1_Proband_F"),
                                  c("Mis0_Proband_M"),
                                  c("Mis0_Proband_F"),
                                  c("Syn_Proband_M"),
                                  c("Syn_Proband_F"))


kaplanis_loop_vars$mutation_rate = c("mu_snp_PTV",
                                     "mu_snp_Mis2",
                                     "mu_snp_Mis1",
                                     "mu_snp_Mis0",
                                     "mu_snp_Syn",
                                     "mu_snp_PTV",
                                     "mu_snp_PTV",
                                     "mu_snp_Mis2",
                                     "mu_snp_Mis2",
                                     "mu_snp_Mis1",
                                     "mu_snp_Mis1",
                                     "mu_snp_Mis0",
                                     "mu_snp_Mis0",
                                     "mu_snp_Syn",
                                     "mu_snp_Syn")


baseprev = 0.01

kaplanis_loop_vars$prevalences = rep(baseprev,15)

N_SPARK = 11052 + 465 + 6543
N_ASC = 1109 + 7291 + 279
N_GeneDX = 11349

weight_SPARK = N_SPARK/(N_SPARK + N_ASC + N_GeneDX)
weight_ASC = N_ASC/(N_SPARK + N_ASC + N_GeneDX)
weight_GeneDX = N_GeneDX/(N_SPARK + N_ASC + N_GeneDX)

prev_weighted = (weight_SPARK * 0.0276) + (weight_ASC * 0.0276) + (weight_GeneDX * 0.01)


#prev_factors = c(0.5,0.75,1,1.25,1.5)
prev_factors = c(1,prev_weighted/0.01)

kaplanis_data <- list(counts = kaplanis_counts,
                      features = kaplanis_features,
                      loop_vars = kaplanis_loop_vars,
                      ptv_scale_factor = kaplanis_ptv_scale_factor,
                      baseprev = baseprev,
                      prev_factors = prev_factors)

ddd_genes_across <-  lapply(1:length(kaplanis_data$loop_vars$subsets),
                            function(i) {

                              processed_input = get_genetic_data(i,kaplanis_data)

                              input_df = processed_input$genetic_data
                              features = processed_input$features

                              return(rownames(input_df))
                            })

ddd_genes_consensus = Reduce("intersect",ddd_genes_across)
kaplanis_data$features = kaplanis_data$features[rownames(kaplanis_data$counts) %in% ddd_genes_consensus,]
kaplanis_data$counts = kaplanis_data$counts[rownames(kaplanis_data$counts) %in% ddd_genes_consensus,]


