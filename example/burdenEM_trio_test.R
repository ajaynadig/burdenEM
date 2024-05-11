gnomad_information <- data.frame(data.table::fread("~/Mirror/robinson_lab/22q/PRS_expression_analyses/data/gnomad.v2.1.1.lof_metrics.by_gene.txt"))

#Load in Fu 2022 SSC/ASC results, to get gene IDs
fu_sscasc <- read.csv("~/Mirror/oconnor_rotation/rare_dn_h2/fu_files/fu_asc_ssc_counts.csv")

#combined_dataset will have counts and mutation rates combined between Fu 2022 cohorts
combined_dataset <- data.frame(gene = fu_sscasc$gene,
                               gene_id = fu_sscasc$gene_id,
                               N = 8028 + 7008)
rownames(combined_dataset) <- combined_dataset$gene_id

#Load in variants from Fu et al, remove indels
fu_variants <- read.csv("~/Mirror/oconnor_rotation/rare_dn_h2/fu_files/fu_variants.csv")
fu_variants_nonindel <- fu_variants[!fu_variants$isIndel,]

#Add PTV info
combined_dataset[,"case_count"] =   sapply(combined_dataset$gene_id,
                                                         function(x) {
                                                           sum(fu_variants_nonindel$Role == "Proband" &
                                                                 fu_variants_nonindel$isPTV &
                                                                 fu_variants_nonindel$gene_id == x,
                                                               na.rm = TRUE)
                                                         })

combined_dataset[,"case_rate"] = gnomad_information$mu_lof[match(combined_dataset$gene_id,gnomad_information$gene_id)]
combined_dataset = combined_dataset[!is.na(combined_dataset$case_count) & !is.na(combined_dataset$case_rate),]

#Get features
gnomad_oe_matched = gnomad_information$oe_lof[match(combined_dataset$gene_id,gnomad_information$gene_id)]
quantiles <- quantile(gnomad_oe_matched, probs = seq(0, 1, by = 0.2))

# Use cut to map each value to its quantile range
quantile_assignments <- cut(gnomad_oe_matched, breaks = quantiles, include.lowest = TRUE,
                            labels = c(1:5))


annot_oe5 <- sapply(1:5,
                    function(x) {
                      as.numeric(quantile_assignments == x)
                    })

rownames(annot_oe5) <- gnomad_information$gene_id[match(combined_dataset$gene_id,gnomad_information$gene_id)]


#Create input data

#Run model
coefs_est = burdenEM_trio(input_data = combined_dataset,
                          features = annot_oe5,
                          prevalence = 0.01)


