
setwd('~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM')

source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/burdenEM_trio.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/EM.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/estimate_heritability.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/io.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/likelihoods.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/model.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/secondary_analysis_functions.R")

source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/example/set_up_asc.R")

output_dir = "~/Mirror/oconnor_rotation/rare_dn_h2/github/outputs/"


load(paste0(output_dir, "models_autism_Apr27_25.Rdata"))
#Get the MLE scaling factors

binomial_LL_ptv_allSPARK = binomial_analysis(model = burdenEM_models_autism[[1]][[5]],
                                             genetic_data_total = get_genetic_data(5,autism_data)$genetic_data,
                                             genetic_data_subsample = get_genetic_data(34,autism_data)$genetic_data)
binomial_LL_ptv_allASC = binomial_analysis(model = burdenEM_models_autism[[1]][[5]],
                                           genetic_data_total = get_genetic_data(5,autism_data)$genetic_data,
                                           genetic_data_subsample = get_genetic_data(31,autism_data)$genetic_data)
binomial_LL_ptv_allGeneDx = binomial_analysis(model = burdenEM_models_autism[[1]][[5]],
                                              genetic_data_total = get_genetic_data(5,autism_data)$genetic_data,
                                              genetic_data_subsample = get_genetic_data(37,autism_data)$genetic_data)


#Forecasting
mean_and_sd <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  return(c(mean = mean_x, sd = sd_x))
}


num_iter = 100
N_new_case_range = seq(1000,38088,length.out = 20)
forecast_df <- data.frame()

for (N_new_case in N_new_case_range) {
  print(N_new_case)

  for (dummy in 1:num_iter){
    SPARK_output <- burdenEM_power_forecasting_mod(get_genetic_data(5,autism_data)$genetic_data,
                                                   burdenEM_models_autism[[1]][[5]],
                                                   N_new_case = N_new_case,
                                                   gamma_scaling_factor = binomial_LL_ptv_allSPARK$MLE)
    SPARK_output_table = SPARK_output$count_table
    SPARK_output_table$dataset = "SPARK"
    SPARK_output_table$N_new_case = N_new_case
    ASC_output <- burdenEM_power_forecasting_mod(get_genetic_data(5,autism_data)$genetic_data,
                                                 burdenEM_models_autism[[1]][[5]],
                                                 N_new_case = N_new_case,
                                                 gamma_scaling_factor = binomial_LL_ptv_allASC$MLE)
    ASC_output_table = ASC_output$count_table
    ASC_output_table$dataset = "ASC"
    ASC_output_table$N_new_case = N_new_case
    GeneDx_output <- burdenEM_power_forecasting_mod(get_genetic_data(5,autism_data)$genetic_data,
                                                    burdenEM_models_autism[[1]][[5]],
                                                    N_new_case = N_new_case,
                                                    gamma_scaling_factor = binomial_LL_ptv_allGeneDx$MLE)
    GeneDx_output_table = GeneDx_output$count_table
    GeneDx_output_table$dataset = "GeneDx"
    GeneDx_output_table$N_new_case = N_new_case

    NoScaling_output <- burdenEM_power_forecasting_mod(get_genetic_data(5,autism_data)$genetic_data,
                                                       burdenEM_models_autism[[1]][[5]],
                                                       N_new_case = N_new_case,
                                                       gamma_scaling_factor = 1)
    NoScaling_output_table = NoScaling_output$count_table
    NoScaling_output_table$dataset = "NoScaling"
    NoScaling_output_table$N_new_case = N_new_case

    iter_df <- rbind(SPARK_output_table,ASC_output_table,GeneDx_output_table)
    forecast_df = rbind(forecast_df, iter_df)
  }
}
forecast_df_average <- aggregate(
  cbind(count_fdr, count_bonferroni) ~ threshold + dataset + N_new_case,
  data = forecast_df,
  FUN = mean_and_sd
)

forecast_df_average_expand <- do.call(data.frame, forecast_df_average)
names(forecast_df_average_expand) <- c(
  "threshold", "dataset", "N_new_case",
  "count_fdr_mean", "count_fdr_sd",
  "count_bonferroni_mean", "count_bonferroni_sd"
)

write.table(forecast_df_average_expand,file = "~/Mirror/oconnor_rotation/rare_dn_h2/outputs/forecasting_output_avgs_May42025.tsv",
            quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
write.table(forecast_df,file = "~/Mirror/oconnor_rotation/rare_dn_h2/outputs/forecasting_output_full_May2025.tsv",
            quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")


