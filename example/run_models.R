path_to_repo = "~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/"
#Source files--unneccessary if package is installed
setwd("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/burdenEM_trio.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/EM.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/estimate_heritability.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/io.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/likelihoods.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/model.R")
source("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/secondary_analysis_functions.R")
Rcpp::sourceCpp("~/Mirror/oconnor_rotation/rare_dn_h2/github/burdenEM/R/EM.cpp")
#set up data structures
source(paste0(path_to_repo,"example/set_up_asc.R"))
source(paste0(path_to_repo,"example/set_up_kaplanis.R"))



# Get the current date in "MonYY" format (e.g., "Mar25")
current_date <- format(Sys.Date(), "%b%d_%y")

output_path = "/Users/anadig/Mirror/oconnor_rotation/rare_dn_h2/github/outputs/"

run_autism = TRUE
run_ddd = TRUE

#Make the bootstrap resamples to keep consistent across models
if (run_autism) {
  autism_ngene = nrow(autism_data$counts)

  bootstrap_samples_autism <- sapply(1:100,
                                     function(dummy) {
                                       sample(1:autism_ngene,replace = TRUE)
                                     })

  burdenEM_models_autism <- lapply(autism_data$prev_factors,
                                   function(prev_factor) {
                                     print("Prevalence Scaling Factor")
                                     print(prev_factor)
                                     lapply(1:length(autism_data$loop_vars$datasets),
                                            function(i) {

                                              processed_input = get_genetic_data(i,autism_data)

                                              input_df = processed_input$genetic_data
                                              features = processed_input$features
                                              print(colSums(features))

                                              model <- burdenEM_trio(input_data = input_df,
                                                                     component_endpoints = seq(0,log(1/0.01),length.out = 10),
                                                                     features = features,
                                                                     null_sim = TRUE,
                                                                   #  n_null = 5,
                                                                     #max_iter = 5,
                                                                     prevalence = autism_data$loop_vars$prevalences[i] * prev_factor,
                                                                     estimate_posteriors = TRUE,
                                                                     bootstrap_samples = bootstrap_samples_autism)
                                              # print(model$delta)
                                              print(model$mutvar_output$total_h2)
                                              print(model$mutvar_output$enrichment)
                                              print(model$mutvar_output$mutvar_CI)
                                              #print(model$mutvar_output$total_h2_p)
                                              return(model)
                                            })
                                   })

  save(burdenEM_models_autism,
       file = paste0(output_path,"models_autism_", current_date, ".Rdata"))
}



if (run_ddd){
  ddd_ngene = nrow(kaplanis_data$counts)

  bootstrap_samples_ddd <- sapply(1:100,
                                  function(dummy) {
                                    sample(1:ddd_ngene,replace = TRUE)
                                  })

  burdenEM_models_DDD <-lapply(kaplanis_data$prev_factors,
                               function(prev_factor) {
                                 print("Prevalence Scaling Factor")
                                 print(prev_factor)
                                 lapply(1:length(kaplanis_data$loop_vars$subsets),
                                        function(i) {

                                          processed_input = get_genetic_data(i,kaplanis_data)

                                          input_df = processed_input$genetic_data
                                          features = processed_input$features
                                          print(colSums(features))
                                          model <- burdenEM_trio(input_data = input_df,
                                                                 component_endpoints = seq(0,log(1/0.01),length.out = 10),
                                                                 features = features,
                                                                 null_sim = FALSE,
                                                                 #num_iter = 5,
                                                                 prevalence = kaplanis_data$loop_vars$prevalences[i] * prev_factor,
                                                                 estimate_posteriors = TRUE,
                                                                 bootstrap_samples = bootstrap_samples_ddd)
                                          #print(model$delta)
                                          print(model$mutvar_output$total_h2)
                                          print(model$mutvar_output$enrichment)
                                          print(model$mutvar_output$mutvar_CI)
                                          # print(model$mutvar_output$total_h2_p)
                                          return(model)
                                        })
                               })
  save(burdenEM_models_DDD,
       file = paste0(output_path,"models_ddd_", current_date, ".Rdata"))


}









