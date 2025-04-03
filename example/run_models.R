source("example/set_up_asc.R")
source("example/set_up_kaplanis.R")


# Get the current date in "MonYY" format (e.g., "Mar25")
current_date <- format(Sys.Date(), "%b%y")

output_path = "/Users/anadig/Mirror/oconnor_rotation/rare_dn_h2/repo/burdenEM/outputs/"

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
                                                                     null_sim = FALSE,
                                                                     #max_iter = 5,
                                                                     prevalence = autism_data$loop_vars$prevalences[i] * prev_factor,
                                                                     estimate_posteriors = TRUE,
                                                                     bootstrap_samples = bootstrap_samples_autism)
                                              # print(model$delta)
                                              print(model$heritability_output$total_h2)
                                              print(model$heritability_output$enrichment)
                                              print(model$heritability_output$heritability_CI)
                                              #print(model$heritability_output$total_h2_p)
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
                                          print(model$heritability_output$total_h2)
                                          print(model$heritability_output$enrichment)
                                          print(model$heritability_output$heritability_CI)
                                          # print(model$heritability_output$total_h2_p)
                                          return(model)
                                        })
                               })
  save(burdenEM_models_DDD,
       file = paste0(output_path,"models_ddd_", current_date, ".Rdata"))


}









