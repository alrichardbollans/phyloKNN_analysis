library(doParallel)
library(foreach)

num_cores <- 12#set to 12 as this function seems to find all threads on cluster node detectCores()  # Use all cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

case = 'ultrametric'
simulation_ev_model= 'APM'
bin_or_cont='binary'
number_of_simulation_iterations = 10

foreach(iter = 1:number_of_simulation_iterations) %dopar% {
  repo_path = Sys.getenv('KEWSCRATCHPATH')
  source(file.path(repo_path, 'phyloKNN_analysis', 'imputation','R_binary_imputation_helper_functions.R'))
  missingness_types = c('mcar')
  
  print(iteration)
  for (missing_type in missingness_types) { 
    run_picante_models(case, simulation_ev_model, iteration, missing_type, bin_or_cont, 'AP')
    run_corHMM_models(case, simulation_ev_model, iteration, missing_type, bin_or_cont, 'AP')
    
  }
}
