library(doParallel)
library(foreach)

num_cores <- 12#set to 32 as this function seems to find all threads on cluster node detectCores()  # Use all cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

number_of_simulation_iterations = 100
foreach(iter = 1:number_of_simulation_iterations) %dopar% {
# for(iter in 1:number_of_simulation_iterations){
  repo_path = Sys.getenv('KEWSCRATCHPATH')
  source(file.path(repo_path, 'phyloKNN_analysis', 'imputation','R_continuous_imputation_helper_functions.R'))
  
  missingness_types = c('mcar', 'phyloNa')
  cases = c('ultrametric', 'with_extinct')
  binary_ev_models = c('ER', 'ARD', 'BiSSE', 'HiSSE', 'bBMT', 'Clonality')
  continuous_ev_models = c('BM', 'OU', 'EB', 'LB', 'BMT', 'Seed Mass')
  print(iter)
  for (missing_type in missingness_types) {
    for(simulation_ev_model in binary_ev_models){# Keep the inner loop sequential
        if(simulation_ev_model == 'Clonality'){
            # run_picante_models('ultrametric', simulation_ev_model, iter, missing_type, 'binary')
            # run_corHMM_models('ultrametric', simulation_ev_model, iter, missing_type, 'binary')
        }
        else{
            for(case in cases){
                # run_picante_models(case, simulation_ev_model, iter, missing_type, 'binary')
                # run_corHMM_models(case, simulation_ev_model, iter, missing_type, 'binary')
             }

        }
    }

     for(simulation_ev_model in continuous_ev_models){# Keep the inner loop sequential
        if(simulation_ev_model == 'Seed Mass'){
            # run_phylopars_models('ultrametric', simulation_ev_model, iter, missing_type, 'continuous')
          run_PEM_models('ultrametric', simulation_ev_model, iter, missing_type, 'continuous')
            # run_picante_models('ultrametric', simulation_ev_model, iter, missing_type, 'continuous')
            
        }
        else{
            for(case in cases){
                # run_phylopars_models(case, simulation_ev_model, iter, missing_type, 'continuous')
              run_PEM_models(case, simulation_ev_model, iter, missing_type, 'continuous')
                # run_picante_models(case, simulation_ev_model, iter, missing_type, 'continuous')
             }

        }

    }
    
  }
}
