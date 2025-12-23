repo_path <- Sys.getenv("KEWSCRATCHPATH")   
source('helpful_phyl_methods.R')
source('helper_simulation_methods.R')
source('add_eigenvectors.R')

for(iteration in 1:number_of_repetitions){
  missingness_types = c('mcar', 'phyloNa')
  binary_ev_models = c('ER', 'ARD', 'BiSSE', 'HiSSE', 'bBMT', 'Clonality')
  continuous_ev_models = c('BM', 'OU', 'EB', 'LB', 'BMT', 'Seed Mass')
  print(iteration)
  for (missing_type in missingness_types) {
    for(simulation_ev_model in c(binary_ev_models, continuous_ev_models)){# Keep the inner loop sequential
      if(simulation_ev_model == 'Clonality' | simulation_ev_model == 'Seed Mass'){
        cases = c('ultrametric')
      }
      else{
        cases = c('ultrametric', 'with_extinct')
        
        
      }
      for(case in cases){
        setup_ = set_up(case, simulation_ev_model, iteration, missing_type)
        output_folder = get_input_data_paths(case, simulation_ev_model, iteration)$value_path
        get_PEMS(setup_, output_folder, missing_type)
      }
    }
  }
}

