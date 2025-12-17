
repo_path = Sys.getenv('KEWSCRATCHPATH')
source(file.path(repo_path, 'phyloKNN_analysis', 'imputation','R_binary_imputation_helper_functions.R'))
possible_phylopars_models = c('BM', 'mvOU','OU', "lambda", "kappa", "delta", "EB", "star") # provided here https://cran.r-project.org/web/packages/Rphylopars/Rphylopars.pdf

format_phylopars <- function(phylopars_predictions,kfold_test_plants, target){
  plants_to_predict_that_are_in_tree = intersect(kfold_test_plants, rownames(phylopars_predictions))
  phylopars_predictions = phylopars_predictions[plants_to_predict_that_are_in_tree,target]
  output_data = data.frame(phylopars_predictions)
  
  
  output_data['estimate'] = output_data['phylopars_predictions']
  output_data = output_data[c('estimate')]
  # make index into column
  output_data <- cbind(accepted_species = rownames(output_data), output_data)
  rownames(output_data) <- 1:nrow(output_data)
  return(output_data)
}

run_phylopars_models <- function(case, simulation_ev_model, iteration, missing_type, bin_or_cont){
  out_path = get_prediction_data_paths(case, simulation_ev_model, iteration, missing_type)
  dir.create(out_path, recursive=TRUE)
  out_file = file.path(out_path, 'phylopars.csv')
  if (file.exists(out_file)){
    print(paste('file already exists for phylopars', case, simulation_ev_model, iteration, missing_type, sep=':'))
  }else{
    print(paste('running phylopars', case, simulation_ev_model, iteration, missing_type, sep=':'))
    setup_ = set_up(case, simulation_ev_model, iteration, missing_type)
    labelled_tree = setup_$labelled_tree
    # if(!ape::is.ultrametric(labelled_tree)){
    #   labelled_tree = phytools::force.ultrametric(labelled_tree) # phylopars needs an ultrametric tree but only for OU model
    # }
    missing_values_with_tree_labels = setup_$missing_values_with_tree_labels
    target = setup_$target
    non_missing_data = setup_$non_missing_data
    skfolds = get_folds(non_missing_data)
    training_tree = setup_$training_tree
    # if(!ape::is.ultrametric(training_tree)){
    #   training_tree = phytools::force.ultrametric(training_tree) # phylopars needs an ultrametric tree but only for OU model
    # }
    
    
    best_mae = -1
    best_ev_model = possible_phylopars_models[1]
    for (ev_model in possible_phylopars_models) {
      if((ev_model == 'OU' & case == 'with_extinct') | (ev_model == 'mvOU' & case == 'with_extinct')){
        print('OU model not currently supported for non-ultrametric trees.')
      }else{
        
      
      mae_for_this_config = 0
      number_of_successful_folds = 0
      for (i in 1:number_of_folds) {
        
        fold_indices <- skfolds[[i]]
        
        kfold_test_plants = non_missing_data[fold_indices,]$accepted_species
        
        test_data_with_tree_labels = data.frame(non_missing_data)
        # set some to unknown using NA
        test_data_with_tree_labels[[target]][test_data_with_tree_labels$label %in% kfold_test_plants] = NA
        
        phylopars_data = data.frame(test_data_with_tree_labels)
        colnames(phylopars_data)[1]  <- "species" #First column name of trait_data MUST be 'species' (all lower case).
        
        phylopars_data = subset(phylopars_data, select = c("species", target))
        ## Catch errors, this can happen for certain ev models (think just kappa)
        ### and OU model not currently supported for non-ultrametric trees.
        try(
          {
            p_v = Rphylopars::phylopars(phylopars_data, training_tree, model = ev_model)
            
            phylopars_predictions = p_v$anc_recon
            
            out = format_phylopars(phylopars_predictions, kfold_test_plants,target)
            validation_data = non_missing_data[non_missing_data$accepted_species %in% kfold_test_plants,]
            
            df_merge <- merge(out,validation_data,by="accepted_species")
            mae_for_this_fold = Metrics::mae(df_merge[[target]], df_merge$estimate)
            mae_for_this_config = mae_for_this_config+mae_for_this_fold
            number_of_successful_folds = number_of_successful_folds+1
          }, silent = FALSE
        )
        
      }
      if(number_of_successful_folds!=0){
        mae_for_this_config = mae_for_this_config/number_of_successful_folds
        if (mae_for_this_config<best_mae || best_mae==-1){
          best_mae=mae_for_this_config
          best_ev_model = ev_model
        }
      }
      }
    }
    
    # Now use best model
    final_test_data_with_tree_labels = data.frame(missing_values_with_tree_labels)
    
    final_phylopars_data = data.frame(final_test_data_with_tree_labels)
    colnames(final_phylopars_data)[1]  <- "species" #First column name of trait_data MUST be 'species' (all lower case).

    final_phylopars_data = subset(final_phylopars_data, select = c("species", target))
    final_p_v = Rphylopars::phylopars(final_phylopars_data, labelled_tree, model = best_ev_model)
    final_phylopars_predictions = final_p_v$anc_recon
    
    plants_to_predict = final_phylopars_data[is.na(final_phylopars_data[[target]]),]$species
    
    final_out = format_phylopars(final_phylopars_predictions, plants_to_predict,target)
    
    
    write.csv(final_out, out_file, row.names = FALSE)
    
    param_df = data.frame(best_ev_model=c(best_ev_model))
    write.csv(param_df, file.path(out_path, 'phylopars_hparams.csv'), row.names = FALSE)
  }
  
  
  
}
