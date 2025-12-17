
number_of_folds=5
possible_corr_ev_models = c('ARD', 'ER')

possible_rate_cats = c(1,2)

repo_path = Sys.getenv('KEWSCRATCHPATH')
source(file.path(repo_path, 'phyloKNN_analysis', 'data','helpful_phyl_methods.R'))
prediction_path = file.path(repo_path, 'phyloKNN_analysis', 'imputation')


get_prediction_data_paths <- function(case, ev_model, iteration, missingness_type) {
  return(file.path(prediction_path, case, ev_model, as.character(iteration), missingness_type))
}

# Methods to format PI outputs to standarise like phyestimatedisc output
format_corhmm <- function(corhmm_output, plant_names_to_predict, ratecat){
  # print('Formatting')
  # Extract tip states from corhmm output
  # This could be much cleaner. See https://github.com/thej022214/corHMM/issues/62
  output_data = data.frame(corhmm_output$tip.states)
  
  # Filter data for plants to predict
  plants_to_predict_that_are_in_tree = intersect(plant_names_to_predict, rownames(output_data))
  output_data = output_data[plants_to_predict_that_are_in_tree,]
  
  
  if (ratecat ==1){
    output_data['0']=output_data$X1
    output_data['1']= output_data$X2
  }
  if (ratecat ==2){
    output_data['0']=output_data$X1+output_data$X3
    output_data['1']= output_data$X2+output_data$X4
  }
  
  output_data['estimated.state'] = lapply(output_data['1'], round)
  
  # make index into column
  output_data <- cbind(accepted_species = rownames(output_data), output_data)
  rownames(output_data) <- 1:nrow(output_data)
  return(output_data)
}

get_folds <- function(non_missing_data){
  if(length(unique(non_missing_data[[target]])) == 1) {
    skfolds =NULL
    unique_target = unique(non_missing_data[[target]])[1]
    print(paste('Unique target found for', iteration, unique_target, sep=':'))
  } else {
    # This does unstratified sampling
    skfolds = caret::createFolds(non_missing_data[[target]], k=number_of_folds)
    unique_target = NULL
  }
  return(skfolds)
}

calculate_brier <- function(y_score,y_true){
  return(mean((y_score - y_true)^2))
}

calculate_inverse_AP <- function(y_score,y_true){
    library(reticulate)
  use_virtualenv("~/Documents/virtual_environments/scratch_interpreter")
  sklearn <- import("sklearn")
  score = -sklearn$metrics$average_precision_score(y_true, y_score)
  print('AP:')
  print(score)
  return(score)
}

run_corHMM_models <- function(case, simulation_ev_model, iteration, missing_type, bin_or_cont, scorer){
    
    outpath = get_prediction_data_paths(case, simulation_ev_model, iteration, missing_type)
    dir.create(outpath, recursive=TRUE)
    out_file = file.path(outpath, 'corHMM.csv')
    if (file.exists(out_file)){
      print(paste('file already exists for corhmm', case, simulation_ev_model, iteration, missing_type, sep=':'))
    }else{
      print(paste('running corhmm', case, simulation_ev_model, iteration, missing_type, sep=':'))
      if (missing(scorer)) {
        scorer = 'brier'
      }
      if (!(scorer %in% c('brier', 'AP'))) {
        stop('Incorrectly specified scorer')
      }
      setup_ = set_up(case, simulation_ev_model, iteration, missing_type)
      labelled_tree = setup_$labelled_tree
      missing_values_with_tree_labels = setup_$missing_values_with_tree_labels
      target = setup_$target
      non_missing_data = setup_$non_missing_data
      skfolds = get_folds(non_missing_data)
      training_tree = setup_$training_tree
      unique_target = setup_$unique_target
      
      if (is.null(unique_target)){
        best_brier_score = 1
        best_rate_cat = possible_rate_cats[1]
        best_ev_model = possible_corr_ev_models[1]
        for (rate_cat in possible_rate_cats) {
          for (ev_model in possible_corr_ev_models) {
            
            brier_score_for_this_config = 0
            number_of_successful_folds = 0
            for (i in 1:number_of_folds) {
              fold_indices <- skfolds[[i]]
              kfold_test_plants = non_missing_data[fold_indices,]$accepted_species
              
              
              #### corHMM
              test_data_with_tree_labels = data.frame(non_missing_data)
              # set unknown using '?' for corhmm
              # corhmm will estimate trait values for all tips in tree that are in trait data with ? value
              test_data_with_tree_labels[[target]][test_data_with_tree_labels$accepted_species %in% kfold_test_plants] = '?'
              cor_trait_data = data.frame(test_data_with_tree_labels)
              cor_trait_data = subset(cor_trait_data, select = c("accepted_species", target))
              ## Catch errors, this can happen where split means target is uniform value
              try(
                {
                  corHMM_predicted_values = corHMM::corHMM(training_tree, cor_trait_data,model=ev_model,
                                                           rate.cat = rate_cat, get.tip.states = TRUE, n.cores = 10)
                  
                  out = format_corhmm(corHMM_predicted_values, kfold_test_plants, rate_cat)
                  
                  validation_data = non_missing_data[non_missing_data$accepted_species %in% kfold_test_plants,]
                  
                  df_merge <- merge(out,validation_data,by="accepted_species") 
                  
                  f_t = df_merge$`1`
                  o_t = as.numeric(df_merge[[target]])
                  if (scorer == 'brier') {
                    brier_score_for_this_fold = calculate_brier(f_t,o_t)
                  } else if( scorer == 'AP'){
                    brier_score_for_this_fold = calculate_inverse_AP(f_t,o_t)
                  }
                  
                  if(!is.na(brier_score_for_this_fold)){
                    brier_score_for_this_config = brier_score_for_this_config + brier_score_for_this_fold
                    number_of_successful_folds = number_of_successful_folds+1
                  }
                  
                }, silent = TRUE)
            }
            if(number_of_successful_folds!=0){
              brier_score_for_this_config = brier_score_for_this_config/number_of_successful_folds
              if (brier_score_for_this_config<best_brier_score){
                best_brier_score=brier_score_for_this_config
                best_rate_cat = rate_cat
                best_ev_model = ev_model
              }
            }
            
          }
        }
        
        final_test_data_with_tree_labels = data.frame(missing_values_with_tree_labels)
        plants_to_predict = final_test_data_with_tree_labels[is.na(final_test_data_with_tree_labels[[target]]),]$accepted_species
        # set unknown using '?' for corhmm
        # corhmm will estimate trait values for all tips in tree that are in trait data with ? value
        final_test_data_with_tree_labels[[target]][is.na(final_test_data_with_tree_labels[[target]])] = '?'
        final_cor_trait_data = data.frame(final_test_data_with_tree_labels)
        final_cor_trait_data = subset(final_cor_trait_data, select = c("accepted_species", target))
        
        final_corHMM_predicted_values = corHMM::corHMM(labelled_tree, final_cor_trait_data,model=best_ev_model,
                                                       rate.cat = best_rate_cat, get.tip.states = TRUE, n.cores = 10)
        
        
        final_out = format_corhmm(final_corHMM_predicted_values, plants_to_predict, best_rate_cat)
        
        final_out = subset(final_out, select = c("accepted_species", '0','1'))
      }else{
        # Where the target value is unique in the training data, just assign that value to predictions
        final_test_data_with_tree_labels = data.frame(missing_values_with_tree_labels)
        final_out = final_test_data_with_tree_labels[is.na(final_test_data_with_tree_labels[[target]]),]
        
        if(unique_target==1){
          final_out[['0']] <- 0
          final_out[['1']] <- 1
          
        }else if(unique_target==0){
          final_out[['0']] <- 1
          final_out[['1']] <- 0
          
        }else {
          stop(paste("Unknown unique target for binary case", bin_or_cont))
        }
        
        final_out = subset(final_out, select = c("accepted_species", '0','1'))
        best_ev_model =NULL
        best_rate_cat =NULL
      }
      
      
      
      write.csv(final_out, out_file, row.names = FALSE)
      param_df = data.frame(best_ev_model=c(best_ev_model), best_rate_cat=c(best_rate_cat))
      write.csv(param_df, file.path(outpath, 'corHMM_hparams.csv'), row.names = FALSE)
    }
    
    
}

run_picante_instance<- function(bin_or_cont,trait_data,target,plants_to_predict,tree, ev_model, method){
  picante_data = data.frame(trait_data)
  picante_data = subset(picante_data, select = c("accepted_species", target))
  
  # print('Drop test samples')
  picante_train_data = subset(picante_data, !(accepted_species %in% plants_to_predict))
  
  named_vector <- as.character(picante_train_data[[target]])
  names(named_vector) <- picante_train_data$accepted_species
  
  ## Set best.state to false following https://github.com/skembel/picante/issues/33

  # See https://rdrr.io/cran/ape/man/ace.html for model parameters
  ### Following molinavenegas_how_2024 check the assign random state to begin with and rerun if there is warnings
  ## Specifically following appendix 1
  if(bin_or_cont=='binary'){
    should_break =FALSE
    max_iter =100
    for (k in 1:max_iter) {
      # print(paste("Iteration", k))  # Check if loop is progressing
      # print(should_break)
      if (should_break){break}
      
      ip=runif(n = 1, min = .01, max = .99)
      tryCatch({
        # print("Calling phyEstimateDisc...")  # Debugging message to confirm function call
        phyEstimate_predicted_values = picante::phyEstimateDisc(phy=tree, method=method,
                                                                model = ev_model,
                                                                trait = named_vector,
                                                                best.state=FALSE,
                                                                ip = ip)
        # # If no error or warning, break the loop after successful assignment
        # print("phyEstimateDisc executed successfully.")  # Debugging message
        should_break <- TRUE
        # print(paste("Successfully stored value in iteration", k))
      },
      warning = function(w){
        # print('not breaking')
        # print(w)
        should_break=FALSE
        
      })
      if(k==max_iter){
        print('Max iter reached for picante')
        phyEstimate_predicted_values = picante::phyEstimateDisc(phy=tree, 
                                                                model = ev_model,
                                                                trait = named_vector,
                                                                best.state=FALSE,
                                                                ip = ip)
      }
      
    }
  } else if(bin_or_cont=='continuous'){
    phyEstimate_predicted_values = picante::phyEstimate(phy=tree, trait = named_vector,
                                            method = method, model=ev_model)
  }
  
  # make index into column
  output_data <- cbind(accepted_species = rownames(phyEstimate_predicted_values), phyEstimate_predicted_values)
  rownames(output_data) <- 1:nrow(output_data)
  return(output_data)
}

run_picante_models <- function(case, simulation_ev_model, iteration, missing_type, bin_or_cont, scorer){
  
  outpath = get_prediction_data_paths(case, simulation_ev_model, iteration, missing_type)
  dir.create(outpath, recursive=TRUE)
  out_file = file.path(outpath, 'picante.csv')
  if (file.exists(out_file)){
    print(paste('file already exists for picante', case, simulation_ev_model, iteration, missing_type, sep=':'))
  } else{
    print(paste('running picante', case, simulation_ev_model, iteration, missing_type, sep=':'))
    if (missing(scorer)) {
      scorer = 'brier'
    }
    
    
    setup_ = set_up(case, simulation_ev_model, iteration, missing_type)
    labelled_tree = setup_$labelled_tree
    if(!ape::is.binary(labelled_tree)){
      labelled_tree=ape::multi2di(labelled_tree)
    }
    missing_values_with_tree_labels = setup_$missing_values_with_tree_labels
    target = setup_$target
    non_missing_data = setup_$non_missing_data
    skfolds = get_folds(non_missing_data)
    training_tree = setup_$training_tree
    if(!ape::is.binary(training_tree)){
      training_tree=ape::multi2di(training_tree)
    }
    if (bin_or_cont == "binary") {
      possible_picante_ev_models = c('ARD', 'ER')# 'ER' and SYM same for binary traits
      possible_picante_methods = c('ML') ## Only ML is available for discrete characters. See https://rdrr.io/cran/ape/man/ace.html
    } else if  (bin_or_cont == "continuous"){
      possible_picante_ev_models = c('BM')#, 'mvOU','OU', "lambda", "kappa", "delta", "EB", "star") #check https://github.com/emmanuelparadis/ape/issues/139
      possible_picante_methods = c('REML','ML', 'pic')
    }
    
    unique_target = setup_$unique_target
    if (is.null(unique_target)){
      best_score = NULL
      best_ev_model = possible_picante_ev_models[1]
      best_method = possible_picante_methods[1]
      for (ev_model in possible_picante_ev_models) {
        for (method in possible_picante_methods){
          
          if (!(method=='pic' && ev_model != "BM")){ #the "pic" method can be used only with model = "BM"
            score_for_this_config = 0
            number_of_successful_folds = 0
            for (i in 1:number_of_folds) {
              fold_indices <- skfolds[[i]]
              kfold_test_plants = non_missing_data[fold_indices,]$accepted_species
              
              
              #### picante
              try({
                output_data = run_picante_instance(bin_or_cont, non_missing_data,target,kfold_test_plants,training_tree,ev_model, method)
                
                validation_data = non_missing_data[non_missing_data$accepted_species %in% kfold_test_plants,]
                
                df_merge <- merge(output_data,validation_data,by="accepted_species") 
                if (bin_or_cont == "binary") {
                  f_t = as.numeric(df_merge$`1`)
                  o_t = as.numeric(df_merge[[target]])
                  if (scorer == 'brier') {
                    brier_score_for_this_fold = calculate_brier(f_t,o_t)
                  } else if( scorer == 'AP'){
                    brier_score_for_this_fold = calculate_inverse_AP(f_t,o_t)
                  } else{
                    stop('Incorrectly specified scorer')
                  }
                  if(!is.na(brier_score_for_this_fold)){
                    score_for_this_config = score_for_this_config + brier_score_for_this_fold
                    number_of_successful_folds = number_of_successful_folds+1
                  }} else if (bin_or_cont == "continuous"){
                    mae_for_this_fold = Metrics::mae(df_merge[[target]], df_merge$estimate)
                    if(!is.na(mae_for_this_fold)){
                      score_for_this_config = score_for_this_config + mae_for_this_fold
                      number_of_successful_folds = number_of_successful_folds+1
                    }
                  }
              }, silent = TRUE)
              
            }
            if(number_of_successful_folds!=0){
              score_for_this_config = score_for_this_config/number_of_successful_folds
              if (is.null(best_score) || score_for_this_config<best_score){
                best_score=score_for_this_config
                best_ev_model = ev_model
                best_method = method
              }
            }
          }
          
          
        }
        
        
      }
      
      plants_to_predict = missing_values_with_tree_labels[is.na(missing_values_with_tree_labels[[target]]),]$accepted_species
      
      final_out = run_picante_instance(bin_or_cont,missing_values_with_tree_labels,target,plants_to_predict,labelled_tree,
                                       best_ev_model, best_method)
      
      # Reorder matrix based on the species_order
      final_out <- final_out[match(plants_to_predict, final_out[, "accepted_species"]), ]
      
      if (bin_or_cont == "binary") {
        final_out = subset(final_out, select = c("accepted_species", '0','1'))
      } else if (bin_or_cont == "continuous") {
        final_out = subset(final_out, select = c("accepted_species", 'estimate'))
      }
    } else{
      # Where the target value is unique in the training data, just assign that value to predictions
      if (bin_or_cont == "binary") {
        final_test_data_with_tree_labels = data.frame(missing_values_with_tree_labels)
        final_out = final_test_data_with_tree_labels[is.na(final_test_data_with_tree_labels[[target]]),]
        
        if(unique_target==1){
          final_out[['0']] <- 0
          final_out[['1']] <- 1
          
        }else if(unique_target==0){
          final_out[['0']] <- 1
          final_out[['1']] <- 0
          
        }else {
          stop(paste("Unknown unique target for binary case", bin_or_cont))
        }
        
        final_out = subset(final_out, select = c("accepted_species", '0','1'))
        best_ev_model = NULL
        best_method = NULL
      } else if (bin_or_cont == "continuous") {
        stop(paste("Found unique target for continuous case", bin_or_cont))
      }
    }
    
    write.csv(final_out, out_file, row.names = FALSE)
    param_df = data.frame(best_ev_model=c(best_ev_model), best_method=c(best_method))
    write.csv(param_df, file.path(outpath, 'picante_hparams.csv'), row.names = FALSE)
  }
  
  
  
}


