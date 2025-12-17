library(dplyr)
repo_path = Sys.getenv('KEWSCRATCHPATH')
input_data_path = file.path(repo_path, 'phyloKNN_analysis', 'data')

replace_space_with_underscore_in_name<- function(given_name){
  new = gsub(" ", "_",given_name)
  return(new)
}

replace_underscore_with_space_in_name<- function(given_tree_name){
  new = gsub("_", " ",given_tree_name)
  return(new)
}

set_labels_on_tree_to_acc_name<- function(tree){
  tree$tip.label = as.character(lapply(tree$tip.label,replace_underscore_with_space_in_name))
  return(tree)
}

get_matching_labels <- function(tree,data){
  # Gets data which appears in tree and appends 'label' column
  # First match by accepted names
  accepted_label_matches <-
    data %>%
    mutate(label = accepted_species) %>%  # Create a copy of 'accepted_species' as 'label'
    filter(label %in% tree$tip.label)     # Keep only rows where 'label' is in tree$tip.label
  
  # matching_labels = accepted_label_matches$label
  # Then drop any NaNs just in case
  data_with_tree_labels_no_nan = tidyr::drop_na(accepted_label_matches,'label')
  
  return(data_with_tree_labels_no_nan)
}


subset_tree <- function(tree, node_list) {
  drop_list <- tree$tip.label[! tree$tip.label %in% node_list]
  
  return(ape::drop.tip(tree, drop_list))
}

get_subset_of_tree_from_names <- function(tree, names_to_include){
  
  lab_data = data.frame(accepted_species=names_to_include)
  
  return(get_subset_of_tree_from_data(lab_data, tree))
}

get_subset_of_tree_from_data <- function(data, tree){
  
  lab_data = data.frame(data)
  
  # print(lab_data)
  labels = get_matching_labels(tree,lab_data)$label
  
  # drop all tips we haven't found matches for
  f_tree <- subset_tree(tree, labels)
  
  return(f_tree)
}

get_iteration_path_from_base <- function(base, case, ev_model, iteration) {
  if (ev_model == "BiSSE" || ev_model == "HiSSE") {
    basepath <- file.path(base,'simulations', case, ev_model)
  } else if(ev_model %in% c('Clonality', 'Seed Mass')){
    basepath = file.path(base, 'real_data', case, ev_model)
    
  }else if(ev_model %in% c('APM')){
    basepath = file.path(base, 'my_apm_data', case, ev_model)
    
  }else {
    basepath = file.path(base, 'simulations', case, 'standard')
  }
  
  treepath = file.path(basepath, as.character(iteration))
  value_path = file.path(basepath, as.character(iteration), ev_model)
  
  return(list(treepath=treepath,value_path=value_path ))
}

get_input_data_paths <- function(case, ev_model, iteration) {
  return(get_iteration_path_from_base(input_data_path, case, ev_model, iteration))
}

set_up <- function(case, ev_model, iteration, missing_type){
  
  input_data = get_input_data_paths(case, ev_model, iteration)
  
  missing_values = read.csv(file.path(input_data$value_path, paste(missing_type,'_values.csv',sep='')))
  
  
  prepared_tree = ape::read.tree(file.path(input_data$treepath, 'tree.tre'))
  
  # prepared_tree = set_labels_on_tree_to_acc_name(tree)
  labelled_tree = get_subset_of_tree_from_data(missing_values,prepared_tree)
  missing_values_with_tree_labels=get_matching_labels(labelled_tree,missing_values)
  
  target = names(missing_values)[2]
  non_missing_data = missing_values_with_tree_labels[!is.na(missing_values_with_tree_labels[[target]]),]
  
  if(length(unique(non_missing_data[[target]])) == 1) {
    unique_target = unique(non_missing_data[[target]])[1]
    print(paste('Unique target found for', iteration, unique_target, sep=':'))
  } else {

    unique_target = NULL
  }
  
  training_tree = get_subset_of_tree_from_data(subset(non_missing_data, select = c("accepted_species")),labelled_tree)
  
  ground_truth = read.csv(file.path(input_data$value_path, 'ground_truth.csv'))
  
  return(list(labelled_tree=labelled_tree, missing_values_with_tree_labels=missing_values_with_tree_labels,
              target=target, non_missing_data=non_missing_data, training_tree=training_tree, unique_target=unique_target, ground_truth=ground_truth))
}