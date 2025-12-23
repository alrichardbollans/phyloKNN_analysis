do_broken_stick <- function(sim_folder, eigenvalues, missing_type){
  eigenValue.df <- data.frame(eigenvalues)
  
  # As in Bachman, Steven Philip, et al. "Extinction risk predictions for the world's flowering plants to support their conservation." bioRxiv (2023): 2023-08.
  # use the Broken Stick method to select the eigenvectors
  # Discussed in Diniz-Filho, J. A. F., Sant'Ana, C. E. R. d., & Bini, L. M. (1998). An Eigenvector Method for Estimating Phylogenetic Inertia. Evolution, 52(5), 1247-1262. https://doi.org/10.2307/2411294
  # Original citation Jackson, Donald A. "Stopping rules in principal components analysis: a comparison of heuristical and statistical approaches." Ecology 74.8 (1993): 2204-2214.
  # The code calculates Broken Stick values cumulatively in reverse order (e.g., for j=1, it sums the last eigenvalue's expected variance). Sorting these values in decreasing order coincidentally aligns them with the correct order.

  number_eigenvalues <- length(eigenvalues)
  broken_stick_model_df <- data.frame(j=seq(1:number_eigenvalues), p=0)
  # The first value of p is set to 1/n, which is the expected proportion of variance explained by the first eigenvector under the broken stick model.
  broken_stick_model_df$p[1] <- 1/number_eigenvalues 
  #The subsequent values of p are calculated iteratively by adding 1/(n+1-i) to the previous value, following the broken stick distribution.
  for (i in 2:number_eigenvalues)
  {
    broken_stick_model_df$p[i] = broken_stick_model_df$p[i-1] + (1/(number_eigenvalues + 1 - i))
  }
  #Finally, the p values are scaled to percentages by multiplying by 100 and dividing by n.
  broken_stick_model_df$p <- 100*broken_stick_model_df$p/number_eigenvalues
  
  # get a data frame with selected vectors
  eigenValue.df$percent.expected <- 100*eigenvalues/sum(eigenvalues) # Percent expected explained variance by the eigenvector according to eigenvalue
  eigenValue.df$broken.stick <- sort(broken_stick_model_df$p, decreasing = T) # Expected explained according to p value
  eigenValue.df$diff <- eigenValue.df$percent.expected - eigenValue.df$broken.stick
  selected_eigenValue.df <- eigenValue.df[eigenValue.df$diff > 0,] # Use those vectors where the variance explained according to eigenvalue is > according to p value
  num_selected = nrow(selected_eigenValue.df) # first X vectors selected
  
  explained_var= sum(selected_eigenValue.df$percent.expected) # X% of the variance
  
  param_df = data.frame(broken_stick_number=c(num_selected),explained_var=c(explained_var))
  
  if (missing(missing_type)){
    write.csv(param_df,file.path(sim_folder,"broken_stick_parameters.csv"))
  } else{
    write.csv(param_df,file.path(sim_folder,paste(missing_type,'PEM_broken_stick_parameters.csv', sep='_')))
  }
  
  return(param_df)
}

decompose_tree <- function(sim_folder){
  ground_tree = ape::read.tree(file.path(sim_folder, 'tree.tre'))
  pvrd <- PVR::PVRdecomp(ground_tree)
  # get the vectors and their names into a data frame
  eigenVec.df <- data.frame(pvrd@Eigen$vectors,row.names=pvrd@phylo$tip.label)
  write.csv(eigenVec.df, file = file.path(sim_folder,'all_eigenvectors.csv'))
  do_broken_stick(sim_folder, pvrd@Eigen$values)
}

get_PEMS <- function(setup_, output_folder, missing_type){
  
  labelled_tree = setup_$labelled_tree
  non_missing_data = setup_$non_missing_data
  target = setup_$target
  
  # Follows package tutorial https://cran.r-project.org/web/packages/MPSEM/vignettes/REM_with_MPSEM.html#cross-validating-pem-predictions
  if(length(unique(non_missing_data[[target]])) == 1){
    
    # In this case L-BFGS-B breaks, so use defaults
    MPSEM::Phylo2DirectedGraph(
      labelled_tree
    ) -> all.pgraph
    
    MPSEM::PEM.build(all.pgraph,
      d = "distance",
      sp="species"
    ) -> trained.PEM_opt
    as.data.frame(trained.PEM_opt) -> all_df
  } else{
    # Estimate weighting parameters empirically
    # print('Estimate weighting parameters empirically')
    
    training_tree = setup_$training_tree
    
    
    MPSEM::Phylo2DirectedGraph(
      training_tree
    ) -> train.pgraph
    
    MPSEM::PEM.fitSimple(
      y = non_missing_data[,target],
      x = NULL,
      w = train.pgraph,
      d = "distance",
      sp="species",
      lower = 0,
      upper = 1
    ) -> trained.PEM_opt
    as.data.frame(trained.PEM_opt) -> trained.PEM_df
    
    # Get values for test tips
    missing_values_with_tree_labels = setup_$missing_values_with_tree_labels
    nan_tips = missing_values_with_tree_labels[is.na(missing_values_with_tree_labels[target]), ]$accepted_species
    MPSEM::getGraphLocations(
      labelled_tree,
      targets = nan_tips
    ) -> targets.loc
    test_PEM_df = data.frame(MPSEM::Locations2PEMscores(trained.PEM_opt, targets.loc)$scores)
    all_df = rbind(trained.PEM_df,test_PEM_df)
  }
  
  
  write.csv(all_df, file = file.path(output_folder,paste(missing_type,'all_PEMS.csv', sep='_')))
  # Looking at the source code, it makes sense the d values are the eigenvalues
  # https://github.com/guenardg/MPSEM_dev/blob/f8b4a384001808951814f5bd38636cf1427bce92/Old/PEM-functions.R#L359C3-L359C6
  bs_values = do_broken_stick(output_folder,trained.PEM_opt$d, missing_type)
  
  
  ## In e.g. https://pmc.ncbi.nlm.nih.gov/articles/PMC4612606/#B11 a forward stepwise procedure is usd to pick eigenvectors
  ## and this is a utilitiy in the package for linear models:
  # MPSEM::lmforwardsequentialAICc(
  #   y = non_missing_data[,target],
  #   object = trained.PEM_opt
  # ) -> lm2
}

