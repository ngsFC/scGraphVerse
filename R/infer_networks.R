infer_networks <- function(count_matrices_list, method = "GENIE3") {
  if (!method %in% c("GENIE3", "GRNBoost2", "JRF")) {
    stop("Invalid method. Choose either 'GENIE3', 'GRNBoost2', or 'JRF'.")
  }
  
  # Detect the number of available cores
  nCores <- parallel::detectCores()-1
  network_results <- list()
  
  # Loop over each matrix in the list
  for (j in seq_along(count_matrices_list)) {
    
    if (method == "GENIE3") {
      # Apply GENIE3 using dynamic nCores
      regulatory_network <- GENIE3(t(count_matrices_list[[j]]), nCores = nCores)
      regulatory_network <- getLinkList(regulatory_network)
      netout <- regulatory_network
      
    } else if (method == "GRNBoost2") {
      # Apply GRNBoost2
      count_matrix_df <- as.data.frame(count_matrices_list[[j]])
      genes <- colnames(count_matrix_df)
      
      df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), 
                                    columns = genes, 
                                    index = rownames(count_matrix_df))
      
      netout <- arboreto$grnboost2(df_pandas, gene_names = genes)
      
    } else if (method == "JRF") {
      # Normalize the matrix
      normalized_matrix <- t(apply(count_matrices_list[[j]], 1, function(x) {
        (x - mean(x)) / sd(x)
      }))
      
      jrf_matrices <- lapply(count_matrices_list, t)
      
      # Run JRF
      netout <- JRF(X = jrf_matrices, 
                    genes.name = rownames(jrf_matrices[[1]]), 
                    ntree = 500, 
                    mtry = round(sqrt(length(rownames(jrf_matrices[[1]])) - 1)))
    }
    
    network_results[[j]] <- netout
  }
  
  return(network_results)
}