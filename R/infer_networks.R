infer_networks <- function(count_matrices_list, method = "GENIE3") {
  # Check if the method is valid
  if (!method %in% c("GENIE3", "GRNBoost2")) {
    stop("Invalid method. Choose either 'GENIE3' or 'GRNBoost2'.")
  }
  
  # Prepare the result list
  network_results <- list()
  
  # Loop over each matrix in the list
  for (j in seq_along(count_matrices_list)) {
    network_j <- list()  # To store results for each count matrix
    
    if (method == "GENIE3") {
      # Apply GENIE3 on the j-th count matrix
      regulatory_network <- GENIE3(t(count_matrices_list[[j]]))
      netout <- getLinkList(regulatory_network)
      network_j[[paste0("Matrix_", j)]] <- netout
      
    } else if (method == "GRNBoost2") {
      # Convert the R matrix to a pandas DataFrame
      count_matrix_df <- as.data.frame(count_matrices_list[[j]])
      genes <- colnames(count_matrix_df)
      
      # Convert to pandas DataFrame using reticulate
      df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), 
                                    columns = genes, 
                                    index = rownames(count_matrix_df))
      
      # Apply GRNBoost2 on the pandas DataFrame
      netout <- arboreto$grnboost2(df_pandas, gene_names = genes)
      network_j[[paste0("Matrix_", j)]] <- netout
    }
    
    # Store results for the j-th matrix
    network_results[[paste0("Adjacency_", j)]] <- network_j
  }
  
  # Return the full list of network results
  return(network_results)
}

