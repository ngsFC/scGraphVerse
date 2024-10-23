infer_networks <- function(count_matrices_list, method = "GENIE3") {
  # Check if the method is valid
  if (!method %in% c("GENIE3", "GRNBoost2", "JRF")) {
    stop("Invalid method. Choose either 'GENIE3', 'GRNBoost2', or 'JRF'.")
  }
  
  # Prepare the result list
  network_results <- list()
  
  # Loop over each matrix in the list
  for (j in seq_along(count_matrices_list)) {
    
    if (method == "GENIE3") {
      # Apply GENIE3 on the j-th count matrix
      regulatory_network <- GENIE3(t(count_matrices_list[[j]]))
      netout <- getLinkList(regulatory_network)
      network_results[[j]] <- netout
      
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
      network_results[[j]] <- netout
    
    } else if (method == "JRF") {
      # Apply JRF method on the list of matrices
      # JRF expects a list of transposed matrices
      jrf_matrices <- lapply(count_matrices_list, t)
      
      # Run JRF with the transposed matrices
      result <- JRF(X = jrf_matrices, 
                    genes.name = rownames(jrf_matrices[[1]]), 
                    ntree = 500, 
                    mtry = round(sqrt(length(rownames(jrf_matrices[[1]])) - 1)))
      
      first_two_columns <- result[, 1:2]
      
      for (i in 1:length(jrf_matrices)) {
        specific_column <- result[, 2 + i]
        netout <- cbind(first_two_columns, specific_column)
        network_results[[i]] <- netout
      }
    }
  }
  
  # Return the full list of network results
  return(network_results)
}

