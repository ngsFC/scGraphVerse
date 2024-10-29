generate_adjacency <- function(df_list) {
  # Initialize a list to store the adjacency matrices
  adj_matrix_list <- list()
  
  # Loop over each data frame in the input list
  for (k in seq_along(df_list)) {
    data <- df_list[[k]]
    
    # Extract unique genes to determine matrix dimensions
    unique_genes <- sort(unique(c(data$Gene1, data$Gene2)))
    p <- length(unique_genes)
    
    # Initialize an empty matrix with zeroes
    adj_matrix <- matrix(0, nrow = p, ncol = p)
    rownames(adj_matrix) <- colnames(adj_matrix) <- unique_genes
    
    # Fill the matrix with weights
    for (i in 1:nrow(data)) {
      gene1 <- data$Gene1[i]
      gene2 <- data$Gene2[i]
      weight <- data$weight[i]
      
      # Find the row and column indices for each gene
      row_index <- which(rownames(adj_matrix) == gene1)
      col_index <- which(colnames(adj_matrix) == gene2)
      
      # Assign the weight to both [gene1, gene2] and [gene2, gene1] to ensure symmetry
      adj_matrix[row_index, col_index] <- weight
      adj_matrix[col_index, row_index] <- weight
    }
    
    # Add the matrix to the list
    adj_matrix_list[[k]] <- adj_matrix
  }
  
  return(adj_matrix_list)
}
