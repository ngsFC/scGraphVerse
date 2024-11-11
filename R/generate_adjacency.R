generate_adjacency <- function(df_list, ground.truth) {
  adjacency_matrix_list <- list()
  
  # Loop over each data frame in the input list
  for (k in seq_along(df_list)) {
    data <- df_list[[k]]
    
    # Initialize an empty adjacency matrix with the same dimensions as ground.truth
    adjacency_matrix <- matrix(0, nrow = nrow(ground.truth), ncol = ncol(ground.truth))
    rownames(adjacency_matrix) <- rownames(ground.truth)
    colnames(adjacency_matrix) <- colnames(ground.truth)
    
    # Populate the adjacency matrix based on the first, second, and third columns of the data frame
    for (i in 1:nrow(data)) {
      gene1 <- as.character(data[i, 1])  # First column: Gene1
      gene2 <- as.character(data[i, 2])  # Second column: Gene2
      weight <- data[i, 3]               # Third column: Weight
      
      # Place the weight in the adjacency matrix in the correct position
      if (gene1 %in% rownames(adjacency_matrix) && gene2 %in% colnames(adjacency_matrix)) {
        adjacency_matrix[gene1, gene2] <- weight
      }
    }
    
    # Ensure the diagonal is set to 0, matching ground.truth
    diag(adjacency_matrix) <- 0
    
    # Add the adjacency matrix to the list
    adjacency_matrix_list[[k]] <- adjacency_matrix
  }
  
  return(adjacency_matrix_list)
}