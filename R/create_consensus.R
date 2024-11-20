create_consensus <- function(adj_matrix_list, method = "vote", weighted_list = NULL, theta = 0.04, threshold = 0.5, ncores = 1) {
  # Sum all adjacency matrices to get consensus values (default for 'vote' and 'union')
  consensus_matrix <- Reduce("+", adj_matrix_list)
  
  if (method == "vote") {
    # Define threshold (present in at least 75% of the adjacency matrices)
    threshold_value <- round(length(adj_matrix_list) * 0.75)
    
    # Apply threshold: 1 if >= threshold, 0 if < threshold
    consensus_matrix[consensus_matrix < threshold_value] <- 0
    consensus_matrix[consensus_matrix >= threshold_value] <- 1
  } else if (method == "union") {
    # Union method: any position with a value >= 1 in any matrix is set to 1
    consensus_matrix[consensus_matrix >= 1] <- 1
  } else if (method == "INet") {
    # Ensure weighted list is provided
    if (is.null(weighted_list)) {
      stop("For method 'INet', a weighted_list must be provided.")
    }
    
    # Element-wise replacement of 1s in binary matrices with corresponding weights
    list_weighted <- mapply(function(binary_mat, weighted_mat) {
      if (!is.matrix(binary_mat) || !is.matrix(weighted_mat)) {
        stop("Both binary and weighted elements must be matrices.")
      }
      if (!all(dim(binary_mat) == dim(weighted_mat))) {
        stop("Dimension mismatch between binary adjacency matrix and weighted matrix.")
      }
      # Replace 1s in binary_mat with weights from weighted_mat
      newmat <- binary_mat * weighted_mat
      return(newmat)
    }, adj_matrix_list, weighted_list, SIMPLIFY = FALSE)
    
    # Normalize the weighted matrices
    list_genie3_norm <- lapply(list_weighted, function(mat) mat / max(mat, na.rm = TRUE))
    
    # Generate consensus network using INet
    Con <- consensusNet(list_genie3_norm, theta = theta, ncores = ncores, threshold = threshold)
    
    # Convert the consensus network graph to an adjacency matrix
    consensus_matrix <- as.matrix(as_adjacency_matrix(Con$graphConsensus))
  } else {
    stop("Invalid method. Choose 'vote', 'union', or 'INet'.")
  }
  
  return(consensus_matrix)
}
