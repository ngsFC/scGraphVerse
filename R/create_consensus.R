create_consensus <- function(adj_matrix_list, method = "vote") {
  # Sum all adjacency matrices to get consensus values
  consensus_matrix <- Reduce("+", adj_matrix_list)
  
  if (method == "vote") {
    # Define threshold (present in at least 75% of the adjacency matrices)
    threshold <- round(length(adj_matrix_list) * 0.75)
    
    # Apply threshold: 1 if >= threshold, 0 if < threshold
    consensus_matrix[consensus_matrix < threshold] <- 0
    consensus_matrix[consensus_matrix >= threshold] <- 1
  } else if (method == "union") {
    # Union method: any position with a value >= 1 in any matrix is set to 1
    consensus_matrix[consensus_matrix >= 1] <- 1
  } else {
    stop("Invalid method. Choose 'vote' or 'union'.")
  }
  
  return(consensus_matrix)
}

