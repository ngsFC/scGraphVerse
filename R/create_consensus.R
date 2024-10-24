create_consensus <- function(adj_matrix_list) {
  # Sum all adjacency matrices to get consensus values
  consensus_matrix <- Reduce("+", adj_matrix_list)
  
  # Define threshold (present in at least 75% of the adjacency matrices)
  threshold <- round(length(adj_matrix_list) * 0.75)
  
  # Modify the consensus matrix in place: 1 if >= threshold, 0 if < threshold
  consensus_matrix[consensus_matrix < threshold] <- 0
  consensus_matrix[consensus_matrix >= threshold] <- 1
  
  return(consensus_matrix)
}
