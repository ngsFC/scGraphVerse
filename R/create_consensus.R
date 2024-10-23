create_consensus <- function(adj_matrix_list) {
  # Sum all adjacency matrices to get consensus values
  consensus_matrix <- Reduce("+", adj_matrix_list)
  
  # Define threshold (present in at least 75% of the adjacency matrices)
  threshold <- round(length(adj_matrix_list) * 0.75)
  
  # Create binary consensus matrix (1 if present in >= threshold matrices)
  consensus_matrix_binary <- consensus_matrix >= threshold
  return(consensus_matrix_binary)
}
