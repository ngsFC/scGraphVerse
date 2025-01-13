BioGRID_Adj <- function(top_genes, biogrid_data) {
  
  filtered_data <- biogrid_data %>%
    filter(Interactor_A %in% top_genes | Interactor_B %in% top_genes) %>%
    filter(Score != "-") %>%  # Ensure valid scores
    mutate(Score = as.numeric(Score))  # Convert scores to numeric
  
  if (nrow(filtered_data) == 0) {
    stop("No interactions found for the provided genes.")
  }
  
  unique_proteins <- unique(c(filtered_data$Interactor_A, filtered_data$Interactor_B))
  
  weighted_adj_matrix <- matrix(0, 
                                nrow = length(unique_proteins), 
                                ncol = length(unique_proteins),
                                dimnames = list(unique_proteins, unique_proteins))
  
  for (i in seq_len(nrow(filtered_data))) {
    a <- filtered_data$Interactor_A[i]
    b <- filtered_data$Interactor_B[i]
    score <- filtered_data$Score[i]
    weighted_adj_matrix[a, b] <- score
    weighted_adj_matrix[b, a] <- score
  }
  
  binary_adj_matrix <- ifelse(weighted_adj_matrix > 0, 1, 0)
  
  connected_indices <- which(rowSums(weighted_adj_matrix) > 0 & colSums(weighted_adj_matrix) > 0)
  weighted_adj_matrix <- weighted_adj_matrix[connected_indices, connected_indices, drop = FALSE]
  binary_adj_matrix <- binary_adj_matrix[connected_indices, connected_indices, drop = FALSE]
  
  if (length(connected_indices) == 0) {
    stop("No connected nodes remain after filtering.")
  }
  
  return(list(weighted = weighted_adj_matrix, binary = binary_adj_matrix))
}
