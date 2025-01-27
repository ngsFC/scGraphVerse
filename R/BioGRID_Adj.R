#' Generate Weighted and Binary Adjacency Matrices from BioGRID Data
#'
#' This function creates weighted and binary adjacency matrices based on interactions 
#' from BioGRID data for the provided list of genes.
#'
#' @param top_genes A character vector of gene identifiers to filter the interactions.
#' @param biogrid_data A data frame containing BioGRID interaction data. The data should have 
#'   at least the columns `Interactor_A`, `Interactor_B`, and `Score`.
#' @return A list containing two adjacency matrices:
#'   \item{weighted}{A weighted adjacency matrix with interaction scores.}
#'   \item{binary}{A binary adjacency matrix indicating the presence of interactions (1 or 0).}
#' @examples
#' \dontrun{
#' adj_matrices <- BioGRID_Adj(top_genes = c("GENE1", "GENE2"), biogrid_data = my_biogrid_data)
#' }
#' @export
BioGRID_Adj <- function(top_genes, biogrid_data) {
  
  # Filter the data to include only interactions involving the top genes
  filtered_data <- biogrid_data %>%
    filter(Interactor_A %in% top_genes | Interactor_B %in% top_genes) %>%
    filter(Score != "-") %>%  # Ensure valid scores
    mutate(Score = as.numeric(Score))  # Convert scores to numeric
  
  if (nrow(filtered_data) == 0) {
    stop("No interactions found for the provided genes.")
  }
  
  # Get unique proteins involved in the interactions
  unique_proteins <- unique(c(filtered_data$Interactor_A, filtered_data$Interactor_B))
  
  # Create the weighted adjacency matrix
  weighted_adj_matrix <- matrix(0, 
                                nrow = length(unique_proteins), 
                                ncol = length(unique_proteins),
                                dimnames = list(unique_proteins, unique_proteins))
  
  # Populate the weighted adjacency matrix with interaction scores
  for (i in seq_len(nrow(filtered_data))) {
    a <- filtered_data$Interactor_A[i]
    b <- filtered_data$Interactor_B[i]
    score <- filtered_data$Score[i]
    weighted_adj_matrix[a, b] <- score
    weighted_adj_matrix[b, a] <- score
  }
  
  # Create the binary adjacency matrix based on interaction presence
  binary_adj_matrix <- ifelse(weighted_adj_matrix > 0, 1, 0)
  
  # Remove disconnected nodes by filtering rows and columns with no interactions
  connected_indices <- which(rowSums(weighted_adj_matrix) > 0 & colSums(weighted_adj_matrix) > 0)
  weighted_adj_matrix <- weighted_adj_matrix[connected_indices, connected_indices, drop = FALSE]
  binary_adj_matrix <- binary_adj_matrix[connected_indices, connected_indices, drop = FALSE]
  
  if (length(connected_indices) == 0) {
    stop("No connected nodes remain after filtering.")
  }
  
  return(list(weighted = weighted_adj_matrix, binary = binary_adj_matrix))
}

