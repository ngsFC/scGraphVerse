generate_adjacency <- function(link_list) {
  
  # Initialize an empty list to store adjacency matrices with weights
  adj_matrix_weight_list <- list()
  
  # Iterate over each matrix in the input list
  for (name in names(link_list)) {
    # Get the current link list
    link_list_genie3 <- link_list[[name]]
    
    # Get all unique gene names
    gene_names <- unique(c(link_list_genie3$Gene1, link_list_genie3$Gene2))
    
    # Create an empty adjacency matrix
    adj_matrix_genie3 <- matrix(0, nrow = length(gene_names), ncol = length(gene_names))
    rownames(adj_matrix_genie3) <- colnames(adj_matrix_genie3) <- gene_names
    
    # Fill in the adjacency matrix with the weight values
    for (i in 1:nrow(link_list_genie3)) {
      gene1 <- link_list_genie3$Gene1[i]
      gene2 <- link_list_genie3$Gene2[i]
      weight <- link_list_genie3$weight[i]
      
      adj_matrix_genie3[gene1, gene2] <- weight
      adj_matrix_genie3[gene2, gene1] <- weight  # Ensure symmetry
    }
    
    # Store the adjacency matrix in the list
    adj_matrix_weight_list[[name]] <- adj_matrix_genie3
  }
  
  return(adj_matrix_weight_list)
}

