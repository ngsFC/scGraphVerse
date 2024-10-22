generate_adjacency <- function(link_list) {
  adj_matrix_weight_list <- list()
  
  for (i in seq_along(link_list)) {
    link_df <- as.data.frame(link_list[[i]])
    gene_names <- unique(c(link_df$Gene1, link_df$Gene2))
    
    adj_matrix <- matrix(0, nrow = length(gene_names), ncol = length(gene_names))
    rownames(adj_matrix) <- colnames(adj_matrix) <- gene_names
    
    for (j in 1:nrow(link_df)) {
      gene1 <- link_df$Gene1[j]
      gene2 <- link_df$Gene2[j]
      w <- link_df$weight[j]
      
      adj_matrix[gene1, gene2] <- w
      adj_matrix[gene2, gene1] <- w  # Ensure symmetry
    }
    
    adj_matrix_weight_list[[i]] <- adj_matrix
  }
  
  return(adj_matrix_weight_list)
}

