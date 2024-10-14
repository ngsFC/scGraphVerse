generate_adjacency <- function(link_list) {
  adj_matrix_weight_list <- list()
  
  for (name in seq_along(link_list)) {
    link_list_genie3 <- link_list[[name]]
    gene_names <- unique(c(link_list_genie3$Gene1, link_list_genie3$Gene2))
    
    adj_matrix_genie3 <- matrix(0, nrow = length(gene_names), ncol = length(gene_names))
    rownames(adj_matrix_genie3) <- colnames(adj_matrix_genie3) <- gene_names
    
    for (i in 1:nrow(link_list_genie3)) {
      gene1 <- link_list_genie3$Gene1[i]
      gene2 <- link_list_genie3$Gene2[i]
      w <- link_list_genie3$weight[i]
      
      adj_matrix_genie3[gene1, gene2] <- w
      adj_matrix_genie3[gene2, gene1] <- w  # Ensure symmetry
    }
    
    adj_matrix_weight_list[[name]] <- adj_matrix_genie3
  }
  
  return(adj_matrix_weight_list)
}

