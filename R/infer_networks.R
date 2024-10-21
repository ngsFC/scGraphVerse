infer_networks <- function(count_matrices_list, method = "GENIE3") {
  
  if (!method %in% c("GENIE3", "GRNBoost2")) {
    stop("Invalid method. Choose either 'GENIE3' or 'GRNBoost2'.")
  }
  
  network_results <- list()
  
  for (j in seq_along(count_matrices_list)) {
    matrices_j <- count_matrices_list[[j]]
    network_j <- list()
    
    if (method == "GENIE3") {
      # Apply GENIE3 to each count matrix
      for (i in seq_along(matrices_j)) {
        regulatory_network_genie3 <- GENIE3(t(matrices_j))  # Transpose for GENIE3
        genie3out <- getLinkList(regulatory_network_genie3)
        network_j[[paste0("Matrix_", i)]] <- genie3out
      }
      
    } else if (method == "GRNBoost2") {
      # Apply GRNBoost2 to each count matrix
      for (i in seq_along(matrices_j)) {
        count_matrix_df <- as.data.frame(matrices_j)
        genes <- colnames(count_matrix_df)
        
        df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), columns = genes, index = rownames(count_matrix_df))
        grn_links <- arboreto$grnboost2(df_pandas, gene_names = genes)
        network_j[[paste0("Matrix_", i)]] <- grn_links
      }
    }
    
    network_results[[paste0("Adjacency_", j)]] <- network_j
  }
  
  return(network_results)
}
