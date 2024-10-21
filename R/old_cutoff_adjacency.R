cutoff_adjacency <- function(count_matrices, weighted_adjm_list, n, method = "GRNBoost2", weight_function = "mean") {
  
  shuffle_rows <- function(matrix, seed_vector) {
    shuffled_matrix <- matrix
    for (i in 1:nrow(matrix)) {
      set.seed(seed_vector[i])
      shuffled_matrix[i, ] <- sample(matrix[i, ])
    }
    return(shuffled_matrix)
  }
  
  create_shuffled_matrices <- function(original_matrix, n) {
    shuffled_matrices <- list()
    for (i in 1:n) {
      seed_vector <- sample(1:10000, nrow(original_matrix))
      shuffled_matrix <- shuffle_rows(original_matrix, seed_vector)
      shuffled_matrices[[i]] <- shuffled_matrix
    }
    return(shuffled_matrices)
  }
  
  all_percentile_values <- list()
  binary_adjm_list <- list()
  
  for (mat_index in seq_along(count_matrices)) {
    original_matrix <- count_matrices[[mat_index]]
    
    shuffled_matrices_list <- create_shuffled_matrices(original_matrix, n)
    
    all_grn_links <- list()
    
    if (method == "GRNBoost2") {
      # Run GRNBoost2 on each shuffled matrix
      for (i in 1:length(shuffled_matrices_list)) {
        count_matrix_df <- as.data.frame(shuffled_matrices_list[[i]])
        genes <- colnames(count_matrix_df)
        
        # Convert to pandas DataFrame (assuming Python interface setup for GRNBoost2)
        df_pandas <- pandas$DataFrame(data = count_matrix_df, columns = genes, index = rownames(count_matrix_df))
        grn_links <- arboreto$grnboost2(df_pandas, gene_names = genes)
        all_grn_links[[i]] <- grn_links
      }
      
    } else if (method == "GENIE3") {
      # Run GENIE3 on each shuffled matrix
      for (i in seq_along(shuffled_matrices_list)) {
        count_matrix_transposed <- t(shuffled_matrices_list[[i]])
        regulatory_network_genie3 <- GENIE3(count_matrix_transposed, nCores=16)
        genie3out <- getLinkList(regulatory_network_genie3)
        all_grn_links[[i]] <- genie3out
      }
    } else {
      stop("Invalid method. Please choose either 'GRNBoost2' or 'GENIE3'.")
    }
    
    all_grn_links <- simmetric(all_grn_links, weight_function = weight_function)
    
    percentile_values <- sapply(all_grn_links, function(grn_result) {
      ordered_weights <- sort(grn_result$weight)
      threshold <- quantile(ordered_weights, 0.95)
      return(threshold)
    })
    
    all_percentile_values[[mat_index]] <- percentile_values
    
    mean_value <- mean(percentile_values)
    
    weighted_adjm <- weighted_adjm_list[[mat_index]]
    binary_adjm <- ifelse(weighted_adjm > mean_value, 1, 0)
    
    binary_adjm_list[[mat_index]] <- binary_adjm
  }
  
  combined_percentile_values <- unlist(all_percentile_values)
  overall_mean_value <- mean(combined_percentile_values)
  
  hist(combined_percentile_values, main = "Histogram of Percentile Values", 
       xlab = "Percentile Values", col = "grey")
  abline(v = overall_mean_value, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("Mean:", round(overall_mean_value, 3)), col = "red", lwd = 2, lty = 2)
  
  return(list(binary_adjm_list = binary_adjm_list, percentile_values = all_percentile_values, overall_mean_value = overall_mean_value))
}
