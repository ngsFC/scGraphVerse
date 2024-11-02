cutoff_adjacency <- function(count_matrices, weighted_adjm_list, ground.truth, n, method = "GRNBoost2", weight_function = "mean") {
  
  # Function to shuffle rows of a matrix
  shuffle_rows <- function(matrix, seed_vector) {
    shuffled_matrix <- matrix
    for (i in 1:nrow(matrix)) {
      set.seed(seed_vector[i])
      shuffled_matrix[i, ] <- sample(matrix[i, ])
    }
    return(shuffled_matrix)
  }
  
  # Function to create n shuffled matrices
  create_shuffled_matrices <- function(original_matrix, n) {
    shuffled_matrices <- list()
    for (i in 1:n) {
      seed_vector <- sample(1:10000, nrow(original_matrix))
      shuffled_matrix <- shuffle_rows(original_matrix, seed_vector)
      shuffled_matrices[[i]] <- shuffled_matrix
    }
    return(shuffled_matrices)
  }
  
  # Initialize lists to store results
  binary_adjm_list <- list()
  
  # Main loop through count matrices
  for (mat_index in seq_along(count_matrices)) {
    original_matrix <- count_matrices[[mat_index]]
    
    # Initialize a list for storing percentiles for this specific matrix
    all_percentile_values <- list()
    
    # Create shuffled matrices
    shuffled_matrices_list <- create_shuffled_matrices(original_matrix, n)
    
    # Variable to store inferred network results
    if (method == "JRF") {
      # Apply JRF directly to the list of shuffled matrices
      jrf_mat <- infer_networks(shuffled_matrices_list, method = "JRF")
      
      # Transform JRF output into jrf_list format
      jrf_list <- list()
      importance_columns <- grep("importance", names(jrf_mat[[1]]), value = TRUE)
      
      for (i in seq_along(importance_columns)) {
        # Select the 'gene1', 'gene2', and the current 'importance' column
        df <- jrf_mat[[1]][, c("gene1", "gene2", importance_columns[i])]
        
        # Rename the importance column to its original name (e.g., importance1, importance2, etc.)
        names(df)[3] <- importance_columns[i]
        
        # Add the data frame to the output list
        jrf_list[[i]] <- df
      }
      
      # Generate adjacency matrices and calculate percentiles for each matrix in jrf_list
      for (k in seq_along(jrf_list)) {
        network_results_adjm <- generate_adjacency(list(jrf_list[[k]]), ground.truth)
        symmetric_network <- symmetrize(network_results_adjm, weight_function = weight_function)
        symmetric_network <- symmetric_network[[1]]
        
        # Extract and order weights from the symmetric network
        upper_triangle_weights <- symmetric_network[upper.tri(symmetric_network)]
        ordered_weights <- sort(upper_triangle_weights, decreasing = TRUE)
        
        # Calculate the 95th percentile and add to the list
        percentile_95 <- quantile(ordered_weights, 0.95)
        all_percentile_values <- c(all_percentile_values, percentile_95)
      }
    } else {
      # For other methods (GENIE3, GRNBoost2), infer networks for each shuffled matrix
      for (shuffled_matrix in shuffled_matrices_list) {
        network_results <- infer_networks(list(shuffled_matrix), method = method)
        
        # Make the result symmetric using the symmetrize function
        network_results_adjm <- generate_adjacency(network_results, ground.truth)
        symmetric_network <- symmetrize(network_results_adjm, weight_function = weight_function)
        symmetric_network <- symmetric_network[[1]]
        
        # Extract and order weights from the symmetric network
        upper_triangle_weights <- symmetric_network[upper.tri(symmetric_network)]
        ordered_weights <- sort(upper_triangle_weights, decreasing = TRUE)
        
        # Calculate the 95th percentile and add to the list
        percentile_95 <- quantile(ordered_weights, 0.95)
        all_percentile_values <- c(all_percentile_values, percentile_95)
      }
    }
    
    # Compute mean of 95th percentiles
    mean_percentile <- mean(unlist(all_percentile_values))
    
    # Apply cutoff to the corresponding weighted adjacency matrix
    weighted_adjm <- weighted_adjm_list[[mat_index]]
    binary_adjm <- ifelse(weighted_adjm >= mean_percentile, 1, 0)
    binary_adjm_list[[mat_index]] <- binary_adjm
    
    # Print cutoff value for each matrix
    cat("Matrix", mat_index, "Mean 95th Percentile Cutoff:", mean_percentile, "\n")
  }
  return(binary_adjm_list)
}
