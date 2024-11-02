cutoff_adjacency <- function(count_matrices, weighted_adjm_list, ground.truth, n, method = "GENIE3", weight_function = "mean") {
  
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
    
    # Create shuffled matrices for the current original matrix
    shuffled_matrices_list <- create_shuffled_matrices(original_matrix, n)
    
    # If method is "JRF", infer network across all shuffled matrices at once
    if (method == "JRF") {
      jrf_mat <- infer_networks(shuffled_matrices_list, method = method)
      network_results <- list()
      
      importance_columns <- grep("importance", names(jrf_mat[[1]]), value = TRUE)
      
      for (i in seq_along(importance_columns)) {
        df <- jrf_mat[[1]][, c("gene1", "gene2", importance_columns[i])]
        names(df)[3] <- importance_columns[i]
        network_results[[i]] <- df
      }
      
      # Process the JRF output by iterating over each inferred network in the results
      for (i in seq_along(network_results)) {
        # Make each result symmetric and calculate the 95th percentile
        network_results_adjm <- generate_adjacency(network_results[[i]], ground.truth)
        symmetric_network <- symmetrize(network_results_adjm, weight_function = weight_function)
        symmetric_network <- symmetric_network[[1]]
        
        # Extract weights and compute the 95th percentile for the symmetric network
        upper_triangle_weights <- symmetric_network[upper.tri(symmetric_network)]
        ordered_weights <- sort(upper_triangle_weights, decreasing = TRUE)
        percentile_95 <- quantile(ordered_weights, 0.95)
        
        all_percentile_values[[i]] <- percentile_95
      }
      
    } else {
      # For other methods, run network inference within the loop for each shuffled matrix
      for (i in seq_along(shuffled_matrices_list)) {
        shuffled_matrix <- shuffled_matrices_list[[i]]
        network_results <- infer_networks(list(shuffled_matrix), method = method)
        
        # Make the result symmetric using the symmetrize function
        network_results_adjm <- generate_adjacency(network_results, ground.truth)
        symmetric_network <- symmetrize(network_results_adjm, weight_function = weight_function)
        symmetric_network <- symmetric_network[[1]]
        
        # Extract weights and compute the 95th percentile for the symmetric network
        upper_triangle_weights <- symmetric_network[upper.tri(symmetric_network)]
        ordered_weights <- sort(upper_triangle_weights, decreasing = TRUE)
        percentile_95 <- quantile(ordered_weights, 0.95)
        
        all_percentile_values[[i]] <- percentile_95
      }
    }
    
    # Compute mean of 95th percentiles (rounded to 3 decimal places)
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