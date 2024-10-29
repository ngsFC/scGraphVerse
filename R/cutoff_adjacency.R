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
    
    all_grn_links <- list()
    
    # Compute 95th percentile for each shuffled matrix
    for (i in seq_along(shuffled_matrices_list)) {
      shuffled_matrix <- shuffled_matrices_list[[i]]
      
      # Infer network using the specified method
      network_results <- infer_networks(list(shuffled_matrix), method = method)
      
      # Handle the output for "JRF"
      if (method == "JRF") {
        # JRF produces a list of n networks, each with 2+ columns
        # We need to process the first network result, which will have 2 columns + the third column for weights
        network_results <- network_results[[1]]
      }
      
      # Make the result symmetric using the symmetrize function
      network_results_adjm <- generate_adjacency(network_results, ground.truth)
      symmetric_network <- symmetrize(network_results_adjm, weight_function = weight_function)
      symmetric_network <- symmetric_network[[1]]
      # Get weights from the symmetric network and order them
      upper_triangle_weights <- symmetric_network[upper.tri(symmetric_network)]
      ordered_weights <- sort(upper_triangle_weights, decreasing = TRUE)
      
      # Calculate the 95th percentile of the ordered weights
      percentile_95 <- quantile(ordered_weights, 0.95)
      all_percentile_values[[i]] <- percentile_95
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
 