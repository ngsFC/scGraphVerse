cutoff_adjacency <- function(count_matrices, weighted_adjm_list, n, method = "GRNBoost2", weight_function = "mean") {
  
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
  all_percentile_values <- list()
  binary_adjm_list <- list()
  
  # Main loop through count matrices
  for (mat_index in seq_along(count_matrices)) {
    original_matrix <- count_matrices[[mat_index]]
    
    # Create shuffled matrices
    shuffled_matrices_list <- create_shuffled_matrices(original_matrix, n)
    
    all_grn_links <- list()
    
    # Compute 95th percentile for each shuffled matrix
    for (i in seq_along(shuffled_matrices_list)) {
      shuffled_matrix <- shuffled_matrices_list[[i]]
      
      # Infer network using the specified method
      network_results <- infer_networks(list(shuffled_matrix), method = method)
      
      # Make the result symmetric using the simmetric function
      symmetric_network <- simmetric(network_results, weight_function = weight_function)[[1]]
      
      # Get weights from the symmetric network and order them
      ordered_weights <- sort(symmetric_network$weight)
      
      # Calculate the 95th percentile of the ordered weights
      percentile_95 <- quantile(ordered_weights, 0.95)
      all_percentile_values[[mat_index]] <- c(all_percentile_values[[mat_index]], percentile_95)
    }
    
    # Compute mean of 95th percentiles (rounded to 3 decimal places)
    mean_percentile <- round(mean(unlist(all_percentile_values[[mat_index]])), 3)
    
    # Apply cutoff to the corresponding weighted adjacency matrix
    weighted_adjm <- weighted_adjm_list[[mat_index]]
    binary_adjm <- ifelse(weighted_adjm >= mean_percentile, 1, 0)
    binary_adjm_list[[mat_index]] <- binary_adjm
    
    # Print cutoff value for each matrix
    cat("Matrix", mat_index, "Mean 95th Percentile Cutoff:", mean_percentile, "\n")
  }
  
  # Plot histograms for cutoff values
  for (mat_index in seq_along(all_percentile_values)) {
    hist(unlist(all_percentile_values[[mat_index]]), main = paste("Histogram of Cutoff Values for Matrix", mat_index),
         xlab = "Cutoff Values", ylab = "Frequency", col = "lightblue", breaks = 20)
  }
  
  return(binary_adjm_list)
}

