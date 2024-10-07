adj_cutoff <- function(adj_matrix_list, weight_quantile = 0.70) {
  
  # Initialize an empty list to store adjacency matrices after applying cutoff
  adj_matrix_binary_list <- list()
  plot_list <- list()
  
  # Iterate over each adjacency matrix
  for (name in names(adj_matrix_list)) {
    # Get the current adjacency matrix
    adj_matrix_genie3 <- adj_matrix_list[[name]]
    
    # Get the weights from the upper triangle of the matrix (excluding diagonal)
    weights <- as.vector(adj_matrix_genie3[upper.tri(adj_matrix_genie3)])
    
    # Calculate the weight cutoff based on the specified quantile
    weight_cutoff <- quantile(weights, weight_quantile, na.rm = TRUE)
    
    # Create the weight distribution plot with cutoff line
    weight_distribution_plot <- ggplot(data.frame(weight = weights), aes(x = weight)) +
      geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
      geom_vline(aes(xintercept = weight_cutoff), color = "red", linetype = "dashed", size = 1) +
      theme_minimal() +
      labs(title = paste("Distribution of Weights for", name),
           x = "Weight",
           y = "Frequency",
           caption = paste("Cut-off point:", round(weight_cutoff, 3))) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Store the plot in the plot list
    plot_list[[name]] <- weight_distribution_plot
    
    # Create a binary matrix by applying the cutoff
    adj_matrix_binary <- adj_matrix_genie3 >= weight_cutoff
    diag(adj_matrix_binary) <- 1  # Optional: Set diagonal to 1 to indicate self-connections
    
    # Convert logical matrix to numeric (0 and 1)
    adj_matrix_binary <- as.numeric(adj_matrix_binary)
    dim(adj_matrix_binary) <- dim(adj_matrix_genie3)
    rownames(adj_matrix_binary) <- rownames(adj_matrix_genie3)
    colnames(adj_matrix_binary) <- colnames(adj_matrix_genie3)
    
    # Store the binary adjacency matrix in the list
    adj_matrix_binary_list[[name]] <- adj_matrix_binary
  }
  
  return(list(binary_matrices = adj_matrix_binary_list, plots = plot_list))
}

