generate_adjacency <- function(link_list, weight_quantile = 0.70, output_rds_file = "genie3_adjacency_matrices.rds") {
  
  # Initialize an empty list to store adjacency matrices
  adj_matrix_list <- list()
  # Initialize an empty list to store plots
  plot_list <- list()
  
  # Iterate over each matrix in the input list
  for (name in names(link_list)) {
    # Get the current link list
    link_list_genie3 <- link_list[[name]]
    
    # Get all unique gene names
    gene_names <- unique(c(link_list_genie3$Gene1, link_list_genie3$Gene2))
    
    # Create an empty adjacency matrix
    adj_matrix_genie3 <- matrix(0, nrow = length(gene_names), ncol = length(gene_names))
    rownames(adj_matrix_genie3) <- colnames(adj_matrix_genie3) <- gene_names
    
    # Calculate the weight cutoff based on the specified quantile
    weight_cutoff <- quantile(link_list_genie3$weight, weight_quantile)
    
    # Create the weight distribution plot with cutoff line
    weight_distribution_plot <- ggplot(data = link_list_genie3, aes(x = weight)) +
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
    
    # Fill in the adjacency matrix based on the weight cutoff
    for (i in 1:nrow(link_list_genie3)) {
      gene1 <- link_list_genie3$Gene1[i]
      gene2 <- link_list_genie3$Gene2[i]
      weight <- link_list_genie3$weight[i]
      
      if (weight >= weight_cutoff) {
        adj_matrix_genie3[gene1, gene2] <- 1
        adj_matrix_genie3[gene2, gene1] <- 1  # Ensure symmetry
      }
    }
    
    # Set diagonal to 1 (optional, depending on your specific requirements)
    diag(adj_matrix_genie3) <- 1
    
    # Store the adjacency matrix in the list
    adj_matrix_list[[name]] <- adj_matrix_genie3
  }
 
  return(adj_matrix_list)

  # Save the list of adjacency matrices to an RDS file
  saveRDS(adj_matrix_list, output_rds_file)
  
  # Arrange all plots in a grid layout
  grid.arrange(grobs = plot_list, ncol = 2)  # Adjust 'ncol' to arrange plots in desired number of columns
  
  cat("Adjacency matrices saved to:", output_rds_file, "\n")
}

