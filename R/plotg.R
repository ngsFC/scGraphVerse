plotg <- function(adj_matrix_list) {
  # Load required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required but not installed.")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("The 'ggraph' package is required but not installed.")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("The 'gridExtra' package is required but not installed.")
  }
  
  library(igraph)
  library(ggraph)
  library(gridExtra)
  
  # List to store individual plots
  plot_list <- list()
  
  # Loop through each adjacency matrix and create plots
  for (i in seq_along(adj_matrix_list)) {
    adj_matrix <- adj_matrix_list[[i]]
    
    # Convert adjacency matrix to igraph object
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = NULL, diag = FALSE)
    
    # Create the plot using ggraph
    plot_title <- paste("Graph", i, "\nNodes:", vcount(g), "Edges:", ecount(g))
    
    p <- ggraph(g, layout = "fr") + 
      geom_edge_link(color = "gray", width = 0.5) +  # Set edge color and width
      geom_node_point(color = "steelblue", size = 0.7) +  # Set node color and size
      labs(title = plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"  # Remove legends to avoid conflicts
      )
    
    # Add the plot to the list
    plot_list[[i]] <- p
  }
  
  # Combine all plots into a grid
  final_plot <- grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  
  # Print the combined plot
  print(final_plot)
}
