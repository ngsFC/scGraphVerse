#' Plot Multiple Graphs from Adjacency Matrices
#'
#' This function generates and visualizes graphs from a list of adjacency matrices.
#' Each adjacency matrix is converted into an `igraph` object and plotted using the
#' `ggraph` package with a force-directed layout. The generated plots are then 
#' arranged in a grid for easy comparison. 
#'
#' @param adj_matrix_list A list of adjacency matrices, where each matrix represents
#'   an undirected graph. The function assumes that the matrices are symmetric and
#'   have zeros on the diagonal (no self-loops).
#' 
#' @details
#' Each adjacency matrix in `adj_matrix_list` is processed individually, converting it
#' into an `igraph` object. The graphs are then visualized using the `ggraph` package, which
#' employs a force-directed layout for graph visualization. The title of each plot displays the
#' graph index and basic graph statistics (number of nodes and edges). All individual graphs are
#' arranged into a grid layout, with the number of columns dynamically determined by the square root
#' of the total number of graphs.
#'
#' @note
#' This function requires the following packages:
#'   - `igraph`
#'   - `ggraph`
#'   - `gridExtra`
#'
#' If any of these packages are not installed, the function will stop and display an error message.
#'
#' @return A combined plot grid of all graphs generated from the adjacency matrices.
#' Each graph is visualized with nodes and edges, and the title includes the graph's index,
#' number of nodes, and number of edges.
#'
#' @examples
#' # Example adjacency matrices (binary, undirected)
#' adj_matrix_1 <- matrix(c(0, 1, 0, 1, 1, 0), nrow = 3, byrow = TRUE)
#' adj_matrix_2 <- matrix(c(0, 1, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
#'
#' # Create a list of adjacency matrices
#' adj_matrices <- list(adj_matrix_1, adj_matrix_2)
#'
#' # Call the plotg function to visualize the graphs
#' plotg(adj_matrices)
#'
#' @import igraph
#' @import ggraph
#' @import gridExtra
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

