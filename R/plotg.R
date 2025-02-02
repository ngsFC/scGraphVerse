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
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom gridExtra grid.arrange
#'
#' @export

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
  
  # Validate input
  if (!is.list(adj_matrix_list)) {
    stop("Input must be a list of adjacency matrices.")
  }
  
  plot_list <- list()
  
  for (i in seq_along(adj_matrix_list)) {
    adj_matrix <- adj_matrix_list[[i]]
    
    # Validate matrix
    if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
      warning(paste("Skipping graph", i, ": Adjacency matrix must be square."))
      next
    }
    
    if (!all(adj_matrix == t(adj_matrix))) {
      warning(paste("Skipping graph", i, ": Adjacency matrix must be symmetric."))
      next
    }
    
    # Convert adjacency matrix to igraph object
    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = NULL, diag = FALSE)
    
    # Remove disconnected nodes (degree 0)
    g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
    
    # Ensure the graph has edges before plotting
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
      warning(paste("Skipping graph", i, ": No edges remaining after removing disconnected nodes."))
      next
    }
    
    # Graph title
    plot_title <- paste("Graph", i, "\nNodes:", igraph::vcount(g), "Edges:", igraph::ecount(g))
    
    # Create the plot
    p <- ggraph::ggraph(g, layout = "fr") + 
      ggraph::geom_edge_link(color = "gray", width = 0.5) +  
      ggraph::geom_node_point(color = "steelblue", size = 3) +  
      ggplot2::labs(title = plot_title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  if (length(plot_list) == 0) {
    stop("No valid graphs were generated.")
  }
  
  # Combine all plots into a grid
  final_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  
  return(final_plot)
}

