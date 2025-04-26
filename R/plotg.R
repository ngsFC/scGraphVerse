#' Visualize Graphs from Adjacency Matrices
#'
#' Generates and arranges multiple graph visualizations from a list 
#' of adjacency matrices. Each matrix is converted into an undirected 
#' \pkg{igraph} object and visualized using a force-directed layout via \pkg{ggraph}.
#'
#' @param adj_matrix_list A list of square, symmetric adjacency matrices
#'   with zeros on the diagonal (no self-loops). Each matrix represents an undirected graph.
#'
#' @return
#' A grid of plots displaying all valid graphs in the input list.
#'
#' @details
#' Each adjacency matrix is validated to ensure it is square and symmetric.
#' Disconnected nodes (degree zero) are removed prior to visualization.
#' Graphs are visualized with a force-directed layout using \pkg{ggraph},
#' and multiple plots are arranged into a grid with \pkg{gridExtra}.
#'
#' Each subplot title includes the graph index, number of nodes, and number of edges.
#'
#' @note
#' This function requires the following packages to be installed:
#' \pkg{igraph}, \pkg{ggraph}, and \pkg{gridExtra}.
#' If any are missing, an informative error will be thrown.
#'
#' @importFrom igraph graph_from_adjacency_matrix delete_vertices V degree vcount ecount
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @export
#'
#' @examples
#' # Create two simple adjacency matrices
#' adj1 <- matrix(c(0,1,0,1,0,1,0,1,0), nrow = 3)
#' adj2 <- matrix(c(0,1,1,1,0,0,1,0,0), nrow = 3)
#'
#' # Visualize the graphs
#' plotg(list(adj1, adj2))

plotg <- function(adj_matrix_list) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required but not installed.")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("The 'ggraph' package is required but not installed.")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("The 'gridExtra' package is required but not installed.")
  }
  
  if (!is.list(adj_matrix_list)) {
    stop("Input must be a list of adjacency matrices.")
  }
  
  plot_list <- list()
  
  for (i in seq_along(adj_matrix_list)) {
    adj_matrix <- adj_matrix_list[[i]]
    
    if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
      warning(paste("Skipping graph", i, ": Adjacency matrix must be square."))
      next
    }
    
    if (!all(adj_matrix == t(adj_matrix))) {
      warning(paste("Skipping graph", i, ": Adjacency matrix must be symmetric."))
      next
    }
    
    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = NULL, diag = FALSE)
    
    # Remove disconnected nodes (degree 0)
    g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
    
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
      warning(paste("Skipping graph", i, ": No edges remaining after removing disconnected nodes."))
      next
    }
    
    plot_title <- paste("Graph", i, "\nNodes:", igraph::vcount(g), "Edges:", igraph::ecount(g))
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
  
  final_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  
  return(final_plot)
}

