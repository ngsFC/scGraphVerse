#' Compare Consensus and Original Graphs
#'
#' This function compares a consensus adjacency matrix with an original adjacency matrix, visualizing the true 
#' positives (TP), false negatives (FN), and false positives (FP) between the two graphs. It uses the 
#' `ggraph` package to create visualizations that show the agreement and disagreement between the two graphs.
#' The function generates two plots:
#' 1. A plot showing the edges present in both the ground truth and consensus matrix (TP) and the edges 
#'    present only in the ground truth (FN).
#' 2. A plot showing the edges present in the consensus matrix but absent in the ground truth (FP).
#'
#' @param consensus_matrix A binary adjacency matrix representing the consensus graph.
#' @param original_matrix A binary adjacency matrix representing the original or ground truth graph.
#'
#' @return This function does not return any values. It generates two plots:
#'         - One for true positives and false positives.
#'         - Another for false negatives.
#'
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom gridExtra grid.arrange
#'
#' @export
compare_consensus <- function(consensus_matrix, original_matrix) {
  
  # Helper function to create graphs from adjacency matrices
  create_graph <- function(adj_matrix) {
    graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
    return(graph)
  }
  
  # Helper function for ggraph plot
  plot_ggraph <- function(graph, plot_title, edge_colors) {
    ggraph(graph, layout = "fr") +
      geom_edge_link(aes(color = edge_colors), width = 0.5) +
      geom_node_point(color = "steelblue", size = 0.7) +
      scale_color_manual(values = c("red" = "red", "blue" = "blue")) +  # Map red and blue
      labs(title = plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
  }
  
  # Create igraph objects from adjacency matrices
  graph_original <- create_graph(original_matrix)
  graph_consensus <- create_graph(consensus_matrix)
  
  # Extract edge lists
  original_edges <- as_edgelist(graph_original)
  consensus_edges <- as_edgelist(graph_consensus)
  
  # Convert edge lists to strings for easy comparison
  original_edges_str <- apply(original_edges, 1, function(x) paste(sort(x), collapse = "-"))
  consensus_edges_str <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse = "-"))
  
  # Determine edge colors for TP and FP
  edge_colors <- ifelse(original_edges_str %in% consensus_edges_str, "red", "blue")  # Red = TP, Blue = FP
  
  # True Positives (TP) and False Positives (FP)
  TP_count <- sum(edge_colors == "red")
  FP_count <- sum(edge_colors == "blue")
  
  # Plot ground truth with TP (red) and FP (blue)
  set.seed(1234)
  plot_comparison <- plot_ggraph(graph_original, 
                                 plot_title = paste("Ground Truth & Consensus\nTP:", TP_count, "FP:", FP_count), 
                                 edge_colors = edge_colors)
  
  # False Negatives (FN) edges
  FN_edges <- setdiff(original_edges_str, consensus_edges_str)
  FN_count <- length(FN_edges)
  if (FN_count > 0) {
    FN_edges_matrix <- matrix(unlist(strsplit(FN_edges, "-")), ncol = 2, byrow = TRUE)
    graph_fn <- graph_from_edgelist(FN_edges_matrix, directed = FALSE)
    
    # Plot False Negatives
    set.seed(1234)
    plot_fn <- ggraph(graph_fn, layout = "fr") +
      geom_edge_link(color = "blue", width = 0.5) +  # FN are blue
      geom_node_point(color = "steelblue", size = 0.7) +
      labs(title = paste("False Negatives\nFN:", FN_count)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
  } else {
    plot_fn <- NULL  # No FN to plot
  }
  
  # Combine plots using gridExtra
  if (!is.null(plot_fn)) {
    grid.arrange(plot_comparison, plot_fn, nrow = 1)
  } else {
    grid.arrange(plot_comparison)
  }
}
