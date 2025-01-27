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
#'         - One for true positives and false negatives.
#'         - Another for false positives.
#'
#' @details
#' The function compares two graphs: the consensus graph and the ground truth graph. It computes:
#' - True Positives (TP): Edges that are present in both graphs.
#' - False Negatives (FN): Edges that are present in the original graph but not in the consensus graph.
#' - False Positives (FP): Edges that are present in the consensus graph but not in the original graph.
#' 
#' The true positives (TP) and false negatives (FN) are plotted together, while the false positives (FP) are plotted separately.
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
  plot_ggraph <- function(graph, plot_title, edge_colors = NULL) {
    ggraph(graph, layout = "fr") +
      geom_edge_link(aes(color = edge_colors), width = 0.5) +
      geom_node_point(color = "steelblue", size = 0.7) +
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
  
  # Edge comparison
  original_edges <- as_edgelist(graph_original)
  consensus_edges <- as_edgelist(graph_consensus)
  
  original_edges_str <- apply(original_edges, 1, function(x) paste(sort(x), collapse = "-"))
  consensus_edges_str <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse = "-"))
  
  edge_colors_ground_truth <- ifelse(original_edges_str %in% consensus_edges_str, "red", "blue")
  TP_count <- sum(edge_colors_ground_truth == "red")
  FN_count <- sum(edge_colors_ground_truth == "blue")
  
  # Plot ground truth with TP and FN
  set.seed(1234)
  plot_ground_truth <- plot_ggraph(graph_original, plot_title = paste("Ground Truth\nTP:", TP_count, "FN:", FN_count), edge_colors = edge_colors_ground_truth)
  
  # False Positive Graph
  FP_edges <- setdiff(consensus_edges_str, original_edges_str)
  FP_count <- length(FP_edges)
  FP_edges_matrix <- matrix(unlist(strsplit(FP_edges, "-")), ncol = 2, byrow = TRUE)
  graph_fp <- graph_from_edgelist(FP_edges_matrix, directed = FALSE)
  
  # Add purple color to FP edges
  edge_colors_fp <- rep("purple", nrow(FP_edges_matrix))
  set.seed(1234)
  plot_fp <- ggraph(graph_fp, layout = "fr") +
    geom_edge_link(color = "purple", width = 0.5) +
    geom_node_point(color = "steelblue", size = 0.7) +
    labs(title = paste("False Positives\nFP:", FP_count)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    )
  
  # Combine plots using gridExtra
  grid.arrange(plot_ground_truth, plot_fp, nrow = 1)
}

