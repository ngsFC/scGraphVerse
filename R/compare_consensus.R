#' Compare Consensus and Original Graphs
#'
#' This function compares a consensus adjacency matrix with an original adjacency matrix, visualizing edges as:
#' - True Positives (TP): edges that are 1 in both (red).
#' - False Positives (FP): edges that are 1 in the original but 0 in the consensus (blue).
#' - False Negatives (FN): edges that are 0 in the original but 1 in the consensus (purple).
#'
#' It generates two plots:
#' 1. A plot of the original (ground truth) graph with TP edges in red and FP edges in blue.
#' 2. A plot showing only the FN edges (purple).
#'
#' @param consensus_matrix A binary adjacency matrix representing the consensus graph.
#' @param original_matrix A binary adjacency matrix representing the original or ground truth graph.
#'
#' @return This function does not return any values. It generates two plots:
#'         - The first for True Positives (red) and False Positives (blue).
#'         - The second for False Negatives (purple).
#'
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist graph_from_edgelist
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom gridExtra grid.arrange
#'
#' @export
compare_consensus <- function(consensus_matrix, original_matrix) {
  
  # Helper function to create an igraph from an adjacency matrix
  create_graph <- function(adj_matrix) {
    graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  }
  
  # Helper function for a ggraph plot
  plot_ggraph <- function(graph, plot_title, edge_colors) {
    ggraph(graph, layout = "fr") +
      geom_edge_link(aes(color = edge_colors), width = 0.5) +
      geom_node_point(color = "steelblue", size = 0.7) +
      scale_color_manual(values = c("red" = "red", "blue" = "blue")) +
      labs(title = plot_title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
  }
  
  # Create igraph objects
  graph_original  <- create_graph(original_matrix)
  graph_consensus <- create_graph(consensus_matrix)
  
  # Extract edge lists (each row is a pair of vertices)
  original_edges  <- as_edgelist(graph_original)
  consensus_edges <- as_edgelist(graph_consensus)
  
  # Convert each edge to a "string" for easy set comparisons
  original_edges_str  <- apply(original_edges,  1, function(x) paste(sort(x), collapse = "-"))
  consensus_edges_str <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse = "-"))
  
  # Color edges of the original graph:
  #   red  = TP (edge also in consensus)
  #   blue = FP (edge not in consensus)
  edge_colors <- ifelse(original_edges_str %in% consensus_edges_str, "red", "blue")
  
  # Count TPs and FPs
  TP_count <- sum(edge_colors == "red")
  FP_count <- sum(edge_colors == "blue")
  
  # First plot: original graph with TP (red) and FP (blue)
  set.seed(1234)
  plot_comparison <- plot_ggraph(
    graph_original,
    plot_title = paste("Ground Truth & Consensus\nTP:", TP_count, " FP:", FP_count),
    edge_colors = edge_colors
  )
  
  # Identify edges that are in consensus but not in original => False Negatives (FN)
  FN_edges_str <- setdiff(consensus_edges_str, original_edges_str)
  FN_count     <- length(FN_edges_str)
  
  # Build a separate graph from those FN edges
  if (FN_count > 0) {
    # Turn the "a-b" strings back into a two-column edge list
    FN_edges_matrix <- do.call(
      rbind,
      lapply(strsplit(FN_edges_str, "-"), function(x) as.numeric(x))
    )
    graph_fn <- graph_from_edgelist(FN_edges_matrix, directed = FALSE)
    
    # Plot these FN edges (purple)
    set.seed(1234)
    plot_fn <- ggraph(graph_fn, layout = "fr") +
      geom_edge_link(color = "purple", width = 0.5) +
      geom_node_point(color = "steelblue", size = 0.7) +
      labs(title = paste("False Negatives\nFN:", FN_count)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
  } else {
    plot_fn <- NULL
  }
  
  # Combine the two plots
  if (!is.null(plot_fn)) {
    grid.arrange(plot_comparison, plot_fn, nrow = 1)
  } else {
    grid.arrange(plot_comparison)
  }
}
