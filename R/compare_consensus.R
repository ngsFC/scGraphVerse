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
  
  # Helper function
  create_graph <- function(adj_matrix) {
    graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  }
  
  # Make igraph objects
  graph_original  <- create_graph(original_matrix)
  graph_consensus <- create_graph(consensus_matrix)
  
  # Edge lists
  original_edges  <- as_edgelist(graph_original)
  consensus_edges <- as_edgelist(graph_consensus)
  
  # "Stringify" edges for set comparison
  original_edges_str  <- apply(original_edges,  1, function(x) paste(sort(x), collapse = "-"))
  consensus_edges_str <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse = "-"))
  
  # Color in the original graph:
  #   red = TP (in both)
  #   blue = FP (in original only)
  edge_colors <- ifelse(original_edges_str %in% consensus_edges_str, "red", "blue")
  
  # Count TPs and FPs
  TP_count <- sum(edge_colors == "red")
  FP_count <- sum(edge_colors == "blue")
  
  # Plot 1: Original graph (TP=red, FP=blue)
  plot_1 <- ggraph(graph_original, layout = "fr") +
    geom_edge_link(aes(color = edge_colors), width = 0.5) +
    geom_node_point(color = "steelblue", size = 0.7) +
    scale_color_manual(values = c("red" = "red", "blue" = "blue")) +
    labs(title = paste("Ground Truth & Consensus\nTP:", TP_count, " FP:", FP_count)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    )
  
  # FN edges: in consensus only
  FN_edges_str <- setdiff(consensus_edges_str, original_edges_str)
  FN_count     <- length(FN_edges_str)
  
  if (FN_count > 0) {
    # NO as.numeric() here!
    FN_edges_matrix <- do.call(
      rbind,
      lapply(strsplit(FN_edges_str, "-"), function(x) x)
    )
    graph_fn <- graph_from_edgelist(FN_edges_matrix, directed = FALSE)
    
    # Plot 2: FN edges in purple
    plot_2 <- ggraph(graph_fn, layout = "fr") +
      geom_edge_link(color = "purple", width = 0.5) +
      geom_node_point(color = "steelblue", size = 0.7) +
      labs(title = paste("False Negatives\nFN:", FN_count)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
    
    grid.arrange(plot_1, plot_2, nrow = 1)
  } else {
    # Just show the first plot if no FN
    grid.arrange(plot_1)
  }
}
