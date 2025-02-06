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
  
  # 1) Build igraphs from adjacency matrices
  graph_original  <- graph_from_adjacency_matrix(original_matrix,  mode="undirected", diag=FALSE)
  graph_consensus <- graph_from_adjacency_matrix(consensus_matrix, mode="undirected", diag=FALSE)
  
  # 2) Extract edges as two-column matrices
  original_edges  <- as_edgelist(graph_original)
  consensus_edges <- as_edgelist(graph_consensus)
  
  # "Stringify" each edge so we can compare sets easily
  original_edges_str  <- apply(original_edges,  1, function(x) paste(sort(x), collapse="-"))
  consensus_edges_str <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse="-"))
  
  # 3) In the original (gold standard) graph:
  #    Mark edges that are also in consensus as TP (red) and those missing as FN (blue)
  E(graph_original)$color <- ifelse(original_edges_str %in% consensus_edges_str, "red", "blue")
  
  # Count TPs and FNs for the original graph:
  TP_count <- sum(E(graph_original)$color == "red")
  FN_count <- sum(E(graph_original)$color == "blue")
  
  # 4) Remove isolated vertices (degree=0) from the original graph before plotting
  graph_original_no_isolates <- delete_vertices(graph_original,
                                                V(graph_original)[degree(graph_original) == 0])
  
  # 5) First plot: Original graph with edges colored red (TP) or blue (FN)
  plot_1 <- ggraph(graph_original_no_isolates, layout="fr") +
    geom_edge_link(aes(color = I(color)), width = 0.7) +
    geom_node_point(color="steelblue", size = 2) +
    #labs(title = paste("Ground Truth\nTP:", TP_count, "FN:", FN_count)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position="none"
    )
  
  # 6) Build the FP graph: edges that are in the consensus (predicted) network but not in the original (gold standard)
  FP_edges_str <- setdiff(consensus_edges_str, original_edges_str)
  FP_count     <- length(FP_edges_str)
  
  # Debugging: Check if FP_edges_str contains valid edges
  print(FP_edges_str)  # Uncomment for debugging
  
  if (FP_count > 0) {
    FP_edges_list <- strsplit(FP_edges_str, "-")
    
    # Ensure all entries have exactly two elements
    FP_edges_list <- Filter(function(x) length(x) == 2, FP_edges_list)
    
    if (length(FP_edges_list) > 0) {
      FP_edges_mat <- do.call(rbind, FP_edges_list)
      graph_fp     <- graph_from_edgelist(FP_edges_mat, directed=FALSE)
      
      # Remove isolated vertices from the FP subgraph
      graph_fp_no_isolates <- delete_vertices(graph_fp, V(graph_fp)[degree(graph_fp) == 0])
      
      # Second plot: all FP edges in purple
      plot_2 <- ggraph(graph_fp_no_isolates, layout="fr") +
        geom_edge_link(color = "purple", width = 1) +
        geom_node_point(color = "steelblue", size = 2) +
        #labs(title = paste("False Positives\nFP:", FP_count)) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust=0.5, size=14, face="bold"),
          legend.position="none"
        )
      
      # 7) Arrange both plots side-by-side
      grid.arrange(plot_1, plot_2, nrow=1)
    } else {
      # No valid FP edges
      grid.arrange(plot_1)
    }
  } else {
    # No FP edges
    grid.arrange(plot_1)
  }
}
