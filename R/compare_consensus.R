#' Compare Consensus and Original Graphs (or STRINGdb Network)
#'
#' This function compares a consensus adjacency matrix with either:
#' - A provided reference adjacency matrix.
#' - A STRINGdb-derived adjacency matrix if no reference is provided.
#'
#' It generates plots visualizing edges:
#' - **True Positives (TP or CE)**: Edges that are 1 in both.
#' - **False Negatives (FN or ME)**: Edges that are 1 in the reference but 0 in the consensus.
#' - **False Positives (FP or EE)** (optional): Edges that are 1 in the consensus but 0 in the reference.
#'
#' @param consensus_matrix A binary adjacency matrix representing the consensus graph.
#' @param reference_matrix (Optional) A binary adjacency matrix as a reference (ground truth). If `NULL`, STRINGdb is used.
#' @param false_plot Logical; if TRUE, plot False Positives (FP) (default = FALSE).
#'
#' @return A ggplot object (or a combined plot if `false_plot = TRUE`).
#'
#' @importFrom igraph V graph_from_adjacency_matrix as_edgelist graph_from_edgelist delete_vertices
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix Matrix
#' @importFrom STRINGdb STRINGdb
#' @export
#'
#' @examples
#' set.seed(42)
#' original <- matrix(sample(0:1, 25, replace = TRUE, prob = c(0.8, 0.2)), 5, 5)
#' consensus <- matrix(sample(0:1, 25, replace = TRUE, prob = c(0.8, 0.2)), 5, 5)
#' diag(original) <- diag(consensus) <- 0
#' compare_consensus(consensus, reference_matrix = original, false_plot = TRUE)
#'
#' # Using STRINGdb
#' compare_consensus(consensus)

compare_consensus <- function(consensus_matrix, reference_matrix = NULL, false_plot = FALSE) {
  
  # Validate consensus_matrix
  if (!is.matrix(consensus_matrix)) {
    stop("consensus_matrix must be a binary adjacency matrix.")
  }
  
  # If reference_matrix is NULL, use STRINGdb to get adjacency matrix
  use_STRINGdb <- is.null(reference_matrix)
  
  if (use_STRINGdb) {
    if (is.null(rownames(consensus_matrix))) {
      stop("consensus_matrix must have row names to query STRINGdb.")
    }
    
    adj <- stringdb_adjacency(
      genes = rownames(consensus_matrix),
      species = 9606,
      required_score = 900,
      keep_all_genes = TRUE
    )$binary
    
    adj <- adj[rownames(consensus_matrix), rownames(consensus_matrix)]
    adj <- symmetrize(list(adj), weight_function = "mean")[[1]]
    
    reference_matrix <- adj  # Use STRINGdb adjacency as reference
  }
  
  # Validate reference_matrix
  if (!is.matrix(reference_matrix)) {
    stop("reference_matrix must be a binary adjacency matrix.")
  }
  
  if (!identical(dim(consensus_matrix), dim(reference_matrix))) {
    stop("Both matrices must have the same dimensions.")
  }
  
  # Create graph objects
  graph_reference <- igraph::graph_from_adjacency_matrix(reference_matrix, mode = "undirected", diag = FALSE)
  graph_consensus <- igraph::graph_from_adjacency_matrix(consensus_matrix, mode = "undirected", diag = FALSE)
  
  # Extract edges as sorted strings
  edge_to_string <- function(edge_list) {
    apply(edge_list, 1, function(x) paste(sort(x), collapse = "-"))
  }
  
  reference_edges_str <- edge_to_string(as_edgelist(graph_reference))
  consensus_edges_str <- edge_to_string(as_edgelist(graph_consensus))
  
  # Adjust labels based on data source
  if (use_STRINGdb) {
    TP_label <- "CE (Confirmed Edges)"
    FN_label <- "ME (Missing Edges)"
    FP_label <- "EE (Extra Edges)"
  } else {
    TP_label <- "TP (True Positives)"
    FN_label <- "FN (False Negatives)"
    FP_label <- "FP (False Positives)"
  }
  
  # Identify TP and FN edges
  edge_colors <- ifelse(reference_edges_str %in% consensus_edges_str, "red", "blue")
  
  # Remove isolated nodes
  graph_reference_no_isolates <- delete_vertices(graph_reference, V(graph_reference)[degree(graph_reference) == 0])
  
  # Plot reference graph with TP and FN edges
  TP_count <- sum(edge_colors == "red")
  FN_count <- sum(edge_colors == "blue")
  
  set.seed(1234)
  plot_tp_fn <- ggraph(graph_reference_no_isolates, layout = "fr") +
    geom_edge_link(aes(color = I(edge_colors)), width = 0.7) +
    geom_node_point(color = "steelblue", size = 1.5) +
    labs(title = paste("Reference Graph\n", TP_label, ":", TP_count, FN_label, ":", FN_count)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Handle FP plot conditionally
  if (false_plot) {
    FP_edges_str <- setdiff(consensus_edges_str, reference_edges_str)
    FP_count <- length(FP_edges_str)
    
    if (FP_count > 0) {
      FP_edges_list <- strsplit(FP_edges_str, "-")
      FP_edges_list <- Filter(function(x) length(x) == 2, FP_edges_list)
      
      if (length(FP_edges_list) > 0) {
        FP_edges_mat <- do.call(rbind, FP_edges_list)
        graph_fp <- igraph::graph_from_edgelist(FP_edges_mat, directed = FALSE)
        
        # Remove isolated nodes
        graph_fp_no_isolates <- delete_vertices(graph_fp, V(graph_fp)[degree(graph_fp) == 0])
        
        # Plot FP edges
        plot_fp <- ggraph(graph_fp_no_isolates, layout = "fr") +
          geom_edge_link(color = "purple", width = 1) +
          geom_node_point(color = "steelblue", size = 2) +
          labs(title = paste(FP_label, ":", FP_count)) +
          theme_minimal() +
          theme(legend.position = "none")
        
        # Combine plots
        return(wrap_plots(plot_tp_fn, plot_fp, nrow = 1))
      }
    }
  }
  
  # Return single plot if FP plot not needed
  return(plot_tp_fn)
}
