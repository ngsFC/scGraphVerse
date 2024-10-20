compare_consensus <- function(adj_matrix_list, original_matrix) {
  
  create_consensus_matrix <- function(adj_matrix_list) {
    # Sum all adjacency matrices to get consensus values
    consensus_matrix <- Reduce("+", adj_matrix_list)
    # Define threshold (present in at least half of the adjacency matrices)
    threshold <- (length(adj_matrix_list)*3) / 4
    # Create binary consensus matrix (1 if present in >= threshold matrices)
    consensus_matrix_binary <- consensus_matrix >= threshold
    #diag(consensus_matrix_binary) <- 1
    return(consensus_matrix_binary)
  }
  
  plot_non_isolated_consensus <- function(consensus_matrix, title = "Consensus Graph") {
    graph <- graph_from_adjacency_matrix(consensus_matrix, mode = "undirected", diag = FALSE)
    non_isolated_vertices <- V(graph)[degree(graph) > 0]
    subgraph <- induced_subgraph(graph, non_isolated_vertices)
    plot(subgraph, 
         main = title, 
         vertex.label.color = "black",
         vertex.size = 5, 
         edge.width = 2, 
         vertex.label.cex = 0.8,
         layout = layout_with_fr)
  }
  
  plot_original_with_highlight <- function(original_matrix, consensus_matrix) {
    graph_original <- graph_from_adjacency_matrix(original_matrix, mode = "undirected", diag = FALSE)
    graph_consensus <- graph_from_adjacency_matrix(consensus_matrix, mode = "undirected", diag = FALSE)
    
    original_edges <- as_edgelist(graph_original)
    consensus_edges <- as_edgelist(graph_consensus)
    
    # Create sets of edges for comparison
    original_edges_set <- apply(original_edges, 1, function(x) paste(sort(x), collapse = "-"))
    consensus_edges_set <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse = "-"))
    
    # Color edges red if they are in the consensus matrix, blue otherwise
    edge_colors <- ifelse(original_edges_set %in% consensus_edges_set, "red", "blue")
    
    # Plot the original graph with highlighted edges
    plot(graph_original, 
         edge.color = edge_colors, 
         main = "Original Graph with Consensus Highlight",
         vertex.label.color = "black",
         vertex.size = 5, 
         edge.width = 2, 
         vertex.label.cex = 0.8,
         layout = layout_with_fr)
  }

  consensus_matrix <- create_consensus_matrix(adj_matrix_list)
  
  par(mfrow = c(1, 2)) 
  plot_non_isolated_consensus(consensus_matrix, title = "Consensus Graph")
  plot_non_isolated_consensus(original_matrix, title = "Original Graph")
  par(mfrow = c(1, 1))  # Reset plotting area
  
  plot_original_with_highlight(original_matrix, consensus_matrix)
  
  return(list("Consensus_Matrix" = consensus_matrix))
}

