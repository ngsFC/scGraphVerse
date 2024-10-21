compare_consensus <- function(adj_matrix_list, original_matrix) {
  
  create_consensus_matrix <- function(adj_matrix_list) {
    # Sum all adjacency matrices to get consensus values
    consensus_matrix <- Reduce("+", adj_matrix_list)
    # Define threshold (present in at least 75% of the adjacency matrices)
    threshold <- round(length(adj_matrix_list) * 0.75)
    # Create binary consensus matrix (1 if present in >= threshold matrices)
    consensus_matrix_binary <- consensus_matrix >= threshold
    return(consensus_matrix_binary)
  }
  
  plot_non_isolated_consensus <- function(consensus_matrix, title = "Consensus Graph") {
    graph <- graph_from_adjacency_matrix(consensus_matrix, mode = "undirected", diag = FALSE)
    non_isolated_vertices <- V(graph)[degree(graph) > 0]
    subgraph <- induced_subgraph(graph, non_isolated_vertices)
    
    # Calculate metrics for plot title
    num_nodes <- igraph::gorder(subgraph)  # Number of nodes
    num_edges <- igraph::gsize(subgraph)  # Number of edges
    
    # Plot the subgraph with customized settings
    plot(subgraph, 
         main = paste(title, "\nNodes:", num_nodes, "Edges:", num_edges),
         vertex.label = NA,  # Remove labels
         vertex.size = 6, 
         edge.width = 2, 
         vertex.color = "orange",  # Change node color
         layout = igraph::layout_with_fr)
  }
  
  plot_original_matrix <- function(original_matrix, title = "Original Graph") {
    graph <- graph_from_adjacency_matrix(original_matrix, mode = "undirected", diag = FALSE)
    non_isolated_vertices <- V(graph)[degree(graph) > 0]
    subgraph <- induced_subgraph(graph, non_isolated_vertices)
    
    # Calculate metrics for plot title
    num_nodes <- igraph::gorder(subgraph)  # Number of nodes
    num_edges <- igraph::gsize(subgraph)  # Number of edges
    
    # Plot the subgraph for the original matrix
    plot(subgraph, 
         main = paste(title, "\nNodes:", num_nodes, "Edges:", num_edges),
         vertex.label = NA,  # Remove labels
         vertex.size = 6, 
         edge.width = 2, 
         vertex.color = "lightblue",  # Change node color
         layout = igraph::layout_with_fr)
  }
  
  plot_comparison_graph <- function(original_matrix, consensus_matrix, title = "Original vs Consensus Edges") {
    graph_original <- graph_from_adjacency_matrix(original_matrix, mode = "undirected", diag = FALSE)
    graph_consensus <- graph_from_adjacency_matrix(consensus_matrix, mode = "undirected", diag = FALSE)
    
    # Create edge colors based on whether they are in the consensus or not
    edge_colors <- ifelse(E(graph_original) %in% E(graph_consensus), "red", "blue")
    
    # Plot the original graph with edges colored based on consensus
    plot(graph_original, 
         main = title,
         vertex.label = NA,  # Remove labels
         vertex.size = 6, 
         edge.width = 2, 
         edge.color = edge_colors,  # Color edges
         vertex.color = "lightgreen",  # Original nodes in light green
         layout = igraph::layout_with_fr)
  }

  # Create the consensus matrix from the list of adjacency matrices
  consensus_matrix <- create_consensus_matrix(adj_matrix_list)
  
  # Plot the consensus graph
  par(mfrow = c(1, 2))  # Set up for side by side plotting
  plot_non_isolated_consensus(consensus_matrix, title = "Final Consensus Graph")
  
  # Plot the original graph with node/edge counts in the title
  plot_original_matrix(original_matrix, title = "Original Matrix Graph")
  
  # Plot comparison of original matrix with consensus highlighted
  par(mfrow = c(1, 1))  # Return to single plot layout
  plot_comparison_graph(original_matrix, consensus_matrix, title = "Original with Consensus Edges")
  
  return(consensus_matrix)
}

