compare_consensus <- function(consensus_matrix, original_matrix) {
  
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
    
    set.seed(1234)
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
    # Create graphs for the original and consensus matrices
    graph_original <- graph_from_adjacency_matrix(original_matrix, mode = "undirected", diag = FALSE)
    graph_consensus <- graph_from_adjacency_matrix(consensus_matrix, mode = "undirected", diag = FALSE)
    
    # Ensure all nodes from the original graph are included in the consensus graph
    all_nodes <- union(V(graph_original)$name, V(graph_consensus)$name)
    graph_original <- igraph::induced_subgraph(graph_original, all_nodes)
    graph_consensus <- igraph::induced_subgraph(graph_consensus, all_nodes)
    
    # Explicitly compare edges to identify True Positives (TP) and False Negatives (FN)
    original_edges <- as_edgelist(graph_original)
    consensus_edges <- as_edgelist(graph_consensus)
    
    # Convert edge lists to strings for easy comparison
    original_edges_str <- apply(original_edges, 1, function(x) paste(sort(x), collapse = "-"))
    consensus_edges_str <- apply(consensus_edges, 1, function(x) paste(sort(x), collapse = "-"))
    
    # Edge colors: red for TP (in both), blue for FN (in original but not in consensus)
    edge_colors <- ifelse(original_edges_str %in% consensus_edges_str, "red", "blue")
    
    # Count TP and FN
    TP_count <- sum(edge_colors == "red")
    FN_count <- sum(edge_colors == "blue")
    
    set.seed(1234)
    # Plot the original graph with edges color-coded
    plot(graph_original, 
         main = paste(title, "\nTP:", TP_count, "FN:", FN_count),
         vertex.label = NA,  # Remove labels
         vertex.size = 6, 
         edge.width = 2, 
         edge.color = edge_colors,  # Color edges
         vertex.color = "lightgreen",  # Nodes in light green
         layout = igraph::layout_with_fr)
    
    # Add a legend to explain the edge colors
    legend("topright", 
           legend = c("True Positive (TP)", "False Negative (FN)"), 
           col = c("red", "blue"), 
           lty = 1,  # Line type
           lwd = 2,  # Line width
           cex = 0.8)  # Text size
  }
  
  # Plot the consensus graph
  par(mfrow = c(1, 2))  # Set up for side-by-side plotting
  plot_non_isolated_consensus(consensus_matrix, title = "Final Consensus Graph")
  
  # Plot the original graph with node/edge counts in the title
  plot_original_matrix(original_matrix, title = "Original Matrix Graph")
  
  # Plot comparison of original matrix with consensus highlighted
  par(mfrow = c(1, 1))  # Return to single plot layout
  plot_comparison_graph(original_matrix, consensus_matrix, title = "Original with Consensus Edges")
}

# Example usage:
# Load or define your consensus_matrix and original_matrix here as adjacency matrices.
# consensus_matrix <- ...
# original_matrix <- ...

# Call the function with your matrices
# compare_consensus(consensus_matrix, original_matrix)
