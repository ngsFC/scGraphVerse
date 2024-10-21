plotg <- function(adj_matrix_list) {
  # Load igraph package (if not already loaded)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required but not installed.")
  }
  library(igraph)
  
  # Determine the number of matrices
  num_matrices <- length(adj_matrix_list)
  
  # Set up dynamic layout based on the number of matrices
  num_cols <- ceiling(sqrt(num_matrices))  # Number of columns
  num_rows <- ceiling(num_matrices / num_cols)  # Number of rows
  par(mfrow = c(num_rows, num_cols))
  
  # Initialize an empty list to store graph metrics data frames
  graph_metrics_list <- list()
  
  for (i in seq_along(adj_matrix_list)) {
    # Create graph from the adjacency matrix
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix_list[[i]], mode = "undirected", diag = FALSE)
    
    # Remove isolated vertices (nodes with no connections)
    non_isolated_vertices <- igraph::V(graph)[igraph::degree(graph) > 0]
    subgraph <- igraph::induced_subgraph(graph, non_isolated_vertices)
    
    # Calculate metrics for plot title
    num_nodes <- igraph::gorder(subgraph)  # Number of nodes
    num_edges <- igraph::gsize(subgraph)  # Number of edges
    
    # Plot the subgraph with customized settings
    plot(subgraph, 
         main = paste("Graph", i, "
Nodes:", num_nodes, "Edges:", num_edges),
         vertex.label = NA,  # Remove labels
         vertex.size = 6, 
         edge.width = 2, 
         vertex.color = "orange",  # Change node color
         layout = igraph::layout_with_fr)
    
    # Calculate and store graph metrics
    avg_degree <- mean(igraph::degree(subgraph))  # Average degree
    density <- igraph::graph.density(subgraph)  # Graph density
    diameter <- if (num_nodes > 1) igraph::diameter(subgraph) else NA  # Diameter (NA if only 1 node)
    avg_path_length <- if (num_nodes > 1) igraph::mean_distance(subgraph) else NA  # Average path length (NA if only 1 node)
    num_components <- igraph::components(subgraph)$no  # Number of connected components
    
    # Create a data frame to hold the metrics for this graph
    graph_metrics <- data.frame(
      Matrix_Index = i,
      Num_Nodes = num_nodes,
      Num_Edges = num_edges,
      Avg_Degree = avg_degree,
      Density = density,
      Diameter = diameter,
      Avg_Path_Length = avg_path_length,
      Num_Components = num_components,
      stringsAsFactors = FALSE
    )
    
    # Append the data frame to the list
    graph_metrics_list[[i]] <- graph_metrics
  }
  
  # Reset plotting area to default
  par(mfrow = c(1, 1))
  
  # Combine all data frames into a single data frame for easier viewing
  combined_metrics_df <- do.call(rbind, graph_metrics_list)
  
  # Return the combined data frame and individual graph metrics
  return(list("Graph_Metrics_List" = graph_metrics_list, 
              "Combined_Metrics" = combined_metrics_df))
}
