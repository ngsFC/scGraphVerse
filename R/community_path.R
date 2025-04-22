#' Community Detection and Pathway Enrichment
#'
#' Detects communities in a gene network using one or two methods and performs pathway enrichment.
#'
#' @param adj_matrix Square adjacency matrix with gene names as row/col names.
#' @param methods A character vector: either 1 method (just apply it) or 2 methods (compare and choose best). Default = "louvain".
#' @param pathway_db "KEGG" or "Reactome" (default: "KEGG").
#' @param genes_path Minimum community size for enrichment (default: 5).
#' @param plot Logical; if TRUE, plot the community structure (default: TRUE).
#' @param verbose Logical; if TRUE, show progress messages (default: TRUE).
#'
#' @return A list containing community assignments, enrichment results, and the graph object.
#' @importFrom igraph graph_from_adjacency_matrix V degree induced_subgraph vcount ecount
#' @importFrom robin membershipCommunities robinCompare robinAUC
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom ReactomePA enrichPathway
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom ggplot2 aes scale_color_manual labs theme_minimal theme element_text
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
community_path <- function(adj_matrix,
                           methods = "louvain",
                           pathway_db = "KEGG",
                           genes_path = 5,
                           plot = TRUE,
                           verbose = TRUE) {
  
  if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("Input adjacency matrix must be a square matrix.")
  }
  if (is.null(rownames(adj_matrix))) {
    stop("Adjacency matrix must have row names (gene names).")
  }
  
  gene_names <- rownames(adj_matrix)
  graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  igraph::V(graph)$name <- gene_names
  
  if (verbose) message("Detecting communities...")
  
  if (length(methods) == 1) {
    best_method <- methods[1]
    best_communities <- robin::membershipCommunities(graph, method = best_method)
  } else if (length(methods) == 2) {
    res <- tryCatch({
      robin::robinCompare(graph, method1 = methods[1], method2 = methods[2])
    }, error = function(e) stop("robinCompare failed: ", conditionMessage(e)))
    
    auc <- robin::robinAUC(res)
    if (!is.numeric(auc) || length(auc) != 2) stop("Unexpected AUC result from robinCompare.")
    
    if (auc[1] < auc[2]) {
      best_method <- methods[1]
      best_communities <- res$Communities1
    } else {
      best_method <- methods[2]
      best_communities <- res$Communities2
    }
  } else {
    stop("methods must be a character vector of length 1 or 2.")
  }
  
  igraph::V(graph)$community <- as.factor(best_communities)
  
  # Remove isolated nodes for plotting
  non_isolated_nodes <- igraph::degree(graph) > 0
  
  if (plot) {
    plot_graph <- igraph::induced_subgraph(graph, vids = igraph::V(graph)[non_isolated_nodes])
    num_communities <- length(unique(igraph::V(plot_graph)$community))
    
    colors <- if (num_communities <= 12) {
      RColorBrewer::brewer.pal(num_communities, "Set3")
    } else {
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(num_communities)
    }
    
    plot_title <- paste0("Community Structure (", best_method, ")\nNodes: ",
                         igraph::vcount(plot_graph), " Edges: ", igraph::ecount(plot_graph))
    
    g <- ggraph::ggraph(plot_graph, layout = "fr") +
      ggraph::geom_edge_link(color = "gray", width = 0.5) +
      ggraph::geom_node_point(ggplot2::aes(color = community), size = 3) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(title = plot_title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
      )
    print(g)
  }
  
  if (verbose) message("Running pathway enrichment...")
  
  pathway_results <- list()
  non_isolated_genes <- igraph::V(graph)$name[non_isolated_nodes]
  
  for (comm in unique(igraph::V(graph)$community)) {
    genes <- intersect(
      igraph::V(graph)$name[igraph::V(graph)$community == comm],
      non_isolated_genes
    )
    
    if (length(genes) < genes_path) {
      pathway_results[[as.character(comm)]] <- NULL
      next
    }
    
    entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = genes,
                                    column = "ENTREZID",
                                    keytype = "SYMBOL",
                                    multiVals = "first")
    entrez <- na.omit(entrez)
    
    if (length(entrez) >= genes_path) {
      enrich <- tryCatch({
        switch(pathway_db,
               "KEGG" = clusterProfiler::enrichKEGG(gene = entrez, organism = "hsa", keyType = "kegg"),
               "Reactome" = ReactomePA::enrichPathway(gene = entrez, organism = "human"),
               NULL)
      }, error = function(e) {
        warning("Enrichment failed for community ", comm, ": ", conditionMessage(e))
        NULL
      })
      
      if (!is.null(enrich) && nrow(enrich@result) > 0) {
        pathway_results[[as.character(comm)]] <- enrich
      } else {
        pathway_results[[as.character(comm)]] <- NULL
      }
    } else {
      pathway_results[[as.character(comm)]] <- NULL
    }
  }
  
  return(list(
    communities = list(best_method = best_method, membership = best_communities),
    pathways = pathway_results,
    graph = graph
  ))
}
