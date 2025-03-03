community_path <- function(adj_matrix, methods = c("louvain"), pathway_db = "KEGG") {
  required_pkgs <- c("igraph", "ggraph", "ggplot2", "RColorBrewer", "clusterProfiler", "org.Hs.eg.db", "ReactomePA")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), 
         "\nInstall them using BiocManager::install() or install.packages().")
  }
  
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  
  
  if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("Error: Input adjacency matrix must be square.")
  }
  
  if (is.null(rownames(adj_matrix))) {
    stop("Error: Adjacency matrix must have row names corresponding to gene names.")
  }
  gene_names <- rownames(adj_matrix)
  
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  V(graph)$name <- gene_names
  
  community_results <- list()
  
  if ("louvain" %in% methods) {
    community_results$louvain <- cluster_louvain(graph)$membership
  }
  if ("walktrap" %in% methods) {
    community_results$walktrap <- cluster_walktrap(graph)$membership
  }
  if ("infomap" %in% methods) {
    community_results$infomap <- cluster_infomap(graph)$membership
  }
  
  method_to_plot <- methods[1]
  V(graph)$community <- as.factor(community_results[[method_to_plot]])
  
  num_communities <- length(unique(V(graph)$community))
  community_colors <- if (num_communities <= 12) {
    brewer.pal(num_communities, "Set3")
  } else {
    colorRampPalette(brewer.pal(12, "Set3"))(num_communities)
  }
  
  set.seed(1234)
  plot_title <- paste("Community Structure (", method_to_plot, ")\nNodes:", vcount(graph), "Edges:", ecount(graph), sep = "")
  
  p <- ggraph(graph, layout = "fr") + 
    geom_edge_link(color = "gray", width = 0.5) +
    geom_node_point(aes(color = community), size = 3) +
    scale_color_manual(values = community_colors) +
    labs(title = plot_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  print(p)
  
  pathway_results <- list()
  
  for (comm in unique(V(graph)$community)) {
    genes_in_community <- V(graph)$name[V(graph)$community == comm]
    
    entrez_ids <- mapIds(org.Hs.eg.db, keys = genes_in_community, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    entrez_ids <- na.omit(entrez_ids)  # Remove missing values
    
    if (length(entrez_ids) >= 10) {  # Require minimum genes for enrichment
      if (pathway_db == "KEGG") {
        enrich_res <- enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "kegg")
      } else if (pathway_db == "Reactome") {
        enrich_res <- enrichPathway(gene = entrez_ids, organism = "human")
      } else {
        warning("Invalid pathway_db argument. Use 'KEGG' or 'Reactome'. Skipping pathway analysis.")
        enrich_res <- NULL
      }
      
      if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
        pathway_results[[as.character(comm)]] <- enrich_res
      } else {
        pathway_results[[as.character(comm)]] <- NULL
      }
    } else {
      pathway_results[[as.character(comm)]] <- NULL
    }
  }
  
  return(list(communities = community_results, pathways = pathway_results))
}

