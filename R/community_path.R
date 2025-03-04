community_path <- function(adj_matrix, methods = NULL, pathway_db = "KEGG", genes_path = 5) {
  required_pkgs <- c("robin", "igraph", "ggraph", "ggplot2", "RColorBrewer", "clusterProfiler", "org.Hs.eg.db", "ReactomePA")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), 
         "\nInstall them using BiocManager::install() or install.packages().")
  }
  
  library(robin)
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
  
  # Create an igraph object directly from the adjacency matrix
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  V(graph)$name <- gene_names
  
  # Community detection
  if (is.null(methods)) {
    # Compare two methods
    robin_result <- robinCompare(graph, method1 = "louvain", method2 = "fastGreedy")
    auc_values <- robinAUC(robin_result)  
    
    if (is.numeric(auc_values) && length(auc_values) == 2) {
      best_method <- ifelse(auc_values[1] < auc_values[2], "louvain", "fastGreedy")
    } else {
      stop("Unexpected output from robinAUC(). Please check the function output.")
    }
    
    community_assignments <- membershipCommunities(graph, method = best_method)
  } else if (length(methods) == 1) {
    best_method <- methods[1]
    community_assignments <- membershipCommunities(graph, method = best_method)
  } else {
    robin_result <- robinCompare(graph, method1 = methods[1], method2 = methods[2])
    auc_values <- robinAUC(robin_result)  
    
    if (is.numeric(auc_values) && length(auc_values) == 2) {
      best_method <- ifelse(auc_values[1] < auc_values[2], methods[1], methods[2])
    } else {
      stop("Unexpected output from robinAUC(). Please check the function output.")
    }
    
    community_assignments <- membershipCommunities(graph, method = best_method)
  }
  
  # Assign communities to nodes
  V(graph)$community <- as.factor(community_assignments)
  
  num_communities <- length(unique(V(graph)$community))
  community_colors <- if (num_communities <= 12) {
    brewer.pal(num_communities, "Set3")
  } else {
    colorRampPalette(brewer.pal(12, "Set3"))(num_communities)
  }
  
  set.seed(1234)
  plot_title <- paste("Community Structure (", best_method, ")\nNodes:", vcount(graph), "Edges:", ecount(graph), sep = "")
  
  p <- ggraph(graph, layout = "fr") + 
    geom_edge_link(color = "gray", width = 0.5) +
    geom_node_point(aes(color = community), size = 3) +
    scale_color_manual(values = community_colors) +
    labs(title = plot_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  print(p)
  
  # Pathway enrichment analysis
  pathway_results <- list()
  
  for (comm in unique(V(graph)$community)) {
    genes_in_community <- V(graph)$name[V(graph)$community == comm]
    
    # Check if community meets minimum size requirement
    if (length(genes_in_community) < genes_path) {
      pathway_results[[as.character(comm)]] <- NULL  # Skip small communities
      next
    }
    
    entrez_ids <- mapIds(org.Hs.eg.db, keys = genes_in_community, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    entrez_ids <- na.omit(entrez_ids)  # Remove missing values
    
    if (length(entrez_ids) >= genes_path) {  # Ensure sufficient genes for pathway analysis
      if (pathway_db == "KEGG") {
        tryCatch({
          enrich_res <- enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "kegg")
          if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
            pathway_results[[as.character(comm)]] <- enrich_res
          } else {
            pathway_results[[as.character(comm)]] <- NULL
          }
        }, error = function(e) {
          warning("KEGG enrichment failed for community ", comm, ": ", conditionMessage(e))
          pathway_results[[as.character(comm)]] <- NULL
        })
      } else if (pathway_db == "Reactome") {
        tryCatch({
          enrich_res <- enrichPathway(gene = entrez_ids, organism = "human")
          if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
            pathway_results[[as.character(comm)]] <- enrich_res
          } else {
            pathway_results[[as.character(comm)]] <- NULL
          }
        }, error = function(e) {
          warning("Reactome enrichment failed for community ", comm, ": ", conditionMessage(e))
          pathway_results[[as.character(comm)]] <- NULL
        })
      } else {
        warning("Invalid pathway_db argument. Use 'KEGG' or 'Reactome'. Skipping pathway analysis.")
        pathway_results[[as.character(comm)]] <- NULL
      }
    } else {
      pathway_results[[as.character(comm)]] <- NULL
    }
  }
  
  return(list(communities = list(best_method = best_method, membership = community_assignments), pathways = pathway_results))
}
