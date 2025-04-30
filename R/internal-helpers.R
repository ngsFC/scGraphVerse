.detect_communities <- function(graph, methods) {
  if (length(methods) == 1) {
    best_method <- methods[1]
    best_communities <- robin::membershipCommunities(graph, method = best_method)
  } else if (length(methods) == 2) {
    res <- tryCatch(
      robin::robinCompare(graph, method1 = methods[1], method2 = methods[2]),
      error = function(e) stop("robinCompare failed: ", conditionMessage(e))
    )

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

  return(list(best_method = best_method, best_communities = best_communities))
}

.plot_communities <- function(graph, best_method) {
  non_isolated_nodes <- igraph::degree(graph) > 0
  plot_graph <- igraph::induced_subgraph(graph, vids = igraph::V(graph)[non_isolated_nodes])
  num_communities <- length(unique(igraph::V(plot_graph)$community))

  colors <- if (num_communities <= 12) {
    RColorBrewer::brewer.pal(num_communities, "Set3")
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(num_communities)
  }

  plot_title <- paste0(
    "Community Structure (", best_method, ")\nNodes: ",
    igraph::vcount(plot_graph), " Edges: ", igraph::ecount(plot_graph)
  )

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

  plot(g)
}

.enrich_communities <- function(graph, non_isolated_nodes, pathway_db, genes_path) {
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
      multiVals = "first"
    )
    entrez <- na.omit(entrez)

    if (length(entrez) >= genes_path) {
      enrich <- tryCatch(
        {
          switch(pathway_db,
            "KEGG" = clusterProfiler::enrichKEGG(gene = entrez, organism = "hsa", keyType = "kegg"),
            "Reactome" = ReactomePA::enrichPathway(gene = entrez, organism = "human"),
            NULL
          )
        },
        error = function(e) {
          warning("Enrichment failed for community ", comm, ": ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(enrich) && nrow(enrich@result) > 0) {
        pathway_results[[as.character(comm)]] <- enrich
      } else {
        pathway_results[[as.character(comm)]] <- NULL
      }
    } else {
      pathway_results[[as.character(comm)]] <- NULL
    }
  }

  return(pathway_results)
}

