#' Community Detection and Pathway Enrichment Analysis
#'
#' Detects gene communities within an adjacency network using one or two
#' community detection methods, and performs pathway enrichment for each
#' detected community.
#'
#' @param adj_matrix A square adjacency matrix. Row and column names must
#'   correspond to gene symbols.
#' @param methods A character vector of one or two community detection
#'   methods supported by \pkg{robin}. If two are given, performance is
#'   compared and the best is selected. Default: \code{"louvain"}.
#' @param pathway_db Character string specifying the pathway database to use:
#'   \code{"KEGG"} or \code{"Reactome"}. Default: \code{"KEGG"}.
#' @param genes_path Integer. Minimum number of genes per community to run
#'   enrichment analysis. Default: \code{5}.
#' @param plot Logical. If \code{TRUE}, generates a plot of detected
#'   communities. Default: \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, shows progress messages.
#'   Default: \code{TRUE}.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{communities}: List with \code{best_method} and a named
#'       vector of community membership per gene.
#'     \item \code{pathways}: List of enrichment results per community
#'       (only for communities meeting size threshold).
#'     \item \code{graph}: The \pkg{igraph} object with community
#'       annotations.
#'   }
#'
#' @details If two methods are provided, the function uses
#'   \code{robinCompare} and selects the method with higher AUC. Pathway
#'   enrichment is done via \pkg{clusterProfiler} (KEGG) or via
#'   \pkg{ReactomePA} (Reactome). Communities smaller than
#'   \code{genes_path} are excluded.
#'
#' @importFrom igraph graph_from_adjacency_matrix V degree
#'   induced_subgraph vcount ecount
#' @importFrom robin membershipCommunities robinCompare robinAUC
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom ReactomePA enrichPathway
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom ggplot2 aes scale_color_manual labs theme_minimal
#'   theme element_text
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' data(count_matrices)
#' data(adj_truth)
#' networks <- infer_networks(
#'     count_matrices_list = count_matrices,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#' head(networks[[1]])
#'
#' wadj_list <- generate_adjacency(networks)
#' swadj_list <- symmetrize(wadj_list, weight_function = "mean")
#'
#' binary_listj <- cutoff_adjacency(
#'     count_matrices = count_matrices,
#'     weighted_adjm_list = swadj_list,
#'     n = 2,
#'     method = "GENIE3",
#'     quantile_threshold = 0.99,
#'     nCores = 1,
#'     debug = TRUE
#' )
#' head(binary_listj[[1]])
#'
#' consensus <- create_consensus(binary_listj, method = "vote")
#' comm_cons <- community_path(consensus)
#'
community_path <- function(
    adj_matrix,
    methods = "louvain",
    pathway_db = "KEGG",
    genes_path = 5,
    plot = TRUE,
    verbose = TRUE) {
    if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
        stop("adj_matrix must be a square matrix.")
    }
    if (is.null(rownames(adj_matrix))) {
        stop("adj_matrix must have row names (gene symbols).")
    }

    gene_names <- rownames(adj_matrix)
    graph <- igraph::graph_from_adjacency_matrix(
        adj_matrix,
        mode = "undirected", diag = FALSE
    )
    igraph::V(graph)$name <- gene_names

    if (verbose) message("Detecting communities...")
    comm_res <- .detect_communities(graph, methods)
    best_method <- comm_res$best_method
    best_communities <- comm_res$best_communities
    igraph::V(graph)$community <- as.factor(best_communities)

    non_iso <- igraph::degree(graph) > 0
    if (plot) {
        .plot_communities(graph, best_method)
    }

    if (verbose) message("Running pathway enrichment...")
    pathway_results <- .enrich_communities(
        graph, non_iso, pathway_db, genes_path
    )

    list(
        communities = list(
            best_method = best_method,
            membership  = best_communities
        ),
        pathways = pathway_results,
        graph = graph
    )
}
