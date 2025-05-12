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
#' genes <- c(
#'     "TP53",   "BRCA1", "EGFR",   "MYC",     "CDKN1A",
#'     "BCL2",   "MDM2",  "PTEN",   "AKT1",    "MAPK1"
#' )
#'
#' adj <- matrix(
#'     0,
#'     nrow = length(genes),
#'     ncol = length(genes),
#'     dimnames = list(genes, genes)
#' )
#'
#' edge_list <- list(
#'     c("TP53", "MDM2"),
#'     c("TP53", "CDKN1A"),
#'     c("BRCA1", "BCL2"),
#'     c("PTEN", "AKT1"),
#'     c("EGFR", "MAPK1"),
#'     c("MYC", "CDKN1A"),
#'     c("MYC", "BCL2"),
#'     c("AKT1", "MAPK1")
#' )
#'
#' for (e in edge_list) {
#'     adj[e[1], e[2]] <- 1
#'     adj[e[2], e[1]] <- 1
#' }
#'
#' diag(adj) <- 0
#' result <- community_path(
#'     adj_matrix = adj,
#'     methods    = "louvain",
#'     pathway_db = "KEGG",
#'     genes_path = 2,
#'     plot       = TRUE,
#'     verbose    = TRUE
#' )
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
