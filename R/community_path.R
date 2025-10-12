#' Community Detection and Pathway Enrichment Analysis
#'
#' Detects gene communities within an adjacency network using one or two
#' community detection methods, and performs pathway enrichment for each
#' detected community.
#'
#' @param adj_matrix A square adjacency matrix or a
#'   \linkS4class{SummarizedExperiment} object containing a single adjacency
#'   matrix. Row and column names must correspond to gene symbols.
#' @param methods A character vector of one or two community detection
#'   methods supported by \pkg{robin}. If two are given, performance is
#'   compared and the best is selected. Default: \code{"louvain"}.
#' @param pathway_db Character string specifying the pathway database to use:
#'   \code{"KEGG"} or \code{"Reactome"}. Default: \code{"KEGG"}.
#' @param organism Character string specifying the organism: \code{"human"} or
#'   \code{"mouse"}. Default: \code{"human"}.
#' @param genes_path Integer. Minimum number of genes per community to run
#'   enrichment analysis. Default: \code{5}.
#' @param plot Logical. If \code{TRUE}, generates a plot of detected
#'   communities. Default: \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, shows progress messages.
#'   Default: \code{TRUE}.
#' @param method_params List of parameters for community detection methods.
#'   Common parameters include:
#'   \itemize{
#'     \item \code{resolution}: Resolution parameter for Louvain/Leiden
#'     (default: 1)
#'     \item \code{steps}: Number of steps for Walktrap (default: 4)
#'     \item \code{spins}: Number of spins for Spinglass (default: 25)
#'     \item \code{nb.trials}: Number of trials for Infomap (default: 10)
#'   }
#' @param comparison_params List of parameters for robin comparison:
#'   \itemize{
#'     \item \code{measure}: Stability measure ("vi", "nmi", "split.join",
#'      "adjusted.rand"). Default: "vi"
#'     \item \code{type}: Robin construction type ("dependent", "independent").
#'      Default: "independent"
#'     \item \code{rewire.w.type}: Rewiring strategy for weighted graphs.
#'      Default: "Rewire"
#'   }
#' @param BPPARAM BiocParallel backend for parallel processing.
#' Default: \code{BiocParallel::bpparam()}.
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
#' @importFrom grDevices colorRampPalette
#' @importFrom BiocParallel bpparam
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @examples
#' data(toy_counts)
#'
#'
#' # Infer networks (toy_counts is already a MultiAssayExperiment)
#' networks <- infer_networks(
#'     count_matrices_list = toy_counts,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#' head(networks[[1]])
#'
#' # Generate and symmetrize adjacency matrices (returns SummarizedExperiment)
#' wadj_se <- generate_adjacency(networks)
#' swadj_se <- symmetrize(wadj_se, weight_function = "mean")
#'
#' # Apply cutoff (returns SummarizedExperiment)
#' binary_se <- cutoff_adjacency(
#'     count_matrices = toy_counts,
#'     weighted_adjm_list = swadj_se,
#'     n = 1,
#'     method = "GENIE3",
#'     quantile_threshold = 0.95,
#'     nCores = 1,
#'     debug = TRUE
#' )
#'
#' # Create consensus (returns SummarizedExperiment)
#' consensus <- create_consensus(binary_se, method = "union")
#'
#' # community_path now accepts SummarizedExperiment objects directly
#' comm_cons <- community_path(consensus)
#'
community_path <- function(
    adj_matrix,
    methods = "louvain",
    pathway_db = "KEGG",
    organism = c("human", "mouse"),
    genes_path = 5,
    plot = TRUE,
    verbose = TRUE,
    method_params = list(),
    comparison_params = list(),
    BPPARAM = BiocParallel::bpparam()) {

    organism <- match.arg(organism)
    # Check for optional packages
    if (plot) {
        if (!requireNamespace("ggraph", quietly = TRUE)) {
            stop(
                "Package 'ggraph' is required for plotting. ",
                "Install it with: BiocManager::install('ggraph') ",
                "or set plot = FALSE"
            )
        }
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
            stop(
                "Package 'ggplot2' is required for plotting. ",
                "Install it with: BiocManager::install('ggplot2') ",
                "or set plot = FALSE"
            )
        }
    }

    if (length(methods) > 1 && !requireNamespace("robin", quietly = TRUE)) {
        stop(
            "Package 'robin' is required for comparing multiple methods. ",
            "Install it with: BiocManager::install('robin') ",
            "or use a single method"
        )
    }

    if (pathway_db != "none") {
        if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
            stop(
                "Package 'clusterProfiler' is required . ",
                "Install it with: BiocManager::install('clusterProfiler') ",
                "or set pathway_db = 'none'"
            )
        }
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            stop(
                "Package 'org.Hs.eg.db' is required for pathway enrichment. ",
                "Install it with: BiocManager::install('org.Hs.eg.db') ",
                "or set pathway_db = 'none'"
            )
        }
    }

    # Handle SummarizedExperiment input
    if (inherits(adj_matrix, "SummarizedExperiment")) {
        adj_matrix <- SummarizedExperiment::assay(adj_matrix, 1)
    }

    # Validate input
    if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
        stop(
            "adj_matrix must be a square matrix or a SummarizedExperiment ",
            "containing a square matrix."
        )
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
    comm_res <- .detect_communities(
        graph, methods, method_params, comparison_params, BPPARAM
    )
    best_method <- comm_res$best_method
    best_communities <- comm_res$best_communities
    igraph::V(graph)$community <- as.factor(best_communities)

    non_iso <- igraph::degree(graph) > 0
    if (plot) {
        .plot_communities(graph, best_method)
    }

    if (verbose) message("Running pathway enrichment...")
    pathway_results <- .enrich_communities(
        graph, non_iso, pathway_db, organism, genes_path
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
