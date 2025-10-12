#' Create a SummarizedExperiment for Network Storage
#'
#' Constructs a SummarizedExperiment container for multiple gene regulatory
#' network adjacency matrices with shared gene space (p×p matrices).
#'
#' @param networks A list of adjacency matrices (all must have same dimensions)
#' @param networkData A DataFrame with metadata for each network
#' @param geneData Optional. A DataFrame with gene-level annotations
#' @param metadata Optional. List of global metadata
#'
#' @return A SummarizedExperiment object where each assay is a p×p network
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' # Example 1: Building SE from a list of adjacency matrices
#' data("toy_adj_matrix")
#'
#' # Create a list of network matrices
#' net_list <- list(
#'     network1 = toy_adj_matrix,
#'     network2 = toy_adj_matrix
#' )
#'
#' # Build SummarizedExperiment
#' network_se <- build_network_se(net_list)
#' network_se
#'
#' # Example 2: Using with inferred networks
#' data("toy_counts")
#'
#' networks <- infer_networks(
#'     count_matrices_list = toy_counts,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#'
#' # generate_adjacency() internally uses build_network_se()
#' wadj_se <- generate_adjacency(networks)
#' wadj_se
build_network_se <- function(networks,
    networkData = NULL,
    geneData = NULL,
    metadata = list()) {
    gene_names <- .validate_network_list(networks)

    if (is.null(networkData)) {
        net_names <- names(networks)
        if (is.null(net_names)) {
            net_names <- paste0("network_", seq_along(networks))
            names(networks) <- net_names
        }
        networkData <- S4Vectors::DataFrame(
            network = net_names,
            n_edges = vapply(networks, sum, numeric(1)),
            row.names = net_names
        )
    }

    if (is.null(geneData)) {
        geneData <- S4Vectors::DataFrame(
            gene = gene_names,
            row.names = gene_names
        )
    }

    assay_list <- networks
    names(assay_list) <- rownames(networkData)

    SummarizedExperiment::SummarizedExperiment(
        assays = assay_list,
        rowData = geneData,
        colData = geneData,
        metadata = c(metadata, list(object_type = "network_collection"))
    )
}
