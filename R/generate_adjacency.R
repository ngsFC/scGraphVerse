#' Generate Adjacency Matrices from Gene Interaction Tables
#'
#' Constructs adjacency matrices from a list of data frames (network edge lists)
#' and returns them in a \linkS4class{SummarizedExperiment} object.
#'
#' @param df_list A list of data frames. Each data frame must have three
#'   columns:
#'   \describe{
#'     \item{Gene1}{Character. First gene in the interaction.}
#'     \item{Gene2}{Character. Second gene in the interaction.}
#'     \item{Weight}{Numeric. Weight or strength of the interaction from
#'       \code{Gene1} to \code{Gene2}.}
#'   }
#' @param nCores Integer. Number of CPU cores to use for parallel
#'   processing. Defaults to the number of available workers from the
#'   current \pkg{BiocParallel} backend.
#'
#' @return A \linkS4class{SummarizedExperiment} object where each assay is a
#'   square numeric adjacency matrix (p×p genes). Diagonal entries are set to
#'   zero (no self-interactions).
#'
#' @details The function first identifies all unique genes across all data
#'   frames to define the matrix dimensions. For each interaction table,
#'   it places the corresponding weights at the appropriate gene-pair
#'   positions. Parallelization is handled by \pkg{BiocParallel} for
#'   improved performance on multiple datasets.
#'
#'   Missing weights (\code{NA}) are ignored during construction. Only
#'   gene pairs matching the global gene list are inserted.
#'
#' @importFrom BiocParallel bplapply bpworkers bpparam MulticoreParam
#' @export
#'
#' @examples
#' data("toy_counts")
#'
#' # Infer networks (toy_counts is already a MultiAssayExperiment)
#' networks <- infer_networks(
#'     count_matrices_list = toy_counts,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#' head(networks[[1]])
#'
#' # Generate adjacency matrices
#' wadj_se <- generate_adjacency(networks) # returns SummarizedExperiment
#' head(wadj_se[[1]])
generate_adjacency <- function(
    df_list,
    nCores = 1) {
    if (!is.list(df_list) ||
        !all(vapply(df_list, is.data.frame, logical(1)))) {
        stop("df_list must be a list of data frames")
    }

    BPPARAM <- BiocParallel::MulticoreParam(nCores)

    all_genes <- sort(unique(unlist(
        BiocParallel::bplapply(
            df_list,
            function(data) {
                unique(
                    c(
                        as.character(data[[1]]),
                        as.character(data[[2]])
                    )
                )
            },
            BPPARAM = BPPARAM
        )
    )))

    template_matrix <- matrix(
        0,
        nrow     = length(all_genes),
        ncol     = length(all_genes),
        dimnames = list(all_genes, all_genes)
    )

    adjacency_matrix_list <- BiocParallel::bplapply(
        df_list,
        function(data) {
            if (ncol(data) < 3) {
                stop(
                    "Each data frame must have at least 3 columns ",
                    "(source, target, weight)"
                )
            }

            adjacency_matrix <- template_matrix

            for (i in seq_len(nrow(data))) {
                gene1 <- as.character(data[[1]][i])
                gene2 <- as.character(data[[2]][i])
                weight <- as.numeric(data[[3]][i])

                if (!is.na(gene1) &&
                    !is.na(gene2) &&
                    !is.na(weight) &&
                    gene1 %in% rownames(adjacency_matrix) &&
                    gene2 %in% colnames(adjacency_matrix)
                ) {
                    adjacency_matrix[gene1, gene2] <- weight
                }
            }

            diag(adjacency_matrix) <- 0
            adjacency_matrix
        },
        BPPARAM = BPPARAM
    )

    if (is.null(names(adjacency_matrix_list))) {
        names(adjacency_matrix_list) <- paste0(
            "network_",
            seq_along(adjacency_matrix_list)
        )
    }

    build_network_se(
        networks = adjacency_matrix_list,
        networkData = S4Vectors::DataFrame(
            network = names(adjacency_matrix_list),
            n_edges = vapply(
                adjacency_matrix_list, function(x) sum(x > 0),
                numeric(1)
            ),
            row.names = names(adjacency_matrix_list)
        ),
        metadata = list(type = "weighted")
    )
}
