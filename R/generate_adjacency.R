#' Generate Adjacency Matrices from Gene Interaction Tables
#'
#' Constructs adjacency matrices from a list of data frames, each
#' representing weighted gene–gene interactions.
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
#' @return A list of square numeric adjacency matrices. Each matrix has
#'   genes as row and column names, and weights as values. Diagonal
#'   entries are set to zero (no self-interactions).
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
#' data("count_matrices")
#'
#' networks <- infer_networks(
#'     count_matrices_list = count_matrices,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#' head(networks[[1]])
#'
#' wadj_list <- generate_adjacency(networks)
#' head(wadj_list[[1]])
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

    adjacency_matrix_list
}
