#' Toy MultiAssayExperiment for Network Inference
#'
#' A simulated dataset generated using zinb_simdata() for demonstrating
#' network inference functions in the scGraphVerse package.
#'
#' @format A \linkS4class{MultiAssayExperiment} object containing 3
#' SingleCellExperiment objects (experiments), each with 10 genes x 40 cells.
#' The data was simulated using a known ground truth network (toy_adj_matrix)
#' with zero-inflated negative binomial distributions.
#'
#' @usage data(toy_counts)
#' @examples
#' data(toy_counts)
#' toy_counts
#'
#' # Access individual experiments
#' MultiAssayExperiment::experiments(toy_counts)
#'
#' # Use directly with infer_networks
#' networks <- infer_networks(
#'     count_matrices_list = toy_counts,
#'     method = "GENIE3",
#'     nCores = 1
#' )
"toy_counts"
