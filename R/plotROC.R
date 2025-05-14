#' Plot ROC Curves for Inferred Networks
#'
#' Computes and visualizes Receiver Operating Characteristic (ROC) curves for
#' a list of predicted adjacency matrices compared against a binary ground
#' truth network.
#'
#' @param matrices_list A list of square matrices representing predicted
#'   interactions. Each must share dimnames with \code{ground_truth}; entries
#'   may be binary (0/1) or continuous weights.
#' @param ground_truth A square binary matrix indicating true interactions (1)
#'   in the upper triangle. Must match dims and names of each element of
#'   \code{matrices_list}.
#' @param plot_title Character string. Title for the ROC plot.
#' @param is_binary Logical. If \code{TRUE}, treat matrices as binary
#'   predictions. Default \code{FALSE} for weighted predictions.
#'
#' @return A list with:\cr
#'   \code{auc}: data frame of AUC per matrix.\cr
#'   \code{plot}: the ROC plot (via ggplot2).
#'
#' @details For binary matrices, a single TPR/FPR point is computed per
#'   matrix. For weighted ones, a full ROC curve is computed from continuous
#'   scores. Diagonals are ignored; symmetry is not enforced.
#'
#' @import ggplot2
#' @importFrom pROC roc
#' @importFrom scales hue_pal
#' @importFrom dplyr bind_rows arrange
#' @export
#'
#' @examples
#' data(count_matrices)
#' data(adj_truth)
#'
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
#' roc_res <- plotROC(
#'     swadj_list,
#'     adj_truth,
#'     plot_title = "ROC Curve: GENIE3",
#'     is_binary = FALSE
#' )
#' roc_res$plot
#' auc_joint <- roc_res$auc
plotROC <- function(
    matrices_list,
    ground_truth,
    plot_title,
    is_binary = FALSE) {
    truth_vec <- as.vector(
        ground_truth[upper.tri(ground_truth)]
    )

    roc_data <- data.frame()
    auc_values <- data.frame(
        Matrix = character(),
        AUC    = numeric()
    )

    if (is_binary) {
        binary_points <- data.frame(
            FPR = numeric(),
            TPR = numeric()
        )

        for (i in seq_along(matrices_list)) {
            pred_vec <- .prepare_prediction_vectors(
                matrices_list[[i]],
                ground_truth
            )
            coords <- .compute_binary_roc_point(
                pred_vec,
                truth_vec
            )
            binary_points <- dplyr::bind_rows(
                binary_points,
                data.frame(
                    FPR = coords$FPR,
                    TPR = coords$TPR
                )
            )
        }

        binary_points <- dplyr::arrange(
            binary_points,
            FPR
        )
        auc <- sum(
            diff(c(0, binary_points$FPR, 1)) *
                (c(0, binary_points$TPR) +
                    c(binary_points$TPR, 1)) / 2
        )

        roc_data <- data.frame(
            FPR = c(0, binary_points$FPR, 1),
            TPR = c(0, binary_points$TPR, 1),
            Matrix = paste(
                "Binary Matrices (AUC=", round(auc, 2), ")",
                sep = ""
            )
        )
        auc_values <- data.frame(
            Matrix = "Binary Matrices",
            AUC    = auc
        )
        p <- .plot_roc_curve(
            roc_data,
            binary_points,
            plot_title
        )
    } else {
        for (i in seq_along(matrices_list)) {
            pred_vec <- .prepare_prediction_vectors(
                matrices_list[[i]],
                ground_truth
            )
            res <- .compute_weighted_roc_curve(
                pred_vec,
                truth_vec,
                paste("Weighted Matrix", i)
            )
            roc_data <- dplyr::bind_rows(roc_data, res$df)
            auc_values <- dplyr::bind_rows(
                auc_values,
                data.frame(
                    Matrix = paste("Weighted Matrix", i),
                    AUC    = res$auc
                )
            )
        }
        p <- .plot_roc_curve(roc_data, NULL, plot_title)
    }

    list(
        auc  = auc_values,
        plot = p
    )
}
