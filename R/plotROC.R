#' Plot ROC Curves for Inferred Networks
#'
#' Computes and visualizes Receiver Operating Characteristic (ROC) curves for
#' predicted adjacency matrices stored in a \linkS4class{SummarizedExperiment}
#' object, compared against a binary ground truth network.
#'
#' @param predicted_se A \linkS4class{SummarizedExperiment} object containing
#'   predicted adjacency matrices as assays. Each matrix must share dimnames
#'   with \code{ground_truth};entries may be binary (0/1) or continuous weights.
#' @param ground_truth A square binary matrix indicating true interactions (1)
#'   in the upper triangle. Must match dims and names of each assay in
#'   \code{predicted_se}.
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
#' @importFrom dplyr bind_rows arrange
#' @importFrom SummarizedExperiment assays
#' @export
#'
#' @examples
#' data(toy_counts)
#' data(toy_adj_matrix)
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
#' # plotROC now accepts SummarizedExperiment directly
#' roc_res <- plotROC(
#'     swadj_se,
#'     toy_adj_matrix,
#'     plot_title = "ROC Curve: GENIE3",
#'     is_binary = FALSE
#' )
#' roc_res$plot
#' auc_joint <- roc_res$auc
plotROC <- function(
    predicted_se,
    ground_truth,
    plot_title,
    is_binary = FALSE) {
    # Check for plotting packages
    BiocBaseUtils::checkInstalled("ggplot2")
    if (!is_binary) {
        BiocBaseUtils::checkInstalled("pROC")
    }

    # Accept SummarizedExperiment and extract assays
    if (inherits(predicted_se, "SummarizedExperiment")) {
        matrices_list <- as.list(SummarizedExperiment::assays(predicted_se))
    } else {
        stop("predicted_se must be a SummarizedExperiment object")
    }

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
