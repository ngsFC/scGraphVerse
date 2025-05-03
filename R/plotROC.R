#' Plot ROC Curves for Inferred Networks
#'
#' Computes and visualizes Receiver Operating Characteristic (ROC) curves
#' for a list of predicted adjacency matrices compared against a binary ground truth network.
#'
#' @param matrices_list A list of square matrices representing predicted interactions.
#'   Each matrix must have the same dimensions and row/column names as \code{ground_truth}.
#'   Entries may be binary (0/1) or continuous weights.
#' @param ground_truth A square binary matrix indicating true interactions (1) in the upper triangle.
#'   Must match the dimensions and names of each matrix in \code{matrices_list}.
#' @param plot_title Character string. Title for the ROC plot.
#' @param is_binary Logical. If \code{TRUE}, matrices are treated as binary predictions.
#'   Default is \code{FALSE}, for weighted predictions.
#'
#' @return
#' A data frame containing the Area Under the Curve (AUC) for each matrix.
#' The ROC plot is displayed automatically.
#'
#' @details
#' For binary matrices, a single TPR/FPR point is computed per matrix.
#' For weighted matrices, a full ROC curve is calculated based on continuous prediction scores.
#' Diagonal entries are ignored. Symmetry between rows and columns is not enforced.
#'
#' @import ggplot2
#' @importFrom pROC roc
#' @importFrom scales hue_pal
#' @importFrom dplyr bind_rows arrange
#' @export
#'
#' @examples
#' mat1 <- matrix(runif(100), nrow = 10)
#' mat2 <- matrix(runif(100), nrow = 10)
#' ground_truth <- matrix(sample(c(0, 1), 100, replace = TRUE), nrow = 10)
#' diag(ground_truth) <- 0
#' ground_truth[lower.tri(ground_truth)] <- 0
#'
#' plotROC(
#'   matrices_list = list(mat1, mat2),
#'   ground_truth = ground_truth,
#'   plot_title = "ROC for Network Inference",
#'   is_binary = FALSE
#' )
plotROC <- function(matrices_list, ground_truth, plot_title, is_binary = FALSE) {
  truth_vec <- as.vector(ground_truth[upper.tri(ground_truth)])
  
  roc_data <- data.frame()
  auc_values <- data.frame(Matrix = character(), AUC = numeric())
  
  if (is_binary) {
    binary_points <- data.frame(FPR = numeric(), TPR = numeric())
    
    for (i in seq_along(matrices_list)) {
      pred_vec <- .prepare_prediction_vectors(matrices_list[[i]], ground_truth)
      coords <- .compute_binary_roc_point(pred_vec, truth_vec)
      binary_points <- dplyr::bind_rows(binary_points, data.frame(FPR = coords$FPR, TPR = coords$TPR))
    }
    
    binary_points <- dplyr::arrange(binary_points, FPR)
    auc <- sum(diff(c(0, binary_points$FPR, 1)) * (c(0, binary_points$TPR) + c(binary_points$TPR, 1)) / 2)
    
    roc_data <- data.frame(
      FPR = c(0, binary_points$FPR, 1),
      TPR = c(0, binary_points$TPR, 1),
      Matrix = paste("Binary Matrices (AUC=", round(auc, 2), ")", sep = "")
    )
    auc_values <- data.frame(Matrix = "Binary Matrices", AUC = auc)
    p <- .plot_roc_curve(roc_data, binary_points, plot_title)
  } else {
    for (i in seq_along(matrices_list)) {
      pred_vec <- .prepare_prediction_vectors(matrices_list[[i]], ground_truth)
      res <- .compute_weighted_roc_curve(pred_vec, truth_vec, paste("Weighted Matrix", i))
      roc_data <- dplyr::bind_rows(roc_data, res$df)
      auc_values <- dplyr::bind_rows(auc_values, data.frame(Matrix = paste("Weighted Matrix", i), AUC = res$auc))
    }
    p <- .plot_roc_curve(roc_data, NULL, plot_title)
  }
  
  return(list(
    auc = auc_values,
    plot = p
  ))
}

