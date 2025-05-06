#' Compute Performance Scores for Predicted Adjacency Matrices
#'
#' Computes classification metrics by comparing predicted adjacency matrices
#' to a ground truth binary network. Visualizes the performance scores using a radar (spider) plot.
#'
#' @param ground_truth A square binary adjacency matrix representing the ground truth network.
#'   Values must be 0 or 1. Only the upper triangle is used for evaluation.
#' @param predicted_list A list of predicted adjacency matrices to evaluate.
#'   Each must have the same dimensions and row/column names as \code{ground_truth}.
#' @param zero_diag Logical. If \code{TRUE} (default), sets the diagonal of \code{ground_truth} to zero
#'   before evaluation, removing self-loops.
#'
#' @return
#' A list with one element:
#' \itemize{
#'   \item \code{Statistics}: A data frame containing evaluation metrics (TP, TN, FP, FN, TPR, FPR, Precision, F1, MCC)
#'   for each predicted matrix.
#' }
#'
#' @details
#' For each predicted matrix, the confusion matrix is computed using the upper triangle
#' (non-self edges). Metrics including True Positive Rate (TPR), False Positive Rate (FPR),
#' Precision, F1-score, and Matthews Correlation Coefficient (MCC) are calculated.
#'
#' A radar plot is automatically generated summarizing the key scores across matrices.
#'
#' @note
#' Requires the \pkg{fmsb}, \pkg{dplyr}, and \pkg{tidyr} packages.
#'
#' @importFrom fmsb radarchart
#' @importFrom dplyr select
#' @importFrom tidyr all_of
#' @export
#'
#' @examples
#' # Simulate ground truth and predictions
#' ground_truth <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#' diag(ground_truth) <- 0
#' pred1 <- ground_truth
#' pred2 <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
#'
#' # Compute scores and generate radar plot
#' result <- pscores(ground_truth, list(pred1, pred2))
#' result$Statistics
pscores <- function(ground_truth, predicted_list, zero_diag = TRUE) {
  ground_truth <- as.matrix(ground_truth)
  if (zero_diag) diag(ground_truth) <- 0
  gt_u <- ground_truth[upper.tri(ground_truth)]
  
  stats_list <- lapply(seq_along(predicted_list), function(i) {
    pred <- as.matrix(predicted_list[[i]])
    pu   <- as.numeric(pred[upper.tri(pred)] > 0)
    .compute_confusion_metrics(pu, gt_u, i)
  })
  stats_df <- do.call(rbind, stats_list)
  
  all_metrics <- c("TPR", "FPR", "Precision", "F1", "MCC")
  radar_all   <- .plot_metrics_radar(stats_df, all_metrics)
  
  list(
    Statistics = stats_df,
    Radar      = radar_all
  )
}
