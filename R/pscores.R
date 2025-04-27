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
  if (!is.matrix(ground_truth) || nrow(ground_truth) != ncol(ground_truth)) {
    stop("`ground_truth` must be a square matrix.")
  }

  if (!all(ground_truth %in% c(0, 1))) {
    stop("`ground_truth` must contain only binary values (0/1).")
  }

  if (!is.list(predicted_list) || !all(vapply(predicted_list, is.matrix))) {
    stop("`predicted_list` must be a list of matrices.")
  }

  ground_truth <- as.matrix(ground_truth)
  if (zero_diag) diag(ground_truth) <- 0

  # Extract upper triangle values for ground truth
  get_upper_tri <- function(mat) mat[upper.tri(mat)]
  ground_truth_upper <- get_upper_tri(ground_truth)

  metrics <- c("TP", "TN", "FP", "FN", "TPR", "FPR", "Precision", "F1", "MCC")
  stat_rows <- lapply(seq_along(predicted_list), function(i) {
    pred <- predicted_list[[i]]
    if (!all(dim(pred) == dim(ground_truth))) {
      stop(sprintf("Predicted matrix %d has mismatched dimensions.", i))
    }
    pred_upper <- as.numeric(get_upper_tri(pred) > 0)
    gt_upper <- ground_truth_upper

    TP <- sum(pred_upper == 1 & gt_upper == 1)
    TN <- sum(pred_upper == 0 & gt_upper == 0)
    FP <- sum(pred_upper == 1 & gt_upper == 0)
    FN <- sum(pred_upper == 0 & gt_upper == 1)

    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    F1 <- ifelse((Precision + TPR) > 0, 2 * (Precision * TPR) / (Precision + TPR), 0)
    denominator <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
      as.numeric(TN + FP) * as.numeric(TN + FN))
    MCC <- ifelse(denominator > 0, (TP * TN - FP * FN) / denominator, 0)

    data.frame(
      Predicted_Matrix = paste("Matrix", i),
      TP, TN, FP, FN, TPR, FPR, Precision, F1, MCC
    )
  })

  stats_df <- do.call(rbind, stat_rows)

  # Radar chart
  radar_metrics <- c("TPR", "FPR", "Precision", "F1", "MCC")
  radar_data <- stats_df %>%
    dplyr::select(Predicted_Matrix, dplyr::all_of(radar_metrics))

  radar_scaled <- rbind(rep(1, length(radar_metrics)), rep(0, length(radar_metrics)), radar_data[, -1])
  colors <- rainbow(nrow(stats_df))

  par(mar = c(2, 2, 2, 2))
  fmsb::radarchart(
    radar_scaled,
    axistype = 2,
    pcol = colors,
    plty = 1,
    plwd = 2,
    cglcol = "grey",
    caxislabels = seq(0, 1, 0.2),
    vlcex = 1.1
  )
  legend("topright", legend = radar_data$Predicted_Matrix, col = colors, lty = 1, lwd = 2)

  return(list(Statistics = stats_df))
}
