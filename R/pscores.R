#' Performance Scores for Predicted Matrices
#' 
#' This function evaluates the performance of predicted adjacency matrices against a ground truth matrix.
#' It computes various classification metrics such as True Positive Rate (TPR), False Positive Rate (FPR),
#' Precision, F1-score, and Matthews Correlation Coefficient (MCC). Additionally, it visualizes the results
#' using a radar chart.
#'
#' @param ground_truth A square adjacency matrix representing the ground truth network.
#' @param predicted_list A list of square adjacency matrices to be evaluated against the ground truth.
#' @param zero_diag Logical; if TRUE, sets the diagonal of the ground truth matrix to zero.
#' 
#' @return A list containing a data frame with computed statistics for each predicted matrix.
#' 
#' @importFrom dplyr select
#' @importFrom tidyr all_of
#' @importFrom fmsb radarchart
#' 
#' @export
pscores <- function(ground_truth, predicted_list, zero_diag = TRUE) {
  
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("fmsb", quietly = TRUE)  # fmsb for radar plot
  
  ground_truth <- as.matrix(ground_truth)
  
  if (zero_diag) diag(ground_truth) <- 0  
  
  all_matrices <- c(list(ground_truth), predicted_list)
  num_matrices <- length(all_matrices)
  
  stats_df <- data.frame(
    Predicted_Matrix = c("Ground Truth", paste("Matrix", seq_along(predicted_list))),
    TP = integer(num_matrices),
    TN = integer(num_matrices),
    FP = integer(num_matrices),
    FN = integer(num_matrices),
    TPR = numeric(num_matrices),
    FPR = numeric(num_matrices),
    Precision = numeric(num_matrices),
    F1 = numeric(num_matrices),
    MCC = numeric(num_matrices),
    stringsAsFactors = FALSE
  )
  
  get_upper_tri <- function(mat) {
    mat[upper.tri(mat)]
  }
  
  ground_truth_upper <- as.numeric(get_upper_tri(ground_truth) > 0)
  
  stats_df[-1, ] <- do.call(rbind, lapply(seq_along(predicted_list), function(i) {
    matrix_i <- as.matrix(predicted_list[[i]])
    binary_i <- as.numeric(get_upper_tri(matrix_i) > 0)
    
    TP <- sum(binary_i == 1 & ground_truth_upper == 1)
    TN <- sum(binary_i == 0 & ground_truth_upper == 0)
    FP <- sum(binary_i == 1 & ground_truth_upper == 0)
    FN <- sum(binary_i == 0 & ground_truth_upper == 1)
    
    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    F1 <- ifelse((Precision + TPR) > 0, 2 * (Precision * TPR) / (Precision + TPR), 0)
    
    denominator <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
                          as.numeric(TN + FP) * as.numeric(TN + FN))
    MCC <- ifelse(denominator > 0, (TP * TN - FP * FN) / denominator, 0)
    
    c(Predicted_Matrix = paste("Matrix", i), TP, TN, FP, FN, TPR, FPR, Precision, F1, MCC)
  }))
  
  stats_df[, 2:ncol(stats_df)] <- lapply(stats_df[, 2:ncol(stats_df)], as.numeric)
  stats_df_filtered <- stats_df[-1, ]
  
  metrics <- c("TPR", "FPR", "Precision", "F1", "MCC")
  radar_data <- stats_df_filtered %>%
    dplyr::select(Predicted_Matrix, all_of(metrics))
  
  max_vals <- rep(1, length(metrics))  # Max values set to 1
  min_vals <- rep(0, length(metrics))  # Min values set to 0
  
  radar_data_scaled <- rbind(max_vals, min_vals, radar_data[, -1])
  
  num_matrices <- nrow(stats_df_filtered)
  colors <- rainbow(num_matrices)  # Generates distinct colors for any number of matrices
  
  par(mar = c(2, 2, 2, 2))
  radarchart(
    radar_data_scaled, 
    axistype = 2, 
    pcol = colors,
    pfcol = NA, 
    plwd = 2, 
    plty = 1, 
    #title = "Pentagon Spider Chart of Performance Metrics",
    cglcol = "grey", cglty = 1, axislabcol = "black",
    caxislabels = seq(0, 1, 0.2),
    vlcex = 1.1 
  )
  
  legend("topright", legend = radar_data$Predicted_Matrix, col = colors, lty = 1, lwd = 2)
  
  return(list(Statistics = stats_df_filtered))
}
