#' Calculate performance metrics for predicted adjacency matrices against ground truth
#'
#' This function computes a variety of performance metrics (such as precision, recall, F1 score, 
#' accuracy, true positive rate, false positive rate) for a list of predicted adjacency matrices 
#' by comparing them with a provided ground truth matrix. The function calculates metrics 
#' for the upper triangular part of the adjacency matrices, excluding the diagonal.
#' It then creates a bar plot comparing the performance of each predicted matrix.
#'
#' @param ground_truth A square matrix representing the ground truth adjacency matrix. The diagonal 
#'        elements will be set to zero.
#' @param predicted_list A list of square matrices representing predicted adjacency matrices that 
#'        will be compared to the ground truth.
#'
#' @return A list containing:
#'   - `Statistics`: A data frame with the performance metrics (TP, TN, FP, FN, TPR, FPR, Precision, 
#'      F1, Accuracy) for each matrix in the `predicted_list`.
#'   - A bar plot displaying the comparison of metrics (TPR, FPR, F1, Precision, Accuracy) 
#'      across all matrices.
#' 
#' @details 
#' This function computes the following performance metrics:
#'   - **True Positive Rate (TPR)**: The proportion of actual positive edges that were correctly 
#'     predicted (also known as recall).
#'   - **False Positive Rate (FPR)**: The proportion of actual negative edges that were incorrectly 
#'     predicted as positive.
#'   - **Precision**: The proportion of predicted positive edges that were actually positive.
#'   - **F1 Score**: The harmonic mean of Precision and TPR.
#'   - **Accuracy**: The proportion of correct predictions (TP + TN) to all predictions (TP + TN + FP + FN).
#' 
#' The function also excludes diagonal elements from both the ground truth and predicted matrices,
#' and operates only on the upper triangular part of the adjacency matrices (which represents the 
#' directed edges).
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text position_dodge labs theme_minimal
#' @importFrom tidyr pivot_longer
#' 
#' @export
pscores <- function(ground_truth, predicted_list) {
  
  # Ensure ground truth is a matrix and zero out the diagonal
  ground_truth <- as.matrix(ground_truth)
  diag(ground_truth) <- 0  # Set diagonal to zero
  
  # Combine ground truth with predicted list
  all_matrices <- c(list(ground_truth), predicted_list)
  num_matrices <- length(all_matrices)
  
  # Create a dataframe to store additional statistics
  stats_df <- data.frame(
    Matrix = c("Ground Truth", paste("Matrix", seq_along(predicted_list))),
    Edges = integer(num_matrices),
    Nodes = integer(num_matrices),
    TP = integer(num_matrices),
    TN = integer(num_matrices),
    FP = integer(num_matrices),
    FN = integer(num_matrices),
    TPR = numeric(num_matrices),
    FPR = numeric(num_matrices),
    Precision = numeric(num_matrices),
    F1 = numeric(num_matrices),
    Accuracy = numeric(num_matrices),
    stringsAsFactors = FALSE
  )
  
  # Function to get the upper triangular part of a matrix (excluding diagonal)
  get_upper_tri <- function(mat) {
    mat[upper.tri(mat)]  # Extracts the upper triangular elements
  }
  
  # Binary version of the ground truth upper triangle
  ground_truth_upper <- ifelse(get_upper_tri(ground_truth) > 0, 1, 0)
  
  # Calculate statistics for each predicted matrix
  for (i in 1:num_matrices) {
    matrix_i <- as.matrix(all_matrices[[i]])
    upper_i <- get_upper_tri(matrix_i)
    binary_i <- ifelse(upper_i > 0, 1, 0)
    
    # Calculate the number of edges (non-zero values) and nodes (size of the matrix)
    stats_df$Edges[i] <- sum(binary_i)  # Count non-zero values in the upper triangle
    stats_df$Nodes[i] <- nrow(matrix_i)  # Number of nodes is the number of rows (or cols)
    
    # Calculate TP, TN, FP, FN
    TP <- sum(binary_i == 1 & ground_truth_upper == 1)
    TN <- sum(binary_i == 0 & ground_truth_upper == 0)
    FP <- sum(binary_i == 1 & ground_truth_upper == 0)
    FN <- sum(binary_i == 0 & ground_truth_upper == 1)
    
    # Store the counts
    stats_df$TP[i] <- TP
    stats_df$TN[i] <- TN
    stats_df$FP[i] <- FP
    stats_df$FN[i] <- FN
    
    # Calculate metrics
    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)  # True Positive Rate (Recall)
    FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)  # False Positive Rate
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)  # Precision
    F1 <- ifelse((Precision + TPR) > 0, 2 * (Precision * TPR) / (Precision + TPR), 0)  # F1 Score
    Accuracy <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), 0)  # Accuracy
    
    # Store the metrics
    stats_df$TPR[i] <- TPR
    stats_df$FPR[i] <- FPR
    stats_df$Precision[i] <- Precision
    stats_df$F1[i] <- F1
    stats_df$Accuracy[i] <- Accuracy
  }
  
  # Remove ground truth from stats_df for plotting
  stats_df_filtered <- stats_df[-1, ]  # Exclude the first row (Ground Truth)
  
  # Reorder the columns in the desired metric order
  ordered_metrics <- c("TPR", "FPR", "F1", "Precision", "Accuracy")
  
  long_stats_df <- stats_df_filtered %>%
    tidyr::pivot_longer(cols = all_of(ordered_metrics), names_to = "Metric", values_to = "Value")
  
  # Reorder the 'Metric' factor levels to ensure the desired order
  long_stats_df$Metric <- factor(long_stats_df$Metric, levels = ordered_metrics)
  
  # Plot the metrics comparison across matrices
  plot <- ggplot(long_stats_df, aes(x = Matrix, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5) +
    labs(title = "Metrics Comparison Across Matrices", x = "Matrix Index", y = "Value") +
    theme_minimal()
  
  print(plot)
  
  # Return the statistics data frame (filtered for plotting)
  list(Statistics = stats_df_filtered)
}

