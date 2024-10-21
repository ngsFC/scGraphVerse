pscores <- function(ground_truth, predicted_list) {
  
  # Ensure ground truth is a matrix and zero out the diagonal
  ground_truth <- as.matrix(ground_truth)
  diag(ground_truth) <- 0  # Set diagonal to zero
  
  # Combine ground truth with predicted list
  all_matrices <- c(list(ground_truth), predicted_list)
  num_matrices <- length(all_matrices)
  
  # Prepare matrices to store results
  jaccard_matrix <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- c("Ground Truth", paste("Matrix", seq_along(predicted_list)))
  
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
    
    # Calculate Jaccard Index with Ground Truth
    if (i > 1) {  # Skip the ground truth comparison with itself
      intersection <- sum(binary_i & ground_truth_upper)
      union <- sum(binary_i | ground_truth_upper)
      jaccard_matrix[i, 1] <- ifelse(union > 0, intersection / union, 0)
    }
  }
  
  # Plot the metrics using ggplot2
  library(ggplot2)
  long_stats_df <- stats_df %>%
    tidyr::pivot_longer(cols = c(TPR, FPR, Precision, F1, Accuracy), names_to = "Metric", values_to = "Value")
  
  plot <- ggplot(long_stats_df, aes(x = Matrix, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5) +
    labs(title = "Metrics Comparison Across Matrices", x = "Matrix Index", y = "Value") +
    theme_minimal()
  
  # Print the plot
  print(plot)
  
  # Return both the Jaccard matrix and stats dataframe
  list(Jaccard_Matrix = jaccard_matrix, Statistics = stats_df)
}
