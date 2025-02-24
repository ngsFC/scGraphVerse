pscores <- function(ground_truth, predicted_list) {
  
  # Ensure ground truth is a matrix and zero out the diagonal
  ground_truth <- as.matrix(ground_truth)
  diag(ground_truth) <- 0  # Set diagonal to zero
  
  # Combine ground truth with predicted list
  all_matrices <- c(list(ground_truth), predicted_list)
  num_matrices <- length(all_matrices)
  
  # Create a dataframe to store statistics (remove Accuracy, add MCC)
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
    MCC = numeric(num_matrices),
    stringsAsFactors = FALSE
  )
  
  # Helper function: get upper triangular (excluding diagonal)
  get_upper_tri <- function(mat) {
    mat[upper.tri(mat)]
  }
  
  # Binarize the ground truth upper triangle (edge = 1)
  ground_truth_upper <- ifelse(get_upper_tri(ground_truth) > 0, 1, 0)
  
  # Loop over each matrix (including ground truth)
  for (i in 1:num_matrices) {
    matrix_i <- as.matrix(all_matrices[[i]])
    upper_i <- get_upper_tri(matrix_i)
    binary_i <- ifelse(upper_i > 0, 1, 0)
    
    # Count edges and nodes
    stats_df$Edges[i] <- sum(binary_i)
    stats_df$Nodes[i] <- nrow(matrix_i)
    
    # Compute confusion matrix components
    TP <- sum(binary_i == 1 & ground_truth_upper == 1)
    TN <- sum(binary_i == 0 & ground_truth_upper == 0)
    FP <- sum(binary_i == 1 & ground_truth_upper == 0)
    FN <- sum(binary_i == 0 & ground_truth_upper == 1)
    
    stats_df$TP[i] <- TP
    stats_df$TN[i] <- TN
    stats_df$FP[i] <- FP
    stats_df$FN[i] <- FN
    
    # Calculate performance metrics
    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)           # Recall
    FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    F1 <- ifelse((Precision + TPR) > 0, 2 * (Precision * TPR) / (Precision + TPR), 0)
    
    # Compute MCC (Matthews Correlation Coefficient)
    denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    MCC <- ifelse(denominator > 0, (TP * TN - FP * FN) / denominator, 0)
    
    stats_df$TPR[i] <- TPR
    stats_df$FPR[i] <- FPR
    stats_df$Precision[i] <- Precision
    stats_df$F1[i] <- F1
    stats_df$MCC[i] <- MCC
  }
  
  # Remove the ground truth row for plotting purposes
  stats_df_filtered <- stats_df[-1, ]
  
  # Prepare data for the radar plot: select the desired metrics
  ordered_metrics <- c("TPR", "FPR", "F1", "Precision", "MCC")
  radar_data <- stats_df_filtered[, ordered_metrics]
  rownames(radar_data) <- stats_df_filtered$Matrix
  
  # fmsb requires the first two rows to be max and min values.
  # Here we set max for TPR, FPR, F1, Precision to 1 and MCC to 1;
  # for the minimum, we use 0 (or -1 for MCC)
  max_values <- c(TPR = 1, FPR = 1, F1 = 1, Precision = 1, MCC = 1)
  min_values <- c(TPR = 0, FPR = 0, F1 = 0, Precision = 0, MCC = -1)
  radar_data <- rbind(max_values, min_values, radar_data)
  
  # Plot the radar chart using the fmsb package
  library(fmsb)
  # Generate colors for the predicted matrices
  num_pred <- nrow(radar_data) - 2
  colors_border <- rainbow(num_pred)
  
  radarchart(radar_data, axistype = 1,
             # Customize polygon appearance
             pcol = colors_border, pfcol = scales::alpha(colors_border, 0.3), plwd = 2,
             # Customize grid appearance
             cglcol = "grey", cglty = 1, axislabcol = "grey",
             # Adjust axis labels (accounting for MCC range)
             caxislabels = seq(-1, 1, 0.5),
             cglwd = 0.8,
             title = "Radar Plot of Performance Metrics")
  
  legend("topright", legend = rownames(radar_data[-c(1,2),]), bty = "n",
         pch = 20, col = colors_border, text.col = "black", cex = 0.8)
  
  # Return the statistics data frame for the predicted matrices
  list(Statistics = stats_df_filtered)
}
