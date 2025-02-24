pscores <- function(ground_truth, predicted_list, zero_diag = TRUE) {
  
  # Load required packages
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("fmsb", quietly = TRUE)  # fmsb for radar plot
  
  # Ensure ground truth is a matrix
  ground_truth <- as.matrix(ground_truth)
  
  # Optionally zero out the diagonal
  if (zero_diag) diag(ground_truth) <- 0  
  
  # Combine matrices
  all_matrices <- c(list(ground_truth), predicted_list)
  num_matrices <- length(all_matrices)
  
  # Create results dataframe
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
  
  # Helper function: extract upper triangle (excluding diagonal)
  get_upper_tri <- function(mat) {
    mat[upper.tri(mat)]
  }
  
  # Convert ground truth to binary
  ground_truth_upper <- as.numeric(get_upper_tri(ground_truth) > 0)
  
  # Compute statistics for each matrix
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
  
  # Convert numeric columns
  stats_df[, 2:ncol(stats_df)] <- lapply(stats_df[, 2:ncol(stats_df)], as.numeric)
  
  # Remove ground truth row
  stats_df_filtered <- stats_df[-1, ]
  
  # ---- 🚀 **Pentagon Spider Plot with fmsb (Infinite Colors & Scale Labels)** ----
  library(fmsb)
  
  # Select only the five performance metrics
  metrics <- c("TPR", "FPR", "Precision", "F1", "MCC")
  radar_data <- stats_df_filtered %>%
    select(Predicted_Matrix, all_of(metrics))
  
  # Ensure fixed axis limits (0 to 1)
  max_vals <- rep(1, length(metrics))  # Max values set to 1
  min_vals <- rep(0, length(metrics))  # Min values set to 0
  
  # Prepare data for `fmsb` (first row: max, second row: min, then actual values)
  radar_data_scaled <- rbind(max_vals, min_vals, radar_data[, -1])
  
  # Generate an infinite set of colors dynamically based on the number of matrices
  num_matrices <- nrow(stats_df_filtered)
  colors <- rainbow(num_matrices)  # Generates distinct colors for any number of matrices
  
  # Plot the radar chart
  par(mar = c(2, 2, 2, 2))  # Adjust margins
  radarchart(
    radar_data_scaled, 
    axistype = 2, 
    pcol = colors,  # Infinite colors
    pfcol = NA,  # No fill color (only lines)
    plwd = 2, 
    plty = 1, 
    #title = "Pentagon Spider Chart of Performance Metrics",
    cglcol = "grey", cglty = 1, axislabcol = "black",
    caxislabels = seq(0, 1, 0.2),  # More axis labels (0, 0.2, 0.4, ..., 1)
    vlcex = 1.1  # Increase text size
  )
  
  # Add legend with auto-generated colors
  legend("topright", legend = radar_data$Predicted_Matrix, col = colors, lty = 1, lwd = 2)
  
  # Return structured output
  return(list(Statistics = stats_df_filtered))
}
