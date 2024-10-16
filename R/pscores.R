pscores <- function(ground_truth, predicted_list) {
  
  # Ensure that ground_truth is a matrix
  ground_truth <- as.matrix(ground_truth)
  
  # Step 1: Calculate the Jaccard Index between each pair of matrices (including ground truth and predicted)
  all_matrices <- c(list(ground_truth), predicted_list)
  num_matrices <- length(all_matrices)
  
  # Create an empty matrix to store the Jaccard Index values
  jaccard_matrix <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- c("Ground Truth", paste("Matrix", seq_along(predicted_list)))
  
  # Calculate Jaccard Index for each pair of matrices
  for (i in 1:num_matrices) {
    for (j in 1:num_matrices) {
      if (i == j) {
        jaccard_matrix[i, j] <- 1
      } else {
        matrix_i <- as.matrix(all_matrices[[i]])
        matrix_j <- as.matrix(all_matrices[[j]])
        binary_i <- ifelse(matrix_i > 0, 1, 0)
        binary_j <- ifelse(matrix_j > 0, 1, 0)
        
        intersection <- sum(binary_i & binary_j)
        union <- sum(binary_i | binary_j)
        jaccard_index <- ifelse(union > 0, intersection / union, 0)
        jaccard_matrix[i, j] <- jaccard_index
      }
    }
  }
  
  # Step 2: Create a heatmap for the Jaccard Index among all matrices
  jaccard_df <- as.data.frame(as.table(jaccard_matrix))
  colnames(jaccard_df) <- c("Matrix1", "Matrix2", "Jaccard_Index")
  
  jaccard_heatmap <- ggplot(jaccard_df, aes(x = Matrix1, y = Matrix2, fill = Jaccard_Index)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0.5, space = "Lab", name = "Jaccard Index") +
    geom_text(aes(label = round(Jaccard_Index, 2)), color = "black", size = 4) +
    theme_minimal() +
    ggtitle("Jaccard Index Heatmap Among All Matrices") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  # Step 3: Calculate metrics for each predicted matrix against the ground truth
  metrics_results <- data.frame(
    Matrix_Index = character(),
    TP = numeric(),
    TN = numeric(),
    FP = numeric(),
    FN = numeric(),
    TPR = numeric(),      # True Positive Rate (Recall)
    FPR = numeric(),      # False Positive Rate
    Precision = numeric(), 
    Recall = numeric(),   # Same as TPR
    F1_Score = numeric(),
    Accuracy = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(predicted_list)) {
    predicted <- as.matrix(predicted_list[[i]])
    
    if (!all(dim(predicted) == dim(ground_truth))) {
      stop(paste("Matrix dimensions do not match for matrix at index", i))
    }
    
    predicted_binary <- ifelse(predicted > 0, 1, 0)
    
    TP <- sum((ground_truth == 1) & (predicted_binary == 1))
    TN <- sum((ground_truth == 0) & (predicted_binary == 0))
    FP <- sum((ground_truth == 0) & (predicted_binary == 1))
    FN <- sum((ground_truth == 1) & (predicted_binary == 0))
    
    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)  # True Positive Rate / Recall
    FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)  # False Positive Rate
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)  # Precision
    Recall <- TPR  # Recall is the same as TPR
    F1_Score <- ifelse((Precision + Recall) > 0, 2 * (Precision * Recall) / (Precision + Recall), 0)  # F1 Score
    Accuracy <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), 0)  # Accuracy
    
    # Store the metrics
    metrics_results <- rbind(
      metrics_results,
      data.frame(
        Matrix_Index = paste("Matrix", i),
        TP = TP,
        TN = TN,
        FP = FP,
        FN = FN,
        TPR = TPR,
        FPR = FPR,
        Precision = Precision,
        Recall = Recall,
        F1_Score = F1_Score,
        Accuracy = Accuracy,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Step 4: Create a bar plot for all metrics
  # Melt the metrics data for plotting
  melted_metrics <- melt(metrics_results, id.vars = "Matrix_Index", 
                         measure.vars = c("TPR", "FPR", "Precision", "Recall", "F1_Score", "Accuracy"))
  
  # Create the bar plot
  metrics_barplot <- ggplot(melted_metrics, aes(x = factor(Matrix_Index), y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(value, 2)), vjust = -0.3, position = position_dodge(width = 0.9), size = 3.5) +
    scale_fill_brewer(palette = "Set3", name = "Metrics") +
    ggtitle("Metrics Comparison Across Matrices") +
    theme_minimal() +
    labs(x = "Matrix Index", y = "Value") +
    theme(legend.position = "right", legend.title = element_text(size = 12))
  
  # Return both results and plots
  return(list("Jaccard_Heatmap" = jaccard_heatmap, 
              "Metrics_Barplot" = metrics_barplot, 
              "Metrics_Results" = metrics_results))
}
