plotROC <- function(weighted_matrices_list, ground_truth, plot_title) {
  
  # Convert the upper triangle of the ground truth matrix to a vector (excluding the diagonal)
  truth_vec <- as.vector(ground_truth[upper.tri(ground_truth)])
  
  roc_data <- data.frame()
  auc_values <- vector("character", length(weighted_matrices_list))
  
  # Loop through each weighted adjacency matrix in the input list
  for (i in seq_along(weighted_matrices_list)) {
    # Ensure the matrix dimensions match the ground truth
    weight_matrix <- as.data.frame(weighted_matrices_list[[i]])
    weight_matrix <- weight_matrix[rownames(ground_truth), colnames(ground_truth)]
    
    # Convert the upper triangle of the matrix to a vector (excluding the diagonal)
    pred_vec <- as.vector(as.matrix(weight_matrix)[upper.tri(weight_matrix)])
    
    # Calculate the ROC curve
    roc_obj <- roc(truth_vec, pred_vec, plot=FALSE)
    
    # Extract TPR and FPR values for ggplot and add them to the data frame
    roc_df <- data.frame(
      FPR = rev(roc_obj$specificities), 
      TPR = rev(roc_obj$sensitivities),
      Matrix = paste("Matrix", i)
    )
    
    # Append the data for each matrix
    roc_data <- bind_rows(roc_data, roc_df)
    
    # Store AUC values for the legend
    auc_values[i] <- paste("Matrix", i, "(AUC=", round(roc_obj$auc, 2), ")", sep = "")
  }
  
  color_count <- length(unique(roc_data$Matrix))
  colors <- brewer.pal(min(color_count, 9), "Set1")
  
  ggplot(roc_data, aes(x=FPR, y=TPR, color=Matrix)) +
    geom_line(size=1.2) +
    labs(title=plot_title,
         x="False Positive Rate (1 - Specificity)",
         y="True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
      legend.position = "bottom") +
    scale_color_manual(values=colors, labels=auc_values) +  
    scale_x_reverse()
}
