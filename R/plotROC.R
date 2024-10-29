# Function to calculate and plot ROC curves for multiple weighted adjacency matrices using ggplot
plotROC <- function(weighted_matrices_list, ground_truth, plot_title) {
  
  # Convert the ground truth matrix to a vector
  truth_vec <- as.vector(ground_truth)
  
  # Initialize a data frame to store ROC data for ggplot
  roc_data <- data.frame()
  auc_values <- vector("character", length(weighted_matrices_list))
  
  # Loop through each weighted adjacency matrix in the input list
  for (i in seq_along(weighted_matrices_list)) {
    # Ensure the matrix dimensions match the ground truth
    weight_matrix <- as.data.frame(weighted_matrices_list[[i]])
    weight_matrix <- weight_matrix[rownames(ground_truth), colnames(ground_truth)]
    
    # Convert the matrix to a vector
    pred_vec <- as.vector(as.matrix(weight_matrix))
    
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
  
  # Generate a dynamic color palette for an arbitrary number of matrices
  color_count <- length(unique(roc_data$Matrix))
  colors <- brewer.pal(min(color_count, 9), "Set1")
  
  # Plotting with ggplot2
  ggplot(roc_data, aes(x=FPR, y=TPR, color=Matrix)) +
    geom_line(size=1.2) +
    labs(title=plot_title,
         x="False Positive Rate (1 - Specificity)",
         y="True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"  # Move legend to the bottom
    ) +
    scale_color_manual(values=colors, labels=auc_values) +  # Display AUC values in the legend
    scale_x_reverse()  # Reverse the x-axis
}