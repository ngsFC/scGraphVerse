plotROC <- function(matrices_list, ground_truth, plot_title, is_binary = FALSE) {
  truth_vec <- as.vector(ground_truth[upper.tri(ground_truth)])
  
  roc_data <- data.frame()
  auc_values <- data.frame(Matrix = character(), AUC = numeric())
  
  if (is_binary) {
    # For binary matrices, calculate the ROC curve from aggregated points
    binary_points <- data.frame(FPR = numeric(), TPR = numeric())
    for (i in seq_along(matrices_list)) {
      binary_matrix <- as.data.frame(matrices_list[[i]])
      binary_matrix <- binary_matrix[rownames(ground_truth), colnames(ground_truth)]
      pred_vec <- as.vector(as.matrix(binary_matrix)[upper.tri(binary_matrix)])
      
      # Calculate TPR and FPR for this binary matrix
      tp <- sum(pred_vec == 1 & truth_vec == 1)
      fp <- sum(pred_vec == 1 & truth_vec == 0)
      fn <- sum(pred_vec == 0 & truth_vec == 1)
      tn <- sum(pred_vec == 0 & truth_vec == 0)
      
      TPR <- tp / (tp + fn)  # Sensitivity
      FPR <- fp / (fp + tn)  # 1 - Specificity
      
      binary_points <- bind_rows(binary_points, data.frame(
        FPR = FPR,
        TPR = TPR
      ))
    }
    
    # Sort the points by FPR to create a proper ROC curve
    binary_points <- binary_points %>% arrange(FPR)
    
    # Compute the AUC manually from the binary points
    auc_value <- sum(diff(c(0, binary_points$FPR, 1)) * (c(0, binary_points$TPR) + c(binary_points$TPR, 1)) / 2)
    auc_values <- data.frame(Matrix = "Binary Matrices", AUC = auc_value)
    
    # Prepare data for plotting the ROC curve
    roc_data <- data.frame(
      FPR = c(0, binary_points$FPR, 1),
      TPR = c(0, binary_points$TPR, 1),
      Matrix = paste("Binary Matrices (AUC=", round(auc_value, 2), ")", sep = "")
    )
    
  } else {
    # For weighted matrices, calculate the ROC curve as usual
    for (i in seq_along(matrices_list)) {
      weight_matrix <- as.data.frame(matrices_list[[i]])
      weight_matrix <- weight_matrix[rownames(ground_truth), colnames(ground_truth)]
      pred_vec <- as.vector(as.matrix(weight_matrix)[upper.tri(weight_matrix)])
      
      roc_obj <- roc(truth_vec, pred_vec, plot=FALSE)
      auc_value <- round(roc_obj$auc, 2)
      
      # Store AUC values in a data frame
      auc_values <- bind_rows(auc_values, data.frame(Matrix = paste("Weighted Matrix", i), AUC = auc_value))
      
      roc_df <- data.frame(
        FPR = 1 - roc_obj$specificities, 
        TPR = roc_obj$sensitivities,
        Matrix = paste("Weighted Matrix", i, "(AUC=", auc_value, ")", sep = "")
      )
      
      roc_data <- bind_rows(roc_data, roc_df)
    }
  }
  
  # Generate dynamic colors for all matrices
  total_matrices <- length(unique(roc_data$Matrix))
  colors <- hue_pal()(total_matrices)
  
  # Plot the ROC curve
  p <- ggplot() +
    labs(title=plot_title,
         x="False Positive Rate (1 - Specificity)",
         y="True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey")
  
  # Add lines for weighted matrices or binary matrix ROC curve
  if (!is_binary) {
    p <- p +
      geom_line(data=roc_data, aes(x=FPR, y=TPR, color=Matrix), size=1.2)
  } else {
    # Add the line connecting binary dots and draw the ROC curve
    p <- p +
      geom_line(data=roc_data, aes(x=FPR, y=TPR, color=Matrix), size=1.2) +
      geom_point(data=binary_points, aes(x=FPR, y=TPR), size=2, color="blue")
  }
  
  # Add dynamic color mapping
  p <- p + scale_color_manual(values=colors)
  
  print(p)
  
  return(auc_values)
}
