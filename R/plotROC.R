#' Plot ROC Curves for Inferred Networks
#'
#' This function computes and plots Receiver Operating Characteristic (ROC) curves 
#' for a list of predicted adjacency matrices against a binary ground truth network. 
#' It supports both binary (0/1) and weighted prediction matrices.
#'
#' @param matrices_list A list of square matrices representing predicted interactions. 
#'   Each matrix must have the same dimensions and row/column names as `ground_truth`. 
#'   Entries may be binary (0/1) or continuous weights.
#' @param ground_truth A square binary matrix indicating true interactions (1) 
#'   in the upper triangle. Must have the same dimensions and names as each matrix in `matrices_list`.
#' @param plot_title Character. Title of the ROC plot.
#' @param is_binary Logical. If `TRUE`, matrices are treated as binary predictions. 
#'   Defaults to `FALSE` for weighted prediction matrices.
#'
#' @return A data frame containing the Area Under the Curve (AUC) for each matrix.
#'
#' @details 
#' For binary matrices, the ROC curve is computed using a single TPR/FPR point per matrix. 
#' For weighted matrices, the full ROC curve is calculated using the continuous prediction scores. 
#' The resulting plot includes one ROC curve per matrix, with AUC values shown in the legend.
#'
#' Diagonal elements are ignored in all comparisons. The function assumes symmetry is not required.
#'
#' @examples
#' \dontrun{
#' mat1 <- matrix(runif(100), nrow = 10)
#' mat2 <- matrix(runif(100), nrow = 10)
#' ground_truth <- matrix(sample(c(0, 1), 100, replace = TRUE), nrow = 10)
#' diag(ground_truth) <- 0
#' ground_truth[lower.tri(ground_truth)] <- 0
#'
#' plotROC(matrices_list = list(mat1, mat2), 
#'         ground_truth = ground_truth,
#'         plot_title = "ROC for Network Inference", 
#'         is_binary = FALSE)
#' }
#'
#' @import ggplot2
#' @importFrom pROC roc
#' @importFrom scales hue_pal
#' @importFrom dplyr bind_rows arrange
#' @export

plotROC <- function(matrices_list, ground_truth, plot_title, is_binary = FALSE) {
  
  # Convert the ground truth matrix to a vector of values for the upper triangle
  truth_vec <- as.vector(ground_truth[upper.tri(ground_truth)])
  
  # Initialize empty data frames to store ROC data and AUC values
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
      
      binary_points <- dplyr::bind_rows(binary_points, data.frame(
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
      
      # Compute ROC for the weighted matrix
      roc_obj <- roc(truth_vec, pred_vec, plot=FALSE)
      auc_value <- round(roc_obj$auc, 2)
      
      # Store AUC values in a data frame
      auc_values <- dplyr::bind_rows(auc_values, data.frame(Matrix = paste("Weighted Matrix", i), AUC = auc_value))
      
      # Prepare data for plotting the ROC curve
      roc_df <- data.frame(
        FPR = 1 - roc_obj$specificities, 
        TPR = roc_obj$sensitivities,
        Matrix = paste("Weighted Matrix", i, "(AUC=", auc_value, ")", sep = "")
      )
      
      roc_data <- dplyr::bind_rows(roc_data, roc_df)
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
    p <- p + geom_line(data=roc_data, aes(x=FPR, y=TPR, color=Matrix), size=1.2)
  } else {
    # Add the line connecting binary dots and draw the ROC curve
    p <- p + geom_line(data=roc_data, aes(x=FPR, y=TPR, color=Matrix), size=1.2) +
      geom_point(data=binary_points, aes(x=FPR, y=TPR), size=2, color="blue")
  }
  
  # Add dynamic color mapping
  p <- p + scale_color_manual(values=colors)
  
  # Print the plot
  print(p)
  
  return(auc_values)
}

