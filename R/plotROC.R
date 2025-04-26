#' Plot ROC Curves for Inferred Networks
#'
#' Computes and visualizes Receiver Operating Characteristic (ROC) curves
#' for a list of predicted adjacency matrices compared against a binary ground truth network.
#'
#' @param matrices_list A list of square matrices representing predicted interactions.
#'   Each matrix must have the same dimensions and row/column names as \code{ground_truth}.
#'   Entries may be binary (0/1) or continuous weights.
#' @param ground_truth A square binary matrix indicating true interactions (1) in the upper triangle.
#'   Must match the dimensions and names of each matrix in \code{matrices_list}.
#' @param plot_title Character string. Title for the ROC plot.
#' @param is_binary Logical. If \code{TRUE}, matrices are treated as binary predictions.
#'   Default is \code{FALSE}, for weighted predictions.
#'
#' @return
#' A data frame containing the Area Under the Curve (AUC) for each matrix.
#' The ROC plot is displayed automatically.
#'
#' @details
#' For binary matrices, a single TPR/FPR point is computed per matrix.
#' For weighted matrices, a full ROC curve is calculated based on continuous prediction scores.
#' Diagonal entries are ignored. Symmetry between rows and columns is not enforced.
#'
#' @import ggplot2
#' @importFrom pROC roc
#' @importFrom scales hue_pal
#' @importFrom dplyr bind_rows arrange
#' @export
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

plotROC <- function(matrices_list, ground_truth, plot_title, is_binary = FALSE) {
  
  truth_vec <- as.vector(ground_truth[upper.tri(ground_truth)])
  
  roc_data <- data.frame()
  auc_values <- data.frame(Matrix = character(), AUC = numeric())
  
  if (is_binary) {
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
    
    binary_points <- binary_points %>% arrange(FPR)
    auc_value <- sum(diff(c(0, binary_points$FPR, 1)) * (c(0, binary_points$TPR) + c(binary_points$TPR, 1)) / 2)
    auc_values <- data.frame(Matrix = "Binary Matrices", AUC = auc_value)
    
    roc_data <- data.frame(
      FPR = c(0, binary_points$FPR, 1),
      TPR = c(0, binary_points$TPR, 1),
      Matrix = paste("Binary Matrices (AUC=", round(auc_value, 2), ")", sep = "")
    )
    
  } else {
    for (i in seq_along(matrices_list)) {
      weight_matrix <- as.data.frame(matrices_list[[i]])
      weight_matrix <- weight_matrix[rownames(ground_truth), colnames(ground_truth)]
      pred_vec <- as.vector(as.matrix(weight_matrix)[upper.tri(weight_matrix)])
      
      roc_obj <- roc(truth_vec, pred_vec, plot=FALSE)
      auc_value <- round(roc_obj$auc, 2)
      
      auc_values <- dplyr::bind_rows(auc_values, data.frame(Matrix = paste("Weighted Matrix", i), AUC = auc_value))
      
      roc_df <- data.frame(
        FPR = 1 - roc_obj$specificities, 
        TPR = roc_obj$sensitivities,
        Matrix = paste("Weighted Matrix", i, "(AUC=", auc_value, ")", sep = "")
      )
      
      roc_data <- dplyr::bind_rows(roc_data, roc_df)
    }
  }
  
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
  
  if (!is_binary) {
    p <- p + geom_line(data=roc_data, aes(x=FPR, y=TPR, color=Matrix), size=1.2)
  } else {
    p <- p + geom_line(data=roc_data, aes(x=FPR, y=TPR, color=Matrix), size=1.2) +
      geom_point(data=binary_points, aes(x=FPR, y=TPR), size=2, color="blue")
  }
  
  p <- p + scale_color_manual(values=colors)
  
  # Print the plot
  #print(p)
  
  return(auc_values)
}

