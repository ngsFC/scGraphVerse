plotROC <- function(wlist, gtruth, plot_title = "ROC Curve") {
  
  roc_data <- data.frame()
  auc_values <- data.frame(Matrix = character(), AUC = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(wlist)) {
    # For each element in the list, create the ground truth
    roclist <- wlist
    roclist[[i]]$ground_truth <- mapply(function(reg, tgt) {
      return(ifelse(gtruth[reg, tgt] == 1, 1, 0))
    }, roclist[[i]]$Gene1, roclist[[i]]$Gene2)
    
    # Generate ROC curve
    roc_curve <- roc(roclist[[i]]$ground_truth, roclist[[i]]$weight)
    
    # Create data for plotting the ROC curve
    roc_df <- data.frame(
      FPR = 1 - roc_curve$specificities,
      TPR = roc_curve$sensitivities,
      Matrix = paste("Matrix", i)
    )
    
    # Store the AUC values
    auc_values <- rbind(auc_values, data.frame(Matrix = paste("Matrix", i), AUC = round(auc(roc_curve), 2)))
    
    # Append to the full ROC data
    roc_data <- rbind(roc_data, roc_df)
  }
  
  # Create the ROC plot
  roc_plot <- ggplot(roc_data, aes(x = 1 - FPR, y = TPR, color = Matrix)) +
    geom_line(size = 1) +
    labs(title = plot_title) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    scale_x_reverse() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = " "))
  
  # Add AUC values as text to the plot
  roc_plot <- roc_plot + 
    geom_text(data = auc_values, aes(x = 0.10, y = 0.2 - 0.05 * as.numeric(gsub("Matrix ", "", Matrix)),
                                     label = paste(Matrix, "- AUC:", AUC)),
              color = "black", size = 4, show.legend = FALSE)
  
  return(roc_plot)
}
