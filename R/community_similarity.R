community_similarity <- function(control_output, predicted_list) {
  required_pkgs <- c("igraph", "fmsb")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), 
         "\nInstall them using BiocManager::install() or install.packages().")
  }
  
  library(igraph)
  library(fmsb)
  
  # Extract control community membership
  control_comm <- control_output$communities$membership
  
  # Initialize results dataframe
  metrics_list <- list()
  
  # Compute similarity metrics for each predicted community structure
  for (i in seq_along(predicted_list)) {
    predicted_comm <- predicted_list[[i]]$communities$membership
    
    vi_value <- compare(control_comm, predicted_comm, method = "vi")  # ✅ Variation of Information (VI)
    nmi_value <- compare(control_comm, predicted_comm, method = "nmi")  # ✅ Normalized Mutual Information (NMI)
    ari_value <- compare(control_comm, predicted_comm, method = "adjusted.rand")  # ✅ Adjusted Rand Index (ARI)
    split_join_value <- compare(control_comm, predicted_comm, method = "split.join")  # ✅ Split-Join Distance
    
    metrics_list[[paste0("Predicted_", i)]] <- c(vi_value, nmi_value, ari_value, split_join_value)
  }
  
  # Convert list to dataframe
  metrics_df <- as.data.frame(do.call(rbind, metrics_list))
  colnames(metrics_df) <- c("VI", "NMI", "ARI", "Split-Join")
  rownames(metrics_df) <- names(metrics_list)
  
  # Normalize values for radar chart (except NMI and ARI which are naturally 0-1)
  normalized_df <- metrics_df
  normalized_df$VI <- (max(metrics_df$VI) - metrics_df$VI) / max(metrics_df$VI)  # VI (inverted scale)
  normalized_df$`Split-Join` <- (max(metrics_df$`Split-Join`) - metrics_df$`Split-Join`) / max(metrics_df$`Split-Join`)  # Split-Join (inverted scale)
  
  # Prepare data for radar chart
  max_vals <- rep(1, ncol(normalized_df))
  min_vals <- rep(0, ncol(normalized_df))
  radar_data <- rbind(max_vals, min_vals, normalized_df)
  
  # Define colors for multiple lines
  colors <- rainbow(nrow(metrics_df))
  
  # Plot radar chart with multiple lines
  par(mar = c(2, 2, 2, 2))
  radarchart(
    radar_data,
    axistype = 2,
    pcol = colors,
    pfcol = NA,
    plwd = 2,
    plty = 1,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = seq(0, 1, 0.2),
    vlcex = 1.1
  )
  
  # Add legend for multiple predictions
  legend("topright", legend = rownames(metrics_df), col = colors, lty = 1, lwd = 2)
  
  # Return dataframe with metrics
  return(metrics_df)
}
